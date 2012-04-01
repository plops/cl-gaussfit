;; the paper I read for this is 1990bentley_kdtree-c++.pdf and
;; friedman 1977, a kd-tree works for other distance measures than
;; euklidean (as opposed to delaunay)
#.(require :alexandria)
(declaim (optimize (speed 0) (safety 3) (debug 3)))

(defpackage :kdtree
  (:use :cl :alexandria)
  (:export #:build-new-tree
	   #:nearest-neighbour
	   #:kd-tree-points
	   #:perm
	   #:get-tree-point
	   #:locate-points-in-circle-around-target))

(in-package :kdtree)

(defconstant +dim+ 2) ;; you can set this to three but then the postscript output doesn't work
(declaim (type fixnum +dim+))

(deftype vec ()
  `(simple-array single-float (#.+dim+)))

(deftype axis ()
  `(member ,@(loop for i below +dim+ collect i)))

(deftype array-index-t ()
  `(unsigned-byte 29))



(defstruct node
  (cutdim (required-argument :cutdim) :type axis)
  (cutval (required-argument :cutval) :type single-float)
  (loson (required-argument :loson) :type (or node leaf))
  (hison (required-argument :hison) :type (or node leaf))
  (father nil :type (or null node))
  (bounds- (make-array +dim+ :element-type 'single-float)
	   :type (simple-array single-float 1))
  (bounds+ (make-array +dim+ :element-type 'single-float)
	   :type (simple-array single-float 1)))

(defstruct leaf
  (lopt (required-argument :lopt) :type array-index-t)
  (hipt (required-argument :hipt) :type array-index-t)
  (father nil :type (or null node)))

(defstruct kd-tree
  (root nil :type (or null node))
  (points (required-argument :points) :type (simple-array vec 1))
  (bucketptr nil
	     :type (or null (simple-array leaf 1)))
  (perm (required-argument :perm) :type (simple-array array-index-t 1)))

(defmacro px (i j)
  ;; i selects point and j the coordinate
  `(aref (aref points (aref perm ,i)) ,j))

(defun get-tree-point (kd-tree i)
  (declare (type kd-tree kd-tree)
	   (type fixnum i)
	   (values vec &optional))
  (with-slots (points
	       perm)
      kd-tree
    (aref points (aref perm i))))

(defun max-spread-in-dim (tree dim l u)
  (declare (type kd-tree tree)
	   (type array-index-t l u)
	   (type axis dim)
	   (values single-float &optional))
  (with-slots (points perm) tree
   (let ((mi (px l dim))
	 (ma (px l dim)))
     (loop for i from (1+ l) upto u do
	  (let ((v (px i dim)))
	    (when (< v mi)
	      (setf mi v))
	    (when (< ma v)
	      (setf ma v))))
     (- ma mi))))

(defun find-max-spread-dimension (tree l u)
  (declare (type kd-tree tree)
	   (type array-index-t l u)
	   (values axis &optional))
  (let ((spread 0f0)
	(dim (the axis 0)))
    (declare (single-float spread))
    (dotimes (j +dim+)
      (let ((v (max-spread-in-dim tree j l u)))
	(when (< spread v)
	  (setf spread v
		dim j))))
    dim))

(defun select (tree l u m cutdim)
  (declare (type kd-tree tree)
	   (type array-index-t l u m)
	   (type axis cutdim)
	   (values fixnum &optional))
  (with-slots (points perm) tree
    (labels ((p (l) ;; accessor
	       (declare (fixnum l)
			(values single-float &optional))
	       (px l cutdim))
	     (swap (a b)
	       (declare (type array-index-t a b))
	       (with-slots (perm) tree
		 (rotatef (aref perm a) (aref perm b)))
	       (values)))
     (loop
	(if (<= u (1+ l))
	    (progn
	      (when (and (= u (1+ l))
			 (< (p u) (p l)))
		(swap u l))
	      (return-from select (aref perm m)))
	    (let ((mid (floor (+ u l) 2)))
	      (swap mid (1+ l))
	      (labels ((enforce< (l u)
			 (declare (type fixnum l u))
			 (when (< (p u) (p l)) (swap u l))
			 (values)))
		;; rearrange for p[l]<=p[l+1]<=p[u]
		(enforce< l u)
		(enforce< (1+ l) u)
		(enforce< l (1+ l)))
	      (let* ((i (1+ l))
		     (j u)
		     (ia (aref perm i))
		     (a (p i)))
		(declare (type fixnum i j))
		(loop ;; scan up and down, put big elements right, small
		   ;; ones left
		   (loop do (incf i) while (< (p i) a))
		   (loop do (decf j) while (< a (p j)))
		   (when (< j i) (return))
		   (swap i j))
		(setf (aref perm (1+ l)) (aref perm j)
		      (aref perm j) ia)
		(when (<= m j) (setf u (1- j)))
		(when (<= j m) (setf l i)))))))))

(defun build (tree l u)
  (declare (type kd-tree tree)
	   (type array-index-t l u)
	   (values (or node leaf) &optional))
  (with-slots (points perm) tree
   (let ((points-in-bucket 14)) ;; change this back to somewhere like 14 for better performance
     (if (<= (1+ (- u l)) (1- points-in-bucket))
	 (make-leaf :lopt l
		    :hipt u)
	 (let ((cutdim (find-max-spread-dimension tree l u))
	       (m (floor (+ l u) 2)))
	   (select tree l u m cutdim)
	   (make-node :cutdim cutdim 
		      :cutval (px m cutdim)
		      :loson (build tree l m)
		      ;; the next select will only shuffle inside l .. m
		      :hison (build tree (1+ m) u)))))))

(defun fill-bucketptr (tree)
  (declare (type kd-tree tree))
  (with-slots (root points perm) tree
   (let* ((n (length points))
	  (noleaf (make-leaf :lopt 0 :hipt 0))
	  (bp (make-array n :element-type 'leaf
			  :initial-element noleaf)))
     (labels ((rec (node)
		(declare (type (or leaf node) node))
		(etypecase node
		  (leaf (with-slots (lopt hipt) node
			  (loop for i from lopt upto hipt do
			       (setf (aref bp (aref perm i)) node))
			  (return-from rec nil)))
		  (node (with-slots (cutval cutdim loson hison) node
			  (rec loson)
			  (rec hison))))))
       (rec root))
     bp)))

(defun fill-fathers (tree)
  (declare (type kd-tree tree))
  (with-slots (root points perm) tree
    (labels ((rec (node father)
	       (declare (type (or leaf node) node)
			(type (or null node) father))
	       (etypecase node
		 (leaf (setf (leaf-father node) father))
		 (node (setf (node-father node) father)
		       (with-slots (loson hison) node
			 (rec loson node)
			 (rec hison node))))))
      (rec root nil))
    tree))

(defun fill-bounds (tree)
  (declare (type kd-tree tree))
  (with-slots (root points perm) tree
    (labels ((make-bounds (ls)
	       (make-array +dim+ :element-type 'single-float
			   :initial-contents ls))
	     (rec (node b- b+)
	       (declare (type (or leaf node) node)
			(type (simple-array single-float 1) b- b+))
	       (etypecase node
		 (leaf nil)
		 (node (setf (node-bounds+ node) b+
			     (node-bounds- node) b-)
		       (with-slots (loson hison) node
			 (rec loson
			      b- 
			      (make-bounds 
			       (loop for k below +dim+ collect
				    (if (= k (node-cutdim node))
					(node-cutval node)
					(aref (node-bounds+ node) k)))))
			 (rec hison
			      (make-bounds 
			       (loop for k below +dim+ collect
				    (if (= k (node-cutdim node))
					(node-cutval node)
					(aref (node-bounds- node) k))))
			      b+))))))
      (rec root 
	   (make-bounds (loop for k below +dim+ collect -1f20))
	   (make-bounds (loop for k below +dim+ collect 1f20))))
    tree))

(defun build-new-tree (points)
  (declare (type (simple-array vec 1) points)
	   (values kd-tree &optional))
  (let* ((n (length points))
	 (perm (make-array n
			   :element-type 'array-index-t
			     :initial-contents 
			     (loop for i below n collect i)))
	 (tree (make-kd-tree :points points
			     :perm perm
			     :root nil))
	 (root (build tree 0 (1- n))))
    (setf (kd-tree-root tree) root
	  (kd-tree-bucketptr tree) (fill-bucketptr tree))
    (fill-bounds tree)
    (fill-fathers tree) 
    tree))

(defun distance2 (tree i j)
  (declare (type kd-tree tree)
	   (type array-index-t i j)
	   (values single-float &optional))
  (with-slots (points perm) tree
   (let ((sum 0f0))
     (declare (type (single-float 0f0) sum))
     (dotimes (k +dim+)
       (let* ((v (- (px i k) (px j k)))
	      (v2 (* v v)))
	 (declare (type (single-float 0f0) v2))
	 (incf sum v2)))
     sum)))

(defmacro with-rnn (&body body)
  ;; nndist2, nntarget and nnptnum are defined outside
  `(labels ((rnn (node) 
	      (declare (type (or leaf node) node))
	      (etypecase node
		;; sequential search in the bucket of the leaf
		(leaf (with-slots (lopt hipt)
			  node
			(loop for i from lopt upto hipt do
			 (let ((d (distance2 tree i nntarget)))
			   (when (< d nndist2)
			     (setf nndist2 d
				   nnptnum (aref perm i)))))
		    (return-from rnn nil)))
	   	(node (with-slots (cutval cutdim loson hison)
       node
     (let* ((x (px nntarget cutdim))
	    (diff (- x cutval))
	    (diff2 (expt diff 2)))
       (if (< diff 0f0)
	   (progn (rnn loson)
		  ;; only check on the other side, when
		  ;; ball of current nndist centered at
		  ;; query record overlaps into the
		  ;; other side
		  (when (< diff2 nndist2)
		    (rnn hison)))
	    ;; search in the closer son, if nothing there,
	    ;; search in farther son as well
	   (progn (rnn hison)
		  (when (< diff2 nndist2)
		    (rnn loson))))))))))
     ,@body))

(defun nearest-neighbour-top-down (tree target)
  (declare (type array-index-t target)
	   (type kd-tree tree)
	   (values array-index-t single-float &optional))
  (with-slots (perm points root) tree
    (let ((nntarget target)
	  (nndist2 1f20)
	  (nnptnum 0))
      (with-rnn
	(rnn root))
      (values nnptnum nndist2))))

(defun nearest-neighbour-bottom-up (tree target)
  (declare (type array-index-t target)
	   (type kd-tree tree)
	   (values array-index-t single-float &optional))
  (let ((nntarget target)
	(nndist2 1f20)
	(nnptnum 0))
    (with-slots (perm points root bucketptr) tree
      (labels ((ball-in-bounds (node targ rad)
		 (declare (type node node)
			  (type array-index-t targ)
			  (type single-float rad))
		 (with-slots (bounds- bounds+) node
		   (loop for k below +dim+ do
			(let ((x (px targ k)))
			  (when (or (<= (- x (aref bounds- k)) rad)
				    (<= (- (aref bounds+ k) x) rad))
			    (return-from ball-in-bounds nil))))
		   t)))
	(with-rnn
	 (let* ((p (aref bucketptr target))
		(old-p p))
	   (rnn p)
	   (loop named trav do
		(setf old-p p
		      p (slot-value p 'father))
		(unless p
		  (return-from trav))
		(with-slots (cutdim cutval loson hison) p
		  (let ((diff2 (expt (- (px target cutdim)
					cutval)
				     2)))
		    (when (<= diff2 nndist2)
		      (if (eq old-p loson)
			  (rnn hison)
			  (rnn loson)))
		    (when (ball-in-bounds p target nndist2)
		      (return-from trav)))))
	   (values nnptnum nndist2)))))))



(defun locate-points-in-circle-around-target (target radius tree)
  (declare (type array-index-t target)
	   (type single-float radius)
	   (type kd-tree tree))
  (with-slots (perm points root bucketptr) tree
   (let* ((nndist radius)
	  (nndist2 (expt nndist 2))
	  (nntarget target)
	  (res nil))
     (labels ((ball-in-bounds (node targ rad)
		(declare (type node node)
			 (type array-index-t targ)
			 (type single-float rad))
		(with-slots (bounds- bounds+) node
		  (loop for k below +dim+ do
		       (let ((x (px targ k)))
			 (when (or (<= (- x (aref bounds- k)) rad)
				   (<= (- (aref bounds+ k) x) rad))
			   (return-from ball-in-bounds nil))))
		  t)) 
	      (rfrnn (node)
		(declare (type (or leaf node) node))
		(etypecase node
		  (leaf (with-slots (lopt hipt) node
			  (loop for i from lopt upto hipt do
			       (when (< (distance2 tree i nntarget)
					nndist2)
				 (push (aref perm i) res)))
			  (return-from rfrnn nil)))
		  (node (with-slots (cutval cutdim loson hison) node
			  (let* ((x (px nntarget cutdim))
				 (diff (- x cutval)))
			    (if (< diff 0f0)
				(progn (rfrnn loson)
				       (when (<= (- diff) nndist) 
					 (rfrnn hison)))
				(progn (rfrnn hison)
				       (when (<= diff nndist)
					 (rfrnn loson))))))))))
       (let* ((p (aref bucketptr nntarget))
	      (old-p p))
	 (rfrnn p)
	 (loop named trav do
	      (setf old-p p
		    p (slot-value p 'father))
	      (unless p
		(return-from trav))
	      (with-slots (cutdim cutval loson hison) p
		(let ((diff (- (px target cutdim)
			       cutval)))
		  (if (eq old-p loson)
		      (when (<= (- diff) nndist)
			(rfrnn hison))
		      (when (<= diff nndist)
			(rfrnn loson)))
		    (when (ball-in-bounds p target nndist2)
		      (return-from trav))))))
       (reverse res)))))
#+nil
(let* ((n 30000))
  (time (defparameter *tree* (build-new-tree (make-random-points n))))
  (time (dotimes (i n) (nearest-neighbour i *tree*))))

;; 10s to find all nn for 30000 points when 2, 9 or 14 points in bucket
;; 11.5s for 20 points in bucket
;; 14s for 50 points in bucket

(progn
;; for debugging I draw the points and the rectangles into an eps
;; file.  There is also a function that writes the tree in a format
;; that the program dot (from the graphviz) package can render in a
;; postscript image of the trees graph structure.

(defun readable-float (x)
  "Prints only a few digits of a float X into a string, to make the
output more readable."
  (read-from-string 
   (format nil "~2,4f" x)))


(defun draw-tree-boxes (root)
  (declare (node root))
  (let ((boxes nil))
    (labels ((rec (node)
	      (when (typep node 'node)
		(with-slots (loson hison bounds- bounds+) node
		  (when (or (typep loson 'leaf)
			    (typep hison 'leaf))
		    (push (list bounds- bounds+) boxes))
		  (rec loson)
		  (rec hison)))))
      (rec root))
    (reverse boxes)))

#+nil
(defparameter *boxes*
  (draw-tree-boxes (kd-tree-root gauss::*tree*)))

(defun eps-moveto (x y)
  (declare (type single-float x y))
  (format nil "~f ~f moveto~%" x y))

(defun eps-lineto (x y)
  (declare (type single-float x y))
  (format nil "~f ~f lineto~%" x y))

(defun eps-rectangle (b- b+)
  (declare (type (simple-array single-float 1) b- b+))
  (let ((x0 (aref b- 0))
	(y0 (aref b- 1))
	(x (aref b+ 0))
	(y (aref b+ 1)))
    (format nil "newpath~%~a~a~a~a~astroke~%"
	   (eps-moveto x0 y0)
	   (eps-lineto x y0)
	   (eps-lineto x y)
	   (eps-lineto x0 y)
	   (eps-lineto x0 y0))))

(defun eps-point (x y)
  (declare (single-float x y))
  (format nil "newpath ~f ~f 0.002 0 360 arc closepath fill~%" x y))

(defun eps-tree (fn boxes points)
  (declare (list boxes)
	   ((simple-array vec 1) points)
	   (string fn))
  (with-open-file (s fn :direction :output
		     :if-exists :supersede)
    (format s "%!PS-Adobe-3.0
%%Pages: 1
%%BoundingBox: 0 0 700 500
%%EndComments
10 10 translate
8.0 8.0 scale
0.002 setlinewidth
0 setgray~%")
   #+nil 
   (loop for (b- b+) in boxes do
	(format s "~a" (eps-rectangle b- b+)))
  
   (loop for p in (locate-points-in-circle-around-target
		   31100 3f0 gauss::*tree*) do
	(format s "~a" (eps-point (aref (aref points p) 0)
				  (aref (aref points p) 1))))
   #+nil
   (loop for p across points do
	(format s "~a" (eps-point (aref p 0) (aref p 1))))
   (format s "%%EOF"))))

#+nil
(eps-tree "/dev/shm/o.ps"
	  (draw-tree-boxes (kd-tree-root gauss::*tree*)) 
	  (kd-tree-points gauss::*tree*))
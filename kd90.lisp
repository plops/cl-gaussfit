;; the paper I read for this is 1990bentley_kdtree-c++.pdf
#.(require :alexandria)
(declaim (optimize (speed 3) (safety 3) (debug 3)))

(defpackage :kdtree
  (:use :cl :alexandria)
  (:export #:build-new-tree
	   #:nearest-neighbour
	   #:kd-tree-points
	   #:perm
	   #:get-tree-point))

(in-package :kdtree)

(defconstant +dim+ 2) ;; you can set this to three but then the postscript output doesn't work
(declaim (type fixnum +dim+))

(deftype vec ()
  `(simple-array single-float (#.+dim+)))

(deftype axis ()
  `(member ,@(loop for i below +dim+ collect i)))

(deftype array-index-t ()
  `(unsigned-byte 29))

(defstruct leaf
  (lopt (required-argument :lopt) :type array-index-t)
  (hipt (required-argument :hipt) :type array-index-t))

(defstruct node
  (cutdim (required-argument :cutdim) :type axis)
  (cutval (required-argument :cutval) :type single-float)
  (father (required-argument :father) :type node)
  (loson (required-argument :loson) :type (or node leaf))
  (hison (required-argument :hison) :type (or node leaf)))

(defstruct kd-tree
  (root nil :type (or null node))
  (points (required-argument :points) :type (simple-array vec 1))
  (bucketptr (required-argument :bucketptr) :type (simple-array leaf 1))
  (perm (required-argument :perm) :type (simple-array array-index-t 1)))

(defparameter *points* (make-array 0 :element-type 'vec))
(defparameter *perm* (make-array 0 :element-type 'array-index-t))
(declaim (type (simple-array vec 1) *points*)
	 (type (simple-array array-index-t 1) *perm*))
(declaim (inline px))
(defun px (i j)
  (declare (type array-index-t i j)
	   (values single-float &optional)) 
  ;; i selects point and j the coordinate
  (aref (aref *points* (aref *perm* i)) j))

(defun get-tree-point (kd-tree i)
  (declare (type kd-tree kd-tree)
	   (type fixnum i)
	   (values vec &optional))
  (with-slots (points
	       perm)
      kd-tree
    (aref points (aref perm i))))

(defun max-spread-in-dim (dim l u)
  (declare (type array-index-t l u)
	   (type axis dim)
	   (values single-float &optional))
  (let ((mi (px l dim))
	(ma (px l dim)))
    (loop for i from (1+ l) upto u do
	 (let ((v (px i dim)))
	   (when (< v mi)
	     (setf mi v))
	   (when (< ma v)
	     (setf ma v))))
    (- ma mi)))

(defun find-max-spread-dimension (l u)
  (declare (type array-index-t l u)
	   (values axis &optional))
  (let ((spread 0f0)
	(dim (the axis 0)))
    (declare (single-float spread))
    (dotimes (j +dim+)
      (let ((v (max-spread-in-dim j l u)))
	(when (< spread v)
	  (setf spread v
		dim j))))
    dim))

(defun swap (a b)
  (declare (type array-index-t a b)
	   (values null &optional))
  (let ((h (aref *perm* a)))
    (setf (aref *perm* a) (aref *perm* b)
	  (aref *perm* b) h))
  nil)

(defun select (l u m cutdim)
  (declare (type array-index-t l u m)
	   (type axis cutdim)
	   (values fixnum &optional))
  (labels ((p (l) ;; accessor
	     (declare (fixnum l)
		      (values single-float &optional))
	     (px l cutdim)))
     (loop
       (if (<= u (1+ l))
	   (progn
	     (when (and (= u (1+ l))
			(< (p u) (p l)))
	       (swap u l))
	     (return-from select (aref *perm* m)))
	   (let ((mid (floor (+ u l) 2)))
	     (swap mid (1+ l))
	     (labels ((enforce< (l u)
			(declare (type fixnum l u)
				 (values null &optional))
			(when (< (p u) (p l)) (swap u l))
			nil))
	       ;; rearrange for p[l]<=p[l+1]<=p[u]
	       (enforce< l u)
	       (enforce< (1+ l) u)
	       (enforce< l (1+ l)))
	     (let* ((i (1+ l))
		    (j u)
		    (ia (aref *perm* i))
		    (a (p i)))
	       (declare (type fixnum i j))
	       (loop ;; scan up and down, put big elements right, small
		  ;; ones left
		  (loop do (incf i) while (< (p i) a))
		  (loop do (decf j) while (< a (p j)))
		  (when (< j i) (return))
		  (swap i j))
	       (setf (aref *perm* (1+ l)) (aref *perm* j)
		     (aref *perm* j) ia)
	       (when (<= m j) (setf u (1- j)))
	       (when (<= j m) (setf l i))))))))

(defun build (l u)
  (declare (type array-index-t l u)
	   (values (or node leaf) &optional))
  (let ((points-in-bucket 2)) ;; change this back to somewhere like 14 for better performance
   (if (<= (1+ (- u l)) (1- points-in-bucket))
       (make-leaf :lopt l
		  :hipt u)
       (let ((cutdim (find-max-spread-dimension l u))
	     (m (floor (+ l u) 2)))
	 (select l u m cutdim)
	 (make-node :cutdim cutdim 
		    :cutval (px m cutdim)
		    :loson (build l m) ;; the next select will only shuffle inside l .. m
		    :hison (build (1+ m) u))))))

(defun build-new-tree (points)
  (declare (type (simple-array vec 1) points)
	   (values kd-tree &optional))
  (let* ((n (length points))
	 (*points* points)
	 (*perm* (make-array n
			     :element-type 'array-index-t
			     :initial-contents 
			     (loop for i below n collect i)))
	 (root (build 0 (1- n))))
    (make-kd-tree :points *points*
		  :perm *perm*
		  :root root)))
#+nil
(let* ((n 300))
  (defparameter *tree* (build-new-tree (make-random-points n))))

(defun distance2 (i j)
  (declare (type array-index-t i j)
	   (values single-float &optional))
  (let ((sum 0f0))
    (declare (type (single-float 0f0) sum))
    (dotimes (k +dim+)
      (let* ((v (- (px i k) (px j k)))
	     (v2 (* v v)))
	(declare (type (single-float 0f0) v2))
	(incf sum v2)))
    sum))

(defun nearest-neighbour (target kd-tree)
   (declare (type array-index-t target)
	    (type kd-tree kd-tree)
	    (values array-index-t single-float &optional))
   (with-slots (perm points root)
       kd-tree
     (let* ((nndist2 1f20)
	    (nearest 0)
	    (*points* points)
	    (*perm* perm))
       (labels ((rec (node)
		  (declare (type (or leaf node) node))
		  (etypecase node
		    ;; sequential search in the bucket of the leaf
		    (leaf (with-slots (lopt hipt)
			      node
			    (loop for i from lopt upto hipt do
				 (let ((d (distance2 (aref perm i) target)))
				   (when (< d nndist2)
				     (setf nndist2 d
					   nearest (aref perm i)))))
			    (return-from rec nil)))
		    ;; search in the closer son, if nothing there,
		    ;; search in farther son as well
		    (node (with-slots (cutval cutdim loson hison)
			      node
			    (let* ((x (px target cutdim))
				   (diff (- x cutval)))
			      (if (< diff 0f0)
				  (progn (rec loson)
					 (when (< (expt diff 2) nndist2) ;(< cutval (+ x dist))
					   (rec hison)))
				  (progn (rec hison)
					 (when (< (expt diff 2) nndist2) ;(< (- x dist) cutval)
					   (rec loson))))))))))
	 (rec root))
       (values nearest nndist2))))

#+nil
(defun locate-points-in-circle (center radius tree)
  (declare (type (simple-array single-float 1) center)
	   (type single-float radius)
	   (type kd-tree tree))
  (with-slots ((*perm* perm) (*points* points) root) tree
    (labels ((dist (a-i b)
	       (declare 
		(type array-index-t a-i)
		(type (simple-array single-float 1) b))
	       (let ((sum 0f0))
		 (dotimes (i (length b))
		   (incf sum (* (px a-i i) (aref b i))))
		 sum))
	     (rec (node)
	       (declare (type (or leaf node) node))
	       (etypecase node
		 (leaf (with-slots (lopt hipt) node
			 (loop for i from lopt upto hipt do
			      (let ((d (dist i center)))))))))))))
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

(defun split-box (dim box lo-p cut) ;; its actually not box but rectangle
  "A box is a list 4 values (the coordinates of the min and max
point). This function splits it along the direction DIM (which is
either 0 or 1). When LO-P is true the part with the lower values is
returned."
  (declare ((member 0 1) dim)
	   (list box)
	   (boolean lo-p)
	   (single-float cut)
	   (values list &optional))
  (destructuring-bind (px py qx qy)
      box
    (if (eq dim 0)
	 (if lo-p
	     (list px py cut qy)
	     (list cut py qx qy))
	 (if lo-p
	     (list px py qx cut)
	     (list px cut qx qy)))))
#+nil
(split-box 0 (list 0 0 1 1) t .5)

(defun readable-float (x)
  "Prints only a few digits of a float X into a string, to make the
output more readable."
  (read-from-string 
   (format nil "~2,4f" x)))


(defun draw-tree-boxes (root)
  (declare (node root))
  (let ((boxes nil))
    ;; add boxes of only those nodes at the bottom of the tree, that
    ;; have a leaf in loson or hison
    (labels ((rec (node box)
	      (when (typep node 'node)
		(with-slots (loson hison cutdim cutval)
		    node
		  (when (or (typep loson 'leaf)
			    (typep hison 'leaf))
		    (push box boxes)
		    #+nil (format t "~a~%" (mapcar #'readable-float box)))
		  (rec loson (split-box cutdim box t cutval))
		  (rec hison (split-box cutdim box nil cutval))))))
      (rec root (list 0f0 0f0 1f0 1f0)))
    boxes))

#+nil
(defparameter *boxes*
  (draw-tree-boxes (kd-tree-root *tree*)))


#+nil
(with-open-file (s "/home/martin/tmp/tree.dot" :direction :output
		   :if-exists :supersede)
 (dot-draw-tree s (kd-tree-root *tree*)))

(defun eps-moveto (x y)
  (declare (single-float x y))
  (format nil "~f ~f moveto~%" x y))

(defun eps-lineto (x y)
  (declare (single-float x y))
  (format nil "~f ~f lineto~%" x y))

(defun eps-rectangle (box)
  (declare (list box))
  (destructuring-bind (x0 y0 x y)
      box
    (declare (single-float x0 y0 x y))
    (format nil "newpath~%~a~a~a~a~astroke~%"
	   (eps-moveto x0 y0)
	   (eps-lineto x y0)
	   (eps-lineto x y)
	   (eps-lineto x0 y)
	   (eps-lineto x0 y0))))

(defun eps-point (x y)
  (declare (single-float x y))
  (format nil "newpath ~f ~f 0.04 0 360 arc closepath fill~%" x y))


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
   (loop for b in boxes do
	(format s "~a" (eps-rectangle b)))
   (loop for p across points do
	(format s "~a" (eps-point (aref p 0) (aref p 1))))
   (format s "%%EOF"))))

#+nil
(eps-tree "/dev/shm/o.ps"
	  (draw-tree-boxes (kd-tree-root gauss::*tree*)) 
	  (kd-tree-points gauss::*tree*))
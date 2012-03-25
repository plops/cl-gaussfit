(declaim (optimize (speed 2) (safety 3) (debug 3)))

(defpackage :kdtree
  (:use :cl :math.vector))
(in-package :kdtree)

(defmacro with-arrays (arrays &body body)
  "Provides a corresponding accessor for each array as a local macro,
so that (ARRAY ...) corresponds to (AREF ARRAY ...)."
  `(macrolet ,(mapcar (lambda (array)
                        `(,array (&rest indices) `(aref ,',array ,@indices)))
                      arrays)
     ,@body))

(defconstant +dim+ 2)

(defstruct (kd-node (:include box)
		    (:conc-name "KD-")
		    (:predicate nil)
		    (:copier nil))
  (mom 0 :type fixnum)
  (left 0 :type fixnum)
  (right 0 :type fixnum)
  (point-min 0 :type fixnum)
  (point-max 0 :type fixnum))

#+nil
(let ((a (make-kd-node :min (v .1 .2))))
  (list (sb-sys:vector-sap (vec-coord (box-min a))) 
	(sb-sys:vector-sap (vec-coord (box-min (copy-box a))))))

(defun copy-kd-node (node)
  (let ((a))))


(defstruct (kd-tree (:predicate nil)
		    (:copier nil))
  (points (make-array 0 :element-type 'vec
		      :initial-element (v))
	  :type (simple-array vec 1))
  (boxes (make-array 0 :element-type 'kd-node
		     :initial-element (make-kd-node))
	 :type (simple-array kd-node 1))
  (point-index (make-array 0 :element-type 'fixnum)
	       :type (simple-array fixnum 1))
  (reverse-point-index (make-array 0 :element-type 'fixnum)
		       :type (simple-array fixnum 1))
  (coords (make-array 0 :element-type 'double-float)
	  :type (simple-array double-float 1)))


(defun copy-array (a)
  (declare ((array * *) a)
	   (values (array * *) &optional))
  (let* ((result (make-array (array-dimensions a)
			     :element-type (array-element-type a)))
	 (result1 (sb-ext:array-storage-vector result))
	 (a1 (sb-ext:array-storage-vector a)))
    (dotimes (i (length result1))
      (setf (aref result1 i) (aref a1 i)))
    result))
#+ni
(sb-impl::%fun-type #'select-index)
#+nil
(sb-kernel:specifier-type (sb-impl::%fun-type #'select-index))

(defun select-index (in ar k)
  "Given K in [0..n-1] modifies the index array IN to point to the K+1
  smallest values AR[IN[k]] in the front (in arbitrary order). Return
  value is IN[K] and therefore the K+1-smallest value."
  (declare ((array double-float 1) ar)
	   ((array fixnum 1) in)
	   (fixnum k)
	   (values fixnum &optional))
  (let* ((n (length in))
	 (l 0) ;; index to leftmost element
	 (ir (1- n))) ;; index to rightmost element
    (format t "select-index ~a~%" (list 'n n (when (< n 4) in)))
    (labels ((swap (p q)
	       (unless (and (array-in-bounds-p in p)
			    (array-in-bounds-p in q))
		 (error "out of bounds for index array."))
	       (let ((h (aref in p)))
		 (setf (aref in p) (aref in q)
		       (aref in q) h))))
      (loop
	 (with-arrays (in ar)
	   (if (<= ir (1+ l))
	       ;; active partition contains 1 or 2 el
	       (progn (when (and (eq ir (1+ l))
				 (< (ar (in ir)) (ar (in l))))
			(swap l ir))
		      (return-from select-index (in k)))
	       ;; choose median of left center and right as partitioning element
	       (let ((mid (floor (+ l ir) 2)))
		 (swap mid (1+ l))
		 ;; rearrange for ar[l]<=ar[l+1]<=ar[ir]
		 (when (< (ar (in ir)) (ar (in l)))      (swap ir l))
		 (when (< (ar (in ir)) (ar (in (1+ l)))) (swap ir (1+ l)))
		 (when (< (ar (in (1+ l))) (ar (in l)))  (swap l (1+ l)))
		 (let* ((i (1+ l))
			(j ir)
			(ia (in i))
			(a (ar ia)))
		   (loop ;; scan up and down, put big elements right,
			 ;; small ones left
		     (loop do (incf i) while (< (ar (in i)) a))
		     (loop do (decf j) while (< a (ar (in j))))
		     (when (< j i) (return))
		     (swap i j))
		  (setf (in (1+ l)) (in j)
			(aref in j) ia)
		  (when (<= k j) (setf ir (1- j)))
		  (when (<= j k) (setf l i))))))))))

#+nil
(let* ((ls '(.1 .9 .1 .7 .2 .4 .6 .3 .8 .1 .5))
       (n (length ls))
       (ar (make-array n :element-type 'double-float
		       :initial-contents (mapcar #'(lambda (x) (* 1d0 x)) ls)))
       (in (make-array n :element-type 'fixnum
		       :initial-contents (loop for i below n collect i))))
 (let ((ret (select-index in ar (floor n 2))))
   (list ret 
	 (read-from-string (format nil "~2,1f" (aref ar ret)))
	 (loop for i below n collect 
	      (read-from-string (format nil "~2,1f" (aref ar (aref in i))))))))

;; p. 1121 nr3.pdf
;; strategy for setting up kd-tree POINT-INDEX is an array of integers
;; that index N points, copy the point coordinates into COORDS with
;; all x-coordinates contiguous followed by y and then z. use
;; SELECT-INDEX to partition INDEX-POINTS by their x-coordinates with
;; half the number of points on each side of the partition. These
;; halfes are now to be partitioned again, this time by the value of
;; the y-coordinate. Every daughter gets a pointer to its mother and
;; two indices to the beginning and end of their points.

(deftype axis ()
  `(member ,@(loop for i below +dim+ collect i)))

(defun next-axis (axis)
  (declare (axis axis))
  (mod (1+ axis) +dim+))

(defun kd-tree (in-points)
  (declare ((simple-array vec 1) in-points)
	   (values kd-tree &optional))
  (let* ((n (length in-points))
	 (m (expt 2 (ceiling (log n 2))))
	 (point-index (make-array 
		       n :element-type 'fixnum
		       :initial-contents (loop for i below n collect i)))
	 (nb (1- (min m (- (* 2 n) (floor m 2)))))
	 (boxes (make-array nb :element-type 'kd-node
			    :initial-contents
			    (loop for i below nb collect (make-kd-node))))
	 (coords (make-array 
		  (* +dim+ n) :element-type 'double-float
		  :initial-contents 
		  (nconc 
		   (loop for i below n collect (vec-x (aref in-points i)))
		   (when (<= 2 +dim+)
		     (loop for i below n collect (vec-y (aref in-points i))))
		   (when (<= 3 +dim+) 
		     (loop for i below n collect (vec-z (aref in-points i)))))))
	 (big 12d0)
	 (ma (v big big big))
	 (mi (v* ma -1d0))
	 (task-mom (make-array 50 :element-type 'fixnum))
	 (task-axis (make-array 50 :element-type 'axis))
	 (nowtask 1)
	 (j-box 0))
    (with-arrays (task-mom task-axis boxes mi ma coords)
      (setf (boxes 0) (make-kd-node :min mi :max ma
				    :mom 0 :left 0 :right 0
				    :point-min 0
				    :point-max (1- n))
	    (task-mom 1) 0
	    (task-axis 1) 0)
      (loop while (< 0 nowtask) do
	  (let* ((mom (task-mom nowtask))
		 (axis (task-axis nowtask))
		 (ptmin (kd-point-min (boxes mom)))
		 (ptmax (kd-point-max (boxes mom)))
		 (number-of-points-in-subdivision (1+ (- ptmax ptmin)))
		 (point-index-sub (make-array number-of-points-in-subdivision
					      :element-type 'fixnum
					      :displaced-to point-index
					      :displaced-index-offset ptmin))
		 (coord-sub (make-array n
					:element-type 'double-float
					:displaced-to coords
					:displaced-index-offset (* axis n)))
		 (boundary-point-left 
		  (floor (1- number-of-points-in-subdivision) 2)))
	    (decf nowtask)
	    (select-index point-index-sub coord-sub boundary-point-left)
	    (setf ma (copy-array (kd-max (boxes mom)))
		  mi (copy-array (kd-min (boxes mom)))
		  (ma axis) (coords (+ (aref point-index-sub boundary-point-left)
					(* n axis)))
		  (mi axis) (ma axis)
		  (boxes (incf j-box)) 
		  (make-kd-node :min (copy-array (kd-min (boxes mom)))
				:max ma
				:mom mom :left 0 :right 0
				:point-min ptmin
				:point-max (+ ptmin boundary-point-left))
		  (boxes (incf j-box)) 
		  (make-kd-node :min mi
				:max (copy-array (kd-max (boxes mom)))
				:mom mom :left 0 :right 0
				:point-min (+ ptmin boundary-point-left 1)
				:point-max ptmax)
		  (kd-left (boxes mom)) (1- j-box)
		  (kd-right (boxes mom)) j-box)
	    (format t "~a~%" (list 'left (1- j-box) 'right j-box 'nb nb))
	    
	    (when (< 1 boundary-point-left)
	      (setf (task-mom (incf nowtask)) (1- j-box)
		    (task-axis nowtask) (next-axis axis)))
	    (when (< 2 (- number-of-points-in-subdivision boundary-point-left))
	      (setf (task-mom (incf nowtask)) j-box
		    (task-axis nowtask) (next-axis axis))))))
    (let ((rpoint-index (make-array n :element-type 'fixnum)))
      (loop for i below n do
	   (setf (aref rpoint-index (aref point-index i)) i))
      (make-kd-tree :points in-points
		    :point-index point-index
		    :reverse-point-index rpoint-index
		    :boxes boxes
		    :coords coords))))
#+nil
(progn 
  (defparameter dkaf
    (let* ((n 8) 
	   (points (make-array n :element-type 'vec
			       :initial-contents 
			      (loop for i below n collect
				   (v (- (random 20d0) 10d0)
				      (- (random 20d0) 10d0))))))
      (kd-tree points)))
#+nil  (print-kd-tree-to-eps "/home/martin/tmp/o.eps" dkaf :scale 10d0))
#+nil
(defparameter kldsa
 (with-slots (boxes points)
     dkaf
   (with-arrays (boxes)
    (loop for i below (length boxes) collect
	 (with-slots (left right)
	     (boxes i)
	   (if (and (eq left 0)
		    (eq right 0))
	       (boxes i)))))))

#|
%!PS-Adobe-3.0
%%Pages: 1
%%BoundingBox: 0 0 226.77167 140.15259
%%EndComments
2.8346457 2.8346457 scale
2 2 translate
0.35 setlinewidth 0 setgray newpath
200 204 moveto
100000 0 lineto
0 100000 lineto
200 204 lineto
stroke
%%EOF
|#

(defun eps-moveto (vec)
  (declare (vec vec)
	   (values string &optional))
  (format nil "~f ~f moveto~%" (vec-x vec) (vec-y vec)))

(defun eps-lineto (vec)
  (declare (vec vec)
	   (values string &optional))
  (format nil "~f ~f lineto~%" (vec-x vec) (vec-y vec)))

(defun eps-rectangle (min max)
  (declare (vec min max)
	   (values string &optional))
  (let* ((delta (v- max min))
	 (dx (v (vec-x delta)))
	 (dy (v 0d0 (vec-y delta)))
	 (offset (v 15d0 15d0)))
    (format nil "0 setgray~%newpath~%~a~a~a~a~astroke~%"
	    (eps-moveto (v+ offset min))
	    (eps-lineto (v+ offset (v+ min dx)))
	    (eps-lineto (v+ offset max))
	    (eps-lineto (v+ offset (v+ min dy)))
	    (eps-lineto (v+ offset min)))))

(defun print-kd-tree-to-eps (fn kd-tree &key (scale 1d0))
  (declare (string fn)
	   (kd-tree kd-tree)
	   ((double-float 0d0) scale)
	   (values null &optional))
  (with-open-file (s fn :direction :output :if-exists :supersede)
    (format s "%!PS-Adobe-3.0
%%Pages: 1
%%BoundingBox: 0 0 700 500
%%EndComments~%")
    (format s "~f ~f scale~%" scale scale)
    (with-slots (boxes points)
	kd-tree
      (with-arrays (boxes points)
       (let ((dz (v 1.d0 1.d0)))
	 (format t "~a" (list 'print 'length-boxes (length boxes)))
	 (dotimes (i (length boxes))
	   (with-slots (max min left right)
	       (aref boxes i)
	     (format t "~a~%" (list i left right))
	     (unless (and (eq left 0)
			  (eq right 0))
	       (format s (eps-linewidth (* .02d0 i)))
	       (format s (eps-rectangle (v+ (v* dz (* .04d0 i)) min)
					(v+ (v* dz (* .04d0 i)) max))))))
	 (format s (eps-linewidth .9d0))
	 (dotimes (i (length points))
	   (format s (eps-point (v+ (v 15d0 15d0) (points i))))))))
    (format s "%%EOF~%"))
  nil)
#+nil
(print-kd-tree-to-eps "/home/martin/tmp/o.eps" dkaf :scale 10d0)
(defun eps-linewidth (s)
  (declare (double-float s)
	   (values string &optional))
  (format nil "~f setlinewidth~%" s))
(defun eps-point (vec)
  (declare (vec vec)
	   (values string &optional))
  (format nil "~f ~f stroke [] 0 setdash gsave 1 setlinecap moveto 0 0 rlineto stroke grestore~%" (vec-x vec) (vec-y vec)))

(defun print-kd-tree (kd-tree)
  (declare (kd-tree kd-tree)
	   (values null &optional))
  (with-slots (boxes)
      kd-tree 
      (dotimes (i (length boxes))
	(format t "~a~%" (list i (aref boxes i))))))

#+nil
(print-kd-tree dkaf)

(defun kd-tree-point-distance (tree p q)
  (declare (kd-tree tree)
	   (fixnum p q)
	   (values double-float &optional))
  (with-slots (points)
      tree
    (if (eq p q)
	1d20
	(norm (v- (aref points p) (aref points q))))))

#+nil
(kd-tree-point-distance dkaf 1 4)

(defun locate-point-in-tree (point kd-tree)
  (declare (vec point)
	   (kd-tree kd-tree)
	   (values fixnum &optional))
  (with-slots (boxes)
      kd-tree
    (let ((d1 0)
	  (nb 0)
	  (jdim 0))
      (loop while (< 0 (kd-left (aref boxes nb))) do
	   (setf
	    d1 (kd-left (aref boxes nb))
	    nb (if (<= (aref point jdim) (aref (kd-max (aref boxes d1)) jdim))
		   d1
		   (kd-right (aref boxes nb)))
	    jdim (mod (1+ jdim) 3)))
      nb)))
#+nil
(locate-point-in-tree (v 1d0 1d0 1d0) dkaf)
#+nil
(aref (kd-tree-boxes dkaf)
      (locate-point-in-tree (v 1d0 1d0 1d0) dkaf))

(defun locate-point-index-in-tree (point-index tree)
  "Return the index of the box that contains the point given by
POINT-INDEX."
  (declare (fixnum point-index)
	   (kd-tree tree)
	   (values fixnum &optional))
  (with-slots (boxes reverse-point-index)
      tree
    (with-arrays (boxes)
     (let ((jh (aref reverse-point-index point-index))
	   (nb 0))
       (loop while (< 0 (kd-left (boxes nb))) do
	    (let ((left (kd-left (boxes nb))))
	      (setf nb (if (<= jh (kd-point-max (boxes left)))
			   left
			   (kd-right (boxes nb))))))
       nb))))

#+nil
(locate-point-index-in-tree 9 dkaf)
#.(require :alexandria)

;; Implementation of kd-trees, a datastructure that stores points of a
;; k-dimensional space and enables quick nearest neighbour search.  I
;; write an application that captures z-stacks of c. elegans embryos
;; (small worms) with a microscope. I hope I can follow the centers of
;; single cells as the embyro splits into more and more cells.  This
;; file contains code for debugging output showing the 2d rectangle
;; filled with random points and the corresponding boxes that make up
;; the tree. There is also code to output the final tree as input for
;; dot (apt-get install graphviz) to render a postscript image of the
;; tree. You can see screenshots of those on
;; http://imgur.com/dyzhC&eGPb9l

;; 2010-07-21 kielhorn.martin@googlemail.com


;; the paper I read for this is 1990bentley_kdtree-c++.pdf

(declaim (optimize (speed 2) (safety 3) (debug 3)))

(defpackage :kdtree
  (:use :cl :alexandria)
  (:export #:build-new-tree
	   #:nearest-neighbour
	   #:kd-tree-points
	   #:perm
	   #:get-tree-point))

(in-package :kdtree)

(defconstant +dim+ 2) ;; you can set this to three but then the postscript output doesn't work
(declaim (fixnum +dim+))

(deftype vec ()
  `(simple-array double-float (#.+dim+)))

(deftype axis ()
  `(member ,@(loop for i below +dim+ collect i)))

(deftype array-index ()
  `(unsigned-byte 29))

(defstruct leaf
  (lopt (required-argument :lopt) :type array-index)
  (hipt (required-argument :hipt) :type array-index))

(defstruct node
  (cutdim (required-argument :cutdim) :type axis)
  (cutval (required-argument :cutval) :type double-float)
  (loson (required-argument :loson) :type (or node leaf))
  (hison (required-argument :hison) :type (or node leaf)))

(defstruct kd-tree
  (root nil :type (or null node))
  (points (required-argument :points) :type (simple-array vec 1))
  (perm (required-argument :perm) :type (simple-array array-index 1)))

(defparameter *points* (make-array 0 :element-type 'vec))
(defparameter *perm* (make-array 0 :element-type 'array-index))
(declaim ((simple-array vec 1) *points*)
	 ((simple-array array-index 1) *perm*))
(declaim (inline px))
(defun px (i j)
  (declare (array-index i j)
	   (values double-float &optional))
  (aref (aref *points* (aref *perm* i)) j))

(defun get-tree-point (kd-tree i)
  (declare (kd-tree kd-tree)
	   (fixnum i)
	   (values vec &optional))
  (with-slots (points
	       perm)
      kd-tree
    (aref points (aref perm i))))

(defun max-spread-in-dim (dim l u)
  (declare (array-index l u)
	   (axis dim)
	   (values double-float &optional))
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
  (declare (array-index l u)
	   (values axis &optional))
  (let ((spread 0d0)
	(dim (the axis 0)))
    (declare (double-float spread))
    (dotimes (j +dim+)
      (let ((v (max-spread-in-dim j l u)))
	(when (< spread v)
	  (setf spread v
		dim j))))
    dim))

(defun make-random-vec-of-dim (dim)
  (declare ((integer 0) dim)
	   (values (simple-array double-float 1) &optional))
  (make-array dim
	      :element-type 'double-float
	      :initial-contents (loop repeat dim collect
				     (random 1d0))))

(defun make-random-points (n)
  (declare (array-index n)
	   (values (simple-array vec 1) &optional))
  (make-array n :element-type '(array * (#.+dim+))
	      :initial-contents
	      (loop for i below n collect
		   (make-random-vec-of-dim +dim+))))

#+nil ;; check that all dimensions appear for random points
(loop repeat 30 collect
     (let* ((n 1000)
	    (kd (build-new-tree (make-random-points n)))
	    (*perm* (kd-tree-perm kd))
	    (*points* (kd-tree-points kd)))
       (find-max-spread-dimension 0 (1- n))))

(defun swap (a b)
  (declare (array-index a b)
	   (values null &optional))
  (let ((h (aref *perm* a)))
    (setf (aref *perm* a) (aref *perm* b)
	  (aref *perm* b) h))
  nil)

(defun select (l u m cutdim)
  (declare (array-index l u m)
	   (axis cutdim)
	   (values fixnum &optional))
  (labels ((p (l) ;; accessor
	     (declare (fixnum l)
		      (values double-float &optional))
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
	     (labels ((enforce-order (l u)
			(declare (fixnum l u)
				 (values null &optional))
			(when (< (p u) (p l)) (swap u l))
			nil))
	       ;; rearrange for p[l]<=p[l+1]<=p[u]
	       (enforce-order l u)
	       (enforce-order (1+ l) u)
	       (enforce-order l (1+ l)))
	     (let* ((i (1+ l))
		    (j u)
		    (ia (aref *perm* i))
		    (a (p i)))
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

#+nil 
(let* ((n 35)
       (*perm* (make-perm n))
       (*points* (make-random-points n)))
  (select 0 (1- n) (floor (+ 0 (1- n)) 2)
	  0)
  (loop for i below n collect (read-from-string
			       (format nil "~d" (px i 0)))))

(defun build (l u)
  (declare (array-index l u)
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
		    :loson (build l m)
		    :hison (build (1+ m) u))))))

(defun build-new-tree (points)
  (declare ((simple-array vec 1) points)
	   (values kd-tree &optional))
  (let* ((n (length points))
	 (*points* points)
	 (*perm* (make-array n
			     :element-type 'array-index
			     :initial-contents 
			     (loop for i below n collect i)))
	 (root (build 0 (1- n))))
    (make-kd-tree :points *points*
		  :perm *perm*
		  :root root)))
#+nil
(let* ((n 300))
  (defparameter *tree* (build-new-tree (make-random-points n))))

(defun distance (i j)
  (declare (array-index i j)
	   (values double-float &optional))
  (let ((sum 0d0))
    (dotimes (k +dim+)
      (let ((v (- (px i k) (px j k))))
	(incf sum (* v v))))
    (sqrt sum)))

(defun nearest-neighbour (target kd-tree)
   (declare (array-index target)
	    (kd-tree kd-tree)
	    (values array-index &optional))
   (with-slots (perm points root)
       kd-tree
     (let* ((dist 1d20)
	    (nearest 0)
	    (*points* points)
	    (*perm* perm))
       (labels ((rec (node)
		  (declare ((or leaf node) node))
		  (etypecase node
		    ;; sequential search in the bucket of the leaf
		    (leaf (with-slots (lopt hipt)
			      node
			    (loop for i from lopt upto hipt do
				 (let ((d (distance (aref perm i) target)))
				   (when (< d dist)
				     (setf dist d
					   nearest (aref perm i)))))
			    (return-from rec nil)))
		    ;; search in the closer son, if nothing there,
		    ;; search in farther son as well
		    (node (with-slots (cutval cutdim loson hison)
			      node
			    (let ((x (px target cutdim)))
			      (if (< x cutval)
				  (progn (rec loson)
					 (when (< cutval (+ x dist))
					   (rec hison)))
				  (progn (rec hison)
					 (when (< (- x dist) cutval)
					   (rec loson))))))))))
	 (rec root))
       nearest)))

#+nil
(let* ((n 30000))
  (time (defparameter *tree* (build-new-tree (make-random-points n))))
  (time (dotimes (i n) (nearest-neighbour i *tree*))))

;; 10s to find all nn for 30000 points when 2, 9 or 14 points in bucket
;; 11.5s for 20 points in bucket
;; 14s for 50 points in bucket


;; the following functions generate debugging output but only work for +dim+=2
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
	   (double-float cut)
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
      (rec root (list 0d0 0d0 1d0 1d0)))
    boxes))

#+nil
(defparameter *boxes*
  (draw-tree-boxes (kd-tree-root *tree*)))

(defun dot-draw-tree (stream root)
  "Convert the tree into input for the program dot. The output can be
converted with dot -Tps -o test.ps file.dot."
  (declare (node root))
  (labels ((transform-dim (cutdim)
	     (ecase cutdim
	       (0 'x)
	       (1 'y)))
	   (rec (node)
	     (declare ((or node leaf) node))
	     (when (typep node 'node)
	       (with-slots (hison loson cutdim cutval)
		   node
		 (labels ((print-daughter (node)
			    (declare ((or node leaf) node))
			    (etypecase node
			      (leaf (list (leaf-lopt node)
					  (leaf-hipt node)))
			      (node (list (transform-dim (node-cutdim node))
					(readable-float (node-cutval node)))))))
		   (let ((here (list (transform-dim cutdim)
				     (readable-float cutval)))
			 (left (print-daughter loson))
			 (right (print-daughter hison)))
		     (format stream "\"~a\" -> \"~a\";~%"
			     here left)
		     (format stream "\"~a\" -> \"~a\";~%"
			     here right)))
		 (rec loson)
		 (rec hison)))))
    (format stream "digraph {~%")
    (rec root)
    (format stream "}~%")))

#+nil
(with-open-file (s "/home/martin/tmp/tree.dot" :direction :output
		   :if-exists :supersede)
 (dot-draw-tree s (kd-tree-root *tree*)))

(defun eps-moveto (x y)
  (declare (double-float x y))
  (format nil "~f ~f moveto~%" x y))

(defun eps-lineto (x y)
  (declare (double-float x y))
  (format nil "~f ~f lineto~%" x y))

(defun eps-rectangle (box)
  (declare (list box))
  (destructuring-bind (x0 y0 x y)
      box
    (declare (double-float x0 y0 x y))
    (format nil "newpath~%~a~a~a~a~astroke~%"
	   (eps-moveto x0 y0)
	   (eps-lineto x y0)
	   (eps-lineto x y)
	   (eps-lineto x0 y)
	   (eps-lineto x0 y0))))

(defun eps-point (x y)
  (declare (double-float x y))
  (format nil "newpath ~f ~f 0.004 0 360 arc closepath fill~%" x y))


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
400.0 400.0 scale
0.002 setlinewidth
0 setgray~%")
   (loop for b in boxes do
	(format s "~a" (eps-rectangle b)))
   (loop for p across points do
	(format s "~a" (eps-point (aref p 0) (aref p 1))))
   (format s "%%EOF")))

#+nil
(progn
  (let* ((n 300))
    (defparameter *tree* (build-new-tree (make-random-points n))))
  (defparameter *boxes* (draw-tree-boxes (kd-tree-root *tree*)))
  (eps-tree "/home/martin/tmp/tree.eps" *boxes*
	    (kd-tree-points *tree*)))
)
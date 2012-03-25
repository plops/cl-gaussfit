;; the paper I read for this is 1990bentley_kdtree-c++.pdf
#.(require :alexandria)
(declaim (optimize (speed 0) (safety 3) (debug 3)))

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
  (loson (required-argument :loson) :type (or node leaf))
  (hison (required-argument :hison) :type (or node leaf)))

(defstruct kd-tree
  (root nil :type (or null node))
  (points (required-argument :points) :type (simple-array vec 1))
  (perm (required-argument :perm) :type (simple-array array-index-t 1)))

(defparameter *points* (make-array 0 :element-type 'vec))
(defparameter *perm* (make-array 0 :element-type 'array-index-t))
(declaim ((simple-array vec 1) *points*)
	 ((simple-array array-index-t 1) *perm*))
(declaim (inline px))
(defun px (i j)
  (declare (array-index-t i j)
	   (values single-float &optional))
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
  (declare (array-index-t l u)
	   (axis dim)
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
  (declare (array-index-t l u)
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

(defun make-random-vec-of-dim (dim)
  (declare ((integer 0) dim)
	   (values (simple-array single-float 1) &optional))
  (make-array dim
	      :element-type 'single-float
	      :initial-contents (loop repeat dim collect
				     (random 1f0))))

(defun make-random-points (n)
  (declare (array-index-t n)
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
  (declare (array-index-t a b)
	   (values null &optional))
  (let ((h (aref *perm* a)))
    (setf (aref *perm* a) (aref *perm* b)
	  (aref *perm* b) h))
  nil)

(defun select (l u m cutdim)
  (declare (array-index-t l u m)
	   (axis cutdim)
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
  (declare (array-index-t l u)
	   (values (or node leaf) &optional))
  (let ((points-in-bucket 14)) ;; change this back to somewhere like 14 for better performance
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

(defun distance (i j)
  (declare (array-index-t i j)
	   (values single-float &optional))
  (let ((sum 0f0))
    (dotimes (k +dim+)
      (let ((v (- (px i k) (px j k))))
	(incf sum (* v v))))
    (sqrt sum)))

(defun nearest-neighbour (target kd-tree)
   (declare (array-index-t target)
	    (kd-tree kd-tree)
	    (values array-index-t &optional))
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

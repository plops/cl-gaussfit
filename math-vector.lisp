(defpackage #:math.vector
  (:use :cl)
  (:export #:+dim+
	   #:v
	   #:copy-vec
	   #:vec-x
	   #:vec-y
	   #:vec-z
	   #:v.
	   #:norm
	   #:v+
	   #:v-
	   #:s*
	   #:make-box
	   #:copy-box
	   #:box-min
	   #:box-max
	   #:distance-point-box
	   ))

(in-package #:math.vector)

(declaim (optimize (speed 3) (safety 1) (debug 1)))

(defconstant +dim+ 2)

(defstruct (vec (:constructor %make-vec)
		(:copier %copy-vec)) 
  (coord (make-array +dim+ 
		     :element-type 'single-float
		     :initial-element 0f0)
	 :type (simple-array single-float #.(list +dim+))))


(defmacro def-coord-access ()
  "define vec-x, vec-y, ..., (setf vec-x), ..., (v .2 .4) and 
copy-v."
  (labels ((name (axis)
	     (intern (format nil "VEC-~a" axis))))
   (append '(progn)
	   (loop for a in '(x y z w) and i below +dim+ collect
		`(progn
		   (defun ,(name a) (v)
		     (declare (type vec v)
			      (values single-float &optional))
*		     (aref (vec-coord v) ,i))
		   (defun (setf ,(name a)) (new v)
		     (declare (type vec v)
			      (type single-float new)
			      (values single-float &optional))
		     (setf (aref (vec-coord v) ,i) new))))
	   `((defun v ,(append '(&optional)
			        (loop for a in '(x y z w) and i below +dim+ collect
				    `(,a 0f0)))
	       (declare ,(append '(type single-float)
				 (loop for a in '(x y z w) and i below +dim+
				    collect a)) 
			(values vec &optional))
	       ,(append '(let* ((ve (%make-vec))))
			(loop for a in '(x y z w) and i below +dim+ collect
			     `(setf (,(name a) ve) ,a))
			'(ve)))
	     (defun copy-vec (v)
	       (declare (type vec v)
			(values vec &optional))
	       (let ((b (%copy-vec v)))
		 (setf (vec-coord b) (copy-seq (vec-coord v)))
		 b))))))

(def-coord-access)

#+nil
(v 1 2)
#+nil
(v .1 .2)

#+nil
(let ((a (v .1 .2)))
  (list (copy-v a)
	(setf (vec-x a) .3)
	(vec-x a)
	a))

(defun v. (a b)
  (declare (type vec a b)
	   (values (single-float -1f0 1f0) &optional))
  (let ((sum 0f0)
	(ac (vec-coord a))
	(bc (vec-coord b)))
   (dotimes (i (length ac))
     (incf sum (* (aref ac i) (aref bc i))))
   sum))

(defun norm (v)
  (declare (type vec v)
	   (values (single-float 0f0 *) &optional))
  (let ((a (v. v v))) 
    (if (<= a 0f0)
	0f0
	(sqrt a))))

(defmacro def-vec-op (op)
  `(defun ,(intern (format nil "V~a" op)) (a b)
     (declare (type vec a b)
	      (values vec &optional))
     (let* ((c (v))
	    (cc (vec-coord c))
	    (ac (vec-coord a))
	    (bc (vec-coord b)))
       (dotimes (i (length ac))
	 (setf (aref cc i) (,op (aref bc i) (aref ac i))))
       c)))

(def-vec-op +)
(def-vec-op -)



#+nil
(v+ (v 1 2) (v 3 4))
#+nil
(reduce #'v+ (list (v 1 2)
		   (v 0 2)
		   (v 1 0)
		   (v 0 1)))

(defun s* (scalar vec) 
  (declare (type single-float scalar)
	   (type vec vec))
  (let* ((r (v))
	 (rc (vec-coord r))
	 (vc (vec-coord vec)))
    (dotimes (i (length vc))
      (setf (aref rc i) (* scalar (aref vc i))))
    r))

(defstruct (box (:copier %copy-box))
  (min (v) :type vec)
  (max (v) :type vec))

(defun copy-box (b)
  (declare (type box b)
	   (values box &optional))
  (let ((r (%copy-box b)))
    (setf (box-min r) (copy-vec (box-min b))
	  (box-max r) (copy-vec (box-max b)))))



;; when the point has a coordinate that is greater than hi or smaller
;; than lo it contributes to sum, otherwise not.
(defun distance-point-box (vec box)
  "If VEC is a point outside the BOX, its distance to the nearest
point on BOX is returned. If VEC is inside or on the surface 0.0 is
returned."
  (declare (type box box)
	   (type vec vec)
	   (values (single-float 0f0) &optional))
  (let ((sum 0f0))
    (declare (type (single-float 0f0) sum))
    (dotimes (i +dim+)
      (let ((p (aref (vec-coord vec) i))
	    (l (aref (vec-coord (box-min box)) i))
	    (h (aref (vec-coord (box-max box)) i)))
	(when (< p l)
	  (incf sum (expt (- l p) 2)))
	(when (< h p)
	  (incf sum (expt (- p h) 2)))))
    (sqrt sum)))

#+nil
(distance-point-box (v .1f0)
		    (make-box :min (v) :max (v 1f0 1f0)))

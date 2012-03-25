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
	   #:s*))

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
		     (aref (vec-coord v) ,i))
		   (defun (setf ,(name a)) (new v)
		     (declare (type vec v)
			      (type single-float new)
			      (values single-float &optional))
		     (setf (aref (vec-coord v) ,i) new))))
	   `((defun %v ,(append '(&optional)
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
	     (define-compiler-macro v (&optional (x 0f0) (y 0f0))
	       (if ,(append '(and)
			    (loop for a in '(x y z w) and i below +dim+ collect
				 `(typep ,a 'single-float)))
		   ,(append '(%v)
			   (loop for a in '(x y z w) and i below +dim+ collect
			      a))
		   ,(append '(%v)
			   (loop for a in '(x y z w) and i below +dim+ collect
			      `(* 1f0 ,a)))))
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
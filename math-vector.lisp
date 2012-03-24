(defpackage #:math.vector
  (:use :cl)
  (:export #:v
	   #:vec-x
	   #:vec-y
	   #:vec-z
	   #:v.
	   #:norm
	   #:v+
	   #:v-
	   #:s*
	   #:cross
	   #:rotation-matrix))

(in-package #:math.vector)

(declaim (optimize (speed 3) (safety 1) (debug 1)))

(deftype vec3 ()
    '(simple-array single-float (3)))


(defun v (&optional (x 0s0) (y 0s0) (z 0s0))
  (if (and (typep x 'single-float)
	   (typep y 'single-float)
	   (typep z 'single-float)) 
      (make-array 3 :element-type 'single-float ;; .04e-6s 
		  :initial-contents (list x y z))
      (make-array 3 :element-type 'single-float ;; .5e-6s
		  :initial-contents (mapcar #'(lambda (x) (* 1s0 x))
					    (list x y z)))))
#+nil
(time
 (dotimes (i 10000000)
   (let ((a (v .2))))))

(declaim (inline vec-x vec-y vec-z))
(declaim (ftype (function (vec3) single-float) vec-x vec-y vec-z))
(defun vec-x (v)
  (aref v 0))
(defun vec-y (v)
  (aref v 1))
(defun vec-z (v)
  (aref v 2))

(defun v. (a b)
  (declare (type vec3 a b))
  (let ((sum 0s0))
   (dotimes (i (length a))
     (incf sum (* (aref a i) (aref b i))))
   sum)) 

(defun norm (v)
  (declare (type vec3 v))
  (let ((a (v. v v))) 
    (if (< a 0s0)
	0s0
	(sqrt a))))

(defmacro def-vec-op (op)
  `(defun ,(intern (format nil "V~a" op)) (a b)
     (declare (type vec3 a b))
     (let ((c (v)))
       (dotimes (i (length a))
	 (setf (aref c i) (,op (aref b i) (aref a i))))
       c)))

(def-vec-op +)
(def-vec-op -)

#+nil
(v+ (v 1 2 3) (v 3 4 5))
#+nil
(reduce #'v+ (list (v 1 2 3)
		   (v 0 0 1)
		   (v 1 0 0)
		   (v 0 1 0)))

(defun s* (scalar vec) ;; .04e-6 s
  (declare (type single-float scalar)
	   (type vec3 vec))
  (let ((r (v)))
    (dotimes (i (length vec))
      (setf (aref r i) (* scalar (aref vec i))))
    r))

(defun cross (a b)
  (v (- (* (vec-y a) (vec-z b))
	(* (vec-z a) (vec-y b)))
     (- (* (vec-z a) (vec-x b))
	(* (vec-x a) (vec-z b)))
     (- (* (vec-x a) (vec-y b))
	(* (vec-y a) (vec-x b)))))

#+nil
(cross (v 0 0 1) (v 1 0 0))
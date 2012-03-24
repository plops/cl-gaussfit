(defpackage #:kielhorn.martin.math.vector
  (:use :cl)
  (:export #:v
	   #:vec-x
	   #:vec-y
	   #:vec-z
	   #:v.
	   #:v+
	   #:v-
	   #:s*
	   #:cross
	   #:rotation-matrix))

(in-package #:kielhorn.martin.math.vector)

(declaim (optimize (speed 3) (safety 1) (debug 1)))

(deftype vec3 ()
    '(simple-array single-float (3)))


(defun v (&optional (x 0 xp) (y 0 yp) (z 0 zp))
  (if (or xp yp zp)
      (make-array 3 :element-type 'single-float 
		  :initial-contents (mapcar #'(lambda (x) (* 1s0 x))
					    (list x y z)))
      (make-array 3 :element-type 'single-float 
		  :initial-element 0s0)))

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

(defun s*2 (scalar vec) ;; .3e-6 s
  (declare (type single-float scalar)
	   (type vec3 vec))
  (map '(simple-array single-float 1)
       #'(lambda (x) (* scalar x)) vec))

(defun s*1 (scalar vec) ;; .04e-6 s
  (declare (type single-float scalar)
	   (type vec3 vec))
  (let ((r (v)))
    (dotimes (i (length vec))
      (setf (aref r i) (* scalar (aref vec i))))
    r))
#+nil
(let ((a (v 1 2 3)))
  (time
   (dotimes (i 1000000)
     (s*2 .3 a))))

#+nil
(type-of (s*1 3s2 (v 0 2)))

(defun cross (a b)
  (v (- (* (vec-y a) (vec-z b))
	(* (vec-z a) (vec-y b)))
     (- (* (vec-z a) (vec-x b))
	(* (vec-x a) (vec-z b)))
     (- (* (vec-x a) (vec-y b))
	(* (vec-y a) (vec-x b)))))


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
    '(simple-array double-float (3)))

(defun v (x y &optional (z 0))
  (make-array 3 :element-type 'double-float 
	      :initial-contents (mapcar #'(lambda (x) (* 1d0 x))
					(list x y z))))

(defun vec-x (v)
  (aref v 0))
(defun vec-y (v)
  (aref v 1))
(defun vec-z (v)
  (aref v 2))

(defun v.3 (a b) ;; .57 s
  (declare (type vec3 a b))
  (reduce #'+ (map '(simple-array double-float 1) #'* a b)))

(defun v.2 (a b) ;; .2 s
  (declare (type vec3 a b))
  (loop for p across a and q across b sum (* p q)))

(defun v. (a b) ;; .025 s
  (declare (type vec3 a b)
	   (values double-float &optional))
  (let ((sum 0d0))
   (dotimes (i (length a))
     (incf sum (* (aref a i) (aref b i))))
   sum)) 

#+nil
(time
 (let ((a (v 1 2 3))
       (b (v 4 5 2)))
   (dotimes (i 1000000)
     (v. a b))))

(defmacro def-vec-op (op)
  `(defun ,(intern (format nil "V~a" op)) (a b) 
     (let ((c (make-array (array-dimensions a)
			  :element-type (array-element-type a))))
       (loop for p across a and q across b and i from 0 do
	    (setf (aref c i) (,op p q)))
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

(defun s* (scalar vec)
  (map 'vec3
       #'(lambda (x) (* scalar x)) vec))

#+nil
(type-of (s* 3 (v 0 2)))

(defun cross (a b)
  (v (- (* (vec-y a) (vec-z b))
	(* (vec-z a) (vec-y b)))
     (- (* (vec-z a) (vec-x b))
	(* (vec-x a) (vec-z b)))
     (- (* (vec-x a) (vec-y b))
	(* (vec-y a) (vec-x b)))))

(defun rotation-matrix (angle vect)
  "Create matrix that rotates by ANGLE radians around the direction
 VECT. VECT must be normalized."
  (let* ((u (vec-x vect))
	 (v (vec-y vect))
	 (w (vec-z vect))
	 (c (cos angle))
	 (s (sin angle))
	 (1-c (- 1 c))
	 (su (* s u))
	 (sv (* s v))
	 (sw (* s w)))
    (m (+ c (* 1-c u u))
       (+ (* 1-c u v) sw)
       (- (* 1-c u w) sv)

       (- (* 1-c u v) sw)
       (+ c (* 1-c v v))
       (+ (* 1-c v w) su)

       (+ (* 1-c u w) sv)
       (- (* 1-c v w) su)
       (+ c (* 1-c w w)))))

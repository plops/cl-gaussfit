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

(declaim (optimize (speed 0) (safety 3) (debug 3)))

(deftype vec3 ()
    '(simple-array double-float 3))

(defun v (x y &optional (z 0))
  (make-array 3 :element-type 'double-float :initial-contents (list x y z)))

(defun vec-x (v)
  (aref v 0))
(defun vec-y (v)
  (aref v 1))
(defun vec-z (v)
  (aref v 2))

(defun v. (a b)
  (reduce #'+ (map 'vector #'* a b))) 

(defun iformat (control-string &rest format-arguments)
  (intern 
   (string-upcase 
    (funcall #'format (append (list nil control-string)) 
				  format-arguments))))

#+nil
(iformat "v~a" '+)

(defmacro def-vec-op (op)
  (let ((mname (iformat "binary~s" op)))
   `(progn
      (defgeneric ,mname (a b))
      (defun ,(iformat "v~a" op) (&rest args) 
	(unless (null args)
	  (reduce #',mname args)))
      (defmethod ,mname ((x double-float) (y double-float))
	(+ x y))
      (defmethod ,mname ((x sequence) (y sequence))
	(map (type-of x) #',mname x y)))))

(def-vec-op +)

(defgeneric binary-subtract (a b))

(defun v- (&rest args) 
  "Subtract vectors element-wise."
  (unless (null args)
    (reduce #'binary-subtract args)))

(defmethod binary-subtract ((x number) (y number))
  (- x y))

(defmethod binary-subtract ((x sequence) (y sequence))
  (map (type-of x) #'binary-subtract x y))

(defun v* (scalar vec)
  (map 'vector #'(lambda (x) (* scalar x)) vec))


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

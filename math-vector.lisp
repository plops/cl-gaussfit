(defpackage #:kielhorn.martin.math.vector
  (:use :cl)
  (:export #:v
	   #:vec-x
	   #:vec-y
	   #:vec-z
	   #:vec-i
	   #:vec-conjugate
	   #:vec-realpart
	   #:dot
	   #:norm
	   #:v+
	   #:v-
	   #:v*
	   #:normalized
	   #:ray
	   #:ray-start
	   #:ray-direction
	   #:make-ray
	   #:dc
	   #:rad->deg
	   #:deg->rad
	   #:point-on-ray
	   #:cross
	   #:m
	   #:m*
	   #:mm*
	   #:transpose
	   #:mm+
	   #:ms*
	   #:determinant
	   #:rotation-matrix
	   #:square))

(in-package #:kielhorn.martin.math.vector)

(declaim (optimize (speed 0) (safety 3) (debug 3)))

(deftype vec3 ()
    '(simple-array number 3))

(defun v (x y &optional (z 0))
  (make-array 3 :element-type 'number :initial-contents (list x y z)))

(defun vec-x (v)
  (aref v 0))
(defun vec-y (v)
  (aref v 1))
(defun vec-z (v)
  (aref v 2))
(defun vec-i (v i)
  (aref v i))

(defun dot (a b)
  (reduce #'+ (map 'vector #'* a b))) 

(defun vec-conjugate (a)
  (map 'vector #'conjugate a))
(defun vec-realpart (a)
  (map 'vector #'realpart a))

(defun norm (a)
  (sqrt (dot a (vec-conjugate a))))

(defgeneric binary-add (a b))

(defun v+ (&rest args) 
  "Element-wise add vectors. E.g. (v+ (v 0 1) (v 0 2))"
  (unless (null args)
      (reduce #'binary-add args)))

(defmethod binary-add ((x number) (y number))
  (+ x y))

(defmethod binary-add ((x sequence) (y sequence))
  (map (type-of x) #'binary-add x y))

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

(defun normalized (v)
  (let ((len (norm v)))
    (assert (< 0 len) (len) "Trying to normalize zero vector.")
    (unless (= 0 len) 
      (v* (/ len) v))))

(defun dc (&key (x 0 x-p) (y 0 y-p) (z 1 z-p))
  "Construct direction cosinuses."
  (let ((res 
	 (flet ((circle (x y)
		  (let ((r2 (+ (* x x) (* y y))))
		    (assert (<= r2 1) 
			    (r2) 
			    "Coordinates ~a are not inside of a circle." 
			    (list x y))
		    (sqrt (- 1 r2))))) 
	   (cond ((and x-p y-p z-p)
		  (let ((vec (v x y z)))
		    (assert (= 1 (dot vec vec))
			    (vec)
			    "The given direction cosinuses ~a are not allowed."
			    vec)
		    (v x y z)))
		 ((and x-p z-p)
		  (v x (circle x z) z))
		 ((and x-p y-p)
		 (v x y (circle x y)))
		((and y-p z-p)
		 (v (circle y z) y z))
		((or x-p y-p z-p)
		 (v (if x-p 1 0)
		    (if y-p 1 0)
		    (if z-p 1 0)))
		((not (or x-p y-p z-p))
		 (v 0 0 1))))))
    (assert (and (<= -1 (vec-x res) 1) 
		 (<= -1 (vec-y res) 1)
		 (<= -1 (vec-z res) 1))
	    (res)
	    "The values ~a are out of the range [-1:1]." res)
    res)) 


(defstruct ray
  (start (v 0 0 0))
  (direction (dc :x 0 :y 0 :z 1)))

(defun rad->deg (x)
  (* 180 (/ pi) x))

(defun deg->rad (x)
  (* pi (/ 180) x))

(defun point-on-ray (ray nu)
  (v+ (ray-start ray) (v* nu (ray-direction ray))))


(defun cross (a b)
  (v (- (* (vec-y a) (vec-z b))
	(* (vec-z a) (vec-y b)))
     (- (* (vec-z a) (vec-x b))
	(* (vec-x a) (vec-z b)))
     (- (* (vec-x a) (vec-y b))
	(* (vec-y a) (vec-x b)))))

(defun m (a b c d e f g h i)
  (make-array '(3 3)
              :element-type 'number
              :initial-contents (list (list a b c) (list d e f) (list g h i))
              :adjustable nil))

(defun m* (matrix vect)
  (let ((res (v 0 0 0)))
    (loop for i below 3 do
	  (loop for j below 3 do
		(setf (aref res i) 
		      (+ (aref res i) (* (aref matrix i j) (aref vect j))))))
    res))

(defun mm* (matrix-a matrix-b)
  (let ((res (m 0 0 0
		0 0 0
		0 0 0)))
    (loop for i below 3 do
	  (loop for j below 3 do
		(loop for k below 3 do
		      (setf (aref res i j)
			    (+ (aref res i j)
			       (* (aref matrix-a k j)
				  (aref matrix-b i k)))))))
    res))

(defun ms* (s matrix)		     
  (let ((res (m 0 0 0 0 0 0 0 0 0)))
   (loop for i below 3 do
	 (loop for j below 3 do
	       (setf (aref res i j)
		     (* s (aref matrix i j)))))
   res))

(defun mm+ (matrix-a matrix-b)
  (let ((res (m 0 0 0 0 0 0 0 0 0)))
    (loop for i below 3 do
	  (loop for j below 3 do
		(setf (aref res i j)
		      (+ (aref matrix-a i j) (aref matrix-b i j)))))
    res))

(defun transpose (a)
  (let ((res (m 0 0 0 0 0 0 0 0 0)))
    (loop for i below 3 do
	  (loop for j below 3 do
		(setf (aref res i j)
		      (aref a j i))))
    res))

(defun determinant (m)
  (let ((a (aref m 0 0))
	(b (aref m 0 1))
	(c (aref m 0 2))
	(d (aref m 1 0))
	(e (aref m 1 1))
	(f (aref m 1 2))
	(g (aref m 2 0))
	(h (aref m 2 1))
	(i (aref m 2 2)))
    (+ (* a (- (* e i) (* f h)))
       (* b (- (* f g) (* d i)))
       (* c (- (* d h) (* e g))))))

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

(defun square (x)
  (* x x))
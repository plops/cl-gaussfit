(defpackage :utils
  (:use #:cl)
  (:export #:make-displaced-array
	   #:copy-img
	   #:extract-frame
	   #:img-op
	   #:img-apply-mask
	   #:find-local-maxima
	   #:draw-points-into-mask
	   #:extract-square
	   #:insert-rectangle
	   #:img-mul))

(in-package :utils)

(defun make-displaced-array (a)
  (make-array (array-total-size a)
	      :element-type (array-element-type a)
	      :displaced-to a))

(defun copy-img (img)
  (let* ((i1 (make-displaced-array img))
	 (d1 (make-array (array-dimensions i1)
			 :element-type (array-element-type i1)
			 :displaced-to i1))
	 (c1 (subseq d1 0))
	 (c (make-array (array-dimensions img)
			:element-type (array-element-type img)
			:displaced-to c1)))
    c))


#+nil
(time
 (write-fits "/dev/shm/hist.fits" (scale-log (copy-img *hists*))))

(defun extract-frame (img k)
  (destructuring-bind (z y x) (array-dimensions img)
    (declare (ignore z))
    (let* ((start (* x y k))
	   (i1 (make-array (* x y) :element-type (array-element-type img)
			   :displaced-to img
			   :displaced-index-offset start))
	   (c1 (subseq i1 0)))
      (make-array (list y x) 
		  :element-type (array-element-type img)
		  :displaced-to c1))))



(defun img-op (op a b)
  (let* ((a1 (make-displaced-array a))
	 (b1 (make-displaced-array b))
	 (c (make-array (array-dimensions a)
			:element-type 'single-float))
	 (c1 (make-displaced-array c)))
    (dotimes (i (length a1))
      (setf (aref c1 i) (* 1s0 (funcall op (aref a1 i) (aref b1 i)))))
    c))


(defun img-mul (im &optional (factor .001d0))
  (let* ((a (make-array (array-dimensions im)
			:element-type (array-element-type im)))
	 (i1 (make-displaced-array im))
	 (a1 (make-displaced-array a)))
    (dotimes (i (length a1))
      (setf (aref a1 i) (* factor (aref i1 i))))
    a))

(defun draw-points-into-mask (img points &key (mask-border 0))
  (destructuring-bind (h w) (array-dimensions img)
   (let ((a (make-array (array-dimensions img)
			:element-type 'boolean
			:initial-element nil)))
     (loop for j from mask-border below (- h mask-border) do
	  (loop for i from mask-border below (- w mask-border) do
	       (setf (aref a j i) t)))
     (dolist (p points)
       (destructuring-bind (y x val) p
	 (declare (ignore val))
	 (loop for j from -5 upto 5 do
	   (loop for i from -5 upto 5 do
		(let ((xx (+ x i))
		      (yy (+ y j)))
		  (when (and (<= 0 xx (1- w))
			     (<= 0 yy (1- h)))
		    (setf (aref a yy xx) nil)))))))
     a)))

(defun img-apply-mask (img mask &key (background 0))
  (unless (= (array-total-size img)
	     (array-total-size mask))
    (error "mask ~d and image ~d don't have the same length"
	   (array-total-size mask)
	   (array-total-size img)))
  (let* ((a1 (make-displaced-array img))
	 (b1 (make-displaced-array mask))
	 (c (make-array (array-dimensions img)
			:initial-element (coerce background
						 (array-element-type img))
			:element-type (array-element-type img)))
	 (c1 (make-displaced-array c)))
    (dotimes (i (length a1))
      (when (aref b1 i)
	(setf (aref c1 i) (aref a1 i))))
    c))

(defun find-local-maxima (im)
  "returns j i val"
  (destructuring-bind (h w) (array-dimensions im)
    (macrolet 
	((compare-expand (&rest pos)
	   `(and ,@(remove-if 
		    #'null
		    (loop for e in pos and p in '((-1 1) (0 1) (1 1)
						  (-1 0) (0 0) (1 0)
						  (-1 -1) (0 -1) (1 -1)) collect
			 (when (= e 1)
			   (destructuring-bind (x y) p
			     (unless (= 0 x y)
			       ;; j is upside down
			       `(<= (aref im (- j ,y) (+ i ,x)) 
				   (aref im j i))))))))))
      (let ((res ())
	    (m 1)) 
	;; 4 edges
	(let ((i 0) (j 0)) ;; up left
	  (when (compare-expand 0 0 0
				0 1 1
				0 1 1)
	    (push (list j i (aref im j i)) res)))
	(let ((i (1- w)) (j 0)) ;; up right
	  (when (compare-expand 0 0 0
				1 1 0
				1 1 0)
	    (push (list j i (aref im j i)) res)))
	(let ((i 0) (j (1- h))) ;; down left
	  (when (compare-expand 0 1 1
				0 1 1
				0 0 0)
	    (push (list j i (aref im j i)) res)))
	(let ((i (1- w)) (j (1- h))) ;; down right
	  (when (compare-expand 1 1 0
				1 1 0
				0 0 0)
	    (push (list j i (aref im j i)) res)))
	;; left and right column
	(loop for j from m below (- h m 1) do
	     (let ((i 0))
	       (when (compare-expand 0 1 1
				     0 1 1
				     0 1 1)
		 (push (list j i (aref im j i)) res)))
	     (let ((i (1- w)))
	       (when (compare-expand 1 1 0
				     1 1 0
				     1 1 0)
		 (push (list j i (aref im j i)) res))))
	;; top and bottom row
	(loop for i from m below (- w m 1) do
	     (let ((j 0))
	       (when (compare-expand 0 0 0
				     1 1 1
				     1 1 1)
		 (push (list j i (aref im j i)) res)))
	     (let ((j (1- h)))
	       (when (compare-expand 1 1 1
				     1 1 1
				     0 0 0)
		 (push (list j i (aref im j i)) res))))
	(loop for j from m below (- h m 1) do
	     (loop for i from m below (- w m 1) do
		  (when (compare-expand 1 1 1
					1 1 1
					1 1 1)
		    (push (list j i (aref im j i)) res))))
	res))))

;; example for 5px neighbourhood
;; 0 0 0 0 0
;; 0 0 0 0 0 
;; 0 0 x 0 0 
;; 0 0 0 0 0 
;; 0 0 0 0 0
;; (floor 5 2) = 2
(defun extract-square (img y x &key (n 9))
  (destructuring-bind (h w) (array-dimensions img)
    (let* ((a (make-array (list n n)
			  :element-type (array-element-type img)))
	   (offset (floor n 2)))
      ;; make sure that image contains full neighbourhood
      (when (and (<= offset x (- w offset 1)) 
		 (<= offset y (- h offset 1)))
	(dotimes (j n)
	  (dotimes (i n)
	    (setf (aref a j i)
		  (aref img
			(+ y j (- offset))
			(+ x i (- offset))))))
	a))))

(defun insert-rectangle (img kern y x)
  (destructuring-bind (hh ww) (array-dimensions kern)
    (let ((oh (floor hh 2))
	  (ow (floor ww 2)))
     (dotimes (j hh)
       (dotimes (i ww)
	 (setf (aref img (+ y (- j oh)) (+ x (- i ow)))
	       (aref kern j i)))))))



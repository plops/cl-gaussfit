(defpackage :utils
  (:use #:cl)
  (:export #:make-displaced-array
	   #:copy-img
	   #:extract-frame
	   #:img-op
	   #:img-apply-mask
	   #:find-local-maxima))

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
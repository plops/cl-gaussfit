;; nr3.pdf p. 433 selection
(declaim (optimize (speed 3)))
(defun select (ar k)
  "Given K in [0..n-1] modifies the array AR to contain the K+1
  smallest values in the front (in arbitrary order). Return value is
  AR[K] and therefore the K+1-smallest value."
  (declare (type (simple-array single-float 1) ar)
	   (fixnum k)
	   (values single-float &optional))
  (let* ((n (length ar))
	 (a 0d0)
	 (i 0)
	 (j 0)
	 (mid 0)
	 (l 0) ;; index to leftmost element
	 (ir (1- n))) ;; index to rightmost element
    (labels ((swap (p q)
	       (let ((h (aref ar p)))
		 (setf (aref ar p) (aref ar q)
		       (aref ar q) h))))
     (loop
	(if (<= ir (1+ l))
	    ;; active partition contains 1 or 2 el
	    (progn (when (and (eq ir (1+ l))
			      (< (aref ar ir) (aref ar l)))
		     (swap l ir))
		   (return-from select (aref ar k)))
	    ;; choose median of left center and right as partitioning element
	    (progn 
	      (setf mid (floor (+ l ir) 2))
	      (swap mid (1+ l))
	      ;; rearrange for ar[l]<=ar[l+1]<=ar[ir]
	      (when (< (aref ar ir) (aref ar l))
		(swap ir l))
	      (when (< (aref ar ir) (aref ar (1+ l)))
		(swap ir (1+ l)))
	      (when (< (aref ar (1+ l)) (aref ar l))
		(swap l (1+ l)))
	      (setf i (1+ l)
		    j ir
		    a (aref ar i))
	      (loop ;; scan up and down, put big elements right, small ones left
		 (loop do (incf i)
		      while (< (aref ar i) a))
		 (loop do (decf j)
		      while (< a (aref ar j)))
		 (when (< j i)
		   (return))
		 (swap i j))
	      (setf (aref ar (1+ l)) (aref ar j)
		    (aref ar j) a)
	      (when (<= k j)
		(setf ir (1- j)))
	      (when (<= j k)
		(setf l i))))))))

#+nil
(let* ((n 10000000) ;; .2s for 10e6
       (nh (floor n 2))
       (a (make-array n :element-type 'single-float)))
  (dotimes (i n)
    (setf (aref a i) (random .5)))
  (time 
   (select a nh)))
#+nil
(let ((k 4))
  (let ((a (let ((ls '(1 3.2 12 9 8 7 6 4 5 3 2 1)))
	     (make-array (length ls)
			 :element-type 'single-float
			:initial-contents (mapcar #'(lambda (x) (* 1s0 x)) ls)))))
    (select a k)
    a))
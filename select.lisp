;; nr3.pdf p. 433 selection
(declaim (optimize (speed 3)))
(defun select (ar k)
  "Given K in [0..n-1] modifies the array AR to contain the K+1
  smallest values in the front (in arbitrary order). Return value is
  AR[K] and therefore the K+1-smallest value."
  (declare (type (simple-array single-float 1) ar)
	   (type fixnum k)
	   (values single-float &optional))
  (let* ((n (length ar))
	 (a 0s0)
	 (i 0)
	 (j 0)
	 (mid 0)
	 (l 0) ;; index to leftmost element
	 (ir (1- n))) ;; index to rightmost element
    (labels ((swap (p q)
	       (let ((h (aref ar p)))
		 (setf (aref ar p) (aref ar q)
		       (aref ar q) h)
		 (values)))
	     (ensure< (p q)
	       (when (< (aref ar p) (aref ar q))
		 (swap q p))
	       (values)))
     (loop
	(if (<= ir (1+ l))
	    (progn ;; indices of both sides met, finished
	      (when (and (eq ir (1+ l)) ;; active partition contains 1 or 2 el
			 (< (aref ar ir) (aref ar l)))
		(swap l ir))
	      (return-from select (aref ar k)))
	    (progn 
	      ;; choose median of left center and right as
	      ;; partitioning element
	      (setf mid (floor (+ l ir) 2))
	      (swap mid (1+ l))
	      ;; rearrange for ar[l]<=ar[l+1]<=ar[ir]
	      (ensure< ir l)
	      (ensure< ir (1+ l))
	      (ensure< (1+ l) l)
	
	      (setf i (1+ l)
		    j ir
		    a (aref ar i))
	      (loop ;; scan up and down, put big elements right, small ones left
		 (loop do (incf i) while (< (aref ar i) a))
		 (loop do (decf j) while (< a (aref ar j)))
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
(time
 (dotimes (i 1000)
   (let* ((n (1+ (random 10900)))
	  (k (random n))
  
	  (a (make-array n
			 :element-type 'single-float
			 :initial-contents (loop for i below n collect
						(* 1s0 (random 1.12))))))
     (select a k)
     ;;(format t "~a~%" (list 'i i 'n n 'k k))
     (when (and (< (1+ k) (- n 1)))
       (let ((ma (reduce #'max
			 (subseq a 0 (1+ k))))
	     (mi (reduce #'min 
			 (subseq a (1+ k)))))
	 (unless (<= ma mi)
	   (error "there is a too big element on the left")))))))
;; ANL8074a.pdf documentation for minpack

;; how to do fitting on the gpu (they consider poisson noise but it
;; doesn't seem to improve precision, their method is very fast but
;; they don't give code)
;; /home/martin/ondrej-papers/papers/Quan2010OSA.pdf

;; fitting precision 2004 ober
;;  /home/martin/ondrej-papers/papers/sdarticle\ \(3\).pdf

(defvar *rand* .9)

(defun g (x y x0 y0 a b s)
  (let* ((dx (- x x0))
	 (dy (- y y0))
	 (s2 (/ (* s s)))
	 (arg (- (* s2 (+ (* dx dx) (* dy dy)))))
	 (e (exp arg)))
    (+ b (* a e))))


(defun fill-img ()
 (defparameter *img*
   (let ((a (make-array '(5 5) :element-type 'double-float)))
     (destructuring-bind (h w) (array-dimensions a)
       (dotimes (j h)
	 (dotimes (i w)
	   (setf (aref a j i) (g i j 2.3d0 2.5d0 .7d0 .005d0 2d0)))))
     a)))

;; load the data and swap endian
(defvar *imgs*
  (with-open-file (s "example-movie_64x64x30000.raw"
		     :direction :input
		     :element-type '(unsigned-byte 16))
    (let* ((z 30000)
	   (x 64)
	   (y 64)
	   (n (* z x y)) 
	   (a (make-array n :element-type '(unsigned-byte 16)))
	   (b (make-array (list z y x) :element-type '(unsigned-byte 16)))
	   (b1 (sb-ext:array-storage-vector b)))
      (read-sequence a s)
      (dotimes (i n)
	(setf (aref b1 i) (+ (ldb (byte 8 8) (aref a i))
			     (* 256 (ldb (byte 8 0) (aref a i))))))
      b)))

(declaim (optimize (speed 0) (safety 3) (debug 3)))

(defun find-max (stack)
  (declare (type (simple-array (unsigned-byte 16) 3) stack))
  (let* ((ma 0)
	 (a1 (array-storage-vector stack))
	 (n (length a1)))
    (declare (type (unsigned-byte 16) ma))
   (dotimes (i n)
     (setf ma (max ma (aref a1 i))))
   ma))

#+nil
(time (find-max *imgs*))

(defun make-histograms (stack)
  (declare (type (simple-array (unsigned-byte 16) 3) stack))
  (destructuring-bind (z y x) (array-dimensions stack)
    (declare (type fixnum z y x))
    (let* ((ma (1+ (find-max stack)))
	   (n 500)
	   (ni 100) ;; combine ni=100 images into histogram
	   (stack1 (array-storage-vector stack))
	   (hist (make-array (list (floor z ni) n) :element-type 'fixnum)))
       (declare (type (simple-array fixnum 2) hist)
		(type fixnum n ni ma))
       (dotimes (k z) 
	 (let ((ko (floor k ni)))
	   (dotimes (j y)
	     (dotimes (i x) 
	       (incf (aref hist ko
			   (floor (* n (aref stack k j i))
				  ma)))))))
       hist)))
#+nil
(time
 (defparameter *hists* (make-histograms *imgs*)))

;; definition of fits file
;; ftp://nssdc.gsfc.nasa.gov/pub/fits
;; 80 byte per header entry
;; from imagej source
;; ~/src/mm/3rdpartypublic/sourceext/ij_source/ij/plugin/FTIS_Reader.java
;; start with SIMPLE = T / comment
;; = separates key value pair
;; / indicates comment
;; END stops
;; BITPIX = 16 short 32 int -32 float, -64 double
;; NAXIS1 = 129 / width, NAXIS2, NAXIS3
;; BSCALE = 1.34 / should be 1.0, then no calibration in imagej
;; BZERO = 32768 / for ushort (expt 2 15)
;; probably 2147483648 = (expt 2 31) for uint
;; count such lines (each 80 bytes long)
;; offset with data at 2880+2880*(((count*80)-1)/2880)
;; origin at bottom left

;; offset in example file #x1680 = 5760 = 2*2880

(defun generate-fits-header (img &key (bits 16))
  (let* ((dims (array-dimensions img))
	 (dat `((SIMPLE T)
		(BITPIX ,bits)
		(NAXIS ,(length dims))
		(NAXIS1 ,(first (last dims)))
		(NAXIS2 ,(first (last (butlast dims))))
		,(when (= (length dims) 3)
		       (list 'NAXIS3 (first dims))))))
    (let* ((head (with-output-to-string (s)
		   (labels ((ensure-80-chars (str)
			      (format s "~80A" str)))
		     (dolist (e dat)
		       (when e 
			 (ensure-80-chars ;; the imagej writer puts = at 9th position and value starts at 11
			  (format nil "~8A= ~A" (first e) (second e)))))
		     (ensure-80-chars "END"))
		   s))
	   (pad (make-array (- 2880 (length head))
			    :element-type 'character
			    :initial-element #\Space)))
      (concatenate 'string head pad))))

#+nil
(format t "~a~%"
 (generate-fits-header (first *blur*) :bits 16))

(defun write-fits (filename img)
  (let ((head (generate-fits-header img
				    :bits (cond
					    ((eq (array-element-type img) 'single-float) -32)
					    ((or
					      (eq (array-element-type img) '(unsigned-byte 16))
					      (eq (array-element-type img) '(signed-byte 16))) 16)))))
    (with-open-file (s filename
		       :direction :output
		       :if-exists :supersede
		       :if-does-not-exist :create)
      (write-sequence head s))
    (let* ((i1 (array-displacement (copy-img img)))
	   (sap (sb-sys:vector-sap i1))
	   (n (array-total-size img))
	   (h1 (make-array n :element-type '(unsigned-byte 32)))
	   (h2 (make-array n :element-type '(unsigned-byte 32))))
      (dotimes (i n)
	(setf (aref h1 i) (sb-sys:sap-ref-32 sap (* 4 i))))
      (cond ;; swap endian
	((eq (array-element-type img) 'single-float)
	 (dotimes (i n) 
	   (setf (ldb (byte 8 0) (aref h2 i)) (ldb (byte 8 24) (aref h1 i))
		 (ldb (byte 8 8) (aref h2 i)) (ldb (byte 8 16) (aref h1 i))
		 (ldb (byte 8 16) (aref h2 i)) (ldb (byte 8 8) (aref h1 i))
		 (ldb (byte 8 24) (aref h2 i)) (ldb (byte 8 0) (aref h1 i)))))
	((or (eq (array-element-type img) '(unsigned-byte 16))
	     (eq (array-element-type img) '(signed-byte 16)))
	 (dotimes (i n) 
	   (setf (ldb (byte 8 0) (aref h2 i)) (ldb (byte 8 8) (aref h1 i))
		 (ldb (byte 8 8) (aref h2 i)) (ldb (byte 8 0) (aref h1 i))
		 (ldb (byte 8 16) (aref h2 i)) (ldb (byte 8 24) (aref h1 i))
		 (ldb (byte 8 24) (aref h2 i)) (ldb (byte 8 16) (aref h1 i))))))
      (with-open-file (s filename
			 :element-type '(unsigned-byte 32)
			 :direction :output
			 :if-exists :append)
	(write-sequence h2 s))))
  (values))

(defun scale-log (img)
  (let* ((n (array-total-size img)) 
	 (i1 (make-array n
			 :element-type (array-element-type img)
			 :displaced-to img))
	(ma (1+ (reduce #'max i1))))
    (dotimes (i n)
      (let ((v (aref i1 i)))
       (setf (aref i1 i) (if (= v 0)
			     0
			     (floor (* (expt 2 15) (log v))
				    (log ma))))))
    img))

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

#+nil
(time
 (defparameter *hists*
   (destructuring-bind (z y x) (array-dimensions *imgs*)
     (let* ((ma (1+ (find-max *imgs*)))
	    (n 400)
	    (ni 30) ;; combine ni=100 images into histogram
	    (hist (make-array (list (ceiling z ni) n) :element-type 'fixnum)))
       (dotimes (k z)
	 (dotimes (j y)
	   (dotimes (i x)
	     (incf (aref hist (floor k ni) 
			 (floor (* n (aref *imgs* k j i)) ma))))))
       (let ((ma (reduce #'max (array-storage-vector hist)))) 
	 (when (< (1- (expt 2 16)) ma)
	  (error "maximum count ~d in histogram doesn't fit into 16 bit" ma)))
       (let* ((a (make-array (array-dimensions hist)
			     :element-type '(unsigned-byte 16)))
	      (a1 (array-storage-vector a))
	      (h1 (array-storage-vector hist)))
	 (dotimes (i (length a1))
	   (setf (aref a1 i) (aref h1 i)))
	 a)))))


;; i now want to separate background and signal. ideally i would just
;; create masks to select background, but i don't see a robust way to
;; do that. instead i will subtract images. i expect the histogram of
;; those to contain a gaussian distribution centered on 0. its sigma
;; should be the same as the background sigma.

;; i hope that the signal's contribution is separated a bit better in
;; the histogram. for this i shouldn't subtract consecutive images but
;; ensure, that molecules are in different positions in both images.

;; eventually i would like to know the sigma and mean of the
;; background, so that i can select events that are brighter than,
;; i.e. 5 sigma


;; calculate average over each image
#+nil
(defparameter *imgs-avg*
  (destructuring-bind (z y x) (array-dimensions *imgs*)
    (loop for k below z collect
	 (let ((sum 0))
	   (dotimes (j y)
	     (dotimes (i x)
	       (incf sum (aref *imgs* k j i))))
	   sum))))

;; accumulate histograms of the averages to determine which images have events
#+nil
(defparameter *imgs-avg-hist*
 (let* ((ma (1+ (reduce #'max *imgs-avg*)))
	(n 10)
	(hist (make-array n :element-type 'fixnum)))
   (dolist (e *imgs-avg*)
     (incf (aref hist (floor (* n e) ma))))
   (loop for j from 0 and i across hist collect
	(list (floor (* ma (/ 1d0 n) j)) i))))

;; select images with bright events
#+nil
(defparameter *imgs-with-event*
  (let ((ma (1+ (reduce #'max *imgs-avg*))))
    (loop for k from 0 and e in *imgs-avg* 
       when (< 49000 e) ;(< .51 (/ e	ma))
       collect k)))

;; show avgs as well
#+nil
(defparameter *im*
 (loop for e in *imgs-more-events* collect
      (list e (elt *imgs-avg* e))))

;; plot data
#+nil
(progn
 (with-open-file (s "/dev/shm/o.dat"
		    :direction :output
		    :if-exists :supersede
		    :if-does-not-exist :create)
   (format s "~{~{~6d ~6d~}~%~}" *im*))
 (with-open-file (s "/dev/shm/o.gp"
		    :direction :output
		    :if-exists :supersede
		    :if-does-not-exist :create)
   (format s "plot \"/dev/shm/o.dat\" u 1:2 w l; pause -1"))
 (run-program "/usr/bin/gnuplot" '("/dev/shm/o.gp")))


(defun extract-frame (img k)
  (destructuring-bind (z y x) (array-dimensions img)
    (let* ((start (* x y k))
	   (i1 (make-array (* x y) :element-type (array-element-type img)
			   :displaced-to img
			   :displaced-index-offset start))
	   (c1 (subseq i1 0)))
      (make-array (list y x) 
		  :element-type (array-element-type img)
		  :displaced-to c1))))

(defun make-displaced-array (a)
  (make-array (array-total-size a)
	      :element-type (array-element-type a)
	      :displaced-to a))

(defun img-op (op a b)
  (let* ((a1 (make-displaced-array a))
	 (b1 (make-displaced-array b))
	 (c (make-array (array-dimensions a)
			:element-type 'single-float))
	 (c1 (make-displaced-array c)))
    (dotimes (i (length a1))
      (setf (aref c1 i) (* 1s0 (funcall op (aref a1 i) (aref b1 i)))))
    c))

(defun img-mask (img mask &key (background 0))
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



;; ~/src/mm/3rdpartypublic/sourceext/ij_source/ij/plugin/filter/GaussianBlur.java


#+nil
(time
 (let ((start 673)
       (n 100))
   (defparameter *raw*
     (loop for k from start upto (+ start n) collect
	  (ub16->single-2 (extract-frame *imgs* k))))
   (defparameter *blur*
     (let ((imgs nil))
       (loop for k from start upto (+ start n) collect
	    (let* ((s 1.3)
		   (s2 1.5)
		   (dog (img-op #'-
				(blur-float
				 (ub16->single-2
				  (extract-frame *imgs* k))
				 s s 1e-4)
				(blur-float
				 (ub16->single-2
				  (extract-frame *imgs* k)) 
				 s2 s2 1e-4))))
	      dog
	      #+nil (mark-points dog (find-local-maxima dog))))))
   (write-fits "/dev/shm/o.fits" (img-list->stack *blur*))))

;; use corrected two pass algorithm 14.1.8 numerical recipes in C
;; 1/(N-1) * (sum_j (x_j-m)^2 - 1/N (sum_j (x_j-m))^2) 
;; second term corrects round-off error

(defun get-statistics (pic-list &key mask-list)
  (if mask-list
      (let* ((mean (let ((sum 0)
			 (n 0))
		     (loop for e in pic-list and m in mask-list do
			  (let ((e1 (make-displaced-array e))
				(m1 (make-displaced-array m)))
			    (loop for g across e1 and mask across m1 
			       when mask
			       do (incf sum g)
				 (incf n))))
		     (* sum (/ 1d0 n))))
	     (s1 (let ((n 0)
		       (sum 0)) 
		   (loop for e in pic-list and m in mask-list do
			(let ((e1 (make-displaced-array e))
			      (m1 (make-displaced-array m)))
			  (loop for g across e1 and mask across m1 
			     when mask
			     do (incf n) 
			       (incf sum (expt (- g mean) 2)))))
		   sum))
	     (n 0)
	     (s2 (let ((sum 0)) 
		    (loop for e in pic-list and m in mask-list do
			 (let ((e1 (make-displaced-array e))
			       (m1 (make-displaced-array m)))
			   (loop for g across e1 and mask across m1 
			      when mask
			      do (incf n) 
				(incf sum (- g mean)))))
		    sum))

	     (var (/ (- s1 (/ (expt s2 2) n))
		     (1- n))))
	
	(format t "~a" (list 'mean mean (* 1d0 mean) 's1 s1 'stddev1 (sqrt (/ s1 (1- n))) 's2 s2 'var var 'stddev (sqrt var)			     ))
	(values mean (sqrt var)))
      (let* ((n (* (length pic-list) (array-total-size (first pic-list))))
	     (mean (/ (loop for e in pic-list sum
			   (reduce #'+ (make-displaced-array e)))
		      n))
	     (s1 (loop for e in pic-list sum 
		      (loop for g across (make-displaced-array e) sum
			   (expt (- g mean) 2))))
	     (s2 (loop for e in pic-list sum
		      (loop for g across (make-displaced-array e) sum
			   (- g mean))))
	     (var (/ (- s1 (/ (expt s2 2) n))
		     (1- n))))
	(format t "~a" (list 'mean mean (* 1d0 mean) 's1 s1 'stddev1 (sqrt (/ s1 (1- n))) 's2 s2 'var var 'stddev (sqrt var)
			     ))
	(values mean (sqrt var)))))

;; find histogram and statistics of difference of gaussian images
#+nil
(time
 (destructuring-bind (h w) (array-dimensions (first *blur*))
   (let* ((ma (loop for e in *blur* maximize
		   (reduce #'max (make-displaced-array e))))
	  (mi (loop for e in *blur* minimize
		   (reduce #'min (make-displaced-array e))))
	  (mp (* 1.1 (max ma (abs mi))))
	  (n 600)
	  (hist (make-array n :element-type 'fixnum))
	  ee q)
     (loop for e in *blur* do
	  (multiple-value-setq (ee q)
	    (calc-hist e :n n :minv (* 1.1 mi) :maxv (* 1.1 ma) :append hist)))
     (defparameter *dog-mean* 0)
     (defparameter *dog-stddev* 0)
     (multiple-value-setq (*dog-mean*
			   *dog-stddev*)
       (get-statistics *blur*))
     (with-open-file (s "/dev/shm/dog-all.dat" :direction :output
			:if-exists :supersede
			:if-does-not-exist :create)
      (format s "~{~{~f ~f~%~}~} # mean = ~f stddev = ~f~%"
	      (loop for i across hist and g in q 
		 collect (list g 
			       (if (= 0 i) 0 i)))
	      *dog-mean*
	      *dog-stddev*)))))

;; remove maxima that are not brighter than .5 stddev
#+nil
(time
  (progn
    (defparameter *blur-big-ma* ())
    (defparameter *blur-ma*
      (loop for e in *blur* collect
	   (let* ((c (copy-img e))
		  (ma (find-local-maxima c))
		  (big-ma (remove-if #'(lambda (e) 
					 (< (third e)
					    (+ *dog-mean* 
					       (* .5 *dog-stddev*))))
				     ma)))
	     (push big-ma *blur-big-ma*)
	     (mark-points c big-ma))))
    (setf *blur-big-ma* (reverse *blur-big-ma*))
    (write-fits "/dev/shm/blur-ma.fits" (img-list->stack *blur-ma*))))

(defun draw-mask (img points &key (mask-border 0))
  (destructuring-bind (h w) (array-dimensions img)
   (let ((a (make-array (array-dimensions img)
			:element-type 'boolean
			:initial-element nil)))
     (loop for j from mask-border below (- h mask-border) do
	  (loop for i from mask-border below (- w mask-border) do
	       (setf (aref a j i) t)))
     (dolist (p points)
       (destructuring-bind (y x val) p
	 (loop for j from -5 upto 5 do
	   (loop for i from -5 upto 5 do
		(let ((xx (+ x i))
		      (yy (+ y j)))
		  (when (and (<= 0 xx (1- w))
			     (<= 0 yy (1- h)))
		    (setf (aref a yy xx) nil)))))))
     a)))

;; mask 10x10 area around bright maxima in dog images
#+nil
(time
 (progn 
   (with-open-file (s "/dev/shm/hists.gp" :direction :output
		     :if-exists :supersede
		     :if-does-not-exist :create)
    (format s "set log y
plot \"raw.dat\" u 1:2 w l, \"raw-all.dat\" u 1:2 w l, \"dog-all.dat\" u 1:2 w l, \"dog.dat\" u 1:2 w l
pause -1
"))
   (let* ((ma (loop for e in *blur* maximize
		   (reduce #'max (make-displaced-array e))))
	  (mi (loop for e in *blur* minimize
		   (reduce #'min (make-displaced-array e))))
	  (n 600)
	  (hist (make-array n :element-type 'fixnum))
	  ee qq
	  (masks
	   (loop for im in *blur* and points in *blur-big-ma* collect
		(let ((m (draw-mask im points :mask-border 0)))
		  (multiple-value-setq (ee qq)
		    (calc-hist im :n n 
			       :minv (* (if (< mi 0) 1.1 .9) mi)
			       :maxv (* (if (< ma 0) 1.1 .9) ma)
			       :append hist :mask m))
		  m)))) 
     
     (multiple-value-bind (mean stddev) (get-statistics *blur* :mask-list masks)
       (with-open-file (s "/dev/shm/dog.dat" :direction :output
			  :if-does-not-exist :create :if-exists :supersede)
	 (format s "~{~{~f ~f~%~}~}~% # mean = ~f stddev = ~f~%"
		 (loop for i across hist and g in qq
		    collect (list g 
				  (if (= 0 i) 0 i)))
		 mean
		 stddev))
       (defparameter *blur-10x10masked* (loop for e in *blur* and m in masks
					   collect
					     (img-mask (copy-img e) m
						       :background mean))))
     (write-fits "/dev/shm/blur-10x10masked.fits" 
		 (img-list->stack *blur-10x10masked*)))))

#+nil
(time
 (format t "~A~%"
	 (progn 
	   (defparameter *dog-mean-2* 0)
	   (defparameter *dog-stddev-2* 0)
	   (let ((mask-list
		  (loop for i below (length *blur*) collect
		       (draw-mask (first *blur*) 
				  (elt *blur-big-ma* i) :mask-border 0))))
	     (multiple-value-setq (*dog-mean-2*
				   *dog-stddev-2*) 
	       (get-statistics *blur* :mask-list mask-list))))))

;; select maxima again, but use better estimate of background stddev
#+nil
(time
  (progn
    (defparameter *blur-big-ma-2* ())
    (defparameter *blur-ma-2*
      (loop for e in *blur* collect
	   (let* ((c (copy-img e))
		  (ma (find-local-maxima c))
		  (big-ma (remove-if #'(lambda (e) 
					 (< (third e)
					    (+ *dog-mean-2* 
					       (* 3 *dog-stddev-2*))))
				     ma)))
	     (push big-ma *blur-big-ma-2*)
	     (mark-points c big-ma))))
    (setf *blur-big-ma-2* (reverse *blur-big-ma-2*))
    (write-fits "/dev/shm/blur-ma-2.fits" (img-list->stack *blur-ma-2*))))

;; get histogram for the full raw data
#+nil
(time
 (let* ((ims (loop for e in *raw* collect (ub16->single-2 e)))
	(ma (loop for e in ims maximize
		 (reduce #'max (make-displaced-array e))))
	(mi (loop for e in ims minimize
		 (reduce #'min (make-displaced-array e))))
	(n 600)
	(hist (make-array n :element-type 'fixnum))
	ee qq
	)
   (loop for e in ims collect
	(multiple-value-setq (ee qq)
	  (calc-hist e :n n :minv (* .9 mi) 
		     :maxv (* 1.1 ma) :append hist)))

   (multiple-value-bind (mean stddev)
       (get-statistics ims)
     (with-open-file (s "/dev/shm/raw-all.dat" :direction :output
			:if-does-not-exist :create :if-exists :supersede)
       (format s "~{~{~f ~f~%~}~} # mean = ~f stddev = ~f~%"
	       (loop for i across hist and g in qq
		  collect (list g 
				(if (= 0 i) 0 i)))
	       mean stddev)))))

;; get the histogram for the masked raw data
#+nil
(time
 (let* ((ims (loop for e in *raw* collect (ub16->single-2 e)))
	(ma (loop for e in ims maximize
		 (reduce #'max (make-displaced-array e))))
	(mi (loop for e in ims minimize
		 (reduce #'min (make-displaced-array e))))
	(n 600)
	(hist (make-array n :element-type 'fixnum))
	ee qq
	(masks (loop for i below (length ims) collect
		    (let ((e (elt ims i))
			  (m (draw-mask (first ims)
					(elt *blur-big-ma-2* i) :mask-border 0)))
		      
		      (multiple-value-setq (ee qq)
			(calc-hist e :n n :minv (* .9 mi) 
				   :maxv (* 1.1 ma) :append hist :mask m))
		      m))))
   
   (multiple-value-bind (mean stddev)
       (get-statistics ims :mask-list masks)
     (let ((masked-images 
	    (loop for e in ims and m in masks collect
		 (img-mask (ub16->single-2 e) m
			   :background mean))))
       (with-open-file (s "/dev/shm/raw.dat" :direction :output
			  :if-does-not-exist :create :if-exists :supersede)
	 (format s "~{~{~f ~f~%~}~} # mean = ~f stddev = ~f~%"
		 (loop for i across hist and g in qq
		    collect (list g 
				  (if (= 0 i) 0 i)))
		 mean stddev))
       (defparameter *raw-10x10masked* masked-images))
     (write-fits "/dev/shm/raw-marked.fits"
		 (img-list->stack 
		  (loop for e in ims and points in *blur-big-ma-2* collect
		       (mark-points (copy-img e) points
				    :value (+ mean (* -5 stddev))))))
     (write-fits "/dev/shm/raw-10x10masked.fits" 
		 (img-list->stack *raw-10x10masked*)))))

;; example for 5px neighbourhood
;; 0 0 0 0 0
;; 0 0 0 0 0 
;; 0 0 x 0 0 
;; 0 0 0 0 0 
;; 0 0 0 0 0
;; (floor 5 2) = 2
(defun extract (img y x &key (n 9))
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

(defun insert (img kern y x)
  (destructuring-bind (hh ww) (array-dimensions kern)
    (let ((oh (floor hh 2))
	  (ow (floor ww 2)))
     (dotimes (j hh)
       (dotimes (i ww)
	 (setf (aref img (+ y (- j oh)) (+ x (- i ow)))
	       (aref kern j i)))))))

(defun img-mul (im &optional (factor .001d0))
  (let* ((a (make-array (array-dimensions im)
			:element-type (array-element-type im)))
	 (i1 (make-displaced-array im))
	 (a1 (make-displaced-array a)))
    (dotimes (i (length a1))
      (setf (aref a1 i) (* factor (aref i1 i))))
    a))

#+nil
(length *all-fits*)

#+nil
(time
 (progn
   (defparameter *all-fits* nil)
   (loop for k from 0 below (length *raw*) do
	(progn
	  (progn ;; extract areas around peaks in an image
	    (defparameter *current-image* k)
	    (defparameter *raw-blobs*
	      (let ((im (elt *raw* *current-image*)))
		(loop for e in (elt *blur-big-ma-2* *current-image*) collect
		     (destructuring-bind (j i val) e
		       (extract im j i :n 5)))))
	    (write-fits "/dev/shm/raw-blobs.fits"
			(img-list->stack (remove-if #'null *raw-blobs*))))

	  (progn ;; run gauss fits on the extracted images
	    (defparameter
		*fits*
	      (loop for e in *raw-blobs* and point in 
		   (elt *blur-big-ma-2* *current-image*) collect
		   (when e
		     (setf *img* (img-mul e .001))
		     (destructuring-bind (h w) (array-dimensions e)
		       (destructuring-bind (j i val) point
			 (multiple-value-bind (x err fnorm)
			     (fit-gaussian :x0 (floor w 2) :y0 (floor h 2)
					   :a 1 :b .49 :sigma .8)
			   (list
			    fnorm val
			    (loop for e across x and r in err collect
				 (list e r)))))))))
	    (loop for num from 0 and (fnorm val x+err) in (remove-if #'null *fits*) do
		 (format t "~3d ~4,2f ~4,0f ~{~{~7,2f ~3,2f~}~}~%"
			 num fnorm val x+err))
	    (setf *all-fits* (append *all-fits* (list *fits*))))
	  
	  #+nil
	  (progn ;; generate images of the fitted functions
	    (defparameter *fit-calc*
	      (loop for e in *fits* collect
		   (when e
		     (destructuring-bind (fnorm val x+err) e
		       (destructuring-bind ((x dx) (y dy) (a da) (b db) (s ds)) x+err
			 (let* ((n 5) 
				(ar (make-array (list n n) :element-type 'single-float)))
			   (dotimes (j n)
			     (dotimes (i n)
			       (setf (aref ar j i) (coerce (g i j x y a b s)
							   'single-float))))
			   ar))))))
	    (let ((a (make-array (array-dimensions (first *raw*))
				 :element-type 'single-float
				 :initial-element .49)))
	      (loop for f in *fit-calc* and pos in (elt *blur-big-ma-2* *current-image*) do
		   (when f
		     (destructuring-bind (j i val) pos
		       (insert a f j i))))
	      (write-fits "/dev/shm/fit-calc.fits" 
			  (img-list->stack
			   (let ((orig (img-mul (elt *raw* *current-image*)
						.001)))
			     (list a
				   orig
				   (img-op #'- orig a)))))))))))



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
#+nil
(let ((ma (find-local-maxima (first *blur*))))
 (format t "~{~{~3d ~3d ~7,4f~%~}~}~% ( ~a )~%"
	 (sort 
	  ma
	  #'> :key #'third)
	 (length ma)))

(defun mark-points (img ls &key (value .1))
  (let ((v (coerce value (array-element-type img))))
   (dolist (e ls)
     (destructuring-bind (j i amp) e
       (setf (aref img j i) v))))
  img)

(defun img-list->stack (ls)
  (destructuring-bind (h w) (array-dimensions (first ls))
   (let* ((z (length ls))
	  (vol (make-array (list z h w)
			  :element-type (array-element-type (first ls)))))
     (dotimes (k z)
       (let ((im (elt ls k)))
	(dotimes (j h)
	  (dotimes (i w)
	    (setf (aref vol k j i) (aref im j i))))))
     vol)))

(defun blur-float (img sigma-x sigma-y accuracy)
  (when (< 0 sigma-x)
    (blur-1-direction img sigma-x accuracy T))
  (when (< 0 sigma-y)
    (blur-1-direction img sigma-y accuracy nil))
  img)

(defun blur-1-direction (img sigma accuracy x-direction)
  (destructuring-bind (h w) (array-dimensions img)
    (let* ((length (if x-direction w h))
	   (pixels (make-displaced-array img))
	   (point-inc (if x-direction 1 w))
	   (line-inc (if x-direction w 1))
	   (line-from 0)
	   (line-to (if x-direction h w))
	   (write-from 0)
	   (write-to (if x-direction w h))
	   (gauss-kernel (make-gaussian-kernel :sigma sigma :accuracy accuracy))
	   (kradius (array-dimension gauss-kernel 1))
	   (cache (make-array length :element-type 'single-float))
	   (pixel0 (* line-from line-inc))
	   (read-from (max 0 (- write-from kradius)))
	   (read-to (min length (+ write-to kradius))))
      (loop for line from line-from below line-to do
	   (let ((p (+ pixel0 (* read-from point-inc))))
	     (loop for i from read-from below read-to do
		  (setf (aref cache i) (aref pixels p))
		  (incf p point-inc))
	     (convolve-line cache pixels gauss-kernel
			    write-from write-to pixel0 point-inc))
	   (incf pixel0 line-inc)))))

(defun convolve-line (in pixels kernel write-from write-to point0 point-inc)
  (declare (type (simple-array single-float 2) kernel))
  (macrolet ((kern (i)
	       `(aref kernel 0 ,i))
	     (kern-sum (i)
	       `(aref kernel 1 ,i))
	     (do-inside ()
	       `(loop for k from 1 below kradius do
		     (let ((v 0))
		       (when (<= 0 (- i k)) 
			 (incf v (aref in (- i k))))
		       (when (< (+ i k) len)
			 (incf v (aref in (+ i k))))
		       (incf res (* (kern k) v)))))
	     (do-updates ()
	       `(progn (setf (aref pixels p) (* 1s0 res))
		       (incf p point-inc)
		       (incf i)))
	     (with-res (&body body)
	       `(let ((res (* (aref in i)
			      (kern 0))))
		  (when (< i kradius)
		    (incf res (* (kern-sum i) first)))
		  (when (<= len (+ i kradius)) ;; i think <= is correct in both loops 
		    (incf res (* (kern-sum (- len i 1))
				 last)))
		  ,@body)))
    (let* ((len (length in))
	   (first (aref in 0)) ;; replace out-of-edge with nearest edge pixel
	   (last (aref in (1- len)))
	   (kradius (array-dimension kernel 1))
	   (first-part (min kradius len))
	   (p (+ point0 (* write-from point-inc)))
	   (i write-from))
      (loop while (< i first-part) ;; while sum includes pixels < 0
	 do 
	   (with-res
	       (do-inside)
	     (do-updates)))
      
      (let ((i-end-inside (min (- len kradius) write-to)))
	
	(loop while (< i i-end-inside) ;; easy case: address only pixels within the line 
	   do
	     (let ((res (* (aref in i) (kern 0))))
	       (loop for k from 1 below kradius do
		    (incf res (* (kern k)
				 (+ (aref in (- i k))
				    (aref in (+ i k))))))
	       (do-updates)))
	(setf i i-end-inside))
      (loop while (< i write-to) ;; sum includes pixels >= len
	 do
	   (with-res
	       (do-inside)
	     (do-updates)))))
  (values))

(defun make-gaussian-kernel (&key (sigma .6) (accuracy .001))
  "sigma .. radius of decay to 1/e in pixels
accuracy .. 1e-3 to 1e-4 for 16bit images. smaller is better"
  (let* ((n (1+ (ceiling (* sigma (sqrt (* -2 (log accuracy)))))))
	 (s (/ 1s0 (* sigma sigma)))
	 (a (make-array (list 2 n) :element-type 'single-float)))
    (dotimes (i n)
      (setf (aref a 0 i) (exp (* -.5 i i s))))
    (when (< 3 n) ;; edge correction to avoid cutoff at finite value
      ;; second order polynomial is zero at first out-of-kernel pixel
      (let* ((sqrt-slope single-float-positive-infinity)
	     (r (loop named slope for r from (1- n) downto (floor n 2) do
		     (let ((v (/ (sqrt (aref a 0 r))
				 (- n r))))
		       (if (< v sqrt-slope)
			   (setf sqrt-slope v)
			   (return-from slope r))))))
	(loop for r1 from (+ r 2) below n do
	     (setf (aref a 0 r1) (expt (* sqrt-slope (- n r1)) 2)))))
    (let* ((sum/ (/ 1s0 (+ (aref a 0 0)
			(* 2 (loop for i from 1 below n sum (aref a 0 i))))))
	   (rsum (+ .5 (* .5 sum/ (aref a 0 0)))))
      (defparameter *rsum* (list 'sum (/ s) 'rsum rsum))
      (loop for i below n do
	   (let ((v (* sum/ (aref a 0 i))))
	    (setf (aref a 0 i) v)
	    (decf rsum v)
	    (setf (aref a 1 i) rsum))))
    a))
#+nil
(make-gaussian-kernel :sigma 1.2 :accuracy 1e-4)
;; i compiled the java code in c and this is the result
;; 0 sum 0.332452 rsum 0.333774
;; 1 sum 0.234927 rsum 0.098847
;; 2 sum 0.082898 rsum 0.015950
;; 3 sum 0.014607 rsum 0.001343
;; 4 sum 0.001285 rsum 0.000058
;; 5 sum 0.000056 rsum 0.000001
;; 6 sum 0.000001 rsum 0.000000

#+nil
(make-gaussian-kernel :sigma 5.4 :accuracy 1e-4)

(defun calc-hist (img &key (n 30) minv maxv (append nil) mask)
  (let* ((i1 (make-displaced-array img))
	 (ma (if maxv maxv (1+ (ceiling (reduce #'max i1)))))
	 (mi (if minv minv (1- (floor (reduce #'min i1)))))
	 (hist (if append
		   append
		   (make-array n :element-type 'fixnum)))
	 (g (loop for x below n collect
		 (+ mi
		    (* x (/ 1d0 n) (- ma mi))))))
    
    ;; x = n (g-mi)/(ma-mi)
    ;; (x (ma-mi) / n)+mi = g
    (if mask
	(let ((mask1 (make-displaced-array mask)))
	 (dotimes (i (length i1))
	   (when (aref mask1 i)
	     (incf (aref hist (floor (* n (- (aref i1 i) mi))
				     (- ma mi)))))))
	(dotimes (i (length i1))
	  (incf (aref hist (floor (* n (- (aref i1 i) mi))
				  (- ma mi))))))
    (values hist g)))

#+nil ;; histogram of subtracted images
(time
 (let* ((n 300)
	(mi (1- (floor (loop for e in *subtr* minimize
			    (reduce #'min (make-displaced-array e))))))
	(ma (1+ (floor (loop for e in *subtr* maximize
			    (reduce #'max (make-displaced-array e))))))
	(mp (max (abs mi) (abs ma)))
	(h (make-array n :element-type 'fixnum))
	(q nil)
	(e nil))
   (dolist (e *subtr*)
     (multiple-value-setq (e q)
       (calc-hist e :n n :append h :minv (- mp) :maxv mp)))
   (format t "~{~{~7,0f ~5,2f~%~}~}"
	   (loop for i across h and g in q 
	      collect (list g 
			    (if (= 0 i) 0 i))))))


;; print out the images
#+nil
(destructuring-bind (z y x) (array-dimensions *imgs*)
  (loop for e in *imgs-more-events* do
       (format t "* ~6,'0d *~%" e)
       (let ((ma (reduce #'max
			 (array-displacement
			  (extract-frame *imgs* e)))))
	 (dotimes  (j y)
	   (dotimes (i x)
	     (format t "~d " (floor (* 9 (aref *imgs* e j i))
				    ma)))
	   (terpri)))
       (terpri)))

(defun find-val-in-image (img val)
  (destructuring-bind (y x) (array-dimensions img)
    (dotimes (j y)
      (dotimes (i x)
	(when (= (aref img j i) val)
	  (return-from find-val-in-image (values j i)))))))

;; frame 673 ff contains a maximum not too close to border at x=2 y=7
(defun locate-max-in-frame (stack frame-index)
 (let* ((frame (extract-frame stack frame-index))
	(ma (reduce #'max (array-displacement frame))))
   (find-val-in-image frame ma)))
#+nil
(locate-max-in-frame *imgs* 673)

(defun ub16->sb16 (img)
  (let* ((a (make-array (array-dimensions img)
			:element-type '(signed-byte 16)))
	 (a1 (array-storage-vector a))
	 (n (length a1))
	 (i1 (array-displacement img)))
    (dotimes (i n)
      (setf (aref a1 i) (coerce (aref i1 i) '(signed-byte 16))))
    a))

(defun ub16->double-2 (img)
  (if (eq (array-element-type img) 'double-float)
      img
      (let* ((a (make-array (array-dimensions img)
			    :element-type 'double-float))
	     (a1 (array-storage-vector a))
	     (n (length a1))
	     (i1 (array-displacement img)))
	(dotimes (i n)
	  (setf (aref a1 i) (* 1d0 (aref i1 i))))
	a)))

(defun ub16->single-2 (img)
  (if (eq (array-element-type img) 'single-float)
      img
      (let* ((a (make-array (array-dimensions img)
			    :element-type 'single-float))
	     (a1 (array-storage-vector a))
	     (n (length a1))
	     (i1 (make-displaced-array img)))
	(dotimes (i n)
	  (setf (aref a1 i) (* 1s0 (aref i1 i))))
	a)))

;; print 2 digits of the sigma sx and the same for the value x
;; if sx>1 print all digits
(defun print-with-error (x sx)
 (let* ((sigma-digits (- (floor (log sx 10))))
	(n (if (<= sigma-digits 0) 
	       0 ;; if sigma > 1 
	       (1+ sigma-digits))))
   (if (<= sigma-digits 0)
       (format nil "~4d ~4d" (floor x) (floor sx))
       (format nil (format nil "~~5,~df ~~5,~df" n n) x sx))))

#+nil
(print-with-error 2.2442 0.13)

#|
jacobian([b+a*exp(-((x-xx)^2+(y-yy)^2)/sigma^2)-f],[xx,yy,a,b,sigma]);
dx = x-xx
dy = y-yy
s2 = 1/s^2
arg = -(dx^2+dy^2) s2
e = exp(arg)
ff = a e
f(x) = a e + b
f_xx = 2 dx s2 ff
f_yy = 2 dy s2 ff
f_a  = e
f_b  = 1
f_s  = - 2 arg/s ff
|#

(defvar *img* nil)

(sb-alien::define-alien-callback fcn2
    void
    ((m (* int)) ; = 8x8 or similar
     (n (* int)) ; = 5
     (x (* double)) ; [xx,yy,a,b,sigma]
     (fvec (* double)) ; m
     (fjac (* double)) ; ldfjac,n
     (ldfjac (* int)) ; 1  ;; leading dimension of fjac
     (iflag (* int))) ; 1 1: fvec, 2: fjac
  (destructuring-bind (h w) (array-dimensions *img*)
    (let* ((xx (deref x 0))
	   (yy (deref x 1))
	   (a (deref x 2))
	   (b (deref x 3))
	   (s (deref x 4))) 
      #+nil (format t "bla ~a ~%" (list (deref iflag 0)
				  (loop for i below 5 collect (deref x i))))
      (ecase (deref iflag 0)
	(1 (dotimes (j h)
	     (dotimes (i w)
	       (let* ((dx (- (* 1d0 i) xx))
		      (dy (- (* 1d0 j) yy))
		      (s2 (/ (* s s))) ;; I just realize that my
				       ;; variable is not sigma:
				       ;; s^2 = 2 sigma^2 -> sigma = s/sqrt(2)
		      (arg (- (* s2 (+ (* dx dx) (* dy dy)))))
		      (e (exp arg)) 
		      (p (+ i (* w j))))
		 (setf (deref fvec p)
		       (+ (- (aref *img* j i)) b (* a e)))))))
	(2 (let ((ww (deref ldfjac 0)))
	     (macrolet ((f (a b)
			  `(deref fjac (+ ,a (* ww ,b)))))
	       (dotimes (j h)
		 (dotimes (i w)
		   ;; f_xx = 2 dx s2 f
		   ;; f_yy = 2 dy s2 f
		   ;; f_a  = e
		   ;; f_b  = 1
		   ;; f_s  = - 2 arg/s f	
		   (let* ((dx (- (* 1d0 i) xx))
			  (dy (- (* 1d0 j) yy))
			  (s2 (/ (* s s)))
			  (arg (- (* s2 (+ (* dx dx) (* dy dy)))))
			  (e (exp arg))
			  (f (* e a))
			  (fxx (* 2 dx s2 f))
			  (fyy (* 2 dy s2 f))
			  (fa e)
			  (fb 1d0)
			  (fs (/ (* -2 arg f)
				 s))
			  (p (+ i (* w j)))) ;; i (width) is the fast index
		     (setf (f p 0) fxx ;; the pixels are the fast index, the variables are slow
			   (f p 1) fyy
			   (f p 2) fa
			   (f p 3) fb
			   (f p 4) fs))))))))))  
  (values))


(load-shared-object "/usr/lib/libminpack.so")

(define-alien-routine lmder1_ void
  (fcn (* int)) ; callback
  (m int :copy) ; 1 input
  (n int :copy) ; 1
  (x (* double))			; n in/out
  (fvec (* double)) ; m out
  (fjac (* double)) ; ldfjac,n out 
  (ldfjac int :copy) ; 1  in
  (tol double :copy) ; 1 in
  (info int :out) ; 1 out
  (ipvt (* int)) ; n out
  (wa (* double)) ; lwa out
  (lwa int :copy)) ; 1 in


(defun norm (ls)
  (sqrt
   (loop for i below (length ls) 
      sum (expt (elt ls i) 2))))

(defun fit-gaussian (&key (x0 2) (y0 2) (a 1) (b 1) (sigma 1))
  (destructuring-bind (h w) (array-dimensions *img*)
   (let* ((m (* h w)) ;; pixels
	  (n 5) ;; variables
	  (x (make-array n :element-type 'double-float
			 :initial-contents
			 (mapcar #'(lambda (x) (* 1d0 x))
				 (list x0 y0 a b sigma))))
	  (ldfjac m)
	  (lwa (+ m (* 5 n)))
	  (wa (make-array lwa :element-type 'double-float
			  :initial-element 0d0))
	  (fvec (make-array m :element-type 'double-float
			    :initial-element 0d0))
	  (fjac (make-array (list n ldfjac) ;; pixels are fast index
			    ;; in fortran this is m x n
			    ;; F'_ij = \partial f_i/\partial x_j
			    ;; 0<i<m, 0<j<n
			    ;; columns (along j) store the gradient
			    :element-type 'double-float
			    :initial-element 0d0))
	  (ipvt (make-array n
			    :element-type '(signed-byte 32)
			    :initial-element 0))
	  (tol .001d0 #+nil (sqrt double-float-epsilon)))
     (sb-sys:with-pinned-objects (x wa fvec fjac ipvt)
       (labels ((a (ar)
		  (sb-sys:vector-sap 
		   (sb-ext:array-storage-vector
		    ar))))
	 (lmder1_ (alien-sap fcn2) m n (a x) (a fvec) (a fjac) ldfjac tol
		  (a ipvt) (a wa) lwa)))
     
     (let ((fnorm (norm fvec))
	   (fjnorm (make-array n :element-type 'double-float))
	   (eps (* 9 .05)))
       
       ;; |x_i - x*_i| <= s_i    with x_j = x*_j for j/=i
       ;; implies that:
       ;; ||F(x)|| <= (1+eps) ||F(x*)||
       ;; ||F'(x*) e_i|| = ||J e_i||
       ;; ||J e_p(j)|| = ||R e_j||
      (dotimes (j n)
	 (let ((l (- (aref ipvt j) 1)))
	   ;; fjac contains upper triangular matrix R
	   (setf (aref fjnorm l) 
		 (norm (loop for i upto j
			  collect (aref fjac j i))))))
      (values x
	      (loop for i below n collect
				    (* (sqrt eps)
				       (/ fnorm
					  (aref fjnorm i))))
	      fnorm)))))


;; p.18 the i-th row of the jacobian is the gradient of the i-th residual

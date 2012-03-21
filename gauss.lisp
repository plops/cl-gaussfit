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
    (+ (random *rand*) b (* a e))))


(defun fill-img ()
 (defparameter *img*
   (let ((a (make-array '(5 5) :element-type 'double-float)))
     (destructuring-bind (h w) (array-dimensions a)
       (dotimes (j h)
	 (dotimes (i w)
	   (setf (aref a j i) (g i j 2.3d0 2.5d0 .7d0 .005d0 2d0)))))
     a)))

;; load the data and swap endian
(defparameter *imgs*
  (with-open-file (s "example-movie_64x64x30000.raw"
		     :direction :input
		     :element-type '(unsigned-byte 16))
    (let* ((z 30000)
	   (x 64)
	   (y 64)
	   (n (* z x y)) ;; the last 1000 images look wrong
	   (a (make-array n :element-type '(unsigned-byte 16)))
	   (b (make-array (list z y x) :element-type '(unsigned-byte 16)))
	   (b1 (sb-ext:array-storage-vector b)))
      (read-sequence a s)
      (dotimes (i n)
	(setf (aref b1 i) (+ (ldb (byte 8 8) (aref a i))
			     (* 256 (ldb (byte 8 0) (aref a i))))))
      b)))

(declaim (optimize (speed 3) (safety 1) (debug 1)))

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

(defun generate-fits-header (&key (bits 16) (w 484) (h 400))
 (let* ((dat `((SIMPLE T)
	       (BITPIX ,bits)
	       (NAXIS 2)
	       (NAXIS1 ,w)
	       (NAXIS2 ,h))))
   (with-output-to-string (s)
     (labels ((ensure-80-chars (str)
		(format s "~80A" str)))
       (dolist (e dat)
	 (ensure-80-chars ;; the imagej writer puts = at 9th position and value starts at 11
	  (format nil "~8A= ~A" (first e) (second e))))
       (ensure-80-chars "END"))
     (format nil "~2880A" s))))

(defun write-fits (filename img)
  (destructuring-bind (h w) (array-dimensions img)
    (let ((head (generate-fits-header :w w :h h)))
      (with-open-file (s filename
			 :direction :output
			 :if-exists :supersede
			 :if-does-not-exist :create)
	(write-sequence head s))
      (let* ((i1 (make-array (array-total-size img)
			     :element-type (array-element-type img)
			     :displaced-to img))
	     (n (length i1))
	     (h1 (make-array n :element-type '(unsigned-byte 16))))
	(dotimes (i n) ;; swap endian
	  (setf (aref h1 i) (+ (ldb (byte 8 8) (aref i1 i))
			       (* 256 (ldb (byte 8 0) (aref i1 i))))))
	(with-open-file (s filename
			   :element-type '(unsigned-byte 16)
			   :direction :output
			   :if-exists :append)
	  (write-sequence h1 s)))))
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
  (let* ((i1 (array-storage-vector img))
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

(defun uniq (ls)
  (let ((old nil))
   (loop for e in ls
      unless (and old
		  (= old e))
      collect (progn (setf old e) e))))

;; select 3 preceding and 3 following images for each event
#+nil
(defparameter *imgs-more-events*
 (let ((r ())
       (d 2))
   (dolist (e *imgs-with-event*)
     (loop for i from (- e d) upto (+ e d) do
	  (when (<= 0 i 28999)
	    (push i r))))
   (uniq (sort r #'<))))

;; show avgs as well
#+nil
(defparameter *im*
 (loop for e in *imgs-more-events* collect
      (list e (elt *imgs-avg* e))))

#+nil
(length *im*)

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
			:element-type 'double-float))
	 (c1 (make-displaced-array c)))
    (dotimes (i (length a1))
      (setf (aref c1 i) (* 1d0 (funcall op (aref a1 i) (aref b1 i)))))
    c))

#+nil
(defparameter *subtr*
  (loop for i below 100 collect
   (img-op #'-
	   (extract-frame *imgs* i)
	   (extract-frame *imgs* (+ i 30)))))

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

(defun uniq (ls)
  (let ((old nil))
   (loop for e in ls
      unless (and old
		  (= old e))
      collect (progn (setf old e) e))))

;; select 3 preceding and 3 following images for each event
#+nil
(defparameter *imgs-more-events*
 (let ((r ())
       (d 2))
   (dolist (e *imgs-with-event*)
     (loop for i from (- e d) upto (+ e d) do
	  (when (<= 0 i 28999)
	    (push i r))))
   (uniq (sort r #'<))))

;; show avgs as well
#+nil
(defparameter *im*
 (loop for e in *imgs-more-events* collect
      (list e (elt *imgs-avg* e))))

#+nil
(length *im*)

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
			:element-type 'double-float))
	 (c1 (make-displaced-array c)))
    (dotimes (i (length a1))
      (setf (aref c1 i) (* 1d0 (funcall op (aref a1 i) (aref b1 i)))))
    c))

#+nil
(defparameter *subtr*
  (loop for i below 100 collect
   (img-op #'-
	   (extract-frame *imgs* i)
	   (extract-frame *imgs* (+ i 30)))))

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

;; unfortunately the above didn't work. so i probably have to convolve
;; the images first

;; ~/src/mm/3rdpartypublic/sourceext/ij_source/ij/plugin/filter/GaussianBlur.java

(defun convolve-line (in out kern read-from read-to 
		      write-from write-to point0 point-inc)
  (let ((kradius (length kern)))
    ))

(defun make-gaussian-kernel (&key (sigma .6d0) (accuracy 1d-3))
  "sigma .. radius of decay to 1/e in pixels
accuracy .. 1e-3 to 1e-4 for 16bit images. smaller is better"
  (let* ((n (1+ (ceiling (* sigma (sqrt (* -2 (log accuracy)))))))
	 (s (/ 1d0 (* sigma sigma)))
	 (a (make-array (list 2 n) :element-type 'double-float)))
    (dotimes (i n)
      (setf (aref a 0 i) (exp (* -.5 i i s))))
    (when (< 3 n) ;; edge correction to avoid cutoff at finite value
      ;; second order polynomial is zero at first out-of-kernel pixel
      (let* ((sqrt-slope double-float-positive-infinity)
	     (r (loop named slope for r from (1- n) downto (floor n 2) do
		     (let ((v (/ (sqrt (aref a 0 r))
				 (- n r))))
		       (if (< v sqrt-slope)
			   (setf sqrt-slope v)
			   (return-from slope r))))))
	(loop for r1 from (+ r 2) below n do
	     (setf (aref a 0 r1) (expt (* sqrt-slope (- n r1)) 2)))))
    (let ((s (/ 1d0 (+ (aref a 0 0)
		       (* 2 (loop for i from 1 below n sum (aref a 0 i))))))
	  (rsum (+ .5 (* .5 s (aref a 0 0)))))
      (format t "rsum = ~a~%" (list rsum (/ s) (aref a 0 0)))
      (loop for i below n do
	   (let ((v (* s (aref a 0 i))))
	    (setf (aref a 0 i) v)
	    (decf rsum v)
	    (setf (aref a 1 i) rsum))))
    a))
#+nil
(make-gaussian-kernel)
#+nil
(format t "~a~%"
 (make-gaussian-kernel :sigma 5.4 :accuracy 1d-2))

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

(defun uniq (ls)
  (let ((old nil))
   (loop for e in ls
      unless (and old
		  (= old e))
      collect (progn (setf old e) e))))

;; select 3 preceding and 3 following images for each event
#+nil
(defparameter *imgs-more-events*
 (let ((r ())
       (d 2))
   (dolist (e *imgs-with-event*)
     (loop for i from (- e d) upto (+ e d) do
	  (when (<= 0 i 28999)
	    (push i r))))
   (uniq (sort r #'<))))

;; show avgs as well
#+nil
(defparameter *im*
 (loop for e in *imgs-more-events* collect
      (list e (elt *imgs-avg* e))))

#+nil
(length *im*)

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
			:element-type 'double-float))
	 (c1 (make-displaced-array c)))
    (dotimes (i (length a1))
      (setf (aref c1 i) (* 1d0 (funcall op (aref a1 i) (aref b1 i)))))
    c))

#+nil
(time
 (defparameter *subtr*
   (loop for i below 28000 collect
	(img-op #'-
		(extract-frame *imgs* i)
		(extract-frame *imgs* (+ i 300))))))


(defun calc-hist (img &key (n 30) minv maxv (append nil))
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
    (dotimes (i (length i1))
      (incf (aref hist (floor (* n (- (aref i1 i) mi))
			      (- ma mi)))))
    (values hist g)))

#+nil
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

(defun ub16->double-2 (img)
  (let* ((a (make-array (array-dimensions img)
			:element-type 'double-float))
	 (a1 (array-storage-vector a))
	 (n (length a1))
	 (i1 (array-displacement img)))
    (dotimes (i n)
      (setf (aref a1 i) (* .001d0 (aref i1 i))))
    a))

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

#+nil
(list
 (print-with-error 10242.1213 123.02)
 (print-with-error 10242.1213 1.02))

#+nil
(time
 (let ((oldx 0)
       (oldy 0))
  (loop for k from 673 upto 712 do
       (defparameter *img* (ub16->double-2 (extract-frame *imgs* k)))
       (multiple-value-bind (y0 x0) (locate-max-in-frame *imgs* k)
	 (multiple-value-bind (theta sigmas fnorm)
	     ;; use old fit position if it isn't too far off from the maximum
	     (if (< (sqrt (+ (expt (- oldx x0) 2) (expt (- oldy y0) 2))) 1)
		 (fit-gaussian :x0 oldx :y0 oldy :a 1 :b .5 :sigma .85)
		 (fit-gaussian :x0 x0 :y0 y0 :a 1 :b .5 :sigma .85))
	   (destructuring-bind (x y a b s) (loop for e across theta collect e)
	     (destructuring-bind (sx sy sa sb ss) sigmas
	       (format t "~a~%" (list k
				      (print-with-error x sx)
				      (print-with-error y sy)
				      (print-with-error a sa)
				      (print-with-error b sb)
				      (print-with-error s ss)
				      (format nil "~5,2f" fnorm))))
	     (setf oldx x
		   oldy y)))))))


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
	  (tol .05d0 #+nil (sqrt double-float-epsilon)))
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
#+NIL
(progn
  (setf *rand* .001)
  (fill-img)
  (time
   (dotimes (i 30)
     (run2))))
#+nil
(progn
  (format t "~{~7s~}~%" '(rand x0 y0 a b s sx sy sa sb ss))
  (let* ((nj 60)
	 (n 20)
	 (h (make-array (list nj n 11) ;; fast: rand x y a b s sx sy sa sb ss 
			:element-type 'double-float)))
    (dotimes (j nj)
      (dotimes (i n)
	(setf *rand* (+ .001 (* j (/ .1d0 12))))
	(fill-img)
	(run2)
	(loop for k below (length *res*) do
	     (setf (aref h j i k) (elt *res* k)))))
    (defparameter *scan4* h)))
#+nil
(with-open-file (s "/dev/shm/o.gp" :direction :output :if-does-not-exist :create
		   :if-exists :supersede)
  (format s "#set yrange [0:1];
 plot ")
  (let* ((j (+ 2 (* 2 4))) 
	 (n (+ j 1)))
    (loop for i from j upto n do
	 (format s "\"/dev/shm/o.dat\" u 1:~d w l~c" i (if (< i n) #\, #\;))))
  (format s "pause -1"))
#+nil
(with-open-file (s "/dev/shm/o.dat" :direction :output :if-does-not-exist :create
		   :if-exists :supersede)
 (format 
  s "~{~{~a ~}~%~}~%"
  (let ((a *scan4*))
    (destructuring-bind (nj n ne) (array-dimensions a)
      (loop for j below nj collect
	   (append (list (format nil "~3f" (aref a j 0 0)))
		   (let ((real-mean '(2.3d0 2.5d0 .7d0 .005d0 2d0)))
		    (flet ((stat (col)
			     (let* ((avg (loop for i below n sum
					      (* (/ n) (aref a j i col))))
				    (stddev (sqrt (loop for i below n
						     sum
						     (/ (expt (- (aref a j i col)
								 (elt real-mean (1- col)))
							      2)
							n)))))
			       (format nil "~5f" stddev)))
			   (st (col)
			     (let* ((avg (loop for i below n sum
					      (* (/ n) (aref a j i col))))
				    (stddev (sqrt (loop for i below n
						     sum
						     (/ (expt (- (aref a j i col)
								 avg)
							      2)
							n)))))
			       (format nil "~5f" (* 3 avg)))))
					; x y a b s
		      (list (stat 1) (st 6) (stat 2) (st 7)
			    (stat 3) (st 8) (stat 4) (st 9) 
			    (stat 5) (st 10))))))))))

;; p.18 the i-th row of the jacobian is the gradient of the i-th residual

#+nil
(let ((jac *fjac2*))
  (terpri)
 (destructuring-bind (h w) (array-dimensions jac)
   (dotimes (j h)
     (dotimes (i w)
       (format t "~5,3f " (aref jac j i)))
     (terpri)))) 
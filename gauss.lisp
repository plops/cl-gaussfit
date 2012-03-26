#+nil
(progn
  (loop for i in '("utils" "gauss-fit" "gauss-blur" "fits-file" "statistics")
       do
       (load (format nil "~a.lisp" i))))

(defpackage :gauss
  (:use :cl :utils :gauss-fit 
	:gauss-blur :fits-file
	:statistics
	:sb-ext))
(in-package :gauss)

(declaim (optimize (speed 0) (safety 3) (debug 3)))

(defvar *rand* .9)

(defun g (x y x0 y0 a b s)
  (let* ((dx (- x x0))
	 (dy (- y y0))
	 (s2 (/ (* s s)))
	 (arg (- (* s2 (+ (* dx dx) (* dy dy)))))
	 (e (exp arg)))
    (+ b (* a e))))

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

(write-fits "/dev/shm/o.fits"
	    (ub16->single-2
			  (extract-frame *imgs* 1))
	    #+nil
	    (img-op #'-
	     (blur-float (ub16->single-2
			  (extract-frame *imgs* 1))
			 1.3 1.3 1e-1)

	     (blur-float (ub16->single-2
			  (extract-frame *imgs* 1))
			 1.5 1.5 1e-1)))

#+nil
(time
 (let ((start 0)
       (n 100)
       (damp (make-array (list 64 64) :element-type 'single-float)))
   (loop for j from 10 below (- 64 10) do
     (loop for i from 10 below (- 64 10) do
	  (setf (aref damp j i) 1f0)))
   (blur-float damp 5f0 5f0 .02)
   (defparameter *blur*
     (loop for k from start below (+ start n) collect
	  (let* ((s 1.3)
		 (s2 1.5)
		 (im (ub16->single-2
		      (extract-frame *imgs* k)))
		 (g1 (blur-float im s s 1e-5))
		 (g2 (blur-float (copy-img im) s2 s2 1e-5))
		 (dog (img-op #'* damp (img-op #'- g1 g2)
		       )))
	    (when (= 0 (mod k 100))
	      (format t "~a~%" (list 'at k 'until (+ start n))))
	    dog
	    #+nil (mark-points dog (find-local-maxima dog)))))
   (write-fits "/dev/shm/o.fits" (img-list->stack *blur*))))


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
    #+nil (write-fits "/dev/shm/blur-ma.fits" (img-list->stack *blur-ma*))))

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
	 (declare (ignore val))
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
					     (img-apply-mask (copy-img e) m
						       :background mean))))
     #+nil (write-fits "/dev/shm/blur-10x10masked.fits" 
		 (img-list->stack *blur-10x10masked*)))))

#+nil
(defparameter *blur-10x10masked* nil)

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
    #+nil (write-fits "/dev/shm/blur-ma-2.fits" (img-list->stack *blur-ma-2*))))

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
			(calc-hist e :n n 
				   :minv (* (if (< mi 0) 1.1 .9) mi)
				   :maxv (* (if (< ma 0) 1.1 .9) ma)
				   :append hist :mask m))
		      m))))
   
   (multiple-value-bind (mean stddev)
       (get-statistics ims :mask-list masks)
     (let ((masked-images 
	    (loop for e in ims and m in masks collect
		 (img-apply-mask (ub16->single-2 e) m
			   :background mean))))
       (with-open-file (s "/dev/shm/raw.dat" :direction :output
			  :if-does-not-exist :create :if-exists :supersede)
	 (format s "~{~{~f ~f~%~}~} # mean = ~f stddev = ~f~%"
		 (loop for i across hist and g in qq
		    collect (list g 
				  (if (= 0 i) 0 i)))
		 mean stddev))
       (defparameter *raw-10x10masked* masked-images))
     #+nil (write-fits "/dev/shm/raw-marked.fits"
		 (img-list->stack 
		  (loop for e in ims and points in *blur-big-ma-2* collect
		       (mark-points (copy-img e) points
				    :value (+ mean (* -5 stddev))))))
     #+nil (write-fits "/dev/shm/raw-10x10masked.fits" 
		 (img-list->stack *raw-10x10masked*)))))
#+nil
(defparameter *raw-10x10masked* masked-images)
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
(dotimes (i 1000)
  (sleep 100)
  (with-open-file (s "/home/martin/0316/cl-gaussfit/store3.lisp"
		     :direction :output
		     :if-exists :append
		     :if-does-not-exist :create)
   (write *all-fits* :stream s)
   nil))

#+nil
(with-open-file (s "/home/martin/0316/cl-gaussfit/store-blur-big-ma-2.lisp"
		     :direction :output
		     :if-exists :append
		     :if-does-not-exist :create)
   (write *blur-big-ma-2* :stream s)
   nil)

#+nil
(time
 (defparameter *loaded*
   (with-open-file (s "/home/martin/0316/cl-gaussfit/store.lisp")
     (read s))))
#+nil
(setf *all-fits* *loaded*)

#+nil
(time 
 (progn ;; run the fit on a few images (100 takes 160 seconds)
   (defparameter *all-fits* nil)
   (loop for k from 0 below (length *raw*) do
	(progn
	  (progn ;; extract areas around peaks in an image
	    (defparameter *current-image* k)
	    (format t "~d~%" k)
	    (defparameter *raw-blobs*
	      (let ((im (elt *raw* *current-image*)))
		(loop for e in (elt *blur-big-ma-2* *current-image*) collect
		     (destructuring-bind (j i val) e
		       (extract im j i :n 5)))))
	  #+nil  (write-fits "/dev/shm/raw-blobs.fits"
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
	    #+nil(loop for num from 0 and (fnorm val x+err) in (remove-if #'null *fits*) do
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

#+nil
(progn ;; create high res image
  (let* ((nn (length *all-fits*))
	 (step (1- nn))
	 (kk 0)
	(ims
	 (progn ;loop for ky from 0 below 6 collect
	  (loop for kk from 0 below (- nn step) by step collect
	       (destructuring-bind (hh ww) (array-dimensions (first *raw*))
		 (let* ((sc 6)
			(h (* sc hh))
			(w (* sc ww))
			(ar (make-array (list h w) :element-type 'single-float)))
		   (loop for e in (subseq *all-fits* kk (+ kk step)) and pos in 
			(subseq *blur-big-ma-2* kk (+ kk step)) do
			(loop for f in e and p in pos do
			     (when f
			       (destructuring-bind (fnorm val x+err) f
				 (destructuring-bind ((x dx) (y dy) 
						      (a da) (b db) (s ds)) x+err
				   (destructuring-bind (j i val) p
				     (when 
					 (and dx
					     #+nil  (< .7 s))
				       (incf 
					(aref ar 
					      (min (1- h) 
						   (max 0 
							(round (* sc (+ j y -2))))) 
					      (min (1- w) 
						   (max 0 
							(round (* sc (+ i x -2))))))))))))))
		   ar))))))
    (write-fits "/dev/shm/high.fits" (img-list->stack ims))))
#+nil
(progn ;; look at the data
  (with-open-file (ss "/dev/shm/dx-center.dat"
		     :direction :output
		     :if-exists :supersede
		     :if-does-not-exist :create)
   (loop for e in *all-fits* and pos in 
	*blur-big-ma-2* do
	(loop for f in e and p in pos do
	     (when f
	       (destructuring-bind (fnorm val x+err) f
		 (destructuring-bind ((x dx) (y dy) 
				      (a da) (b db) (s ds)) x+err
		   (destructuring-bind (j i val) p
		     (when 
			 (and dx)
		       (let ((xx (sqrt (+ (expt (- x 2) 2)
					  (expt (- y 2) 2))) )
			     (dd (sqrt (+ (expt dx 2)
					  (expt dy 2)))))
			 (when a (and (< 0 x 4) (< 0 y 4) (< fnorm .8)
				    (< 0 a 2) (< .4 b .8) (< dd 2) (< .7 s 3.5))
			   (format ss "~5,6f ~5,6f~%" 
				  (+ x i) (+ y j)))))))))))))

;; ensure sigma > .8 so that those peaks, that were cut on the border of the
;; window are rejected
;; many points have dx < .2

(defun mark-points (img ls &key (value .1))
  (let ((v (coerce value (array-element-type img))))
   (dolist (e ls)
     (destructuring-bind (j i amp) e
       (declare (ignore amp))
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


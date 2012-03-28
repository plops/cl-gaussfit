#+nil
(progn
  (loop for i in '("utils" "gauss-fit" "gauss-blur"
		   "fits-file" "statistics" "kd90")
       do
       (load (format nil "~a.lisp" i))))


(defpackage :gauss
  (:use :cl :utils :gauss-fit 
	:gauss-blur :fits-file
	:statistics :kdtree
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
(time
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
       b))))

#+nil
(time
 (let ((start 20000)
       (n 10000))
   (defparameter *blur*
     (loop for k from start below (+ start n) collect
	  (let* ((s 1.2)
		 (s2 1.3)
		 (im (ub16->single-2
		      (extract-frame *imgs* k)))
		 (g1 (blur-float im s s 1e-3))
		 (g2 (blur-float (copy-img im) s2 s2 1e-3))
		 (dog (img-op #'- g1 g2)))
	    (when (= 0 (mod k 100))
	      (format t "~a~%" 
		      (list 'at k 'until (+ start n))))
	    dog
	    #+nil (mark-points dog (find-local-maxima dog)))))
   (write-fits "/dev/shm/o.fits" (img-list->stack *blur*))
   (progn ;; store subset of the raw data
     (destructuring-bind (z h w) (array-dimensions *imgs*)
       (let* ((part (make-array (list n h w)
				:element-type (array-element-type *imgs*)
				:displaced-to 
				(subseq (make-displaced-array *imgs*)
					0 (* n w h)))))
	 (write-fits "/dev/shm/raw.fits" part)
	 (defparameter *imgs-part* part))))))




#+nil
(time
 (progn ;; find histogram and statistics of difference of gaussian images
   (defparameter *dog-mean* 0)
   (defparameter *dog-stddev* 0)
   (multiple-value-setq (*dog-mean*
			 *dog-stddev*)
     (get-statistics *blur*))))

#+nil
(time 
 (progn ;; histogram of difference of gaussian filtered images
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
      (with-open-file (s "/dev/shm/dog-all.dat" :direction :output
			 :if-exists :supersede
			 :if-does-not-exist :create)
	(format s "~{~{~f ~f~%~}~} # mean = ~f stddev = ~f~%"
		(loop for i across hist and g in q 
		   collect (list g 
				 (if (= 0 i) 0 i)))
		*dog-mean*
		*dog-stddev*))))))


#+nil
(time
  (progn ;; remove maxima that are not brighter than .5 stddev
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
	     (mark-points c big-ma :value (- *dog-mean* (* 3 *dog-stddev*))))))
    (setf *blur-big-ma* (reverse *blur-big-ma*))
    (write-fits "/dev/shm/blur-ma.fits" (img-list->stack *blur-ma*))))

#+nil
(defparameter *blur-ma* nil)

#+nil
(time
 (progn ;; mask 10x10 area around bright maxima in dog images
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
		(let ((m (draw-points-into-mask im points :mask-border 0)))
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
     (write-fits "/dev/shm/blur-10x10masked.fits" 
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
		       (draw-points-into-mask (first *blur*) 
				  (elt *blur-big-ma* i) :mask-border 0))))
	     (multiple-value-setq (*dog-mean-2*
				   *dog-stddev-2*) 
	       (get-statistics *blur* :mask-list mask-list))))))

#+nil
(time
  (progn ;; select maxima again, but use better estimate of background stddev
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
	     (mark-points c big-ma :value (- *dog-mean-2*
					     (* 15 *dog-stddev-2*))))))
    (setf *blur-big-ma-2* (reverse *blur-big-ma-2*))
    (write-fits "/dev/shm/blur-ma-2.fits" (img-list->stack *blur-ma-2*))))

#+nil
(defparameter *blur-ma* nil)
#+nil
(defparameter *blur-ma-2* nil)

#+nil
(sb-ext:gc :full t)

#+nil
(defun centroid (k j i &optional (n 2))
  (let* ((x 0)
	 (y 0)
	 (count 0)
	 (bg 490))
    (loop for jj from (- n) upto n do
	 (loop for ii from (- n) upto n do
	      (let ((g (- (aref *imgs-part* 
				k (+ jj j) (+ ii i))
			  bg)))
	       (incf count g)
	       (incf x (* ii g))
	       (incf y (* jj g)))))
    (if (< 0 count)
     (let ((s (/ 1d0 count)))
       (values (* s y) (* s x)))
     (values 0d0 0d0))))

#+nil
(progn ;; store data for later use
  (with-open-file (s "/home/martin/0316/cl-gaussfit/store-centro4c.lisp"
		     :direction :output
		     :if-exists :supersede
		     :if-does-not-exist :create
		     )
    (write `(progn (defparameter *1dog-mean-2* ,*dog-mean-2*)
		   (defparameter *1dog-stddev-2* ,*dog-stddev-2*)
		   (defparameter *1blur-big-ma-2* ',*blur-big-ma-2*)
		   (defparameter *1all-fits* ',*all-fits*)) :stream s))
  nil)

#+nil
(time
 (load "/home/martin/0316/cl-gaussfit/sicher/store-centro4.lisp"))

#+nil
(progn
  (defparameter *dog-mean-2* *1dog-mean-2*)
  (defparameter *dog-stddev-2* *1dog-stddev-2*)
  (defparameter *blur-big-ma-2* *1blur-big-ma-2*)
  (defparameter *all-fits* *1all-fits*))
#+nil
(defparameter *all-fits* *178fits* )


#+nil 
(progn ;; run the fit on a few images (100 takes 160 seconds)
  (defparameter *all-fits* nil)
  (let ((n 5))
    (destructuring-bind (z h w) (array-dimensions *imgs-part*)
      (loop for k from 0 below z do
	   (format t "~a~%" k)
	   (progn
	     (progn ;; run gauss fits on the extracted images
	       (defparameter
		   *fits*
		 (loop for point in (elt *blur-big-ma-2* k) collect
		      (destructuring-bind (j i val) point
			(when (and (< 2 i (- w 2 1))
				   (< 2 j (- h 2 1))
					;(< 100 val)
				   )
			  (multiple-value-bind (x err fnorm)
			      (fit-gaussian 
			       :stack *imgs*
			       :window-w n
			       :center-x i :center-y j
			       :center-slice (+ 20000 k)
			       :x0 (floor n 2) 
			       :y0 (floor n 2)
			       :a 1f0 :b .49 :sigma .8)
			    (list
			     fnorm val (list i j)
			     (loop for e across x and r in err collect
				  (list e r))))))))
	       (progn
		 (with-open-file (str "/dev/shm/centro.dat" :direction :output :if-exists :append
				      :if-does-not-exist :create)
		   (loop for num from 0 and (fnorm val (i j) x+err) in
			(remove-if #'null *fits*) do
			(destructuring-bind ((x dx) (y dy) (a da) (b db) (s ds)) x+err
			  (multiple-value-bind (cy cx) (centroid k j i)
			    (format str "~3d ~3d ~5,2f ~5,2f ~3d ~5,2f ~5,2f ~5,2f ~4,2f ~4,0f ~{~{~7,2f ~3,2f~}~}~%"
				    num  
				    i (+ i cx) (+ i x -2) 
				    j (+ j cy) (+ j y -2)
				    (sqrt (+ (expt (- cx (- x 2)) 2)
					     (expt (- cy (- y 2)) 2)))
				    fnorm val x+err))))))
	       
	       (setf *all-fits* (append *all-fits* (list *fits*))))
	     
	     )))))
(declaim (optimize (debug 3) (safety 3)))




#+nil
(time
 (progn ;; generate images of the fitted functions
   (let ((res nil))
     (loop for if below (length *all-fits*) by 1 collect
	  (let ((a (make-array (cdr (array-dimensions *imgs*))
			       :element-type 'single-float
			       :initial-element .49))
		(fitcalc (loop for e in (elt *all-fits* if) collect
			      (when e
				(destructuring-bind (fnorm val (ii jj) x+err) e
				  (destructuring-bind ((x dx) (y dy) (a da)
						       (b db) (s ds)) x+err
				    (let* ((n 5) 
					   (ar (make-array (list n n)
							   :element-type 'single-float)))
				      (dotimes (j n)
					(dotimes (i n)
					  (setf (aref ar j i) (coerce (g i j x y a b s)
								      'single-float))))
				      ar)))))))
	    (loop for f in fitcalc and pos in (elt *blur-big-ma-2* if) do
		 (when f
		   (destructuring-bind (j i val) pos
		     (insert-rectangle a f j i))))
	    (let ((orig (img-mul (ub16->single-2 (extract-frame *imgs* 
								(+ if 20000)))
				 .001)))
	      (setf res (append res 
				(list ;a
				      ;orig
				       (img-op #'- orig a)))))))
     (write-fits "/dev/shm/fit-calc.fits" (img-list->stack res)))))


#+nil
(progn ;; create high res image
  (let* ((nn (length *all-fits*))
	 (step (1- nn))
	(ims
	 (progn ;loop for ky from 0 below 6 collect
	  (loop for kk from 0 below (- nn step) by step collect
	       (destructuring-bind (zz hh ww) (array-dimensions *imgs*)
		 (let* ((sc 3)
			(h (* sc hh))
			(w (* sc ww))
			(ar (make-array (list h w) :element-type 'single-float
					:initial-element 0f0)))
		   (loop for e in (subseq *all-fits* kk (+ kk step)) and pos in 
			(subseq *blur-big-ma-2* kk (+ kk step)) do
			(loop for f in e and p in pos do
			     (when f
			       (destructuring-bind (fnorm val (ii jj) x+err) f
				 (destructuring-bind ((x dx) (y dy) 
						      (a da) (b db) (s ds)) x+err
				   (destructuring-bind (j i val) p
				     (when 
					 t
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
(with-open-file (file "/dev/shm/o.dat" :direction :output 
		   :if-exists :supersede
		   :if-does-not-exist :create)
 (loop for e in *all-fits* do
      (loop for f in e do
	   (when f
	     (destructuring-bind (fnorm val (i j) x+err) f
	       (destructuring-bind ((x dx) (y dy) 
				    (a da) (b db) (s ds)) x+err
		 (format file "~12f ~12f~%" x s)))))))

;; ensure sigma > .8 so that those peaks, that were cut on the border of the
;; window are rejected
;; many points have dx < .2
#+nil
(time
 (progn ;; build a kdtree of the points
   (let* ((points nil))
     (loop for e in *all-fits* do
	  (loop for f in e do
	       (when f
		 (destructuring-bind (fnorm val (i j) x+err) f
		   (destructuring-bind ((x dx) (y dy) 
					(a da) (b db) (s ds)) x+err
		     (when (< .6 s 1)
		       (push (make-array 2 :element-type 'single-float
					 :initial-contents
					 (mapcar #'(lambda (x) (coerce x 'single-float))
						 (list (+ i x -2) (+ j y -2))))
			     points)))))))
     (let ((point-a (make-array (length points)
				:element-type 'vec
				:initial-contents points)))
       (defparameter *point-a* point-a)
       (time (defparameter *tree* (build-new-tree point-a)))
       (length point-a)))))

#+nil
(time
 (progn ;; 43s ;; find nearest neighbour distances
   (let* ((n (length *point-a*))
	  (dists (make-array n :element-type 'single-float)))
     (dotimes (i n)
       (multiple-value-bind (p d)
	   (nearest-neighbour i *tree*)
	 (setf (aref dists i) d)))
     (defparameter *dists* dists))))
 
#+nil
(progn ;; create histogram of nearest neighbour distances 
  (let* ((n 12800)
	 (hist (make-array n :element-type '(unsigned-byte 64)))
	 (ma 2.5 ;(reduce #'max *dists*)
	   ))
    (loop for d across *dists* do
	 (incf (aref hist (min (1- n) (floor (* n d (/ (* 1. ma))))))))
    ;; figure out x axis of histogram
    ;; i = n d /(1. ma)
    ;; d = 1.1 ma i / n
    (defparameter *dists-hist* hist)
    (with-open-file (s "/dev/shm/dist-hists.dat"
		       :direction :output
		       :if-exists :supersede
		       :if-does-not-exist :create)
      (dotimes (i (length hist))
	(format s "~f ~f~%"
		(/ (* 1 ma i) n) 
		(aref hist i))))
    (with-open-file (s "/dev/shm/dist-hists.gp"
		   :direction :output
		   :if-exists :supersede
		   :if-does-not-exist :create)
      (format s "set xlabel \"distance in px\"; set ylabel \"frequency\";plot \"/dev/shm/dist-hists.dat\" u 1:2 w p pt 0; pause -1"))))



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


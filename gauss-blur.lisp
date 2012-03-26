(defpackage :gauss-blur
  (:use #:cl #:utils)
  (:export #:blur-float
	   #:make-gaussian-kernel))

(in-package :gauss-blur)

;; ~/src/mm/3rdpartypublic/sourceext/ij_source/ij/plugin/filter/GaussianBlur.java

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
      (let* ((sqrt-slope sb-ext:single-float-positive-infinity)
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

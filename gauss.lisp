(defun g (x y x0 y0 a b s)
  (let* ((dx (- x x0))
	 (dy (- y y0))
	 (s2 (/ (* s s)))
	 (arg (- (* s2 (+ (* dx dx) (* dy dy)))))
	 (e (exp arg)))
    (+ b (* a e))))

(defparameter *img*
 (let ((a (make-array '(16 16) :element-type 'double-float)))
   (destructuring-bind (h w) (array-dimensions a)
     (dotimes (j h)
       (dotimes (i w)
	 (setf (aref a j i) (g i j 8.3d0 8.5d0 .7d0 .005d0 3d0)))))
   a))

#|
jacobian([b+a*exp(-((x-xx)^2+(y-yy)^2)/sigma^2)-f],[xx,yy,a,b,sigma]);
dx = x-xx
dy = y-yy
s2 = 1/s^2
arg = -(dx^2+dy^2) s2
e = exp(arg)
f(x) = a e + b
f_xx = 2 dx s2 f
f_yy = 2 dy s2 f
f_a  = e
f_b  = 1
f_s  = 2 arg/s f
|#

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
      (format t "bla ~a ~%" (list (deref iflag 0)
				  (loop for i below 5 collect (deref x i))))
      (ecase (deref iflag 0)
	(1 (dotimes (j h)
	     (dotimes (i w)
	       (let* ((dx (- (* 1d0 i) xx))
		      (dy (- (* 1d0 j) yy))
		      (s2 (/ (* s s)))
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
		   ;; f_s  = 2 arg/s f	
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
			  (p (+ i (* w j))))
		     (setf (f p 0) fxx
			   (f p 1) fyy
			   (f p 2) fa
			   (f p 3) fb
			   (f p 4) fs)
		     (let ((h 1d-5))
		      (let ((dxx 
				 (/ (- (g i j (+ xx h) yy a b s) (g i j xx yy a b s)) h))
			       (dyy
				 (/ (- (g i j xx (+ yy h) a b s) (g i j xx yy a b s)) h))
			       (da 
				 (/ (- (g i j xx yy (+ a h) b s) (g i j xx yy a b s)) h))
			       (db 
				 (/ (- (g i j xx yy a (+ b h) s) (g i j xx yy a b s)) h))
			       (ds 
				 (/ (- (g i j xx yy a b (+ s h)) (g i j xx yy a b s)) h)))
			(format t "checkderiv ~1,12f ~{~2,3f ~}~%" 
				(reduce #'(lambda (x y) (+ x (* y y)))
					(list (- fxx dxx)
					      (- fyy dyy)
					      (- fa da)
					      (- fb db) 
					      (- fs ds)))
				(list (/ (- fxx dxx) dxx)
					      (/ (- fyy dyy) dyy)
					      (/ (- fa da) da)
					      (/ (- fb db) db) 
					      (/ (- fs ds) ds)))))
		     )))))))))  
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

(defun run2 ()
  (destructuring-bind (h w) (array-dimensions *img*)
   (let* ((m (* h w))
	  (n 5)
	  (x (make-array n :element-type 'double-float
			 :initial-contents
			 (mapcar #'(lambda (x) (* 1d0 x))
				 '(8 8 1 .01 3.01))))
	  (ldfjac m)
	  (lwa (+ m (* 5 n)))
	  (wa (make-array lwa :element-type 'double-float
			  :initial-element 0d0))
	  (fvec (make-array m :element-type 'double-float
			    :initial-element 0d0))
	  (fjac (make-array (list ldfjac n) 
			    :element-type 'double-float
			    :initial-element 0d0))
	  (ipvt (make-array n
			    :element-type '(signed-byte 64)
			    :initial-element 0))
	  (tol (sqrt double-float-epsilon)))
     (sb-sys:with-pinned-objects (x wa fvec fjac ipvt)
       (labels ((a (ar)
		  (sb-sys:vector-sap 
		   (sb-ext:array-storage-vector
		    ar))))
	 (lmder1_ (alien-sap fcn2) m n (a x) (a fvec) (a fjac) ldfjac tol
		  (a ipvt) (a wa) lwa)))
     #+nil
      (format t "~a ~%" (reduce #'(lambda (x y) (+ x (* y y))) fvec)))))
#+NIL
(run2)
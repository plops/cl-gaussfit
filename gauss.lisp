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

(defun g (x y x0 y0 a b)
  (+ b
   (* a (/ 2 pi)
      (exp (* -.5 (+ (expt (- x x0) 2)
		     (expt (- y y0) 2)))))))


(load-shared-object "/usr/lib/libminpack.so")

(sb-alien::define-alien-callback fcn
    void
    ((m (* int)) ; 1
     (n (* int)) ; 1
     (x (* double)) ; n
     (fvec (* double)) ; m
     (fjac (* double)) ; ldfjac,n
     (ldfjac (* int)) ; 1  
     (iflag (* int))) ; 1 1: fvec, 2: fjac
  (let ((y (make-array 15 :element-type 'double-float
		       :initial-contents
		       (mapcar #'(lambda (x) (* 1d0 x))
			       '(.14 .18 .22 .25 .29 .32 .35 .39
				 .37 .58 .73 .96 1.34 2.1 4.39)))))
    (format t "bla ~a ~%" (list (deref iflag 0)
				(loop for i below 3 collect (deref x i))))
    (ecase (deref iflag 0)
      (1 (dotimes (i 15)
	   (let* ((t1 i)
		  (t2 (- 15 i))
		  (t3 (if (< i 8) t1 t2)))
	     (setf (deref fvec i) 
		   (- (aref y i)
		      (+ (deref x 0)
			 (/ t1 (+ (* (deref x 1) t2)
				  (* (deref x 2) t3)))))))))
      (2 (let ((w (deref ldfjac 0))) 
	   (dotimes (i 15)
	     (let*  ((t1 i)
		     (t2 (- 15 i))
		     (t3 (if (< i 8) t1 t2))
		     (t4 (expt (+ (* (deref x 1) t2)
				  (* (deref x 2) t3))
			       2)))
	       (macrolet ((f (a b)
			    `(deref fjac (+ ,a (* w ,b)))))
		(setf (f i 0) -1d0
		      (f i 1) (/ (* t1 t2) t4)
		      (f i 2) (/ (* t1 t3) t4)))))))))  
  (values))


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


(defun run ()
 (let* ((m 15)
	(n 3)
       
	(x (make-array n :element-type 'double-float
		       :initial-element 1d0))
	(ldfjac 15)
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
   (sb-sys:without-gcing
    (sb-sys:with-pinned-objects (x  wa fvec fjac ipvt)
      (labels ((a (ar)
		 (sb-sys:vector-sap 
		  (sb-ext:array-storage-vector
		   ar))))
	(lmder1_ (alien-sap fcn) m n (a x) (a fvec) (a fjac) ldfjac tol
		 (a ipvt) (a wa) lwa))))))

#+nil
(time (run))
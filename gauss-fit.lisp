(defpackage :gauss-fit
  (:use #:cl)
  (:export #:fit-gaussian
	   #:print-with-error1))

;; ANL8074a.pdf documentation for minpack

;; note that minpack isn't thread safe

;; how to do fitting on the gpu (they consider poisson noise but it
;; doesn't seem to improve precision, their method is very fast but
;; they don't give code)
;; /home/martin/ondrej-papers/papers/Quan2010OSA.pdf

;; fitting precision 2004 ober
;;  /home/martin/ondrej-papers/papers/sdarticle\ \(3\).pdf


;; jacobian([b+a*exp(-((x-xx)^2+(y-yy)^2)/sigma^2)-f],[xx,yy,a,b,sigma]);
;; dx = x-xx
;; dy = y-yy
;; s2 = 1/s^2
;; arg = -(dx^2+dy^2) s2
;; e = exp(arg)
;; ff = a e
;; f(x) = a e + b
;; f_xx = 2 dx s2 ff
;; f_yy = 2 dy s2 ff
;; f_a  = e
;; f_b  = 1
;; f_s  = - 2 arg/s ff

(defvar *imgs* nil) ;; 16bit 3d stack, note that values are divided by
		    ;; 1000 before use

(defvar *current-center* (list 0 0 0)) ;; integer position of the
				       ;; current maximum
(defvar *current-window-size* (list 5 5))

(sb-alien::define-alien-callback fcn2
    void
    (;(m (* int)) ; = 8x8 or similar
     ;(n (* int)) ; = 5
     (x (* double)) ; [xx,yy,a,b,sigma]
     (fvec (* double)) ; m
     (fjac (* double)) ; ldfjac,n
     (ldfjac (* int)) ; 1  ;; leading dimension of fjac
     (iflag (* int))) ; 1 1: fvec, 2: fjac
  (destructuring-bind (pz py px) *current-center*
    (destructuring-bind (h w) *current-window-size*
      (let* ((xx (deref x 0))
	     (yy (deref x 1))
	     (a (deref x 2))
	     (b (deref x 3))
	     (s (deref x 4))) 
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
			 (+ (* -.001 (aref *imgs*
					   pz 
					   (+ py j (- (floor h 2)))
					   (+ px i (- (floor w 2))))) 
			    b (* a e)))))))
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
			     (f p 4) fs)))))))))))  
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

(defun fit-gaussian (&key (stack *imgs*) 
		     (window-w 5) (window-h window-w)
		     (center-x (floor window-w 2)) (center-y (floor window-h 2))
		     (x0 2) (y0 2) (a 1) (b 1) (sigma 1))
  (setf *imgs* stack)
  (setf *current-center* (list center-y center-x))
  (setf *current-window-size* (list window-h window-w))
  (let* ((m (* window-h window-w)) ;; pixels
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
	 (tol .01d0 #+nil (sqrt double-float-epsilon)))
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
		   (unless (< (abs (aref fjnorm i)) double-float-epsilon)
		     (* (sqrt eps)
			(/ fnorm
			   (aref fjnorm i)))))
	      fnorm))))


;; p.18 the i-th row of the jacobian is the gradient of the i-th residual


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


;; ANL8074a.pdf documentation for minpack

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

(defparameter *imgs*
  (with-open-file (s "example-movie_small5.raw"
		     :direction :input
		     :element-type '(unsigned-byte 16))
    (let* ((n (* 30000 5 5))
	   (a (make-array n :element-type '(unsigned-byte 16)))
	   (b (make-array (list 30000 5 5) :element-type '(unsigned-byte 16)))
	   (b1 (sb-ext:array-storage-vector b)))
      (read-sequence a s)
      (dotimes (i n)
	(setf (aref b1 i) (+ (ldb (byte 8 8) (aref a i))
			     (* 256 (ldb (byte 8 0) (aref a i))))))
      b)))

(defparameter *imgs-avg*
  (destructuring-bind (z y x) (array-dimensions *imgs*)
    (loop for k below z collect
	 (let ((sum 0))
	   (dotimes (j y)
	     (dotimes (i x)
	       (incf sum (aref *imgs* k j i))))
	   sum))))

(defparameter *imgs-avg-hist*
 (let* ((ma (1+ (reduce #'max *imgs-avg*)))
	(n 10)
	(hist (make-array n :element-type 'fixnum)))
   (dolist (e *imgs-avg*)
     (incf (aref hist (floor (* n e) ma))))
   hist))

(time
 (defparameter *imgs-with-event*
   (let ((ma (1+ (reduce #'max *imgs-avg*))))
     (loop for k below (length *imgs-avg*)
	when (< .4 (/ (elt *imgs-avg* k)
		      ma))
	collect
	(list k (elt *imgs-avg* k))))))

(time
 (defparameter *imgs-with-event*
   (let ((ma (1+ (reduce #'max *imgs-avg*))))
     (loop for k from 0 and e in *imgs-avg* 
	when (< .4 (/ e	ma))
	collect
	(list k e)))))

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


(defun calc-sigmas (x fvec &optional (epsilon 0.05))
  (destructuring-bind (h w) (array-dimensions *img*)
    (let* ((m (* h w))
	   (n (length x))
	   (res (make-array n :element-type 'double-float
			    :initial-element 0d0))
	   (xx (aref x 0))
	   (yy (aref x 1))
	   (a (aref x 2))
	   (s (aref x 4)))
      (labels ((icf (a v)
		 (incf (aref res a) (* v v))))
       (dotimes (j h)
	 (dotimes (i w)
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
	     (icf 0 fxx)
	     (icf 1 fyy)
	     (icf 2 fa)
	     (icf 3 fb)
	     (icf 4 fs))))
       
       (dotimes (i n)
	 (setf (aref res i) (* (sqrt epsilon)
			       (/ (sqrt (loop for e across fvec sum (* e e)))
				  (sqrt (aref res i))))))
       res))))



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

(defun run2 ()
  (destructuring-bind (h w) (array-dimensions *img*)
   (let* ((m (* h w)) ;; pixels
	  (n 5) ;; variables
	  (x (make-array n :element-type 'double-float
			 :initial-contents
			 (mapcar #'(lambda (x) (* 1d0 x))
				 '(2 2 1 1 1))))
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
	  (tol .0001d0 #+nil (sqrt double-float-epsilon)))
     (sb-sys:with-pinned-objects (x wa fvec fjac ipvt)
       (labels ((a (ar)
		  (sb-sys:vector-sap 
		   (sb-ext:array-storage-vector
		    ar))))
	 (lmder1_ (alien-sap fcn2) m n (a x) (a fvec) (a fjac) ldfjac tol
		  (a ipvt) (a wa) lwa)))
     
     (let ((fnorm (norm fvec))
	   (fjnorm (make-array n :element-type 'double-float))
	   (eps .05))
       
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

      (format t "~{~9s ~{~6,3f ~}~%~} ~%" 
	      (list 		    'rel-sigma (loop for i below n collect
				    (/ (* (sqrt eps)
					  (/ fnorm
					     (aref fjnorm i)))
				       (aref x i)))
		    'x (loop for e across x collect e)
		    'sigma (loop for i below n collect
				    (* (sqrt eps)
				       (/ fnorm
					  (aref fjnorm i))))
		    'sigma-my (loop for e across (calc-sigmas x fvec)
				   collect e)))
      (defparameter *res* (append (list *rand*)
				  (loop for i across x collect i) 
				  (loop for i below n collect
				       (/ (* (sqrt eps)
					     (/ fnorm
						(aref fjnorm i)))
					  (abs (aref x i))))))
     #+nil (format t "~{~8,4f~}~%" *res*)
      (defparameter *fjac* fjac)))))
#+NIL
(progn
  (setf *rand* .001)
  (fill-img)
  (time
   (dotimes (i 30)
     (run2))))

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

(with-open-file (s "/dev/shm/o.gp" :direction :output :if-does-not-exist :create
		   :if-exists :supersede)
  (format s "#set yrange [0:1];
 plot ")
  (let* ((j (+ 2 (* 2 4))) 
	 (n (+ j 1)))
    (loop for i from j upto n do
	 (format s "\"/dev/shm/o.dat\" u 1:~d w l~c" i (if (< i n) #\, #\;))))
  (format s "pause -1"))

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
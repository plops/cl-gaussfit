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

(sb-alien::define-alien-callback bla sb-alien:double ((x sb-alien:double))
  (* 12d0))

(sb-alien:load-shared-object "/home/martin/0316/cl-gaussfit/lib1.so")

(define-alien-routine fun sb-alien:double (fptr (* int)))

(fun (sb-alien:alien-sap bla))

(load-shared-object "/usr/lib/libminpack.so")

(define-alien-routine lmder1_ void
  (fcn (* int)) ; callback
  (m (* int)) ; 1
  (n (* int)) ; 1
  (x (* double)) ; n
  (fvec (* double)) ; m
  (fjac (* double)) ; ldfjac,n 
  (ldfjac (* int)) ; 1
  (tol (* double)) ; 1
  (info (* int)) ; 1
  (ipvt (* int)) ; n
  (wa (* double)) ; lwa
  (lwa (* int))) ; 1

(sb-alien::define-alien-callback fcn
    void
    ((m (* int)) ; 1
     (n (* int)) ; 1
     (x (* double)) ; n
     (fvec (* double)) ; m
     (fjac (* double)) ; ldfjac,n
     (ldfjac (* int)) ; 1
     (iflag (* int))) ; 1
  (values))


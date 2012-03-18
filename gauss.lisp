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

               [                          2           2 ]
               [                - (y - yy)  - (x - xx)  ]
               [                ----------------------- ]
               [                             2          ]
(%o1)  Col 1 = [                        sigma           ]
               [ 2 a (x - xx) %e                        ]
               [ -------------------------------------- ]
               [                      2                 ]
               [                 sigma                  ]
         [                          2           2 ]
         [                - (y - yy)  - (x - xx)  ]
         [                ----------------------- ]
         [                             2          ]
 Col 2 = [                        sigma           ]
         [ 2 a (y - yy) %e                        ]
         [ -------------------------------------- ]
         [                      2                 ]
         [                 sigma                  ]
         [             2           2 ]
         [   - (y - yy)  - (x - xx)  ]
 Col 3 = [   ----------------------- ] Col 4 = [ 1 ]
         [                2          ]
         [           sigma           ]
         [ %e                        ]
         [                                             2           2 ]
         [                                   - (y - yy)  - (x - xx)  ]
         [                                   ----------------------- ]
         [                                                2          ]
 Col 5 = [                  2           2            sigma           ]
         [   2 a (- (y - yy)  - (x - xx) ) %e                        ]
         [ - ------------------------------------------------------- ]
         [                                3                          ]
         [                           sigma                           ]
|#

(defun g (x y x0 y0 a b)
  (+ b
   (* a (/ 2 pi)
      (exp (* -.5 (+ (expt (- x x0) 2)
		     (expt (- y y0) 2)))))))
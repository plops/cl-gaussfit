(defpackage :utils
  (:use #:cl)
  (:export #:make-displaced-array))

(in-package :utils)

(defun make-displaced-array (a)
  (make-array (array-total-size a)
	      :element-type (array-element-type a)
	      :displaced-to a))
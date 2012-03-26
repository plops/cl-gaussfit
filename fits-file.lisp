(defpackage :fits-file
  (:use #:cl #:utils)
  (:export #:write-fits))

(in-package :fits-file)

;; definition of fits file
;; ftp://nssdc.gsfc.nasa.gov/pub/fits
;; 80 byte per header entry
;; from imagej source
;; ~/src/mm/3rdpartypublic/sourceext/ij_source/ij/plugin/FTIS_Reader.java
;; start with SIMPLE = T / comment
;; = separates key value pair
;; / indicates comment
;; END stops
;; BITPIX = 16 short 32 int -32 float, -64 double
;; NAXIS1 = 129 / width, NAXIS2, NAXIS3
;; BSCALE = 1.34 / should be 1.0, then no calibration in imagej
;; BZERO = 32768 / for ushort (expt 2 15)
;; probably 2147483648 = (expt 2 31) for uint
;; count such lines (each 80 bytes long)
;; offset with data at 2880+2880*(((count*80)-1)/2880)
;; origin at bottom left

;; offset in example file #x1680 = 5760 = 2*2880

(defun generate-fits-header (img &key (bits 16))
  (unless bits
    (error "bits isn't set correctly"))
  (let* ((dims (array-dimensions img))
	 (dat `((SIMPLE T)
		(BITPIX ,bits)
		(NAXIS ,(length dims))
		(NAXIS1 ,(first (last dims)))
		(NAXIS2 ,(first (last (butlast dims))))
		,(when (= (length dims) 3)
		       (list 'NAXIS3 (first dims))))))
    (let* ((head (with-output-to-string (s)
		   (labels ((ensure-80-chars (str)
			      (format s "~80A" str)))
		     (dolist (e dat)
		       (when e 
			 (ensure-80-chars ;; the imagej writer puts = at 9th position and value starts at 11
			  (format nil "~8A= ~A" (first e) (second e)))))
		     (ensure-80-chars "END"))
		   s))
	   (pad (make-array (- 2880 (length head))
			    :element-type 'character
			    :initial-element #\Space)))
      (concatenate 'string head pad))))

#+nil
(format t "~a~%"
 (generate-fits-header (first *blur*) :bits 16))

(defun write-fits (filename img)
  (let ((head 
	 (generate-fits-header 
	  img
	  :bits (cond
		  ((equal (array-element-type img) 'single-float) 
		   -32)
		  ((or
		    (equal (array-element-type img) '(unsigned-byte 16))
		    (equal (array-element-type img) '(signed-byte 16)))
		   16)))))
    (with-open-file (s filename
		       :direction :output
		       :if-exists :supersede
		       :if-does-not-exist :create)
      (write-sequence head s))
        (let* ((i1 (array-displacement (copy-img img)))
		    (sap (sb-sys:vector-sap i1))
		    (n (cond
		  ((equal (array-element-type img) 'single-float) 
		   (array-total-size img))
		  ((or
		    (equal (array-element-type img) '(unsigned-byte 16))
		    (equal (array-element-type img) '(signed-byte 16)))
		   (floor (array-total-size img)
			  2))))
	   (h1 (make-array n :element-type '(unsigned-byte 32)))
	   (h2 (make-array n :element-type '(unsigned-byte 32))))
      (dotimes (i n)
	(setf (aref h1 i) (sb-sys:sap-ref-32 sap (* 4 i))))
      (cond ;; swap endian
	((equal (array-element-type img) 'single-float)
	 (dotimes (i n) 
	   (setf (ldb (byte 8 0) (aref h2 i)) (ldb (byte 8 24) (aref h1 i))
		 (ldb (byte 8 8) (aref h2 i)) (ldb (byte 8 16) (aref h1 i))
		 (ldb (byte 8 16) (aref h2 i)) (ldb (byte 8 8) (aref h1 i))
		 (ldb (byte 8 24) (aref h2 i)) (ldb (byte 8 0) (aref h1 i)))))
	((or (equal (array-element-type img) '(unsigned-byte 16))
	     (equal (array-element-type img) '(signed-byte 16)))
	 (dotimes (i n) 
	   (setf (ldb (byte 8 0) (aref h2 i)) (ldb (byte 8 8) (aref h1 i))
		 (ldb (byte 8 8) (aref h2 i)) (ldb (byte 8 0) (aref h1 i))
		 (ldb (byte 8 16) (aref h2 i)) (ldb (byte 8 24) (aref h1 i))
		 (ldb (byte 8 24) (aref h2 i)) (ldb (byte 8 16) (aref h1 i))))))
      (with-open-file (s filename
			 :element-type '(unsigned-byte 32)
			 :direction :output
			 :if-exists :append)
	(write-sequence h2 s))))
  (values))

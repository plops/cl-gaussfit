(require :asdf)
(require :run)
(defpackage :run
  (:use :cl :kielhorn.martin.math.vector))
(in-package :run)

;; jedes dreieck, dass der algorithmus jemals in betracht zieht wird
;; in ein array geschrieben. jedesmal wenn ein dreieck in 3 kleinere
;; geteilt wird, werden links zu den drei toechtern gesetzt. wenn eine
;; seite getauscht wird, dann werden die zwei neuen dreiecke als
;; toechter der beiden alten eingetragen

(defstruct triel
  "VERTICES enthaelt dreiecksecken als indizes in
array-points. DAUGHTERS enthaelt entweder pointer zu drei enthaltenen
kleineren dreiecken, zwei dreiecken die nur ueberlappen -- durch
seitentausch erzeugt und -1 als letztes element oder dreimal -1 --
bedeutet keine tochter. STATUS ist nicht null, wenn das dreieck gerade
lebt."
  (vertices #(0 0 0) :type '(array integer 3))			
  (daughters #(-1 -1 -1) :type '(array integer 3))
  (status 1 :type integer))

(defparameter array-triel (make-array 1000 :element-type 'triel))

(defun get-triangle (i)
  (triel-vertices (aref array-triel i)))

(defun get-triangle-a (i)
  (aref (get-triangle i) 0))
(defun get-triangle-b (i)
  (aref (get-triangle i) 1))
(defun get-triangle-c (i)
  (aref (get-triangle i) 2))


;; wir brauchen zwei schnelle lookups:

;; 1) gegeben ein punkt und eine gegenueberliegende seite, finde den
;; eckpunkt des gegenueberliegenden dreiecks. dafuer hashen wir fuer
;; jedes dreieck (immer CCW) mit ecken A, B, C den eckpunkt A unter
;; der linie: ht-lines[h(B)-h(C)]=A. (und zyklisch vertauscht)

;; um einen punkt auf der anderen seite von BC zu finden, schauen wir
;; einfach in h(C)-h(B) nach

(defparameter ht-lines (make-hash-table))

;; 2) gegeben drei punkte (in beliebiger reihenfolge), wo ist der
;; eintrag im grossen array? wir speichern den index j unter einer hash
;; die aus den drei geXORten punkten besteht:
;; ht-triangles(h(A)^h(B)^h(C))=j

(defparameter ht-triangles (make-hash-table))



;; basic data structure point

(defparameter array-points (make-array 3000 :element-type 'vec3))
(defparameter n-points 0)
(defun append-point (point)
  (setf (aref array-points n-points) point
	n-points (1+ n-points))
  (1- n-points))
(defun get-point (i)
  (aref array-points i))


(defun vertices (ti)
  "mit gebenenen index in array-triel liefere die koordinaten der
  eckpunkte"
  (mapcar #'get-point
	  (get-triangle i)))

;; geometry stuff

(defclass circle  ()
  ((center :reader center :initarg :center :initform (v 0 0))
   (r^2 :reader r^2 :initarg :r^2 :initform 1)))

(defgeneric r (obj))

(defmethod r ((circ circle))
  (sqrt (r^2 circ)))

(defun centroid (ti)
  "Center of mass of a triangle TI."
  (v* 1/3 (reduce 'v+ (vertices ti))))
 
(defun delta (ti)
  "Calculate the edge directions of the triangle (in CCW order)."
  (let ((result '())
	(old-p))
    (dolist (p (vertices ti))
      (when old-p
	(push (v- p old-p) result))
      (setf old-p p))
    (push (v- (first (vertices ti)) old-p) result)
    (nreverse result)))

(defun incircle (ti)
  "Inscribed circle of the triangle (Numerical Recipes 2007 p. 1112)."
  (let* ((points (vertices ti))
	 (dtri (delta ti))
	 (side-lengths (mapcar #'norm dtri))
	 (perimeter (reduce #'+ side-lengths))
	 (center (v* (/ perimeter) 
		     (reduce #'v+ 
			     (mapcar #'(lambda (p d) 
					 (v (* (vec-x p)  (norm d))
					    (* (vec-y p) (norm d))))
				     points dtri))))
	 (r2 (/ (reduce #'* 
			(mapcar #'(lambda (l) (- perimeter l))
				side-lengths))
		perimeter)))
    (make-instance 'circle :center center :r^2 r2)))

(define-condition det0-error (error)
  ((message :accessor message :initarg :message)))

(defun det0-error (message)
  (error 'det0-error :message message))

(defmethod print-object ((obj det0-error) stream)
  (print-unreadable-object (obj stream :type t :identity t)
    (format stream "~a" (message obj))))

(defun circumcircle (ti)
  "Circumscribed circle of a triangle."
  (let* ((points (vertices ti))
	 (ba (v- (second points)
		 (first points)))
	 (ca (v- (third points)
		 (first points)))
	 (ba2 (dot ba ba))
	 (ca2 (dot ca ca))
	 (ndet (- (* (vec-x ba) (vec-y ca))
		  (* (vec-y ba) (vec-x ca))))
	 (ndet/ (if (= 0 ndet)
		    (det0-error "No circle thru colinear points.")
		    (/ .5 ndet)))
	 (dcenter (v* ndet/ (v (- (* ba2 (vec-y ca))
				  (* ca2 (vec-y ba)))
			       (- (* ca2 (vec-x ba))
				  (* ba2 (vec-x ca))))))
	 (r2 (dot dcenter dcenter)))
    (make-instance 'circle :center (v+ (first points) dcenter) :r^2 r2)))

(defun contains (ti point)
  "Is POINT inside of TRIANGLE? Returns -1 if point is outside, 0 if on border, 1 if inside."
  (let* ((dtri (delta ti))
	 (rs (mapcar #'(lambda (vertex)
			 (v- vertex point))
		     (vertices ti)))
	 ;; measure area between r and each edge vector
	 (areas (mapcar #'(lambda (b r) (vec-z (cross r b)))
			dtri rs)))
    (every #'(lambda (x) (>= x 0))
	   areas)))

(defun contains-point (point-index)
  "Return index of triangle that contains POINT-INDEX. -1 for failure."
  (block block-let
    (let ((j 0)
	  (k 0)
	  (i 0))
      (loop while (<= (stat (get-triangle k)) 0) ;; walk down all triangles which are disabled
	    do
	    (block block-loop 
	      (loop ;; check up to three daughters [running i]
	       (setf j (aref (daughters (aref trilist k)) i))
	       (unless (< j 0)
		 (when (> (contains (aref trilist j) (aref global-points point-index))
			  0)
		   (format t "point is contained in triangle j=~a." j)
		   (return-from block-loop)))
	       (setf i (1+ i))
	       (when (= i 3)
		 (format t "point not found in any daughter.")
		 (return-from block-let -1))))
	    (setf k j) ;; new mother
	    )
      k)))
#+nil
(defun point-in-circumcircle (d a b c)
  "Return positive if point D is inside the circumcircle of triangle ABC. Zero if on border. Negative if outside."
  (let* ((cc (circumcircle (make-triangle a b c)))
	 (vr (v- d (center cc)))
	 (radd (dot vr vr)))
    (- (r^2 cc) radd)))
#+nil
(defmethod erase-triangle ((del delaunay) a b c d0 d1 d2)
  "Erase triangle ABC and inactivate it in TRILIST after setting its daughters."
  (let* ((key (logxor (sxhash a)
		      (sxhash b)
		      (sxhash c)))
	 (trihash (trihash del))
	 (value (gethash key trihash)))
    (if value
	(progn
	  (remhash key trihash)
	  (let* ((tri (aref (trilist del) value))
		 (dau (daughters tri)))
	    (setf (aref dau 0) d0
		  (aref dau 1) d1
		  (aref dau 2) d2
		  (stat tri) 0))
	  (setf (ntri del) (1- (ntri del))))
	(error "Non-existant triangle."))))
#+nil
(defmethod store-triangle ((del delaunay) a b c)
  (with-slots (trilist ntree ntri trihash linehash) del
    (setf (aref (trilist del) ntree) 
	  (make-triangle-element a b c)
	   
	  (gethash (logxor (sxhash a) 
			   (sxhash b) 
			   (sxhash c))
		   trihash)
	  ntree
	   
	  (gethash (- (sxhash b) (sxhash c)) linehash) a
	  (gethash (- (sxhash c) (sxhash a)) linehash) b
	  (gethash (- (sxhash a) (sxhash b)) linehash) c
	  ntree (1+ ntree)
	  ntri (1+ ntri))
    (1- ntree)))
#+nil
(defmethod insert-point ((del delaunay) r)
  "R is an index for global-points."
  (let* ((tno (contains-point del r))
	 (tri (if (< tno 0)
		  (error "Points should be fuzzed.")
		  (aref (trilist del) tno)))
	 (vs (vertices tri))
	 (i (first vs))
	 (j (second vs))
	 (k (third vs))
	 (d0 (store-triangle del r i j))
	 (d1 (store-triangle del r j k))
	 (d2 (store-triangle del r k i))
	 (tasks ()))
    (push (list r i j) tasks)
    (push (list r j k) tasks)
    (push (list r k i) tasks)
    (erase-triangle del i j k d0 d1 d2)
    (loop named task-loop
	  while tasks
	  do
	  (destructuring-bind (s i j)
	      (pop tasks)
	    (let* ((key (- (sxhash j)
			   (sxhash i)))	;; look up fourth point
		   (l (gethash key (linehash del))))
	      (when l
		(when (> (point-in-circumcircle (get-global-point l) ;; needs legalization
						(get-global-point j)
						(get-global-point s)
						(get-global-point i))
			 0)
		  (let* ((d0 (store-triangle del s l j)) ;; create 2 new triangles
			 (d1 (store-triangle del s i l)))
		    (erase-triangle del s i j d0 d1 -1)	;; delete old triangles
		    (erase-triangle del l j i d0 d1 -1)
		    (remhash (- (sxhash i) ;; erase line in both directions
				(sxhash j)) (linehash del))
		    (remhash (- (sxhash j)
				(sxhash i)) (linehash del))
		    (push (list s l j) tasks) ;; two new edges need checking
		    (push (list s i l) tasks)))))))))
  
;;;; hash table stuff

(defun triangle-key (tri)
  "Generate a hash key for a triangle so that it can be found in the hash if the three points are given in any order."
  (logxor (sxhash (get-triangle-a tri))
	  (sxhash (get-triangle-b tri))
	  (sxhash (get-triangle-c tri))))

(defun line-key (ai bi)
  "Generate a hash key for a line (between the points A and B). Used to find the other point of the triangle if only one edge is given."
  (- (sxhash ai)
     (sxhash bi)))


(defun store-triangle (a b c)
  "Store the points A, B and C into ARRAY-POINTS and the triangle into ARRAY-TRIANGLES. A reference of the triangle goes into the triangle hash table HT-TRIANGLES and store a reference to each of the points with their opposing edge in another hash HT-LINES."
  (let* ((ai (append-point a))
	 (bi (append-point b))
	 (ci (append-point c))
	 (ti (append-triangle (triangle ai bi ci))))
    (setf (gethash (triangle-key ti) ht-triangles) ti
	  ;; store ccw line and opposing point, circling points
	  (gethash (line-key bi ci) ht-lines) ai
	  (gethash (line-key ai bi) ht-lines) ci
	  (gethash (line-key ci ai) ht-lines) bi)
    ti))


(defun clear-all ()
  (setf n-triangles 0
	n-points 0)
  (clrhash ht-triangles)
  (clrhash ht-lines))


(defun print-triangles ()
  (loop for val being each hash-value in ht-triangles
	using (hash-key key)
	collecting (list key (get-triangle val))))

(defun print-lines ()
  (loop for val being each hash-value in ht-lines
	using (hash-key key)
	collecting (list key (get-point val))))

(defun print-points ()
  (loop for i below n-points collecting
	(get-point i)))

#|
(clear-all)

(store-triangle (v -1 0) (v 1 0) (v 0 2))
(store-triangle (v .4 .3) (v .6 .22) (v .3 .1))
(store-triangle (v .6 .3) (v .2 .5) (v .3 .2))

(hash-table-count ht-triangles)
(hash-table-count ht-lines)

(print-points)
(print-triangles)
(print-lines)

(centroid 0)
(delta 0)
(incircle 0)
(circumcircle 0)
(contains 0 (center (incircle 0)))

(run)
|#
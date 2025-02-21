;;;; Copyright 2025 Google LLC
;;;; 
;;;; Licensed under the Apache License, Version 2.0 (the "License");
;;;; you may not use this file except in compliance with the License.
;;;; You may obtain a copy of the License at
;;;; 
;;;;     https://www.apache.org/licenses/LICENSE-2.0
;;;; 
;;;; Unless required by applicable law or agreed to in writing, software
;;;; distributed under the License is distributed on an "AS IS" BASIS,
;;;; WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
;;;; See the License for the specific language governing permissions and
;;;; limitations under the License.


;;;; "Life at the Boundary of Chemical Kinetics and Program Execution":
;;;; Core Framework Code


;; Make a new vector by populating its entries with values
;; obtained by mapping a function over its indices.
(define (make-vector-iter n f)
  (let ((vec (make-vector n #f)))
    (n-times (k n) (vector-set! vec k (f k)))
    vec))

;; Sequence-probability-table sp-table
;; #(sp-table v-striding v-probs)
(defobj (sp-table-striding 'sp-table x) (vector-ref x 1))
(defobj (sp-table-probs 'sp-table x) (vector-ref x 2))
(define (sp-table striding probs) (vector 'sp-table striding probs))
(define (sp-table? x)
  (and (vector? x)
       (= (vector-length x) 3)
       (eq? (vector-ref x 0) 'sp-table)))

;;;; === Tape-Views ===
;;
;; The left-list / right-list are such that the farthest-out item comes first.
;; Index zero sits on the right-list.
(defobj (t-view-l-len 't-view tv) (vector-ref tv 1))
(defobj (t-view-r-len 't-view tv) (vector-ref tv 2))
(defobj (t-view-l-list 't-view tv) (vector-ref tv 3))
(defobj (t-view-r-list 't-view tv) (vector-ref tv 4))
(defobj (t-view-visible-size 't-view tv)
  (+ (vector-ref tv 1) (vector-ref tv 2)))

(define (t-view l-len r-len l-list r-list)
  (vector 't-view l-len r-len l-list r-list))
(define (t-view-from-lists l-list r-list)
  (vector t-view (length l-list) (length r-list) l-list r-list))
(define (t-view? x)
  (and (vector? x)
       (= (vector-length x) 5)
       (eq? (vector-ref x 0) 't-view)))

;; The fundamental operation on a tape-view is to perform a lookup
;; for a given index. We have these cases:
;;
;; - We already unfolded the tape to the target-position and
;;   just can determine the value.
;;
;; - We do not have the value, but there is an unfold-operation
;;   that brings us closer to finding it. This is either a left-unfold
;;   or a right-unfold, and we know the prefix, respectively suffix.
;;
;; - Tape-view-lookup is an internal operation that needs some context-helpers.

;; Helper: Given two lists xs, ys, determine the sequence-rank of the substring
;; formed by the first k entries from (reversed (xs ++ (reversed ys))).
;; We know the length of xs and ys in advance.
;;
;; "Generic case" example:
;; We have: size-a = 10, cl-k = 7 (so, need 6 elements),
;;          xs = (2 3 4), ys = (5 6 7 8 9)
;; The tape-sequence hence is (asterisks mark index-zero):
;; ... 5 6 7 8 9 *4* 3 2 ? ...
;; The 6-prefix of the '?' is: 7 8 9 4 3 2, which has rank 789432.
(define (prefix-rank-1 tv size-a cl-k)
  (let ((xs (t-view-r-list tv))
        (ys (t-view-l-list tv))
        (len-xs (t-view-r-len tv))
        (len-ys (t-view-l-len tv)))
    (let loop-xs ((rest-xs xs) (rank-now 0) (scaling 1) (syms-todo (- cl-k 1)))
      (cond
       ((= syms-todo 0) rank-now) ; done.
       ((not (null? rest-xs))  ; we still have items in xs
        (loop-xs (cdr rest-xs)
                 (+ (* (car rest-xs) scaling) rank-now)
                 (* scaling size-a)
                 (- syms-todo 1)))
       (#t
        ;; We still have to scan, but xs have been consumed.
        ;; We need to process the final syms-todo elements of
        ;; the ys-list.
        (let loop-ys ((rest-ys (list-tail ys (- len-ys syms-todo)))
                      (rank-now rank-now)
                      (syms-done-ys 0))
          (if (= syms-done-ys syms-todo)
              rank-now
              (loop-ys (cdr rest-ys)
                       (+ rank-now
                          (* (car rest-ys) (expt size-a (- cl-k 2 syms-done-ys))))
                       (+ 1 syms-done-ys)))))))))


; We also need the counterpart, suffix-rank.
;; Sticking with the generic example above, the tape-sequence (again) is
;; ... ? 5 6 7 8 9 *4* 3 2 ...
;; If we asked for a length-7 suffix of the '?', the rank of that would be 5678943.
(define (suffix-rank-1 tv size-a cl-k)
  (let ((xs (t-view-r-list tv))
        (ys (t-view-l-list tv))
        (len-xs (t-view-r-len tv))
        (len-ys (t-view-l-len tv)))
    (let loop-ys ((rest-ys ys) (rank-now 0) (syms-todo (- cl-k 1)))
      (cond
       ((= syms-todo 0) rank-now)
       ((not (null? rest-ys))
        (loop-ys (cdr rest-ys)
                 (+ (* rank-now size-a) (car rest-ys))
                 (- syms-todo 1)))
       (#t
        ;; ran out of ys, need to get some xs from the tail of that list.
        (let loop-xs ((rest-xs (list-tail xs (- len-xs syms-todo)))
                      (rank-now (* rank-now (expt size-a syms-todo)))
                      (syms-done-xs 0))
          (if (= syms-done-xs syms-todo)
              rank-now
              (loop-xs
               (cdr rest-xs)
               (+ rank-now (* (car rest-xs)
                              (expt size-a syms-done-xs)))
               (+ 1 syms-done-xs)))))))))

;; Return either the entry, or the default if tape was not unfolded
;; far enough yet.
(define (tv-lookup tv index size-a cl-k)
  (check-tag 't-view tv)
  (if (>= index 0)  ; right-lookup
      (let ((r-len (t-view-r-len tv))
            (r-list (t-view-r-list tv)))
        (if (< index r-len)
            ;; we have an entry for this.
            (list-ref r-list (- r-len 1 index))
            '()))
      ;; otherwise, left-lookup.
      (let* ((l-len (t-view-l-len tv))
             (l-list (t-view-l-list tv))
             (offset (+ index l-len)))
        (if (>= offset 0)
            ;; we have an entry for this.
            (list-ref l-list offset)
            '()))))

;; Extend a tape-view. The new entry that becomes available has value `value`.
(define (tv-extend tv to-right? value)
  (if to-right?
      (t-view (t-view-l-len tv) (+ 1 (t-view-r-len tv))
              (t-view-l-list tv) (cons value (t-view-r-list tv)))
      (t-view (+ 1 (t-view-l-len tv)) (t-view-r-len tv)
              (cons value (t-view-l-list tv)) (t-view-r-list tv))))


;; Copy the data from a tape-view to a buffer-vector.
;; (Should be used only in the (slow) "reference"
;; implementation of the algorithm).
(define (tv-copy-to-buffer tv u16-vec)
  (let loop-l
      ((todo (t-view-l-list tv))
       (index 0))
    (if (not (null? todo))
        (begin
          (u16-vector-set! u16-vec index (car todo))
          (loop-l (cdr todo) (+ 1 index)))))
  (let loop-r
      ((todo (t-view-r-list tv))
       (index (+ (t-view-l-len tv) (t-view-r-len tv) -1)))
    (if (not (null? todo))
        (begin
          (u16-vector-set! u16-vec index (car todo))
          (loop-r (cdr todo) (- index 1))))))


(define (tv-as-index tv size-a)
  (let
    ((index-left
      (let loop-l
          ((todo (t-view-l-list tv))
           (current 0))
        (if (null? todo)
            current
            (loop-l (cdr todo)
                    (+ (* size-a current) (car todo))))))
     (num-r-total (t-view-r-len tv)))
    (let loop-r
        ((todo (t-view-r-list tv))
         (i 0)
         (current (* index-left (expt size-a num-r-total))))
      (if (null? todo)
          current
          (loop-r
           (cdr todo) (+ 1 i)
           (+ current (* (car todo) (expt size-a i))))))))
       

;; A tape-view-pair (tvp) gives us the original tape-state and the
;; tape-state due to program-action.
(defobj (tv-pair-orig 'tv-pair tvp) (vector-ref tvp 1))
(defobj (tv-pair-adjusted 'tv-pair tvp) (vector-ref tvp 2))
(define (tv-pair tv-orig tv-adjusted)
  (vector 'tv-pair tv-orig tv-adjusted))
(define (tv-pair? x)
  (and (vector? x)
       (= (vector-length x) 3)
       (eq? (vector-ref x 0) 'tv-pair)))


;; Extending a tape-view pair extends both the 'original-state' tape-view
;; and the 'adjusted-state' tape-view.
(define (tv-pair-extend tvp to-right? value)
  (let ((tv-orig (tv-pair-orig tvp))
        (tv-adjusted (tv-pair-adjusted tvp)))
    (let ((tv-orig-extended (tv-extend tv-orig to-right? value)))
      (if (eq? tv-orig tv-adjusted)
          ;; Optimization: If the original and adjusted view are not
          ;; just equal, but identical, the extended tv-pair will also
          ;; have identical t-view-s. That way, we can later answer the
          ;; "did the tape stay completely untouched?" question
          ;; with a simple eq? check.
          (tv-pair tv-orig-extended tv-orig-extended)
          (tv-pair
           tv-orig-extended
           (tv-extend (tv-pair-adjusted tvp) to-right? value))))))


;; Given a sufficiently-extended tape-view-pair,
;; produce an adjusted one where the value at the target-index
;; is the given value.
(define (tv-pair-adjust-extended tvp index value)
  (let ((tv-orig (tv-pair-orig tvp))
        (tv-adjusted (tv-pair-adjusted tvp)))
    (if (>= index 0)
        ;; Need to adjust right-list.
        (if (< (t-view-r-len tv-adjusted) (+ 1 index))
            (error
             "Assertion A1 violated: Tape should have been right-extended sufficiently far."
             `((tv-pair . ,tvp)
               (index . ,index)))
            (let loop ((r-todo (reverse (t-view-r-list tv-adjusted)))
                       (r-done '())
                       (k 0))
              (if (null? r-todo)
                  (tv-pair tv-orig
                           (t-view
                            (t-view-l-len tv-adjusted)
                            (t-view-r-len tv-adjusted)
                            (t-view-l-list tv-adjusted)
                            r-done))
                  (loop (cdr r-todo)
                        (cons (if (= k index) value (car r-todo)) r-done)
                        (+ 1 k)))))
        ;; Otherwise, need to adjust left-list.
        (if (< (t-view-l-len tv-adjusted) (- index))
            (error
             "Assertion A2 violated: Tape should have been left-extended sufficiently far."
             `((tv-pair . ,tvp)
               (index . ,index)))
            (let loop ((l-todo (reverse (t-view-l-list tv-adjusted)))
                       (l-done '())
                       (k -1))
              (if (null? l-todo)
                  (tv-pair tv-orig
                           (t-view
                            (t-view-l-len tv-adjusted)
                            (t-view-r-len tv-adjusted)
                            l-done
                            (t-view-r-list tv-adjusted)))
                  (loop (cdr l-todo)
                        (cons (if (= k index) value (car l-todo)) l-done)
                        (- k 1))))))))


;; Finding the affected-range after user fn-eval execution is done.
;;
;; (We need to do this in order to determine by how much to further
;; unfold the tape(s) to have all changes covered by our probabilistic
;; tape-content model.)
;;
;; Depending on that affected-range, we may have to expand further
;; so that we have full context-length (cl-k) coverage of every
;; change in both directions.
;; The most straightforward way to get to full coverage is by
;; splitting the multiverse with further lookups, but we also might
;; want to speed up this step by instead building some custom
;; looping, since there no longer is any evaluation happening
;; we would have to return to.
(define (tv-pair-maybe-max-change-index tvp)
  ;; First, process the right-facing list.
  (let loop-r ((todo-orig (t-view-r-list (tv-pair-orig tvp)))
               (todo-adjusted (t-view-r-list (tv-pair-adjusted tvp)))
               (index (- (t-view-r-len (tv-pair-orig tvp)) 1)))
    (if (not (null? todo-orig))
        (if (not (eqv? (car todo-orig) (car todo-adjusted)))
            ;; We have our index.
            index
            ;; Otherwise, they are the same. Recurse.
            (loop-r (cdr todo-orig) (cdr todo-adjusted) (- index 1)))
        ;; We are done processing r-list; it might still be that
        ;; the highest-index change has a negative index.
        ;; Here, we need to find the closest-to-the-end change.
        (let loop-l ((todo-orig (t-view-l-list (tv-pair-orig tvp)))
                     (todo-adjusted (t-view-l-list (tv-pair-adjusted tvp)))
                     (index (- (t-view-l-len (tv-pair-orig tvp))))
                     (candidate-result '()))
          (if (null? todo-orig) candidate-result
              (loop-l (cdr todo-orig) (cdr todo-adjusted) (+ 1 index)
                      (if (eqv? (car todo-orig) (car todo-adjusted))
                          candidate-result
                          index)))))))

;; Finding the smallest index where there has been a tape-adjustment.
(define (tv-pair-maybe-min-change-index tvp)
  (if (eq? (tv-pair-orig tvp) (tv-pair-adjusted tvp))
      '()  ; short-cut: identical objects mean: tape was untouched.
      ;; First, process the left-facing list.
      (let loop-l ((todo-orig (t-view-l-list (tv-pair-orig tvp)))
                   (todo-adjusted (t-view-l-list (tv-pair-adjusted tvp)))
                   (index (- (t-view-l-len (tv-pair-orig tvp)))))
        (if (not (null? todo-orig))
            (if (not (eqv? (car todo-orig) (car todo-adjusted)))
                index
                (loop-l (cdr todo-orig) (cdr todo-adjusted) (+ index 1)))
            ;; We are done processing the l-list and found no change.
            ;; Need to scan the r-list for the smallest-index change.
            (let loop-r ((todo-orig (t-view-r-list (tv-pair-orig tvp)))
                         (todo-adjusted (t-view-r-list (tv-pair-adjusted tvp)))
                         (index (- (t-view-r-len (tv-pair-orig tvp)) 1))
                         (candidate-result '()))
              (if (null? todo-orig)
                  candidate-result
                  (loop-r (cdr todo-orig) (cdr todo-adjusted)
                          (- index 1)
                          (if (eqv? (car todo-orig) (car todo-adjusted))
                              candidate-result index))))))))


;;;; === Probabilistic Tape-composition Model ===

;; We model tape-content via a Markov processes.
;; 
;; We start from k-sequence probabilities, and these give
;; us (k-1)-prefix/suffix Markov transition probabilities.
;; 
;; Size of the alphabet is `size-a`, and `cl-k` is the 'correlation
;; length', i.e. length of subsequences that we track.
;; 
;; For updates, we need to know the margial probabilities that sum over
;; the first, respectively last index of an array of (implicit)
;; shape `[{size-a}]*{cl-k}` that is stored in flattened form.
;; These again are sp-tables, but with one axis less.
(define (sp-table-marginal sp-t sum-over-last-index?)
  (let* ((striding (sp-table-striding sp-t))
         (probs (sp-table-probs sp-t))
         (cl-k (vector-length striding))
         (cl-k-shortened (- cl-k 1))
         (size-a (/ (fp-vector-length (sp-table-probs sp-t))
                    (vector-ref striding 0)))
         (data-size-reduced (vector-ref striding 0))
         (shortened-striding
          (let ((v (make-vector cl-k-shortened 0)))
            (n-times (j cl-k-shortened)
                     (vector-set! v j (vector-ref striding (+ 1 j))))
            v))
         (reduced-probs (make-fp-vector data-size-reduced))
         (co-k-scaling (if sum-over-last-index? size-a 1))
         (k-scaling (if sum-over-last-index? 1 data-size-reduced)))
    (n-times (co-k data-size-reduced)
      (let ((offset (* co-k co-k-scaling)))
        (let loop ((total 0) (k 0))
          (if (= k size-a)
              (fp-vector-set! reduced-probs co-k total)
              (loop (+ total (fp-vector-ref probs (+ offset (* k-scaling k))))
                    (+ k 1))))))
    (sp-table shortened-striding reduced-probs)))

;; subsequence-probability info: sp-table with marginals and alphabet.
;; Here, we use the following properties:
;; - For a Markov process, it does not matter if we get marginals
;;   by summing over the first or the last index.
;; - We can populate the in-use marginals lazily.
(define (sp-info symbols sp-table)
  (let* ((cl-k (vector-length (sp-table-striding sp-table)))
         (result (make-vector (+ 3 cl-k) '())))
    (vector-set! result 0 'sp-info)
    (vector-set! result 1 symbols)
    ;; + 2, since object structure is:
    ;; #(sp-table #(a b ...) {0-marginals} {1-marginals} ... {k-marginals})
    (vector-set! result (+ 2 cl-k) sp-table)
    result))
(defobj (sp-info-symbols 'sp-info x) (vector-ref x 1))
(defobj (sp-info-table 'sp-info x) (vector-ref x (- (vector-length x) 1)))
(defobj (sp-info-cl-k 'sp-info x) (- (vector-length x) 3))
(define (sp-info-marginals the-sp-info subseq-length)
  (check-tag 'sp-info the-sp-info)
  (let ((smallest-so-far
         (let find-smallest-we-have ((k subseq-length))
           (if (not (null? (vector-ref the-sp-info (+ 2 k))))
               k
               (find-smallest-we-have (+ 1 k))))))
    (let fill-in-needed ((k (- smallest-so-far 1)))
      (if (>= k subseq-length)
          (begin
            (vector-set!
             the-sp-info (+ 2 k)
             (sp-table-marginal (vector-ref the-sp-info (+ 3 k)) #t))
            (fill-in-needed (- k 1)))))
    (vector-ref the-sp-info (+ 2 subseq-length))))


;;;; === Multiverse State ===

;; A multiverse-state may have a "cause for existence"
;; which is one of these different possible things:
;;
;; - '() if the cause is missing, i.e. a tape-lookup was performed,
;;   and this already caused a multiverse-split (with a resulting
;;   world-probability reduction), relative to which we now work,
;;   until the multiverse stack is unwound.
;;
;; - a triplet #(data-tape? index continuation), which specifies
;;   a lookup that may require further world-splitting to be satisfied.
;;
;; - A vector #(probs+options maybe-index continuation), which specifies
;;   what continuation to take from a
;;   (choose-1 ... #(...probs...) #(...options...) ...) form.
;;   If maybe-index is '(), this has not been unfolded yet.

(defobj (mv-state-prob 'mv-state mvs) (vector-ref mvs 1))
;; tape-view-pair
(defobj (mv-state-tvp 'mv-state mvs data-tape?)
  (vector-ref mvs (if data-tape? 3 2)))
(defobj (mv-state-cause 'mv-state mvs) (vector-ref mvs 4))
(define (mv-state prob p-tvp d-tvp maybe-cause)
  (vector 'mv-state prob p-tvp d-tvp maybe-cause))
(define (mv-state-neutered x)
  (mv-state (vector-ref x 1) (vector-ref x 2)
            (vector-ref x 3) '()))
(define (mv-state? x)
  (and (vector? x)
       (= (vector-length x) 5)
       (eq? (vector-ref x 0) 'mv-state)))


;; Helper for make-new-mv-state: Map a seq-index to a r-list.
(define (mv-state-r-list cl-k size-a seq-index)
  ;; if size-a=10 and cl-k=5, then seq-index 12345
  ;; will become the r-list '(5 4 3 2 1).
  (let loop ((done '())
             (i-todo seq-index)
             (power-todo (- cl-k 1)))
    (if (< power-todo 0) done
        (let* ((pow (expt size-a power-todo))
               (next-id (floor (/ i-todo pow))))
          (loop (cons next-id done)
                (- i-todo (* next-id pow))
                (- power-todo 1))))))

(define (make-new-mv-state size-a cl-k prob p-seq-index d-seq-index)
  (let* ((p-r-list (mv-state-r-list cl-k size-a p-seq-index))
         (d-r-list (mv-state-r-list cl-k size-a d-seq-index))
         (tv-p (t-view 0 cl-k '() p-r-list))
         (tv-d (t-view 0 cl-k '() d-r-list)))
    ;; Note that initially, the original and adjusted tape-view
    ;; are identically the same object.
    (mv-state prob (tv-pair tv-p tv-p) (tv-pair tv-d tv-d) '())))


;; Unfolding the multiverse-state due to a tape-get.
;; This will be called *only* by mv-stack-process when a
;; need-for-unfolding has been recognized.
(define (mv-state-unfold-for-tape-get
         the-mv-state cause
         the-sp-info
         mv-stack-to-add-to)
  ;; We did find a need to unfold the tape.
  (let* ((data-tape? (vector-ref cause 0))
         (index (vector-ref cause 1))
         (tvp (mv-state-tvp the-mv-state data-tape?))
         (tv-orig (tv-pair-orig tvp))
         (tv-adjusted (tv-pair-adjusted tvp))
         ;; Unfolding will move us one step closer to the goal of
         ;; being able to resolve the lookup encoded in `cause`.
         ;; As such, we need to know:
         ;; - What is the size of the alphabet?
         ;; - Do we need to unfold to the right or to the left?
         ;; - What is the "effective correlation length" at play?
         ;;   (This is the minimum of cl-k and
         ;;   (+ 1 as-many-entries-as-we-have)).
         (size-a (vector-length (sp-info-symbols the-sp-info)))
         (tv-size (t-view-visible-size tv-orig))
         (cl-k-effective (min (sp-info-cl-k the-sp-info) (+ 1 tv-size)))
         (probs-seq
          (sp-table-probs (sp-info-marginals the-sp-info cl-k-effective)))
         (probs-marginals
          (sp-table-probs (sp-info-marginals the-sp-info (- cl-k-effective 1))))
         (subseq-rank
          (if (>= index 0)
              (prefix-rank-1 tv-orig size-a cl-k-effective)
              (suffix-rank-1 tv-orig size-a cl-k-effective)))
         (p-marginal (fp-vector-ref probs-marginals subseq-rank))
         (subseq-rank-offset
          (if (>= index 0)
              (* size-a subseq-rank)  ; subseq is a prefix.
              subseq-rank  ; subseq is a suffix
              ))
         (k-scaling
          (if (>= index 0)
              1  ; subseq is a prefix.
              (fp-vector-length probs-marginals)  ; subseq is a suffix.
              )))
    ;; We need to compute the extended mv-stack where we added all the
    ;; sub-split cases.
    (let loop-add-mv-state
        ((k 0)
         ;; TODO: Decide whether (and in what form) we want to have
         ;; such an internal consistency check enabled. The
         ;; observation is that with such a check, violations get
         ;; detected due to the ODE-integrator's interpolation
         ;; sometimes feeding us invalid probability data - unless we
         ;; constrain the maximum time-step. This however only seems
         ;; to happen in "harmless" situations where the total
         ;; probability-contribution is many orders of magnitude below
         ;; what's relevant. Having such a check would also protect
         ;; against mis-use of (vector-choose ...) - but this is only
         ;; one of many possible mis-uses where we cannot do anything
         ;; about the others. In any case, if we want to do something
         ;; here, we likely should not (error ...) out, but somehow
         ;; communicate upwards that there was a problem.
         ;;
         ;; (p-relative-total 0.0)
         ;; (p-relative-collected '())
         (mv-stack-now mv-stack-to-add-to))
      (if (= k size-a)
          (begin
            ;; (if (not (<= (abs (- p-relative-total 1.0)) 0.001))
            ;;     ;; This would e.g. trigger on a vec-choose where
            ;;     ;; probbabilities do not sum to 1.
            ;;     ;; Check is deliberately somewhat lax.
            ;;     (error
            ;;      "Assertion A3 violated: relative probabilities must add to 1.0 - got:"
            ;;      p-relative-total p-relative-collected)
            ;;     mv-stack-now  ; done, return new mv-stack with additions.
            ;;     )
            mv-stack-now)
          (let* ((seq-rank (+ subseq-rank-offset (* k-scaling k)))
                 (p-here (max 0.0 (fp-vector-ref probs-seq seq-rank)))
                 ;; The max() below is to safeguard against weird
                 ;; numerical rounding issues.
                 (p-relative (if (= p-here 0) 0
                                 ;; (max ...) to handle numerical noise.
                                 (/ p-here (max p-here p-marginal))))
                 (prob (mv-state-prob the-mv-state))
                 (prob-new (* prob p-relative)))
            (if (> prob-new 0)
                (let* ((the-p-tvp (mv-state-tvp the-mv-state #f))
                       (the-d-tvp (mv-state-tvp the-mv-state #t))
                       (new-p-tvp
                        (if data-tape?
                            the-p-tvp
                            (tv-pair-extend the-p-tvp (>= index 0) k)))
                       (new-d-tvp
                        (if data-tape?
                            (tv-pair-extend the-d-tvp (>= index 0) k)
                            the-d-tvp))
                       (mv-state-new
                        (mv-state prob-new new-p-tvp new-d-tvp
                         ;; cause is being inherited under splitting.
                         cause)))
                  (loop-add-mv-state (+ k 1)
                                     ;; (+ p-relative-total p-relative)
                                     ;; (cons p-relative p-relative-collected)
                                     (cons mv-state-new mv-stack-now)))
                ;; otherwise, this universe is "impossible", and we just iterate.
                (loop-add-mv-state (+ k 1)
                                   ;; (+ p-relative-total p-relative)
                                   ;; (cons p-relative p-relative-collected)
                                   mv-stack-now)))))))


;; Unfolding the multiverse-state due to a choose.
;; This will be called *only* by mv-stack-process when a
;; need-for-unfolding has been recognized.
(define (mv-state-unfold-for-choose
         the-mv-state cause
         the-sp-info
         mv-stack-to-add-to)  
  (let* ((probs (car (vector-ref cause 0)))
         (options (cdr (vector-ref cause 0)))
         (num-options (vector-length probs)))
    (let loop-add-mv-state
        ((k 0)
         ;; TODO: As above - decide if we want such a consistency check here.
         ;; (p-relative-total 0.0)
         ;; (p-relative-collected '())
         (mv-stack-now mv-stack-to-add-to))
      (if (= k num-options)
          (begin
            ;; (if (not (<= (abs (- p-relative-total 1.0)) 0.001))
            ;;     (error
            ;;      "Assertion A4 violated: relative probabilities must add to 1.0 - got:"
            ;;      p-relative-total p-relative-collected)
            ;;     mv-stack-now  ; done, return new mv-stack with additions.
            ;;     )
            mv-stack-now)
          (let ((mv-state-new
                 (mv-state (* (max 0.0 (vector-ref probs k))
                              (mv-state-prob the-mv-state))
                           (mv-state-tvp the-mv-state #f)
                           (mv-state-tvp the-mv-state #t)
                           ;; Rebuilding a new cause with the index set to k.
                           (vector (vector-ref cause 0) k (vector-ref cause 2)))))
            (loop-add-mv-state (+ k 1)
                               ;; (+ p-relative-total (max 0.0 (vector-ref probs k)))
                               ;; (cons (max 0.0 (vector-ref probs k)) p-relative-collected)
                               (cons mv-state-new mv-stack-now)))))))


;;;; === Algorithmic-ODE Core Mechanism ===
;;
;; Basic idea: Given some symbol-alphabet (as a vector) and some
;; interpreter-function that is of the form...:
;;
;; (lambda (tape-get tape-set!)
;;   ... execute the program ...)
;;
;; we want to have a function that maps a `sp-table` to the d/dt
;; for that sp-table.


;; Processing the multiverse-stack
;; (Actually, as we will see: Unfolding the top-of-stack)
;;
;; - (!!!) This cleanup gets triggered by fn-eval finishing its work.
;;
;; - If the stack is empty, we are done.
;;
;; - If the top stack entry does not have a cause / requesting-continuation,
;;   it exists because we performed a lookup successfully that became the
;;   basis for a state-of-world relative to which further splitting
;;   and cleaning-up would have happened. But "everything that needed
;;   splitting was processed", so nothing is to happen with this
;;   intermediate checkpoint anymore. We can remove it.
;;
;; - If the top stack entry has a requesting-continuation, ...
;;   see if we can actually do the element-lookup that asked
;;   for this requesting-continuation.
;;   - If yes, replace the top stack entry with a "neutered"
;;     one that no longer has the continuation: this will be
;;     the top-world that is active and seen by the continuation
;;     that has the successful lookup; this will ultimately be
;;     removed-from-stack by termination of fn-eval.
;;     ALSO, leave by activating the continuation.
;;   - If no, remove the element with requesting-continuation from
;;     world-stack and add size-a many entries to world-stack
;;     that have the same requesting-continuation but different
;;     mv-state state-of-the-world, as-unfolded.

(define (mv-stack-process mv-stack-ref the-sp-info)
  (cond
   ((null? (car mv-stack-ref))  ; We are done, nothing to do.
    #t)
   ((null? (mv-state-cause (caar mv-stack-ref)))
    ;; Top-of-mv-stack is a 'resolved checkpoint'.
    ;; This does not ask for further sub-splitting work,
    ;; so we can clean it up by removing it from the mv-stack.
    (set-car! mv-stack-ref (cdar mv-stack-ref))
    (mv-stack-process mv-stack-ref the-sp-info))
   (#t
    ;; Top-of-mv-stack has a cause-continuation.
    (let* ((this-mv-state (caar mv-stack-ref))
           (cause (mv-state-cause this-mv-state))
           (index (vector-ref cause 1))
           (cont (vector-ref cause 2)))
      (if (not (pair? (vector-ref cause 0)))
          ;;;; Cause comes from a tape-get.
          (let* ((data-tape? (vector-ref cause 0))
                 (tvp (mv-state-tvp this-mv-state data-tape?))
                 (size-a (vector-length (sp-info-symbols the-sp-info)))
                 (cl-k (vector-length
                        (sp-table-striding (sp-info-table the-sp-info))))
                 (maybe-lookup-result
                  (tv-lookup (tv-pair-adjusted tvp) index size-a cl-k)))
            (if (not (null? maybe-lookup-result))
                ;; We could perform the lookup, the answer is clear,
                ;; no need to split the world further.
                ;; We need to replace this top-of-world-state with
                ;; "what now is the current checkpoint-state-of-the-world"
                ;; that does not have a cause-continuation.
                (let ((this-mv-state-neutered
                       (mv-state-neutered this-mv-state)))
                  ;; Replace top element.
                  (set-car! mv-stack-ref
                            (cons this-mv-state-neutered (cdar mv-stack-ref)))
                  ;; Actually invoke continuation; further execution happens
                  ;; in a world where current-unfolding-state is a "checkpoint".
                  (cont  ; replaces current-continuation with cont,
                         ; no further execution here.
                   maybe-lookup-result)))
            ;; Otherwise, lookup was not successful.
            ;; We are removing the current mv-state,
            ;; and adding mv-states for all possible unfoldings.
            (begin
              (set-car! mv-stack-ref
                        (mv-state-unfold-for-tape-get
                         this-mv-state cause the-sp-info (cdar mv-stack-ref)))
              ;; Here, we can now recurse into processing the stack
              ;; post-unfolding.
              (mv-stack-process mv-stack-ref the-sp-info)
              ))
          ;;;; Otherwise, cause comes from a (choose ...)
            (if (null? index)
                (begin
                  ;; This (choose ...) was not yet unfolded. Unfold it.
                  (set-car! mv-stack-ref
                            (mv-state-unfold-for-choose
                             this-mv-state cause the-sp-info (cdar mv-stack-ref)))
                  ;; Here, we can now recurse into processing the stack
                  ;; post-unfolding.
                  (mv-stack-process mv-stack-ref the-sp-info))
                ;; Otherwise, we are in a world described by the "neutered"
                ;; 'choose' element.
                ;; Replace top multiverse stack element with "neutered" form
                ;; and invoke continuation.
                (let* ((probs (car (vector-ref cause 0)))
                       (options (cdr (vector-ref cause 0)))
                       (this-mv-state-neutered (mv-state-neutered this-mv-state)))
                  (set-car! mv-stack-ref
                            (cons this-mv-state-neutered (cdar mv-stack-ref)))
                  (cont
                   ;; replaces current-continuation, as above.
                   (vector-ref options index)))))))))

;; Internal tape-get. This performs a tape-lookup which,
;; if successful, simply returns a value, and if not
;; (since the tape has not been unfolded deeply enough),
;; puts a new mv-state that represents the unfold-request
;; on the top of the stack and exits out via `exit-out-cont`,
;; which will process the stack via (mv-stack-process).
(define (tape-get-1 mv-stack-ref the-sp-info data-tape? index exit-out-cont)
  (call-with-current-continuation
   (lambda (tape-get-1-return-cont)
     (let* ((this-mv-state (caar mv-stack-ref))
            (tvp (mv-state-tvp this-mv-state data-tape?))
            (tv-adjusted (tv-pair-adjusted tvp))
            (size-a (vector-length (sp-info-symbols the-sp-info)))
            (cl-k (vector-length (sp-table-striding
                                  (sp-info-table the-sp-info))))
            (maybe-lookup-result (tv-lookup tv-adjusted index size-a cl-k)))
       (if (not (null? maybe-lookup-result))
           ;; Tape-lookup was successful, we just return
           ;; the value on the tape, via our continuation.
           (tape-get-1-return-cont maybe-lookup-result)
           ;; Otherwise, lookup failed. We have to put the
           ;; lookup-request that needs multiverse-unfolding
           ;; onto the multiverse-stack and then exit out
           ;; via the exit-out-cont.
           (let ((new-mv-state
                  (mv-state
                   (mv-state-prob this-mv-state)
                   (mv-state-tvp this-mv-state #f)
                   (mv-state-tvp this-mv-state #t)
                   (vector data-tape? index tape-get-1-return-cont))))
             (set-car! mv-stack-ref (cons new-mv-state (car mv-stack-ref)))
             (exit-out-cont #f)))))))


;; Internal tape-set:
;;  - First, perform a tape-get that ensures that the tape has been unfolded
;;    to the needed depth (this may split the universe).
;;  - Then, build a new state-of-the-world that has the tape-view-pair
;;    of the tape in question modified in the expected place.
;;  - This then becomes a new "checkpoint", so we push it onto
;;    the multiverse-stack.
(define (tape-set-1 mv-stack-ref
                    the-sp-info data-tape? index value exit-out-cont)
  (let ((unused-lookup-value  ; only needed to ensure we are in a world
                                        ; where the tape is unfolded far enough.
         ;; Here, we can directly forward the exit-out-cont, since the adjustment
         ;; of the tape-view on the multiverse-stack happens in the continuation
         ;; of tape-get-1 right here.
         (tape-get-1 mv-stack-ref the-sp-info data-tape? index exit-out-cont)))
    (let* ((the-mv-state (caar mv-stack-ref))
           (tvp-post-adjustment
            (tv-pair-adjust-extended (mv-state-tvp the-mv-state data-tape?)
                                     index value))
           (new-mv-state
            (mv-state
             (mv-state-prob the-mv-state)
             (if data-tape? (mv-state-tvp the-mv-state #f) tvp-post-adjustment)
             (if data-tape? tvp-post-adjustment (mv-state-tvp the-mv-state #t))
             '()  ; no cause-continuation.
             )))
      (if (not (null? (mv-state-cause the-mv-state)))
          (error
           "Assertion A5 violated: tape-set-1 top-of-stack universe must not have unresolved cause."
           the-mv-state)
          (set-car!
           mv-stack-ref
           (cons new-mv-state (car mv-stack-ref)))))))


;; Performance note: `(choose ...) is computationally costly - in the
;; longer run, we might want to give users an alternative approach that
;; roughly looks as follows:
;;
;; (let ((choose-ab (make-choice-fn '((1.0 a) (1.0 b)))))
;;   (register-problem ...
;;     (let ((a-or-b (choose-ab)))
;;        ...)
;;
;; That way, relevant arg-processing can be done at make-choice-fn
;; invocation time, rather than at-choice-time.


;; Internal choose. This performs a nondeterministic choice.
(define (choose-1 mv-stack-ref the-sp-info probs options exit-out-cont)
  (call-with-current-continuation
   (lambda (choose-1-return-cont)
     (let* ((this-mv-state (caar mv-stack-ref))
            (new-mv-state
             (mv-state
              (mv-state-prob this-mv-state)
              (mv-state-tvp this-mv-state #f)
              (mv-state-tvp this-mv-state #t)
              (vector (cons probs options) '() choose-1-return-cont))))
       (set-car! mv-stack-ref (cons new-mv-state (car mv-stack-ref)))
       (exit-out-cont #f)))))


;; Helper, used inside compute-dy/dt, adjusts the dy/dt-accumulator.
;; vec-buffer-ref is #(vec-orig vec-adjusted size-valid-prefix)
(define (dy/dt-collect-into-accumulator!
         the-sp-info world-prob dy/dt-accumulator vec-buffer-ref)
  (let* ((the-sp-table (sp-info-table the-sp-info))
         (size-a (vector-length (sp-info-symbols the-sp-info)))
         (striding (sp-table-striding the-sp-table))
         (cl-k (vector-length striding)))
    (letrec
        ((seq-index
          (lambda (vec vector-offset)
            (let loop ((total 0)
                       (v-i vector-offset)
                       (i 0))
              (if (= i cl-k) total
                  (loop (+ total (* (vector-ref striding i)
                                    (u16-vector-ref vec v-i)))
                        (+ 1 v-i)
                        (+ 1 i))))))
         (process-vec
          (lambda (size)
            (let ((offset-max (- size cl-k)))
              (let loop ((offset 0))
                (if (<= offset offset-max)
                    (let ((seq-index-orig
                           (seq-index (vector-ref vec-buffer-ref 0)
                                      offset))
                          (seq-index-adjusted
                           (seq-index (vector-ref vec-buffer-ref 1)
                                      offset)))
                      (if (not (= seq-index-orig seq-index-adjusted))
                          ;; There actually is a change here.
                          (begin
                            (fp-vector-set!
                             dy/dt-accumulator seq-index-orig
                             (- (fp-vector-ref dy/dt-accumulator seq-index-orig)
                                world-prob))
                            (fp-vector-set!
                             dy/dt-accumulator seq-index-adjusted
                             (+ (fp-vector-ref dy/dt-accumulator
                                               seq-index-adjusted)
                                world-prob))))
                      (loop (+ offset 1)))))))))
    (process-vec (vector-ref vec-buffer-ref 2)))))

;; If we are done with eval-fn evaluation, we still need to unfold tapes
;; so that every observed tape-change has full coverage, and then collect
;; all changes to get us the dy/dt.
;;
;; vec-buffer-ref: see dy/dt-collect-into-accumulator!
;; This mutable buffer's vectors will be replaced with a sufficiently-extended
;; vectors as-needed, and have the orig-sequence and adjusted-sequence unspooled
;; into them, for windowed processing.
(define (finish-fn-eval-reference  ; "slow reference implementation"
         mv-stack-ref the-sp-info exit-out-cont
         vec-buffer-ref dy/dt-accumulator)
  (let*
      ((the-mv-state (caar mv-stack-ref))
       (debug-flag (cdr mv-stack-ref))  ; TODO: Document use of CDR on stack-ref for debug-flag!
       (cl-k (vector-length (sp-table-striding (sp-info-table the-sp-info))))
       (tvp-p (mv-state-tvp the-mv-state #f))
       (tvp-d (mv-state-tvp the-mv-state #t))
       ;; If there was no change at all, -min and -max are off-but-meaningless,
       ;; since then, change-index-[pd]-(min|max) will be '(), and this is what
       ;; we use for making decisions.
       (unfolded-range-p-min (- (t-view-l-len (tv-pair-orig tvp-p))))
       (unfolded-range-p-max (- (t-view-r-len (tv-pair-orig tvp-p)) 1))         
       (unfolded-range-d-min (- (t-view-l-len (tv-pair-orig tvp-d))))
       (unfolded-range-d-max (- (t-view-r-len (tv-pair-orig tvp-d)) 1))
       (change-index-p-min (tv-pair-maybe-min-change-index tvp-p))
       (change-index-p-max (if (null? change-index-p-min)
                               change-index-p-min
                               (tv-pair-maybe-max-change-index tvp-p)))
       (change-index-d-min (tv-pair-maybe-min-change-index tvp-d))
       (change-index-d-max (if (null? change-index-d-min)
                               change-index-d-min
                               (tv-pair-maybe-max-change-index tvp-d)))
       (expand-tape-l-r
        (lambda (data-tape?
                 change-index-min change-index-max
                 unfolded-range-min unfolded-range-max)
          ;; Expanding the tape so that every cl-k wide window that contains
          ;; a change will be visible.
          (if (not (null? change-index-min))
              ;; We have a min-position for a change on the tape.
              ;; This implies we will also have a max-position.
              ;;
              ;; TODO: expanding the world-state so we have full coverage from
              ;; leftmost-touched-window to rightmost-touched-window can create
              ;; a lot of avoidable effort. What we rather need to do:
              ;; - Generic situation on the tape (. = unfolded, X = adjusted):
              ;;   ..X..X...X.
              ;; - If context-window is 4, we need to first take this world
              ;;   and unfold to have 4-wide coverage to the left, and
              ;;   process everything we can here.
              ;; - We then also need to start from the
              ;;   last-processed-window-start *on the original world-situation*,
              ;;   and unfold-towards-right from that.
              (begin
                (let loop-expand-l
                    ((index (- unfolded-range-min 1))
                     ;; Example: cl-k=5, unfolded-min=-4, change=-2:
                     ;; Need to unfold at indices: -2, -3, -4, -5, -6 to cover
                     ;; a 5-sequence that ends with the change, so here have
                     ;; to do 2 more steps (= 5 - 1 - ((-2) - (-4)))
                     (todo (- cl-k 1 (- change-index-min unfolded-range-min))))
                  (if (> todo 0)
                      (begin
                        (tape-get-1 mv-stack-ref
                                    the-sp-info data-tape? index exit-out-cont)
                        (loop-expand-l (- index 1) (- todo 1)))))
                (let loop-expand-r
                    ((index (+ unfolded-range-max 1))
                     ;; Example: cl-k=5, unfolded-range-max=4, change=2:
                     ;; Need to unfold at indices: 2, 3, 4, 5, 6 to cover
                     ;; a 5-sequence that ends with the change. Also here, we have
                     ;; to do 2 more steps.
                     (todo (- cl-k 1 (- unfolded-range-max change-index-max))))
                  (if (> todo 0)
                      (begin
                        (tape-get-1 mv-stack-ref
                                    the-sp-info data-tape? index exit-out-cont)
                        (loop-expand-r (+ index 1) (- todo 1))))))))))
    (expand-tape-l-r #f
                     change-index-p-min change-index-p-max
                     unfolded-range-p-min unfolded-range-p-max)
    (expand-tape-l-r #t
                     change-index-d-min change-index-d-max
                     unfolded-range-d-min unfolded-range-d-max)
    (let ((collect-dy/dt-contributions
           (lambda (data-tape?)
             ;; Here, we are in a world where both tapes have been unfolded
             ;; as far as they need be.
             ;; There no longer is any reason to split the world any further.
             (if (not (null? (if data-tape?
                                 unfolded-range-d-min
                                 unfolded-range-p-min)))
                 ;; There was a change to the tape in question.
                 ;; (It does not matter that we also did some tape-expanding
                 ;; afterwards.)
                 ;;
                 ;; We need to again lookup the top-of-multiverse-stack state,
                 ;; since that typically will have been changed by expanding
                 ;; the tapes.
                 (let* ((the-mv-state (caar mv-stack-ref))
                        (tvp (mv-state-tvp the-mv-state data-tape?))
                        (tv-orig (tv-pair-orig tvp))
                        (view-size (+ (t-view-l-len tv-orig) (t-view-r-len tv-orig)))
                        (world-prob (mv-state-prob the-mv-state)))
                   (vector-set! vec-buffer-ref 2 view-size)            
                   ;; Ensure the vec-buffer has sufficient size.
                   (if (> view-size (u16-vector-length
                                     (vector-ref vec-buffer-ref 0)))
                       ;; We can afford to be generous with re-allocation size.
                       (begin
                         (vector-set! vec-buffer-ref 0
                                      (make-u16-vector
                                       (min 32 (* view-size 2)) 0))
                         (vector-set! vec-buffer-ref 1
                                      (make-u16-vector
                                       (min 32 (* view-size 2)) 0))))
                   ;; Copy over the data.
                   (tv-copy-to-buffer tv-orig (vector-ref vec-buffer-ref 0))
                   (tv-copy-to-buffer (tv-pair-adjusted tvp)
                                      (vector-ref vec-buffer-ref 1))
                   (if debug-flag
                       (let ((seq-old (make-vector-iter
                                       view-size
                                       (lambda (k)
                                         (u16-vector-ref
                                          (vector-ref vec-buffer-ref 0) k))))
                             (seq-new (make-vector-iter
                                       view-size
                                       (lambda (k)
                                         (u16-vector-ref
                                          (vector-ref vec-buffer-ref 1) k)))))
                         (if (not (equal? seq-old seq-new))
                             (let ((prog-unfolded
                                    (tv-pair-orig
                                     (mv-state-tvp the-mv-state #f))))
                               (display
                                `((p-world . ,world-prob)
                                  (prog . ,(reverse
                                            (t-view-r-list prog-unfolded)))
                                  (data-tape? . ,data-tape?)
                                  (seq-old . ,seq-old)
                                  (seq-new . ,seq-new)))
                               (display "\n")))))
                   (dy/dt-collect-into-accumulator!
                    the-sp-info world-prob dy/dt-accumulator vec-buffer-ref))))))
      (collect-dy/dt-contributions #f)
      (collect-dy/dt-contributions #t)
      )))


;;;; BEGIN CODE FOR: faster inner loop
;;;; We no longer expand windows unnecessarily, or use call/cc magic, for the innermost loop.

;; Given vec-buffer-ref == #(u16-vec0 u16-vec1 ...),
;; check whether the entries in index-range [i-min .. i-max] are all equal.
(define (vec-buffer-range-equal? vec-buffer-ref i-min i-max)
  (let ((v0 (vector-ref vec-buffer-ref 0))
        (v1 (vector-ref vec-buffer-ref 1)))
    (let loop ((i i-min))
      (cond
       ((> i i-max) #t)
       ((= (u16-vector-ref v0 i) (u16-vector-ref v1 i))
        (loop (+ 1 i)))
       (#t #f)))))

;; XXX for debugging.
(define (ddd-vec-buffer-ranges vec-buffer-ref i-min i-max)
  (let* ((v0 (vector-ref vec-buffer-ref 0))
         (v1 (vector-ref vec-buffer-ref 1))
         (d0 (make-vector-iter
              (- (+ 1 i-max) i-min)
              (lambda (k) (u16-vector-ref v0 (+ i-min k)))))
         (d1 (make-vector-iter
              (- (+ 1 i-max) i-min)
              (lambda (k) (u16-vector-ref v1 (+ i-min k))))))
    (cons d0 d1)))

;; Extracting a fused-index from a subrange of indices.
;; The fused-index is the rank of the subsequence of indices
;; in the indicated range.
(define (u16-vec-fused-index u16-vec index-step-factor i-min i-max)
  (let loop ((i i-max)
             (scale 1)
             (total 0))
    (if (< i i-min)
        total
        (loop (- i 1)
              (* scale index-step-factor)
              (+ total (* scale (u16-vector-ref u16-vec i)))))))


(define (dy/dt-accumulate-part
         dy/dt-accumulator
         current-probs
         ;; p-relative-factor is (/ p-world p-middle-part-orig);
         ;; we then multiply this with the Markov-chain-table
         ;; probability of the fully expanded sequence.
         ;; Note that if the middle-part equals the full sequence,
         ;; this is just p-world.
         p-relative-factor 
         base
         offset-middle-part-orig
         offset-middle-part-adjusted         
         scaling-left-part
         num-indices-left
         num-indices-right)
  (let ((i-right-end (expt base num-indices-right))
        (i-left-end (expt base num-indices-left)))
    (let loop-left ((i-left 0))  ; left part makes "slow" loop.
      (if (= i-left i-left-end)
          #t
          (let loop-right ((i-right 0))
            (if (= i-right i-right-end)
                ;; Done with this right-loop, continuing
                ;; with next left-loop iteration.
                (loop-left (+ 1 i-left))
                (let* ((index-contrib-left-right
                        (+ (* scaling-left-part i-left) i-right))
                       (index-accumulator-orig
                        (+ index-contrib-left-right
                           offset-middle-part-orig))
                       (index-accumulator-adjusted
                        (+ index-contrib-left-right
                           offset-middle-part-adjusted))
                       (p-full-sequence-orig
                        (fp-vector-ref current-probs index-accumulator-orig))
                       (p-delta (* p-full-sequence-orig p-relative-factor)))
                  (if (and (> p-delta 0.0)
                           (not (= index-accumulator-orig
                                   index-accumulator-adjusted)))
                      (begin
                        (fp-vector-set!
                         dy/dt-accumulator index-accumulator-orig
                         (- (fp-vector-ref dy/dt-accumulator
                                           index-accumulator-orig)
                            p-delta))
                        (fp-vector-set!
                         dy/dt-accumulator index-accumulator-adjusted
                         (+ (fp-vector-ref dy/dt-accumulator
                                           index-accumulator-adjusted)
                            p-delta))))
                  (loop-right (+ 1 i-right)))))))))


;;;; *** Adding contributions to the dy/dt-accumulator ***
;;
;; In the general case, we have some unfolded tape that covers
;; the leftmost and rightmost change, and may extend beyond it.
;; We then have to slide all possible windows of length cl-k
;; that have overlap with [change-index-min .. change-index-max]
;; over the tape, potentially left-extending and maybe also
;; right-extending at the same time if the window covers
;; not-yet-expanded tape-regions.
;;
;; IMPORTANT: when we do handle a window that does not contain
;; the full "revealed context" on the tape, we must make sure that
;; expansion is done in a way that respects the revealed context.
;; This will then further reduce world-probability.
;;
;; The generic situation is as follows - V = Visible,
;; M = Visible and Mutated across -orig and -adjusted.
;; We are considering a context-length of 6 - example EX1:
;;
;;                 ... V V M V V M V ...
;;
;;           ...(? ? ? V V M)V V M V ...          #0
;;             ...(? ? V V M V)V M V ...          #1
;;               ...(? V V M V V)M V ...          #2
;;                 ...(V V M V V M)V ...          #3
;;                 ... V(V M V V M V)...          #4
;;                 ... V V(M V V M V ?)...        #5
;;                 ... V V M(V V M V ? ?)...      #6
;;                 ... V V M V(V M V ? ? ?)...    #7
;;                 ... V V M V V(M V ? ? ? ?)...  #8
;;
;; Here, for those windows that do not have `?`-s on them,
;; i.e. #3 and #4, we can simply read off probability rate-of-change
;; adjustments from the current multiverse world-probability,
;; but when we extend to the left and to the right, we need
;; to do this taking the full current-window context into account.
;;
;; There also are situations where we have to left-extend AND right-extend,
;; as in this example EX2:
;;
;;                       ... M V M ...
;;
;;            ....(? ? ? ? ? M)V M...              #0
;;              ....(? ? ? ? M V)M...              #1
;;                 ...(? ? ? M V M)...             #2
;;                   ...(? ? M V M ?)...           #3
;;                     ...(? M V M ? ?)...         #4
;;                       ...(M V M ? ? ?)...       #5
;;                        ...M(V M ? ? ? ?)...     #6
;;                        ...M V(M ? ? ? ? ?)...   #7
;;
;; General strategy:
;;
;; Broadly speaking, since the tape-content is fixed, we can do what
;; the call/cc mechanism does more compactly, avoiding some very costly
;; bookkeeping.
;;
;; - We represent a pair of orig/adjusted stretches with this data:
;;   <orig index into a flattened shape-[A]*k array>, <same for adjusted>, <k>
;;
;; - We also know the sequence-length cl-k.
;;
;; - We start with int-representations of the revealed orig and
;;   adjusted sequence, and the probability to be in this situation.
;;
;; - For left-shifting, there are two operations:
;;   - Increase window-length
;;   - Shift window
;; - Same for right-shifting.
;;
;; - We start with a potentially-short/potentially-long window.
;;
;;;; XXX-OBSOLETE-REST-OF-COMMENT-BELOW-XXX
;;
;; - We start recursion from the rightmost affected position,
;;   left-shifting the window in every step. We know how far
;;   to right-extend from that position.
;;
;; - We recursively-and-iteratively left-pad by drawing
;;   from the existing length-k and -(k+1) -marginal distributions.
;;
;; - At every left-padding step, we do two things:
;;   - recursively left-pad further.
;;   - also start recursively right-padding.
;;
;; - At the deepest recursion level, the next-step orig/adjusted-index align.
;;   - We also know the probability with which we recursed into this situation.
;;   - We need to update the rate-accumulator at the orig- and adjusted-index.

;; Being able to turn on debuging for accumulator-updates is so useful
;; that we provide code for this.
(define (report-accum-update-1 size-a cl-k
                               p-current
                               i-orig i-adjusted len-i
                               i-orig-here i-adjusted-here
                               accum probs)
  (let* ((vec-from-index
          (lambda (i-x target-len)
            (let ((result (make-vector target-len #f)))
              (let loop ((i 0))
                (if (= i target-len) result
                    (begin
                      (vector-set!
                       result (- target-len 1 i)
                       (remainder (quotient i-x (expt size-a i)) size-a))
                      (loop (+ 1 i))))))))
         (v-orig (vec-from-index i-orig len-i))
         (v-adjusted (vec-from-index i-adjusted len-i)))
    (display `(ddd-a- ,p-current at ,(vec-from-index i-orig-here cl-k) from
                      ,v-orig -> ,v-adjusted / ,i-orig-here ->
                      ,i-adjusted-here len-i ,len-i))
    (display "\n")
    (display `(ddd-a+ ,p-current at ,(vec-from-index i-adjusted-here cl-k) from
                      ,v-orig -> ,v-adjusted / ,i-orig-here ->
                      ,i-adjusted-here len-i ,len-i))
    (display "\n\n")
    ))


(define (lr-rec-extend-1 the-sp-info p-current i-orig i-adjusted len-i accum)
  (let* (;; Setting this to #t enables debugging accumulator-contributions.
         (debug-accum-contributions #f)
         (cl-k (vector-length (sp-table-striding (sp-info-table the-sp-info))))
         (size-a (vector-length (sp-info-symbols the-sp-info)))
         (prefix-modulus (expt size-a (- cl-k 1)))
         (window-modulus (* prefix-modulus size-a)))
    (letrec
      ((get-prob
        (lambda (index len-i)
          (fp-vector-ref
           (sp-table-probs (sp-info-marginals the-sp-info len-i))
           index)))
       ;;
       (get-prob-relative
        (lambda (index-short len-short index-long len-long)
          (let ((p-long (get-prob index-long len-long)))
            (if (= p-long 0.0) 0.0
                (/ p-long (max  ; (max ...) handles numerical noise.
                           p-long
                           (get-prob index-short len-short)))))))
       ;;
       (accumulate-dp/dt
        ;; Also returns whether there was work to be done.
        (lambda (p-current i-orig-maybe-extended
                           i-adjusted-maybe-extended)
          ;; (display `(ddd-accumulate-dp/dt
          ;;            p-current ,p-current
          ;;            i-orig-maybe-extended ,i-orig-maybe-extended
          ;;            i-adjusted-maybe-extended ,i-adjusted-maybe-extended))
          ;; (display "\n")
          (let ((i-orig-here (remainder i-orig-maybe-extended window-modulus))
                (i-adjusted-here
                 (remainder i-adjusted-maybe-extended window-modulus)))
            (if (= i-orig-here i-adjusted-here)
                #f  ; nothing to do
                (let ((d-p-orig/dt-prev (fp-vector-ref accum i-orig-here))
                      (d-p-adjusted/dt-prev
                       (fp-vector-ref accum i-adjusted-here)))
                  (fp-vector-set! accum i-orig-here
                                  (- d-p-orig/dt-prev p-current))
                  (fp-vector-set! accum i-adjusted-here
                                  (+ d-p-adjusted/dt-prev p-current))
                  (if debug-accum-contributions
                      (report-accum-update-1
                       size-a cl-k
                       p-current
                       i-orig i-adjusted len-i
                       i-orig-here i-adjusted-here
                       accum
                       (sp-table-probs (sp-info-marginals the-sp-info cl-k))
                       ))
                  )))))
       ;;
       (extend-ri-from-prefix
        (lambda (p-current i-orig-prefix i-adjusted-prefix)
          ;; (display `(ddd-extend-ri-from-prefix p-current ,p-current i-orig-prefix ,i-orig-prefix
          ;;                                      i-adjusted-prefix ,i-adjusted-prefix))
          ;; (display "\n")
          (if (= i-orig-prefix i-adjusted-prefix)
              #t  ; Same prefix, no change.
              (n-times (i size-a)
                (let* ((i-orig (+ (* i-orig-prefix size-a) i))
                       (i-adjusted (+ (* i-adjusted-prefix size-a) i))
                       (p-relative
                        (get-prob-relative i-orig-prefix (- cl-k 1)
                                           i-orig cl-k)))
                  (if (> p-relative 0.0)
                      (let ((p-new (* p-current p-relative)))
                        (accumulate-dp/dt p-new i-orig i-adjusted)
                        (extend-ri-from-prefix
                         p-new
                         (remainder i-orig prefix-modulus)
                         (remainder i-adjusted prefix-modulus)))))))))
       ;;
       (extend-le
        (lambda (p-current i-orig i-adjusted
                           len-i do-right-shifts-from-here)
          ;; (display `(ddd-extend-le p-current ,p-current i-orig ,i-orig i-adjusted ,i-adjusted
          ;;                          len-i ,len-i do-right-shifts-from-here ,do-right-shifts-from-here))
          ;; (display "\n")
          ;; If there is no change in need of processing anymore, we are done.
          (if (= i-orig i-adjusted)
              ;; In the current window, there is no change. We are done.
              #t
              ;; Otherwise, there is a change.
              ;; We first have to determine if we have a too-short window.
              ;; In that case, we need to left-extend the reading frame,
              ;; but keep its right-position.
              (begin
                (cond
                 ((< len-i cl-k)
                  ;; Try left-extending in all possible ways.
                  (n-times (i size-a)
                    (let* ((i-scaled (* i (expt size-a len-i)))
                           (i-orig-ext (+ i-scaled i-orig))
                           (i-adjusted-ext (+ i-scaled i-adjusted))
                           (len-i-next (+ 1 len-i))
                           (p-relative
                            (get-prob-relative
                             i-orig len-i i-orig-ext len-i-next)))
                      (if (> p-relative 0.0)
                          (extend-le (* p-current p-relative)
                                     i-orig-ext i-adjusted-ext
                                     len-i-next
                                     ;; As we are left-extending, once we reach
                                     ;; full prefix-length, we are at the point
                                     ;; where we do right-extending.
                                     (= len-i-next (- cl-k 1)))))))
                 ((= len-i cl-k)
                  ;; We have a just-full window - AND are guaranteed
                  ;; to have a change.
                  ;; First, perform the probability-adjustment right here.
                  ;; (Using the lowest cl-k many size-a-digits).
                  (accumulate-dp/dt p-current i-orig i-adjusted)
                  (n-times (i size-a)
                    (let* ((i-scaled (* i (expt size-a (- len-i 1))))
                           (i-suffix-orig (quotient i-orig size-a))
                           (i-suffix-adjusted (quotient i-adjusted size-a))
                           (i-orig-next (+ i-scaled i-suffix-orig))
                           (i-adjusted-next (+ i-scaled i-suffix-adjusted))
                           (p-relative
                            (get-prob-relative
                             i-suffix-orig (- len-i 1) i-orig-next len-i)))
                      (if (> p-relative 0.0)
                          (extend-le (* p-current p-relative)
                                     i-orig-next i-adjusted-next
                                     len-i
                                     ;; If we left-shifted a full frame,
                                     ;; we do not right shift from there.
                                     #f)))))
                 (#t  ;  (> len-i cl-k)
                  ;; We have extra digits to the left of the window.
                  ;; Accumulate and shift left.
                  (accumulate-dp/dt p-current i-orig i-adjusted)
                  (extend-le p-current
                             (quotient i-orig size-a)
                             (quotient i-adjusted size-a)
                             (- len-i 1)
                             ;; Never do right-extension from such a
                             ;; left-shifted window.
                             #f)))
                ;; We are done accumulating and possibly left-extending.
                ;; Let's see if we have to right-extend.
                (if do-right-shifts-from-here
                    (extend-ri-from-prefix
                     p-current
                     (remainder i-orig prefix-modulus)
                     (remainder i-adjusted prefix-modulus))))))))
    (extend-le p-current i-orig i-adjusted len-i
               ;; We can do right-extending right here if we have
               ;; at least a full prefix.
               (>= len-i (- cl-k 1))))))


;; If we are done with eval-fn evaluation, we still need to unfold tapes
;; so that every observed tape-change has full coverage, and then collect
;; all changes to get us the dy/dt.
;;
;; Since this function does not use any of the call/cc magic, and implements
;; the innermost loop, it would be a prime candidate for replacement with
;; a more performant c-lambda implementation.
;;
;; vec-buffer-ref: see dy/dt-collect-into-accumulator!
;; This mutable buffer's vectors will be replaced with a sufficiently-extended
;; vectors as-needed, and have the orig-sequence and adjusted-sequence unspooled
;; into them, for windowed processing.
(define (finish-fn-eval-fast-fixed
         mv-stack-ref the-sp-info exit-out-cont
         vec-buffer-ref dy/dt-accumulator)
  (let*
      ((the-mv-state (caar mv-stack-ref))
       (cl-k (vector-length (sp-table-striding (sp-info-table the-sp-info))))
       (size-a (vector-length (sp-info-symbols the-sp-info)))
       (p-current (mv-state-prob the-mv-state))
       (collect-dy/dt-contributions
        (lambda (data-tape?)
          (let* ((tvp (mv-state-tvp the-mv-state data-tape?))
                 (tv-orig (tv-pair-orig tvp))
                 (tv-adjusted (tv-pair-adjusted tvp))
                 (len-visible (+ (t-view-l-len tv-orig)
                                 (t-view-r-len tv-orig))))
            ;; XXX
            ;; (display `(ddd-calling-lr-rec-extend-1 tv-orig ,tv-orig tv-adjusted ,tv-adjusted))
            ;; (display "\n")
            (lr-rec-extend-1
             the-sp-info p-current
             (tv-as-index tv-orig size-a)
             (tv-as-index tv-adjusted size-a)
             len-visible
             dy/dt-accumulator)))))
    ;; (display "DDD collect p-tape\n")
    (collect-dy/dt-contributions #f)
    ;; (display "DDD collect d-tape\n")    
    (collect-dy/dt-contributions #t)))



;; "Reference" implementation uses splitting primitives.
;; (define finish-fn-eval finish-fn-eval-reference)
(define finish-fn-eval finish-fn-eval-fast-fixed)

;;;; END CODE FOR: faster inner loop

;; Computing dy/dt, basic variant.
;; Here, the tape-get and tape-set! functions directly work with
;; alphabet-indices, not symbols. Symbols might be more convenient
;; to work with in some cases.
;;
;; A key issue here is that we have to provide an exit-continuation
;; that gets invoked after an universe-splitting event and proceeds to
;; (mv-stack-process mv-stack-ref sp-info), which in general will
;; invoke continuations, and in this way come back to executing
;; the continuation of a fn-eval invocation.

(define (compute-dy/dt fn-eval symbols the-sp-table debug-flag)
  (let* ((cl-k (vector-length (sp-table-striding the-sp-table)))
         (sp-probs (sp-table-probs the-sp-table))
         (size-a (vector-length symbols))
         (num-params (fp-vector-length (sp-table-probs the-sp-table)))
         (the-sp-info (sp-info symbols the-sp-table))
         (dy/dt-accumulator (make-fp-vector num-params 0.0))
         (the-mv-stack-ref (cons '() debug-flag))
         (vec-buffer-ref (vector (make-u16-vector 32 0)
                                 (make-u16-vector 32 0) 0)))
    ;; (display "XXX-EXPERIMENTAL: Have nonzero limit on world-probability!\n")
    ;; We start with an empty tape-view.
    (let ((the-mv-state (make-new-mv-state size-a 0 1 0 0)))
      (set-car! the-mv-stack-ref (cons the-mv-state '()))
      (call-with-current-continuation
       (lambda (exit-out-cont)
         (let ((tape-get
                (lambda (data-tape? index)
                  (tape-get-1 the-mv-stack-ref the-sp-info
                              data-tape? index exit-out-cont)))
               (tape-set!
                (lambda (data-tape? index value)
                  (tape-set-1 the-mv-stack-ref the-sp-info
                              data-tape? index value exit-out-cont)))
               (vector-choose
                (lambda (v-probs v-options)
                  (choose-1 the-mv-stack-ref the-sp-info
                            v-probs v-options
                            exit-out-cont))))
           ;; TODO: Might want to protect the clean-up with a dynamic-wind.
           (fn-eval tape-get tape-set! vector-choose)
           (finish-fn-eval the-mv-stack-ref the-sp-info
                           exit-out-cont vec-buffer-ref dy/dt-accumulator)
           )))
      ;; this is exit-out-cont: 
      (mv-stack-process the-mv-stack-ref the-sp-info))
    dy/dt-accumulator))

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

;; The gambit-C based framework.

(include "gambit_macros.scm")
(include "tape_multiverse.scm")

;; The framework provides a global problem-registry;
;; this is necessary for having a simple Gambit/Python
;; interface and merely repackages other global information,
;; i.e. the symbol content of the compiled shared object.
(define _problem-registry (make-table))


(define-macro (register-problem problem-key symbols . body)
  (let ((g-symbols (gensym 'g-symbols))
        (g-key (gensym 'g-key)))
    `(let ((,g-symbols ,symbols)
           (,g-key ,problem-key))
       (table-set!
        _problem-registry
        ,g-key  ; used-once, so no need for a (let ((symbol ...)))
        (cons ,g-symbols
              (tape-evaluator ,g-symbols . ,body))))))


;; Integrating user-defined problems into the compiled shared object.
(include "problems.scm")

;; While we could wrap execution into (with-exception-catcher ...),
;; it likely is more useful to surface problems with user code
;; by dropping into a Scheme REPL.

(define populate-f64vector-from-c-data
  (c-lambda
   (scheme-object (pointer double) scheme-object)
   void
   ;; Note: replacing `target[i] = ...` with:
   ;; `___F64VECTORSET(___arg3, i, f64data[i])`;
   ;; fails, i.e. does not populate the vector as expected,
   ;; apparently due to misaligned memory striding.
   "
    int n = ___S32UNBOX(___arg1);
    double *f64data = (double *) ___arg2;
    double* target = &(___F64VECTORREF(___arg3, 0));
    for (int i = 0; i < n; ++i) {
      target[i] = f64data[i];
    }
    ___return;
    "))


(define f64vector-data-to-c
  (c-lambda
   (scheme-object scheme-object (pointer double))
    void
    "
    int n = ___S32UNBOX(___arg2);
    double *dst = (double *) ___arg3;
    double *src = &(___F64VECTORREF(___arg1, 0));
    for (int i = 0; i < n; ++i) {
      dst[i] = src[i];
    }
    ___return;
    "))


(define (c-compute-dy-dt problem-tag cl-k int-debug? c-probs-in c-probs-out)
  (let* ((symbols+fn-eval (table-ref _problem-registry problem-tag))
         (symbols (car symbols+fn-eval))
         (size-a (vector-length symbols))
         (fn-eval (cdr symbols+fn-eval))
         (debug? (not (= int-debug? 0)))
         (total-size (expt size-a cl-k))
         (probs (make-fp-vector total-size 0.0))
         (striding (make-vector-iter
                    cl-k
                    (lambda (k) (expt size-a (- cl-k 1 k)))))
         (the-sp-table (sp-table striding probs)))
    (populate-f64vector-from-c-data total-size c-probs-in probs)
    (let ((dy/dt (compute-dy/dt fn-eval symbols the-sp-table debug?)))
      (f64vector-data-to-c dy/dt total-size c-probs-out))))

;;;; === C Interface ===

;; For debugging, so we can use printf() inside (c-lambda ...) code.
(c-declare "#include <stdio.h>")  

(c-define (c_register_problems n) (int64) int64 "c_register_problems" "extern"
   (begin
     ;; Registering the 'canary' problem that gets exercised at python
     ;; module import time to validate that the Scheme/Python interface works.
     (register-problem
      "__canary_problem_radioactive_decay"
      #(A B)
      (if (eq? (tape-get-sym #t 0) 'B)
          (tape-set-sym! #t 0 'A)))
     ;; Registering the problems from `problems.scm`.
     (register-problems)
     (+ 1 n)))


(c-define (c_compute_dy_dt problem-tag cl-k int-debug? c-probs-in c-probs-out)
          (char-string int64 int64 (pointer double) (pointer double)) void
          "c_compute_dy_dt" "extern"
  (c-compute-dy-dt
   problem-tag
   ;; (+ 0 ...) for type conversion.
   (+ 0 cl-k) (+ 0 int-debug?)
   c-probs-in c-probs-out))

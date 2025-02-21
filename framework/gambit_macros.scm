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


;; Gambit uses its own (define-macro ...), which looks a lot like
;; Common Lisp's DEFMACRO, but does not support define-syntax
;; out of the box.
;;
;; So, we provide alternative macro definitions here.

(define-macro (n-times sym-and-end . body)
  (let ((sym-n (car sym-and-end))
        (end (cadr sym-and-end))
        (gs-end (gensym 'gs-end))
        (gs-loop (gensym 'gs-loop)))
    `(let ((,gs-end ,end))
       (let ,gs-loop ((,sym-n 0))
            (if (< ,sym-n ,gs-end)
                (begin
                  (begin . ,body)
                  (,gs-loop (+ 1 ,sym-n)))
                #t)))))

(define-macro (check-tag sym x)
  (let ((gs-x (gensym 'gs-x)))
    `(let ((,gs-x ,x))
       (if (not (and (vector? ,gs-x)
                     (>= (vector-length ,gs-x) 1)
                     (eq? (vector-ref ,gs-x 0) ,sym)))
           (error "Expected tag:" ,sym ,gs-x)
           #t))))

(define-macro (defobj name-tag-x-rest . body)
  (let ((name (car name-tag-x-rest))
        (tag (cadr name-tag-x-rest))
        (x (caddr name-tag-x-rest))
        (rest (cdddr name-tag-x-rest)))
    `(define (,name ,x . ,rest)
       (check-tag ,tag ,x) . ,body)))


;;;; === Specialized Vector Representations ===
;;
;; For Gambit, we use the specialised f64vector vectors.
(define make-fp-vector make-f64vector)
(define fp-vector f64vector)
(define fp-vector-length f64vector-length)
(define fp-vector-ref f64vector-ref)
(define fp-vector-set! f64vector-set!)

(define make-u16-vector make-u16vector)
(define u16-vector u16vector)
(define u16-vector-length u16vector-length)
(define u16-vector-ref u16vector-ref)
(define u16-vector-set! u16vector-set!)


;; Wrapper for fn-eval definition that provides tape-get-sym / tape-set-sym!
;; and also tape-get / tape-set!
;; We can then do:
;; (define fn-eval (tape-evaluator #(A B) ...body...))

;; We need helpers to get vector-choose args from the more convenient choose args.
(define (probs-from-weights-and-options-1 n list-weight+option)
  (let ((result (make-vector n))
        (total-weight
         (let loop ((total 0) (todo list-weight+option))
           (if (null? todo)
               total
               (loop (+ total (caar todo)) (cdr todo))))))
    (let loop ((k 0) (todo list-weight+option))
      (if (null? todo)
          result
          (begin
            (vector-set! result k (/ (caar todo) total-weight))
            (loop (+ 1 k) (cdr todo)))))))


(define (options-from-weights-and-options-1 n list-weight+option)
  (let ((result (make-vector n)))
    (let loop ((k 0) (todo list-weight+option))
      (if (= k n) result
          (begin
            (vector-set! result k (cadar todo))
            (loop (+ 1 k) (cdr todo)))))))

    
(define-macro (tape-evaluator symbols . body)
  (let ((g-symbols (gensym 'g-symbols))
        (g-index-by-symbol (gensym 'g-index-by-symbol))
        (g-loop (gensym 'loop))
        (g-k (gensym 'k))
        )
    `(let* ((,g-symbols ,symbols)
            (,g-index-by-symbol (make-table)))
       ;; Populating the index-by-symbol table.
       (let ,g-loop ((,g-k (- (vector-length ,g-symbols) 1)))
            (table-set! ,g-index-by-symbol (vector-ref ,g-symbols ,g-k) ,g-k)
            (if (>= ,g-k 1) (,g-loop (- ,g-k 1))))
       (lambda (tape-get tape-set! vector-choose)
         ;; No risk of any symbol clash here.
         (let ((tape-get-sym
                (lambda (data-tape? index)
                  (vector-ref ,g-symbols (tape-get data-tape? index))))
               (tape-set-sym!
                (lambda (data-tape? index sym)
                  (tape-set! data-tape? index (table-ref ,g-index-by-symbol sym))))
               (choose
                (lambda (list-weight+option)
                  (let ((n (length list-weight+option)))
                    (vector-choose
                     (probs-from-weights-and-options-1 n list-weight+option)
                     (options-from-weights-and-options-1 n list-weight+option))))))
           . ,body)))))

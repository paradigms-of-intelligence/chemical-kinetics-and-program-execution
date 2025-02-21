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

;; Loaded by framework.scm and compiled alongside other code into a
;; .so to be dlopen()ed by Python.

(define (register-problems)
  ;;;;
  ;;;; Example 1: "Radioactive decay"
  ;;;;
  (register-problem
   "ex1-radioactive-decay"
   #(A B)
   (if (eq? (tape-get-sym #t 0) 'B)
       (tape-set-sym! #t 0 'A)))
  ;;;;
  ;;;; Example 2: "Ferromagnetic chain"
  ;;;;
  (let ((param-J 1.0)
        (param-h -0.25)
        (beta 1.0)
        )
    (register-problem
     "ex2-ferromagnetic-chain"
     #(D U)
     (let ((p-mid (tape-get-sym #t 0))
           (p-left (tape-get-sym #t -1))
           (p-right (tape-get-sym #t +1)))
       (let* ((energy-J (+ (if (eq? p-left p-mid) 1 -1)
                           (if (eq? p-mid p-right) 1 -1)))
              (factor-a (exp (- (* beta param-J (+ 4 (* 2 energy-J))))))
              (factor-b
               (if (eqv? (> param-h 0) (eq? p-mid 'U))
                   ;; Either field is pushing up, and we are in an
                   ;; up-configuration, or field is pushing down,
                   ;; and we are in a down-configuration.
                   ;; In both cases, flip-rate has to be suppressed.
                   (exp (- (* 2 beta (abs param-h))))
                   1.0))
              (p-flip (* factor-a factor-b))
              (p-stay (- 1 p-flip)))
         (if (choose `((,p-flip #t) (,p-stay #f)))
             (tape-set-sym! #t 0 (if (eq? p-mid 'U) 'D 'U))
             #t)))))
  ;;;;  
  ;;;; Example 3: "Nylon copolymerization"
  ;;;;
  (register-problem
   "ex3-copolymerization"
   #(O A M N)
   (let ((p0 (tape-get-sym #f 0)))
     (if (and (not (eq? p0 'O))
              (eq? (tape-get-sym #f -1) 'O)
              (eq? (tape-get-sym #f +1) 'O))
         ;; We have an isolated monomer on the P-tape.
         (let ((d0 (tape-get-sym #t 0)))
           (if (or (and (eq? p0 'A) (or (eq? d0 'M) (eq? d0 'N)))
                   (and (eq? d0 'A) (or (eq? p0 'M) (eq? p0 'N))))
               ;; We have compatible monomers.
               ;; Need to see if we have an available end to connect to.
               ;; We try both sides with equal probability.
               ;; If there is no open end on one side, we give up.
               (let* ((i (choose '((1.0 -1) (1.0 +1))))
                      (di (tape-get-sym #t i)))
                 (if (and (eq? di 'O)
                          ;; post-addition, we also have an endpoint,
                          ;; i.e. this is not an accidental place where
                          ;; two chain-ends meet.
                          (eq? (tape-get-sym #t (* 2 i)) 'O))
                     (begin
                       ;; Remove the monomer from the P-tape.
                       (tape-set-sym! #f 0 'O)
                       ;; Connect the monomer on the D-tape.
                       (tape-set-sym! #t i p0)))))))))
  ;;
  ;; Variant 1: Preference for alternation.
  ;;
  (register-problem
   "ex3var1-copolymerization"
   #(O A M N)
   (let ((p0 (tape-get-sym #f 0)))
     (if (and (not (eq? p0 'O))
              (eq? (tape-get-sym #f -1) 'O)
              (eq? (tape-get-sym #f +1) 'O))
         ;; We have an isolated monomer on the P-tape.
         (let ((d0 (tape-get-sym #t 0)))
           (if (or (and (eq? p0 'A) (or (eq? d0 'M) (eq? d0 'N)))
                   (and (eq? d0 'A) (or (eq? p0 'M) (eq? p0 'N))))
               ;; We have compatible monomers.
               ;; Need to see if we have an available end to connect to.
               ;; We try both sides with equal probability.
               ;; If there is no open end on one side, we give up.
               (let* ((i (choose '((1.0 -1) (1.0 +1))))
                      (di (tape-get-sym #t i)))
                 (if (and (eq? di 'O)
                          ;; post-addition, we also have an endpoint,
                          ;; i.e. this is not an accidental place where
                          ;; two chain-ends meet.
                          (eq? (tape-get-sym #t (* 2 i)) 'O))
                     ;; if the opposite-side of the D-tape unit is of the same
                     ;; type as the new unit, prevent the reaction with some
                     ;; probability.
                     (if (and
                          ;; New unit is M or N.
                          (not (eq? p0 'A))
                          ;; Same unit already on other side of A.
                          (eq? (tape-get-sym #t (- i)) p0)
                          ;; High rejection rate for this case. - XXX transfer to paper
                          (choose '((75.0 #t) (25.0 #f))))
                         #f  ; do nothing
                         (begin  ; otherwise, react as before.
                           ;; Remove the monomer from the P-tape.
                           (tape-set-sym! #f 0 'O)
                           ;; Connect the monomer on the D-tape.
                           (tape-set-sym! #t i p0))))))))))
  ;;
  ;; Variant 2: Reversible reactions.
  ;;
  (register-problem
   ;; Sequences `O[AMN]?...` and `...[AMN]?O` plus `OOO` can depolymerize,
   ;; but with low relative rate.
   ;; Rate-ratio {polymerization} : {depolymerization} is related to
   ;; the free enthalpy of the reaction, i.e. thermodynamic stability
   ;; of the polymer relative to monomers.
   "ex3var2-copolymerization"
   #(O A M N)
   (let ((p0 (tape-get-sym #f 0)))
     (if
      (eq? p0 'O)
      ;; "program-tape" cell 0 is empty, try dissociation.
      (if (and (eq? (tape-get-sym #f -1) 'O)
               (eq? (tape-get-sym #f +1) 'O))
          ;; We have free space on the P-tape.
          (let ((d0 (tape-get-sym #t 0)))
            (if (not (eq? d0 'O))
                (let ((d1-right (tape-get-sym #t 1))
                      (d1-left (tape-get-sym #t -1)))
                  (if (= 1 (+ (if (eq? d1-left 'O) 0 1)
                              (if (eq? d1-right 'O) 0 1)))
                      ;; We are at the end of a chain.
                      ;; Depolymerization happens at a reduced rate,
                      ;; since we take the polymer to be thermodynamically
                      ;; more stable than the monomers.
                      (if (choose '((1.0 #t) (50.0 #f)))
                          (begin
                            (tape-set-sym! #f 0 d0)
                            (tape-set-sym! #t 0 'O))))))))
      ;; else, "program-tape" cell 0 is not empty, try polycondensation.
      (if (and (eq? (tape-get-sym #f -1) 'O)
               (eq? (tape-get-sym #f +1) 'O))
          ;; We have an isolated monomer on the P-tape.
          (let ((d0 (tape-get-sym #t 0)))
            (if (or (and (eq? p0 'A) (or (eq? d0 'M) (eq? d0 'N)))
                    (and (eq? d0 'A) (or (eq? p0 'M) (eq? p0 'N))))
                ;; We have compatible monomers.
                ;; Need to see if we have an available end to connect to.
                ;; We try both sides with equal probability.
                ;; If there is no open end on one side, we give up.
                (let* ((i (choose '((1.0 -1) (1.0 +1))))
                       (di (tape-get-sym #t i)))
                  (if (and (eq? di 'O)
                           ;; post-addition, we also have an endpoint,
                           ;; i.e. this is not an accidental place where
                           ;; two chain-ends meet.
                           (eq? (tape-get-sym #t (* 2 i)) 'O))
                      (begin
                        ;; Remove the monomer from the P-tape.
                        (tape-set-sym! #f 0 'O)
                        ;; Connect the monomer on the D-tape.
                        (tape-set-sym! #t i p0))))))))))
  ;;;;  
  ;;;; Example 4: "Chemical Turing Machine"
  ;;;;
  ;; We can move some definitions out of problem-registration.
  (let* ((is-io? (lambda (x) (or (eq? x 'I) (eq? x 'O))))
         (px-relative-stability-reverse-suppression-factor 0.05)
         (px-reverse-suppression-choices
          `((,(- 1.0 px-relative-stability-reverse-suppression-factor) #f)
            (,px-relative-stability-reverse-suppression-factor #t))))
    (register-problem
     "ex4-chemical-turing"
     #(A B C D I O P X S)  ; S = Solvent, P = Powered, X = De-Powered.
     (let ((p0 (tape-get-sym #f 0)))
       (cond
        ((and (eq? p0 'P)  ; powered->de-powered
              ;; We need to suppress this by a factor 2, since otherwise,
              ;; back- and forward- reaction constants would not be the same at
              ;; px-relative-stability-reverse-suppression-factor = 0.
              (choose '((1.0 #t) (1.0 #f)))
              )
         (let ((d0 (tape-get-sym #t 0)))
           (cond
            ((and (eq? d0 'A)
                  ;; Can we advance?
                  (is-io? (tape-get-sym #t 1))
                  ;; Post-advancement, we again have to be in a valid state
                  ;; where the cursor-marker is just before an I or O.
                  (is-io? (tape-get-sym #t 2))                  
                  )
             (tape-set-sym! #f 0 'X)             
             (tape-set-sym! #t 0 'I)
             (tape-set-sym! #t 1 'B))
            ((and (eq? d0 'B)
                  (is-io? (tape-get-sym #t 1))
                  (is-io? (tape-get-sym #t 2)))
             (tape-set-sym! #f 0 'X)                
             (tape-set-sym! #t 0 'O)
             (tape-set-sym! #t 1 'C))
            ((and (eq? d0 'C)
                  (is-io? (tape-get-sym #t 1))
                  (is-io? (tape-get-sym #t 2)))
             (tape-set-sym! #f 0 'X)  
             (tape-set-sym! #t 0 'I)
             (tape-set-sym! #t 1 'D)))))
        ((eq? p0 'X)  ; de-powered->powered
         (let ((d0 (tape-get-sym #t 0)))
           (if (and (or (eq? d0 'B) (eq? d0 'C) (eq? d0 'D))
                    (is-io? (tape-get-sym #t -1))  ; Can move back
                    (is-io? (tape-get-sym #t -2))  ; Won't move next to a cursor.
                    ;; Also, the previous symbol needs to be compatible with
                    ;; the forward-reaction end-state.
                    (or (and (eq? d0 'C) (eq? (tape-get-sym #t -1) 'O))
                        (and (not (eq? d0 'C)) (eq? (tape-get-sym #t -1) 'I)))
                    ;; If P is thermodynamically more stable than X, we need to
                    ;; further suppress this reaction.
                    (choose px-reverse-suppression-choices))
               (begin
                 (tape-set-sym! #f 0 'P)
                 (tape-set-sym! #t 0 (choose '((1.0 I) (1.0 O))))
                 (tape-set-sym! #t -1
                                (cond ((eq? d0 'B) 'A)
                                      ((eq? d0 'C) 'B)
                                      ((eq? d0 'D) 'C)))))))
        ))))
  ;; Variant 1: Equal thermodynamic stability for forward- and backward-reaction.
  (let* ((is-io? (lambda (x) (or (eq? x 'I) (eq? x 'O))))
         (px-relative-chemical-stability-reverse-suppression-factor 0.0)  ; no suppression
         (px-reverse-suppression-choices
          `((,(- 1.0 px-relative-chemical-stability-reverse-suppression-factor) #t)
            (,px-relative-chemical-stability-reverse-suppression-factor #f))))
    (register-problem
     "ex4var1-chemical-turing"
     #(A B C D I O P X S)  ; S = Solvent, P = Powered, X = De-Powered.
     (let ((p0 (tape-get-sym #f 0)))
       (cond
        ((and (eq? p0 'P)  ; powered->de-powered
              ;; We need to suppress this by a factor 2, since otherwise,
              ;; back- and forward- reaction constants would not be the same at
              ;; relative-chemical-stability-factor-px = 1.
              (choose '((1.0 #t) (1.0 #f)))
              )
         (let ((d0 (tape-get-sym #t 0)))
           (cond
            ((and (eq? d0 'A)
                  ;; Can we advance?
                  (is-io? (tape-get-sym #t 1))
                  ;; Post-advancement, we again have to be in a valid state
                  ;; where the cursor-marker is just before an I or O.
                  (is-io? (tape-get-sym #t 2))                  
                  )
             (tape-set-sym! #f 0 'X)
             (tape-set-sym! #t 0 'I)
             (tape-set-sym! #t 1 'B))
            ((and (eq? d0 'B)
                  (is-io? (tape-get-sym #t 1))
                  (is-io? (tape-get-sym #t 2)))
             (tape-set-sym! #f 0 'X)                
             (tape-set-sym! #t 0 'O)
             (tape-set-sym! #t 1 'C))
            ((and (eq? d0 'C)
                  (is-io? (tape-get-sym #t 1))
                  (is-io? (tape-get-sym #t 2)))
             (tape-set-sym! #f 0 'X)  
             (tape-set-sym! #t 0 'I)
             (tape-set-sym! #t 1 'D)))))
        ((eq? p0 'X)  ; de-powered->powered
         (let ((d0 (tape-get-sym #t 0)))
           (if (and (or (eq? d0 'B) (eq? d0 'C) (eq? d0 'D))
                    (is-io? (tape-get-sym #t -1))  ; Can move back
                    (is-io? (tape-get-sym #t -2))  ; Won't move next to a cursor.
                    ;; Also, the previous symbol needs to be compatible with the forward-reaction end-state.
                    (or (and (eq? d0 'C) (eq? (tape-get-sym #t -1) 'O))
                        (and (not (eq? d0 'C)) (eq? (tape-get-sym #t -1) 'I)))
                    ;; If P is thermodynamically more stable than X, we need to further suppress
                    ;; this reaction.
                    (choose px-reverse-suppression-choices))
               (begin
                 (tape-set-sym! #f 0 'P)
                 (tape-set-sym! #t 0 (choose '((1.0 I) (1.0 O))))
                 (tape-set-sym! #t -1
                                (cond ((eq? d0 'B) 'A)
                                      ((eq? d0 'C) 'B)
                                      ((eq? d0 'D) 'C)))))))
        ))))
  ;; Variant 2: Detachable Evaluator, Symbol-Under-Cursor.
  (let* ((is-io? (lambda (x) (or (eq? x 'I) (eq? x 'O))))
         (choice-IO '((1.0 I) (1.0 O)))
         (choice-1:1 '((1.0 #t) (1.0 #f)))
         (beta 1.0)  ; Adjustable 1/(k_B T) factor.
         ;; Free enthalpies of formation. A large G-E disfavors
         ;; evaluators in solution.
         (G-P 6.0) (G-X 0.0) (G-E 1.0)
         ;; Final D-state must be thermodynamically less stable
         ;; if we want D to detach more easily to E than A.
         (G-A -1.0) (G-B -1.0) (G-C -1.0) (G-D 1.5)
         ;; With these parameters, the fastest reactions are the A+P->B+X
         ;; type reactions. 
         (Delta-G-fastest (- (+ G-B G-X) (+ G-A G-P)))
         (get-rate-factor
          (lambda (G-left G-right)
            (let ((rate-factor
                   (exp (- (* beta (- G-right G-left Delta-G-fastest))))))
              (if (> rate-factor 1.001)
                  (error
                   "Setup error: Delta-G-fastest not actually fastest.")
                  (min 1.0 rate-factor)))))
         (rate-choices
          (lambda (G-left G-right)
            (let ((r (get-rate-factor G-left G-right)))
              `((,r #t) (,(- 1 r) #f)))))
         (rate-choices-A+P->B+X (rate-choices (+ G-A G-P) (+ G-B G-X)))
         (rate-choices-B+X->A+P (rate-choices (+ G-B G-X) (+ G-A G-P)))
         (rate-choices-B+P->C+X (rate-choices (+ G-B G-P) (+ G-C G-X)))
         (rate-choices-C+X->B+P (rate-choices (+ G-C G-X) (+ G-B G-P)))
         (rate-choices-C+P->D+X (rate-choices (+ G-C G-P) (+ G-D G-X)))
         (rate-choices-D+X->C+P (rate-choices (+ G-D G-X) (+ G-C G-P)))
         (rate-choices-A->E (rate-choices G-A G-E))
         (rate-choices-D->E (rate-choices G-D G-E))                  
         (rate-choices-E->A+D
          (let ((r-A (get-rate-factor G-E G-A))
                (r-D (get-rate-factor G-E G-D)))
            (if (> (+ r-A r-D) 1.0)
                ;; In order to handle this case, we would have to set
                ;; Delta-G-fastest to make this rate fastest.
                (error "E->A+D rates too high to merge, given Delta-G-fastest.")
                `((,r-A A) (,r-D D) (,(- 1.0 r-A r-D) #f))))))
    ;; It can be useful to show rates at problem registration time,
    ;; for visual inspection. This can be done as follows:
    ;;(begin
    ;;  (display `(DEBUG rates
    ;;                   rate-choices-A+P->B+X ,rate-choices-A+P->B+X
    ;;                   rate-choices-B+P->C+X ,rate-choices-B+P->C+X               
    ;;                   rate-choices-C+P->D+X ,rate-choices-C+P->D+X
    ;;                   rate-choices-D+X->C+P ,rate-choices-D+X->C+P
    ;;                   rate-choices-C+X->B+P ,rate-choices-C+X->B+P               
    ;;                   rate-choices-B+X->A+P ,rate-choices-B+X->A+P
    ;;                   rate-choices-A->E ,rate-choices-A->E
    ;;                   rate-choices-D->E ,rate-choices-D->E
    ;;                   rate-choices-E->A+D ,rate-choices-E->A+D))
    ;;  (display "\n"))
    (register-problem
     "ex4var2-chemical-turing"
     ;; S = Solvent, P = Powered, X = De-Powered, E = Detached Evaluator.     
     #(A B C D I O P X S E) 
     (let ((p0 (tape-get-sym #f 0)))
       (cond
        ((and (eq? p0 'P)  ; powered->de-powered
              ;; Data tape is "?[IO][IO]" - so, if "?" is a cursor,
              ;; we can advance to a valid state.
              (is-io? (tape-get-sym #t 1))
              (is-io? (tape-get-sym #t 2))
              ;; We need to suppress this by another factor 2, since
              ;; back-reaction is two different reactions, depending on
              ;; what bit gets written.
              (choose choice-1:1))
         (let ((d0 (tape-get-sym #t 0)))
           (cond
            ((and (eq? d0 'A) (choose rate-choices-A+P->B+X))
             (tape-set-sym! #f 0 'X)
             (tape-set-sym! #t 0 'I)
             (tape-set-sym! #t 1 'B))
            ((and (eq? d0 'B) (choose rate-choices-B+P->C+X))
             (tape-set-sym! #f 0 'X)
             (tape-set-sym! #t 0 'O)
             (tape-set-sym! #t 1 'C))
            ((and (eq? d0 'C) (choose rate-choices-C+P->D+X))
             (tape-set-sym! #f 0 'X)
             (tape-set-sym! #t 0 'I)
             (tape-set-sym! #t 1 'D)))))
        ((and (eq? p0 'X)  ; de-powered->powered
              ;; Data tape is "[IO][IO]?" - so, if "?" is a cursor,
              ;; we can un-advance to a valid state.
              (is-io? (tape-get-sym #t -1))
              (is-io? (tape-get-sym #t -2)))
         (let ((d0 (tape-get-sym #t 0)))
           (cond
            ((and (eq? d0 'B) (choose rate-choices-B+X->A+P))
             (tape-set-sym! #f 0 'P)
             (tape-set-sym! #t 0 (choose choice-IO))
             (tape-set-sym! #t -1 'A))
            ((and (eq? d0 'C) (choose rate-choices-C+X->B+P))
             (tape-set-sym! #f 0 'P)
             (tape-set-sym! #t 0 (choose choice-IO))
             (tape-set-sym! #t -1 'B))
            ((and (eq? d0 'D) (choose rate-choices-D+X->C+P))
             (tape-set-sym! #f 0 'P)
             (tape-set-sym! #t 0 (choose choice-IO))
             (tape-set-sym! #t -1 'C)))))
        ((and (eq? p0 'E)  ; Detached evaluator that can attach.
              (is-io? (tape-get-sym #t 0))
              (is-io? (tape-get-sym #t +1))
              (is-io? (tape-get-sym #t -1))
              (choose choice-1:1))  ; We overwrite one bit.
         (let ((A-D-f (choose rate-choices-E->A+D)))
           (cond
            ((eq? A-D-f 'A)
             (tape-set-sym! #f 0 'S)
             (tape-set-sym! #t 0 'A))
            ((eq? A-D-f 'D)
             (tape-set-sym! #f 0 'S)
             (tape-set-sym! #t 0 'D)))))
        ((and (eq? p0 'S)
              (is-io? (tape-get-sym #t +1))
              (is-io? (tape-get-sym #t -1)))
         (let ((d0 (tape-get-sym #t 0)))
           (cond
            ((and (eq? d0 'A) (choose rate-choices-A->E))
             (tape-set-sym! #f 0 'E)
             (tape-set-sym! #t 0 (choose choice-IO)))
            ((and (eq? d0 'D) (choose rate-choices-D->E))
             (tape-set-sym! #f 0 'E)
             (tape-set-sym! #t 0 (choose choice-IO)))
            )))))))
  ;;;; =================
  ;;;;  
  ;;;; Example 5: "Simple Machine Language (with guaranteed termination)"
  ;;;;
  (let ((single-R-can-execute #f))
    (register-problem
     "ex5-msrtf-machine"
     #(M S R T F)
     (let loop ((Q 4) (Is 0) (Ip 0) (Id 0) (Op #f) (NT 0) (NR 0) (NF 0))
       (let ((op-todo (if (> Q 0) (tape-get-sym #f Ip) Op)))
         (if (= Q 4)
             (cond
              ((eq? op-todo 'S)
               (loop (- Q 1) Is (+ 1 Ip) Id op-todo 0 0 0))
              ((and (eq? op-todo 'R) single-R-can-execute)
               (tape-set! #t Id (modulo (+ 1 (tape-get #t Id)) 5))))
             ;; Otherwise, not-first-op.
             (case op-todo
               ((T)
                (let ((activated? (and (> NT 0) (> NF 0))))
                  (if activated? (tape-set! #t Id (tape-get #f Is)))
                  (if (not (or (= Q 1) (= Q -3)))
                      (loop (- Q 1)
                            (if activated? (+ 1 Is) Is)
                            (if (> Q 0) (+ 1 Ip) Ip)
                            (if activated? (+ 1 Id) Id)
                            op-todo
                            1 NR NF))))
               ((R)
                (if (> NR 0)
                    (tape-set! #t Id (modulo (+ 1 (tape-get #t Id)) 5)))
                (if (not (or (= Q 1) (= Q -3)))
                    (loop (- Q 1) Is
                          (if (> Q 0) (+ 1 Ip) Ip)
                          Id
                          op-todo
                          NT 1 NF)))
               ((F)
                (if (not (or (= Q 1) (= Q -3)))
                    (loop (- Q 1) Is 
                          (if (> Q 0) (+ 1 Ip) Ip)
                          Id
                          op-todo
                          NT NR 1)))
               ((M)
                (if (or (eq? Op 'R) (eq? Op 'T))
                    (loop -1 Is Ip Id Op NT NR NF)))))))))
  ;;
  ;; Variant 1: A single R-operation can execute independently.
  ;;
  (let ((single-R-can-execute #t))
    (register-problem
     "ex5var1-msrtf-machine"
     #(M S R T F)
     (let loop ((Q 4) (Is 0) (Ip 0) (Id 0) (Op #f) (NT 0) (NR 0) (NF 0))
       (let ((op-todo (if (> Q 0) (tape-get-sym #f Ip) Op)))
         (if (= Q 4)
             (cond
              ((eq? op-todo 'S)
               (loop (- Q 1) Is (+ 1 Ip) Id op-todo 0 0 0))
              ((and (eq? op-todo 'R) single-R-can-execute)
               (tape-set! #t Id (modulo (+ 1 (tape-get #t Id)) 5))))
             ;; Otherwise, not-first-op.
             (case op-todo
               ((T)
                (let ((activated? (and (> NT 0) (> NF 0))))
                  (if activated? (tape-set! #t Id (tape-get #f Is)))
                  (if (not (or (= Q 1) (= Q -3)))
                      (loop (- Q 1)
                            (if activated? (+ 1 Is) Is)
                            (if (> Q 0) (+ 1 Ip) Ip)
                            (if activated? (+ 1 Id) Id)
                            op-todo
                            1 NR NF))))
               ((R)
                (if (> NR 0)
                    (tape-set! #t Id (modulo (+ 1 (tape-get #t Id)) 5)))
                (if (not (or (= Q 1) (= Q -3)))
                    (loop (- Q 1) Is
                          (if (> Q 0) (+ 1 Ip) Ip)
                          Id
                          op-todo
                          NT 1 NF)))
               ((F)
                (if (not (or (= Q 1) (= Q -3)))
                    (loop (- Q 1) Is 
                          (if (> Q 0) (+ 1 Ip) Ip)
                          Id
                          op-todo
                          NT NR 1)))
               ((M)
                (if (or (eq? Op 'R) (eq? Op 'T))
                    (loop -1 Is Ip Id Op NT NR NF)))))))))
  ;;;;
  ;;;; Example 6: "Mini-BFF"
  ;;;;
  (let ((alphabet
         #(sym< sym> sym-cl sym-cr  ; cl/cr = curly left/right {}
           sym- sym+ sym-dot sym-comma
           sym-bl sym-br sym0 sym-nop)  ; bl/br = bracket left/right []
         ))
    (register-problem
     "ex6-mini-bff"
     alphabet
     (let loop ((max-num-syms-to-still-read 10)
                (p-offset 0)
                (d0-offset 0)  ; "head 0" offset
                (d1-offset 12)  ; "head 1" offset
                ;; Scan-mode N<0: Look for (-N)-th [-bracket on the left.
                ;; Scan mode N>0: Look for N-th ]-bracket on the right.
                ;; Scan mode 0: Execute local command.
                (scan-mode 0))
       (if (= max-num-syms-to-still-read 0)
           #f  ; Done.
           (let ((op (tape-get-sym #f p-offset)))
             (cond
              ((< scan-mode 0)
               (cond
                ((eq? op 'sym-bl)
                 (if (= scan-mode -1)
                     (loop (- max-num-syms-to-still-read 1)
                           (+ p-offset 1)  ; right after [
                           d0-offset d1-offset 0)
                     (loop (- max-num-syms-to-still-read 1)
                           (+ p-offset -1)
                           d0-offset d1-offset (+ scan-mode 1))))
                ((eq? op 'sym-br)
                 (loop (- max-num-syms-to-still-read 1)
                       (+ p-offset -1)
                       d0-offset d1-offset (+ scan-mode -1)))
                (#t
                 (loop (- max-num-syms-to-still-read 1)
                       (+ p-offset -1) d0-offset d1-offset scan-mode))))
              ((> scan-mode 0)
               (cond
                ((eq? op 'sym-br)
                 (if (= scan-mode 1)
                     (loop (- max-num-syms-to-still-read 1)
                           (+ p-offset 1)  ; right after ]
                           d0-offset d1-offset 0)
                     (loop (- max-num-syms-to-still-read 1)
                           (+ p-offset 1)
                           d0-offset d1-offset (+ scan-mode -1))))
                ((eq? op 'sym-bl)
                 (loop (- max-num-syms-to-still-read 1)
                       (+ p-offset 1)
                       d0-offset d1-offset (+ scan-mode 1)))
                (#t
                 (loop (- max-num-syms-to-still-read 1)
                       (+ p-offset 1) d0-offset d1-offset scan-mode))))
              (#t
               (cond
                ((or (eq? op sym<) (eq? op sym>))
                 (loop (- max-num-syms-to-still-read 1)
                       (+ p-offset 1)
                       (+ d0-offset (if (eq? op sym<) -1 +1))
                       d1-offset 0))
                ((or (eq? op sym-cl) (eq? op sym-cr))
                 (loop (- max-num-syms-to-still-read 1)
                       (+ p-offset 1)
                       d0-offset
                       (+ d1-offset (if (eq? op sym<) -1 +1))
                       0))
                ((or (eq? op sym+) (eq? op sym-))
                 (tape-set! #t d0-offset
                            (modulo
                             (+ (tape-get #t d0-offset 1) (if (eq? op sym+) +1 -1))
                             (vector-length alphabet)))
                 (loop (- max-num-syms-to-still-read 1)
                       (+ p-offset 1) d0-offset d1-offset 0))
                ((eq? op sym-dot)
                 (tape-set! #t d1-offset (tape-get #t d0-offset))
                 (loop (- max-num-syms-to-still-read 1)
                       (+ p-offset 1) d0-offset d1-offset 0))
                ((eq? op sym-comma)
                 (tape-set! #t d0-offset (tape-get #t d1-offset))
                 (loop (- max-num-syms-to-still-read 1)
                       (+ p-offset 1) d0-offset d1-offset 0))
                ((eq? op sym-bl)
                 (loop (- max-num-syms-to-still-read 1)
                       (+ p-offset 1) d0-offset d1-offset
                       (if (eq? sym0 (tape-get-sym #t d0-offset))
                           ;; Either enter scan-mode or make this a no-op.
                           +1 0)))
                ((eq? op sym-br)
                 (if (eq? sym0 (tape-get-sym #t d0-offset))
                     ;; no-op.
                     (loop (- max-num-syms-to-still-read 1)
                           (+ p-offset 1) d0-offset d1-offset 0)
                     ;; otherwise, scan backwards.
                     (loop (- max-num-syms-to-still-read 1)
                           (+ p-offset -1) d0-offset d1-offset -1)))
                (#t  ; no-op.
                 (loop (- max-num-syms-to-still-read 1)
                       (+ p-offset 1) d0-offset d1-offset 0))))))))))
  ;;;; ===  
  (display "=== Registered Problems ===\n")
  (let loop ((todo (table->list _problem-registry)))
    (if (null? todo) #t
        (begin
          (display (caar todo))
          (display "\n")
          (loop (cdr todo)))))
  (display "======\n")
  )

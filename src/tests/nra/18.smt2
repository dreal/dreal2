(set-logic QF_NRA)
(declare-fun x1 () Real)
(declare-fun x2 () Real)

(assert (< x1 2))
(assert (< x2 2))
(assert (> x1 -2))
(assert (> x2 -2))

(assert (<= (+ (^ x1 2) (^ x2 3)) 0))

(check-sat)
(exit)

(set-logic QF_NRA)
(declare-fun x1 () Real)
(declare-fun x2 () Real)

(assert (< x1 2))
(assert (< x2 2))
(assert (> x1 -2))
(assert (> x2 -2))

(assert (< (+ (^ x1 2) (sin x2)) 0))

(check-sat)
(exit)

(set-logic QF_NRA)
(declare-fun skoX () Real)
(declare-fun skoY () Real)
(declare-fun pi () Real)
(assert (and (<= (+ (* 2 skoY) (* (- 1) pi)) 0) (and (<= (* (- 20) skoX) (- 1)) (not (<= (+ (* (- 1) skoX) skoY) 0)))))
(set-info :status sat)
(check-sat)
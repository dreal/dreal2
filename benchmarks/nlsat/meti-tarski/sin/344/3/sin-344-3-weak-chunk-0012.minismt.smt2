(set-logic QF_NRA)
(declare-fun skoY () Real)
(declare-fun skoX () Real)
(declare-fun skoZ () Real)
(assert (and (not (<= (+ skoY skoX skoZ) 0)) (and (not (<= (* (- 10) skoX) (- 11))) (and (not (<= (* (- 10) skoY) (- 11))) (and (not (<= (* (- 10) skoZ) (- 11))) (and (not (<= (* 10 skoX) 1)) (and (not (<= (* 10 skoY) 1)) (not (<= (* 10 skoZ) 1)))))))))
(set-info :status sat)
(check-sat)
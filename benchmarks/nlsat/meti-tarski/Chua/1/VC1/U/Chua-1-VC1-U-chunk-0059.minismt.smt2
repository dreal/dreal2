(set-logic QF_NRA)
(declare-fun skoS () Real)
(declare-fun skoC () Real)
(declare-fun skoX () Real)
(assert (and (not (<= (* 19 skoX) 1000)) (or (not (<= (+ (* (- 689) skoS) (* 1770 skoC)) 0)) (not (<= (+ (* 689 skoS) (* (- 1770) skoC)) 0)))))
(set-info :status sat)
(check-sat)
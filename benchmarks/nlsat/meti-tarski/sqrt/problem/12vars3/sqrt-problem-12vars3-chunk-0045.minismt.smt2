(set-logic QF_NRA)
(declare-fun skoSMX () Real)
(declare-fun skoSX () Real)
(declare-fun skoX () Real)
(assert (and (not (<= (+ (* (- 1) skoSMX) skoSX) 0)) (and (not (<= skoX 0)) (and (<= (* (- 1) skoSMX) 0) (and (<= (* (- 1) skoSX) 0) (and (<= skoX 1) (and (= (+ (* (- 1) skoX) (* skoSX skoSX)) 1) (= (+ skoX (* skoSMX skoSMX)) 1))))))))
(set-info :status sat)
(check-sat)
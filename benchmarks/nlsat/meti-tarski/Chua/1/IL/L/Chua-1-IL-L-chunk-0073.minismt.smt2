(set-logic QF_NRA)
(declare-fun skoX () Real)
(declare-fun skoS () Real)
(declare-fun skoC () Real)
(assert (and (<= (+ (* (- 42) skoS) (* (- 235) skoC)) 0) (and (<= (+ (* (- 5175000) skoX) (* 207 (* skoX skoX))) (- 43125000000)) (and (or (not (<= (+ (* 42 skoS) (* 235 skoC)) 0)) (<= skoX 0)) (and (or (<= (+ (* (- 42) skoS) (* (- 235) skoC)) 0) (<= skoX 0)) (and (= (+ (* skoS skoS) (* skoC skoC)) 1) (and (<= skoX 289) (<= (* (- 1) skoX) 0))))))))
(set-info :status unsat)
(check-sat)
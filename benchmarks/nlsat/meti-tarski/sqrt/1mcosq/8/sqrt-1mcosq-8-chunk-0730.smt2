(set-logic QF_NRA)

(declare-fun skoY () Real)
(declare-fun skoX () Real)
(declare-fun pi () Real)
(assert (and (<= (* skoX (* skoX (/ (- 1.) 4.))) (/ (- 1.) 2.)) (and (<= (* skoX (* skoX (/ 1. 4.))) (/ 1. 2.)) (and (not (<= (* skoY (+ (* skoX (* pi (- 2.))) (* skoY (* skoY (+ (* skoX (* pi (/ 2. 3.))) (* skoY (* skoY (+ (* skoX (* pi (/ (- 1.) 12.))) (* skoY (* skoY (* skoX (* pi (/ 1. 288.))))))))))))) (- 1.))) (and (<= skoY (+ (/ (- 1.) 5.) (* pi (/ 1. 2.)))) (and (not (<= pi (/ 15707963. 5000000.))) (and (not (<= (/ 31415927. 10000000.) pi)) (and (<= (/ 1. 10.) skoX) (not (<= skoY skoX))))))))))
(set-info :status unsat)
(check-sat)
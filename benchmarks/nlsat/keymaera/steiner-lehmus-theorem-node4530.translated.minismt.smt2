(set-logic QF_NRA)
(declare-fun y0 () Real)
(declare-fun x1uscore0dollarsk!3 () Real)
(declare-fun x0 () Real)
(declare-fun y1uscore0dollarsk!2 () Real)
(declare-fun cuscore0dollarsk!5 () Real)
(declare-fun buscore0dollarsk!6 () Real)
(declare-fun xuscore0dollarsk!4 () Real)
(declare-fun x2uscore0dollarsk!1 () Real)
(declare-fun y2uscore0dollarsk!0 () Real)
(assert (not (<= (+ y1uscore0dollarsk!2 (* x0 y1uscore0dollarsk!2) (* (- 1) (* y0 x1uscore0dollarsk!3))) 0)))
(assert (not (<= y0 0)))
(assert (not (<= cuscore0dollarsk!5 0)))
(assert (= (+ (* cuscore0dollarsk!5 cuscore0dollarsk!5) (* (- 1) (* buscore0dollarsk!6 buscore0dollarsk!6))) 1))
(assert (= (+ (* cuscore0dollarsk!5 cuscore0dollarsk!5) (* (- 1) (* buscore0dollarsk!6 buscore0dollarsk!6)) (* (- 1) (* x0 x0)) (* (- 1) (* y0 y0)) (* 2 (* y0 buscore0dollarsk!6))) 0))
(assert (= (+ (* y0 xuscore0dollarsk!4) (* cuscore0dollarsk!5 xuscore0dollarsk!4) (* (- 1) (* buscore0dollarsk!6 xuscore0dollarsk!4)) (* (- 1) (* x0 cuscore0dollarsk!5)) (* x0 buscore0dollarsk!6)) 0))
(assert (= (+ (* (- 1) (* cuscore0dollarsk!5 cuscore0dollarsk!5)) (* buscore0dollarsk!6 buscore0dollarsk!6) (* x1uscore0dollarsk!3 x1uscore0dollarsk!3) (* y1uscore0dollarsk!2 y1uscore0dollarsk!2) (* (- 2) (* y1uscore0dollarsk!2 buscore0dollarsk!6))) 0))
(assert (= (+ (* 2 x1uscore0dollarsk!3) (* (- 1) (* x0 x0)) (* (- 1) (* y0 y0)) (* 2 (* x1uscore0dollarsk!3 x0)) (* 2 (* y0 y1uscore0dollarsk!2))) (- 1)))
(assert (= (+ y1uscore0dollarsk!2 (* (- 1) y2uscore0dollarsk!0) (* x1uscore0dollarsk!3 y2uscore0dollarsk!0) (* (- 1) (* y1uscore0dollarsk!2 x2uscore0dollarsk!1))) 0))
(assert (= (+ (* (- 1) y0) y2uscore0dollarsk!0 (* x0 y2uscore0dollarsk!0) (* (- 1) (* y0 x2uscore0dollarsk!1))) 0))
(assert (not (<= (+ (* (- 2) (* x0 x0)) (* (- 1) (* y0 y0)) (* (- 1) (* x0 x0 x0 x0))) (- 3))))
(assert (not (= (+ (* (- 2) x0) (* x0 x0) (* y0 y0)) 3)))
(assert (not (= (+ (* 2 x0) (* x0 x0) (* y0 y0)) (- 1))))
(assert (not (= x0 (- 1))))
(assert (not (= x0 0)))
(assert (not (= y0 0)))
(assert (<= (+ (* (- 2) x2uscore0dollarsk!1) (* (- 1) (* x0 x0)) (* (- 1) (* y0 y0)) (* y2uscore0dollarsk!0 y2uscore0dollarsk!0) (* x2uscore0dollarsk!1 x2uscore0dollarsk!1) (* 2 (* x0 xuscore0dollarsk!4)) (* (- 1) (* xuscore0dollarsk!4 xuscore0dollarsk!4))) (- 1)))
(check-sat)
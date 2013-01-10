(set-logic QF_NRA)
(set-info :source |
From termination analysis of term rewriting.

Submitted by Harald Roman Zankl <Harald.Zankl@uibk.ac.at>

|)
(set-info :smt-lib-version 2.0)
(set-info :category "industrial")
(set-info :status unknown)
(declare-fun x6 () Real)
(declare-fun x23 () Real)
(declare-fun x40 () Real)
(declare-fun x57 () Real)
(declare-fun x13 () Real)
(declare-fun x30 () Real)
(declare-fun x47 () Real)
(declare-fun x64 () Real)
(declare-fun x3 () Real)
(declare-fun x20 () Real)
(declare-fun x37 () Real)
(declare-fun x54 () Real)
(declare-fun x71 () Real)
(declare-fun x10 () Real)
(declare-fun x27 () Real)
(declare-fun x44 () Real)
(declare-fun x61 () Real)
(declare-fun x0 () Real)
(declare-fun x17 () Real)
(declare-fun x34 () Real)
(declare-fun x51 () Real)
(declare-fun x68 () Real)
(declare-fun x7 () Real)
(declare-fun x24 () Real)
(declare-fun x41 () Real)
(declare-fun x58 () Real)
(declare-fun x14 () Real)
(declare-fun x31 () Real)
(declare-fun x48 () Real)
(declare-fun x65 () Real)
(declare-fun x4 () Real)
(declare-fun x21 () Real)
(declare-fun x38 () Real)
(declare-fun x55 () Real)
(declare-fun x72 () Real)
(declare-fun x11 () Real)
(declare-fun x28 () Real)
(declare-fun x45 () Real)
(declare-fun x62 () Real)
(declare-fun x1 () Real)
(declare-fun x18 () Real)
(declare-fun x35 () Real)
(declare-fun x52 () Real)
(declare-fun x69 () Real)
(declare-fun x8 () Real)
(declare-fun x25 () Real)
(declare-fun x42 () Real)
(declare-fun x59 () Real)
(declare-fun x15 () Real)
(declare-fun x32 () Real)
(declare-fun x49 () Real)
(declare-fun x66 () Real)
(declare-fun x5 () Real)
(declare-fun x22 () Real)
(declare-fun x39 () Real)
(declare-fun x56 () Real)
(declare-fun x73 () Real)
(declare-fun x12 () Real)
(declare-fun x29 () Real)
(declare-fun x46 () Real)
(declare-fun x63 () Real)
(declare-fun x2 () Real)
(declare-fun x19 () Real)
(declare-fun x36 () Real)
(declare-fun x53 () Real)
(declare-fun x70 () Real)
(declare-fun x9 () Real)
(declare-fun x26 () Real)
(declare-fun x43 () Real)
(declare-fun x60 () Real)
(declare-fun x16 () Real)
(declare-fun x33 () Real)
(declare-fun x50 () Real)
(declare-fun x67 () Real)
(assert (>= x6 0))
(assert (>= x23 0))
(assert (>= x40 0))
(assert (>= x57 0))
(assert (>= x13 0))
(assert (>= x30 0))
(assert (>= x47 0))
(assert (>= x64 0))
(assert (>= x3 0))
(assert (>= x20 0))
(assert (>= x37 0))
(assert (>= x54 0))
(assert (>= x71 0))
(assert (>= x10 0))
(assert (>= x27 0))
(assert (>= x44 0))
(assert (>= x61 0))
(assert (>= x0 0))
(assert (>= x17 0))
(assert (>= x34 0))
(assert (>= x51 0))
(assert (>= x68 0))
(assert (>= x7 0))
(assert (>= x24 0))
(assert (>= x41 0))
(assert (>= x58 0))
(assert (>= x14 0))
(assert (>= x31 0))
(assert (>= x48 0))
(assert (>= x65 0))
(assert (>= x4 0))
(assert (>= x21 0))
(assert (>= x38 0))
(assert (>= x55 0))
(assert (>= x72 0))
(assert (>= x11 0))
(assert (>= x28 0))
(assert (>= x45 0))
(assert (>= x62 0))
(assert (>= x1 0))
(assert (>= x18 0))
(assert (>= x35 0))
(assert (>= x52 0))
(assert (>= x69 0))
(assert (>= x8 0))
(assert (>= x25 0))
(assert (>= x42 0))
(assert (>= x59 0))
(assert (>= x15 0))
(assert (>= x32 0))
(assert (>= x49 0))
(assert (>= x66 0))
(assert (>= x5 0))
(assert (>= x22 0))
(assert (>= x39 0))
(assert (>= x56 0))
(assert (>= x73 0))
(assert (>= x12 0))
(assert (>= x29 0))
(assert (>= x46 0))
(assert (>= x63 0))
(assert (>= x2 0))
(assert (>= x19 0))
(assert (>= x36 0))
(assert (>= x53 0))
(assert (>= x70 0))
(assert (>= x9 0))
(assert (>= x26 0))
(assert (>= x43 0))
(assert (>= x60 0))
(assert (>= x16 0))
(assert (>= x33 0))
(assert (>= x50 0))
(assert (>= x67 0))
(assert (let ((?v_0 (+ x0 (+ (+ (+ (* x1 x5) (* x2 x6)) (* x3 x7)) (* x4 x8)))) (?v_2 (+ (+ (+ (* x1 x9) (* x2 x13)) (* x3 x17)) (* x4 x21))) (?v_3 (+ (+ (+ (* x1 x10) (* x2 x14)) (* x3 x18)) (* x4 x22))) (?v_4 (+ (+ (+ (* x1 x11) (* x2 x15)) (* x3 x19)) (* x4 x23))) (?v_5 (+ (+ (+ (* x1 x12) (* x2 x16)) (* x3 x20)) (* x4 x24))) (?v_1 (+ x0 (+ (+ (+ (* x1 x30) (* x2 x31)) (* x3 x32)) (* x4 x33))))) (let ((?v_11 (and (and (and (> ?v_0 x25) (>= ?v_0 x25)) (and (and (and (>= ?v_2 x26) (>= ?v_3 x27)) (>= ?v_4 x28)) (>= ?v_5 x29))) (and (and (> ?v_0 ?v_1) (>= ?v_0 ?v_1)) (and (and (and (>= ?v_2 (+ (+ (+ (* x1 x34) (* x2 x38)) (* x3 x42)) (* x4 x46))) (>= ?v_3 (+ (+ (+ (* x1 x35) (* x2 x39)) (* x3 x43)) (* x4 x47)))) (>= ?v_4 (+ (+ (+ (* x1 x36) (* x2 x40)) (* x3 x44)) (* x4 x48)))) (>= ?v_5 (+ (+ (+ (* x1 x37) (* x2 x41)) (* x3 x45)) (* x4 x49))))))) (?v_7 (+ x50 (+ (+ (+ (* x54 x30) (* x55 x31)) (* x56 x32)) (* x57 x33)))) (?v_6 (+ x50 (+ (+ (+ (* x54 x5) (* x55 x6)) (* x56 x7)) (* x57 x8)))) (?v_8 (+ x51 (+ (+ (+ (* x58 x30) (* x59 x31)) (* x60 x32)) (* x61 x33)))) (?v_9 (+ x52 (+ (+ (+ (* x62 x30) (* x63 x31)) (* x64 x32)) (* x65 x33)))) (?v_10 (+ x53 (+ (+ (+ (* x66 x30) (* x67 x31)) (* x68 x32)) (* x69 x33))))) (and (and (and (and ?v_11 (and (and (> ?v_6 ?v_7) (and (and (and (>= ?v_6 ?v_7) (>= (+ x51 (+ (+ (+ (* x58 x5) (* x59 x6)) (* x60 x7)) (* x61 x8))) ?v_8)) (>= (+ x52 (+ (+ (+ (* x62 x5) (* x63 x6)) (* x64 x7)) (* x65 x8))) ?v_9)) (>= (+ x53 (+ (+ (+ (* x66 x5) (* x67 x6)) (* x68 x7)) (* x69 x8))) ?v_10))) (and (and (and (and (and (and (and (and (and (and (and (and (and (and (and (>= (+ (+ (+ (* x54 x9) (* x55 x13)) (* x56 x17)) (* x57 x21)) (+ (+ (+ (* x54 x34) (* x55 x38)) (* x56 x42)) (* x57 x46))) (>= (+ (+ (+ (* x54 x10) (* x55 x14)) (* x56 x18)) (* x57 x22)) (+ (+ (+ (* x54 x35) (* x55 x39)) (* x56 x43)) (* x57 x47)))) (>= (+ (+ (+ (* x54 x11) (* x55 x15)) (* x56 x19)) (* x57 x23)) (+ (+ (+ (* x54 x36) (* x55 x40)) (* x56 x44)) (* x57 x48)))) (>= (+ (+ (+ (* x54 x12) (* x55 x16)) (* x56 x20)) (* x57 x24)) (+ (+ (+ (* x54 x37) (* x55 x41)) (* x56 x45)) (* x57 x49)))) (>= (+ (+ (+ (* x58 x9) (* x59 x13)) (* x60 x17)) (* x61 x21)) (+ (+ (+ (* x58 x34) (* x59 x38)) (* x60 x42)) (* x61 x46)))) (>= (+ (+ (+ (* x58 x10) (* x59 x14)) (* x60 x18)) (* x61 x22)) (+ (+ (+ (* x58 x35) (* x59 x39)) (* x60 x43)) (* x61 x47)))) (>= (+ (+ (+ (* x58 x11) (* x59 x15)) (* x60 x19)) (* x61 x23)) (+ (+ (+ (* x58 x36) (* x59 x40)) (* x60 x44)) (* x61 x48)))) (>= (+ (+ (+ (* x58 x12) (* x59 x16)) (* x60 x20)) (* x61 x24)) (+ (+ (+ (* x58 x37) (* x59 x41)) (* x60 x45)) (* x61 x49)))) (>= (+ (+ (+ (* x62 x9) (* x63 x13)) (* x64 x17)) (* x65 x21)) (+ (+ (+ (* x62 x34) (* x63 x38)) (* x64 x42)) (* x65 x46)))) (>= (+ (+ (+ (* x62 x10) (* x63 x14)) (* x64 x18)) (* x65 x22)) (+ (+ (+ (* x62 x35) (* x63 x39)) (* x64 x43)) (* x65 x47)))) (>= (+ (+ (+ (* x62 x11) (* x63 x15)) (* x64 x19)) (* x65 x23)) (+ (+ (+ (* x62 x36) (* x63 x40)) (* x64 x44)) (* x65 x48)))) (>= (+ (+ (+ (* x62 x12) (* x63 x16)) (* x64 x20)) (* x65 x24)) (+ (+ (+ (* x62 x37) (* x63 x41)) (* x64 x45)) (* x65 x49)))) (>= (+ (+ (+ (* x66 x9) (* x67 x13)) (* x68 x17)) (* x69 x21)) (+ (+ (+ (* x66 x34) (* x67 x38)) (* x68 x42)) (* x69 x46)))) (>= (+ (+ (+ (* x66 x10) (* x67 x14)) (* x68 x18)) (* x69 x22)) (+ (+ (+ (* x66 x35) (* x67 x39)) (* x68 x43)) (* x69 x47)))) (>= (+ (+ (+ (* x66 x11) (* x67 x15)) (* x68 x19)) (* x69 x23)) (+ (+ (+ (* x66 x36) (* x67 x40)) (* x68 x44)) (* x69 x48)))) (>= (+ (+ (+ (* x66 x12) (* x67 x16)) (* x68 x20)) (* x69 x24)) (+ (+ (+ (* x66 x37) (* x67 x41)) (* x68 x45)) (* x69 x49)))))) (and (> ?v_7 x70) (and (and (and (>= ?v_7 x70) (>= ?v_8 x71)) (>= ?v_9 x72)) (>= ?v_10 x73)))) (and (and (> x30 x5) (and (and (and (>= x30 x5) (>= x31 x6)) (>= x32 x7)) (>= x33 x8))) (and (and (and (and (and (and (and (and (and (and (and (and (and (and (and (>= x34 x9) (>= x35 x10)) (>= x36 x11)) (>= x37 x12)) (>= x38 x13)) (>= x39 x14)) (>= x40 x15)) (>= x41 x16)) (>= x42 x17)) (>= x43 x18)) (>= x44 x19)) (>= x45 x20)) (>= x46 x21)) (>= x47 x22)) (>= x48 x23)) (>= x49 x24)))) ?v_11))))
(check-sat)
(exit)
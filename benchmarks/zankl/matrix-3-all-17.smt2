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
(declare-fun x13 () Real)
(declare-fun x3 () Real)
(declare-fun x20 () Real)
(declare-fun x10 () Real)
(declare-fun x27 () Real)
(declare-fun x0 () Real)
(declare-fun x17 () Real)
(declare-fun x7 () Real)
(declare-fun x24 () Real)
(declare-fun x14 () Real)
(declare-fun x4 () Real)
(declare-fun x21 () Real)
(declare-fun x11 () Real)
(declare-fun x1 () Real)
(declare-fun x18 () Real)
(declare-fun x8 () Real)
(declare-fun x25 () Real)
(declare-fun x15 () Real)
(declare-fun x5 () Real)
(declare-fun x22 () Real)
(declare-fun x12 () Real)
(declare-fun x2 () Real)
(declare-fun x19 () Real)
(declare-fun x9 () Real)
(declare-fun x26 () Real)
(declare-fun x16 () Real)
(assert (>= x6 0))
(assert (>= x23 0))
(assert (>= x13 0))
(assert (>= x3 0))
(assert (>= x20 0))
(assert (>= x10 0))
(assert (>= x27 0))
(assert (>= x0 0))
(assert (>= x17 0))
(assert (>= x7 0))
(assert (>= x24 0))
(assert (>= x14 0))
(assert (>= x4 0))
(assert (>= x21 0))
(assert (>= x11 0))
(assert (>= x1 0))
(assert (>= x18 0))
(assert (>= x8 0))
(assert (>= x25 0))
(assert (>= x15 0))
(assert (>= x5 0))
(assert (>= x22 0))
(assert (>= x12 0))
(assert (>= x2 0))
(assert (>= x19 0))
(assert (>= x9 0))
(assert (>= x26 0))
(assert (>= x16 0))
(assert (let ((?v_0 (+ x16 (+ (+ (* x19 x16) (* x20 x17)) (* x21 x18)))) (?v_1 (+ x17 (+ (+ (* x22 x16) (* x23 x17)) (* x24 x18)))) (?v_2 (+ x18 (+ (+ (* x25 x16) (* x26 x17)) (* x27 x18))))) (let ((?v_18 (+ x4 (+ (+ (* x7 ?v_0) (* x8 ?v_1)) (* x9 ?v_2)))) (?v_19 (+ x5 (+ (+ (* x10 ?v_0) (* x11 ?v_1)) (* x12 ?v_2)))) (?v_20 (+ x6 (+ (+ (* x13 ?v_0) (* x14 ?v_1)) (* x15 ?v_2))))) (let ((?v_3 (+ x0 (+ (+ (* x1 ?v_18) (* x2 ?v_19)) (* x3 ?v_20)))) (?v_4 (+ (+ (* x19 x19) (* x20 x22)) (* x21 x25))) (?v_5 (+ (+ (* x22 x19) (* x23 x22)) (* x24 x25))) (?v_6 (+ (+ (* x25 x19) (* x26 x22)) (* x27 x25)))) (let ((?v_35 (+ (+ (* x7 ?v_4) (* x8 ?v_5)) (* x9 ?v_6))) (?v_36 (+ (+ (* x10 ?v_4) (* x11 ?v_5)) (* x12 ?v_6))) (?v_37 (+ (+ (* x13 ?v_4) (* x14 ?v_5)) (* x15 ?v_6)))) (let ((?v_14 (+ (+ (* x1 ?v_35) (* x2 ?v_36)) (* x3 ?v_37))) (?v_7 (+ (+ (* x19 x20) (* x20 x23)) (* x21 x26))) (?v_8 (+ (+ (* x22 x20) (* x23 x23)) (* x24 x26))) (?v_9 (+ (+ (* x25 x20) (* x26 x23)) (* x27 x26)))) (let ((?v_47 (+ (+ (* x7 ?v_7) (* x8 ?v_8)) (* x9 ?v_9))) (?v_48 (+ (+ (* x10 ?v_7) (* x11 ?v_8)) (* x12 ?v_9))) (?v_49 (+ (+ (* x13 ?v_7) (* x14 ?v_8)) (* x15 ?v_9)))) (let ((?v_15 (+ (+ (* x1 ?v_47) (* x2 ?v_48)) (* x3 ?v_49))) (?v_10 (+ (+ (* x19 x21) (* x20 x24)) (* x21 x27))) (?v_11 (+ (+ (* x22 x21) (* x23 x24)) (* x24 x27))) (?v_12 (+ (+ (* x25 x21) (* x26 x24)) (* x27 x27)))) (let ((?v_59 (+ (+ (* x7 ?v_10) (* x8 ?v_11)) (* x9 ?v_12))) (?v_60 (+ (+ (* x10 ?v_10) (* x11 ?v_11)) (* x12 ?v_12))) (?v_61 (+ (+ (* x13 ?v_10) (* x14 ?v_11)) (* x15 ?v_12)))) (let ((?v_16 (+ (+ (* x1 ?v_59) (* x2 ?v_60)) (* x3 ?v_61))) (?v_13 (+ x0 (+ (+ (* x1 x4) (* x2 x5)) (* x3 x6)))) (?v_21 (+ x4 (+ (+ (* x7 x4) (* x8 x5)) (* x9 x6)))) (?v_22 (+ x5 (+ (+ (* x10 x4) (* x11 x5)) (* x12 x6)))) (?v_23 (+ x6 (+ (+ (* x13 x4) (* x14 x5)) (* x15 x6))))) (let ((?v_17 (+ x0 (+ (+ (* x1 ?v_21) (* x2 ?v_22)) (* x3 ?v_23)))) (?v_38 (+ (+ (* x7 x7) (* x8 x10)) (* x9 x13))) (?v_39 (+ (+ (* x10 x7) (* x11 x10)) (* x12 x13))) (?v_40 (+ (+ (* x13 x7) (* x14 x10)) (* x15 x13))) (?v_50 (+ (+ (* x7 x8) (* x8 x11)) (* x9 x14))) (?v_51 (+ (+ (* x10 x8) (* x11 x11)) (* x12 x14))) (?v_52 (+ (+ (* x13 x8) (* x14 x11)) (* x15 x14))) (?v_62 (+ (+ (* x7 x9) (* x8 x12)) (* x9 x15))) (?v_63 (+ (+ (* x10 x9) (* x11 x12)) (* x12 x15))) (?v_64 (+ (+ (* x13 x9) (* x14 x12)) (* x15 x15)))) (let ((?v_80 (and (and (and (and (> ?v_3 x0) (>= ?v_3 x0)) (and (and (>= ?v_14 x1) (>= ?v_15 x2)) (>= ?v_16 x3))) (and (and (> ?v_3 ?v_13) (>= ?v_3 ?v_13)) (and (and (>= ?v_14 (+ (+ (* x1 x7) (* x2 x10)) (* x3 x13))) (>= ?v_15 (+ (+ (* x1 x8) (* x2 x11)) (* x3 x14)))) (>= ?v_16 (+ (+ (* x1 x9) (* x2 x12)) (* x3 x15)))))) (and (and (> ?v_3 ?v_17) (>= ?v_3 ?v_17)) (and (and (>= ?v_14 (+ (+ (* x1 ?v_38) (* x2 ?v_39)) (* x3 ?v_40))) (>= ?v_15 (+ (+ (* x1 ?v_50) (* x2 ?v_51)) (* x3 ?v_52)))) (>= ?v_16 (+ (+ (* x1 ?v_62) (* x2 ?v_63)) (* x3 ?v_64))))))) (?v_24 (+ x4 (+ (+ (* x7 ?v_21) (* x8 ?v_22)) (* x9 ?v_23)))) (?v_25 (+ x5 (+ (+ (* x10 ?v_21) (* x11 ?v_22)) (* x12 ?v_23)))) (?v_26 (+ x6 (+ (+ (* x13 ?v_21) (* x14 ?v_22)) (* x15 ?v_23))))) (let ((?v_27 (+ x16 (+ (+ (* x19 ?v_24) (* x20 ?v_25)) (* x21 ?v_26)))) (?v_28 (+ x17 (+ (+ (* x22 ?v_24) (* x23 ?v_25)) (* x24 ?v_26)))) (?v_29 (+ x18 (+ (+ (* x25 ?v_24) (* x26 ?v_25)) (* x27 ?v_26))))) (let ((?v_32 (+ x16 (+ (+ (* x19 ?v_27) (* x20 ?v_28)) (* x21 ?v_29)))) (?v_33 (+ x17 (+ (+ (* x22 ?v_27) (* x23 ?v_28)) (* x24 ?v_29)))) (?v_34 (+ x18 (+ (+ (* x25 ?v_27) (* x26 ?v_28)) (* x27 ?v_29))))) (let ((?v_31 (+ x16 (+ (+ (* x19 ?v_32) (* x20 ?v_33)) (* x21 ?v_34)))) (?v_30 (+ x4 (+ (+ (* x7 ?v_18) (* x8 ?v_19)) (* x9 ?v_20)))) (?v_41 (+ (+ (* x7 ?v_38) (* x8 ?v_39)) (* x9 ?v_40))) (?v_42 (+ (+ (* x10 ?v_38) (* x11 ?v_39)) (* x12 ?v_40))) (?v_43 (+ (+ (* x13 ?v_38) (* x14 ?v_39)) (* x15 ?v_40)))) (let ((?v_44 (+ (+ (* x19 ?v_41) (* x20 ?v_42)) (* x21 ?v_43))) (?v_45 (+ (+ (* x22 ?v_41) (* x23 ?v_42)) (* x24 ?v_43))) (?v_46 (+ (+ (* x25 ?v_41) (* x26 ?v_42)) (* x27 ?v_43)))) (let ((?v_71 (+ (+ (* x19 ?v_44) (* x20 ?v_45)) (* x21 ?v_46))) (?v_72 (+ (+ (* x22 ?v_44) (* x23 ?v_45)) (* x24 ?v_46))) (?v_73 (+ (+ (* x25 ?v_44) (* x26 ?v_45)) (* x27 ?v_46))) (?v_53 (+ (+ (* x7 ?v_50) (* x8 ?v_51)) (* x9 ?v_52))) (?v_54 (+ (+ (* x10 ?v_50) (* x11 ?v_51)) (* x12 ?v_52))) (?v_55 (+ (+ (* x13 ?v_50) (* x14 ?v_51)) (* x15 ?v_52)))) (let ((?v_56 (+ (+ (* x19 ?v_53) (* x20 ?v_54)) (* x21 ?v_55))) (?v_57 (+ (+ (* x22 ?v_53) (* x23 ?v_54)) (* x24 ?v_55))) (?v_58 (+ (+ (* x25 ?v_53) (* x26 ?v_54)) (* x27 ?v_55)))) (let ((?v_74 (+ (+ (* x19 ?v_56) (* x20 ?v_57)) (* x21 ?v_58))) (?v_75 (+ (+ (* x22 ?v_56) (* x23 ?v_57)) (* x24 ?v_58))) (?v_76 (+ (+ (* x25 ?v_56) (* x26 ?v_57)) (* x27 ?v_58))) (?v_65 (+ (+ (* x7 ?v_62) (* x8 ?v_63)) (* x9 ?v_64))) (?v_66 (+ (+ (* x10 ?v_62) (* x11 ?v_63)) (* x12 ?v_64))) (?v_67 (+ (+ (* x13 ?v_62) (* x14 ?v_63)) (* x15 ?v_64)))) (let ((?v_68 (+ (+ (* x19 ?v_65) (* x20 ?v_66)) (* x21 ?v_67))) (?v_69 (+ (+ (* x22 ?v_65) (* x23 ?v_66)) (* x24 ?v_67))) (?v_70 (+ (+ (* x25 ?v_65) (* x26 ?v_66)) (* x27 ?v_67)))) (let ((?v_77 (+ (+ (* x19 ?v_68) (* x20 ?v_69)) (* x21 ?v_70))) (?v_78 (+ (+ (* x22 ?v_68) (* x23 ?v_69)) (* x24 ?v_70))) (?v_79 (+ (+ (* x25 ?v_68) (* x26 ?v_69)) (* x27 ?v_70)))) (and (and ?v_80 (and (and (> ?v_30 ?v_31) (and (and (>= ?v_30 ?v_31) (>= (+ x5 (+ (+ (* x10 ?v_18) (* x11 ?v_19)) (* x12 ?v_20))) (+ x17 (+ (+ (* x22 ?v_32) (* x23 ?v_33)) (* x24 ?v_34))))) (>= (+ x6 (+ (+ (* x13 ?v_18) (* x14 ?v_19)) (* x15 ?v_20))) (+ x18 (+ (+ (* x25 ?v_32) (* x26 ?v_33)) (* x27 ?v_34)))))) (and (and (and (and (and (and (and (and (>= (+ (+ (* x7 ?v_35) (* x8 ?v_36)) (* x9 ?v_37)) (+ (+ (* x19 ?v_71) (* x20 ?v_72)) (* x21 ?v_73))) (>= (+ (+ (* x7 ?v_47) (* x8 ?v_48)) (* x9 ?v_49)) (+ (+ (* x19 ?v_74) (* x20 ?v_75)) (* x21 ?v_76)))) (>= (+ (+ (* x7 ?v_59) (* x8 ?v_60)) (* x9 ?v_61)) (+ (+ (* x19 ?v_77) (* x20 ?v_78)) (* x21 ?v_79)))) (>= (+ (+ (* x10 ?v_35) (* x11 ?v_36)) (* x12 ?v_37)) (+ (+ (* x22 ?v_71) (* x23 ?v_72)) (* x24 ?v_73)))) (>= (+ (+ (* x10 ?v_47) (* x11 ?v_48)) (* x12 ?v_49)) (+ (+ (* x22 ?v_74) (* x23 ?v_75)) (* x24 ?v_76)))) (>= (+ (+ (* x10 ?v_59) (* x11 ?v_60)) (* x12 ?v_61)) (+ (+ (* x22 ?v_77) (* x23 ?v_78)) (* x24 ?v_79)))) (>= (+ (+ (* x13 ?v_35) (* x14 ?v_36)) (* x15 ?v_37)) (+ (+ (* x25 ?v_71) (* x26 ?v_72)) (* x27 ?v_73)))) (>= (+ (+ (* x13 ?v_47) (* x14 ?v_48)) (* x15 ?v_49)) (+ (+ (* x25 ?v_74) (* x26 ?v_75)) (* x27 ?v_76)))) (>= (+ (+ (* x13 ?v_59) (* x14 ?v_60)) (* x15 ?v_61)) (+ (+ (* x25 ?v_77) (* x26 ?v_78)) (* x27 ?v_79)))))) ?v_80))))))))))))))))))))))
(check-sat)
(exit)
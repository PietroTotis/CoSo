universe_0(1,15).
universe(X) :- universe_0(X).
universe_1(16, 1).
universe_1(17, 1).
universe_1(18, 1).
universe_1(19, 1).
universe_1(20, 1).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(16, 1).
dom1_0(17, 1).
dom1_0(18, 1).
dom1_0(19, 1).
dom1_0(20, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
dom2_0(16, 1).
dom2_0(17, 1).
dom2_0(18, 1).
dom2_0(19, 1).
dom2_0(20, 1).
dom2(X) :- dom2_0(X).
universe(X) :- dom2(X).
multisubset_1(A) :- universe(A), .
used_1(X,0) :- multisubset_1(X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_2(A,B) :- universe(A), universe(B), A<=B.
used_2(X,0) :- multisubset_2(X, _). 
used_2(X,1) :- multisubset_2(_, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_3(A,B,C) :- universe(A), universe(B), universe(C), A<=B, B<=C.
used_3(X,0) :- multisubset_3(X, _, _). 
used_3(X,1) :- multisubset_3(_, X, _). 
used_3(X,2) :- multisubset_3(_, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_4(A,B,C,D) :- universe(A), universe(B), universe(C), universe(D), A<=B, B<=C, C<=D.
used_4(X,0) :- multisubset_4(X, _, _, _). 
used_4(X,1) :- multisubset_4(_, X, _, _). 
used_4(X,2) :- multisubset_4(_, _, X, _). 
used_4(X,3) :- multisubset_4(_, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_5(A,B,C,D,E) :- universe(A), universe(B), universe(C), universe(D), universe(E), A<=B, B<=C, C<=D, D<=E.
used_5(X,0) :- multisubset_5(X, _, _, _, _). 
used_5(X,1) :- multisubset_5(_, X, _, _, _). 
used_5(X,2) :- multisubset_5(_, _, X, _, _). 
used_5(X,3) :- multisubset_5(_, _, _, X, _). 
used_5(X,4) :- multisubset_5(_, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_6(A,B,C,D,E,F) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), A<=B, B<=C, C<=D, D<=E, E<=F.
used_6(X,0) :- multisubset_6(X, _, _, _, _, _). 
used_6(X,1) :- multisubset_6(_, X, _, _, _, _). 
used_6(X,2) :- multisubset_6(_, _, X, _, _, _). 
used_6(X,3) :- multisubset_6(_, _, _, X, _, _). 
used_6(X,4) :- multisubset_6(_, _, _, _, X, _). 
used_6(X,5) :- multisubset_6(_, _, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_7(A,B,C,D,E,F,G) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G.
used_7(X,0) :- multisubset_7(X, _, _, _, _, _, _). 
used_7(X,1) :- multisubset_7(_, X, _, _, _, _, _). 
used_7(X,2) :- multisubset_7(_, _, X, _, _, _, _). 
used_7(X,3) :- multisubset_7(_, _, _, X, _, _, _). 
used_7(X,4) :- multisubset_7(_, _, _, _, X, _, _). 
used_7(X,5) :- multisubset_7(_, _, _, _, _, X, _). 
used_7(X,6) :- multisubset_7(_, _, _, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_8(A,B,C,D,E,F,G,H) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G, G<=H.
used_8(X,0) :- multisubset_8(X, _, _, _, _, _, _, _). 
used_8(X,1) :- multisubset_8(_, X, _, _, _, _, _, _). 
used_8(X,2) :- multisubset_8(_, _, X, _, _, _, _, _). 
used_8(X,3) :- multisubset_8(_, _, _, X, _, _, _, _). 
used_8(X,4) :- multisubset_8(_, _, _, _, X, _, _, _). 
used_8(X,5) :- multisubset_8(_, _, _, _, _, X, _, _). 
used_8(X,6) :- multisubset_8(_, _, _, _, _, _, X, _). 
used_8(X,7) :- multisubset_8(_, _, _, _, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_9(A,B,C,D,E,F,G,H,I) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G, G<=H, H<=I.
used_9(X,0) :- multisubset_9(X, _, _, _, _, _, _, _, _). 
used_9(X,1) :- multisubset_9(_, X, _, _, _, _, _, _, _). 
used_9(X,2) :- multisubset_9(_, _, X, _, _, _, _, _, _). 
used_9(X,3) :- multisubset_9(_, _, _, X, _, _, _, _, _). 
used_9(X,4) :- multisubset_9(_, _, _, _, X, _, _, _, _). 
used_9(X,5) :- multisubset_9(_, _, _, _, _, X, _, _, _). 
used_9(X,6) :- multisubset_9(_, _, _, _, _, _, X, _, _). 
used_9(X,7) :- multisubset_9(_, _, _, _, _, _, _, X, _). 
used_9(X,8) :- multisubset_9(_, _, _, _, _, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_10(A,B,C,D,E,F,G,H,I,J) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G, G<=H, H<=I, I<=J.
used_10(X,0) :- multisubset_10(X, _, _, _, _, _, _, _, _, _). 
used_10(X,1) :- multisubset_10(_, X, _, _, _, _, _, _, _, _). 
used_10(X,2) :- multisubset_10(_, _, X, _, _, _, _, _, _, _). 
used_10(X,3) :- multisubset_10(_, _, _, X, _, _, _, _, _, _). 
used_10(X,4) :- multisubset_10(_, _, _, _, X, _, _, _, _, _). 
used_10(X,5) :- multisubset_10(_, _, _, _, _, X, _, _, _, _). 
used_10(X,6) :- multisubset_10(_, _, _, _, _, _, X, _, _, _). 
used_10(X,7) :- multisubset_10(_, _, _, _, _, _, _, X, _, _). 
used_10(X,8) :- multisubset_10(_, _, _, _, _, _, _, _, X, _). 
used_10(X,9) :- multisubset_10(_, _, _, _, _, _, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_11(A,B,C,D,E,F,G,H,I,J,K) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G, G<=H, H<=I, I<=J, J<=K.
used_11(X,0) :- multisubset_11(X, _, _, _, _, _, _, _, _, _, _). 
used_11(X,1) :- multisubset_11(_, X, _, _, _, _, _, _, _, _, _). 
used_11(X,2) :- multisubset_11(_, _, X, _, _, _, _, _, _, _, _). 
used_11(X,3) :- multisubset_11(_, _, _, X, _, _, _, _, _, _, _). 
used_11(X,4) :- multisubset_11(_, _, _, _, X, _, _, _, _, _, _). 
used_11(X,5) :- multisubset_11(_, _, _, _, _, X, _, _, _, _, _). 
used_11(X,6) :- multisubset_11(_, _, _, _, _, _, X, _, _, _, _). 
used_11(X,7) :- multisubset_11(_, _, _, _, _, _, _, X, _, _, _). 
used_11(X,8) :- multisubset_11(_, _, _, _, _, _, _, _, X, _, _). 
used_11(X,9) :- multisubset_11(_, _, _, _, _, _, _, _, _, X, _). 
used_11(X,10) :- multisubset_11(_, _, _, _, _, _, _, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_12(A,B,C,D,E,F,G,H,I,J,K,L) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G, G<=H, H<=I, I<=J, J<=K, K<=L.
used_12(X,0) :- multisubset_12(X, _, _, _, _, _, _, _, _, _, _, _). 
used_12(X,1) :- multisubset_12(_, X, _, _, _, _, _, _, _, _, _, _). 
used_12(X,2) :- multisubset_12(_, _, X, _, _, _, _, _, _, _, _, _). 
used_12(X,3) :- multisubset_12(_, _, _, X, _, _, _, _, _, _, _, _). 
used_12(X,4) :- multisubset_12(_, _, _, _, X, _, _, _, _, _, _, _). 
used_12(X,5) :- multisubset_12(_, _, _, _, _, X, _, _, _, _, _, _). 
used_12(X,6) :- multisubset_12(_, _, _, _, _, _, X, _, _, _, _, _). 
used_12(X,7) :- multisubset_12(_, _, _, _, _, _, _, X, _, _, _, _). 
used_12(X,8) :- multisubset_12(_, _, _, _, _, _, _, _, X, _, _, _). 
used_12(X,9) :- multisubset_12(_, _, _, _, _, _, _, _, _, X, _, _). 
used_12(X,10) :- multisubset_12(_, _, _, _, _, _, _, _, _, _, X, _). 
used_12(X,11) :- multisubset_12(_, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_13(A,B,C,D,E,F,G,H,I,J,K,L,M) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G, G<=H, H<=I, I<=J, J<=K, K<=L, L<=M.
used_13(X,0) :- multisubset_13(X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_13(X,1) :- multisubset_13(_, X, _, _, _, _, _, _, _, _, _, _, _). 
used_13(X,2) :- multisubset_13(_, _, X, _, _, _, _, _, _, _, _, _, _). 
used_13(X,3) :- multisubset_13(_, _, _, X, _, _, _, _, _, _, _, _, _). 
used_13(X,4) :- multisubset_13(_, _, _, _, X, _, _, _, _, _, _, _, _). 
used_13(X,5) :- multisubset_13(_, _, _, _, _, X, _, _, _, _, _, _, _). 
used_13(X,6) :- multisubset_13(_, _, _, _, _, _, X, _, _, _, _, _, _). 
used_13(X,7) :- multisubset_13(_, _, _, _, _, _, _, X, _, _, _, _, _). 
used_13(X,8) :- multisubset_13(_, _, _, _, _, _, _, _, X, _, _, _, _). 
used_13(X,9) :- multisubset_13(_, _, _, _, _, _, _, _, _, X, _, _, _). 
used_13(X,10) :- multisubset_13(_, _, _, _, _, _, _, _, _, _, X, _, _). 
used_13(X,11) :- multisubset_13(_, _, _, _, _, _, _, _, _, _, _, X, _). 
used_13(X,12) :- multisubset_13(_, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_14(A,B,C,D,E,F,G,H,I,J,K,L,M,N) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G, G<=H, H<=I, I<=J, J<=K, K<=L, L<=M, M<=N.
used_14(X,0) :- multisubset_14(X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,1) :- multisubset_14(_, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,2) :- multisubset_14(_, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,3) :- multisubset_14(_, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_14(X,4) :- multisubset_14(_, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_14(X,5) :- multisubset_14(_, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_14(X,6) :- multisubset_14(_, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_14(X,7) :- multisubset_14(_, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_14(X,8) :- multisubset_14(_, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_14(X,9) :- multisubset_14(_, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_14(X,10) :- multisubset_14(_, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_14(X,11) :- multisubset_14(_, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_14(X,12) :- multisubset_14(_, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_14(X,13) :- multisubset_14(_, _, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_16(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G, G<=H, H<=I, I<=J, J<=K, K<=L, L<=M, M<=N, N<=O, O<=P.
used_16(X,0) :- multisubset_16(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,1) :- multisubset_16(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,2) :- multisubset_16(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,3) :- multisubset_16(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,4) :- multisubset_16(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,5) :- multisubset_16(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_16(X,6) :- multisubset_16(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_16(X,7) :- multisubset_16(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_16(X,8) :- multisubset_16(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_16(X,9) :- multisubset_16(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_16(X,10) :- multisubset_16(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_16(X,11) :- multisubset_16(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_16(X,12) :- multisubset_16(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_16(X,13) :- multisubset_16(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_16(X,14) :- multisubset_16(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_16(X,15) :- multisubset_16(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_17(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G, G<=H, H<=I, I<=J, J<=K, K<=L, L<=M, M<=N, N<=O, O<=P, P<=Q.
used_17(X,0) :- multisubset_17(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,1) :- multisubset_17(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,2) :- multisubset_17(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,3) :- multisubset_17(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,4) :- multisubset_17(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,5) :- multisubset_17(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,6) :- multisubset_17(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_17(X,7) :- multisubset_17(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_17(X,8) :- multisubset_17(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_17(X,9) :- multisubset_17(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_17(X,10) :- multisubset_17(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_17(X,11) :- multisubset_17(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_17(X,12) :- multisubset_17(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_17(X,13) :- multisubset_17(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_17(X,14) :- multisubset_17(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_17(X,15) :- multisubset_17(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_17(X,16) :- multisubset_17(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_18(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q), universe(R), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G, G<=H, H<=I, I<=J, J<=K, K<=L, L<=M, M<=N, N<=O, O<=P, P<=Q, Q<=R.
used_18(X,0) :- multisubset_18(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,1) :- multisubset_18(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,2) :- multisubset_18(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,3) :- multisubset_18(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,4) :- multisubset_18(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,5) :- multisubset_18(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,6) :- multisubset_18(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,7) :- multisubset_18(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_18(X,8) :- multisubset_18(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_18(X,9) :- multisubset_18(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_18(X,10) :- multisubset_18(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_18(X,11) :- multisubset_18(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_18(X,12) :- multisubset_18(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_18(X,13) :- multisubset_18(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_18(X,14) :- multisubset_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_18(X,15) :- multisubset_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_18(X,16) :- multisubset_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_18(X,17) :- multisubset_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_19(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q), universe(R), universe(S), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G, G<=H, H<=I, I<=J, J<=K, K<=L, L<=M, M<=N, N<=O, O<=P, P<=Q, Q<=R, R<=S.
used_19(X,0) :- multisubset_19(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,1) :- multisubset_19(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,2) :- multisubset_19(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,3) :- multisubset_19(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,4) :- multisubset_19(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,5) :- multisubset_19(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,6) :- multisubset_19(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,7) :- multisubset_19(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,8) :- multisubset_19(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_19(X,9) :- multisubset_19(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_19(X,10) :- multisubset_19(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_19(X,11) :- multisubset_19(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_19(X,12) :- multisubset_19(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_19(X,13) :- multisubset_19(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_19(X,14) :- multisubset_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_19(X,15) :- multisubset_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_19(X,16) :- multisubset_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_19(X,17) :- multisubset_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_19(X,18) :- multisubset_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.
multisubset_20(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q), universe(R), universe(S), universe(T), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G, G<=H, H<=I, I<=J, J<=K, K<=L, L<=M, M<=N, N<=O, O<=P, P<=Q, Q<=R, R<=S, S<=T.
used_20(X,0) :- multisubset_20(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,1) :- multisubset_20(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,2) :- multisubset_20(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,3) :- multisubset_20(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,4) :- multisubset_20(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,5) :- multisubset_20(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,6) :- multisubset_20(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,7) :- multisubset_20(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,8) :- multisubset_20(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,9) :- multisubset_20(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_20(X,10) :- multisubset_20(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_20(X,11) :- multisubset_20(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_20(X,12) :- multisubset_20(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_20(X,13) :- multisubset_20(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_20(X,14) :- multisubset_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_20(X,15) :- multisubset_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_20(X,16) :- multisubset_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_20(X,17) :- multisubset_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_20(X,18) :- multisubset_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_20(X,19) :- multisubset_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=14.
df_3_0(16, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=2.

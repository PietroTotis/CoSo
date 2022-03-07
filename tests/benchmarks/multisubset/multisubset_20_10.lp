universe_0(1, 1).
universe_0(19, 1).
universe_0(20, 1).
universe(X) :- universe_0(X).
universe_1(2,17).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(1, 1).
dom1_0(19, 1).
dom1_0(20, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
dom2_0(1, 1).
dom2_0(19, 1).
dom2_0(20, 1).
dom2(X) :- dom2_0(X).
universe(X) :- dom2(X).
multisubset_1(A) :- universe(A), .
used_1(X,0) :- multisubset_1(X). 
multisubset_2(A,B) :- universe(A), universe(B), A<=B.
used_2(X,0) :- multisubset_2(X, _). 
used_2(X,1) :- multisubset_2(_, X). 
multisubset_3(A,B,C) :- universe(A), universe(B), universe(C), A<=B, B<=C.
used_3(X,0) :- multisubset_3(X, _, _). 
used_3(X,1) :- multisubset_3(_, X, _). 
used_3(X,2) :- multisubset_3(_, _, X). 
multisubset_4(A,B,C,D) :- universe(A), universe(B), universe(C), universe(D), A<=B, B<=C, C<=D.
used_4(X,0) :- multisubset_4(X, _, _, _). 
used_4(X,1) :- multisubset_4(_, X, _, _). 
used_4(X,2) :- multisubset_4(_, _, X, _). 
used_4(X,3) :- multisubset_4(_, _, _, X). 
multisubset_5(A,B,C,D,E) :- universe(A), universe(B), universe(C), universe(D), universe(E), A<=B, B<=C, C<=D, D<=E.
used_5(X,0) :- multisubset_5(X, _, _, _, _). 
used_5(X,1) :- multisubset_5(_, X, _, _, _). 
used_5(X,2) :- multisubset_5(_, _, X, _, _). 
used_5(X,3) :- multisubset_5(_, _, _, X, _). 
used_5(X,4) :- multisubset_5(_, _, _, _, X). 
multisubset_6(A,B,C,D,E,F) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), A<=B, B<=C, C<=D, D<=E, E<=F.
used_6(X,0) :- multisubset_6(X, _, _, _, _, _). 
used_6(X,1) :- multisubset_6(_, X, _, _, _, _). 
used_6(X,2) :- multisubset_6(_, _, X, _, _, _). 
used_6(X,3) :- multisubset_6(_, _, _, X, _, _). 
used_6(X,4) :- multisubset_6(_, _, _, _, X, _). 
used_6(X,5) :- multisubset_6(_, _, _, _, _, X). 
multisubset_7(A,B,C,D,E,F,G) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G.
used_7(X,0) :- multisubset_7(X, _, _, _, _, _, _). 
used_7(X,1) :- multisubset_7(_, X, _, _, _, _, _). 
used_7(X,2) :- multisubset_7(_, _, X, _, _, _, _). 
used_7(X,3) :- multisubset_7(_, _, _, X, _, _, _). 
used_7(X,4) :- multisubset_7(_, _, _, _, X, _, _). 
used_7(X,5) :- multisubset_7(_, _, _, _, _, X, _). 
used_7(X,6) :- multisubset_7(_, _, _, _, _, _, X). 
multisubset_8(A,B,C,D,E,F,G,H) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G, G<=H.
used_8(X,0) :- multisubset_8(X, _, _, _, _, _, _, _). 
used_8(X,1) :- multisubset_8(_, X, _, _, _, _, _, _). 
used_8(X,2) :- multisubset_8(_, _, X, _, _, _, _, _). 
used_8(X,3) :- multisubset_8(_, _, _, X, _, _, _, _). 
used_8(X,4) :- multisubset_8(_, _, _, _, X, _, _, _). 
used_8(X,5) :- multisubset_8(_, _, _, _, _, X, _, _). 
used_8(X,6) :- multisubset_8(_, _, _, _, _, _, X, _). 
used_8(X,7) :- multisubset_8(_, _, _, _, _, _, _, X). 
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

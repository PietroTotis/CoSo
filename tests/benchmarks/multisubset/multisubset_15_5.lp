universe_0(1, 1).
universe_0(2, 1).
universe_0(3, 1).
universe_0(4, 1).
universe_0(8, 1).
universe_0(9, 1).
universe_0(10, 1).
universe_0(11, 1).
universe_0(12, 1).
universe_0(13, 1).
universe_0(14, 1).
universe_0(15, 1).
universe(X) :- universe_0(X).
universe_1(5,3).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(2, 1).
dom1_0(3, 1).
dom1_0(14, 1).
dom1_0(15, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
dom2_0(1, 1).
dom2_0(2, 1).
dom2_0(3, 1).
dom2_0(4, 1).
dom2_0(8, 1).
dom2_0(9, 1).
dom2_0(10, 1).
dom2_0(11, 1).
dom2_0(12, 1).
dom2_0(13, 1).
dom2_0(14, 1).
dom2_0(15, 1).
dom2(X) :- dom2_0(X).
universe(X) :- dom2(X).
dom3_0(3, 1).
dom3_0(4, 1).
dom3(X) :- dom3_0(X).
universe(X) :- dom3(X).
dom4_0(3, 1).
dom4_0(4, 1).
dom4(X) :- dom4_0(X).
universe(X) :- dom4(X).
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
multisubset_15(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), A<=B, B<=C, C<=D, D<=E, E<=F, F<=G, G<=H, H<=I, I<=J, J<=K, K<=L, L<=M, M<=N, N<=O.
used_15(X,0) :- multisubset_15(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,1) :- multisubset_15(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,2) :- multisubset_15(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,3) :- multisubset_15(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,4) :- multisubset_15(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_15(X,5) :- multisubset_15(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_15(X,6) :- multisubset_15(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_15(X,7) :- multisubset_15(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_15(X,8) :- multisubset_15(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_15(X,9) :- multisubset_15(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_15(X,10) :- multisubset_15(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_15(X,11) :- multisubset_15(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_15(X,12) :- multisubset_15(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_15(X,13) :- multisubset_15(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_15(X,14) :- multisubset_15(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 

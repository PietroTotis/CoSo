universe_0(1, 1).
universe_0(2, 1).
universe_0(3, 1).
universe_0(4, 1).
universe_0(5, 1).
universe_0(6, 1).
universe_0(7, 1).
universe_0(8, 1).
universe_0(9, 1).
universe_0(10, 1).
universe_0(11, 1).
universe_0(12, 1).
universe_0(13, 1).
universe_0(14, 1).
universe_0(15, 1).
universe(X) :- universe_0(X).
universe(X) :- universe(X).
dom1_0(1, 1).
dom1_0(2, 1).
dom1_0(3, 1).
dom1_0(4, 1).
dom1_0(5, 1).
dom1_0(6, 1).
dom1_0(7, 1).
dom1_0(8, 1).
dom1_0(9, 1).
dom1_0(10, 1).
dom1_0(11, 1).
dom1_0(12, 1).
dom1_0(13, 1).
dom1_0(14, 1).
dom1_0(15, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
dom2_0(1, 1).
dom2_0(12, 1).
dom2_0(13, 1).
dom2(X) :- dom2_0(X).
universe(X) :- dom2(X).
sequence_10(A,B,C,D,E,F,G,H,I,Jpf_1_0(1, 1).
pf_1_0(2, 1).
pf_1_0(3, 1).
pf_1_0(4, 1).
pf_1_0(5, 1).
pf_1_0(6, 1).
pf_1_0(7, 1).
pf_1_0(8, 1).
pf_1_0(9, 1).
pf_1_0(10, 1).
pf_1_0(11, 1).
pf_1_0(12, 1).
pf_1_0(13, 1).
pf_1_0(14, 1).
pf_1_0(15, 1).
pf_1(X) :- pf_1_0(X).
universe(X) :- pf_1(X).
pf_6(X) :- pf_6_0(X).
universe(X) :- pf_6(X).
pf_8_0(2, 1).
pf_8_0(3, 1).
pf_8_0(4, 1).
pf_8_0(5, 1).
pf_8_0(6, 1).
pf_8_0(7, 1).
pf_8_0(8, 1).
pf_8_0(9, 1).
pf_8_0(10, 1).
pf_8_0(11, 1).
pf_8_0(14, 1).
pf_8_0(15, 1).
pf_8(X) :- pf_8_0(X).
universe(X) :- pf_8(X).
pf_9_0(2, 1).
pf_9_0(3, 1).
pf_9_0(4, 1).
pf_9_0(5, 1).
pf_9_0(6, 1).
pf_9_0(7, 1).
pf_9_0(8, 1).
pf_9_0(9, 1).
pf_9_0(10, 1).
pf_9_0(11, 1).
pf_9_0(14, 1).
pf_9_0(15, 1).
pf_9(X) :- pf_9_0(X).
universe(X) :- pf_9(X).
) :- universe(A), pf_1(B), universe(C), universe(D), universe(E), universe(F), pf_6(G), universe(H), pf_8(I), pf_9(J).
used_10(X,0) :- sequence_10(X, _, _, _, _, _, _, _, _, _). 
used_10(X,1) :- sequence_10(_, X, _, _, _, _, _, _, _, _). 
used_10(X,2) :- sequence_10(_, _, X, _, _, _, _, _, _, _). 
used_10(X,3) :- sequence_10(_, _, _, X, _, _, _, _, _, _). 
used_10(X,4) :- sequence_10(_, _, _, _, X, _, _, _, _, _). 
used_10(X,5) :- sequence_10(_, _, _, _, _, X, _, _, _, _). 
used_10(X,6) :- sequence_10(_, _, _, _, _, _, X, _, _, _). 
used_10(X,7) :- sequence_10(_, _, _, _, _, _, _, X, _, _). 
used_10(X,8) :- sequence_10(_, _, _, _, _, _, _, _, X, _). 
used_10(X,9) :- sequence_10(_, _, _, _, _, _, _, _, _, X). 
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
sequence_11(A,B,C,D,E,F,G,H,I,J,Kpf_1_0(1, 1).
pf_1_0(2, 1).
pf_1_0(3, 1).
pf_1_0(4, 1).
pf_1_0(5, 1).
pf_1_0(6, 1).
pf_1_0(7, 1).
pf_1_0(8, 1).
pf_1_0(9, 1).
pf_1_0(10, 1).
pf_1_0(11, 1).
pf_1_0(12, 1).
pf_1_0(13, 1).
pf_1_0(14, 1).
pf_1_0(15, 1).
pf_1(X) :- pf_1_0(X).
universe(X) :- pf_1(X).
pf_6(X) :- pf_6_0(X).
universe(X) :- pf_6(X).
pf_8_0(2, 1).
pf_8_0(3, 1).
pf_8_0(4, 1).
pf_8_0(5, 1).
pf_8_0(6, 1).
pf_8_0(7, 1).
pf_8_0(8, 1).
pf_8_0(9, 1).
pf_8_0(10, 1).
pf_8_0(11, 1).
pf_8_0(14, 1).
pf_8_0(15, 1).
pf_8(X) :- pf_8_0(X).
universe(X) :- pf_8(X).
pf_9_0(2, 1).
pf_9_0(3, 1).
pf_9_0(4, 1).
pf_9_0(5, 1).
pf_9_0(6, 1).
pf_9_0(7, 1).
pf_9_0(8, 1).
pf_9_0(9, 1).
pf_9_0(10, 1).
pf_9_0(11, 1).
pf_9_0(14, 1).
pf_9_0(15, 1).
pf_9(X) :- pf_9_0(X).
universe(X) :- pf_9(X).
) :- universe(A), pf_1(B), universe(C), universe(D), universe(E), universe(F), pf_6(G), universe(H), pf_8(I), pf_9(J), universe(K).
used_11(X,0) :- sequence_11(X, _, _, _, _, _, _, _, _, _, _). 
used_11(X,1) :- sequence_11(_, X, _, _, _, _, _, _, _, _, _). 
used_11(X,2) :- sequence_11(_, _, X, _, _, _, _, _, _, _, _). 
used_11(X,3) :- sequence_11(_, _, _, X, _, _, _, _, _, _, _). 
used_11(X,4) :- sequence_11(_, _, _, _, X, _, _, _, _, _, _). 
used_11(X,5) :- sequence_11(_, _, _, _, _, X, _, _, _, _, _). 
used_11(X,6) :- sequence_11(_, _, _, _, _, _, X, _, _, _, _). 
used_11(X,7) :- sequence_11(_, _, _, _, _, _, _, X, _, _, _). 
used_11(X,8) :- sequence_11(_, _, _, _, _, _, _, _, X, _, _). 
used_11(X,9) :- sequence_11(_, _, _, _, _, _, _, _, _, X, _). 
used_11(X,10) :- sequence_11(_, _, _, _, _, _, _, _, _, _, X). 
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
sequence_12(A,B,C,D,E,F,G,H,I,J,K,Lpf_1_0(1, 1).
pf_1_0(2, 1).
pf_1_0(3, 1).
pf_1_0(4, 1).
pf_1_0(5, 1).
pf_1_0(6, 1).
pf_1_0(7, 1).
pf_1_0(8, 1).
pf_1_0(9, 1).
pf_1_0(10, 1).
pf_1_0(11, 1).
pf_1_0(12, 1).
pf_1_0(13, 1).
pf_1_0(14, 1).
pf_1_0(15, 1).
pf_1(X) :- pf_1_0(X).
universe(X) :- pf_1(X).
pf_6(X) :- pf_6_0(X).
universe(X) :- pf_6(X).
pf_8_0(2, 1).
pf_8_0(3, 1).
pf_8_0(4, 1).
pf_8_0(5, 1).
pf_8_0(6, 1).
pf_8_0(7, 1).
pf_8_0(8, 1).
pf_8_0(9, 1).
pf_8_0(10, 1).
pf_8_0(11, 1).
pf_8_0(14, 1).
pf_8_0(15, 1).
pf_8(X) :- pf_8_0(X).
universe(X) :- pf_8(X).
pf_9_0(2, 1).
pf_9_0(3, 1).
pf_9_0(4, 1).
pf_9_0(5, 1).
pf_9_0(6, 1).
pf_9_0(7, 1).
pf_9_0(8, 1).
pf_9_0(9, 1).
pf_9_0(10, 1).
pf_9_0(11, 1).
pf_9_0(14, 1).
pf_9_0(15, 1).
pf_9(X) :- pf_9_0(X).
universe(X) :- pf_9(X).
) :- universe(A), pf_1(B), universe(C), universe(D), universe(E), universe(F), pf_6(G), universe(H), pf_8(I), pf_9(J), universe(K), universe(L).
used_12(X,0) :- sequence_12(X, _, _, _, _, _, _, _, _, _, _, _). 
used_12(X,1) :- sequence_12(_, X, _, _, _, _, _, _, _, _, _, _). 
used_12(X,2) :- sequence_12(_, _, X, _, _, _, _, _, _, _, _, _). 
used_12(X,3) :- sequence_12(_, _, _, X, _, _, _, _, _, _, _, _). 
used_12(X,4) :- sequence_12(_, _, _, _, X, _, _, _, _, _, _, _). 
used_12(X,5) :- sequence_12(_, _, _, _, _, X, _, _, _, _, _, _). 
used_12(X,6) :- sequence_12(_, _, _, _, _, _, X, _, _, _, _, _). 
used_12(X,7) :- sequence_12(_, _, _, _, _, _, _, X, _, _, _, _). 
used_12(X,8) :- sequence_12(_, _, _, _, _, _, _, _, X, _, _, _). 
used_12(X,9) :- sequence_12(_, _, _, _, _, _, _, _, _, X, _, _). 
used_12(X,10) :- sequence_12(_, _, _, _, _, _, _, _, _, _, X, _). 
used_12(X,11) :- sequence_12(_, _, _, _, _, _, _, _, _, _, _, X). 
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
sequence_13(A,B,C,D,E,F,G,H,I,J,K,L,Mpf_1_0(1, 1).
pf_1_0(2, 1).
pf_1_0(3, 1).
pf_1_0(4, 1).
pf_1_0(5, 1).
pf_1_0(6, 1).
pf_1_0(7, 1).
pf_1_0(8, 1).
pf_1_0(9, 1).
pf_1_0(10, 1).
pf_1_0(11, 1).
pf_1_0(12, 1).
pf_1_0(13, 1).
pf_1_0(14, 1).
pf_1_0(15, 1).
pf_1(X) :- pf_1_0(X).
universe(X) :- pf_1(X).
pf_6(X) :- pf_6_0(X).
universe(X) :- pf_6(X).
pf_8_0(2, 1).
pf_8_0(3, 1).
pf_8_0(4, 1).
pf_8_0(5, 1).
pf_8_0(6, 1).
pf_8_0(7, 1).
pf_8_0(8, 1).
pf_8_0(9, 1).
pf_8_0(10, 1).
pf_8_0(11, 1).
pf_8_0(14, 1).
pf_8_0(15, 1).
pf_8(X) :- pf_8_0(X).
universe(X) :- pf_8(X).
pf_9_0(2, 1).
pf_9_0(3, 1).
pf_9_0(4, 1).
pf_9_0(5, 1).
pf_9_0(6, 1).
pf_9_0(7, 1).
pf_9_0(8, 1).
pf_9_0(9, 1).
pf_9_0(10, 1).
pf_9_0(11, 1).
pf_9_0(14, 1).
pf_9_0(15, 1).
pf_9(X) :- pf_9_0(X).
universe(X) :- pf_9(X).
) :- universe(A), pf_1(B), universe(C), universe(D), universe(E), universe(F), pf_6(G), universe(H), pf_8(I), pf_9(J), universe(K), universe(L), universe(M).
used_13(X,0) :- sequence_13(X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_13(X,1) :- sequence_13(_, X, _, _, _, _, _, _, _, _, _, _, _). 
used_13(X,2) :- sequence_13(_, _, X, _, _, _, _, _, _, _, _, _, _). 
used_13(X,3) :- sequence_13(_, _, _, X, _, _, _, _, _, _, _, _, _). 
used_13(X,4) :- sequence_13(_, _, _, _, X, _, _, _, _, _, _, _, _). 
used_13(X,5) :- sequence_13(_, _, _, _, _, X, _, _, _, _, _, _, _). 
used_13(X,6) :- sequence_13(_, _, _, _, _, _, X, _, _, _, _, _, _). 
used_13(X,7) :- sequence_13(_, _, _, _, _, _, _, X, _, _, _, _, _). 
used_13(X,8) :- sequence_13(_, _, _, _, _, _, _, _, X, _, _, _, _). 
used_13(X,9) :- sequence_13(_, _, _, _, _, _, _, _, _, X, _, _, _). 
used_13(X,10) :- sequence_13(_, _, _, _, _, _, _, _, _, _, X, _, _). 
used_13(X,11) :- sequence_13(_, _, _, _, _, _, _, _, _, _, _, X, _). 
used_13(X,12) :- sequence_13(_, _, _, _, _, _, _, _, _, _, _, _, X). 
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
sequence_14(A,B,C,D,E,F,G,H,I,J,K,L,M,Npf_1_0(1, 1).
pf_1_0(2, 1).
pf_1_0(3, 1).
pf_1_0(4, 1).
pf_1_0(5, 1).
pf_1_0(6, 1).
pf_1_0(7, 1).
pf_1_0(8, 1).
pf_1_0(9, 1).
pf_1_0(10, 1).
pf_1_0(11, 1).
pf_1_0(12, 1).
pf_1_0(13, 1).
pf_1_0(14, 1).
pf_1_0(15, 1).
pf_1(X) :- pf_1_0(X).
universe(X) :- pf_1(X).
pf_6(X) :- pf_6_0(X).
universe(X) :- pf_6(X).
pf_8_0(2, 1).
pf_8_0(3, 1).
pf_8_0(4, 1).
pf_8_0(5, 1).
pf_8_0(6, 1).
pf_8_0(7, 1).
pf_8_0(8, 1).
pf_8_0(9, 1).
pf_8_0(10, 1).
pf_8_0(11, 1).
pf_8_0(14, 1).
pf_8_0(15, 1).
pf_8(X) :- pf_8_0(X).
universe(X) :- pf_8(X).
pf_9_0(2, 1).
pf_9_0(3, 1).
pf_9_0(4, 1).
pf_9_0(5, 1).
pf_9_0(6, 1).
pf_9_0(7, 1).
pf_9_0(8, 1).
pf_9_0(9, 1).
pf_9_0(10, 1).
pf_9_0(11, 1).
pf_9_0(14, 1).
pf_9_0(15, 1).
pf_9(X) :- pf_9_0(X).
universe(X) :- pf_9(X).
) :- universe(A), pf_1(B), universe(C), universe(D), universe(E), universe(F), pf_6(G), universe(H), pf_8(I), pf_9(J), universe(K), universe(L), universe(M), universe(N).
used_14(X,0) :- sequence_14(X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,1) :- sequence_14(_, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,2) :- sequence_14(_, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,3) :- sequence_14(_, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_14(X,4) :- sequence_14(_, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_14(X,5) :- sequence_14(_, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_14(X,6) :- sequence_14(_, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_14(X,7) :- sequence_14(_, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_14(X,8) :- sequence_14(_, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_14(X,9) :- sequence_14(_, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_14(X,10) :- sequence_14(_, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_14(X,11) :- sequence_14(_, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_14(X,12) :- sequence_14(_, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_14(X,13) :- sequence_14(_, _, _, _, _, _, _, _, _, _, _, _, _, X). 
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
sequence_15(A,B,C,D,E,F,G,H,I,J,K,L,M,N,Opf_1_0(1, 1).
pf_1_0(2, 1).
pf_1_0(3, 1).
pf_1_0(4, 1).
pf_1_0(5, 1).
pf_1_0(6, 1).
pf_1_0(7, 1).
pf_1_0(8, 1).
pf_1_0(9, 1).
pf_1_0(10, 1).
pf_1_0(11, 1).
pf_1_0(12, 1).
pf_1_0(13, 1).
pf_1_0(14, 1).
pf_1_0(15, 1).
pf_1(X) :- pf_1_0(X).
universe(X) :- pf_1(X).
pf_6(X) :- pf_6_0(X).
universe(X) :- pf_6(X).
pf_8_0(2, 1).
pf_8_0(3, 1).
pf_8_0(4, 1).
pf_8_0(5, 1).
pf_8_0(6, 1).
pf_8_0(7, 1).
pf_8_0(8, 1).
pf_8_0(9, 1).
pf_8_0(10, 1).
pf_8_0(11, 1).
pf_8_0(14, 1).
pf_8_0(15, 1).
pf_8(X) :- pf_8_0(X).
universe(X) :- pf_8(X).
pf_9_0(2, 1).
pf_9_0(3, 1).
pf_9_0(4, 1).
pf_9_0(5, 1).
pf_9_0(6, 1).
pf_9_0(7, 1).
pf_9_0(8, 1).
pf_9_0(9, 1).
pf_9_0(10, 1).
pf_9_0(11, 1).
pf_9_0(14, 1).
pf_9_0(15, 1).
pf_9(X) :- pf_9_0(X).
universe(X) :- pf_9(X).
) :- universe(A), pf_1(B), universe(C), universe(D), universe(E), universe(F), pf_6(G), universe(H), pf_8(I), pf_9(J), universe(K), universe(L), universe(M), universe(N), universe(O).
used_15(X,0) :- sequence_15(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,1) :- sequence_15(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,2) :- sequence_15(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,3) :- sequence_15(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,4) :- sequence_15(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_15(X,5) :- sequence_15(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_15(X,6) :- sequence_15(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_15(X,7) :- sequence_15(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_15(X,8) :- sequence_15(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_15(X,9) :- sequence_15(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_15(X,10) :- sequence_15(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_15(X,11) :- sequence_15(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_15(X,12) :- sequence_15(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_15(X,13) :- sequence_15(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_15(X,14) :- sequence_15(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.

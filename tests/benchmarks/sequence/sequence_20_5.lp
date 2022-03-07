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
universe_0(16, 1).
universe_0(17, 1).
universe_0(18, 1).
universe_0(19, 1).
universe_0(20, 1).
universe(X) :- universe_0(X).
universe(X) :- universe(X).
dom1_0(1, 1).
dom1_0(2, 1).
dom1_0(7, 1).
dom1_0(8, 1).
dom1_0(9, 1).
dom1_0(10, 1).
dom1_0(11, 1).
dom1_0(12, 1).
dom1_0(13, 1).
dom1_0(14, 1).
dom1_0(15, 1).
dom1_0(16, 1).
dom1_0(17, 1).
dom1_0(18, 1).
dom1_0(19, 1).
dom1_0(20, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
dom2_0(1, 1).
dom2_0(2, 1).
dom2_0(3, 1).
dom2_0(4, 1).
dom2_0(5, 1).
dom2_0(6, 1).
dom2_0(7, 1).
dom2_0(8, 1).
dom2_0(9, 1).
dom2_0(10, 1).
dom2_0(11, 1).
dom2_0(12, 1).
dom2_0(13, 1).
dom2_0(14, 1).
dom2_0(15, 1).
dom2_0(16, 1).
dom2_0(17, 1).
dom2_0(18, 1).
dom2_0(19, 1).
dom2_0(20, 1).
dom2(X) :- dom2_0(X).
universe(X) :- dom2(X).
dom3_0(14, 1).
dom3_0(15, 1).
dom3_0(16, 1).
dom3_0(17, 1).
dom3_0(18, 1).
dom3_0(19, 1).
dom3_0(20, 1).
dom3(X) :- dom3_0(X).
universe(X) :- dom3(X).
sequence_1(A) :- universe(A).
used_1(X,0) :- sequence_1(X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
sequence_2(A,B) :- universe(A), universe(B).
used_2(X,0) :- sequence_2(X, _). 
used_2(X,1) :- sequence_2(_, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
sequence_3(A,B,C) :- universe(A), universe(B), universe(C).
used_3(X,0) :- sequence_3(X, _, _). 
used_3(X,1) :- sequence_3(_, X, _). 
used_3(X,2) :- sequence_3(_, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
sequence_4(A,B,C,D) :- universe(A), universe(B), universe(C), universe(D).
used_4(X,0) :- sequence_4(X, _, _, _). 
used_4(X,1) :- sequence_4(_, X, _, _). 
used_4(X,2) :- sequence_4(_, _, X, _). 
used_4(X,3) :- sequence_4(_, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
sequence_6(A,B,C,D,E,F) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F).
used_6(X,0) :- sequence_6(X, _, _, _, _, _). 
used_6(X,1) :- sequence_6(_, X, _, _, _, _). 
used_6(X,2) :- sequence_6(_, _, X, _, _, _). 
used_6(X,3) :- sequence_6(_, _, _, X, _, _). 
used_6(X,4) :- sequence_6(_, _, _, _, X, _). 
used_6(X,5) :- sequence_6(_, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
sequence_7(A,B,C,D,E,F,G) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G).
used_7(X,0) :- sequence_7(X, _, _, _, _, _, _). 
used_7(X,1) :- sequence_7(_, X, _, _, _, _, _). 
used_7(X,2) :- sequence_7(_, _, X, _, _, _, _). 
used_7(X,3) :- sequence_7(_, _, _, X, _, _, _). 
used_7(X,4) :- sequence_7(_, _, _, _, X, _, _). 
used_7(X,5) :- sequence_7(_, _, _, _, _, X, _). 
used_7(X,6) :- sequence_7(_, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
sequence_8(A,B,C,D,E,F,G,H) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H).
used_8(X,0) :- sequence_8(X, _, _, _, _, _, _, _). 
used_8(X,1) :- sequence_8(_, X, _, _, _, _, _, _). 
used_8(X,2) :- sequence_8(_, _, X, _, _, _, _, _). 
used_8(X,3) :- sequence_8(_, _, _, X, _, _, _, _). 
used_8(X,4) :- sequence_8(_, _, _, _, X, _, _, _). 
used_8(X,5) :- sequence_8(_, _, _, _, _, X, _, _). 
used_8(X,6) :- sequence_8(_, _, _, _, _, _, X, _). 
used_8(X,7) :- sequence_8(_, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
sequence_9(A,B,C,D,E,F,G,H,I) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I).
used_9(X,0) :- sequence_9(X, _, _, _, _, _, _, _, _). 
used_9(X,1) :- sequence_9(_, X, _, _, _, _, _, _, _). 
used_9(X,2) :- sequence_9(_, _, X, _, _, _, _, _, _). 
used_9(X,3) :- sequence_9(_, _, _, X, _, _, _, _, _). 
used_9(X,4) :- sequence_9(_, _, _, _, X, _, _, _, _). 
used_9(X,5) :- sequence_9(_, _, _, _, _, X, _, _, _). 
used_9(X,6) :- sequence_9(_, _, _, _, _, _, X, _, _). 
used_9(X,7) :- sequence_9(_, _, _, _, _, _, _, X, _). 
used_9(X,8) :- sequence_9(_, _, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
sequence_10(A,B,C,D,E,F,G,H,I,J) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J).
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
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
sequence_11(A,B,C,D,E,F,G,H,I,J,K) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K).
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
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
sequence_12(A,B,C,D,E,F,G,H,I,J,K,L) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L).
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
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
sequence_13(A,B,C,D,E,F,G,H,I,J,K,L,M) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M).
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
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
:- C = #count{N:used(S,N),df_0(S)}, C!=13.
sequence_14(A,B,C,D,E,F,G,H,I,J,K,L,M,N) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N).
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
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
:- C = #count{N:used(S,N),df_0(S)}, C!=13.
:- C = #count{N:used(S,N),df_0(S)}, C!=14.
sequence_15(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O).
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
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
:- C = #count{N:used(S,N),df_0(S)}, C!=13.
:- C = #count{N:used(S,N),df_0(S)}, C!=14.
:- C = #count{N:used(S,N),df_0(S)}, C!=15.
sequence_16(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P).
used_16(X,0) :- sequence_16(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,1) :- sequence_16(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,2) :- sequence_16(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,3) :- sequence_16(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,4) :- sequence_16(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,5) :- sequence_16(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_16(X,6) :- sequence_16(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_16(X,7) :- sequence_16(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_16(X,8) :- sequence_16(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_16(X,9) :- sequence_16(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_16(X,10) :- sequence_16(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_16(X,11) :- sequence_16(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_16(X,12) :- sequence_16(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_16(X,13) :- sequence_16(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_16(X,14) :- sequence_16(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_16(X,15) :- sequence_16(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
:- C = #count{N:used(S,N),df_0(S)}, C!=13.
:- C = #count{N:used(S,N),df_0(S)}, C!=14.
:- C = #count{N:used(S,N),df_0(S)}, C!=15.
:- C = #count{N:used(S,N),df_0(S)}, C!=16.
sequence_17(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q).
used_17(X,0) :- sequence_17(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,1) :- sequence_17(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,2) :- sequence_17(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,3) :- sequence_17(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,4) :- sequence_17(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,5) :- sequence_17(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,6) :- sequence_17(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_17(X,7) :- sequence_17(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_17(X,8) :- sequence_17(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_17(X,9) :- sequence_17(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_17(X,10) :- sequence_17(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_17(X,11) :- sequence_17(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_17(X,12) :- sequence_17(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_17(X,13) :- sequence_17(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_17(X,14) :- sequence_17(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_17(X,15) :- sequence_17(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_17(X,16) :- sequence_17(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
:- C = #count{N:used(S,N),df_0(S)}, C!=13.
:- C = #count{N:used(S,N),df_0(S)}, C!=14.
:- C = #count{N:used(S,N),df_0(S)}, C!=15.
:- C = #count{N:used(S,N),df_0(S)}, C!=16.
:- C = #count{N:used(S,N),df_0(S)}, C!=17.
sequence_18(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q), universe(R).
used_18(X,0) :- sequence_18(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,1) :- sequence_18(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,2) :- sequence_18(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,3) :- sequence_18(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,4) :- sequence_18(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,5) :- sequence_18(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,6) :- sequence_18(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,7) :- sequence_18(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_18(X,8) :- sequence_18(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_18(X,9) :- sequence_18(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_18(X,10) :- sequence_18(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_18(X,11) :- sequence_18(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_18(X,12) :- sequence_18(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_18(X,13) :- sequence_18(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_18(X,14) :- sequence_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_18(X,15) :- sequence_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_18(X,16) :- sequence_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_18(X,17) :- sequence_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
:- C = #count{N:used(S,N),df_0(S)}, C!=13.
:- C = #count{N:used(S,N),df_0(S)}, C!=14.
:- C = #count{N:used(S,N),df_0(S)}, C!=15.
:- C = #count{N:used(S,N),df_0(S)}, C!=16.
:- C = #count{N:used(S,N),df_0(S)}, C!=17.
:- C = #count{N:used(S,N),df_0(S)}, C!=18.
sequence_19(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q), universe(R), universe(S).
used_19(X,0) :- sequence_19(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,1) :- sequence_19(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,2) :- sequence_19(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,3) :- sequence_19(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,4) :- sequence_19(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,5) :- sequence_19(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,6) :- sequence_19(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,7) :- sequence_19(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,8) :- sequence_19(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_19(X,9) :- sequence_19(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_19(X,10) :- sequence_19(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_19(X,11) :- sequence_19(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_19(X,12) :- sequence_19(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_19(X,13) :- sequence_19(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_19(X,14) :- sequence_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_19(X,15) :- sequence_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_19(X,16) :- sequence_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_19(X,17) :- sequence_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_19(X,18) :- sequence_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
:- C = #count{N:used(S,N),df_0(S)}, C!=13.
:- C = #count{N:used(S,N),df_0(S)}, C!=14.
:- C = #count{N:used(S,N),df_0(S)}, C!=15.
:- C = #count{N:used(S,N),df_0(S)}, C!=16.
:- C = #count{N:used(S,N),df_0(S)}, C!=17.
:- C = #count{N:used(S,N),df_0(S)}, C!=18.
:- C = #count{N:used(S,N),df_0(S)}, C!=19.
sequence_20(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q), universe(R), universe(S), universe(T).
used_20(X,0) :- sequence_20(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,1) :- sequence_20(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,2) :- sequence_20(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,3) :- sequence_20(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,4) :- sequence_20(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,5) :- sequence_20(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,6) :- sequence_20(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,7) :- sequence_20(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,8) :- sequence_20(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,9) :- sequence_20(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_20(X,10) :- sequence_20(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_20(X,11) :- sequence_20(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_20(X,12) :- sequence_20(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_20(X,13) :- sequence_20(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_20(X,14) :- sequence_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_20(X,15) :- sequence_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_20(X,16) :- sequence_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_20(X,17) :- sequence_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_20(X,18) :- sequence_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_20(X,19) :- sequence_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
:- C = #count{N:used(S,N),df_0(S)}, C!=13.
:- C = #count{N:used(S,N),df_0(S)}, C!=14.
:- C = #count{N:used(S,N),df_0(S)}, C!=15.
:- C = #count{N:used(S,N),df_0(S)}, C!=16.
:- C = #count{N:used(S,N),df_0(S)}, C!=17.
:- C = #count{N:used(S,N),df_0(S)}, C!=18.
:- C = #count{N:used(S,N),df_0(S)}, C!=19.
:- C = #count{N:used(S,N),df_0(S)}, C!=20.

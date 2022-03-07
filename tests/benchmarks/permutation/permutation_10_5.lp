universe_0(1, 1).
universe_0(2, 1).
universe_0(3, 1).
universe_0(4, 1).
universe_0(5, 1).
universe(X) :- universe_0(X).
universe_1(6,5).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(1, 1).
dom1_0(2, 1).
dom1_0(3, 1).
dom1_0(4, 1).
dom1_0(5, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
permutation_1(A) :- universe(A).
used_1(X,0) :- permutation_1(X). 
:- permutation_1(A), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
permutation_2(A,B) :- universe(A), universe(B).
used_2(X,0) :- permutation_2(X, _). 
used_2(X,1) :- permutation_2(_, X). 
:- permutation_2(A,B), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
permutation_3(A,B,C) :- universe(A), universe(B), universe(C).
used_3(X,0) :- permutation_3(X, _, _). 
used_3(X,1) :- permutation_3(_, X, _). 
used_3(X,2) :- permutation_3(_, _, X). 
:- permutation_3(A,B,C), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
permutation_4(A,B,C,D) :- universe(A), universe(B), universe(C), universe(D).
used_4(X,0) :- permutation_4(X, _, _, _). 
used_4(X,1) :- permutation_4(_, X, _, _). 
used_4(X,2) :- permutation_4(_, _, X, _). 
used_4(X,3) :- permutation_4(_, _, _, X). 
:- permutation_4(A,B,C,D), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
permutation_6(A,B,C,D,E,F) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F).
used_6(X,0) :- permutation_6(X, _, _, _, _, _). 
used_6(X,1) :- permutation_6(_, X, _, _, _, _). 
used_6(X,2) :- permutation_6(_, _, X, _, _, _). 
used_6(X,3) :- permutation_6(_, _, _, X, _, _). 
used_6(X,4) :- permutation_6(_, _, _, _, X, _). 
used_6(X,5) :- permutation_6(_, _, _, _, _, X). 
:- permutation_6(A,B,C,D,E,F), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
permutation_7(A,B,C,D,E,F,G) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G).
used_7(X,0) :- permutation_7(X, _, _, _, _, _, _). 
used_7(X,1) :- permutation_7(_, X, _, _, _, _, _). 
used_7(X,2) :- permutation_7(_, _, X, _, _, _, _). 
used_7(X,3) :- permutation_7(_, _, _, X, _, _, _). 
used_7(X,4) :- permutation_7(_, _, _, _, X, _, _). 
used_7(X,5) :- permutation_7(_, _, _, _, _, X, _). 
used_7(X,6) :- permutation_7(_, _, _, _, _, _, X). 
:- permutation_7(A,B,C,D,E,F,G), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
permutation_8(A,B,C,D,E,F,G,H) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H).
used_8(X,0) :- permutation_8(X, _, _, _, _, _, _, _). 
used_8(X,1) :- permutation_8(_, X, _, _, _, _, _, _). 
used_8(X,2) :- permutation_8(_, _, X, _, _, _, _, _). 
used_8(X,3) :- permutation_8(_, _, _, X, _, _, _, _). 
used_8(X,4) :- permutation_8(_, _, _, _, X, _, _, _). 
used_8(X,5) :- permutation_8(_, _, _, _, _, X, _, _). 
used_8(X,6) :- permutation_8(_, _, _, _, _, _, X, _). 
used_8(X,7) :- permutation_8(_, _, _, _, _, _, _, X). 
:- permutation_8(A,B,C,D,E,F,G,H), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
permutation_9(A,B,C,D,E,F,G,H,I) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I).
used_9(X,0) :- permutation_9(X, _, _, _, _, _, _, _, _). 
used_9(X,1) :- permutation_9(_, X, _, _, _, _, _, _, _). 
used_9(X,2) :- permutation_9(_, _, X, _, _, _, _, _, _). 
used_9(X,3) :- permutation_9(_, _, _, X, _, _, _, _, _). 
used_9(X,4) :- permutation_9(_, _, _, _, X, _, _, _, _). 
used_9(X,5) :- permutation_9(_, _, _, _, _, X, _, _, _). 
used_9(X,6) :- permutation_9(_, _, _, _, _, _, X, _, _). 
used_9(X,7) :- permutation_9(_, _, _, _, _, _, _, X, _). 
used_9(X,8) :- permutation_9(_, _, _, _, _, _, _, _, X). 
:- permutation_9(A,B,C,D,E,F,G,H,I), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
permutation_10(A,B,C,D,E,F,G,H,I,J) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J).
used_10(X,0) :- permutation_10(X, _, _, _, _, _, _, _, _, _). 
used_10(X,1) :- permutation_10(_, X, _, _, _, _, _, _, _, _). 
used_10(X,2) :- permutation_10(_, _, X, _, _, _, _, _, _, _). 
used_10(X,3) :- permutation_10(_, _, _, X, _, _, _, _, _, _). 
used_10(X,4) :- permutation_10(_, _, _, _, X, _, _, _, _, _). 
used_10(X,5) :- permutation_10(_, _, _, _, _, X, _, _, _, _). 
used_10(X,6) :- permutation_10(_, _, _, _, _, _, X, _, _, _). 
used_10(X,7) :- permutation_10(_, _, _, _, _, _, _, X, _, _). 
used_10(X,8) :- permutation_10(_, _, _, _, _, _, _, _, X, _). 
used_10(X,9) :- permutation_10(_, _, _, _, _, _, _, _, _, X). 
:- permutation_10(A,B,C,D,E,F,G,H,I,J), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(4, 1).
df_0_0(5, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.

universe_0(1, 1).
universe_0(2, 1).
universe_0(3, 1).
universe_0(4, 1).
universe_0(5, 1).
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
universe_1(6,5).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(1, 1).
dom1_0(2, 1).
dom1_0(3, 1).
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
permutation_1(A) :- universe(A).
used_1(X,0) :- permutation_1(X). 
:- permutation_1(A), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
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
df_1_0(1, 1).
df_1_0(2, 1).
df_1_0(3, 1).
df_1_0(11, 1).
df_1_0(12, 1).
df_1_0(13, 1).
df_1_0(14, 1).
df_1_0(15, 1).
df_1_0(16, 1).
df_1_0(17, 1).
df_1_0(18, 1).
df_1_0(19, 1).
df_1_0(20, 1).
df_1(X) :- df_1_0(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=0.
:- C = #count{N:used(S,N),df_1(S)}, C!=1.
:- C = #count{N:used(S,N),df_1(S)}, C!=2.
:- C = #count{N:used(S,N),df_1(S)}, C!=3.
:- C = #count{N:used(S,N),df_1(S)}, C!=4.
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
:- C = #count{N:used(S,N),df_1(S)}, C!=6.
:- C = #count{N:used(S,N),df_1(S)}, C!=7.
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.
permutation_2(A,B) :- universe(A), universe(B).
used_2(X,0) :- permutation_2(X, _). 
used_2(X,1) :- permutation_2(_, X). 
:- permutation_2(A,B), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
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
df_1_0(1, 1).
df_1_0(2, 1).
df_1_0(3, 1).
df_1_0(11, 1).
df_1_0(12, 1).
df_1_0(13, 1).
df_1_0(14, 1).
df_1_0(15, 1).
df_1_0(16, 1).
df_1_0(17, 1).
df_1_0(18, 1).
df_1_0(19, 1).
df_1_0(20, 1).
df_1(X) :- df_1_0(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=0.
:- C = #count{N:used(S,N),df_1(S)}, C!=1.
:- C = #count{N:used(S,N),df_1(S)}, C!=2.
:- C = #count{N:used(S,N),df_1(S)}, C!=3.
:- C = #count{N:used(S,N),df_1(S)}, C!=4.
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
:- C = #count{N:used(S,N),df_1(S)}, C!=6.
:- C = #count{N:used(S,N),df_1(S)}, C!=7.
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.
permutation_3(A,B,C) :- universe(A), universe(B), universe(C).
used_3(X,0) :- permutation_3(X, _, _). 
used_3(X,1) :- permutation_3(_, X, _). 
used_3(X,2) :- permutation_3(_, _, X). 
:- permutation_3(A,B,C), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
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
df_1_0(1, 1).
df_1_0(2, 1).
df_1_0(3, 1).
df_1_0(11, 1).
df_1_0(12, 1).
df_1_0(13, 1).
df_1_0(14, 1).
df_1_0(15, 1).
df_1_0(16, 1).
df_1_0(17, 1).
df_1_0(18, 1).
df_1_0(19, 1).
df_1_0(20, 1).
df_1(X) :- df_1_0(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=0.
:- C = #count{N:used(S,N),df_1(S)}, C!=1.
:- C = #count{N:used(S,N),df_1(S)}, C!=2.
:- C = #count{N:used(S,N),df_1(S)}, C!=3.
:- C = #count{N:used(S,N),df_1(S)}, C!=4.
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
:- C = #count{N:used(S,N),df_1(S)}, C!=6.
:- C = #count{N:used(S,N),df_1(S)}, C!=7.
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.
permutation_4(A,B,C,Dpf_3_0(1, 1).
pf_3_0(2, 1).
pf_3_0(3, 1).
pf_3_0(11, 1).
pf_3_0(12, 1).
pf_3_0(13, 1).
pf_3_0(14, 1).
pf_3_0(15, 1).
pf_3_0(16, 1).
pf_3_0(17, 1).
pf_3_0(18, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
) :- universe(A), universe(B), universe(C), pf_3(D).
used_4(X,0) :- permutation_4(X, _, _, _). 
used_4(X,1) :- permutation_4(_, X, _, _). 
used_4(X,2) :- permutation_4(_, _, X, _). 
used_4(X,3) :- permutation_4(_, _, _, X). 
:- permutation_4(A,B,C,D), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
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
df_1_0(1, 1).
df_1_0(2, 1).
df_1_0(3, 1).
df_1_0(11, 1).
df_1_0(12, 1).
df_1_0(13, 1).
df_1_0(14, 1).
df_1_0(15, 1).
df_1_0(16, 1).
df_1_0(17, 1).
df_1_0(18, 1).
df_1_0(19, 1).
df_1_0(20, 1).
df_1(X) :- df_1_0(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=0.
:- C = #count{N:used(S,N),df_1(S)}, C!=1.
:- C = #count{N:used(S,N),df_1(S)}, C!=2.
:- C = #count{N:used(S,N),df_1(S)}, C!=3.
:- C = #count{N:used(S,N),df_1(S)}, C!=4.
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
:- C = #count{N:used(S,N),df_1(S)}, C!=6.
:- C = #count{N:used(S,N),df_1(S)}, C!=7.
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.
permutation_5(A,B,C,D,Epf_3_0(1, 1).
pf_3_0(2, 1).
pf_3_0(3, 1).
pf_3_0(11, 1).
pf_3_0(12, 1).
pf_3_0(13, 1).
pf_3_0(14, 1).
pf_3_0(15, 1).
pf_3_0(16, 1).
pf_3_0(17, 1).
pf_3_0(18, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
) :- universe(A), universe(B), universe(C), pf_3(D), universe(E).
used_5(X,0) :- permutation_5(X, _, _, _, _). 
used_5(X,1) :- permutation_5(_, X, _, _, _). 
used_5(X,2) :- permutation_5(_, _, X, _, _). 
used_5(X,3) :- permutation_5(_, _, _, X, _). 
used_5(X,4) :- permutation_5(_, _, _, _, X). 
:- permutation_5(A,B,C,D,E), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
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
df_1_0(1, 1).
df_1_0(2, 1).
df_1_0(3, 1).
df_1_0(11, 1).
df_1_0(12, 1).
df_1_0(13, 1).
df_1_0(14, 1).
df_1_0(15, 1).
df_1_0(16, 1).
df_1_0(17, 1).
df_1_0(18, 1).
df_1_0(19, 1).
df_1_0(20, 1).
df_1(X) :- df_1_0(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=0.
:- C = #count{N:used(S,N),df_1(S)}, C!=1.
:- C = #count{N:used(S,N),df_1(S)}, C!=2.
:- C = #count{N:used(S,N),df_1(S)}, C!=3.
:- C = #count{N:used(S,N),df_1(S)}, C!=4.
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
:- C = #count{N:used(S,N),df_1(S)}, C!=6.
:- C = #count{N:used(S,N),df_1(S)}, C!=7.
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.
permutation_6(A,B,C,D,E,Fpf_3_0(1, 1).
pf_3_0(2, 1).
pf_3_0(3, 1).
pf_3_0(11, 1).
pf_3_0(12, 1).
pf_3_0(13, 1).
pf_3_0(14, 1).
pf_3_0(15, 1).
pf_3_0(16, 1).
pf_3_0(17, 1).
pf_3_0(18, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
) :- universe(A), universe(B), universe(C), pf_3(D), universe(E), universe(F).
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
df_1_0(1, 1).
df_1_0(2, 1).
df_1_0(3, 1).
df_1_0(11, 1).
df_1_0(12, 1).
df_1_0(13, 1).
df_1_0(14, 1).
df_1_0(15, 1).
df_1_0(16, 1).
df_1_0(17, 1).
df_1_0(18, 1).
df_1_0(19, 1).
df_1_0(20, 1).
df_1(X) :- df_1_0(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=0.
:- C = #count{N:used(S,N),df_1(S)}, C!=1.
:- C = #count{N:used(S,N),df_1(S)}, C!=2.
:- C = #count{N:used(S,N),df_1(S)}, C!=3.
:- C = #count{N:used(S,N),df_1(S)}, C!=4.
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
:- C = #count{N:used(S,N),df_1(S)}, C!=6.
:- C = #count{N:used(S,N),df_1(S)}, C!=7.
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.
permutation_7(A,B,C,D,E,F,Gpf_3_0(1, 1).
pf_3_0(2, 1).
pf_3_0(3, 1).
pf_3_0(11, 1).
pf_3_0(12, 1).
pf_3_0(13, 1).
pf_3_0(14, 1).
pf_3_0(15, 1).
pf_3_0(16, 1).
pf_3_0(17, 1).
pf_3_0(18, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
pf_6_0(1, 1).
pf_6_0(2, 1).
pf_6_0(3, 1).
pf_6_0(11, 1).
pf_6_0(12, 1).
pf_6_0(13, 1).
pf_6_0(14, 1).
pf_6_0(15, 1).
pf_6_0(16, 1).
pf_6_0(17, 1).
pf_6_0(18, 1).
pf_6_0(19, 1).
pf_6_0(20, 1).
pf_6(X) :- pf_6_0(X).
universe(X) :- pf_6(X).
) :- universe(A), universe(B), universe(C), pf_3(D), universe(E), universe(F), pf_6(G).
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
df_1_0(1, 1).
df_1_0(2, 1).
df_1_0(3, 1).
df_1_0(11, 1).
df_1_0(12, 1).
df_1_0(13, 1).
df_1_0(14, 1).
df_1_0(15, 1).
df_1_0(16, 1).
df_1_0(17, 1).
df_1_0(18, 1).
df_1_0(19, 1).
df_1_0(20, 1).
df_1(X) :- df_1_0(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=0.
:- C = #count{N:used(S,N),df_1(S)}, C!=1.
:- C = #count{N:used(S,N),df_1(S)}, C!=2.
:- C = #count{N:used(S,N),df_1(S)}, C!=3.
:- C = #count{N:used(S,N),df_1(S)}, C!=4.
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
:- C = #count{N:used(S,N),df_1(S)}, C!=6.
:- C = #count{N:used(S,N),df_1(S)}, C!=7.
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.
permutation_8(A,B,C,D,E,F,G,Hpf_3_0(1, 1).
pf_3_0(2, 1).
pf_3_0(3, 1).
pf_3_0(11, 1).
pf_3_0(12, 1).
pf_3_0(13, 1).
pf_3_0(14, 1).
pf_3_0(15, 1).
pf_3_0(16, 1).
pf_3_0(17, 1).
pf_3_0(18, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
pf_6_0(1, 1).
pf_6_0(2, 1).
pf_6_0(3, 1).
pf_6_0(11, 1).
pf_6_0(12, 1).
pf_6_0(13, 1).
pf_6_0(14, 1).
pf_6_0(15, 1).
pf_6_0(16, 1).
pf_6_0(17, 1).
pf_6_0(18, 1).
pf_6_0(19, 1).
pf_6_0(20, 1).
pf_6(X) :- pf_6_0(X).
universe(X) :- pf_6(X).
) :- universe(A), universe(B), universe(C), pf_3(D), universe(E), universe(F), pf_6(G), universe(H).
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
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
df_1_0(1, 1).
df_1_0(2, 1).
df_1_0(3, 1).
df_1_0(11, 1).
df_1_0(12, 1).
df_1_0(13, 1).
df_1_0(14, 1).
df_1_0(15, 1).
df_1_0(16, 1).
df_1_0(17, 1).
df_1_0(18, 1).
df_1_0(19, 1).
df_1_0(20, 1).
df_1(X) :- df_1_0(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=0.
:- C = #count{N:used(S,N),df_1(S)}, C!=1.
:- C = #count{N:used(S,N),df_1(S)}, C!=2.
:- C = #count{N:used(S,N),df_1(S)}, C!=3.
:- C = #count{N:used(S,N),df_1(S)}, C!=4.
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
:- C = #count{N:used(S,N),df_1(S)}, C!=6.
:- C = #count{N:used(S,N),df_1(S)}, C!=7.
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.
permutation_9(A,B,C,D,E,F,G,H,Ipf_3_0(1, 1).
pf_3_0(2, 1).
pf_3_0(3, 1).
pf_3_0(11, 1).
pf_3_0(12, 1).
pf_3_0(13, 1).
pf_3_0(14, 1).
pf_3_0(15, 1).
pf_3_0(16, 1).
pf_3_0(17, 1).
pf_3_0(18, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
pf_6_0(1, 1).
pf_6_0(2, 1).
pf_6_0(3, 1).
pf_6_0(11, 1).
pf_6_0(12, 1).
pf_6_0(13, 1).
pf_6_0(14, 1).
pf_6_0(15, 1).
pf_6_0(16, 1).
pf_6_0(17, 1).
pf_6_0(18, 1).
pf_6_0(19, 1).
pf_6_0(20, 1).
pf_6(X) :- pf_6_0(X).
universe(X) :- pf_6(X).
) :- universe(A), universe(B), universe(C), pf_3(D), universe(E), universe(F), pf_6(G), universe(H), universe(I).
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
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
df_1_0(1, 1).
df_1_0(2, 1).
df_1_0(3, 1).
df_1_0(11, 1).
df_1_0(12, 1).
df_1_0(13, 1).
df_1_0(14, 1).
df_1_0(15, 1).
df_1_0(16, 1).
df_1_0(17, 1).
df_1_0(18, 1).
df_1_0(19, 1).
df_1_0(20, 1).
df_1(X) :- df_1_0(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=0.
:- C = #count{N:used(S,N),df_1(S)}, C!=1.
:- C = #count{N:used(S,N),df_1(S)}, C!=2.
:- C = #count{N:used(S,N),df_1(S)}, C!=3.
:- C = #count{N:used(S,N),df_1(S)}, C!=4.
:- C = #count{N:used(S,N),df_1(S)}, C!=5.
:- C = #count{N:used(S,N),df_1(S)}, C!=6.
:- C = #count{N:used(S,N),df_1(S)}, C!=7.
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.

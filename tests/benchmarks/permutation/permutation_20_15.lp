universe_0(1, 1).
universe_0(2, 1).
universe_0(19, 1).
universe_0(20, 1).
universe(X) :- universe_0(X).
universe_1(3,16).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(2, 1).
dom1_0(19, 1).
dom1_0(20, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
permutation_1(A) :- universe(A).
used_1(X,0) :- permutation_1(X). 
:- permutation_1(A), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
df_1_0(2, 1).
df_1_0(19, 1).
df_1_0(20, 1).
df_1(X) :- df_1_0(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=0.
:- C = #count{N:used(S,N),df_1(S)}, C!=1.
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
permutation_2(A,B) :- universe(A), universe(B).
used_2(X,0) :- permutation_2(X, _). 
used_2(X,1) :- permutation_2(_, X). 
:- permutation_2(A,B), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
df_1_0(2, 1).
df_1_0(19, 1).
df_1_0(20, 1).
df_1(X) :- df_1_0(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=0.
:- C = #count{N:used(S,N),df_1(S)}, C!=1.
:- C = #count{N:used(S,N),df_1(S)}, C!=2.
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
permutation_3(A,B,C) :- universe(A), universe(B), universe(C).
used_3(X,0) :- permutation_3(X, _, _). 
used_3(X,1) :- permutation_3(_, X, _). 
used_3(X,2) :- permutation_3(_, _, X). 
:- permutation_3(A,B,C), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
df_1_0(2, 1).
df_1_0(19, 1).
df_1_0(20, 1).
df_1(X) :- df_1_0(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=0.
:- C = #count{N:used(S,N),df_1(S)}, C!=1.
:- C = #count{N:used(S,N),df_1(S)}, C!=2.
:- C = #count{N:used(S,N),df_1(S)}, C!=3.
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
permutation_4(A,B,C,Dpf_3_0(2, 1).
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
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
df_1_0(2, 1).
df_1_0(19, 1).
df_1_0(20, 1).
df_1(X) :- df_1_0(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=0.
:- C = #count{N:used(S,N),df_1(S)}, C!=1.
:- C = #count{N:used(S,N),df_1(S)}, C!=2.
:- C = #count{N:used(S,N),df_1(S)}, C!=3.
:- C = #count{N:used(S,N),df_1(S)}, C!=4.
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
permutation_5(A,B,C,D,Epf_3_0(2, 1).
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
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
df_1_0(2, 1).
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
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
permutation_6(A,B,C,D,E,Fpf_3_0(2, 1).
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
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
df_1_0(2, 1).
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
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
:- C = #count{N:used(S,N),df_4(S)}, C!=6.
permutation_7(A,B,C,D,E,F,Gpf_3_0(2, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
) :- universe(A), universe(B), universe(C), pf_3(D), universe(E), universe(F), universe(G).
used_7(X,0) :- permutation_7(X, _, _, _, _, _, _). 
used_7(X,1) :- permutation_7(_, X, _, _, _, _, _). 
used_7(X,2) :- permutation_7(_, _, X, _, _, _, _). 
used_7(X,3) :- permutation_7(_, _, _, X, _, _, _). 
used_7(X,4) :- permutation_7(_, _, _, _, X, _, _). 
used_7(X,5) :- permutation_7(_, _, _, _, _, X, _). 
used_7(X,6) :- permutation_7(_, _, _, _, _, _, X). 
:- permutation_7(A,B,C,D,E,F,G), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
df_1_0(2, 1).
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
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
:- C = #count{N:used(S,N),df_4(S)}, C!=6.
:- C = #count{N:used(S,N),df_4(S)}, C!=7.
permutation_8(A,B,C,D,E,F,G,Hpf_3_0(2, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
pf_7_0(2, 1).
pf_7_0(19, 1).
pf_7_0(20, 1).
pf_7(X) :- pf_7_0(X).
universe(X) :- pf_7(X).
) :- universe(A), universe(B), universe(C), pf_3(D), universe(E), universe(F), universe(G), pf_7(H).
used_8(X,0) :- permutation_8(X, _, _, _, _, _, _, _). 
used_8(X,1) :- permutation_8(_, X, _, _, _, _, _, _). 
used_8(X,2) :- permutation_8(_, _, X, _, _, _, _, _). 
used_8(X,3) :- permutation_8(_, _, _, X, _, _, _, _). 
used_8(X,4) :- permutation_8(_, _, _, _, X, _, _, _). 
used_8(X,5) :- permutation_8(_, _, _, _, _, X, _, _). 
used_8(X,6) :- permutation_8(_, _, _, _, _, _, X, _). 
used_8(X,7) :- permutation_8(_, _, _, _, _, _, _, X). 
:- permutation_8(A,B,C,D,E,F,G,H), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
df_1_0(2, 1).
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
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
:- C = #count{N:used(S,N),df_4(S)}, C!=6.
:- C = #count{N:used(S,N),df_4(S)}, C!=7.
:- C = #count{N:used(S,N),df_4(S)}, C!=8.
permutation_9(A,B,C,D,E,F,G,H,Ipf_3_0(2, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
pf_7_0(2, 1).
pf_7_0(19, 1).
pf_7_0(20, 1).
pf_7(X) :- pf_7_0(X).
universe(X) :- pf_7(X).
) :- universe(A), universe(B), universe(C), pf_3(D), universe(E), universe(F), universe(G), pf_7(H), universe(I).
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
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
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
df_1_0(2, 1).
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
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
:- C = #count{N:used(S,N),df_4(S)}, C!=6.
:- C = #count{N:used(S,N),df_4(S)}, C!=7.
:- C = #count{N:used(S,N),df_4(S)}, C!=8.
:- C = #count{N:used(S,N),df_4(S)}, C!=9.
permutation_10(A,B,C,D,E,F,G,H,I,Jpf_3_0(2, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
pf_7_0(2, 1).
pf_7_0(19, 1).
pf_7_0(20, 1).
pf_7(X) :- pf_7_0(X).
universe(X) :- pf_7(X).
pf_9_0(2, 1).
pf_9_0(19, 1).
pf_9_0(20, 1).
pf_9(X) :- pf_9_0(X).
universe(X) :- pf_9(X).
) :- universe(A), universe(B), universe(C), pf_3(D), universe(E), universe(F), universe(G), pf_7(H), universe(I), pf_9(J).
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
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
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
df_1_0(2, 1).
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
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
:- C = #count{N:used(S,N),df_4(S)}, C!=6.
:- C = #count{N:used(S,N),df_4(S)}, C!=7.
:- C = #count{N:used(S,N),df_4(S)}, C!=8.
:- C = #count{N:used(S,N),df_4(S)}, C!=9.
:- C = #count{N:used(S,N),df_4(S)}, C!=10.
permutation_11(A,B,C,D,E,F,G,H,I,J,Kpf_3_0(2, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
pf_7_0(2, 1).
pf_7_0(19, 1).
pf_7_0(20, 1).
pf_7(X) :- pf_7_0(X).
universe(X) :- pf_7(X).
pf_9_0(2, 1).
pf_9_0(19, 1).
pf_9_0(20, 1).
pf_9(X) :- pf_9_0(X).
universe(X) :- pf_9(X).
) :- universe(A), universe(B), universe(C), pf_3(D), universe(E), universe(F), universe(G), pf_7(H), universe(I), pf_9(J), universe(K).
used_11(X,0) :- permutation_11(X, _, _, _, _, _, _, _, _, _, _). 
used_11(X,1) :- permutation_11(_, X, _, _, _, _, _, _, _, _, _). 
used_11(X,2) :- permutation_11(_, _, X, _, _, _, _, _, _, _, _). 
used_11(X,3) :- permutation_11(_, _, _, X, _, _, _, _, _, _, _). 
used_11(X,4) :- permutation_11(_, _, _, _, X, _, _, _, _, _, _). 
used_11(X,5) :- permutation_11(_, _, _, _, _, X, _, _, _, _, _). 
used_11(X,6) :- permutation_11(_, _, _, _, _, _, X, _, _, _, _). 
used_11(X,7) :- permutation_11(_, _, _, _, _, _, _, X, _, _, _). 
used_11(X,8) :- permutation_11(_, _, _, _, _, _, _, _, X, _, _). 
used_11(X,9) :- permutation_11(_, _, _, _, _, _, _, _, _, X, _). 
used_11(X,10) :- permutation_11(_, _, _, _, _, _, _, _, _, _, X). 
:- permutation_11(A,B,C,D,E,F,G,H,I,J,K), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
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
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
df_1_0(2, 1).
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
:- C = #count{N:used(S,N),df_1(S)}, C!=11.
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
:- C = #count{N:used(S,N),df_4(S)}, C!=6.
:- C = #count{N:used(S,N),df_4(S)}, C!=7.
:- C = #count{N:used(S,N),df_4(S)}, C!=8.
:- C = #count{N:used(S,N),df_4(S)}, C!=9.
:- C = #count{N:used(S,N),df_4(S)}, C!=10.
:- C = #count{N:used(S,N),df_4(S)}, C!=11.
permutation_12(A,B,C,D,E,F,G,H,I,J,K,Lpf_3_0(2, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
pf_7_0(2, 1).
pf_7_0(19, 1).
pf_7_0(20, 1).
pf_7(X) :- pf_7_0(X).
universe(X) :- pf_7(X).
pf_9_0(2, 1).
pf_9_0(19, 1).
pf_9_0(20, 1).
pf_9(X) :- pf_9_0(X).
universe(X) :- pf_9(X).
pf_11_0(2, 1).
pf_11_0(19, 1).
pf_11_0(20, 1).
pf_11(X) :- pf_11_0(X).
universe(X) :- pf_11(X).
) :- universe(A), universe(B), universe(C), pf_3(D), universe(E), universe(F), universe(G), pf_7(H), universe(I), pf_9(J), universe(K), pf_11(L).
used_12(X,0) :- permutation_12(X, _, _, _, _, _, _, _, _, _, _, _). 
used_12(X,1) :- permutation_12(_, X, _, _, _, _, _, _, _, _, _, _). 
used_12(X,2) :- permutation_12(_, _, X, _, _, _, _, _, _, _, _, _). 
used_12(X,3) :- permutation_12(_, _, _, X, _, _, _, _, _, _, _, _). 
used_12(X,4) :- permutation_12(_, _, _, _, X, _, _, _, _, _, _, _). 
used_12(X,5) :- permutation_12(_, _, _, _, _, X, _, _, _, _, _, _). 
used_12(X,6) :- permutation_12(_, _, _, _, _, _, X, _, _, _, _, _). 
used_12(X,7) :- permutation_12(_, _, _, _, _, _, _, X, _, _, _, _). 
used_12(X,8) :- permutation_12(_, _, _, _, _, _, _, _, X, _, _, _). 
used_12(X,9) :- permutation_12(_, _, _, _, _, _, _, _, _, X, _, _). 
used_12(X,10) :- permutation_12(_, _, _, _, _, _, _, _, _, _, X, _). 
used_12(X,11) :- permutation_12(_, _, _, _, _, _, _, _, _, _, _, X). 
:- permutation_12(A,B,C,D,E,F,G,H,I,J,K,L), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
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
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1_0(2, 1).
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
:- C = #count{N:used(S,N),df_1(S)}, C!=11.
:- C = #count{N:used(S,N),df_1(S)}, C!=12.
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
:- C = #count{N:used(S,N),df_4(S)}, C!=6.
:- C = #count{N:used(S,N),df_4(S)}, C!=7.
:- C = #count{N:used(S,N),df_4(S)}, C!=8.
:- C = #count{N:used(S,N),df_4(S)}, C!=9.
:- C = #count{N:used(S,N),df_4(S)}, C!=10.
:- C = #count{N:used(S,N),df_4(S)}, C!=11.
:- C = #count{N:used(S,N),df_4(S)}, C!=12.
permutation_13(A,B,C,D,E,F,G,H,I,J,K,L,Mpf_3_0(2, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
pf_7_0(2, 1).
pf_7_0(19, 1).
pf_7_0(20, 1).
pf_7(X) :- pf_7_0(X).
universe(X) :- pf_7(X).
pf_9_0(2, 1).
pf_9_0(19, 1).
pf_9_0(20, 1).
pf_9(X) :- pf_9_0(X).
universe(X) :- pf_9(X).
pf_11_0(2, 1).
pf_11_0(19, 1).
pf_11_0(20, 1).
pf_11(X) :- pf_11_0(X).
universe(X) :- pf_11(X).
) :- universe(A), universe(B), universe(C), pf_3(D), universe(E), universe(F), universe(G), pf_7(H), universe(I), pf_9(J), universe(K), pf_11(L), universe(M).
used_13(X,0) :- permutation_13(X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_13(X,1) :- permutation_13(_, X, _, _, _, _, _, _, _, _, _, _, _). 
used_13(X,2) :- permutation_13(_, _, X, _, _, _, _, _, _, _, _, _, _). 
used_13(X,3) :- permutation_13(_, _, _, X, _, _, _, _, _, _, _, _, _). 
used_13(X,4) :- permutation_13(_, _, _, _, X, _, _, _, _, _, _, _, _). 
used_13(X,5) :- permutation_13(_, _, _, _, _, X, _, _, _, _, _, _, _). 
used_13(X,6) :- permutation_13(_, _, _, _, _, _, X, _, _, _, _, _, _). 
used_13(X,7) :- permutation_13(_, _, _, _, _, _, _, X, _, _, _, _, _). 
used_13(X,8) :- permutation_13(_, _, _, _, _, _, _, _, X, _, _, _, _). 
used_13(X,9) :- permutation_13(_, _, _, _, _, _, _, _, _, X, _, _, _). 
used_13(X,10) :- permutation_13(_, _, _, _, _, _, _, _, _, _, X, _, _). 
used_13(X,11) :- permutation_13(_, _, _, _, _, _, _, _, _, _, _, X, _). 
used_13(X,12) :- permutation_13(_, _, _, _, _, _, _, _, _, _, _, _, X). 
:- permutation_13(A,B,C,D,E,F,G,H,I,J,K,L,M), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
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
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
:- C = #count{N:used(S,N),df_0(S)}, C!=13.
df_1_0(2, 1).
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
:- C = #count{N:used(S,N),df_1(S)}, C!=11.
:- C = #count{N:used(S,N),df_1(S)}, C!=12.
:- C = #count{N:used(S,N),df_1(S)}, C!=13.
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
:- C = #count{N:used(S,N),df_4(S)}, C!=6.
:- C = #count{N:used(S,N),df_4(S)}, C!=7.
:- C = #count{N:used(S,N),df_4(S)}, C!=8.
:- C = #count{N:used(S,N),df_4(S)}, C!=9.
:- C = #count{N:used(S,N),df_4(S)}, C!=10.
:- C = #count{N:used(S,N),df_4(S)}, C!=11.
:- C = #count{N:used(S,N),df_4(S)}, C!=12.
:- C = #count{N:used(S,N),df_4(S)}, C!=13.
permutation_14(A,B,C,D,E,F,G,H,I,J,K,L,M,Npf_3_0(2, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
pf_7_0(2, 1).
pf_7_0(19, 1).
pf_7_0(20, 1).
pf_7(X) :- pf_7_0(X).
universe(X) :- pf_7(X).
pf_9_0(2, 1).
pf_9_0(19, 1).
pf_9_0(20, 1).
pf_9(X) :- pf_9_0(X).
universe(X) :- pf_9(X).
pf_11_0(2, 1).
pf_11_0(19, 1).
pf_11_0(20, 1).
pf_11(X) :- pf_11_0(X).
universe(X) :- pf_11(X).
) :- universe(A), universe(B), universe(C), pf_3(D), universe(E), universe(F), universe(G), pf_7(H), universe(I), pf_9(J), universe(K), pf_11(L), universe(M), universe(N).
used_14(X,0) :- permutation_14(X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,1) :- permutation_14(_, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,2) :- permutation_14(_, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,3) :- permutation_14(_, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_14(X,4) :- permutation_14(_, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_14(X,5) :- permutation_14(_, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_14(X,6) :- permutation_14(_, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_14(X,7) :- permutation_14(_, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_14(X,8) :- permutation_14(_, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_14(X,9) :- permutation_14(_, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_14(X,10) :- permutation_14(_, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_14(X,11) :- permutation_14(_, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_14(X,12) :- permutation_14(_, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_14(X,13) :- permutation_14(_, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- permutation_14(A,B,C,D,E,F,G,H,I,J,K,L,M,N), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
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
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
:- C = #count{N:used(S,N),df_0(S)}, C!=13.
:- C = #count{N:used(S,N),df_0(S)}, C!=14.
df_1_0(2, 1).
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
:- C = #count{N:used(S,N),df_1(S)}, C!=11.
:- C = #count{N:used(S,N),df_1(S)}, C!=12.
:- C = #count{N:used(S,N),df_1(S)}, C!=13.
:- C = #count{N:used(S,N),df_1(S)}, C!=14.
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
:- C = #count{N:used(S,N),df_4(S)}, C!=6.
:- C = #count{N:used(S,N),df_4(S)}, C!=7.
:- C = #count{N:used(S,N),df_4(S)}, C!=8.
:- C = #count{N:used(S,N),df_4(S)}, C!=9.
:- C = #count{N:used(S,N),df_4(S)}, C!=10.
:- C = #count{N:used(S,N),df_4(S)}, C!=11.
:- C = #count{N:used(S,N),df_4(S)}, C!=12.
:- C = #count{N:used(S,N),df_4(S)}, C!=14.
permutation_15(A,B,C,D,E,F,G,H,I,J,K,L,M,N,Opf_3_0(2, 1).
pf_3_0(19, 1).
pf_3_0(20, 1).
pf_3(X) :- pf_3_0(X).
universe(X) :- pf_3(X).
pf_7_0(2, 1).
pf_7_0(19, 1).
pf_7_0(20, 1).
pf_7(X) :- pf_7_0(X).
universe(X) :- pf_7(X).
pf_9_0(2, 1).
pf_9_0(19, 1).
pf_9_0(20, 1).
pf_9(X) :- pf_9_0(X).
universe(X) :- pf_9(X).
pf_11_0(2, 1).
pf_11_0(19, 1).
pf_11_0(20, 1).
pf_11(X) :- pf_11_0(X).
universe(X) :- pf_11(X).
) :- universe(A), universe(B), universe(C), pf_3(D), universe(E), universe(F), universe(G), pf_7(H), universe(I), pf_9(J), universe(K), pf_11(L), universe(M), universe(N), universe(O).
used_15(X,0) :- permutation_15(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,1) :- permutation_15(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,2) :- permutation_15(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,3) :- permutation_15(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,4) :- permutation_15(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_15(X,5) :- permutation_15(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_15(X,6) :- permutation_15(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_15(X,7) :- permutation_15(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_15(X,8) :- permutation_15(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_15(X,9) :- permutation_15(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_15(X,10) :- permutation_15(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_15(X,11) :- permutation_15(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_15(X,12) :- permutation_15(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_15(X,13) :- permutation_15(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_15(X,14) :- permutation_15(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- permutation_15(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(2, 1).
df_0_0(19, 1).
df_0_0(20, 1).
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
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
:- C = #count{N:used(S,N),df_0(S)}, C!=13.
:- C = #count{N:used(S,N),df_0(S)}, C!=14.
:- C = #count{N:used(S,N),df_0(S)}, C!=15.
df_1_0(2, 1).
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
:- C = #count{N:used(S,N),df_1(S)}, C!=11.
:- C = #count{N:used(S,N),df_1(S)}, C!=12.
:- C = #count{N:used(S,N),df_1(S)}, C!=13.
:- C = #count{N:used(S,N),df_1(S)}, C!=14.
:- C = #count{N:used(S,N),df_1(S)}, C!=15.
df_2_0(2, 1).
df_2_0(19, 1).
df_2_0(20, 1).
df_2(X) :- df_2_0(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=0.
:- C = #count{N:used(S,N),df_2(S)}, C!=1.
:- C = #count{N:used(S,N),df_2(S)}, C!=2.
:- C = #count{N:used(S,N),df_2(S)}, C!=3.
:- C = #count{N:used(S,N),df_2(S)}, C!=4.
:- C = #count{N:used(S,N),df_2(S)}, C!=5.
:- C = #count{N:used(S,N),df_2(S)}, C!=6.
:- C = #count{N:used(S,N),df_2(S)}, C!=7.
:- C = #count{N:used(S,N),df_2(S)}, C!=8.
:- C = #count{N:used(S,N),df_2(S)}, C!=9.
:- C = #count{N:used(S,N),df_2(S)}, C!=10.
:- C = #count{N:used(S,N),df_2(S)}, C!=11.
:- C = #count{N:used(S,N),df_2(S)}, C!=12.
df_3_0(2, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
universe(X) :- df_3(X).
:- C = #count{N:used(S,N),df_3(S)}, C!=15.
df_4_0(2, 1).
df_4_0(19, 1).
df_4_0(20, 1).
df_4(X) :- df_4_0(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
:- C = #count{N:used(S,N),df_4(S)}, C!=6.
:- C = #count{N:used(S,N),df_4(S)}, C!=7.
:- C = #count{N:used(S,N),df_4(S)}, C!=8.
:- C = #count{N:used(S,N),df_4(S)}, C!=9.
:- C = #count{N:used(S,N),df_4(S)}, C!=10.
:- C = #count{N:used(S,N),df_4(S)}, C!=11.
:- C = #count{N:used(S,N),df_4(S)}, C!=12.
:- C = #count{N:used(S,N),df_4(S)}, C!=14.
:- C = #count{N:used(S,N),df_4(S)}, C!=15.

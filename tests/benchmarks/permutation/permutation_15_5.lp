universe_0(1, 1).
universe_0(2, 1).
universe_0(12, 1).
universe_0(13, 1).
universe_0(14, 1).
universe_0(15, 1).
universe(X) :- universe_0(X).
universe_1(3,9).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(1, 1).
dom1_0(12, 1).
dom1_0(13, 1).
dom1_0(14, 1).
dom1_0(15, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
permutation_1(A) :- universe(A).
used_1(X,0) :- permutation_1(X). 
:- permutation_1(A), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
permutation_2(A,Bpf_1_0(1, 1).
pf_1_0(12, 1).
pf_1_0(13, 1).
pf_1_0(14, 1).
pf_1_0(15, 1).
pf_1(X) :- pf_1_0(X).
universe(X) :- pf_1(X).
) :- universe(A), pf_1(B).
used_2(X,0) :- permutation_2(X, _). 
used_2(X,1) :- permutation_2(_, X). 
:- permutation_2(A,B), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
permutation_3(A,B,Cpf_1_0(1, 1).
pf_1_0(12, 1).
pf_1_0(13, 1).
pf_1_0(14, 1).
pf_1_0(15, 1).
pf_1(X) :- pf_1_0(X).
universe(X) :- pf_1(X).
) :- universe(A), pf_1(B), universe(C).
used_3(X,0) :- permutation_3(X, _, _). 
used_3(X,1) :- permutation_3(_, X, _). 
used_3(X,2) :- permutation_3(_, _, X). 
:- permutation_3(A,B,C), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
permutation_4(A,B,C,Dpf_1_0(1, 1).
pf_1_0(12, 1).
pf_1_0(13, 1).
pf_1_0(14, 1).
pf_1_0(15, 1).
pf_1(X) :- pf_1_0(X).
universe(X) :- pf_1(X).
) :- universe(A), pf_1(B), universe(C), universe(D).
used_4(X,0) :- permutation_4(X, _, _, _). 
used_4(X,1) :- permutation_4(_, X, _, _). 
used_4(X,2) :- permutation_4(_, _, X, _). 
used_4(X,3) :- permutation_4(_, _, _, X). 
:- permutation_4(A,B,C,D), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.

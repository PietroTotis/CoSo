universe_0("e1",5).
universe(X) :- universe_0(X, _).
universe_1("e2",4).
universe(X) :- universe_1(X, _).
universe_2("e3",4).
universe(X) :- universe_2(X, _).
universe_3("e4", 1).
universe(X) :- universe_3(X, _).
universe_4("e5", 1).
universe(X) :- universe_4(X, _).
universe(X) :- universe(X).
dom1_0("e4",1).
dom1(X) :- dom1_0(X, _).
dom1_1("e3",4).
dom1(X) :- dom1_1(X, _).
universe(X) :- dom1(X).
df_0_0("e3",4).
df_0(X) :- df_0_0(X, _).
df_0_1("e4", 1).
df_0(X) :- df_0_1(X, _).
universe(X) :- df_0(X).
permutation_guess_2(A,B) :- pf_0(A), universe(B).
pf_0_0("e4",1).
pf_0(X) :- pf_0_0(X, _).
pf_0_1("e3",4).
pf_0(X) :- pf_0_1(X, _).
universe(X) :- pf_0(X).
1{permutation_2(A,B):permutation_guess_2(A,B)}1.
used_2(X,0) :- permutation_2(X, _). 
used_2(X,1) :- permutation_2(_, X). 
:- universe_0(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- universe_1(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- universe_2(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- universe_3(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- universe_4(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- dom1_0(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- dom1_1(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- C = #count{N:used_2(S,N),df_0(S)}, C=1.

#show permutation_2/2.
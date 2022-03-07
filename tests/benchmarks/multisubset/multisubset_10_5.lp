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
universe(X) :- universe_0(X).
universe(X) :- universe(X).
dom1_0(1, 1).
dom1_0(9, 1).
dom1_0(10, 1).
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
dom2(X) :- dom2_0(X).
universe(X) :- dom2(X).
multisubset_5(A,B,C,D,E) :- universe(A), universe(B), universe(C), universe(D), universe(E), A<=B, B<=C, C<=D, D<=E.
used_5(X,0) :- multisubset_5(X, _, _, _, _). 
used_5(X,1) :- multisubset_5(_, X, _, _, _). 
used_5(X,2) :- multisubset_5(_, _, X, _, _). 
used_5(X,3) :- multisubset_5(_, _, _, X, _). 
used_5(X,4) :- multisubset_5(_, _, _, _, X). 
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
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.

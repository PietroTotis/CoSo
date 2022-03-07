universe_0(1,1).
universe_0(2,8).
universe_0(10,4).
universe_0(14,2).
universe(X) :- universe_0(X, _).
universe(X) :- universe(X).
dom1_0(2,8).
dom1_0(10,4).
dom1_0(14,2).
dom1(X) :- dom1_0(X, _).
universe(X) :- dom1(X).
df_0_0(2,8).
df_0_0(10,4).
df_0_0(14,2).
df_0(X) :- df_0_0(X, _).
universe(X) :- df_0(X).
sequence_guess_5(A,B,C,D,E) :- universe(A), universe(B), universe(C), universe(D), universe(E).
1{sequence_5(A,B,C,D,E):sequence_guess_5(A,B,C,D,E)}1.
used_5(X,0) :- sequence_5(X, _, _, _, _). 
used_5(X,1) :- sequence_5(_, X, _, _, _). 
used_5(X,2) :- sequence_5(_, _, X, _, _). 
used_5(X,3) :- sequence_5(_, _, _, X, _). 
used_5(X,4) :- sequence_5(_, _, _, _, X). 
:- C = #count{N:used_5(S,N),df_0(S)}, C!=3.

#show sequence_5/5.
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
sequence_1(A) :- universe(A).
used_1(X,0) :- sequence_1(X). 
df_0_0(2,8).
df_0_0(10,4).
df_0_0(14,2).
df_0(X) :- df_0_0(X, _).
universe(X) :- df_0(X).
:- C = #count{N:used_1(S,N),df_0(S)}, C!=3.
sequence_2(A,B) :- universe(A), universe(B).
used_2(X,0) :- sequence_2(X, _). 
used_2(X,1) :- sequence_2(_, X). 
:- C = #count{N:used_2(S,N),df_0(S)}, C!=3.
sequence_3(A,B,C) :- universe(A), universe(B), universe(C).
used_3(X,0) :- sequence_3(X, _, _). 
used_3(X,1) :- sequence_3(_, X, _). 
used_3(X,2) :- sequence_3(_, _, X). 
:- C = #count{N:used_3(S,N),df_0(S)}, C!=3.
sequence_4(A,B,C,D) :- universe(A), universe(B), universe(C), universe(D).
used_4(X,0) :- sequence_4(X, _, _, _). 
used_4(X,1) :- sequence_4(_, X, _, _). 
used_4(X,2) :- sequence_4(_, _, X, _). 
used_4(X,3) :- sequence_4(_, _, _, X). 
:- C = #count{N:used_4(S,N),df_0(S)}, C!=3.
sequence_5(A,B,C,D,E) :- universe(A), universe(B), universe(C), universe(D), universe(E).
used_5(X,0) :- sequence_5(X, _, _, _, _). 
used_5(X,1) :- sequence_5(_, X, _, _, _). 
used_5(X,2) :- sequence_5(_, _, X, _, _). 
used_5(X,3) :- sequence_5(_, _, _, X, _). 
used_5(X,4) :- sequence_5(_, _, _, _, X). 
:- C = #count{N:used_5(S,N),df_0(S)}, C!=3.
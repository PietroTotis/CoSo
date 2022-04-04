universe_0("e1",8).
universe(X) :- universe_0(X, _).
universe_1("e2",4).
universe(X) :- universe_1(X, _).
universe_2("e3",2).
universe(X) :- universe_2(X, _).
universe_3("e4", 1).
universe(X) :- universe_3(X, _).
universe(X) :- universe(X).
dom1_0("e1",8).
dom1(X) :- dom1_0(X, _).
dom1_1("e2",4).
dom1(X) :- dom1_1(X, _).
dom1_2("e3",2).
dom1(X) :- dom1_2(X, _).
universe(X) :- dom1(X).
df_0_0("e1",8).
df_0(X) :- df_0_0(X, _).
df_0_1("e2",4).
df_0(X) :- df_0_1(X, _).
df_0_2("e3",2).
df_0(X) :- df_0_2(X, _).
universe(X) :- df_0(X).
sequence_guess_3(A,B,C) :- universe(A), universe(B), universe(C).
1{sequence_3(A,B,C):sequence_guess_3(A,B,C)}1.
used_3(X,0) :- sequence_3(X, _, _). 
used_3(X,1) :- sequence_3(_, X, _). 
used_3(X,2) :- sequence_3(_, _, X). 
:- C = #count{N:used_3(S,N),df_0(S)}, C=0.
:- C = #count{N:used_3(S,N),df_0(S)}, C=1.
:- C = #count{N:used_3(S,N),df_0(S)}, C=2.

#show sequence_3/3.
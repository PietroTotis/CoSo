universe_0("e1",3).
universe(X) :- universe_0(X, _).
universe_1("e2",4).
universe(X) :- universe_1(X, _).
universe_2("e3", 1).
universe(X) :- universe_2(X, _).
universe_3("e3", 1).
universe(X) :- universe_3(X, _).
universe_4("e3", 1).
universe(X) :- universe_4(X, _).
universe(X) :- universe(X).
dom1_0("e3",1).
dom1(X) :- dom1_0(X, _).
dom1_1("e2",4).
dom1(X) :- dom1_1(X, _).
universe(X) :- dom1(X).
subset_guess_6(A,B,C,D,E,F) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), A<=B, B<=C, C<=D, D<=E, E<=F.
1{subset_6(A,B,C,D,E,F):subset_guess_6(A,B,C,D,E,F)}1.
used_6(X,0) :- subset_6(X, _, _, _, _, _). 
used_6(X,1) :- subset_6(_, X, _, _, _, _). 
used_6(X,2) :- subset_6(_, _, X, _, _, _). 
used_6(X,3) :- subset_6(_, _, _, X, _, _). 
used_6(X,4) :- subset_6(_, _, _, _, X, _). 
used_6(X,5) :- subset_6(_, _, _, _, _, X). 
:- universe_0(S,SN), C = #count{N:used_6(S,N)}, C>SN.
:- universe_1(S,SN), C = #count{N:used_6(S,N)}, C>SN.
:- universe_2(S,SN), C = #count{N:used_6(S,N)}, C>SN.
:- universe_3(S,SN), C = #count{N:used_6(S,N)}, C>SN.
:- universe_4(S,SN), C = #count{N:used_6(S,N)}, C>SN.
:- dom1_0(S,SN), C = #count{N:used_6(S,N)}, C>SN.
:- dom1_1(S,SN), C = #count{N:used_6(S,N)}, C>SN.
:- C = #count{N:used_6(S,N),dom1(S)}, C=0.
:- C = #count{N:used_6(S,N),dom1(S)}, C=1.
:- C = #count{N:used_6(S,N),dom1(S)}, C=2.
:- C = #count{N:used_6(S,N),dom1(S)}, C=3.
:- C = #count{N:used_6(S,N),dom1(S)}, C=4.

#show subset_6/6.
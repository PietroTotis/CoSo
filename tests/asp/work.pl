universe_0(4,3).
universe(X) :- universe_0(X, _).
universe_1(7,4).
universe(X) :- universe_1(X, _).
universe_2(1, 1).
universe(X) :- universe_2(X, _).
universe_3(2, 1).
universe(X) :- universe_3(X, _).
universe_4(3, 1).
universe(X) :- universe_4(X, _).
universe(X) :- universe(X).
dom1_0(1,1).
dom1(X) :- dom1_0(X, _).
dom1_1(7,4).
dom1(X) :- dom1_1(X, _).
universe(X) :- dom1(X).
subset_guess_5(A,B,C,D,E) :- universe(A), universe(B), universe(C), universe(D), universe(E), A<=B, B<=C, C<=D, D<=E.
1{subset_5(A,B,C,D,E):subset_guess_5(A,B,C,D,E)}1.
used_5(X,0) :- subset_5(X, _, _, _, _). 
used_5(X,1) :- subset_5(_, X, _, _, _). 
used_5(X,2) :- subset_5(_, _, X, _, _). 
used_5(X,3) :- subset_5(_, _, _, X, _). 
used_5(X,4) :- subset_5(_, _, _, _, X). 
:- universe_0(S,SN), C = #count{N:used_5(S,N)}, C>SN.
:- universe_1(S,SN), C = #count{N:used_5(S,N)}, C>SN.
:- universe_2(S,SN), C = #count{N:used_5(S,N)}, C>SN.
:- universe_3(S,SN), C = #count{N:used_5(S,N)}, C>SN.
:- universe_4(S,SN), C = #count{N:used_5(S,N)}, C>SN.
:- dom1_0(S,SN), C = #count{N:used_5(S,N)}, C>SN.
:- dom1_1(S,SN), C = #count{N:used_5(S,N)}, C>SN.
count(S,C) :- universe_0(S,SN), C = #count{N:used_5(S,N)}.
count(S,C) :- universe_1(S,SN), C = #count{N:used_5(S,N)}.
count(S,C) :- universe_2(S,SN), C = #count{N:used_5(S,N)}.
count(S,C) :- universe_3(S,SN), C = #count{N:used_5(S,N)}.
count(S,C) :- universe_4(S,SN), C = #count{N:used_5(S,N)}.
count(S,C) :- dom1_0(S,SN), C = #count{N:used_5(S,N)}.
count(S,C) :- dom1_1(S,SN), C = #count{N:used_5(S,N)}.
:- C = #count{N:used_5(S,N),dom1(S)}, C=0.
:- C = #count{N:used_5(S,N),dom1(S)}, C=1.
:- C = #count{N:used_5(S,N),dom1(S)}, C=2.
:- C = #count{N:used_5(S,N),dom1(S)}, C=3.
:- C = #count{N:used_5(S,N),dom1(S)}, C=4.

#show subset_5/5.
#show count/2.
#show used_5/2.
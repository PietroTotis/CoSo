universe_0(3,5).
universe(X) :- universe_0(X, _).
universe_1(8,4).
universe(X) :- universe_1(X, _).
universe_2(12,4).
universe(X) :- universe_2(X, _).
universe_3(1, 1).
universe(X) :- universe_3(X, _).
universe_4(2, 1).
universe(X) :- universe_4(X, _).
universe(X) :- universe(X).
dom1_0(1,1).
dom1(X) :- dom1_0(X, _).
dom1_1(12,4).
dom1(X) :- dom1_1(X, _).
universe(X) :- dom1(X).
permutation_guess_2(A,B) :- universe(A), dom1(B).
1{permutation_2(A,B):permutation_guess_2(A,B)}1.
used_2(X,0) :- permutation_2(X, _). 
used_2(X,1) :- permutation_2(_, X). 
:- permutation_guess_2(A,B), universe_0(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- permutation_guess_2(A,B), universe_1(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- permutation_guess_2(A,B), universe_2(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- permutation_guess_2(A,B), universe_3(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- permutation_guess_2(A,B), universe_4(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- permutation_guess_2(A,B), dom1_0(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- permutation_guess_2(A,B), dom1_1(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- C = #count{N:used_2(S,N),dom1(S)}, C=1.

#show permutation_2/2.
bs_0("bs",2).
bs(X) :- bs_0(X, _).
universe(X) :- bs(X).
as_0("as",6).
as(X) :- as_0(X, _).
universe(X) :- as(X).
ns_0("ns",4).
ns(X) :- ns_0(X, _).
universe(X) :- ns(X).
permutation_guess_11(A,B,C,D,E,F,G,H,I,J,K) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K).
1{permutation_11(A,B,C,D,E,F,G,H,I,J,K):permutation_guess_11(A,B,C,D,E,F,G,H,I,J,K)}1.
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
:- bs_0(S,SN), C = #count{N:used_11(S,N)}, C>SN.
:- as_0(S,SN), C = #count{N:used_11(S,N)}, C>SN.
:- ns_0(S,SN), C = #count{N:used_11(S,N)}, C>SN.

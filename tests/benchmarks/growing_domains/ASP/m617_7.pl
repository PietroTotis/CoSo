bs_0("bs",2).
bs(X) :- bs_0(X, _).
universe(X) :- bs(X).
as_0("as",6).
as(X) :- as_0(X, _).
universe(X) :- as(X).
ns_0("ns",4).
ns(X) :- ns_0(X, _).
universe(X) :- ns(X).
permutation_guess_13(A,B,C,D,E,F,G,H,I,J,K,L,M) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M).
1{permutation_13(A,B,C,D,E,F,G,H,I,J,K,L,M):permutation_guess_13(A,B,C,D,E,F,G,H,I,J,K,L,M)}1.
used_13(X,0) :- permutation_13(X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_13(X,1) :- permutation_13(_, X, _, _, _, _, _, _, _, _, _, _, _). 
used_13(X,2) :- permutation_13(_, _, X, _, _, _, _, _, _, _, _, _, _). 
used_13(X,3) :- permutation_13(_, _, _, X, _, _, _, _, _, _, _, _, _). 
used_13(X,4) :- permutation_13(_, _, _, _, X, _, _, _, _, _, _, _, _). 
used_13(X,5) :- permutation_13(_, _, _, _, _, X, _, _, _, _, _, _, _). 
used_13(X,6) :- permutation_13(_, _, _, _, _, _, X, _, _, _, _, _, _). 
used_13(X,7) :- permutation_13(_, _, _, _, _, _, _, X, _, _, _, _, _). 
used_13(X,8) :- permutation_13(_, _, _, _, _, _, _, _, X, _, _, _, _). 
used_13(X,9) :- permutation_13(_, _, _, _, _, _, _, _, _, X, _, _, _). 
used_13(X,10) :- permutation_13(_, _, _, _, _, _, _, _, _, _, X, _, _). 
used_13(X,11) :- permutation_13(_, _, _, _, _, _, _, _, _, _, _, X, _). 
used_13(X,12) :- permutation_13(_, _, _, _, _, _, _, _, _, _, _, _, X). 
:- bs_0(S,SN), C = #count{N:used_13(S,N)}, C>SN.
:- as_0(S,SN), C = #count{N:used_13(S,N)}, C>SN.
:- ns_0(S,SN), C = #count{N:used_13(S,N)}, C>SN.

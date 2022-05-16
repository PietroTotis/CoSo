bs_0("bs",2).
bs(X) :- bs_0(X, _).
universe(X) :- bs(X).
as_0("as",6).
as(X) :- as_0(X, _).
universe(X) :- as(X).
ns_0("ns",4).
ns(X) :- ns_0(X, _).
universe(X) :- ns(X).
permutation_guess_12(A,B,C,D,E,F,G,H,I,J,K,L) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L).
1{permutation_12(A,B,C,D,E,F,G,H,I,J,K,L):permutation_guess_12(A,B,C,D,E,F,G,H,I,J,K,L)}1.
used_12(X,0) :- permutation_12(X, _, _, _, _, _, _, _, _, _, _, _). 
used_12(X,1) :- permutation_12(_, X, _, _, _, _, _, _, _, _, _, _). 
used_12(X,2) :- permutation_12(_, _, X, _, _, _, _, _, _, _, _, _). 
used_12(X,3) :- permutation_12(_, _, _, X, _, _, _, _, _, _, _, _). 
used_12(X,4) :- permutation_12(_, _, _, _, X, _, _, _, _, _, _, _). 
used_12(X,5) :- permutation_12(_, _, _, _, _, X, _, _, _, _, _, _). 
used_12(X,6) :- permutation_12(_, _, _, _, _, _, X, _, _, _, _, _). 
used_12(X,7) :- permutation_12(_, _, _, _, _, _, _, X, _, _, _, _). 
used_12(X,8) :- permutation_12(_, _, _, _, _, _, _, _, X, _, _, _). 
used_12(X,9) :- permutation_12(_, _, _, _, _, _, _, _, _, X, _, _). 
used_12(X,10) :- permutation_12(_, _, _, _, _, _, _, _, _, _, X, _). 
used_12(X,11) :- permutation_12(_, _, _, _, _, _, _, _, _, _, _, X). 
:- bs_0(S,SN), C = #count{N:used_12(S,N)}, C>SN.
:- as_0(S,SN), C = #count{N:used_12(S,N)}, C>SN.
:- ns_0(S,SN), C = #count{N:used_12(S,N)}, C>SN.

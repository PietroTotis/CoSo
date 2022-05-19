bs_0("bs",2).
bs(X) :- bs_0(X, _).
universe(X) :- bs(X).
as_0("as",6).
as(X) :- as_0(X, _).
universe(X) :- as(X).
ns_0("ns",4).
ns(X) :- ns_0(X, _).
universe(X) :- ns(X).
permutation_guess_10(A,B,C,D,E,F,G,H,I,J) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J).
1{permutation_10(A,B,C,D,E,F,G,H,I,J):permutation_guess_10(A,B,C,D,E,F,G,H,I,J)}1.
used_10(X,0) :- permutation_10(X, _, _, _, _, _, _, _, _, _). 
used_10(X,1) :- permutation_10(_, X, _, _, _, _, _, _, _, _). 
used_10(X,2) :- permutation_10(_, _, X, _, _, _, _, _, _, _). 
used_10(X,3) :- permutation_10(_, _, _, X, _, _, _, _, _, _). 
used_10(X,4) :- permutation_10(_, _, _, _, X, _, _, _, _, _). 
used_10(X,5) :- permutation_10(_, _, _, _, _, X, _, _, _, _). 
used_10(X,6) :- permutation_10(_, _, _, _, _, _, X, _, _, _). 
used_10(X,7) :- permutation_10(_, _, _, _, _, _, _, X, _, _). 
used_10(X,8) :- permutation_10(_, _, _, _, _, _, _, _, X, _). 
used_10(X,9) :- permutation_10(_, _, _, _, _, _, _, _, _, X). 
:- bs_0(S,SN), C = #count{N:used_10(S,N)}, C>SN.
:- as_0(S,SN), C = #count{N:used_10(S,N)}, C>SN.
:- ns_0(S,SN), C = #count{N:used_10(S,N)}, C>SN.

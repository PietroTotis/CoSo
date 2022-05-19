bs_0("bs",2).
bs(X) :- bs_0(X, _).
universe(X) :- bs(X).
as_0("as",6).
as(X) :- as_0(X, _).
universe(X) :- as(X).
ns_0("ns",4).
ns(X) :- ns_0(X, _).
universe(X) :- ns(X).
permutation_guess_8(A,B,C,D,E,F,G,H) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H).
1{permutation_8(A,B,C,D,E,F,G,H):permutation_guess_8(A,B,C,D,E,F,G,H)}1.
used_8(X,0) :- permutation_8(X, _, _, _, _, _, _, _). 
used_8(X,1) :- permutation_8(_, X, _, _, _, _, _, _). 
used_8(X,2) :- permutation_8(_, _, X, _, _, _, _, _). 
used_8(X,3) :- permutation_8(_, _, _, X, _, _, _, _). 
used_8(X,4) :- permutation_8(_, _, _, _, X, _, _, _). 
used_8(X,5) :- permutation_8(_, _, _, _, _, X, _, _). 
used_8(X,6) :- permutation_8(_, _, _, _, _, _, X, _). 
used_8(X,7) :- permutation_8(_, _, _, _, _, _, _, X). 
:- bs_0(S,SN), C = #count{N:used_8(S,N)}, C>SN.
:- as_0(S,SN), C = #count{N:used_8(S,N)}, C>SN.
:- ns_0(S,SN), C = #count{N:used_8(S,N)}, C>SN.

bs_0("bs",2).
bs(X) :- bs_0(X, _).
universe(X) :- bs(X).
as_0("as",6).
as(X) :- as_0(X, _).
universe(X) :- as(X).
ns_0("ns",4).
ns(X) :- ns_0(X, _).
universe(X) :- ns(X).
permutation_guess_5(A,B,C,D,E) :- universe(A), universe(B), universe(C), universe(D), universe(E).
1{permutation_5(A,B,C,D,E):permutation_guess_5(A,B,C,D,E)}1.
used_5(X,0) :- permutation_5(X, _, _, _, _). 
used_5(X,1) :- permutation_5(_, X, _, _, _). 
used_5(X,2) :- permutation_5(_, _, X, _, _). 
used_5(X,3) :- permutation_5(_, _, _, X, _). 
used_5(X,4) :- permutation_5(_, _, _, _, X). 
:- bs_0(S,SN), C = #count{N:used_5(S,N)}, C>SN.
:- as_0(S,SN), C = #count{N:used_5(S,N)}, C>SN.
:- ns_0(S,SN), C = #count{N:used_5(S,N)}, C>SN.

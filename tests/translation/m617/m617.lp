bs_0("b", 1).
bs(X) :- bs_0(X, _).
universe(X) :- bs(X).
as_0("a",3).
as(X) :- as_0(X, _).
universe(X) :- as(X).
ns_0("n",2).
ns(X) :- ns_0(X, _).
universe(X) :- ns(X).
permutation_guess_6(A,B,C,D,E,F) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F).
1{permutation_6(A,B,C,D,E,F):permutation_guess_6(A,B,C,D,E,F)}1.
used_6(X,0) :- permutation_6(X, _, _, _, _, _). 
used_6(X,1) :- permutation_6(_, X, _, _, _, _). 
used_6(X,2) :- permutation_6(_, _, X, _, _, _). 
used_6(X,3) :- permutation_6(_, _, _, X, _, _). 
used_6(X,4) :- permutation_6(_, _, _, _, X, _). 
used_6(X,5) :- permutation_6(_, _, _, _, _, X). 
:- bs_0(S,SN), C = #count{N:used_6(S,N)}, C>SN.
:- as_0(S,SN), C = #count{N:used_6(S,N)}, C>SN.
:- ns_0(S,SN), C = #count{N:used_6(S,N)}, C>SN.

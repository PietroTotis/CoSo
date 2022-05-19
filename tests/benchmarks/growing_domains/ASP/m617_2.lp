bs_0("bs",2).
bs(X) :- bs_0(X, _).
universe(X) :- bs(X).
as_0("as",6).
as(X) :- as_0(X, _).
universe(X) :- as(X).
ns_0("ns",4).
ns(X) :- ns_0(X, _).
universe(X) :- ns(X).
permutation_guess_4(A,B,C,D) :- universe(A), universe(B), universe(C), universe(D).
1{permutation_4(A,B,C,D):permutation_guess_4(A,B,C,D)}1.
used_4(X,0) :- permutation_4(X, _, _, _). 
used_4(X,1) :- permutation_4(_, X, _, _). 
used_4(X,2) :- permutation_4(_, _, X, _). 
used_4(X,3) :- permutation_4(_, _, _, X). 
:- bs_0(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- as_0(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- ns_0(S,SN), C = #count{N:used_4(S,N)}, C>SN.

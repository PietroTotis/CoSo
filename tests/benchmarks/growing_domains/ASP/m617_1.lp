bs_0("bs",2).
bs(X) :- bs_0(X, _).
universe(X) :- bs(X).
as_0("as",6).
as(X) :- as_0(X, _).
universe(X) :- as(X).
ns_0("ns",4).
ns(X) :- ns_0(X, _).
universe(X) :- ns(X).
permutation_guess_3(A,B,C) :- universe(A), universe(B), universe(C).
1{permutation_3(A,B,C):permutation_guess_3(A,B,C)}1.
used_3(X,0) :- permutation_3(X, _, _). 
used_3(X,1) :- permutation_3(_, X, _). 
used_3(X,2) :- permutation_3(_, _, X). 
:- bs_0(S,SN), C = #count{N:used_3(S,N)}, C>SN.
:- as_0(S,SN), C = #count{N:used_3(S,N)}, C>SN.
:- ns_0(S,SN), C = #count{N:used_3(S,N)}, C>SN.

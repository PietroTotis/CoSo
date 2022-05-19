bs_0("bs",2).
bs(X) :- bs_0(X, _).
universe(X) :- bs(X).
as_0("as",6).
as(X) :- as_0(X, _).
universe(X) :- as(X).
ns_0("ns",4).
ns(X) :- ns_0(X, _).
universe(X) :- ns(X).
permutation_guess_2(A,B) :- universe(A), universe(B).
1{permutation_2(A,B):permutation_guess_2(A,B)}1.
used_2(X,0) :- permutation_2(X, _). 
used_2(X,1) :- permutation_2(_, X). 
:- bs_0(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- as_0(S,SN), C = #count{N:used_2(S,N)}, C>SN.
:- ns_0(S,SN), C = #count{N:used_2(S,N)}, C>SN.

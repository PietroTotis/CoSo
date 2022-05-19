bs_0("bs",2).
bs(X) :- bs_0(X, _).
universe(X) :- bs(X).
as_0("as",6).
as(X) :- as_0(X, _).
universe(X) :- as(X).
ns_0("ns",4).
ns(X) :- ns_0(X, _).
universe(X) :- ns(X).
permutation_guess_9(A,B,C,D,E,F,G,H,I) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I).
1{permutation_9(A,B,C,D,E,F,G,H,I):permutation_guess_9(A,B,C,D,E,F,G,H,I)}1.
used_9(X,0) :- permutation_9(X, _, _, _, _, _, _, _, _). 
used_9(X,1) :- permutation_9(_, X, _, _, _, _, _, _, _). 
used_9(X,2) :- permutation_9(_, _, X, _, _, _, _, _, _). 
used_9(X,3) :- permutation_9(_, _, _, X, _, _, _, _, _). 
used_9(X,4) :- permutation_9(_, _, _, _, X, _, _, _, _). 
used_9(X,5) :- permutation_9(_, _, _, _, _, X, _, _, _). 
used_9(X,6) :- permutation_9(_, _, _, _, _, _, X, _, _). 
used_9(X,7) :- permutation_9(_, _, _, _, _, _, _, X, _). 
used_9(X,8) :- permutation_9(_, _, _, _, _, _, _, _, X). 
:- bs_0(S,SN), C = #count{N:used_9(S,N)}, C>SN.
:- as_0(S,SN), C = #count{N:used_9(S,N)}, C>SN.
:- ns_0(S,SN), C = #count{N:used_9(S,N)}, C>SN.

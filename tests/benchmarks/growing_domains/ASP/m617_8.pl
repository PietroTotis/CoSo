bs_0("bs",2).
bs(X) :- bs_0(X, _).
universe(X) :- bs(X).
as_0("as",6).
as(X) :- as_0(X, _).
universe(X) :- as(X).
ns_0("ns",4).
ns(X) :- ns_0(X, _).
universe(X) :- ns(X).
permutation_guess_14(A,B,C,D,E,F,G,H,I,J,K,L,M,N) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N).
1{permutation_14(A,B,C,D,E,F,G,H,I,J,K,L,M,N):permutation_guess_14(A,B,C,D,E,F,G,H,I,J,K,L,M,N)}1.
used_14(X,0) :- permutation_14(X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,1) :- permutation_14(_, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,2) :- permutation_14(_, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,3) :- permutation_14(_, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_14(X,4) :- permutation_14(_, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_14(X,5) :- permutation_14(_, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_14(X,6) :- permutation_14(_, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_14(X,7) :- permutation_14(_, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_14(X,8) :- permutation_14(_, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_14(X,9) :- permutation_14(_, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_14(X,10) :- permutation_14(_, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_14(X,11) :- permutation_14(_, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_14(X,12) :- permutation_14(_, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_14(X,13) :- permutation_14(_, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- bs_0(S,SN), C = #count{N:used_14(S,N)}, C>SN.
:- as_0(S,SN), C = #count{N:used_14(S,N)}, C>SN.
:- ns_0(S,SN), C = #count{N:used_14(S,N)}, C>SN.

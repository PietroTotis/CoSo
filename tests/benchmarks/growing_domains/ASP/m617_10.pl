bs_0("bs",2).
bs(X) :- bs_0(X, _).
universe(X) :- bs(X).
as_0("as",6).
as(X) :- as_0(X, _).
universe(X) :- as(X).
ns_0("ns",4).
ns(X) :- ns_0(X, _).
universe(X) :- ns(X).
permutation_guess_16(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P).
1{permutation_16(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P):permutation_guess_16(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P)}1.
used_16(X,0) :- permutation_16(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,1) :- permutation_16(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,2) :- permutation_16(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,3) :- permutation_16(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,4) :- permutation_16(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,5) :- permutation_16(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_16(X,6) :- permutation_16(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_16(X,7) :- permutation_16(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_16(X,8) :- permutation_16(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_16(X,9) :- permutation_16(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_16(X,10) :- permutation_16(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_16(X,11) :- permutation_16(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_16(X,12) :- permutation_16(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_16(X,13) :- permutation_16(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_16(X,14) :- permutation_16(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_16(X,15) :- permutation_16(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- bs_0(S,SN), C = #count{N:used_16(S,N)}, C>SN.
:- as_0(S,SN), C = #count{N:used_16(S,N)}, C>SN.
:- ns_0(S,SN), C = #count{N:used_16(S,N)}, C>SN.

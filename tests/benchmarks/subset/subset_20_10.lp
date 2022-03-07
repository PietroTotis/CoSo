universe_0(1, 1).
universe_0(2, 1).
universe_0(3, 1).
universe_0(4, 1).
universe_0(5, 1).
universe_0(6, 1).
universe_0(7, 1).
universe_0(8, 1).
universe_0(9, 1).
universe_0(10, 1).
universe_0(11, 1).
universe_0(12, 1).
universe_0(13, 1).
universe_0(14, 1).
universe_0(15, 1).
universe_0(16, 1).
universe_0(17, 1).
universe_0(18, 1).
universe_0(19, 1).
universe_0(20, 1).
universe(X) :- universe_0(X).
universe(X) :- universe(X).
dom1_0(1, 1).
dom1_0(2, 1).
dom1_0(18, 1).
dom1_0(19, 1).
dom1_0(20, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
dom2_0(1, 1).
dom2_0(13, 1).
dom2_0(14, 1).
dom2_0(15, 1).
dom2_0(16, 1).
dom2_0(17, 1).
dom2(X) :- dom2_0(X).
universe(X) :- dom2(X).
dom3_0(4, 1).
dom3_0(5, 1).
dom3_0(6, 1).
dom3_0(7, 1).
dom3_0(8, 1).
dom3_0(9, 1).
dom3_0(10, 1).
dom3_0(11, 1).
dom3_0(12, 1).
dom3_0(13, 1).
dom3_0(14, 1).
dom3_0(15, 1).
dom3_0(16, 1).
dom3_0(17, 1).
dom3(X) :- dom3_0(X).
universe(X) :- dom3(X).
subset_10(A,B,C,D,E,F,G,H,I,J) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), A<B, B<C, C<D, D<E, E<F, F<G, G<H, H<I, I<J.
used_10(X,0) :- subset_10(X, _, _, _, _, _, _, _, _, _). 
used_10(X,1) :- subset_10(_, X, _, _, _, _, _, _, _, _). 
used_10(X,2) :- subset_10(_, _, X, _, _, _, _, _, _, _). 
used_10(X,3) :- subset_10(_, _, _, X, _, _, _, _, _, _). 
used_10(X,4) :- subset_10(_, _, _, _, X, _, _, _, _, _). 
used_10(X,5) :- subset_10(_, _, _, _, _, X, _, _, _, _). 
used_10(X,6) :- subset_10(_, _, _, _, _, _, X, _, _, _). 
used_10(X,7) :- subset_10(_, _, _, _, _, _, _, X, _, _). 
used_10(X,8) :- subset_10(_, _, _, _, _, _, _, _, X, _). 
used_10(X,9) :- subset_10(_, _, _, _, _, _, _, _, _, X). 
:- subset_10(A,B,C,D,E,F,G,H,I,J), set(S,SN), C = #count{N:used(S,N)}, C>SN.

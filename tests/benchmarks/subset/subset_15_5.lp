universe_0(1, 1).
universe_0(9, 1).
universe_0(10, 1).
universe_0(11, 1).
universe_0(12, 1).
universe_0(13, 1).
universe_0(14, 1).
universe_0(15, 1).
universe(X) :- universe_0(X).
universe_1(2,7).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(9, 1).
dom1_0(10, 1).
dom1_0(11, 1).
dom1_0(12, 1).
dom1_0(13, 1).
dom1_0(14, 1).
dom1_0(15, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
subset_5(A,B,C,D,E) :- universe(A), universe(B), universe(C), universe(D), universe(E), A<B, B<C, C<D, D<E.
used_5(X,0) :- subset_5(X, _, _, _, _). 
used_5(X,1) :- subset_5(_, X, _, _, _). 
used_5(X,2) :- subset_5(_, _, X, _, _). 
used_5(X,3) :- subset_5(_, _, _, X, _). 
used_5(X,4) :- subset_5(_, _, _, _, X). 
:- subset_5(A,B,C,D,E), set(S,SN), C = #count{N:used(S,N)}, C>SN.
subset_6(A,B,C,D,E,F) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), A<B, B<C, C<D, D<E, E<F.
used_6(X,0) :- subset_6(X, _, _, _, _, _). 
used_6(X,1) :- subset_6(_, X, _, _, _, _). 
used_6(X,2) :- subset_6(_, _, X, _, _, _). 
used_6(X,3) :- subset_6(_, _, _, X, _, _). 
used_6(X,4) :- subset_6(_, _, _, _, X, _). 
used_6(X,5) :- subset_6(_, _, _, _, _, X). 
:- subset_6(A,B,C,D,E,F), set(S,SN), C = #count{N:used(S,N)}, C>SN.
subset_7(A,B,C,D,E,F,G) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), A<B, B<C, C<D, D<E, E<F, F<G.
used_7(X,0) :- subset_7(X, _, _, _, _, _, _). 
used_7(X,1) :- subset_7(_, X, _, _, _, _, _). 
used_7(X,2) :- subset_7(_, _, X, _, _, _, _). 
used_7(X,3) :- subset_7(_, _, _, X, _, _, _). 
used_7(X,4) :- subset_7(_, _, _, _, X, _, _). 
used_7(X,5) :- subset_7(_, _, _, _, _, X, _). 
used_7(X,6) :- subset_7(_, _, _, _, _, _, X). 
:- subset_7(A,B,C,D,E,F,G), set(S,SN), C = #count{N:used(S,N)}, C>SN.
subset_8(A,B,C,D,E,F,G,H) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), A<B, B<C, C<D, D<E, E<F, F<G, G<H.
used_8(X,0) :- subset_8(X, _, _, _, _, _, _, _). 
used_8(X,1) :- subset_8(_, X, _, _, _, _, _, _). 
used_8(X,2) :- subset_8(_, _, X, _, _, _, _, _). 
used_8(X,3) :- subset_8(_, _, _, X, _, _, _, _). 
used_8(X,4) :- subset_8(_, _, _, _, X, _, _, _). 
used_8(X,5) :- subset_8(_, _, _, _, _, X, _, _). 
used_8(X,6) :- subset_8(_, _, _, _, _, _, X, _). 
used_8(X,7) :- subset_8(_, _, _, _, _, _, _, X). 
:- subset_8(A,B,C,D,E,F,G,H), set(S,SN), C = #count{N:used(S,N)}, C>SN.
subset_9(A,B,C,D,E,F,G,H,I) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), A<B, B<C, C<D, D<E, E<F, F<G, G<H, H<I.
used_9(X,0) :- subset_9(X, _, _, _, _, _, _, _, _). 
used_9(X,1) :- subset_9(_, X, _, _, _, _, _, _, _). 
used_9(X,2) :- subset_9(_, _, X, _, _, _, _, _, _). 
used_9(X,3) :- subset_9(_, _, _, X, _, _, _, _, _). 
used_9(X,4) :- subset_9(_, _, _, _, X, _, _, _, _). 
used_9(X,5) :- subset_9(_, _, _, _, _, X, _, _, _). 
used_9(X,6) :- subset_9(_, _, _, _, _, _, X, _, _). 
used_9(X,7) :- subset_9(_, _, _, _, _, _, _, X, _). 
used_9(X,8) :- subset_9(_, _, _, _, _, _, _, _, X). 
:- subset_9(A,B,C,D,E,F,G,H,I), set(S,SN), C = #count{N:used(S,N)}, C>SN.
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
subset_11(A,B,C,D,E,F,G,H,I,J,K) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), A<B, B<C, C<D, D<E, E<F, F<G, G<H, H<I, I<J, J<K.
used_11(X,0) :- subset_11(X, _, _, _, _, _, _, _, _, _, _). 
used_11(X,1) :- subset_11(_, X, _, _, _, _, _, _, _, _, _). 
used_11(X,2) :- subset_11(_, _, X, _, _, _, _, _, _, _, _). 
used_11(X,3) :- subset_11(_, _, _, X, _, _, _, _, _, _, _). 
used_11(X,4) :- subset_11(_, _, _, _, X, _, _, _, _, _, _). 
used_11(X,5) :- subset_11(_, _, _, _, _, X, _, _, _, _, _). 
used_11(X,6) :- subset_11(_, _, _, _, _, _, X, _, _, _, _). 
used_11(X,7) :- subset_11(_, _, _, _, _, _, _, X, _, _, _). 
used_11(X,8) :- subset_11(_, _, _, _, _, _, _, _, X, _, _). 
used_11(X,9) :- subset_11(_, _, _, _, _, _, _, _, _, X, _). 
used_11(X,10) :- subset_11(_, _, _, _, _, _, _, _, _, _, X). 
:- subset_11(A,B,C,D,E,F,G,H,I,J,K), set(S,SN), C = #count{N:used(S,N)}, C>SN.
subset_12(A,B,C,D,E,F,G,H,I,J,K,L) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), A<B, B<C, C<D, D<E, E<F, F<G, G<H, H<I, I<J, J<K, K<L.
used_12(X,0) :- subset_12(X, _, _, _, _, _, _, _, _, _, _, _). 
used_12(X,1) :- subset_12(_, X, _, _, _, _, _, _, _, _, _, _). 
used_12(X,2) :- subset_12(_, _, X, _, _, _, _, _, _, _, _, _). 
used_12(X,3) :- subset_12(_, _, _, X, _, _, _, _, _, _, _, _). 
used_12(X,4) :- subset_12(_, _, _, _, X, _, _, _, _, _, _, _). 
used_12(X,5) :- subset_12(_, _, _, _, _, X, _, _, _, _, _, _). 
used_12(X,6) :- subset_12(_, _, _, _, _, _, X, _, _, _, _, _). 
used_12(X,7) :- subset_12(_, _, _, _, _, _, _, X, _, _, _, _). 
used_12(X,8) :- subset_12(_, _, _, _, _, _, _, _, X, _, _, _). 
used_12(X,9) :- subset_12(_, _, _, _, _, _, _, _, _, X, _, _). 
used_12(X,10) :- subset_12(_, _, _, _, _, _, _, _, _, _, X, _). 
used_12(X,11) :- subset_12(_, _, _, _, _, _, _, _, _, _, _, X). 
:- subset_12(A,B,C,D,E,F,G,H,I,J,K,L), set(S,SN), C = #count{N:used(S,N)}, C>SN.
subset_13(A,B,C,D,E,F,G,H,I,J,K,L,M) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), A<B, B<C, C<D, D<E, E<F, F<G, G<H, H<I, I<J, J<K, K<L, L<M.
used_13(X,0) :- subset_13(X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_13(X,1) :- subset_13(_, X, _, _, _, _, _, _, _, _, _, _, _). 
used_13(X,2) :- subset_13(_, _, X, _, _, _, _, _, _, _, _, _, _). 
used_13(X,3) :- subset_13(_, _, _, X, _, _, _, _, _, _, _, _, _). 
used_13(X,4) :- subset_13(_, _, _, _, X, _, _, _, _, _, _, _, _). 
used_13(X,5) :- subset_13(_, _, _, _, _, X, _, _, _, _, _, _, _). 
used_13(X,6) :- subset_13(_, _, _, _, _, _, X, _, _, _, _, _, _). 
used_13(X,7) :- subset_13(_, _, _, _, _, _, _, X, _, _, _, _, _). 
used_13(X,8) :- subset_13(_, _, _, _, _, _, _, _, X, _, _, _, _). 
used_13(X,9) :- subset_13(_, _, _, _, _, _, _, _, _, X, _, _, _). 
used_13(X,10) :- subset_13(_, _, _, _, _, _, _, _, _, _, X, _, _). 
used_13(X,11) :- subset_13(_, _, _, _, _, _, _, _, _, _, _, X, _). 
used_13(X,12) :- subset_13(_, _, _, _, _, _, _, _, _, _, _, _, X). 
:- subset_13(A,B,C,D,E,F,G,H,I,J,K,L,M), set(S,SN), C = #count{N:used(S,N)}, C>SN.
subset_14(A,B,C,D,E,F,G,H,I,J,K,L,M,N) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), A<B, B<C, C<D, D<E, E<F, F<G, G<H, H<I, I<J, J<K, K<L, L<M, M<N.
used_14(X,0) :- subset_14(X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,1) :- subset_14(_, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,2) :- subset_14(_, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,3) :- subset_14(_, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_14(X,4) :- subset_14(_, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_14(X,5) :- subset_14(_, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_14(X,6) :- subset_14(_, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_14(X,7) :- subset_14(_, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_14(X,8) :- subset_14(_, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_14(X,9) :- subset_14(_, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_14(X,10) :- subset_14(_, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_14(X,11) :- subset_14(_, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_14(X,12) :- subset_14(_, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_14(X,13) :- subset_14(_, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- subset_14(A,B,C,D,E,F,G,H,I,J,K,L,M,N), set(S,SN), C = #count{N:used(S,N)}, C>SN.
subset_15(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), A<B, B<C, C<D, D<E, E<F, F<G, G<H, H<I, I<J, J<K, K<L, L<M, M<N, N<O.
used_15(X,0) :- subset_15(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,1) :- subset_15(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,2) :- subset_15(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,3) :- subset_15(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,4) :- subset_15(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_15(X,5) :- subset_15(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_15(X,6) :- subset_15(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_15(X,7) :- subset_15(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_15(X,8) :- subset_15(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_15(X,9) :- subset_15(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_15(X,10) :- subset_15(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_15(X,11) :- subset_15(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_15(X,12) :- subset_15(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_15(X,13) :- subset_15(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_15(X,14) :- subset_15(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- subset_15(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O), set(S,SN), C = #count{N:used(S,N)}, C>SN.

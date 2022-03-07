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
dom1_0(15, 1).
dom1_0(16, 1).
dom1_0(17, 1).
dom1_0(18, 1).
dom1_0(19, 1).
dom1_0(20, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
dom2_0(9, 1).
dom2_0(10, 1).
dom2_0(11, 1).
dom2_0(12, 1).
dom2_0(13, 1).
dom2_0(14, 1).
dom2(X) :- dom2_0(X).
universe(X) :- dom2(X).
dom3_0(1, 1).
dom3_0(2, 1).
dom3_0(3, 1).
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
dom3(X) :- dom3_0(X).
universe(X) :- dom3(X).
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
subset_16(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), A<B, B<C, C<D, D<E, E<F, F<G, G<H, H<I, I<J, J<K, K<L, L<M, M<N, N<O, O<P.
used_16(X,0) :- subset_16(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,1) :- subset_16(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,2) :- subset_16(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,3) :- subset_16(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,4) :- subset_16(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,5) :- subset_16(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_16(X,6) :- subset_16(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_16(X,7) :- subset_16(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_16(X,8) :- subset_16(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_16(X,9) :- subset_16(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_16(X,10) :- subset_16(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_16(X,11) :- subset_16(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_16(X,12) :- subset_16(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_16(X,13) :- subset_16(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_16(X,14) :- subset_16(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_16(X,15) :- subset_16(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- subset_16(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P), set(S,SN), C = #count{N:used(S,N)}, C>SN.
subset_17(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q), A<B, B<C, C<D, D<E, E<F, F<G, G<H, H<I, I<J, J<K, K<L, L<M, M<N, N<O, O<P, P<Q.
used_17(X,0) :- subset_17(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,1) :- subset_17(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,2) :- subset_17(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,3) :- subset_17(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,4) :- subset_17(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,5) :- subset_17(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,6) :- subset_17(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_17(X,7) :- subset_17(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_17(X,8) :- subset_17(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_17(X,9) :- subset_17(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_17(X,10) :- subset_17(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_17(X,11) :- subset_17(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_17(X,12) :- subset_17(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_17(X,13) :- subset_17(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_17(X,14) :- subset_17(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_17(X,15) :- subset_17(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_17(X,16) :- subset_17(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- subset_17(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q), set(S,SN), C = #count{N:used(S,N)}, C>SN.
subset_18(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q), universe(R), A<B, B<C, C<D, D<E, E<F, F<G, G<H, H<I, I<J, J<K, K<L, L<M, M<N, N<O, O<P, P<Q, Q<R.
used_18(X,0) :- subset_18(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,1) :- subset_18(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,2) :- subset_18(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,3) :- subset_18(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,4) :- subset_18(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,5) :- subset_18(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,6) :- subset_18(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,7) :- subset_18(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_18(X,8) :- subset_18(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_18(X,9) :- subset_18(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_18(X,10) :- subset_18(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_18(X,11) :- subset_18(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_18(X,12) :- subset_18(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_18(X,13) :- subset_18(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_18(X,14) :- subset_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_18(X,15) :- subset_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_18(X,16) :- subset_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_18(X,17) :- subset_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- subset_18(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R), set(S,SN), C = #count{N:used(S,N)}, C>SN.
subset_19(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q), universe(R), universe(S), A<B, B<C, C<D, D<E, E<F, F<G, G<H, H<I, I<J, J<K, K<L, L<M, M<N, N<O, O<P, P<Q, Q<R, R<S.
used_19(X,0) :- subset_19(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,1) :- subset_19(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,2) :- subset_19(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,3) :- subset_19(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,4) :- subset_19(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,5) :- subset_19(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,6) :- subset_19(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,7) :- subset_19(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,8) :- subset_19(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_19(X,9) :- subset_19(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_19(X,10) :- subset_19(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_19(X,11) :- subset_19(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_19(X,12) :- subset_19(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_19(X,13) :- subset_19(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_19(X,14) :- subset_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_19(X,15) :- subset_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_19(X,16) :- subset_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_19(X,17) :- subset_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_19(X,18) :- subset_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- subset_19(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S), set(S,SN), C = #count{N:used(S,N)}, C>SN.
subset_20(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q), universe(R), universe(S), universe(T), A<B, B<C, C<D, D<E, E<F, F<G, G<H, H<I, I<J, J<K, K<L, L<M, M<N, N<O, O<P, P<Q, Q<R, R<S, S<T.
used_20(X,0) :- subset_20(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,1) :- subset_20(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,2) :- subset_20(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,3) :- subset_20(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,4) :- subset_20(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,5) :- subset_20(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,6) :- subset_20(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,7) :- subset_20(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,8) :- subset_20(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,9) :- subset_20(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_20(X,10) :- subset_20(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_20(X,11) :- subset_20(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_20(X,12) :- subset_20(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_20(X,13) :- subset_20(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_20(X,14) :- subset_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_20(X,15) :- subset_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_20(X,16) :- subset_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_20(X,17) :- subset_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_20(X,18) :- subset_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_20(X,19) :- subset_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- subset_20(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T), set(S,SN), C = #count{N:used(S,N)}, C>SN.

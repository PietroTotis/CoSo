universe_0(1, 1).
universe_0(2, 1).
universe_0(3, 1).
universe_0(7, 1).
universe_0(8, 1).
universe_0(9, 1).
universe_0(10, 1).
universe(X) :- universe_0(X).
universe_1(4,3).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(1, 1).
dom1_0(7, 1).
dom1_0(8, 1).
dom1_0(9, 1).
dom1_0(10, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
subset_5(A,B,C,D,E) :- universe(A), universe(B), universe(C), universe(D), universe(E), A<B, B<C, C<D, D<E.
used_5(X,0) :- subset_5(X, _, _, _, _). 
used_5(X,1) :- subset_5(_, X, _, _, _). 
used_5(X,2) :- subset_5(_, _, X, _, _). 
used_5(X,3) :- subset_5(_, _, _, X, _). 
used_5(X,4) :- subset_5(_, _, _, _, X). 
:- subset_5(A,B,C,D,E), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
subset_6(A,B,C,D,E,F) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), A<B, B<C, C<D, D<E, E<F.
used_6(X,0) :- subset_6(X, _, _, _, _, _). 
used_6(X,1) :- subset_6(_, X, _, _, _, _). 
used_6(X,2) :- subset_6(_, _, X, _, _, _). 
used_6(X,3) :- subset_6(_, _, _, X, _, _). 
used_6(X,4) :- subset_6(_, _, _, _, X, _). 
used_6(X,5) :- subset_6(_, _, _, _, _, X). 
:- subset_6(A,B,C,D,E,F), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
subset_7(A,B,C,D,E,F,G) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), A<B, B<C, C<D, D<E, E<F, F<G.
used_7(X,0) :- subset_7(X, _, _, _, _, _, _). 
used_7(X,1) :- subset_7(_, X, _, _, _, _, _). 
used_7(X,2) :- subset_7(_, _, X, _, _, _, _). 
used_7(X,3) :- subset_7(_, _, _, X, _, _, _). 
used_7(X,4) :- subset_7(_, _, _, _, X, _, _). 
used_7(X,5) :- subset_7(_, _, _, _, _, X, _). 
used_7(X,6) :- subset_7(_, _, _, _, _, _, X). 
:- subset_7(A,B,C,D,E,F,G), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(1, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
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
df_0_0(1, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
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
df_0_0(1, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
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
df_0_0(1, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.

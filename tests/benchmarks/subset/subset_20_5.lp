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
universe_0(17, 1).
universe_0(18, 1).
universe_0(19, 1).
universe_0(20, 1).
universe(X) :- universe_0(X).
universe_1(14,3).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(17, 1).
dom1_0(18, 1).
dom1_0(19, 1).
dom1_0(20, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
dom2_0(1, 1).
dom2_0(2, 1).
dom2_0(3, 1).
dom2_0(4, 1).
dom2_0(5, 1).
dom2_0(6, 1).
dom2_0(7, 1).
dom2_0(8, 1).
dom2_0(9, 1).
dom2_0(10, 1).
dom2_0(11, 1).
dom2_0(12, 1).
dom2_0(13, 1).
dom2(X) :- dom2_0(X).
universe(X) :- dom2(X).
subset_5(A,B,C,D,E) :- universe(A), universe(B), universe(C), universe(D), universe(E), A<B, B<C, C<D, D<E.
used_5(X,0) :- subset_5(X, _, _, _, _). 
used_5(X,1) :- subset_5(_, X, _, _, _). 
used_5(X,2) :- subset_5(_, _, X, _, _). 
used_5(X,3) :- subset_5(_, _, _, X, _). 
used_5(X,4) :- subset_5(_, _, _, _, X). 
:- subset_5(A,B,C,D,E), set(S,SN), C = #count{N:used(S,N)}, C>SN.
df_0_0(14,3).
df_0(X) :- df_0_0(X).
df_0_1(17, 1).
df_0_1(18, 1).
df_0_1(19, 1).
df_0_1(20, 1).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.

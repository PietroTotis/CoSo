tvs(1..3).
tvs(4..12).
universe(X) :- tvs(X).
defective(1..3).
universe(X) :- defective(X).
df_0(1..3).
universe(X) :- df_0(X).
subset_guess_5(A,B,C,D,E) :- universe(A), universe(B), universe(C), universe(D), universe(E), A<B, B<C, C<D, D<E.
1{subset_5(A,B,C,D,E):subset_guess_5(A,B,C,D,E)}1.
used_5(X,0) :- subset_5(X, _, _, _, _). 
used_5(X,1) :- subset_5(_, X, _, _, _). 
used_5(X,2) :- subset_5(_, _, X, _, _). 
used_5(X,3) :- subset_5(_, _, _, X, _). 
used_5(X,4) :- subset_5(_, _, _, _, X). 
:- C = #count{N:used_5(S,N),df_0(S)}, C=0.
:- C = #count{N:used_5(S,N),df_0(S)}, C=1.

tvs(1..16).
defective(1..6).
subset_guess_5(A,B,C,D,E) :- tvs(A), tvs(B), tvs(C), tvs(D), tvs(E), A<B, B<C, C<D, D<E.
1{subset_5(A,B,C,D,E):subset_guess_5(A,B,C,D,E)}1.
used_5(X,0) :- subset_5(X, _, _, _, _). 
used_5(X,1) :- subset_5(_, X, _, _, _). 
used_5(X,2) :- subset_5(_, _, X, _, _). 
used_5(X,3) :- subset_5(_, _, _, X, _). 
used_5(X,4) :- subset_5(_, _, _, _, X). 
:- C = #count{N:used_5(S,N),defective(S)}, C=0.
:- C = #count{N:used_5(S,N),defective(S)}, C=1.

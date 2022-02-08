set("e1", 1).
set("e2", 1).
set("e4", 1).

1{permutation(X,Y,Z):seq(X,Y,Z)}1.
seq(X,Y,Z) :- set(X,_), set(Y,_), set(Z,_).

used(X,1) :- permutation(X,_,_).
used(X,2) :- permutation(_,_,X).
used(X,3) :- permutation(_,X,_).

% test(X,Y,Z,S,C):-permutation(X,Y,Z), set(S,SN), C = #count{N:used(S,N)}.
:- permutation(X,Y,Z), set(S,SN), C = #count{N:used(S,N)}, C>SN.

#show permutation/3.
% #show test/5.
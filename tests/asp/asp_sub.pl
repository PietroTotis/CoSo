set(1, 1).
set(2, 3).
set(3, 2).

% one set per model
1{sub(X,Y,Z):seq(X,Y,Z)}1.

% domains plus breaking permutations
seq(X,Y,Z) :- set(X,_), set(Y,_), set(Z,_), X<Y, Y<Z .

% entity-position bind
used(X,1) :- sub(X,_,_).
used(X,2) :- sub(_,_,X).
used(X,3) :- sub(_,X,_).

% test(X,Y,Z,S,C):-permutation(X,Y,Z), set(S,SN), C = #count{N:used(S,N)}.

% do not use more copies than available
:- sub(X,Y,Z), set(S,SN), C = #count{N:used(S,N)}, C>SN.

#show sub/3.
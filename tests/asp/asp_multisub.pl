set(1, 1).
set(2, 3).
set(3, 2).

% one multisubset per model
1{multisub(X,Y,Z):seq(X,Y,Z)}1.

% domains plus breaking permutations
seq(X,Y,Z) :- set(X,_), set(Y,_), set(Z,_), X<=Y, Y<=Z .

% entity-position bind
used(X,1) :- multisub(X,_,_).
used(X,2) :- multisub(_,_,X).
used(X,3) :- multisub(_,X,_).

% test(X,Y,Z,S,C):-permutation(X,Y,Z), set(S,SN), C = #count{N:used(S,N)}.
% :- multisub(X,Y,Z), set(S,SN), C = #count{N:used(S,N)}, C>SN.

#show multisub/3.
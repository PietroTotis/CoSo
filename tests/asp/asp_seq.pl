set("e1", 1).
set("e2", 3).
set("e3", 2).

special("e1").
special("e3").

% one sequence per model
1{sequence(X,Y,Z):seq(X,Y,Z)}1.

% domains
seq(X,Y,Z) :- set(X,_), set(Y,_), set(Z,_).

% domains with positional constraints

% seq(X,"e2",Z) :- set(X,_), set(Z,_).

% counting constraints

used(X,1) :- sequence(X,_,_).
used(X,2) :- sequence(_,_,X).
used(X,3) :- sequence(_,X,_).

:- C = #count{1:used(S,_),special(S)}, C!=1.
:- sequence(X,Y,Z), set(S,SN), C = #count{N:used(S,N)}, C>SN.

#show sequence/3.
% #show test/5.
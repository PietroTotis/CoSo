set("e1", 1).
set("e2", 3).
set("e3", 2).

1{sequence(X,Y,Z):seq(X,Y,Z)}1.
seq(X,Y,Z) :- set(X,_), set(Y,_), set(Z,_).

#show sequence/3.
% #show test/5.
set("e1", 1).
set("e2", 2).
% set("e3", 1).

int(0..3).

part("p1").
part("p2").
% part("p3").

1{put(E,N,P): int(N), N<=EN} 1 :- set(E, EN), part(P).

:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN. 
:- part(P), #count{E,N:put(E,N,P), N>0}==0.
% sum_p(P,S) :- part(P), S=#sum{N,E:put(E,N,P)}.
% :- sum_p(P,0), part(P).


% part(P, A, B, C) :- put("e1",A,P),put("e2",B,P),put("e3",C,P).
part(P, A, B) :- put("e1",A,P),put("e2",B,P).
% empty("p3").

% 1{partition(X,Y,Z):seq(X,Y,Z)}1.
% seq(X,Y,Z) :- set(X,_), set(Y,_), set(Z,_).

% used(X,1) :- partition(X,_,_).
% used(X,2) :- partition(_,_,X).
% used(X,3) :- partition(_,X,_).

% test(X,Y,Z,S,C):-permutation(X,Y,Z), set(S,SN), C = #count{N:used(S,N)}.
% :- partition(X,Y,Z), set(S,SN), C = #count{N:used(S,N)}, C>SN.

% #show put/3.
#show part/3.
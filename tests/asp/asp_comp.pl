set("e1", 1).
set("e2", 2).
% set("e3", 1).

int(0..3).

part("p1").
part("p2").
% part("p3").

% for each part choose how many entities from set E to put there (less than or equal available)
1{put(E,N,P): int(N), N<=EN} 1 :- set(E, EN), part(P).

% the sum of entities E across parts is the number of available E
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN. 
% all parts are non-empty
:- part(P), #count{E,N:put(E,N,P), N>0}==0.

% output for each part how many entities from each 
part(P, A, B) :- put("e1",A,P),put("e2",B,P).

% #show put/3.
#show part/3.
set("e1", 4).
set("e2", 2).
set("e3", 1).

df1("e2").
df1("e3").

int(0..4).

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
part(P, A, B, C) :- put("e1",A,P),put("e2",B,P),put("e3",C,P).

% size constraints 
:-  C=#sum{N,E:put(E,N,"p2")}, C>2.

% count constraints
% the number of parts with 2 elements from df1 cannot be 1 

count_2(P,S) :-  S=#sum{N,E:put(E,N,P), df1(E)}, part(P),S=2.

count(C) :- C=#count{P:count_2(P,2)}.

:- count(1).


% #show put/3.
#show part/4.
% #show count/1.
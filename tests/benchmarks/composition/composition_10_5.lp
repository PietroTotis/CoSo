universe_0(1, 1).
universe_0(2, 1).
universe_0(9, 1).
universe_0(10, 1).
universe(X) :- universe_0(X).
universe_1(3,6).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(1, 1).
dom1_0(9, 1).
dom1_0(10, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
part(0).
part(1).
part(2).
part(3).
part(4).
1{put(E,N,P): int(N), N<=EN} 1 :- universe_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_1(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom1_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.
part(0).
part(1).
part(2).
part(3).
part(4).
part(5).
1{put(E,N,P): int(N), N<=EN} 1 :- universe_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_1(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom1_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.
part(0).
part(1).
part(2).
part(3).
part(4).
part(5).
part(6).
1{put(E,N,P): int(N), N<=EN} 1 :- universe_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_1(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom1_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.
part(0).
part(1).
part(2).
part(3).
part(4).
part(5).
part(6).
part(7).
1{put(E,N,P): int(N), N<=EN} 1 :- universe_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_1(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom1_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.
part(0).
part(1).
part(2).
part(3).
part(4).
part(5).
part(6).
part(7).
part(8).
1{put(E,N,P): int(N), N<=EN} 1 :- universe_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_1(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom1_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.
part(0).
part(1).
part(2).
part(3).
part(4).
part(5).
part(6).
part(7).
part(8).
part(9).
1{put(E,N,P): int(N), N<=EN} 1 :- universe_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_1(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom1_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.

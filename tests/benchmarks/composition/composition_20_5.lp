universe_0(1, 1).
universe_0(2, 1).
universe_0(3, 1).
universe_0(4, 1).
universe_0(5, 1).
universe_0(6, 1).
universe_0(7, 1).
universe_0(8, 1).
universe_0(9, 1).
universe_0(10, 1).
universe_0(11, 1).
universe_0(12, 1).
universe_0(13, 1).
universe_0(14, 1).
universe_0(15, 1).
universe_0(16, 1).
universe_0(17, 1).
universe_0(18, 1).
universe_0(19, 1).
universe_0(20, 1).
universe(X) :- universe_0(X).
universe(X) :- universe(X).
dom1_0(3, 1).
dom1_0(4, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
dom2_0(1, 1).
dom2_0(2, 1).
dom2_0(15, 1).
dom2_0(16, 1).
dom2_0(17, 1).
dom2_0(18, 1).
dom2_0(19, 1).
dom2_0(20, 1).
dom2(X) :- dom2_0(X).
universe(X) :- dom2(X).
dom3_0(5, 1).
dom3_0(6, 1).
dom3_0(7, 1).
dom3_0(8, 1).
dom3_0(9, 1).
dom3_0(10, 1).
dom3_0(11, 1).
dom3_0(12, 1).
dom3_0(13, 1).
dom3_0(14, 1).
dom3(X) :- dom3_0(X).
universe(X) :- dom3(X).
part(0).
part(1).
part(2).
part(3).
part(4).
1{put(E,N,P): int(N), N<=EN} 1 :- universe_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom1_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom2_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom3_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.

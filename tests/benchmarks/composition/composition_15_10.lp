universe_0(1, 1).
universe_0(2, 1).
universe_0(10, 1).
universe_0(11, 1).
universe_0(12, 1).
universe_0(13, 1).
universe_0(14, 1).
universe_0(15, 1).
universe(X) :- universe_0(X).
universe_1(3,7).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(1, 1).
dom1_0(2, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
dom2_0(2, 1).
dom2_0(14, 1).
dom2_0(15, 1).
dom2(X) :- dom2_0(X).
universe(X) :- dom2(X).
dom3_0(1, 1).
dom3_0(10, 1).
dom3_0(11, 1).
dom3_0(12, 1).
dom3_0(13, 1).
dom3(X) :- dom3_0(X).
universe(X) :- dom3(X).
dom4_0(1, 1).
dom4_0(10, 1).
dom4_0(11, 1).
dom4_0(12, 1).
dom4_0(13, 1).
dom4(X) :- dom4_0(X).
universe(X) :- dom4(X).
part(0).
1{put(E,N,P): int(N), N<=EN} 1 :- universe_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_1(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom1_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom2_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom3_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom4_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.
part(0).
part(1).
1{put(E,N,P): int(N), N<=EN} 1 :- universe_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_1(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom1_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom2_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom3_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom4_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.
part(0).
part(1).
part(2).
1{put(E,N,P): int(N), N<=EN} 1 :- universe_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_1(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom1_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom2_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom3_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom4_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.
part(0).
part(1).
part(2).
part(3).
1{put(E,N,P): int(N), N<=EN} 1 :- universe_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_1(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom1_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom2_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom3_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom4_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.
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
1{put(E,N,P): int(N), N<=EN} 1 :- dom2_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom3_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom4_0(E, EN), part(P).
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
1{put(E,N,P): int(N), N<=EN} 1 :- dom2_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom3_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom4_0(E, EN), part(P).
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
1{put(E,N,P): int(N), N<=EN} 1 :- dom2_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom3_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom4_0(E, EN), part(P).
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
1{put(E,N,P): int(N), N<=EN} 1 :- dom2_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom3_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom4_0(E, EN), part(P).
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
1{put(E,N,P): int(N), N<=EN} 1 :- dom2_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom3_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom4_0(E, EN), part(P).
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
1{put(E,N,P): int(N), N<=EN} 1 :- dom2_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom3_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom4_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.

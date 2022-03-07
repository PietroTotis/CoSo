universe_0(1,20).
universe(X) :- universe_0(X).
universe_1(6, 1).
universe_1(7, 1).
universe_1(8, 1).
universe_1(9, 1).
universe_1(10, 1).
universe_1(11, 1).
universe_1(12, 1).
universe_1(13, 1).
universe_1(14, 1).
universe_1(15, 1).
universe_1(16, 1).
universe_1(17, 1).
universe_1(18, 1).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(6, 1).
dom1_0(7, 1).
dom1_0(8, 1).
dom1_0(9, 1).
dom1_0(10, 1).
dom1_0(11, 1).
dom1_0(12, 1).
dom1_0(13, 1).
dom1_0(14, 1).
dom1_0(15, 1).
dom1_0(16, 1).
dom1_0(17, 1).
dom1_0(18, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
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
part(10).
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
part(10).
part(11).
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
part(10).
part(11).
part(12).
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
part(10).
part(11).
part(12).
part(13).
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
part(10).
part(11).
part(12).
part(13).
part(14).
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
part(10).
part(11).
part(12).
part(13).
part(14).
part(15).
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
part(10).
part(11).
part(12).
part(13).
part(14).
part(15).
part(16).
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
part(10).
part(11).
part(12).
part(13).
part(14).
part(15).
part(16).
part(17).
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
part(10).
part(11).
part(12).
part(13).
part(14).
part(15).
part(16).
part(17).
part(18).
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
part(10).
part(11).
part(12).
part(13).
part(14).
part(15).
part(16).
part(17).
part(18).
part(19).
1{put(E,N,P): int(N), N<=EN} 1 :- universe_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_1(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom1_0(E, EN), part(P).
:- set(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.

universe_0("e2",2).
universe(X) :- universe_0(X, _).
universe_1("e3",2).
universe(X) :- universe_1(X, _).
universe_2("e4", 1).
universe(X) :- universe_2(X, _).
universe(X) :- universe(X).
dom1_0("e4",1).
dom1(X) :- dom1_0(X, _).
dom1_1("e3",2).
dom1(X) :- dom1_1(X, _).
universe(X) :- dom1(X).
int(0..5).
df_0_0("e4",1).
df_0(X) :- df_0_0(X, _).
df_0_1("e3",2).
df_0(X) :- df_0_1(X, _).
universe(X) :- df_0(X).
part(0).
part(1).
part(2).
1{put(E,N,P): int(N), N<=EN} 1 :- universe_0(E, EN), part(P).
:- universe_0(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_1(E, EN), part(P).
:- universe_1(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_2(E, EN), part(P).
:- universe_2(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom1_0(E, EN), part(P).
:- dom1_0(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- dom1_1(E, EN), part(P).
:- dom1_1(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.
cf_0_3(P,S) :-  S=#sum{N,E:put(E,N,P), df_0(E)}, part(P), S=3.
cf_0_4(P,S) :-  S=#sum{N,E:put(E,N,P), df_0(E)}, part(P), S=4.
count_0(C) :- C1=#count{P:cf_0_4(P,4)}, C2=#count{P:cf_0_3(P,3)}, C=C1+C2.
:- count_0(0).

part(P, A, B, C) :- put("e2",A,P),put("e3",B,P),put("e4",C,P).
#show part/4.
#show count_0/1. 
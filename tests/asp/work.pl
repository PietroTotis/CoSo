universe_0("e1",4).
universe(X) :- universe_0(X, _).
universe_1("e2",2).
universe(X) :- universe_1(X, _).
universe_2("e3", 1).
universe(X) :- universe_2(X, _).
universe(X) :- universe(X).
df_0("e3",1).
df(X) :- df_0(X, _).
df_1("e2",2).
df(X) :- df_1(X, _).
universe(X) :- df(X).
int(0..7).
df_0_0("e2",2).
df_0(X) :- df_0_0(X, _).
df_0_1("e3", 1).
df_0(X) :- df_0_1(X, _).
universe(X) :- df_0(X).
part(0).
part(1).
1{put(E,N,P): int(N), N<=EN} 1 :- universe_0(E, EN), part(P).
:- universe_0(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_1(E, EN), part(P).
:- universe_1(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- universe_2(E, EN), part(P).
:- universe_2(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- df_0(E, EN), part(P).
:- df_0(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
1{put(E,N,P): int(N), N<=EN} 1 :- df_1(E, EN), part(P).
:- df_1(E,EN), #sum{N,P:put(E,N,P),part(P)}!=EN.
:- part(P), #count{E,N:put(E,N,P), N>0}==0.
:-  C=#sum{N,E:put(E,N,1)}, C=0.
:-  C=#sum{N,E:put(E,N,1)}, C=1.
:-  C=#sum{N,E:put(E,N,1)}, C=2.


cf_0_2(P,S) :-  S=#sum{N,E:put(E,N,P), df_0(E)}, part(P), S=2.
count_0(C) :- C2=#count{P:cf_0_2(P,2)}, C=C2.
:- count_0(0).
:- count_0(2).


show_part(P,A,B,C) :- put("e1",A,P), put("e2",B,P), put("e3",C,P).

#show show_part/4.
#show cf_0_2/2.
workers(1..15).
part(0).
part(1).
part(2).
1{put(E,P):part(P)} 1 :- workers(E).
:-  C=#count{E:put(E,0)}, C!=8.
:-  C=#count{E:put(E,1)}, C!=5.
:-  C=#count{E:put(E,2)}, C!=2.

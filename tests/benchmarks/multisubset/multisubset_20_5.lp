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
universe(X) :- universe_0(X).
universe_1(17,4).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(4, 1).
dom1_0(5, 1).
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
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
multisubset_1(A) :- universe(A), .
used_1(X,0) :- multisubset_1(X). 
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
multisubset_2(A,B) :- universe(A), universe(B), A<=B.
used_2(X,0) :- multisubset_2(X, _). 
used_2(X,1) :- multisubset_2(_, X). 
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
multisubset_3(A,B,C) :- universe(A), universe(B), universe(C), A<=B, B<=C.
used_3(X,0) :- multisubset_3(X, _, _). 
used_3(X,1) :- multisubset_3(_, X, _). 
used_3(X,2) :- multisubset_3(_, _, X). 
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
multisubset_4(A,B,C,D) :- universe(A), universe(B), universe(C), universe(D), A<=B, B<=C, C<=D.
used_4(X,0) :- multisubset_4(X, _, _, _). 
used_4(X,1) :- multisubset_4(_, X, _, _). 
used_4(X,2) :- multisubset_4(_, _, X, _). 
used_4(X,3) :- multisubset_4(_, _, _, X). 
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
multisubset_5(A,B,C,D,E) :- universe(A), universe(B), universe(C), universe(D), universe(E), A<=B, B<=C, C<=D, D<=E.
used_5(X,0) :- multisubset_5(X, _, _, _, _). 
used_5(X,1) :- multisubset_5(_, X, _, _, _). 
used_5(X,2) :- multisubset_5(_, _, X, _, _). 
used_5(X,3) :- multisubset_5(_, _, _, X, _). 
used_5(X,4) :- multisubset_5(_, _, _, _, X). 
df_0_0(4, 1).
df_0_0(5, 1).
df_0_0(6, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(10, 1).
df_0_0(11, 1).
df_0_0(12, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0(X) :- df_0_0(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.

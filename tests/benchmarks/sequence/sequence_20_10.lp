universe_0(1, 1).
universe_0(2, 1).
universe_0(13, 1).
universe_0(14, 1).
universe_0(15, 1).
universe_0(16, 1).
universe_0(17, 1).
universe_0(18, 1).
universe_0(19, 1).
universe_0(20, 1).
universe(X) :- universe_0(X).
universe_1(3,10).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(2, 1).
dom1_0(19, 1).
dom1_0(20, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
dom2_0(13, 1).
dom2_0(14, 1).
dom2_0(15, 1).
dom2_0(16, 1).
dom2_0(17, 1).
dom2_0(18, 1).
dom2(X) :- dom2_0(X).
universe(X) :- dom2(X).
dom3_0(16, 1).
dom3_0(17, 1).
dom3_0(18, 1).
dom3_0(19, 1).
dom3_0(20, 1).
dom3(X) :- dom3_0(X).
universe(X) :- dom3(X).
sequence_10(A,B,C,D,E,F,G,H,I,Jpf_7_0(13, 1).
pf_7_0(14, 1).
pf_7_0(15, 1).
pf_7_0(16, 1).
pf_7_0(17, 1).
pf_7_0(18, 1).
pf_7(X) :- pf_7_0(X).
universe(X) :- pf_7(X).
pf_8_0(13, 1).
pf_8_0(14, 1).
pf_8_0(15, 1).
pf_8_0(16, 1).
pf_8_0(17, 1).
pf_8_0(18, 1).
pf_8(X) :- pf_8_0(X).
universe(X) :- pf_8(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), pf_7(H), pf_8(I), universe(J).
used_10(X,0) :- sequence_10(X, _, _, _, _, _, _, _, _, _). 
used_10(X,1) :- sequence_10(_, X, _, _, _, _, _, _, _, _). 
used_10(X,2) :- sequence_10(_, _, X, _, _, _, _, _, _, _). 
used_10(X,3) :- sequence_10(_, _, _, X, _, _, _, _, _, _). 
used_10(X,4) :- sequence_10(_, _, _, _, X, _, _, _, _, _). 
used_10(X,5) :- sequence_10(_, _, _, _, _, X, _, _, _, _). 
used_10(X,6) :- sequence_10(_, _, _, _, _, _, X, _, _, _). 
used_10(X,7) :- sequence_10(_, _, _, _, _, _, _, X, _, _). 
used_10(X,8) :- sequence_10(_, _, _, _, _, _, _, _, X, _). 
used_10(X,9) :- sequence_10(_, _, _, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(13, 1).
df_0_0(14, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0(X) :- df_0_0(X).
df_0_1(3,10).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.

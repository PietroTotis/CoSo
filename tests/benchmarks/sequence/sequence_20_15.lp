universe_0(1, 1).
universe_0(2, 1).
universe_0(3, 1).
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
universe_1(4,3).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(3, 1).
dom1_0(17, 1).
dom1_0(18, 1).
dom1_0(19, 1).
dom1_0(20, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
dom2_0(1, 1).
dom2_0(2, 1).
dom2_0(3, 1).
dom2_0(7, 1).
dom2_0(8, 1).
dom2_0(9, 1).
dom2_0(10, 1).
dom2_0(11, 1).
dom2_0(12, 1).
dom2_0(13, 1).
dom2_0(14, 1).
dom2_0(15, 1).
dom2_0(16, 1).
dom2_0(17, 1).
dom2_0(18, 1).
dom2_0(19, 1).
dom2_0(20, 1).
dom2(X) :- dom2_0(X).
universe(X) :- dom2(X).
dom3_0(10, 1).
dom3_0(11, 1).
dom3_0(12, 1).
dom3_0(13, 1).
dom3_0(14, 1).
dom3_0(15, 1).
dom3_0(16, 1).
dom3(X) :- dom3_0(X).
universe(X) :- dom3(X).
dom4_0(3, 1).
dom4_0(19, 1).
dom4_0(20, 1).
dom4(X) :- dom4_0(X).
universe(X) :- dom4(X).
dom5_0(2, 1).
dom5_0(15, 1).
dom5_0(16, 1).
dom5(X) :- dom5_0(X).
universe(X) :- dom5(X).
sequence_1(A) :- universe(A).
used_1(X,0) :- sequence_1(X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
sequence_2(A,Bpf_1(X) :- pf_1_0(X).
pf_1_1(4,3).
pf_1(X) :- pf_1_1(X).
universe(X) :- pf_1(X).
) :- universe(A), pf_1(B).
used_2(X,0) :- sequence_2(X, _). 
used_2(X,1) :- sequence_2(_, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
sequence_3(A,B,Cpf_1(X) :- pf_1_0(X).
pf_1_1(4,3).
pf_1(X) :- pf_1_1(X).
universe(X) :- pf_1(X).
) :- universe(A), pf_1(B), universe(C).
used_3(X,0) :- sequence_3(X, _, _). 
used_3(X,1) :- sequence_3(_, X, _). 
used_3(X,2) :- sequence_3(_, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
sequence_4(A,B,C,Dpf_1(X) :- pf_1_0(X).
pf_1_1(4,3).
pf_1(X) :- pf_1_1(X).
universe(X) :- pf_1(X).
) :- universe(A), pf_1(B), universe(C), universe(D).
used_4(X,0) :- sequence_4(X, _, _, _). 
used_4(X,1) :- sequence_4(_, X, _, _). 
used_4(X,2) :- sequence_4(_, _, X, _). 
used_4(X,3) :- sequence_4(_, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
:- C = #count{N:used(S,N),df_5(S)}, C!=4.
sequence_5(A,B,C,D,Epf_1(X) :- pf_1_0(X).
pf_1_1(4,3).
pf_1(X) :- pf_1_1(X).
universe(X) :- pf_1(X).
pf_4_0(1, 1).
pf_4_0(2, 1).
pf_4_0(3, 1).
pf_4_0(7, 1).
pf_4_0(8, 1).
pf_4_0(9, 1).
pf_4_0(17, 1).
pf_4_0(18, 1).
pf_4_0(19, 1).
pf_4_0(20, 1).
pf_4(X) :- pf_4_0(X).
pf_4_1(4,3).
pf_4(X) :- pf_4_1(X).
universe(X) :- pf_4(X).
) :- universe(A), pf_1(B), universe(C), universe(D), pf_4(E).
used_5(X,0) :- sequence_5(X, _, _, _, _). 
used_5(X,1) :- sequence_5(_, X, _, _, _). 
used_5(X,2) :- sequence_5(_, _, X, _, _). 
used_5(X,3) :- sequence_5(_, _, _, X, _). 
used_5(X,4) :- sequence_5(_, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
:- C = #count{N:used(S,N),df_5(S)}, C!=4.
:- C = #count{N:used(S,N),df_5(S)}, C!=5.
sequence_6(A,B,C,D,E,Fpf_1(X) :- pf_1_0(X).
pf_1_1(4,3).
pf_1(X) :- pf_1_1(X).
universe(X) :- pf_1(X).
pf_4_0(1, 1).
pf_4_0(2, 1).
pf_4_0(3, 1).
pf_4_0(7, 1).
pf_4_0(8, 1).
pf_4_0(9, 1).
pf_4_0(17, 1).
pf_4_0(18, 1).
pf_4_0(19, 1).
pf_4_0(20, 1).
pf_4(X) :- pf_4_0(X).
pf_4_1(4,3).
pf_4(X) :- pf_4_1(X).
universe(X) :- pf_4(X).
) :- universe(A), pf_1(B), universe(C), universe(D), pf_4(E), universe(F).
used_6(X,0) :- sequence_6(X, _, _, _, _, _). 
used_6(X,1) :- sequence_6(_, X, _, _, _, _). 
used_6(X,2) :- sequence_6(_, _, X, _, _, _). 
used_6(X,3) :- sequence_6(_, _, _, X, _, _). 
used_6(X,4) :- sequence_6(_, _, _, _, X, _). 
used_6(X,5) :- sequence_6(_, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
:- C = #count{N:used(S,N),df_5(S)}, C!=4.
:- C = #count{N:used(S,N),df_5(S)}, C!=5.
:- C = #count{N:used(S,N),df_5(S)}, C!=6.
sequence_7(A,B,C,D,E,F,Gpf_1(X) :- pf_1_0(X).
pf_1_1(4,3).
pf_1(X) :- pf_1_1(X).
universe(X) :- pf_1(X).
pf_4_0(1, 1).
pf_4_0(2, 1).
pf_4_0(3, 1).
pf_4_0(7, 1).
pf_4_0(8, 1).
pf_4_0(9, 1).
pf_4_0(17, 1).
pf_4_0(18, 1).
pf_4_0(19, 1).
pf_4_0(20, 1).
pf_4(X) :- pf_4_0(X).
pf_4_1(4,3).
pf_4(X) :- pf_4_1(X).
universe(X) :- pf_4(X).
) :- universe(A), pf_1(B), universe(C), universe(D), pf_4(E), universe(F), universe(G).
used_7(X,0) :- sequence_7(X, _, _, _, _, _, _). 
used_7(X,1) :- sequence_7(_, X, _, _, _, _, _). 
used_7(X,2) :- sequence_7(_, _, X, _, _, _, _). 
used_7(X,3) :- sequence_7(_, _, _, X, _, _, _). 
used_7(X,4) :- sequence_7(_, _, _, _, X, _, _). 
used_7(X,5) :- sequence_7(_, _, _, _, _, X, _). 
used_7(X,6) :- sequence_7(_, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
:- C = #count{N:used(S,N),df_5(S)}, C!=4.
:- C = #count{N:used(S,N),df_5(S)}, C!=5.
:- C = #count{N:used(S,N),df_5(S)}, C!=6.
:- C = #count{N:used(S,N),df_5(S)}, C!=7.
sequence_8(A,B,C,D,E,F,G,Hpf_1(X) :- pf_1_0(X).
pf_1_1(4,3).
pf_1(X) :- pf_1_1(X).
universe(X) :- pf_1(X).
pf_4_0(1, 1).
pf_4_0(2, 1).
pf_4_0(3, 1).
pf_4_0(7, 1).
pf_4_0(8, 1).
pf_4_0(9, 1).
pf_4_0(17, 1).
pf_4_0(18, 1).
pf_4_0(19, 1).
pf_4_0(20, 1).
pf_4(X) :- pf_4_0(X).
pf_4_1(4,3).
pf_4(X) :- pf_4_1(X).
universe(X) :- pf_4(X).
) :- universe(A), pf_1(B), universe(C), universe(D), pf_4(E), universe(F), universe(G), universe(H).
used_8(X,0) :- sequence_8(X, _, _, _, _, _, _, _). 
used_8(X,1) :- sequence_8(_, X, _, _, _, _, _, _). 
used_8(X,2) :- sequence_8(_, _, X, _, _, _, _, _). 
used_8(X,3) :- sequence_8(_, _, _, X, _, _, _, _). 
used_8(X,4) :- sequence_8(_, _, _, _, X, _, _, _). 
used_8(X,5) :- sequence_8(_, _, _, _, _, X, _, _). 
used_8(X,6) :- sequence_8(_, _, _, _, _, _, X, _). 
used_8(X,7) :- sequence_8(_, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
:- C = #count{N:used(S,N),df_5(S)}, C!=4.
:- C = #count{N:used(S,N),df_5(S)}, C!=5.
:- C = #count{N:used(S,N),df_5(S)}, C!=6.
:- C = #count{N:used(S,N),df_5(S)}, C!=7.
:- C = #count{N:used(S,N),df_5(S)}, C!=8.
sequence_9(A,B,C,D,E,F,G,H,Ipf_1(X) :- pf_1_0(X).
pf_1_1(4,3).
pf_1(X) :- pf_1_1(X).
universe(X) :- pf_1(X).
pf_4_0(1, 1).
pf_4_0(2, 1).
pf_4_0(3, 1).
pf_4_0(7, 1).
pf_4_0(8, 1).
pf_4_0(9, 1).
pf_4_0(17, 1).
pf_4_0(18, 1).
pf_4_0(19, 1).
pf_4_0(20, 1).
pf_4(X) :- pf_4_0(X).
pf_4_1(4,3).
pf_4(X) :- pf_4_1(X).
universe(X) :- pf_4(X).
) :- universe(A), pf_1(B), universe(C), universe(D), pf_4(E), universe(F), universe(G), universe(H), universe(I).
used_9(X,0) :- sequence_9(X, _, _, _, _, _, _, _, _). 
used_9(X,1) :- sequence_9(_, X, _, _, _, _, _, _, _). 
used_9(X,2) :- sequence_9(_, _, X, _, _, _, _, _, _). 
used_9(X,3) :- sequence_9(_, _, _, X, _, _, _, _, _). 
used_9(X,4) :- sequence_9(_, _, _, _, X, _, _, _, _). 
used_9(X,5) :- sequence_9(_, _, _, _, _, X, _, _, _). 
used_9(X,6) :- sequence_9(_, _, _, _, _, _, X, _, _). 
used_9(X,7) :- sequence_9(_, _, _, _, _, _, _, X, _). 
used_9(X,8) :- sequence_9(_, _, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
:- C = #count{N:used(S,N),df_5(S)}, C!=4.
:- C = #count{N:used(S,N),df_5(S)}, C!=5.
:- C = #count{N:used(S,N),df_5(S)}, C!=6.
:- C = #count{N:used(S,N),df_5(S)}, C!=7.
:- C = #count{N:used(S,N),df_5(S)}, C!=8.
:- C = #count{N:used(S,N),df_5(S)}, C!=9.
sequence_10(A,B,C,D,E,F,G,H,I,Jpf_1(X) :- pf_1_0(X).
pf_1_1(4,3).
pf_1(X) :- pf_1_1(X).
universe(X) :- pf_1(X).
pf_4_0(1, 1).
pf_4_0(2, 1).
pf_4_0(3, 1).
pf_4_0(7, 1).
pf_4_0(8, 1).
pf_4_0(9, 1).
pf_4_0(17, 1).
pf_4_0(18, 1).
pf_4_0(19, 1).
pf_4_0(20, 1).
pf_4(X) :- pf_4_0(X).
pf_4_1(4,3).
pf_4(X) :- pf_4_1(X).
universe(X) :- pf_4(X).
) :- universe(A), pf_1(B), universe(C), universe(D), pf_4(E), universe(F), universe(G), universe(H), universe(I), universe(J).
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
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
:- C = #count{N:used(S,N),df_5(S)}, C!=4.
:- C = #count{N:used(S,N),df_5(S)}, C!=5.
:- C = #count{N:used(S,N),df_5(S)}, C!=6.
:- C = #count{N:used(S,N),df_5(S)}, C!=7.
:- C = #count{N:used(S,N),df_5(S)}, C!=8.
:- C = #count{N:used(S,N),df_5(S)}, C!=9.
:- C = #count{N:used(S,N),df_5(S)}, C!=10.
sequence_11(A,B,C,D,E,F,G,H,I,J,Kpf_1(X) :- pf_1_0(X).
pf_1_1(4,3).
pf_1(X) :- pf_1_1(X).
universe(X) :- pf_1(X).
pf_4_0(1, 1).
pf_4_0(2, 1).
pf_4_0(3, 1).
pf_4_0(7, 1).
pf_4_0(8, 1).
pf_4_0(9, 1).
pf_4_0(17, 1).
pf_4_0(18, 1).
pf_4_0(19, 1).
pf_4_0(20, 1).
pf_4(X) :- pf_4_0(X).
pf_4_1(4,3).
pf_4(X) :- pf_4_1(X).
universe(X) :- pf_4(X).
) :- universe(A), pf_1(B), universe(C), universe(D), pf_4(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K).
used_11(X,0) :- sequence_11(X, _, _, _, _, _, _, _, _, _, _). 
used_11(X,1) :- sequence_11(_, X, _, _, _, _, _, _, _, _, _). 
used_11(X,2) :- sequence_11(_, _, X, _, _, _, _, _, _, _, _). 
used_11(X,3) :- sequence_11(_, _, _, X, _, _, _, _, _, _, _). 
used_11(X,4) :- sequence_11(_, _, _, _, X, _, _, _, _, _, _). 
used_11(X,5) :- sequence_11(_, _, _, _, _, X, _, _, _, _, _). 
used_11(X,6) :- sequence_11(_, _, _, _, _, _, X, _, _, _, _). 
used_11(X,7) :- sequence_11(_, _, _, _, _, _, _, X, _, _, _). 
used_11(X,8) :- sequence_11(_, _, _, _, _, _, _, _, X, _, _). 
used_11(X,9) :- sequence_11(_, _, _, _, _, _, _, _, _, X, _). 
used_11(X,10) :- sequence_11(_, _, _, _, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.
:- C = #count{N:used(S,N),df_1(S)}, C!=11.
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
:- C = #count{N:used(S,N),df_5(S)}, C!=4.
:- C = #count{N:used(S,N),df_5(S)}, C!=5.
:- C = #count{N:used(S,N),df_5(S)}, C!=6.
:- C = #count{N:used(S,N),df_5(S)}, C!=7.
:- C = #count{N:used(S,N),df_5(S)}, C!=8.
:- C = #count{N:used(S,N),df_5(S)}, C!=9.
:- C = #count{N:used(S,N),df_5(S)}, C!=10.
:- C = #count{N:used(S,N),df_5(S)}, C!=11.
sequence_12(A,B,C,D,E,F,G,H,I,J,K,Lpf_1(X) :- pf_1_0(X).
pf_1_1(4,3).
pf_1(X) :- pf_1_1(X).
universe(X) :- pf_1(X).
pf_4_0(1, 1).
pf_4_0(2, 1).
pf_4_0(3, 1).
pf_4_0(7, 1).
pf_4_0(8, 1).
pf_4_0(9, 1).
pf_4_0(17, 1).
pf_4_0(18, 1).
pf_4_0(19, 1).
pf_4_0(20, 1).
pf_4(X) :- pf_4_0(X).
pf_4_1(4,3).
pf_4(X) :- pf_4_1(X).
universe(X) :- pf_4(X).
pf_11_0(3, 1).
pf_11_0(17, 1).
pf_11_0(18, 1).
pf_11_0(19, 1).
pf_11_0(20, 1).
pf_11(X) :- pf_11_0(X).
universe(X) :- pf_11(X).
) :- universe(A), pf_1(B), universe(C), universe(D), pf_4(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), pf_11(L).
used_12(X,0) :- sequence_12(X, _, _, _, _, _, _, _, _, _, _, _). 
used_12(X,1) :- sequence_12(_, X, _, _, _, _, _, _, _, _, _, _). 
used_12(X,2) :- sequence_12(_, _, X, _, _, _, _, _, _, _, _, _). 
used_12(X,3) :- sequence_12(_, _, _, X, _, _, _, _, _, _, _, _). 
used_12(X,4) :- sequence_12(_, _, _, _, X, _, _, _, _, _, _, _). 
used_12(X,5) :- sequence_12(_, _, _, _, _, X, _, _, _, _, _, _). 
used_12(X,6) :- sequence_12(_, _, _, _, _, _, X, _, _, _, _, _). 
used_12(X,7) :- sequence_12(_, _, _, _, _, _, _, X, _, _, _, _). 
used_12(X,8) :- sequence_12(_, _, _, _, _, _, _, _, X, _, _, _). 
used_12(X,9) :- sequence_12(_, _, _, _, _, _, _, _, _, X, _, _). 
used_12(X,10) :- sequence_12(_, _, _, _, _, _, _, _, _, _, X, _). 
used_12(X,11) :- sequence_12(_, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.
:- C = #count{N:used(S,N),df_1(S)}, C!=11.
:- C = #count{N:used(S,N),df_1(S)}, C!=12.
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
:- C = #count{N:used(S,N),df_5(S)}, C!=4.
:- C = #count{N:used(S,N),df_5(S)}, C!=5.
:- C = #count{N:used(S,N),df_5(S)}, C!=6.
:- C = #count{N:used(S,N),df_5(S)}, C!=7.
:- C = #count{N:used(S,N),df_5(S)}, C!=8.
:- C = #count{N:used(S,N),df_5(S)}, C!=9.
:- C = #count{N:used(S,N),df_5(S)}, C!=10.
:- C = #count{N:used(S,N),df_5(S)}, C!=11.
:- C = #count{N:used(S,N),df_5(S)}, C!=12.
sequence_13(A,B,C,D,E,F,G,H,I,J,K,L,Mpf_1(X) :- pf_1_0(X).
pf_1_1(4,3).
pf_1(X) :- pf_1_1(X).
universe(X) :- pf_1(X).
pf_4_0(1, 1).
pf_4_0(2, 1).
pf_4_0(3, 1).
pf_4_0(7, 1).
pf_4_0(8, 1).
pf_4_0(9, 1).
pf_4_0(17, 1).
pf_4_0(18, 1).
pf_4_0(19, 1).
pf_4_0(20, 1).
pf_4(X) :- pf_4_0(X).
pf_4_1(4,3).
pf_4(X) :- pf_4_1(X).
universe(X) :- pf_4(X).
pf_11_0(3, 1).
pf_11_0(17, 1).
pf_11_0(18, 1).
pf_11_0(19, 1).
pf_11_0(20, 1).
pf_11(X) :- pf_11_0(X).
universe(X) :- pf_11(X).
pf_12_0(2, 1).
pf_12_0(15, 1).
pf_12_0(16, 1).
pf_12(X) :- pf_12_0(X).
universe(X) :- pf_12(X).
) :- universe(A), pf_1(B), universe(C), universe(D), pf_4(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), pf_11(L), pf_12(M).
used_13(X,0) :- sequence_13(X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_13(X,1) :- sequence_13(_, X, _, _, _, _, _, _, _, _, _, _, _). 
used_13(X,2) :- sequence_13(_, _, X, _, _, _, _, _, _, _, _, _, _). 
used_13(X,3) :- sequence_13(_, _, _, X, _, _, _, _, _, _, _, _, _). 
used_13(X,4) :- sequence_13(_, _, _, _, X, _, _, _, _, _, _, _, _). 
used_13(X,5) :- sequence_13(_, _, _, _, _, X, _, _, _, _, _, _, _). 
used_13(X,6) :- sequence_13(_, _, _, _, _, _, X, _, _, _, _, _, _). 
used_13(X,7) :- sequence_13(_, _, _, _, _, _, _, X, _, _, _, _, _). 
used_13(X,8) :- sequence_13(_, _, _, _, _, _, _, _, X, _, _, _, _). 
used_13(X,9) :- sequence_13(_, _, _, _, _, _, _, _, _, X, _, _, _). 
used_13(X,10) :- sequence_13(_, _, _, _, _, _, _, _, _, _, X, _, _). 
used_13(X,11) :- sequence_13(_, _, _, _, _, _, _, _, _, _, _, X, _). 
used_13(X,12) :- sequence_13(_, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.
:- C = #count{N:used(S,N),df_1(S)}, C!=11.
:- C = #count{N:used(S,N),df_1(S)}, C!=12.
:- C = #count{N:used(S,N),df_1(S)}, C!=13.
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
:- C = #count{N:used(S,N),df_5(S)}, C!=4.
:- C = #count{N:used(S,N),df_5(S)}, C!=5.
:- C = #count{N:used(S,N),df_5(S)}, C!=6.
:- C = #count{N:used(S,N),df_5(S)}, C!=7.
:- C = #count{N:used(S,N),df_5(S)}, C!=8.
:- C = #count{N:used(S,N),df_5(S)}, C!=9.
:- C = #count{N:used(S,N),df_5(S)}, C!=10.
:- C = #count{N:used(S,N),df_5(S)}, C!=11.
:- C = #count{N:used(S,N),df_5(S)}, C!=12.
:- C = #count{N:used(S,N),df_5(S)}, C!=13.
sequence_14(A,B,C,D,E,F,G,H,I,J,K,L,M,Npf_1(X) :- pf_1_0(X).
pf_1_1(4,3).
pf_1(X) :- pf_1_1(X).
universe(X) :- pf_1(X).
pf_4_0(1, 1).
pf_4_0(2, 1).
pf_4_0(3, 1).
pf_4_0(7, 1).
pf_4_0(8, 1).
pf_4_0(9, 1).
pf_4_0(17, 1).
pf_4_0(18, 1).
pf_4_0(19, 1).
pf_4_0(20, 1).
pf_4(X) :- pf_4_0(X).
pf_4_1(4,3).
pf_4(X) :- pf_4_1(X).
universe(X) :- pf_4(X).
pf_11_0(3, 1).
pf_11_0(17, 1).
pf_11_0(18, 1).
pf_11_0(19, 1).
pf_11_0(20, 1).
pf_11(X) :- pf_11_0(X).
universe(X) :- pf_11(X).
pf_12_0(2, 1).
pf_12_0(15, 1).
pf_12_0(16, 1).
pf_12(X) :- pf_12_0(X).
universe(X) :- pf_12(X).
pf_13_0(1, 1).
pf_13_0(7, 1).
pf_13_0(8, 1).
pf_13_0(9, 1).
pf_13_0(10, 1).
pf_13_0(11, 1).
pf_13_0(12, 1).
pf_13_0(13, 1).
pf_13_0(14, 1).
pf_13_0(17, 1).
pf_13_0(18, 1).
pf_13(X) :- pf_13_0(X).
pf_13_1(4,3).
pf_13(X) :- pf_13_1(X).
universe(X) :- pf_13(X).
) :- universe(A), pf_1(B), universe(C), universe(D), pf_4(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), pf_11(L), pf_12(M), pf_13(N).
used_14(X,0) :- sequence_14(X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,1) :- sequence_14(_, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,2) :- sequence_14(_, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,3) :- sequence_14(_, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_14(X,4) :- sequence_14(_, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_14(X,5) :- sequence_14(_, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_14(X,6) :- sequence_14(_, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_14(X,7) :- sequence_14(_, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_14(X,8) :- sequence_14(_, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_14(X,9) :- sequence_14(_, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_14(X,10) :- sequence_14(_, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_14(X,11) :- sequence_14(_, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_14(X,12) :- sequence_14(_, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_14(X,13) :- sequence_14(_, _, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.
:- C = #count{N:used(S,N),df_1(S)}, C!=11.
:- C = #count{N:used(S,N),df_1(S)}, C!=12.
:- C = #count{N:used(S,N),df_1(S)}, C!=13.
:- C = #count{N:used(S,N),df_1(S)}, C!=14.
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
:- C = #count{N:used(S,N),df_5(S)}, C!=4.
:- C = #count{N:used(S,N),df_5(S)}, C!=5.
:- C = #count{N:used(S,N),df_5(S)}, C!=6.
:- C = #count{N:used(S,N),df_5(S)}, C!=7.
:- C = #count{N:used(S,N),df_5(S)}, C!=8.
:- C = #count{N:used(S,N),df_5(S)}, C!=9.
:- C = #count{N:used(S,N),df_5(S)}, C!=10.
:- C = #count{N:used(S,N),df_5(S)}, C!=11.
:- C = #count{N:used(S,N),df_5(S)}, C!=12.
:- C = #count{N:used(S,N),df_5(S)}, C!=13.
:- C = #count{N:used(S,N),df_5(S)}, C!=14.
sequence_15(A,B,C,D,E,F,G,H,I,J,K,L,M,N,Opf_1(X) :- pf_1_0(X).
pf_1_1(4,3).
pf_1(X) :- pf_1_1(X).
universe(X) :- pf_1(X).
pf_4_0(1, 1).
pf_4_0(2, 1).
pf_4_0(3, 1).
pf_4_0(7, 1).
pf_4_0(8, 1).
pf_4_0(9, 1).
pf_4_0(17, 1).
pf_4_0(18, 1).
pf_4_0(19, 1).
pf_4_0(20, 1).
pf_4(X) :- pf_4_0(X).
pf_4_1(4,3).
pf_4(X) :- pf_4_1(X).
universe(X) :- pf_4(X).
pf_11_0(3, 1).
pf_11_0(17, 1).
pf_11_0(18, 1).
pf_11_0(19, 1).
pf_11_0(20, 1).
pf_11(X) :- pf_11_0(X).
universe(X) :- pf_11(X).
pf_12_0(2, 1).
pf_12_0(15, 1).
pf_12_0(16, 1).
pf_12(X) :- pf_12_0(X).
universe(X) :- pf_12(X).
pf_13_0(1, 1).
pf_13_0(7, 1).
pf_13_0(8, 1).
pf_13_0(9, 1).
pf_13_0(10, 1).
pf_13_0(11, 1).
pf_13_0(12, 1).
pf_13_0(13, 1).
pf_13_0(14, 1).
pf_13_0(17, 1).
pf_13_0(18, 1).
pf_13(X) :- pf_13_0(X).
pf_13_1(4,3).
pf_13(X) :- pf_13_1(X).
universe(X) :- pf_13(X).
) :- universe(A), pf_1(B), universe(C), universe(D), pf_4(E), universe(F), universe(G), universe(H), universe(I), universe(J), universe(K), pf_11(L), pf_12(M), pf_13(N), universe(O).
used_15(X,0) :- sequence_15(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,1) :- sequence_15(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,2) :- sequence_15(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,3) :- sequence_15(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,4) :- sequence_15(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_15(X,5) :- sequence_15(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_15(X,6) :- sequence_15(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_15(X,7) :- sequence_15(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_15(X,8) :- sequence_15(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_15(X,9) :- sequence_15(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_15(X,10) :- sequence_15(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_15(X,11) :- sequence_15(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_15(X,12) :- sequence_15(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_15(X,13) :- sequence_15(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_15(X,14) :- sequence_15(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
df_0_0(1, 1).
df_0_0(2, 1).
df_0_0(3, 1).
df_0_0(7, 1).
df_0_0(8, 1).
df_0_0(9, 1).
df_0_0(15, 1).
df_0_0(16, 1).
df_0_0(17, 1).
df_0_0(18, 1).
df_0_0(19, 1).
df_0_0(20, 1).
df_0(X) :- df_0_0(X).
df_0_1(4,3).
df_0(X) :- df_0_1(X).
universe(X) :- df_0(X).
:- C = #count{N:used(S,N),df_0(S)}, C!=0.
:- C = #count{N:used(S,N),df_0(S)}, C!=1.
:- C = #count{N:used(S,N),df_0(S)}, C!=2.
:- C = #count{N:used(S,N),df_0(S)}, C!=3.
:- C = #count{N:used(S,N),df_0(S)}, C!=4.
:- C = #count{N:used(S,N),df_0(S)}, C!=5.
:- C = #count{N:used(S,N),df_0(S)}, C!=6.
:- C = #count{N:used(S,N),df_0(S)}, C!=7.
:- C = #count{N:used(S,N),df_0(S)}, C!=8.
:- C = #count{N:used(S,N),df_0(S)}, C!=9.
:- C = #count{N:used(S,N),df_0(S)}, C!=10.
:- C = #count{N:used(S,N),df_0(S)}, C!=11.
:- C = #count{N:used(S,N),df_0(S)}, C!=12.
df_1(X) :- df_1_0(X).
df_1_1(4,3).
df_1(X) :- df_1_1(X).
universe(X) :- df_1(X).
:- C = #count{N:used(S,N),df_1(S)}, C!=8.
:- C = #count{N:used(S,N),df_1(S)}, C!=9.
:- C = #count{N:used(S,N),df_1(S)}, C!=10.
:- C = #count{N:used(S,N),df_1(S)}, C!=11.
:- C = #count{N:used(S,N),df_1(S)}, C!=12.
:- C = #count{N:used(S,N),df_1(S)}, C!=13.
:- C = #count{N:used(S,N),df_1(S)}, C!=14.
:- C = #count{N:used(S,N),df_1(S)}, C!=15.
df_2(X) :- df_2_0(X).
df_2_1(4,3).
df_2(X) :- df_2_1(X).
universe(X) :- df_2(X).
:- C = #count{N:used(S,N),df_2(S)}, C!=15.
df_3_0(1, 1).
df_3_0(3, 1).
df_3_0(7, 1).
df_3_0(8, 1).
df_3_0(9, 1).
df_3_0(10, 1).
df_3_0(11, 1).
df_3_0(12, 1).
df_3_0(13, 1).
df_3_0(14, 1).
df_3_0(17, 1).
df_3_0(18, 1).
df_3_0(19, 1).
df_3_0(20, 1).
df_3(X) :- df_3_0(X).
df_3_1(4,3).
df_3(X) :- df_3_1(X).
universe(X) :- df_3(X).
df_4_0(1, 1).
df_4_0(2, 1).
df_4_0(7, 1).
df_4_0(8, 1).
df_4_0(9, 1).
df_4_0(10, 1).
df_4_0(11, 1).
df_4_0(12, 1).
df_4_0(13, 1).
df_4_0(14, 1).
df_4_0(15, 1).
df_4_0(16, 1).
df_4_0(17, 1).
df_4_0(18, 1).
df_4(X) :- df_4_0(X).
df_4_1(4,3).
df_4(X) :- df_4_1(X).
universe(X) :- df_4(X).
:- C = #count{N:used(S,N),df_4(S)}, C!=0.
:- C = #count{N:used(S,N),df_4(S)}, C!=1.
:- C = #count{N:used(S,N),df_4(S)}, C!=2.
:- C = #count{N:used(S,N),df_4(S)}, C!=3.
:- C = #count{N:used(S,N),df_4(S)}, C!=4.
:- C = #count{N:used(S,N),df_4(S)}, C!=5.
df_5_0(3, 1).
df_5_0(19, 1).
df_5_0(20, 1).
df_5(X) :- df_5_0(X).
universe(X) :- df_5(X).
:- C = #count{N:used(S,N),df_5(S)}, C!=4.
:- C = #count{N:used(S,N),df_5(S)}, C!=5.
:- C = #count{N:used(S,N),df_5(S)}, C!=6.
:- C = #count{N:used(S,N),df_5(S)}, C!=7.
:- C = #count{N:used(S,N),df_5(S)}, C!=8.
:- C = #count{N:used(S,N),df_5(S)}, C!=9.
:- C = #count{N:used(S,N),df_5(S)}, C!=10.
:- C = #count{N:used(S,N),df_5(S)}, C!=11.
:- C = #count{N:used(S,N),df_5(S)}, C!=12.
:- C = #count{N:used(S,N),df_5(S)}, C!=13.
:- C = #count{N:used(S,N),df_5(S)}, C!=14.
:- C = #count{N:used(S,N),df_5(S)}, C!=15.

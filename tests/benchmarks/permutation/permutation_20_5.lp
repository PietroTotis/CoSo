universe_0(1, 1).
universe_0(2, 1).
universe_0(3, 1).
universe_0(12, 1).
universe_0(13, 1).
universe_0(14, 1).
universe_0(15, 1).
universe_0(16, 1).
universe_0(17, 1).
universe_0(18, 1).
universe(X) :- universe_0(X).
universe_1(4,17).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(1, 1).
dom1_0(12, 1).
dom1_0(13, 1).
dom1_0(14, 1).
dom1_0(15, 1).
dom1_0(16, 1).
dom1_0(17, 1).
dom1_0(18, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
permutation_1(A) :- universe(A).
used_1(X,0) :- permutation_1(X). 
:- permutation_1(A), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_2(A,B) :- universe(A), universe(B).
used_2(X,0) :- permutation_2(X, _). 
used_2(X,1) :- permutation_2(_, X). 
:- permutation_2(A,B), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_3(A,B,C) :- universe(A), universe(B), universe(C).
used_3(X,0) :- permutation_3(X, _, _). 
used_3(X,1) :- permutation_3(_, X, _). 
used_3(X,2) :- permutation_3(_, _, X). 
:- permutation_3(A,B,C), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_4(A,B,C,D) :- universe(A), universe(B), universe(C), universe(D).
used_4(X,0) :- permutation_4(X, _, _, _). 
used_4(X,1) :- permutation_4(_, X, _, _). 
used_4(X,2) :- permutation_4(_, _, X, _). 
used_4(X,3) :- permutation_4(_, _, _, X). 
:- permutation_4(A,B,C,D), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_6(A,B,C,D,E,Fpf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F).
used_6(X,0) :- permutation_6(X, _, _, _, _, _). 
used_6(X,1) :- permutation_6(_, X, _, _, _, _). 
used_6(X,2) :- permutation_6(_, _, X, _, _, _). 
used_6(X,3) :- permutation_6(_, _, _, X, _, _). 
used_6(X,4) :- permutation_6(_, _, _, _, X, _). 
used_6(X,5) :- permutation_6(_, _, _, _, _, X). 
:- permutation_6(A,B,C,D,E,F), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_7(A,B,C,D,E,F,Gpf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F), universe(G).
used_7(X,0) :- permutation_7(X, _, _, _, _, _, _). 
used_7(X,1) :- permutation_7(_, X, _, _, _, _, _). 
used_7(X,2) :- permutation_7(_, _, X, _, _, _, _). 
used_7(X,3) :- permutation_7(_, _, _, X, _, _, _). 
used_7(X,4) :- permutation_7(_, _, _, _, X, _, _). 
used_7(X,5) :- permutation_7(_, _, _, _, _, X, _). 
used_7(X,6) :- permutation_7(_, _, _, _, _, _, X). 
:- permutation_7(A,B,C,D,E,F,G), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_8(A,B,C,D,E,F,G,Hpf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F), universe(G), universe(H).
used_8(X,0) :- permutation_8(X, _, _, _, _, _, _, _). 
used_8(X,1) :- permutation_8(_, X, _, _, _, _, _, _). 
used_8(X,2) :- permutation_8(_, _, X, _, _, _, _, _). 
used_8(X,3) :- permutation_8(_, _, _, X, _, _, _, _). 
used_8(X,4) :- permutation_8(_, _, _, _, X, _, _, _). 
used_8(X,5) :- permutation_8(_, _, _, _, _, X, _, _). 
used_8(X,6) :- permutation_8(_, _, _, _, _, _, X, _). 
used_8(X,7) :- permutation_8(_, _, _, _, _, _, _, X). 
:- permutation_8(A,B,C,D,E,F,G,H), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_9(A,B,C,D,E,F,G,H,Ipf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F), universe(G), universe(H), universe(I).
used_9(X,0) :- permutation_9(X, _, _, _, _, _, _, _, _). 
used_9(X,1) :- permutation_9(_, X, _, _, _, _, _, _, _). 
used_9(X,2) :- permutation_9(_, _, X, _, _, _, _, _, _). 
used_9(X,3) :- permutation_9(_, _, _, X, _, _, _, _, _). 
used_9(X,4) :- permutation_9(_, _, _, _, X, _, _, _, _). 
used_9(X,5) :- permutation_9(_, _, _, _, _, X, _, _, _). 
used_9(X,6) :- permutation_9(_, _, _, _, _, _, X, _, _). 
used_9(X,7) :- permutation_9(_, _, _, _, _, _, _, X, _). 
used_9(X,8) :- permutation_9(_, _, _, _, _, _, _, _, X). 
:- permutation_9(A,B,C,D,E,F,G,H,I), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_10(A,B,C,D,E,F,G,H,I,Jpf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F), universe(G), universe(H), universe(I), universe(J).
used_10(X,0) :- permutation_10(X, _, _, _, _, _, _, _, _, _). 
used_10(X,1) :- permutation_10(_, X, _, _, _, _, _, _, _, _). 
used_10(X,2) :- permutation_10(_, _, X, _, _, _, _, _, _, _). 
used_10(X,3) :- permutation_10(_, _, _, X, _, _, _, _, _, _). 
used_10(X,4) :- permutation_10(_, _, _, _, X, _, _, _, _, _). 
used_10(X,5) :- permutation_10(_, _, _, _, _, X, _, _, _, _). 
used_10(X,6) :- permutation_10(_, _, _, _, _, _, X, _, _, _). 
used_10(X,7) :- permutation_10(_, _, _, _, _, _, _, X, _, _). 
used_10(X,8) :- permutation_10(_, _, _, _, _, _, _, _, X, _). 
used_10(X,9) :- permutation_10(_, _, _, _, _, _, _, _, _, X). 
:- permutation_10(A,B,C,D,E,F,G,H,I,J), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_11(A,B,C,D,E,F,G,H,I,J,Kpf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F), universe(G), universe(H), universe(I), universe(J), universe(K).
used_11(X,0) :- permutation_11(X, _, _, _, _, _, _, _, _, _, _). 
used_11(X,1) :- permutation_11(_, X, _, _, _, _, _, _, _, _, _). 
used_11(X,2) :- permutation_11(_, _, X, _, _, _, _, _, _, _, _). 
used_11(X,3) :- permutation_11(_, _, _, X, _, _, _, _, _, _, _). 
used_11(X,4) :- permutation_11(_, _, _, _, X, _, _, _, _, _, _). 
used_11(X,5) :- permutation_11(_, _, _, _, _, X, _, _, _, _, _). 
used_11(X,6) :- permutation_11(_, _, _, _, _, _, X, _, _, _, _). 
used_11(X,7) :- permutation_11(_, _, _, _, _, _, _, X, _, _, _). 
used_11(X,8) :- permutation_11(_, _, _, _, _, _, _, _, X, _, _). 
used_11(X,9) :- permutation_11(_, _, _, _, _, _, _, _, _, X, _). 
used_11(X,10) :- permutation_11(_, _, _, _, _, _, _, _, _, _, X). 
:- permutation_11(A,B,C,D,E,F,G,H,I,J,K), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_12(A,B,C,D,E,F,G,H,I,J,K,Lpf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L).
used_12(X,0) :- permutation_12(X, _, _, _, _, _, _, _, _, _, _, _). 
used_12(X,1) :- permutation_12(_, X, _, _, _, _, _, _, _, _, _, _). 
used_12(X,2) :- permutation_12(_, _, X, _, _, _, _, _, _, _, _, _). 
used_12(X,3) :- permutation_12(_, _, _, X, _, _, _, _, _, _, _, _). 
used_12(X,4) :- permutation_12(_, _, _, _, X, _, _, _, _, _, _, _). 
used_12(X,5) :- permutation_12(_, _, _, _, _, X, _, _, _, _, _, _). 
used_12(X,6) :- permutation_12(_, _, _, _, _, _, X, _, _, _, _, _). 
used_12(X,7) :- permutation_12(_, _, _, _, _, _, _, X, _, _, _, _). 
used_12(X,8) :- permutation_12(_, _, _, _, _, _, _, _, X, _, _, _). 
used_12(X,9) :- permutation_12(_, _, _, _, _, _, _, _, _, X, _, _). 
used_12(X,10) :- permutation_12(_, _, _, _, _, _, _, _, _, _, X, _). 
used_12(X,11) :- permutation_12(_, _, _, _, _, _, _, _, _, _, _, X). 
:- permutation_12(A,B,C,D,E,F,G,H,I,J,K,L), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_13(A,B,C,D,E,F,G,H,I,J,K,L,Mpf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M).
used_13(X,0) :- permutation_13(X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_13(X,1) :- permutation_13(_, X, _, _, _, _, _, _, _, _, _, _, _). 
used_13(X,2) :- permutation_13(_, _, X, _, _, _, _, _, _, _, _, _, _). 
used_13(X,3) :- permutation_13(_, _, _, X, _, _, _, _, _, _, _, _, _). 
used_13(X,4) :- permutation_13(_, _, _, _, X, _, _, _, _, _, _, _, _). 
used_13(X,5) :- permutation_13(_, _, _, _, _, X, _, _, _, _, _, _, _). 
used_13(X,6) :- permutation_13(_, _, _, _, _, _, X, _, _, _, _, _, _). 
used_13(X,7) :- permutation_13(_, _, _, _, _, _, _, X, _, _, _, _, _). 
used_13(X,8) :- permutation_13(_, _, _, _, _, _, _, _, X, _, _, _, _). 
used_13(X,9) :- permutation_13(_, _, _, _, _, _, _, _, _, X, _, _, _). 
used_13(X,10) :- permutation_13(_, _, _, _, _, _, _, _, _, _, X, _, _). 
used_13(X,11) :- permutation_13(_, _, _, _, _, _, _, _, _, _, _, X, _). 
used_13(X,12) :- permutation_13(_, _, _, _, _, _, _, _, _, _, _, _, X). 
:- permutation_13(A,B,C,D,E,F,G,H,I,J,K,L,M), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_14(A,B,C,D,E,F,G,H,I,J,K,L,M,Npf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N).
used_14(X,0) :- permutation_14(X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,1) :- permutation_14(_, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,2) :- permutation_14(_, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_14(X,3) :- permutation_14(_, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_14(X,4) :- permutation_14(_, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_14(X,5) :- permutation_14(_, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_14(X,6) :- permutation_14(_, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_14(X,7) :- permutation_14(_, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_14(X,8) :- permutation_14(_, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_14(X,9) :- permutation_14(_, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_14(X,10) :- permutation_14(_, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_14(X,11) :- permutation_14(_, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_14(X,12) :- permutation_14(_, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_14(X,13) :- permutation_14(_, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- permutation_14(A,B,C,D,E,F,G,H,I,J,K,L,M,N), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_15(A,B,C,D,E,F,G,H,I,J,K,L,M,N,Opf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O).
used_15(X,0) :- permutation_15(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,1) :- permutation_15(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,2) :- permutation_15(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,3) :- permutation_15(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_15(X,4) :- permutation_15(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_15(X,5) :- permutation_15(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_15(X,6) :- permutation_15(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_15(X,7) :- permutation_15(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_15(X,8) :- permutation_15(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_15(X,9) :- permutation_15(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_15(X,10) :- permutation_15(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_15(X,11) :- permutation_15(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_15(X,12) :- permutation_15(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_15(X,13) :- permutation_15(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_15(X,14) :- permutation_15(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- permutation_15(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_16(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,Ppf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P).
used_16(X,0) :- permutation_16(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,1) :- permutation_16(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,2) :- permutation_16(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,3) :- permutation_16(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,4) :- permutation_16(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_16(X,5) :- permutation_16(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_16(X,6) :- permutation_16(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_16(X,7) :- permutation_16(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_16(X,8) :- permutation_16(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_16(X,9) :- permutation_16(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_16(X,10) :- permutation_16(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_16(X,11) :- permutation_16(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_16(X,12) :- permutation_16(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_16(X,13) :- permutation_16(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_16(X,14) :- permutation_16(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_16(X,15) :- permutation_16(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- permutation_16(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_17(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Qpf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q).
used_17(X,0) :- permutation_17(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,1) :- permutation_17(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,2) :- permutation_17(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,3) :- permutation_17(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,4) :- permutation_17(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,5) :- permutation_17(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_17(X,6) :- permutation_17(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_17(X,7) :- permutation_17(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_17(X,8) :- permutation_17(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_17(X,9) :- permutation_17(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_17(X,10) :- permutation_17(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_17(X,11) :- permutation_17(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_17(X,12) :- permutation_17(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_17(X,13) :- permutation_17(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_17(X,14) :- permutation_17(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_17(X,15) :- permutation_17(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_17(X,16) :- permutation_17(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- permutation_17(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_18(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,Rpf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q), universe(R).
used_18(X,0) :- permutation_18(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,1) :- permutation_18(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,2) :- permutation_18(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,3) :- permutation_18(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,4) :- permutation_18(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,5) :- permutation_18(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,6) :- permutation_18(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_18(X,7) :- permutation_18(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_18(X,8) :- permutation_18(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_18(X,9) :- permutation_18(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_18(X,10) :- permutation_18(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_18(X,11) :- permutation_18(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_18(X,12) :- permutation_18(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_18(X,13) :- permutation_18(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_18(X,14) :- permutation_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_18(X,15) :- permutation_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_18(X,16) :- permutation_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_18(X,17) :- permutation_18(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- permutation_18(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_19(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,Spf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q), universe(R), universe(S).
used_19(X,0) :- permutation_19(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,1) :- permutation_19(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,2) :- permutation_19(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,3) :- permutation_19(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,4) :- permutation_19(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,5) :- permutation_19(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,6) :- permutation_19(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,7) :- permutation_19(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_19(X,8) :- permutation_19(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_19(X,9) :- permutation_19(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_19(X,10) :- permutation_19(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_19(X,11) :- permutation_19(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_19(X,12) :- permutation_19(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_19(X,13) :- permutation_19(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_19(X,14) :- permutation_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_19(X,15) :- permutation_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_19(X,16) :- permutation_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_19(X,17) :- permutation_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_19(X,18) :- permutation_19(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- permutation_19(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S), set(S,SN), C = #count{N:used(S,N)}, C>SN.
permutation_20(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,Tpf_5_0(1, 1).
pf_5_0(12, 1).
pf_5_0(13, 1).
pf_5_0(14, 1).
pf_5_0(15, 1).
pf_5_0(16, 1).
pf_5_0(17, 1).
pf_5_0(18, 1).
pf_5(X) :- pf_5_0(X).
universe(X) :- pf_5(X).
) :- universe(A), universe(B), universe(C), universe(D), universe(E), pf_5(F), universe(G), universe(H), universe(I), universe(J), universe(K), universe(L), universe(M), universe(N), universe(O), universe(P), universe(Q), universe(R), universe(S), universe(T).
used_20(X,0) :- permutation_20(X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,1) :- permutation_20(_, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,2) :- permutation_20(_, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,3) :- permutation_20(_, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,4) :- permutation_20(_, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,5) :- permutation_20(_, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,6) :- permutation_20(_, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,7) :- permutation_20(_, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,8) :- permutation_20(_, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _, _). 
used_20(X,9) :- permutation_20(_, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _, _). 
used_20(X,10) :- permutation_20(_, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _, _). 
used_20(X,11) :- permutation_20(_, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _, _). 
used_20(X,12) :- permutation_20(_, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _, _). 
used_20(X,13) :- permutation_20(_, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _, _). 
used_20(X,14) :- permutation_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _, _). 
used_20(X,15) :- permutation_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _, _). 
used_20(X,16) :- permutation_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _, _). 
used_20(X,17) :- permutation_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _, _). 
used_20(X,18) :- permutation_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X, _). 
used_20(X,19) :- permutation_20(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, X). 
:- permutation_20(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T), set(S,SN), C = #count{N:used(S,N)}, C>SN.

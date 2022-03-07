universe_0(1, 1).
universe_0(2, 1).
universe_0(3, 1).
universe_0(7, 1).
universe_0(8, 1).
universe_0(9, 1).
universe_0(10, 1).
universe(X) :- universe_0(X).
universe_1(4,3).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(3, 1).
dom1_0(7, 1).
dom1_0(8, 1).
dom1_0(9, 1).
dom1_0(10, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
sequence_5(A,B,C,D,E) :- universe(A), universe(B), universe(C), universe(D), universe(E).
used_5(X,0) :- sequence_5(X, _, _, _, _). 
used_5(X,1) :- sequence_5(_, X, _, _, _). 
used_5(X,2) :- sequence_5(_, _, X, _, _). 
used_5(X,3) :- sequence_5(_, _, _, X, _). 
used_5(X,4) :- sequence_5(_, _, _, _, X). 
sequence_6(A,B,C,D,E,F) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F).
used_6(X,0) :- sequence_6(X, _, _, _, _, _). 
used_6(X,1) :- sequence_6(_, X, _, _, _, _). 
used_6(X,2) :- sequence_6(_, _, X, _, _, _). 
used_6(X,3) :- sequence_6(_, _, _, X, _, _). 
used_6(X,4) :- sequence_6(_, _, _, _, X, _). 
used_6(X,5) :- sequence_6(_, _, _, _, _, X). 
sequence_7(A,B,C,D,E,F,G) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G).
used_7(X,0) :- sequence_7(X, _, _, _, _, _, _). 
used_7(X,1) :- sequence_7(_, X, _, _, _, _, _). 
used_7(X,2) :- sequence_7(_, _, X, _, _, _, _). 
used_7(X,3) :- sequence_7(_, _, _, X, _, _, _). 
used_7(X,4) :- sequence_7(_, _, _, _, X, _, _). 
used_7(X,5) :- sequence_7(_, _, _, _, _, X, _). 
used_7(X,6) :- sequence_7(_, _, _, _, _, _, X). 
sequence_8(A,B,C,D,E,F,G,H) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H).
used_8(X,0) :- sequence_8(X, _, _, _, _, _, _, _). 
used_8(X,1) :- sequence_8(_, X, _, _, _, _, _, _). 
used_8(X,2) :- sequence_8(_, _, X, _, _, _, _, _). 
used_8(X,3) :- sequence_8(_, _, _, X, _, _, _, _). 
used_8(X,4) :- sequence_8(_, _, _, _, X, _, _, _). 
used_8(X,5) :- sequence_8(_, _, _, _, _, X, _, _). 
used_8(X,6) :- sequence_8(_, _, _, _, _, _, X, _). 
used_8(X,7) :- sequence_8(_, _, _, _, _, _, _, X). 
sequence_9(A,B,C,D,E,F,G,H,I) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I).
used_9(X,0) :- sequence_9(X, _, _, _, _, _, _, _, _). 
used_9(X,1) :- sequence_9(_, X, _, _, _, _, _, _, _). 
used_9(X,2) :- sequence_9(_, _, X, _, _, _, _, _, _). 
used_9(X,3) :- sequence_9(_, _, _, X, _, _, _, _, _). 
used_9(X,4) :- sequence_9(_, _, _, _, X, _, _, _, _). 
used_9(X,5) :- sequence_9(_, _, _, _, _, X, _, _, _). 
used_9(X,6) :- sequence_9(_, _, _, _, _, _, X, _, _). 
used_9(X,7) :- sequence_9(_, _, _, _, _, _, _, X, _). 
used_9(X,8) :- sequence_9(_, _, _, _, _, _, _, _, X). 
sequence_10(A,B,C,D,E,F,G,H,I,J) :- universe(A), universe(B), universe(C), universe(D), universe(E), universe(F), universe(G), universe(H), universe(I), universe(J).
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

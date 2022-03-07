universe_0(1, 1).
universe_0(2, 1).
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
universe(X) :- universe_0(X).
universe_1(3,3).
universe(X) :- universe_1(X).
universe(X) :- universe(X).
dom1_0(8, 1).
dom1_0(9, 1).
dom1_0(10, 1).
dom1_0(11, 1).
dom1_0(12, 1).
dom1_0(13, 1).
dom1_0(14, 1).
dom1_0(15, 1).
dom1(X) :- dom1_0(X).
universe(X) :- dom1(X).
dom2_0(1, 1).
dom2_0(6, 1).
dom2_0(7, 1).
dom2_0(8, 1).
dom2_0(9, 1).
dom2_0(10, 1).
dom2_0(11, 1).
dom2_0(12, 1).
dom2(X) :- dom2_0(X).
universe(X) :- dom2(X).
dom3_0(11, 1).
dom3_0(12, 1).
dom3_0(13, 1).
dom3_0(14, 1).
dom3_0(15, 1).
dom3(X) :- dom3_0(X).
universe(X) :- dom3(X).

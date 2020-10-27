:- use_module(library(aproblog)).
            :- use_semiring(
                sr_plus,   % addition (arity 3)
                sr_times,  % multiplication (arity 3)
                sr_zero,   % neutral element of addition
                sr_one,    % neutral element of multiplication
                sr_pos,
                sr_neg,    % negation of fact label
                sr_negate, 
                true,      % requires solving disjoint sum problem?
                true).    % requires solving neutral sum problem?
            sr_zero(0.0).
            sr_one(1.0).
            sr_plus(A, B, C) :- C is A + B.
            sr_times(A, B, C) :- C is A * B.
            sr_pos(A, B, A).
            sr_neg(A, B, 0).
            sr_negate(A, 1.0).
1::v1(1);1::v1(2);1::v1(3);1::v1(4);1::v1(5);1::v1(6);1::v1(7);1::v1(8);1::v1(9);1::v1(10);1::v1(11);1::v1(12);1::v1(13);1::v1(14).
1::v2(1);1::v2(2);1::v2(3);1::v2(4);1::v2(5);1::v2(6);1::v2(7);1::v2(8);1::v2(9);1::v2(10);1::v2(11);1::v2(12);1::v2(13);1::v2(14).
1::v3(1);1::v3(2);1::v3(3);1::v3(4);1::v3(5);1::v3(6);1::v3(7);1::v3(8);1::v3(9);1::v3(10);1::v3(11);1::v3(12);1::v3(13);1::v3(14).
1::v4(1);1::v4(2);1::v4(3);1::v4(4);1::v4(5);1::v4(6);1::v4(7);1::v4(8);1::v4(9);1::v4(10);1::v4(11);1::v4(12);1::v4(13);1::v4(14).
sequences :- v1(V1),v2(V2),v3(V3),v4(V4), V1 \= V2 ,V1 \= V3 ,V1 \= V4 ,V2 \= V3 ,V2 \= V4 ,V3 \= V4 .
query(sequences).

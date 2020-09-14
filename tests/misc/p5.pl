
students(8,[s1,s2,s3,s4,s5,s6,s7,s8]).
dutch(5,[s1,s2,s3,s4,s5,s6]).
french(4, [s3,s4,s5,s6,s7,s8]).
% students([1,13]).
% dutch([1,6]).
% french([4,11]).

structure(seq, sequence, true, students).
size(seq, 4).
% pos(seq,1,french).
% pos(seq,2,french).
% pos(seq,3,dutch).
pos(seq,3,dutch).
% pos(seq,3,not(french)).
% pos(seq,3,not(french)).
% pos(seq,8,not(dutch)).
% in(seq, s2).
% count(seq, french==2).
% count(seq, dutch<2).
count(seq, dutch>=2).
% pos(seq, 1, inter(dutch,french)).
% pos(seq, 2, not(dutch)).
% pos(seq,2,dutch).
% count(seq, inter(dutch,not(french))>1).
% count(seq, dutch>=2).
% count(seq, dutch>3).
% count(seq, dutch<5).
% pos(seq,1,s5) :- pos(seq,3,s1).

% query(seq).
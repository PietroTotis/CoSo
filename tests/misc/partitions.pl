
students(8,[s1,s2,s3,s4,s5,s6,s7,s8]).
% dutch(5,[s1,s2,s3,s4,s5,s6]).
% french(4, [s3,s4,s5,s6,s7,s8]).
% students([1,13]).
% dutch([1,6]).
% french([4,11]).

structure(part, partition, false, students).
size(part,==,4).
count(partition, size(partition,==,2) == 2).
count(partition, count(partition, dutch, >, 0), ==, 2).


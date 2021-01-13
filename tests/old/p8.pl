students([1,5]).
p1([1,3]).
p2([2,4]).

structure(part, partition, false, students).
size(part,==,3).
count(partition, count(partition, p1, >, 0), ==, 2).
% count(partition, count(partition, p1, <, 4), ==, 3).
% count(partition, count(partition, p2, <, 3), ==, 3).
% count(partition, count(partition, p2, >, 0), ==, 3).
count(partition, count(partition, p2, >, 0), ==, 1).
% count(partition, count(partition, p2, >, 0), ==, 3).
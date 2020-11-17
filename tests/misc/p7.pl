students([1,8]).
p1([1,4]).
p2([6,10]).

structure(part, partition, false, students).
size(part,==,3).
count(partition, size(partition,>=,2), ==, 3).
count(partition, size(partition,=<,3), ==, 3).
count(partition, count(partition, p1, >, 0), ==, 3).
count(partition, count(partition, p1, <, 4), ==, 3).
% count(partition, count(partition, p2, <, 3), ==, 3).
% count(partition, count(partition, p2, >, 0), ==, 3).

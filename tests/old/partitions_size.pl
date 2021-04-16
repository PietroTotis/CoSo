
% students = {s1,s2,s3,s4,s5,s6,s7,s8};
dutch = {s1,s2,s3,s4,s5,s6};
d = {s1,s2,s3};
% french = {s3,s4,s5,s6,s7,s8};
% dist = {s1,s2};
% indist rest = {s7, s8};
% students = [1,13].
% dutch = [1,6]}.
% french = [4,11]}.

parts in partitions(universe);
#parts = 2;
#{#p=3 | p in parts} = 2;
#{#d>0 | p in parts} = 2;
% #parts = 4;
% #{#p<=2 | p in parts} = 2;
% #{#dutch > 0 | p in parts} = 2;
% #{#french=0 | p in parts} > 0;


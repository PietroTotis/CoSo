include "globals.mzn";
set of int: universe = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
set of int: dom1 = {10,11,12,13,14,15};
set of int: dom2 = {11,12,13,14};
int: n = 8;
array[1..n] of var universe: sequence;

set of int: df1= universe diff dom2;
constraint sequence[7] in df1;

set of int: df2= universe diff dom2;
constraint sequence[8] in df2;
constraint sum(i in 1..n)(bool2int(sequence[i] in dom2)) <= 7;
set of int: df3= universe diff dom2;
set of int: df4= dom1 intersect df3;
constraint sum(i in 1..n)(bool2int(sequence[i] in df4)) <= 3;
constraint sum(i in 1..n)(bool2int(sequence[i] in dom2)) == 2;
solve satisfy;
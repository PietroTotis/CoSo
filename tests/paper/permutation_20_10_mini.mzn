include "globals.mzn";
set of int: uni = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
set of int: dom1 = {8,9,10,11,12,13};
set of int: dom2 = {4,5,6,7,8,9,10,11,12,13,14,15,16,17};
int: n = 10;
array[1..n] of var uni: sequence;
constraint alldifferent(sequence);
constraint sequence[3] in dom2;

set of int: df1= uni diff dom1;

set of int: df2= dom2 union df1;
constraint sequence[1] in df2;

set of int: df3= uni diff dom2;
constraint sequence[5] in df3;
constraint sum(i in 1..n)(bool2int(sequence[i] in dom2)) >= 0 /\ sum(i in 1..n)(bool2int(sequence[i] in dom2)) <= 4;
set of int: df4= uni diff dom2;
constraint sum(i in 1..n)(bool2int(sequence[i] in df4)) >= 0 /\ sum(i in 1..n)(bool2int(sequence[i] in df4)) <= 4;
solve satisfy;
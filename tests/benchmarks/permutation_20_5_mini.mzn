include "globals.mzn";
set of int: uni = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
set of int: dom1 = {8,9,10,11,12,13};
set of int: dom2 = {4,5,6,7,8,9,10,11,12,13,14,15,16,17};
int: n = 5;
array[1..n] of var uni: sequence;
constraint alldifferent(sequence);
constraint sequence[2] in dom2;
solve satisfy;
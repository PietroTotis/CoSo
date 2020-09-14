include "globals.mzn";
set of int: universe = {1,2,3,4,5,6,7,8,9,10};
set of int: dom1 = {8,9,10};
set of int: dom2 = {4,5,6,7};
int: n = 8;
array[1..n] of var universe: sequence;
constraint alldifferent(sequence);
solve satisfy;
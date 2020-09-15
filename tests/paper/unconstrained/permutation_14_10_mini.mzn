include "globals.mzn";
set of int: universe = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
set of int: dom1 = {5,6,7};
set of int: dom2 = {13,14};
int: n = 10;
array[1..n] of var universe: sequence;
constraint alldifferent(sequence);
solve satisfy;
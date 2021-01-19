include "globals.mzn";
set of int: uni = {1,2,3,4,5,6,7,8,9,10};
set of int: dom1 = {4,5,6,7};
set of int: dom2 = {2,3,4,5,6,7,8,9};
int: n = 5;
array[1..n] of var uni: sequence;
constraint alldifferent(sequence);
constraint sequence[2] in dom2;
solve satisfy;
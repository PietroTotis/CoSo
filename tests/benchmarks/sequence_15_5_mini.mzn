include "globals.mzn";
set of int: uni = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
set of int: dom1 = {4,5,6,7,8,9};
set of int: dom2 = {2,3,4,5,6,7,8,9,10,11,12,13,14};
int: n = 5;
array[1..n] of var uni: sequence;

set of int: df1= uni diff dom2 ([2,14]);
constraint sequence[4] in df1;
solve satisfy;
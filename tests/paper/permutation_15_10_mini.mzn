include "globals.mzn";
set of int: uni = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
set of int: dom1 = {4,5,6,7,8,9};
set of int: dom2 = {2,3,4,5,6,7,8,9,10,11,12,13,14};
int: n = 10;
array[1..n] of var uni: sequence;
constraint alldifferent(sequence);

set of int: df1= uni diff dom2;
constraint sequence[8] in df1;
constraint sequence[2] in dom2;

set of int: df2= dom1 intersect dom2;

set of int: df3= uni diff df2;
constraint sequence[9] in df3;
solve satisfy;
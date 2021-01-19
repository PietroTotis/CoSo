include "globals.mzn";
set of int: uni = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
set of int: dom1 = {8,9,10,11,12,13};
set of int: dom2 = {4,5,6,7,8,9,10,11,12,13,14,15,16,17};
int: n = 5;
var set of uni: sub;
constraint card(sub) == n;
constraint sum(i in sub)(bool2int(i in dom2)) >= 3 /\ sum(i in sub)(bool2int(i in dom2)) <= 6;
solve satisfy;
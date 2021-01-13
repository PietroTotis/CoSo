include "globals.mzn";       

set of int: student = {1,2,3,4,5,6,7,8};
set of int: dutch = {1,2,3,4,5};
set of int: french = {1,4,5,7,8};      
int: n = 3;    
var set of student: sets; 

% constraint sum(i in sets)(bool2int(i in dutch)) >= 3;
% constraint sum(i in sets)(bool2int(i in french)) == 2;
constraint 2 in sets;
constraint card(sets) == n;
solve satisfy;
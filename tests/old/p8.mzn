% Use this editor as a MiniZinc scratch book
set of int: dom = {1,2,3,4,5};
set of int: p1 = {1,2,3};
set of int: p2 = {2,3,4};

array[1..3] of var set of dom: parts;
% var set of dom: p1;
% var set of dom: p2;
% var set of dom: p3;
% var set of dom: p4;

% constraint (p1 union p2 union p3 union p4) = dom;
% constraint (p1 union p2) = dom;
constraint (parts[1] union parts[2] union parts[3]) = dom;
% constraint (p1 union p2) = dom;
constraint parts[1] intersect parts[2] = {};
constraint parts[1] intersect parts[3] = {};
% constraint p1 intersect p4 = {};
constraint parts[2] intersect parts[3] = {};
% constraint p2 intersect p4 = {};
% constraint p3 intersect p4 = {};

constraint parts[1] != {};
constraint parts[2] != {};
constraint parts[3] != {};

constraint(sum(p in 1..3)(bool2int(sum(i in parts[p])(bool2int(i in p1)) > 0)) == 2); 
constraint(sum(p in 1..3)(bool2int(sum(i in parts[p])(bool2int(i in p2)) > 0)) == 1); 

% constraint(sum(i in p1)(bool2int(i in s2)) > 0);
% constraint(sum(i in p2)(bool2int(i in s2)) > 0);
% constraint(sum(i in p3)(bool2int(i in s2)) > 0);
% constraint(sum(i in p1)(bool2int(i in s2)) < 3);
% constraint(sum(i in p2)(bool2int(i in s2)) < 3);
% constraint(sum(i in p3)(bool2int(i in s2)) < 3);

solve satisfy;
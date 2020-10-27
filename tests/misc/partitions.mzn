% Use this editor as a MiniZinc scratch book
set of int: dom = {1,2,3,4,5,6,7,8};

var set of dom: p1;
var set of dom: p2;
var set of dom: p3;
var set of dom: p4;

constraint (p1 union p2 union p3 union p4) = dom;
constraint p1 intersect p2 = {};
constraint p1 intersect p3 = {};
constraint p1 intersect p4 = {};
constraint p2 intersect p3 = {};
constraint p2 intersect p4 = {};
constraint p4 intersect p3 = {};

constraint p1 != {};
constraint p2 != {};
constraint p3 != {};
constraint p4 != {};

constraint p1 < p2 /\ p2 < p3 /\ p3 <p4;
constraint(sum(i in p1)(1) == 2);
constraint(sum(i in p2)(1) == 2);

solve satisfy;
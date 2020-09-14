% Use this editor as a MiniZinc scratch book
set of int: dom = {1,2,3,4,5,6,7,8};

var set of dom: p1;
var set of dom: p2;
var set of dom: p3;

constraint (p1 union p2 union p3) = dom;
constraint p1 intersect p2 = {};
constraint p1 intersect p3 = {};
constraint p2 intersect p2 = {};

constraint p1 != {};
constraint p2 != {};
constraint p3 != {};

solve satisfy;
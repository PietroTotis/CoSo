
set indist russians = {r1,r2};
set americans = {a1,a2,a3};
% spaniards = {s1};
selection in [| russians+americans];
% selection in {| russians+americans+spaniards};
#selection = 3;
#(russians & selection) > 0;
% #(spaniards & selection) > 0;
#(americans & selection) > 0;

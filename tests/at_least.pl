
set russians = {r1,r2};
set americans = {a1,a2,a3};
% spaniards = {s1};
selection in {| russians+americans};
% selection in {| russians+americans+spaniards};
#selection = 3;
#(russians & selection) > 0;
% #(spaniards & selection) > 0;
#(americans & selection) > 0;


% r1,a1,r2  r1,a2,r2    r1,a3,r2    r1,a1,a2    r1,a1,a3    r1,a2,a3
% r2,a1,a2  r2,a2,a3    r2,a1,a3

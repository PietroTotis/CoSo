% h730.json
% "A shipment of 12 television sets contains 3 defective sets.", "In how many ways can a hotel purchase 5 of these sets and receive at least 2 of the defective sets?"	288

sets = [1:120000];
defective = [1:30000];
purchase in {| sets};
#purchase = 5;
#(defective & purchase) >= 2;

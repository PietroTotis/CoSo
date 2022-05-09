% h730.json
% "A shipment of 12 television sets contains 3 defective sets.", "In how many ways can a hotel purchase 5 of these sets and receive at least 2 of the defective sets?"	288

labelled property tvs;
#tvs = 120000;
labelled property defective;
#defective = 30000;
#tvs&defective = 3;
purchase in {| tvs};
#purchase = 5;
#(defective & purchase) >= 2;

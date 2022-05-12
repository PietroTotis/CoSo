% m722.json
% "Fourteen construction workers are to be assigned to three different tasks.", "Seven workers are needed for mixing cement, five for laying bricks, and two for carrying the bricks to the brick layers.", "In how many different ways can the workers be assigned to these tasks?"	72072

labelled property workers;
#workers = 14;
assignments in compositions(workers);
#assignments=3;
#assignments[1] = 7;
#assignments[2] = 5;
#assignments[3] = 2;

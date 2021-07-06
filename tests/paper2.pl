squares = {r2,b2};
indist triangles = {r1,b1,g1,g2,g3};
indist red = {r1,r2};
indist blue = {b1,b2};
indist green = {g1,g2,g3};

perm in [|universe];
#perm = 4;
perm[2] = green;
#squares&perm =2;
squares = {s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12};
indist triangles = {t1,t2,t3,t4,t5,t6,t7};
indist red = {s1,s2,s3,t1,t2,t3,t4};
indist blue = {s4,s5,s6,s7};
indist green = {s8,s9,s10,s11,s12,t5,t6,t6};

perm in [|universe];
#perm = 4;
perm[2] = green;
#squares&perm =2;
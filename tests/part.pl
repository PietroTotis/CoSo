squares = {r2,b2};
indist triangles = {r1,b1,g1,g2,g3};
indist red = {r1,r2};
indist blue = {b1,b2};
indist green = {g1,g2,g3};

parts in partitions(red+green+blue);

#parts = 3;
#{ #part & green = 3 } =1;
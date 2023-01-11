property squares = {r2,b2};
property triangles = {r1,b1,g1,g2,g3};
property red = {r1,r2};
property blue = {b1,b2};
property green = {g1,g2,g3};

parts in {{red+green+blue}};

#parts = 3;
#( #part & green = 3 ) = 1;
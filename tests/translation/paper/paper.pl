universe shapes = {sq_r , sq_b ,tr_r ,tr_b ,tr_g ,tr_g ,tr_g };
property red = {tr_r , sq_r };
property blue = {tr_b , sq_b };
property green = {tr_g };
property triangle = {tr_r ,tr_b ,tr_g };
property squares = {sq_b , sq_r };

perm in [|universe];
#perm = 4;
perm[2] = green;
#squares&perm =2;
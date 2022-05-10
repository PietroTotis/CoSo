universe u = {e1, e1, e1, e1, e1, e1, e2, e2, e3, e4};
property dom1 = {e2, e3};
a in compositions(u);
#a <= 5;
#a[4] & dom1 != 1 ;
#{#dom1 = 1} < 2 ;

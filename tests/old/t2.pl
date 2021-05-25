stuff = [1:10];
f = [1:4];
seq in compositions(stuff);
#seq = 3;
#seq[1] = 3 ;
#{#p=3 | p in seq} = 2;
#{#p&f =2 | p in seq} = 1;
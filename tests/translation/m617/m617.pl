% m617.json
% "In how many distinguishable ways can the letters in B A N A N A be written?"	60

property bs = {b};
property as = {a,a,a};
property ns = {n,n};
words in [| bs+as+ns];
#words = 6;

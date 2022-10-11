shapes_0("5",3).
shapes(X) :- shapes_0(X, _).
shapes_1("1", 1).
shapes(X) :- shapes_1(X, _).
shapes_2("2", 1).
shapes(X) :- shapes_2(X, _).
shapes_3("3", 1).
shapes(X) :- shapes_3(X, _).
shapes_4("4", 1).
shapes(X) :- shapes_4(X, _).
universe(X) :- shapes(X).
red_0("1", 1).
red(X) :- red_0(X, _).
red_1("3", 1).
red(X) :- red_1(X, _).
universe(X) :- red(X).
blue_0("2", 1).
blue(X) :- blue_0(X, _).
blue_1("4", 1).
blue(X) :- blue_1(X, _).
universe(X) :- blue(X).
green_0("5",3).
green(X) :- green_0(X, _).
universe(X) :- green(X).
triangle_0("5",3).
triangle(X) :- triangle_0(X, _).
triangle_1("3", 1).
triangle(X) :- triangle_1(X, _).
triangle_2("4", 1).
triangle(X) :- triangle_2(X, _).
universe(X) :- triangle(X).
squares_0("1", 1).
squares(X) :- squares_0(X, _).
squares_1("2", 1).
squares(X) :- squares_1(X, _).
universe(X) :- squares(X).
df_0_0("1", 1).
df_0(X) :- df_0_0(X, _).
df_0_1("2", 1).
df_0(X) :- df_0_1(X, _).
universe(X) :- df_0(X).
permutation_guess_4(A,B,C,D) :- universe(A), pf_1(B), universe(C), universe(D).
pf_1_0("5",3).
pf_1(X) :- pf_1_0(X, _).
universe(X) :- pf_1(X).
1{permutation_4(A,B,C,D):permutation_guess_4(A,B,C,D)}1.
used_4(X,0) :- permutation_4(X, _, _, _). 
used_4(X,1) :- permutation_4(_, X, _, _). 
used_4(X,2) :- permutation_4(_, _, X, _). 
used_4(X,3) :- permutation_4(_, _, _, X). 
:- shapes_0(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- shapes_1(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- shapes_2(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- shapes_3(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- shapes_4(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- red_0(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- red_1(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- blue_0(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- blue_1(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- green_0(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- triangle_0(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- triangle_1(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- triangle_2(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- squares_0(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- squares_1(S,SN), C = #count{N:used_4(S,N)}, C>SN.
:- C = #count{N:used_4(S,N),df_0(S)}, C=0.
:- C = #count{N:used_4(S,N),df_0(S)}, C=1.
:- C = #count{N:used_4(S,N),df_0(S)}, C=3.
:- C = #count{N:used_4(S,N),df_0(S)}, C=4.

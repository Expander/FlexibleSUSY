Block MODSEL                 # Select model
    1    0                   # mSUGRA
    3    1                   # NMSSM
    6   0                    # flavour violation
   12   1500                 #

Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   0                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   2                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                   # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)

Block SOFTSUSY               # SOFTSUSY specific inputs
    1   1.000000000e-04      # tolerance
    2   2.000000000e+00      # up-quark mixing (=1) or down (=2)
    5   1.000000000E+00      # 2-loop running
    3   1.000000000E+00      # printout
        4   1500                #
    7   2.0             
   18   1.000000000E+00      # use soft Higgs masses as EWSB output

BLOCK SMINPUTS	 	 # 
1	128.962	 	 # 
2	0.000011663900000000002	 	 # 
3	0.1172	 	 # 
4	91.1876	 	 # 
5	4.2	 	 # 
6	172.9	 	 # 
7	1.777	 	 # 
9	80.385	 	 # 
11	0.00051099891	 	 # 
13	1.10565836	 	 # 
21	0.00495	 	 # 
22	0.0025	 	 # 
23	0.1	 	 # 
24	1.42	 	 # 

BLOCK EXTPAR	 	 # 
0	1500.	 	 # 
1	500.	 	 # 
2	1000.	 	 # 
3	3000.	 	 # 
11	0.	 	 # At
12	0.	 	 # 
13	0.	 	 # 
25	10.0	 	 # 
31	1500.1	 	 # 
32	1500.2	 	 # 
33	1500.3	 	 # 
34	1500.4	 	 # 
35	1500.5	 	 # 
36	1500.6	 	 # 
41	1500.7	 	 # 
42	1500.8	 	 # 
43	1500.	 	 # 
44	1500.01	 	 # 
45	1500.02	 	 # 
46	1500.0	 	 # 
47	1500.04	 	 # 
48	1500.05	 	 # 
49	1500.06	 	 # 
61	0.1	 	 # 
62	0.1	 	 # 
63	-10.	 	 # ALambda
64	-10.0	 	 # AKappa
65	900.	 	 # 

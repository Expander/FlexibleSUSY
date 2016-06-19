Needs["TestSuite`", "TestSuite.m"];
Needs["ReadSLHA`", "ReadSLHA.m"];

slha = "
Block MODSEL                 # Select model
#   12    1000                # parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   1                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top quark 2-loop corrections QCD
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   17   0                    # calculate parametric uncertainties
   17   1                    # calculate parametric uncertainties
Block FlexibleSUSYInput
    0   0.00729735           # alpha_em(0)
    1   125.09               # Mh pole
Block SMINPUTS               # Standard Model inputs
    1   1.279340000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166370000e-05      # G_Fermi
    3   1.176000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.200000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.404               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
Block MINPAR                 # Input parameters
    1   0.192                # Lambda(Qin)
Block EXTPAR                 # Input parameters
    0   1000                 # input scale Qin
    1   173.34               # scale QEWSB
Block EXTPAR
    0   1.1
    1   2.2
Block V
    1   10.1
    2   20.2
    3   a
Block M
    1 1  10.1
    1 2  1
    2 2  20.2
    3 3  30.3
Block T3
    1 1 1  10.1
    2 2 2  20.2
    3 3 3  30.3
Block T4
    1 1 1 1  10.1
    2 2 2 2  20.2
    3 3 3 3  30.3
";

pars = {
    {Qin  , {0}            , {EXTPAR, 0}},
    {QEWSB, {0}            , {EXTPAR, 1}},
    {V    , {3}            , V},
    {M    , {3, 3}         , M},
    {T3   , {3, 3, 3}      , T3},
    {T4   , {3, 3, 3, 3}   , T4}
};

values = ReadSLHAString[slha, pars];

TestEquality[Qin /. values, 1.1];
TestEquality[QEWSB /. values, 2.2];
TestEquality[V /. values, {10.1, 20.2, 0}];
TestEquality[M /. values, {{10.1, 1, 0}, {0, 20.2, 0}, {0, 0, 30.3}}];
TestEquality[T3 /. values, {{{10.1, 0, 0}, {0, 0, 0}, {0, 0, 0}},
                            {{0, 0, 0}, {0, 20.2, 0}, {0, 0, 0}},
                            {{0, 0, 0}, {0, 0, 0}, {0, 0, 30.3}}}];

t4 = Table[0, {i,1,3}, {j,1,3}, {k,1,3}, {l,1,3}];
t4[[1,1,1,1]] = 10.1;
t4[[2,2,2,2]] = 20.2;
t4[[3,3,3,3]] = 30.3;

TestEquality[T4 /. values, t4];

Print[""];

PrintTestSummary[];

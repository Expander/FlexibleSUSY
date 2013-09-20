Needs["TestSuite`", "TestSuite.m"];
Needs["EWSB`", "EWSB.m"];

Print["testing FindFreePhasesInEWSB[] ..."];

mssmEwsbEqs = {
    {mu^2 + x^2 + x y + z + 5},
    {Bmu  - x^2 + x y + z + 5}
};

mssmEwsbOutputParameters = { mu, Bmu };

Parameters`AddRealParameter[mssmEwsbOutputParameters];

TestEquality[EWSB`FindFreePhasesInEWSB[mssmEwsbEqs, mssmEwsbOutputParameters],
             {FlexibleSUSY`Sign[mu]}];

PrintTestSummary[];

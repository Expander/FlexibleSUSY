Needs["TestSuite`", "TestSuite.m"];
Needs["EWSB`", "EWSB.m"];

Print["testing CheckEWSBEquations[] ..."];

mssmEwsbEqs = {
    {mu^2 + x^2 + x y + z + 5},
    {Bmu  - x^2 + x y + z + 5}
};

mssmEwsbOutputParameters = { mu, Bmu };

Parameters`AddRealParameter[mssmEwsbOutputParameters];

TestEquality[EWSB`CheckEWSBEquations[mssmEwsbEqs, mssmEwsbOutputParameters],
             {FlexibleSUSY`Sign[mu]}];

PrintTestSummary[];

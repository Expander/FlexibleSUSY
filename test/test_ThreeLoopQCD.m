Needs["SARAH`"];
Needs["TestSuite`", "TestSuite.m"];
Needs["TwoLoop`", "TwoLoop.m"];
Needs["ThreeLoopQCD`", "ThreeLoopQCD.m"];

Start["SM"];

Print["Testing 1L arxiv:hep-ph/9912391 vs. hep-ph/9803493 ..."];

m1 = M/(1 + 
    h GetDeltaMOverMQCDOneLoopMSbar[TopQuark, FlexibleSUSY`M[TopQuark]]);

m1 = Simplify[Normal[Series[m1, {h, 0, 1}]] /. h -> 1];

m2 = M (GetMTopMSbarOverMTopPole[{1, 0, 0, 0}] +
        GetMTopMSbarOverMTopPole[{0, 1, 0, 0}])

TestEquality[Simplify[m1 - m2], 0];

Print["Testing 2L arxiv:hep-ph/9912391 vs. hep-ph/9803493 ..."];

m1 = M/(1 + 
    h GetDeltaMOverMQCDOneLoopMSbar[TopQuark, FlexibleSUSY`M[TopQuark]] + 
    h^2 GetDeltaMOverMQCDTwoLoopMSbar[TopQuark, FlexibleSUSY`M[TopQuark]]);

m1 = Simplify[Normal[Series[m1, {h, 0, 2}]] /. h -> 1];

m2 = M (GetMTopMSbarOverMTopPole[{1, 0, 0, 0}] +
        GetMTopMSbarOverMTopPole[{0, 1, 0, 0}] +
        GetMTopMSbarOverMTopPole[{0, 0, 1, 0}])

TestEquality[Simplify[m1 - m2], 0];

PrintTestSummary[];

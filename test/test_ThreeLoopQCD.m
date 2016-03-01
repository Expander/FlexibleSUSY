Needs["SARAH`"];
Needs["TestSuite`", "TestSuite.m"];
Needs["TwoLoopQCD`", "TwoLoopQCD.m"];
Needs["ThreeLoopQCD`", "ThreeLoopQCD.m"];

Start["SM"];

Print["Testing 1L arxiv:hep-ph/9912391 vs. arxiv:hep-ph/9803493 ..."];

m1 = M/(1 + 
    h GetDeltaMOverMQCDOneLoopMSbar[TopQuark, Q]);

m1 = Simplify[Normal[Series[m1, {h, 0, 1}]] /. h -> 1];

m2 = M (GetMTopMSbarOverMTopPole[{1, 0, 0, 0}] +
        GetMTopMSbarOverMTopPole[{0, 1, 0, 0}])

TestEquality[Simplify[m1 - m2], 0];

Print["Testing 2L arxiv:hep-ph/9912391 vs. arxiv:hep-ph/9803493 ..."];

m1 = M/(1 + 
    h GetDeltaMOverMQCDOneLoopMSbar[TopQuark, Q] + 
    h^2 GetDeltaMOverMQCDTwoLoopMSbar[TopQuark, Q]);

m1 = Simplify[Normal[Series[m1, {h, 0, 2}]] /. h -> 1];

m2 = M (GetMTopMSbarOverMTopPole[{1, 0, 0, 0}] +
        GetMTopMSbarOverMTopPole[{0, 1, 0, 0}] +
        GetMTopMSbarOverMTopPole[{0, 0, 1, 0}])

TestEquality[Simplify[m1 - m2], 0];

Print["Testing 2L arxiv:hep-ph/9803493 vs. arxiv:hep-ph/9912391 ..."];

M1 = m (1 +
        h   GetDeltaMOverMQCDOneLoopMSbar[TopQuark, Q] +
        h^2 GetDeltaMOverMQCDTwoLoopMSbar[TopQuark, Q]);

M2 = m / (GetMTopMSbarOverMTopPole[{1, 0, 0  , 0}] +
          GetMTopMSbarOverMTopPole[{0, h, 0  , 0}] +
          GetMTopMSbarOverMTopPole[{0, 0, h^2, 0}])

M2 = Simplify[Normal[Series[M2, {h, 0, 2}]]];

TestEquality[Simplify[(M1 - M2) /. h -> 1], 0];

PrintTestSummary[];

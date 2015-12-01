Needs["TestSuite`", "TestSuite.m"];
Needs["ThreeLoopMSSM`", "ThreeLoopMSSM.m"];

FlexibleSUSY`$flexiblesusyMetaDir = FileNameJoin[{Directory[], "meta"}];

Print["testing BetaMSSM[] ..."];

TestEquality[Length[ThreeLoopMSSM`BetaMSSM[SARAH`hyperchargeCoupling]], 3];
TestEquality[Length[ThreeLoopMSSM`BetaMSSM[SARAH`leftCoupling       ]], 3];
TestEquality[Length[ThreeLoopMSSM`BetaMSSM[SARAH`strongCoupling     ]], 3];
TestEquality[Length[ThreeLoopMSSM`BetaMSSM[SARAH`UpYukawa           ]], 3];
TestEquality[Length[ThreeLoopMSSM`BetaMSSM[SARAH`DownYukawa         ]], 3];
TestEquality[Length[ThreeLoopMSSM`BetaMSSM[SARAH`ElectronYukawa     ]], 3];

Print["checking BetaMSSM[] and MSSM[] ..."];

SARAH`SARAH[OutputDirectory] = FileNameJoin[{Directory[], "Output"}];
SARAH`Start["MSSM"];
SARAH`CalcRGEs[ReadLists -> True, TwoLoop -> True,
               NoMatrixMultiplication -> False];

FindBetaFunction[lst_List, c_] :=
    Cases[lst, {c | c[__], beta__} :> {beta}][[1]];

UniformTraces[] := {
    trace[Adj[x_], y_] :> trace[y, Adj[x]],
    trace[a_, Adj[b_], x_, Adj[d_]] :> trace[Adj[b], x, Adj[d], a]
};

TestBetaEquality[lst_, c_, loop_] :=
    Module[{sa, sm},
           sm = ThreeLoopMSSM`BetaMSSM[c][[loop]] /. UniformTraces[];
           sa = FindBetaFunction[lst, c][[loop]] /. a_[i1,i2] :> a /. UniformTraces[];
           TestEquality[Expand[sm], Expand[sa]]
          ];

For[l = 1, l <= 2, l++,
    TestBetaEquality[SARAH`BetaGauge, SARAH`hyperchargeCoupling, l];
    TestBetaEquality[SARAH`BetaGauge, SARAH`leftCoupling       , l];
    TestBetaEquality[SARAH`BetaGauge, SARAH`strongCoupling     , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`UpYukawa           , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`DownYukawa         , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`ElectronYukawa     , l];
   ];

PrintTestSummary[];

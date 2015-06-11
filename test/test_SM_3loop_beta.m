Needs["TestSuite`", "TestSuite.m"];
Needs["ThreeLoopSM`", "ThreeLoopSM.m"];

$flexiblesusyMetaDir = FileNameJoin[{Directory[], "meta"}];

Print["testing BetaSM[] ..."];

TestEquality[Length[ThreeLoopSM`BetaSM[SARAH`hyperchargeCoupling]], 3];
TestEquality[Length[ThreeLoopSM`BetaSM[SARAH`leftCoupling       ]], 3];
TestEquality[Length[ThreeLoopSM`BetaSM[SARAH`strongCoupling     ]], 3];
TestEquality[Length[ThreeLoopSM`BetaSM[SARAH`UpYukawa           ]], 3];
TestEquality[Length[ThreeLoopSM`BetaSM[SARAH`DownYukawa         ]], 3];
TestEquality[Length[ThreeLoopSM`BetaSM[SARAH`ElectronYukawa     ]], 3];
TestEquality[Length[ThreeLoopSM`BetaSM[\[Lambda]                ]], 3];

Print["checking BetaSM[] and SM[] ..."];

SARAH`SARAH[OutputDirectory] = FileNameJoin[{Directory[], "Output"}];
SARAH`Start["SM"];
SARAH`CalcRGEs[ReadLists -> True, TwoLoop -> True,
               NoMatrixMultiplication -> False];

FindBetaFunction[lst_List, c_] :=
    Cases[lst, {c | c[__], beta__} :> {beta}][[1]];

TestBetaEquality[lst_, c_, loop_] :=
    Module[{sa, sm},
           sm = ThreeLoopSM`BetaSM[c][[loop]];
           sa = FindBetaFunction[lst, c][[loop]] //. {
               SARAH`UpYukawa       -> gt,
               SARAH`DownYukawa     -> gb,
               SARAH`ElectronYukawa -> g\[Tau],
               SARAH`Adj[a_]        :> a,
               a_[Susyno`LieGroups`i1,SARAH`i2] :> a
           } /. {
               SARAH`MatMul[a__]    :> Times[a],
               SARAH`trace[a__]     :> Times[a],
               SARAH`trace[a_]      :> a
           };
           TestEquality[Expand[sm], Expand[sa]]
          ];

For[l = 1, l <= 2, l++,
    TestBetaEquality[SARAH`BetaGauge, SARAH`hyperchargeCoupling, l];
    TestBetaEquality[SARAH`BetaGauge, SARAH`leftCoupling       , l];
    TestBetaEquality[SARAH`BetaGauge, SARAH`strongCoupling     , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`UpYukawa           , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`DownYukawa         , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`ElectronYukawa     , l];
    TestBetaEquality[SARAH`BetaLijkl, \[Lambda]                , l];
   ];

PrintTestSummary[];

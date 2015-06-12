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

TestBetaEquality[lst_, c_, csm_, loop_] :=
    Module[{sa, sm},
           sm = ThreeLoopSM`BetaSM[csm][[loop]];
           sa = FindBetaFunction[lst, c][[loop]] /. {
               SARAH`MatMul[a__][_,_] :> Times[a],
               SARAH`MatMul[a__]    :> Times[a],
               SARAH`trace[a__]     :> Times[a],
               SARAH`trace[a_]      :> a
           } //. {
               SARAH`Adj[a_]        :> a,
               SARAH`UpYukawa       -> gt,
               SARAH`DownYukawa     -> gb,
               SARAH`ElectronYukawa -> g\[Tau],
               a_[Susyno`LieGroups`i1,SARAH`i2] :> a
           } /. {
               gt -> SARAH`UpYukawa[3,3],
               gb -> SARAH`DownYukawa[3,3],
               g\[Tau] -> SARAH`ElectronYukawa[3,3]
           };
           TestEquality[Expand[sm], Expand[sa]]
          ];

For[l = 1, l <= 2, l++,
    TestBetaEquality[SARAH`BetaGauge, SARAH`hyperchargeCoupling, SARAH`hyperchargeCoupling, l];
    TestBetaEquality[SARAH`BetaGauge, SARAH`leftCoupling       , SARAH`leftCoupling       , l];
    TestBetaEquality[SARAH`BetaGauge, SARAH`strongCoupling     , SARAH`strongCoupling     , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`UpYukawa           , SARAH`UpYukawa           , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`DownYukawa         , SARAH`DownYukawa         , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`ElectronYukawa     , SARAH`ElectronYukawa     , l];
    TestBetaEquality[SARAH`BetaLijkl, \[Lambda]                , \[Lambda]                , l];
    TestBetaEquality[SARAH`BetaBij  , mu2                      , m2                       , l];
   ];

PrintTestSummary[];

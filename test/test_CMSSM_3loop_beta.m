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
    trace[a_, Adj[b_], x_, Adj[d_]] :> trace[Adj[b], x, Adj[d], a],
    trace[args__ /; !FreeQ[{args}, SARAH`Tp[_]]] :> trace[Sequence @@ (Reverse[SARAH`Tp /@ {args}])],
    trace[Adj[Yd], T[Yd], Adj[Yu], Yu] -> trace[Adj[Yu], Yu, Adj[Yd], T[Yd]],
    trace[md2, Yd, Adj[Yd]] -> trace[Yd, Adj[Yd], md2],
    trace[me2, Ye, Adj[Ye]] -> trace[Ye, Adj[Ye], me2],
    trace[mq2, Adj[Yd], Yd] -> trace[Adj[Yd], Yd, mq2],
    trace[mq2, Adj[Yu], Yu] -> trace[Adj[Yu], Yu, mq2],
    trace[mu2, Yu, Adj[Yu]] -> trace[Yu, Adj[Yu], mu2],
    trace[ml2, Adj[Ye], Ye] -> trace[Adj[Ye], Ye, ml2],
    trace[Adj[Yd], Yd, Adj[Yd], T[Yd]] -> trace[Adj[Yd], T[Yd], Adj[Yd], Yd],
    trace[Adj[Ye], Ye, Adj[Ye], T[Ye]] -> trace[Adj[Ye], T[Ye], Adj[Ye], Ye],
    trace[Adj[Yu], Yu, Adj[Yu], T[Yu]] -> trace[Adj[Yu], T[Yu], Adj[Yu], Yu]
};

CalcDifference[a_, b_] := Simplify[Expand[a] - Expand[b]];

TestBetaEquality[lst_, c_, loop_] :=
    Module[{sa, sm},
           sm = ThreeLoopMSSM`BetaMSSM[c][[loop]] /. UniformTraces[];
           sa = FindBetaFunction[lst, c][[loop]] /. {
               Kronecker[_,_] :> 1,
               a_[i1,i2] :> a,
               conj[a_] :> a,
               Conj[a_] :> a,
               Tr1[1] -> 0, (* why is this term missing in the JJ result? *)
               Tr2[2] -> (mHd2 + mHu2 + trace[ml2] + 3*trace[mq2])/2,
               Tr2[3] -> (trace[md2] + 2*trace[mq2] + trace[mu2])/2,
               Tr3[1] -> 0, (* why is this term missing in the JJ result? *)
               Tr2U1[1,1] -> (g1^2*(3*mHd2 + 3*mHu2 + 2*trace[md2] + 6*trace[me2] +
                                    3*trace[ml2] + trace[mq2] + 8*trace[mu2]))/10
           } /. UniformTraces[] /. SARAH`Tp -> SARAH`Adj /. UniformTraces[];
           TestEquality[CalcDifference[sm, sa], 0];
          ];

For[l = 1, l <= 2, l++,
    TestBetaEquality[SARAH`BetaGauge, SARAH`hyperchargeCoupling, l];
    TestBetaEquality[SARAH`BetaGauge, SARAH`leftCoupling       , l];
    TestBetaEquality[SARAH`BetaGauge, SARAH`strongCoupling     , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`UpYukawa           , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`DownYukawa         , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`ElectronYukawa     , l];
    TestBetaEquality[SARAH`BetaMi   , MassB                    , l];
    TestBetaEquality[SARAH`BetaMi   , MassWB                   , l];
    TestBetaEquality[SARAH`BetaMi   , MassG                    , l];
    TestBetaEquality[SARAH`BetaMuij , \[Mu]                    , l];
    TestBetaEquality[SARAH`BetaBij  , B[\[Mu]]                 , l];
    TestBetaEquality[SARAH`BetaTijk , SARAH`TrilinearUp        , l];
    TestBetaEquality[SARAH`BetaTijk , SARAH`TrilinearDown      , l];
    TestBetaEquality[SARAH`BetaTijk , SARAH`TrilinearLepton    , l];
    TestBetaEquality[SARAH`Betam2ij , SARAH`SoftSquark         , l];
    TestBetaEquality[SARAH`Betam2ij , SARAH`SoftUp             , l];
    TestBetaEquality[SARAH`Betam2ij , SARAH`SoftDown           , l];
    TestBetaEquality[SARAH`Betam2ij , SARAH`SoftLeftLepton     , l];
    TestBetaEquality[SARAH`Betam2ij , SARAH`SoftRightLepton    , l];
   ];

PrintTestSummary[];

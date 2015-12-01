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

Roatations[trace[lst__]] :=
    Table[trace @@ RotateLeft[{lst},i], {i,1,Length[{lst}]}];

(* calc difference, taking equivalent trace expressions into account *)
CalcDifference[a_, b_] :=
    Module[{diff, traces, differences = {}, repls = {}, i},
           diff = Simplify[Expand[a] - Expand[b]];
           traces = Cases[diff, trace[__], {0, Infinity}];
           For[i = 1, i <= Length[traces], i++,
               AppendTo[repls, (traces[[i]] -> #)& /@ Roatations[traces[[i]]]];
              ];
           repls = Outer[List, repls];
           differences = (Sequence @@ #)& /@ ((Simplify[Expand[a /. #] - Expand[b]])& /@ repls);
           If[MemberQ[differences, 0], 0, diff]
          ];

TestBetaEquality[lst_, c_, loop_] :=
    Module[{sa, sm, difference},
           sm = ThreeLoopMSSM`BetaMSSM[c][[loop]] /. UniformTraces[];
           sa = FindBetaFunction[lst, c][[loop]] /. {
               Kronecker[_,_] :> 1,
               a_[i1,i2] :> a,
               conj[a_] :> a,
               Tr1[1] -> 0 (* why is this term missing in the JJ result? *)
           } /. UniformTraces[];
           difference = CalcDifference[sm, sa];
           TestEquality[difference, 0];
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

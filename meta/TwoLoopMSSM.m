BeginPackage["TwoLoopMSSM`"];
EndPackage[];

GetMSSMCPEvenHiggsLoopMassMatrix::usage = "Returns the
 loop-corrections to the CP-even Higgs mass matrix in the
 CP-conserving MSSM.  The loop-corrections are taken from
 arxiv:hep-ph/0105096.

Note: The return value contains the contributions from tadpole
 diagrams.

Note: The sign of the mu parameter is opposite to the one in SARAH.
 In order to switch to the SARAH convention, set $signMu = -1;

Usage: GetMSSMCPEvenHiggsLoopMassMatrix[loopOrder -> {1,1,1}, corrections -> {1}]

Parameters:

- loopOrder: List of factors multiplied by each loop order.
  #1: tree-level
  #2: 1-loop level
  #3: 2-loop level
  (default: {1,1,1})

- corrections: List of factors multiplied by each correction.
  (default: {1})
  #1: alpha_t * alpha_s
";

ReplaceStopMasses::usage = "Returns list of replacetment rules which
 replace mst1, mst2 and sin2Theta by DR-bar parameters.";

(* DR-bar parameters *)
{ ht, Mu, mt, At, mQ33, mU33, TanBeta, g3, Q, M3 };

(* Stop mass parameters *)
{ sin2Theta, mst1, mst2 };

$signMu::usage = "Sign of the mu parameter (default: 1).
In order to switch to the SARAH convention, set $signMu = -1";

$signMu = 1;

Begin["TwoLoopMSSM`Private`"];

(* Eqs. (17) of arxiv:hep-ph/0105096 *)
CalculateMStop2[signMu_] :=
    Module[{mst2 = {0,0}, mL2, mR2, Xt, s2t},
           mL2 = mQ33^2 + mt^2;
           mR2 = mU33^2 + mt^2;
           Xt = mt Abs[At + signMu Mu / TanBeta]; (* Eq. (14) and below *)
           (* Eq. (17) *)
           mst2[[1]] = 1/2 ((mL2 + mR2) + Sqrt[(mL2 - mR2)^2 + 4 Abs[Xt]^2]);
           mst2[[2]] = 1/2 ((mL2 + mR2) - Sqrt[(mL2 - mR2)^2 + 4 Abs[Xt]^2]);
           s2t = 2 Abs[Xt] / (mst2[[1]] - mst2[[2]]); (* Eq. (18) *)
           {mst2[[1]], mst2[[2]], s2t}
          ];

(* Eqs. (25)-(27) of arxiv:hep-ph/0105096 *)
CreateMassMatrixCPEven[F1_, F2_, F3_, DeltaF2_, DeltaF3_] :=
    Module[{mm = {{0,0},{0,0}}},
           mm[[1,1]] = 1/2 ht^2 Mu^2 sin2Theta^2 F3;
           mm[[1,2]] = (ht^2 $signMu Mu mt sin2Theta F2
                        + 1/2 ht^2 At $signMu Mu sin2Theta^2 (F3 + DeltaF3));
           mm[[2,1]] = mm[[1,2]];
           mm[[2,2]] = (2 ht^2 mt^2 F1 + 2 ht^2 At mt sin2Theta (F2 + DeltaF2)
                        + 1/2 ht^2 At^2 sin2Theta^2 (F3 + 2 DeltaF3));
           mm
          ];

(* Eq. (31) of arxiv:hep-ph/0105096 *)
GetMSSMCPEvenHiggsLoopMassMatrix1LAlphaTAlphaS[] :=
    Module[{F1, F2, F3, Nc = 3, h = 1/(16 Pi^2)},
           F1 = Nc h Log[mst1^2 mst2^2 / mt^4];
           F2 = Nc h Log[mst1^2 / mst2^2];
           F3 = Nc h (2 - (mst1^2 + mst2^2)/(mst1^2 - mst2^2) Log[mst1^2 / mst2^2]);
           CreateMassMatrixCPEven[F1, F2, F3, 0, 0]
          ];

(* Eq. (A4) of arxiv:hep-ph/0105096 *)
Phi[x_, y_, z_] :=
    Module[{u = x/z, v = y/z, lambda, xp, xm},
           lambda = Sqrt[(1 - u - v)^2 - 4 u v];
           xp = 1/2 (1 + (u - v) - lambda);
           xm = 1/2 (1 - (u - v) - lambda);
           1/lambda (2 Log[xp] Log[xm] - Log[u] Log[v]
                     - 2 (PolyLog[2,xp] + PolyLog[2,xm]) + Pi^2/3)
          ];

(* Eq. (A1) of arxiv:hep-ph/0105096 *)
f1[mt_, mg_, msqu_, msqd_, s2t_, Q_] :=
    Module[{delta},
           delta = mg^4 + mt^4 + msqu^2 - 2 (mg^2 mt^2 + mg^2 msqu + mt^2 msqu);
           (
               4 (mt^2 + mg^2 - mg mt s2t)/(msqu) (1 - Log[mg^2/Q^2])
               + 4 Log[mt^2/mg^2]
               - 2 Log[msqu/mg^2]
               + 2 / delta (
                   + 4 mg^4 Log[msqu/mg^2]
                   + (mg^4 - msqu^2 +mt^2 (10 mg^2 + 3 mt^2 + 2 (mt^2 mg^2 - mt^4)/msqu)) Log[mt^2/mg^2]
                 )
               + (2 mg s2t)/mt (Log[msqu/Q^2]^2 + 2 Log[mt^2/Q^2] Log[msqu/Q^2])
               + (4 mg s2t)/(mt delta) (
                   + mg^2 (msqu - mt^2 - mg^2) Log[msqu/mg^2]
                   + mt^2 (msqu - 3 mg^2 - 2 mt^2 - (mt^2 mg^2 - mt^4)/msqu) Log[mt^2/mg^2]
                 )
               + (4 mg^2 (mt^2 + mg^2 - msqu - 2 mg mt s2t) / delta
                  - (4 mg s2t)/mt) Phi[mt^2, msqu, mg^2]
           )
          ];

(* Eq. (A2) of arxiv:hep-ph/0105096 *)
f2[mt_, mg_, msqu_, msqd_, s2t_, Q_] :=
    Module[{delta, diff},
           delta = mg^4 + mt^4 + msqu^2 - 2 (mg^2 mt^2 + mg^2 msqu + mt^2 msqu);
           diff = msqu - msqd;
           (
               4 (mt^2 + mg^2)/msqu
               - (4 mg s2t)/(mt diff) (
                   + 3 msqu - (mt^2 msqd)/(msqu)
                 )
               + (2 mg s2t)/(mt diff) (
                   + (4 mt^2 + 5 msqu + msqd) Log[msqu/Q^2]
                   - 2 (mt^2 msqd)/msqu Log[mg^2/Q^2]
                 )
               - 4 (mg^2 + mt^2)/msqu Log[mg^2/Q^2] 
               - 2 Log[msqu/mg^2]
               + 2/delta (
                   2 mg^2 (mg^2 + mt^2 - msqu) Log[msqu/mg^2]
                   + 2 mt^2 (
                       + 3 mg^2 + 2 mt^2 - msqu 
                       + (mg^2 mt^2-mt^4)/(msqu)
                     ) Log[mt^2/mg^2]
                 )
               - (4 mg mt s2t)/(msqu delta) (
                   2 msqu mg^2 Log[msqu/mg^2]
                   - ((mt^2 - msqu)^2 - mg^2 (mt^2 + msqu)) Log[mt^2/mg^2]
                 )
               - (8 mg mt)/(s2t diff) (
                   Log[msqu/Q^2] - Log[mt^2/Q^2] Log[msqu/Q^2]
                 )
               - (mg s2t)/(mt diff) (
                   (msqu + msqd) Log[msqu/Q^2]^2
                   + (10 mt^2 - 2 mg^2 + msqu + msqd) Log[mt^2/Q^2] Log[msqu/Q^2] 
                   + (2 mg^2 - 2 mt^2 + msqu + msqd) Log[msqu/Q^2] Log[mg^2/Q^2]
                 )
               + ((8 mg^2 mt^2)/delta
                  - (8 mg mt)/(s2t diff)
                  + (2 s2t (4 mt^2  mg^2 - delta))/(mg mt diff)
                  + (s2t (msqu - mg^2 - mt^2)^3)/(mg mt delta)
                 ) Phi[mt^2, msqu, mg^2]
           )
          ];

(* Eq. (A3) of arxiv:hep-ph/0105096 *)
f3[mt_, mg_, msqu_, msqd_, s2t_, Q_] :=
    Module[{delta, diff},
           delta = mg^4 + mt^4 + msqu^2 - 2 (mg^2 mt^2 + mg^2 msqu + mt^2 msqu);
           diff = msqu - msqd;
           (
               - 4 (msqd (mg^2 + mt^2))/(msqu diff)
               + (4 mg mt s2t)/(diff^2) (
                   21 msqu - msqd^2/msqu
                 )
               + 4/diff ( 
                   (mg^2msqd)/(msqu)Log[mg^2/Q^2]
                   - 2 (mt^2 + mg^2) Log[msqu/Q^2]
                 )
               - (24 mg mt s2t (3 msqu + msqd))/(diff^2)Log[msqu/Q^2]
               + (4 mt^2)/(msqu delta) (
                   2 mg^2 msqu Log[msqu/Q^2] 
                   -mg^2 (mg^2 - mt^2 +msqu) Log[mg^2/Q^2]  
                   - ((mt^2 - msqu)^2-mg^2 (mt^2+msqu)) Log[mt^2/Q^2] 
                 )
               - (4 mg mt s2t)/(msqu delta)(
                   mt^2 (mg^2 -mt^2+msqu) Log[mt^2/Q^2] 
                   - mg^2 (mg^2 -mt^2-msqu) Log[mg^2/Q^2] 
                   + msqu (mg^2 +mt^2-msqu) Log[msqu/Q^2]
                 )
               + 2 (2 mg^2 + 2 mt^2 - msqu - msqd)/diff Log[mt^2 mg^2/Q^4]Log[msqu/Q^2]
               + (12 mg mt s2t)/(diff^2)(
                   2 (mg^2 - mt^2) Log[mg^2/mt^2] Log[msqu/Q^2]
                   + (msqu + msqd) Log[mt^2 mg^2/Q^4] Log[msqu/Q^2]
                 )
               + (8 mg mt)/(s2t diff^2) (
                   - 8 msqu + 2 (3 msqu + msqd)Log[msqu/Q^2] 
                   - 2 (mg^2 - mt^2) Log[mg^2/mt^2] Log[msqu/Q^2]
                   - (msqu + msqd) Log[mt^2 mg^2/Q^4] Log[msqu/Q^2]
                 )
               - (
                   (8/s2t - 12 s2t) (mt (2 delta + (mg^2 + mt^2 - msqu) diff))/(mg diff^2)
                   + (4 delta + 8 mg^2 mt^2)/(mg^2 diff)
                   + (2 (mg^2 + mt^2 - msqu))/(mg^2)
                   - (4 mt^2 (mg^2 + mt^2 - msqu - 2 mg mt s2t))/delta
                 ) Phi[mt^2, msqu, mg^2]
           )
          ];

(* Eqs. (32)-(34) of arxiv:hep-ph/0105096 *)
GetMSSMCPEvenHiggsLoopMassMatrix2LAlphaTAlphaS[] :=
    Module[{unit, CF = 4/3, Nc = 3, h = 1/(16 Pi^2),
            mg, cos2Theta2, F1, F2, F3, F31L, DeltaF2, DeltaF3},
           unit = g3^2 CF Nc h^2;
           cos2Theta2 = 1 - sin2Theta^2;
           mg = Abs[M3];
           (* Eq. (32) *)
           F1 = unit (
               - 6 (1 - Log[mt^2 / Q^2])
               + 5 Log[mst1^2 mst2^2 / mt^4]
               + Log[mst1^2 mst2^2 / mt^4]^2
               + 8 Log[mt^2 / Q^2]^2
               - 4 (Log[mst1^2 / Q^2]^2 + Log[mst2^2 / Q^2]^2)
               - cos2Theta2 (2 - Log[mst1^2 mst2^2 / Q^4] - Log[mst1^2 / mst2^2]^2)
               - sin2Theta^2 (
                   + mst1^2 / mst2^2 (1 - Log[mst1^2 / Q^2])
                   + mst2^2 / mst1^2 (1 - Log[mst2^2 / Q^2])
                 )
               + f1[mt, mg, mst1^2, mst2^2, sin2Theta, Q]
               + f1[mt, mg, mst2^2, mst1^2, -sin2Theta, Q]
           );
           (* Eq. (33) *)
           F2 = unit (
               + 5 Log[mst1^2 / mst2^2]
               - 3 (Log[mst1^2 / Q^2]^2 - Log[mst2^2 / Q^2]^2)
               + cos2Theta2 (
                   + 5 Log[mst1^2 / mst2^2]
                   - (mst1^2 + mst2^2)/(mst1^2 - mst2^2) Log[mst1^2 / mst2^2]^2
                   - 2 / (mst1^2 - mst2^2) (
                       + mst1^2 Log[mst1^2 / Q^2]
                       - mst2^2 Log[mst2^2 / Q^2]
                     ) Log[mst1^2 / mst2^2]
                 )
               + sin2Theta^2 (
                   + mst1^2 / mst2^2 (1 - Log[mst1^2 / Q^2])
                   - mst2^2 / mst1^2 (1 - Log[mst2^2 / Q^2])
                 )
               + f2[mt, mg, mst1^2, mst2^2, sin2Theta, Q]
               - f2[mt, mg, mst2^2, mst1^2, -sin2Theta, Q]
           );
           (* Eq. (33) *)
           F31L = Nc h (2 - (mst1^2 + mst2^2)/(mst1^2 - mst2^2) Log[mst1^2 / mst2^2]);
           F3 = unit (
               + F31L/(h Nc) (3 + 9 cos2Theta2)
               + 4
               - (3 + 13 cos2Theta2) / (mst1^2 - mst2^2) (
                   + mst1^2 Log[mst1^2 / Q^2]
                   - mst2^2 Log[mst2^2 / Q^2]
                 )
               + 3 (mst1^2 + mst2^2)/(mst1^2 - mst2^2) (
                   + Log[mst1^2 / Q^2]^2 - Log[mst2^2 / Q^2]^2
                 )
               - cos2Theta2 (
                   + 4
                   - ((mst1^2 + mst2^2)/(mst1^2 - mst2^2))^2 Log[mst1^2 / mst2^2]^2
                   - 6 (mst1^2 + mst2^2)/(mst1^2 - mst2^2)^2 (
                       + mst1^2 Log[mst1^2 / Q^2]
                       - mst2^2 Log[mst2^2 / Q^2]
                     ) Log[mst1^2 / mst2^2]
                 )
               - sin2Theta^2 (
                   + mst1^2 / mst2^2
                   + mst2^2 / mst1^2
                   + 2 Log[mst1^2 mst2^2 / Q^4]
                   - mst1^2^2 / (mst2^2 (mst1^2 - mst2^2)) Log[mst1^2 / Q^2]
                   + mst2^2^2 / (mst1^2 (mst1^2 - mst2^2)) Log[mst2^2 / Q^2]
                 )
               + f3[mt, mg, mst1^2, mst2^2, sin2Theta, Q]
               + f3[mt, mg, mst2^2, mst1^2, -sin2Theta, Q]
           );
           (* Eq. (35) *)
           DeltaF2 = unit 2 mg / At (
               + Log[mst2^2 / Q^2]^2
               - Log[mst1^2 / Q^2]^2
           );
           (* Eq. (36) *)
           DeltaF3 = unit mg / At (
               + 8
               - 2 (mst1^2 + mst2^2)/(mst1^2 - mst2^2) (
                   + Log[mst2^2 / Q^2]^2
                   - Log[mst1^2 / Q^2]^2
                 )
               + 8 / (mst1^2 - mst2^2) (
                   + mst2^2 Log[mst2^2 / Q^2]
                   - mst1^2 Log[mst1^2 / Q^2]
                 )
           );
           CreateMassMatrixCPEven[F1, F2, F3, DeltaF2, DeltaF3]
          ];

Options[GetMSSMCPEvenHiggsLoopMassMatrix] = {
    loopOrder -> {1,1,1},
    corrections -> {1}
};

GetMSSMCPEvenHiggsLoopMassMatrix[OptionsPattern[]] :=
    (
        OptionValue[loopOrder][[1]] {{0,0}, {0,0}} +
        OptionValue[loopOrder][[2]] OptionValue[corrections][[1]] GetMSSMCPEvenHiggsLoopMassMatrix1LAlphaTAlphaS[] +
        OptionValue[loopOrder][[3]] OptionValue[corrections][[1]] GetMSSMCPEvenHiggsLoopMassMatrix2LAlphaTAlphaS[]
    );

ReplaceStopMasses[] :=
    Module[{mstop = CalculateMStop2[$signMu]},
           {
               mst1 -> Sqrt[mstop[[1]]],
               mst2 -> Sqrt[mstop[[2]]],
               sin2Theta -> mstop[[3]]
           }
          ];

End[];

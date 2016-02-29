BeginPackage["ThreeLoopQCD`", {"SARAH`"}];

GetMTopMSbarOverMTopPole::usage = "Returns the ratio Eq. (10) of
 arxiv:hep-ph/9912391 at the scale Q = Mt_pole.

Example call:

   mTopMSbar = MtPole GetMTopMSbarOverMTopPole[{1,1,1,1}]

The function argument is a list of numbers, which are multiplied by
the tree-level, 1-loop, 2-loop and 3-loop contributions, respectively.
I.e. by setting factors to 0, different loop orders can be disabled.
";

Begin["`Private`"];

GetMTopMSbarOverMTopPole0L[] := 1;

GetMTopMSbarOverMTopPole1L[] :=
    Module[{CF, colorPosition, as},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[SARAH`TopQuark, colorPosition];
           as = SARAH`strongCoupling^2 / (4 Pi);
           -CF as / Pi
          ];

(* 2-loop coefficients *)
d12 = 7/128 - 3/4 Zeta[3] + 1/2 Pi^2 Log[2] - 5/16 Pi^2;
d22 = -1111/384 + 3/8 Zeta[3] - 1/4 Pi^2 Log[2] + Pi^2/12;
d32 = 71/96 + Pi^2/12;
d42 = 143/96 - 1/6 Pi^2;

(* 3-loop coefficients *)
a4 = PolyLog[4, 1/2];
d13 = (
    -2969/768 - 1/16 Pi^2 Zeta[3] 81/16 Zeta[3] + 5/8 Zeta[5] + 29/4 Pi^2 Log[2]
    +1/2 Pi^2 Log[2]^2 - 613/192 Pi^2 - 1/48 Pi^4 - 1/2 Log[2]^4 - 12 a4
);
d23 = (
    13189/4608 - 19/16 Pi^2 Zeta[3] - 773/96 Zeta[3] + 45/16 Zeta[5] - 31/72 Pi^2 Log[2]
    -31/36 Pi^2 Log[2]^2 + 509/576 Pi^2 + 65/432 Pi^4 - 1/18 Log[2]^4 - 4/3 a4
);
d33 = (
    -1322545/124416 + 51/64 Pi^2 Zeta[3] + 1343/288 Zeta[3] - 65/32 Zeta[5]
    -115/72 Pi^2 Log[2] + 11/36 Pi^2 Log[2]^2 - 1955/3456 Pi^2 - 179/3456 Pi^4
    +11/72 Log[2]^4 + 11/3 a4
);
d43 = (
    1283/576 + 55/24 Zeta[3] - 11/9 Pi^2 Log[2] + 2/9 Pi^2 Log[2]^2 + 13/18 Pi^2
    -119/2160 Pi^4 + 1/9 Log[2]^4 + 8/3 a4
);
d53 = (
    1067/576 - 53/24 Zeta[3] + 8/9 Pi^2 Log[2] - 1/9 Pi^2 Log[2]^2 - 85/108 Pi^2
    +91/2160 Pi^4 + 1/9 Log[2]^4 + 8/3 a4
);
d63 = (
    70763/15552 + 89/144 Zeta[3] + 11/18 Pi^2 Log[2] - 1/9 Pi^2 Log[2]^2
    +175/432 Pi^2 + 19/2160 Pi^4 - 1/18 Log[2]^4 - 4/3 a4
);
d73 = (
    144959/15552 + 1/8 Pi^2 Zeta[3] - 109/144 Zeta[3] - 5/8 Zeta[5] + 32/9 Pi^2 Log[2]
    +1/18 Pi^2 Log[2]^2 - 449/144 Pi^2 - 43/1080 Pi^4 - 1/18 Log[2]^4 - 4/3 a4
);
d83 = -5917/3888 + 2/9 Zeta[3] + 13/108 Pi^2;
d93 = -9481/7776 + 11/18 Zeta[3] + 4/135 Pi^2;
d103 = -2353/7776 - 7/18 Zeta[3] - 13/108 Pi^2;

GetMTopMSbarOverMTopPole2L[quark_:SARAH`TopQuark, NH_:1, NL_:5] :=
    Module[{CF, CA, TR, colorPosition, as},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           CA = SA`Casimir[SARAH`VectorG, colorPosition];
           TR = 1/2;
           as = SARAH`strongCoupling^2 / (4 Pi);
           CF (as/Pi)^2 (CF d12 + CA d22 + TR NL d32 + TR NH d42)
          ];

GetMTopMSbarOverMTopPole3L[quark_:SARAH`TopQuark, NH_:1, NL_:5] :=
    Module[{CF, CA, TR, colorPosition, as},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           CA = SA`Casimir[SARAH`VectorG, colorPosition];
           TR = 1/2;
           as = SARAH`strongCoupling^2 / (4 Pi);
           CF (as/Pi)^3 (
               CF^2 d13 + CF CA d23 + CA^2 d33 +
               CF TR NL d43 + CF TR NH d53 + CA TR NL d63 +
               CA TR NH d73 + TR^2 NL NH d83 + TR^2 NH^2 d93 +
               TR^2 NL^2 d103
           )
          ];

GetMTopMSbarOverMTopPole[loopOrder_List:{1,1,1,1}] :=
    Module[{},
           Assert[Length[loopOrder] == 4];
           (
               loopOrder[[1]] GetMTopMSbarOverMTopPole0L[] +
               loopOrder[[2]] GetMTopMSbarOverMTopPole1L[] +
               loopOrder[[3]] GetMTopMSbarOverMTopPole2L[] +
               loopOrder[[4]] GetMTopMSbarOverMTopPole3L[]
           )
          ];

End[];

EndPackage[];


BeginPackage["TwoLoop`", {"SARAH`"}];

GetDeltaMQCD::usage="Returns two-loop QCD contributions to \Delta m_f
in the DRbar scheme.  Taken from hep-ph/0210258 Eq. (60)-(61).";

Begin["Private`"];

GetDeltaMQCD[quark_, renScale_] :=
    Module[{},
           Print["Warning: two-loop mass correction for ", quark, " not available!"];
           Return[0];
          ];

GetDeltaMQCD[quark:SARAH`TopQuark, renScale_] :=
    Module[{CF, CA, colorPosition, alphaStrong, mf, log},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           CA = SA`Casimir[SARAH`VectorG, colorPosition];
           alphaStrong = SARAH`strongCoupling^2 / (4 Pi);
           mf = FlexibleSUSY`M[quark];
           log = Log[(mf^2) / (renScale^2)];
           mf CF (alphaStrong / (4 Pi))^2 (-43 - 12 Zeta[2] + 26 log - 6 log^2
                                           + CF (-59/8 + 30 Zeta[2] - 48 Log[2] Zeta[2]
                                                 + 12 Zeta[3] + (3/2) log + (9/2) log^2)
                                           + CA (1093/24 - 8 Zeta[2] + 24 Log[2] Zeta[2]
                                                 - 6 Zeta[3] - (179/6) log + 11/2 log^2)
                                          )
          ];

GetDeltaMQCD[quark:SARAH`BottomQuark, renScale_] :=
    Module[{CF, CA, colorPosition, alphaStrong, mf, mt, log, logMt},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           CA = SA`Casimir[SARAH`VectorG, colorPosition];
           alphaStrong = SARAH`strongCoupling^2 / (4 Pi);
           mf = FlexibleSUSY`M[quark];
           mt = FlexibleSUSY`M[SARAH`TopQuark];
           log = Log[(mf^2) / (mt^2)];
           logMt = Log[(mf^2) / (^2)]
           mf CF (alphaStrong / (4 Pi))^2 (-623/18 - 8 Zeta[2] + 26 log - 6 log^2
                                           + CF (-59/8 + 30 Zeta[2] - 48 Log[2] Zeta[2]
                                                 + 12 Zeta[3] + (3/2) log + (9/2) log^2)
                                           + CA (1093/24 - 8 Zeta[2] + 24 Log[2] Zeta[2]
                                                 - 6 Zeta[3] - (179/6) log + 11/2 log^2)
                                           - 13/3 logMt + logMt^2
                                          )
          ];

End[];

EndPackage[];

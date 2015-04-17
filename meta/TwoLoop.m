
BeginPackage["TwoLoop`", {"SARAH`"}];

GetDeltaMOverMQCDOneLoop::usage="Returns one-loop QCD contributions to
Delta m_f/m_f in the given renormalization scheme.";

GetDeltaMOverMQCDOneLoopDRbar::usage="Returns one-loop QCD contributions to
Delta m_f/m_f in the DR-bar scheme.  Taken from hep-ph/0210258
Eq. (58)";

GetDeltaMOverMQCDOneLoopMSbar::usage="Returns one-loop QCD contributions to
Delta m_f/m_f in the MS-bar scheme.  Taken from hep-ph/9803493 Eq. (16)";

GetDeltaMOverMQCDTwoLoop::usage="Returns two-loop QCD contributions to
Delta m_f/m_f in the given renormalization scheme.";

GetDeltaMOverMQCDTwoLoopDRbar::usage="Returns two-loop QCD contributions to
Delta m_f/m_f in the DR-bar scheme.  Taken from hep-ph/0210258
Eq. (60)-(61).";

GetDeltaMOverMQCDTwoLoopMSbar::usage="Returns two-loop QCD contributions to
Delta m_f/m_f in the MS-bar scheme.  Taken from hep-ph/9803493 Eq. (17).";

Begin["`Private`"];

SelectRenormalizationScheme::UnknownRenormalizationScheme = "Unknown\
 renormalization scheme `1`.";

(* ******* two-loop ******* *)

GetDeltaMOverMQCDTwoLoop[quark_, renScale_, renScheme_] :=
    Switch[renScheme,
           FlexibleSUSY`DRbar, GetDeltaMOverMQCDTwoLoopDRbar[quark,renScale],
           FlexibleSUSY`MSbar, GetDeltaMOverMQCDTwoLoopMSbar[quark,renScale],
           _, Message[SelectRenormalizationScheme::UnknownRenormalizationScheme, renScheme];
              Quit[1];
          ];

GetDeltaMOverMQCDTwoLoopDRbar[quark_ /; quark === SARAH`TopQuark, renScale_] :=
    Module[{CF, CA, colorPosition, alphaStrong, mf, log, result},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           CA = SA`Casimir[SARAH`VectorG, colorPosition];
           alphaStrong = SARAH`strongCoupling^2 / (4 Pi);
           mf = FlexibleSUSY`M[quark];
           log = Log[(mf / renScale)^2];
           result = CF (alphaStrong / (4 Pi))^2 (
               -43 - 12 Zeta[2] + 26 log - 6 log^2
               + CF (-59/8 + 30 Zeta[2] - 48 Log[2] Zeta[2]
                     + 12 Zeta[3] + (3/2) log + (9/2) log^2)
               + CA (1093/24 - 8 Zeta[2] + 24 Log[2] Zeta[2]
                     - 6 Zeta[3] - (179/6) log + 11/2 log^2));
           Simplify[N[result]]
          ];

GetDeltaMOverMQCDTwoLoopDRbar[quark_ /; quark === SARAH`BottomQuark, renScale_] :=
    Module[{CF, CA, colorPosition, alphaStrong, mf, mt, log, logMt, result},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           CA = SA`Casimir[SARAH`VectorG, colorPosition];
           alphaStrong = SARAH`strongCoupling^2 / (4 Pi);
           mf = FlexibleSUSY`M[quark];
           mt = FlexibleSUSY`M[SARAH`TopQuark];
           log = Log[(mf / renScale)^2];
           logMt = Log[(mf / mt)^2];
           result = CF (alphaStrong / (4 Pi))^2 (
               -623/18 - 8 Zeta[2] + 26 log - 6 log^2
               + CF (-59/8 + 30 Zeta[2] - 48 Log[2] Zeta[2]
                     + 12 Zeta[3] + (3/2) log + (9/2) log^2)
               + CA (1093/24 - 8 Zeta[2] + 24 Log[2] Zeta[2]
                     - 6 Zeta[3] - (179/6) log + 11/2 log^2)
               - 13/3 logMt + logMt^2);
           Simplify[N[result]]
          ];

GetDeltaMOverMQCDTwoLoopDRbar[_, _] := 0;

GetDeltaMOverMQCDTwoLoopMSbar[quark_, renScale_] :=
    Module[{CF, CA, colorPosition, alphaStrong, mf, log,
            t = 1/2, nf = 6, I31, result},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           CA = SA`Casimir[SARAH`VectorG, colorPosition];
           alphaStrong = SARAH`strongCoupling^2 / (4 Pi);
           mf = FlexibleSUSY`M[quark];
           log = Log[(renScale / mf)^2];
           I31 = 3/2 Zeta[3] - 6 Zeta[2] Log[2];
           result = (alphaStrong / (4 Pi))^2 (
               CF CA (1111/24 - 8 Zeta[2] - 4 I31 + 185/6 log + 11/2 log^2)
               - CF t nf (71/6 + 8 Zeta[2] + 26/3 log + 2 log^2)
               + CF^2 (121/8 + 30 Zeta[2] + 8 I31 + 27/2 log + 9/2 log^2)
               - 12 CF t (1 - 2 Zeta[2]));
           Simplify[N[result]]
          ];

(* ******* one-loop ******* *)

GetDeltaMOverMQCDOneLoop[quark_, renScale_, renScheme_] :=
    Switch[renScheme,
           FlexibleSUSY`DRbar, GetDeltaMOverMQCDOneLoopDRbar[quark,renScale],
           FlexibleSUSY`MSbar, GetDeltaMOverMQCDOneLoopMSbar[quark,renScale],
           _, Message[SelectRenormalizationScheme::UnknownRenormalizationScheme, renScheme];
              Quit[1];
          ];

GetDeltaMOverMQCDOneLoopDRbar[quark_, renScale_] :=
    Module[{CF, colorPosition, alphaStrong, mf, result},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           alphaStrong = SARAH`strongCoupling^2 / (4 Pi);
           mf = FlexibleSUSY`M[quark];
           result = CF alphaStrong / (4 Pi) (5 - 3 Log[(mf / renScale)^2]);
           Simplify[N[result]]
          ];

GetDeltaMOverMQCDOneLoopMSbar[quark_, renScale_] :=
    Module[{CF, colorPosition, alphaStrong, mf, result},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           alphaStrong = SARAH`strongCoupling^2 / (4 Pi);
           mf = FlexibleSUSY`M[quark];
           result = CF alphaStrong / (4 Pi) (4 - 3 Log[(mf / renScale)^2]);
           Simplify[N[result]]
          ];

End[];

EndPackage[];

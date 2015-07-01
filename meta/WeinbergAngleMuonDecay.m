(* ::Package:: *)

BeginPackage["WeinbergAngleMuonDecay`", {"SARAH`", "CConversion`", "Parameters`", "TextFormatting`"}];

deltaRho2LoopSM::usage="";

Begin["`Private`"];

deltaRho2LoopSM[]:=Module[{gY, alphaDRbar, xt, gfermi, mt, hmix, expr, sinThetaW, mw, mz, RHO2, result},
                          gY = g1 Sqrt[0.6];
                          alphaDRbar = gY^2 Susyno`LieGroups`g2^2 / (4 Pi (gY^2 + Susyno`LieGroups`g2^2));
                          xt = 3.0 / (8 Pi^2 Sqrt[2]) gfermi mt^2;
                          hmix = (FlexibleSUSY`ZH[0,1] / Sin[ArcTan[FlexibleSUSY`vu/FlexibleSUSY`vd]])^2;
                          expr = alphaDRbar g3^2/(16 Pi^3 sinThetaW^2)(-2.145 mt^2/mw^2 + 1.262 Log[mt/mz] - 2.24 - 0.85 mz^2/mt^2) + xt^2 hmix RHO2[FlexibleSUSY`M[hh[0]]/mt] / 3;
                          result = Parameters`CreateLocalConstRefs[expr];
                          result = result <> "\n" <> "\n";
                          result = result <> CConversion`RValueToCFormString[expr];
                          Return[result];
                          ];

End[];

EndPackage[];

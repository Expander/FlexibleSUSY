(* ::Package:: *)

BeginPackage["WeinbergAngleMuonDecay`", {"SARAH`", "CConversion`", "Parameters`", "TextFormatting`"}];

deltaRho2LoopSM::usage="";
deltaR2LoopSM::usage="";

Begin["`Private`"];

extPars={sinThetaW, rhohat, gfermi, mw, mz, mt, RHO2, deltaR1Loop};
Do[Format[extPars[[i]],CForm]=Format[ToString[extPars[[i]]],OutputForm],{i,Length[extPars]}];

deltaRho2LoopSM[]:=Module[{gY, alphaDRbar, xt, hmix, expr, result},
                          gY = SARAH`g1 Sqrt[0.6];
                          alphaDRbar = gY^2 Susyno`LieGroups`g2^2 / (4 Pi (gY^2 + Susyno`LieGroups`g2^2));
                          xt = 3 / (8 Pi^2 Sqrt[2]) gfermi mt^2;
                          hmix = (FlexibleSUSY`ZH[0,1] / Sin[ArcTan[FlexibleSUSY`vu/FlexibleSUSY`vd]])^2;
                          expr = alphaDRbar SARAH`g3^2/(16 Pi^3 sinThetaW^2)(-2.145 mt^2/mw^2 + 1.262 Log[mt/mz] - 2.24 - 0.85 mz^2/mt^2) + xt^2 hmix RHO2[FlexibleSUSY`M[hh[0]]/mt] / 3;
                          result = Parameters`CreateLocalConstRefs[expr] <> "\n";
                          result = result <> "deltaRho2LoopSM = ";
                          result = result <> CConversion`RValueToCFormString[expr] <> ";";
                          Return[result];
                          ];

deltaR2LoopSM[]:=Module[{gY, alphaDRbar, xt, hmix, expr, result},
                        gY = SARAH`g1 Sqrt[0.6];
                        alphaDRbar = gY^2 Susyno`LieGroups`g2^2 / (4 Pi (gY^2 + Susyno`LieGroups`g2^2));
                        xt = 3 / (8 Pi^2 Sqrt[2]) gfermi mt^2;
                        hmix = (FlexibleSUSY`ZH[0,1] / Sin[ArcTan[FlexibleSUSY`vu/FlexibleSUSY`vd]])^2;
                        expr = alphaDRbar SARAH`g3^2/(16 Pi^3 sinThetaW^2 (1 - sinThetaW^2))(2.145 mt^2/mz^2 + 0.575 Log[mt/mz] - 0.224 - 0.144 mz^2/mt^2) - xt^2 hmix RHO2[FlexibleSUSY`M[hh[0]]/mt] (1 - deltaR1Loop) rhohat / 3;
                        result = Parameters`CreateLocalConstRefs[expr] <> "\n";
                        result = result <> "deltaR2LoopSM = ";
                        result = result <> CConversion`RValueToCFormString[expr] <> ";";
                        Return[result];
                        ];

End[];

EndPackage[];

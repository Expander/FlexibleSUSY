(* ::Package:: *)

BeginPackage["WeinbergAngleMuonDecay`", {"SARAH`", "CConversion`", "Parameters`", "TextFormatting`"}];

deltaRhoHat2LoopSM::usage="";
deltaRHat2LoopSM::usage="";

Begin["`Private`"];

extPars={sinThetaW, rhohat, gfermi, mw, mz, mt, RHO2, deltaRHat1Loop};
Do[Format[extPars[[i]],CForm]=Format[ToString[extPars[[i]]],OutputForm],{i,Length[extPars]}];

(*formula according to (C.6) from hep-ph/9606211*)
deltaRhoHat2LoopSM[]:=Module[{gY, alphaDRbar, xt, hmix, expr, result},
                             gY = SARAH`hyperchargeCoupling FlexibleSUSY`GUTNormalization[SARAH`hyperchargeCoupling];
                             alphaDRbar = gY^2 SARAH`leftCoupling^2 / (4 Pi (gY^2 + SARAH`leftCoupling^2));
                             xt = 3 / (8 Pi^2 Sqrt[2]) gfermi mt^2;
                             hmix = (SARAH`HiggsMixingMatrix[0,1] / Sin[ArcTan[FlexibleSUSY`vu/FlexibleSUSY`vd]])^2;
                             expr = alphaDRbar SARAH`strongCoupling^2/(16 Pi^3 sinThetaW^2)(-2.145 mt^2/mw^2 + 1.262 Log[mt/mz] - 2.24 - 0.85 mz^2/mt^2) + xt^2 hmix RHO2[FlexibleSUSY`M[SARAH`HiggsBoson[0]]/mt] / 3;
                             result = Parameters`CreateLocalConstRefs[expr] <> "\n";
                             result = result <> "deltaRhoHat2LoopSM = ";
                             result = result <> CConversion`RValueToCFormString[expr] <> ";";
                             Return[result];
                             ];

(*formula according to (C.5) from hep-ph/9606211*)
deltaRHat2LoopSM[]:=Module[{gY, alphaDRbar, xt, hmix, expr, result},
                           gY = SARAH`hyperchargeCoupling FlexibleSUSY`GUTNormalization[SARAH`hyperchargeCoupling];
                           alphaDRbar = gY^2 SARAH`leftCoupling^2 / (4 Pi (gY^2 + SARAH`leftCoupling^2));
                           xt = 3 / (8 Pi^2 Sqrt[2]) gfermi mt^2;
                           hmix = (SARAH`HiggsMixingMatrix[0,1] / Sin[ArcTan[FlexibleSUSY`vu/FlexibleSUSY`vd]])^2;
                           expr = alphaDRbar SARAH`strongCoupling^2/(16 Pi^3 sinThetaW^2 (1 - sinThetaW^2))(2.145 mt^2/mz^2 + 0.575 Log[mt/mz] - 0.224 - 0.144 mz^2/mt^2) - xt^2 hmix RHO2[FlexibleSUSY`M[SARAH`HiggsBoson[0]]/mt] (1 - deltaRHat1Loop) rhohat / 3;
                           result = Parameters`CreateLocalConstRefs[expr] <> "\n";
                           result = result <> "deltaRHat2LoopSM = ";
                           result = result <> CConversion`RValueToCFormString[expr] <> ";";
                           Return[result];
                           ];

End[];

EndPackage[];

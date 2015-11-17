(* ::Package:: *)

BeginPackage["WeinbergAngleMuonDecay`", {"SARAH`", "CConversion`", "Parameters`", "TextFormatting`", "TreeMasses`", "WeinbergAngle`"}];

deltaRhoHat2LoopSM::usage="";
deltaRHat2LoopSM::usage="";
RhoHatTree::usage="";

Begin["`Private`"];

(*function taken from WeinbergAngle.m*)
FindMass[massMatrices_List, particle_] :=
    Module[{massExpr},
           massExpr = Cases[massMatrices, TreeMasses`FSMassMatrix[{mass_}, particle, ___] :> mass];
           If[Head[massExpr] =!= List || massExpr === {},
              Print["Error: Could not find mass of ", particle,
                    " in masses list."];
              Return[0];
             ];
           Return[massExpr[[1]] /. SARAH`Weinberg[] -> SARAH`Weinberg];
          ];

extPars={SINTHETAW, RHOHATRATIO, GFERMI, MW, MZ, MT, RHO2, DELTARHAT1LOOP, PIZZTMZ};
Do[Format[extPars[[i]],CForm]=Format[ToString[extPars[[i]]],OutputForm],{i,Length[extPars]}];

(*formula according to (C.6) from hep-ph/9606211*)
deltaRhoHat2LoopSM[]:=Module[{gY, alphaDRbar, xt, hmix, expr, result},
                             gY = SARAH`hyperchargeCoupling FlexibleSUSY`GUTNormalization[SARAH`hyperchargeCoupling];
                             alphaDRbar = gY^2 SARAH`leftCoupling^2 / (4 Pi (gY^2 + SARAH`leftCoupling^2));
                             xt = 3 / (8 Pi^2 Sqrt[2]) GFERMI MT^2;
                             hmix = (SARAH`HiggsMixingMatrix[0,1] / Sin[ArcTan[FlexibleSUSY`vu/FlexibleSUSY`vd]])^2;
                             expr = (alphaDRbar SARAH`strongCoupling^2/(16 Pi^3 SINTHETAW^2)(-2.145 MT^2/MW^2 + 1.262 Log[MT/MZ] - 2.24 - 0.85 MZ^2/MT^2) + xt^2 hmix RHO2[FlexibleSUSY`M[SARAH`HiggsBoson[0]]/MT] / 3) / (1 + PIZZTMZ / MZ^2);
                             result = Parameters`CreateLocalConstRefs[expr] <> "\n";
                             result = result <> "deltaRhoHat2LoopSM = ";
                             result = result <> CConversion`RValueToCFormString[expr] <> ";";
                             Return[result];
                             ];

(*formula according to (C.5) from hep-ph/9606211*)
deltaRHat2LoopSM[]:=Module[{gY, alphaDRbar, xt, hmix, expr, result},
                           gY = SARAH`hyperchargeCoupling FlexibleSUSY`GUTNormalization[SARAH`hyperchargeCoupling];
                           alphaDRbar = gY^2 SARAH`leftCoupling^2 / (4 Pi (gY^2 + SARAH`leftCoupling^2));
                           xt = 3 / (8 Pi^2 Sqrt[2]) GFERMI MT^2;
                           hmix = (SARAH`HiggsMixingMatrix[0,1] / Sin[ArcTan[FlexibleSUSY`vu/FlexibleSUSY`vd]])^2;
                           expr = alphaDRbar SARAH`strongCoupling^2/(16 Pi^3 SINTHETAW^2 (1 - SINTHETAW^2))(2.145 MT^2/MZ^2 + 0.575 Log[MT/MZ] - 0.224 - 0.144 MZ^2/MT^2) - xt^2 hmix RHO2[FlexibleSUSY`M[SARAH`HiggsBoson[0]]/MT] (1 - DELTARHAT1LOOP) RHOHATRATIO / 3;
                           result = Parameters`CreateLocalConstRefs[expr] <> "\n";
                           result = result <> "deltaRHat2LoopSM = ";
                           result = result <> CConversion`RValueToCFormString[expr] <> ";";
                           Return[result];
                           ];

(*calculates tree-level value of rhohat parameter from Z and W tree-level DRbar masses*)
RhoHatTree[massMatrices_List]:=Module[{Ztreemass, Wtreemass, expr, result},
                                       Ztreemass = FindMass[massMatrices,SARAH`VectorZ];
                                       Wtreemass = FindMass[massMatrices,SARAH`VectorW];
                                       expr = Simplify[Wtreemass / Ztreemass / Cos[SARAH`Weinberg]^2/.SARAH`Weinberg->WeinbergAngle`ExpressWeinbergAngleInTermsOfGaugeCouplings[massMatrices]];
                                       result = Parameters`CreateLocalConstRefs[expr] <> "\n";
                                       result = result <> "rhohat_tree = ";
                                       result = result <> CConversion`RValueToCFormString[expr] <> ";";
                                       Return[result];
                                       ];

End[];

EndPackage[];

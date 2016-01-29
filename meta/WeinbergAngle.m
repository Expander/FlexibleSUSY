
BeginPackage["WeinbergAngle`", {"SARAH`", "CConversion`", "Parameters`", "TreeMasses`"}];

ExpressWeinbergAngleInTermsOfGaugeCouplings::usage="";
deltaRhoHat2LoopSM::usage="";
deltaRHat2LoopSM::usage="";
RhoHatTree::usage="";

Begin["`Private`"];

FindMassW2[masses_List] :=
    FindMass2[masses, SARAH`VectorW];

FindMassZ2[masses_List] :=
    FindMass2[masses, SARAH`VectorZ];

FindMass2[masses_List, particle_] :=
    Module[{massExpr},
           massExpr = Cases[masses, TreeMasses`FSMassMatrix[{mass_}, particle, ___] :> mass];
           If[Head[massExpr] =!= List || massExpr === {},
              Print["Error: Could not find mass of ", particle,
                    " in masses list."];
              Return[0];
             ];
           Return[massExpr[[1]] /. SARAH`Weinberg[] -> SARAH`Weinberg];
          ];

(*extracts squared tree-level mass of Z boson before mixing with additional Z`*)
UnmixedZMass2[] :=
    Module[{ZMassMatrix, extraGaugeCouplings, submatrixList, submatrix, mass2Eigenvalues},
           ZMassMatrix = SARAH`MassMatrix[SARAH`VectorZ];
           Assert[MatrixQ[ZMassMatrix]];
           extraGaugeCouplings = Cases[SARAH`Gauge, x_ /; FreeQ[x, SARAH`hypercharge] && FreeQ[x, SARAH`left] && FreeQ[x, SARAH`color] :> x[[4]]];
           submatrixList = ZMassMatrix[[#, #]] & /@ Flatten[Table[{i, j}, {i, 1, Length[ZMassMatrix]}, {j, i + 1, Length[ZMassMatrix]}], 1];
           submatrix = Cases[submatrixList, x_ /; And @@ (FreeQ[x, #] & /@ extraGaugeCouplings)];
           If[Length[submatrix] != 1, Print["Error: Photon-Z mass matrix could not be identified"]; Return[0]];
           mass2Eigenvalues = Eigenvalues[submatrix];
           If[Length[mass2Eigenvalues] != 2 || !MemberQ[mass2Eigenvalues, 0], Print["Error: Determination of UnmixedZMass2 failed"]; Return[0]];
           Return[Select[mass2Eigenvalues, # =!= 0 &][[1]] /. Parameters`ApplyGUTNormalization[]];
          ];

(*calculates rho_0 from SU(2)_L representations of the Higgs multipletts as in (16) from 0801.1345 [hep-ph]*)
RhoZero[] :=
    Module[{hyperchargePos, leftPos, vevlist},
           hyperchargePos = Position[SARAH`Gauge, x_ /; !FreeQ[x, SARAH`hypercharge], {1}][[1, 1]];
           leftPos = Position[SARAH`Gauge, x_ /; !FreeQ[x, SARAH`left], {1}][[1, 1]];
           (* what about non Higgs vevs? non scalar vevs possible? *)
           vevlist = SARAH`DEFINITION[SARAH`EWSB][SARAH`VEVs][[All, 1;;2]];
           vevlist = vevlist /. {Sfieldname_Symbol, vevinfo_List} :> {ToExpression[StringReplace[ToString[Sfieldname], StartOfString ~~ "S" -> ""]], vevinfo};
           (* extract isospin from SU(2)_left representation and its third component from Gell-Mann-Nishijima formula with given hypercharge and electric charge = 0 *)
           vevlist = vevlist /. {fieldname_Symbol, vevinfo_List} :> Flatten[{vevinfo, Cases[SARAH`SuperFields, x_ /; !FreeQ[x[[3]], fieldname] :> {(x[[3 + leftPos]] - 1) / 2, -x[[3 + hyperchargePos]]}]}];
           Return[Simplify[Plus @@ ((#[[3]]^2 - #[[4]]^2 + #[[3]]) (#[[1]] #[[2]] Sqrt[2])^2 & /@ vevlist) / Plus @@ (2 #[[4]]^2 (#[[1]] #[[2]] Sqrt[2])^2 & /@ vevlist)]];
          ];

ExpressWeinbergAngleInTermsOfGaugeCouplings[masses_List] :=
    Module[{solution},
           Print["Expressing Weinberg angle in terms of model parameters ..."];
           solution = ArcCos[Sqrt[FindMassW2[masses] / UnmixedZMass2[] / RhoZero[]]];
           Return[Simplify[solution]];
          ];

extPars={SINTHETAW, RHOHATRATIO, GFERMI, MW, MZ, MT, RHO2, DELTARHAT1LOOP, PIZZTMZ};
Do[Format[extPars[[i]],CForm]=Format[ToString[extPars[[i]]],OutputForm],{i,Length[extPars]}];

(*formula according to (C.6) from hep-ph/9606211*)
deltaRhoHat2LoopSM[]:=
    Module[{gY, alphaDRbar, xt, hmix, expr, result},
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
deltaRHat2LoopSM[]:=
    Module[{gY, alphaDRbar, xt, hmix, expr, result},
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

(*calculates tree-level value of rhohat parameter from umixed and mixed Z mass as well as RhoZero*)
RhoHatTree[massMatrices_List]:=
    Module[{Zmass2unmixed, Zmass2mixed, expr, result},
           Zmass2unmixed = UnmixedZMass2[];
           Zmass2mixed = FindMassZ2[massMatrices];
           expr = Simplify[RhoZero[] Zmass2unmixed / Zmass2mixed /. SARAH`Weinberg -> ExpressWeinbergAngleInTermsOfGaugeCouplings[massMatrices], SARAH`hyperchargeCoupling > 0 && SARAH`leftCoupling > 0];
           result = Parameters`CreateLocalConstRefs[expr] <> "\n";
           result = result <> "rhohat_tree = ";
           result = result <> CConversion`RValueToCFormString[expr] <> ";";
           Return[result];
          ];

End[];

EndPackage[];

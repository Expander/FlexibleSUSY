
BeginPackage["WeinbergAngle`", {"SARAH`", "CConversion`", "Parameters`", "TreeMasses`"}];

ExpressWeinbergAngleInTermsOfGaugeCouplings::usage="";
deltaRhoHat2LoopSM::usage="";
deltaRHat2LoopSM::usage="";
RhoHatTree::usage="";

Begin["`Private`"];

InputFormOfNonStrings[a_String] := a;
InputFormOfNonStrings[a_] := InputForm[a];

DebugPrint[msg___] :=
    If[FlexibleSUSY`FSDebugOutput,
       Print["Debug<WeinbergAngle>: ", Sequence @@ InputFormOfNonStrings /@ {msg}]];

FindWeinbergAngleDef[] :=
    Parameters`FindSymbolDef[SARAH`Weinberg];

FindMassW[masses_List] :=
    FindMass[masses, SARAH`VectorW];

FindMassZ[masses_List] :=
    FindMass[masses, SARAH`VectorZ];

FindMass[masses_List, particle_] :=
    Module[{massExpr},
           massExpr = Cases[masses, TreeMasses`FSMassMatrix[{mass_}, particle, ___] :> mass];
           If[Head[massExpr] =!= List || massExpr === {},
              Print["Error: Could not find mass of ", particle,
                    " in masses list."];
              Return[0];
             ];
           Return[massExpr[[1]] /. SARAH`Weinberg[] -> SARAH`Weinberg];
          ];

SolvesWeinbergEq[eq_, expr_] :=
    Module[{insertedEq},
           insertedEq = TimeConstrained[
               Simplify[eq /. SARAH`Weinberg -> expr],
               FlexibleSUSY`FSSolveWeinbergAngleTimeConstraint,
               False];
           insertedEq === True
          ];

ExpressWeinbergAngleInTermsOfGaugeCouplings::nonSMDef = "Warning: The \
Standard Model definition of the Weinberg angle `1` does not solve the \
system of equations `2`.";

ExpressWeinbergAngleInTermsOfGaugeCouplings::noSolution = "Error: \
could not express the Weinberg angle `1` in terms of the gauge \
couplings.  The equation which could not be solved is: `2`.  I'll use \
the Standard Model definition `3` instead.";

ExpressWeinbergAngleInTermsOfGaugeCouplings[masses_List] :=
    Module[{weinbergDef, eqs, reducedEq, solution, smValue},
           Print["Expressing Weinberg angle in terms of model parameters ..."];
           (* SM value of the Weinberg angle *)
           smValue = Simplify[
               ArcTan[SARAH`hyperchargeCoupling / SARAH`leftCoupling] /.
               Parameters`ApplyGUTNormalization[]];
           (* search weinberg angle definition *)
           weinbergDef = FindWeinbergAngleDef[] /. Parameters`ApplyGUTNormalization[];
           If[weinbergDef === 0,
              Print["Warning: Weinberg angle is defined to be zero."];
              Return[weinbergDef];
             ];
           (* Assumption: the masses contain GUT normalized gauge
              couplings, but the Weinberg angle does not (because we
              take it directly from SARAH`ParameterDefinitions, where
              no GUT normalization has been applied so far)*)
           eqs = {SARAH`Weinberg == weinbergDef,
                  SARAH`Mass[SARAH`VectorW]^2 == FindMassW[masses],
                  SARAH`Mass[SARAH`VectorZ]^2 == FindMassZ[masses]
                 };
           DebugPrint["The 3 equations to determine the Weinberg angle are: ",
                      eqs];
           reducedEq = Eliminate[eqs, {SARAH`Mass[SARAH`VectorW],
                                       SARAH`Mass[SARAH`VectorZ]}];
           DebugPrint["Elimination of ",
                      {SARAH`Mass[SARAH`VectorW], SARAH`Mass[SARAH`VectorZ]},
                      " yields: ", reducedEq];
           (* Try Standard Model definition first *)
           If[SolvesWeinbergEq[reducedEq, smValue],
              Return[smValue];,
              Message[ExpressWeinbergAngleInTermsOfGaugeCouplings::nonSMDef,
                      InputForm[SARAH`Weinberg == smValue],
                      InputForm[eqs]];
             ];
           DebugPrint["Solving equation for ", SARAH`Weinberg,
                      " using a time constraint of ",
                      FlexibleSUSY`FSSolveWeinbergAngleTimeConstraint];
           Off[Solve::ifun];
           solution = TimeConstrained[Solve[reducedEq, SARAH`Weinberg, Reals],
                                      FlexibleSUSY`FSSolveWeinbergAngleTimeConstraint,
                                      {}];
           On[Solve::ifun];
           If[Head[solution] =!= List || solution === {},
              Message[ExpressWeinbergAngleInTermsOfGaugeCouplings::noSolution,
                      InputForm[SARAH`Weinberg],
                      InputForm[reducedEq],
                      InputForm[SARAH`Weinberg == smValue]];
              Return[smValue];
             ];
           solution = Simplify[#[[1,2]],
                               Assumptions ->
                               SARAH`Weinberg > 0 && Element[SARAH`Weinberg, Reals] &&
                               SARAH`hyperchargeCoupling > 0 && Element[SARAH`hyperchargeCoupling, Reals] &&
                               SARAH`leftCoupling > 0 && Element[SARAH`leftCoupling, Reals]]& /@ solution;
           solution = Select[solution, !NumericQ[#]&];
           If[solution === {},
              Message[ExpressWeinbergAngleInTermsOfGaugeCouplings::noSolution,
                      InputForm[SARAH`Weinberg],
                      InputForm[reducedEq],
                      InputForm[SARAH`Weinberg == smValue]];
              Return[smValue];
             ];
           Sort[solution, ByteCount[#1] < ByteCount[#2]&][[1]]
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
                                       Ztreemass = FindMassZ[massMatrices];
                                       Wtreemass = FindMassW[massMatrices];
                                       expr = Apart[Simplify[Wtreemass / Ztreemass / Cos[SARAH`Weinberg]^2/.SARAH`Weinberg->ExpressWeinbergAngleInTermsOfGaugeCouplings[massMatrices]]];
                                       result = Parameters`CreateLocalConstRefs[expr] <> "\n";
                                       result = result <> "rhohat_tree = ";
                                       result = result <> CConversion`RValueToCFormString[expr] <> ";";
                                       Return[result];
                                       ];

End[];

EndPackage[];

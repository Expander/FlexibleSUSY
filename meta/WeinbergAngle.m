
BeginPackage["WeinbergAngle`", {"SARAH`", "CConversion`", "Parameters`", "SelfEnergies`", "TextFormatting`", "ThresholdCorrections`", "TreeMasses`", "Vertices`"}];

GetBottomMass::usage="";
GetTopMass::usage="";

ExpressWeinbergAngleInTermsOfGaugeCouplings::usage="";
DeltaRhoHat2LoopSM::usage="";
DeltaRHat2LoopSM::usage="";
RhoHatTree::usage="";
InitGenerationOfDiagrams::usage="";
DeltaVBwave::usage="";
CreateDeltaVBContributions::usage="";

Begin["`Private`"];

GetBottomMass[] := ThresholdCorrections`GetParameter[TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]]];
GetTopMass[] := ThresholdCorrections`GetParameter[TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]]];

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
           Return[TreeMasses`ReplaceDependenciesReverse[massExpr[[1]]]];
          ];

(*extracts squared tree-level mass of Z boson before mixing with additional Z`*)
UnmixedZMass2[] :=
    Module[{ZMassMatrix, extraGaugeCouplings, submatrixList, submatrix, mass2Eigenvalues},
           ZMassMatrix = SARAH`MassMatrix[SARAH`VectorZ];
           Assert[MatrixQ[ZMassMatrix]];
           extraGaugeCouplings = Cases[SARAH`Gauge, x_ /; FreeQ[x, SARAH`hypercharge] && FreeQ[x, SARAH`left] && FreeQ[x, SARAH`color] :> x[[4]]];
           submatrixList = ZMassMatrix[[#, #]] & /@ Flatten[Table[{i, j}, {i, 1, Length[ZMassMatrix]}, {j, i + 1, Length[ZMassMatrix]}], 1];
           submatrix = Cases[submatrixList, x_ /; And @@ (FreeQ[x, #] & /@ extraGaugeCouplings)];
           If[Length[submatrix] != 1, Print["Error: Photon-Z mass matrix could not be identified"]; Return[0];];
           mass2Eigenvalues = Eigenvalues[submatrix];
           If[Length[mass2Eigenvalues] != 2 || !MemberQ[mass2Eigenvalues, 0], Print["Error: Determination of UnmixedZMass2 failed"]; Return[0];];
           Return[Select[mass2Eigenvalues, # =!= 0 &][[1]] /. Parameters`ApplyGUTNormalization[]];
          ];

(*extracts squared tree-level mass of W boson before mixing with additional W`*)
UnmixedWMass2[] :=
    Module[{WMassMatrix, extraGaugeCouplings, submatrixList, submatrix, mass2Eigenvalues},
           WMassMatrix = SARAH`MassMatrix[SARAH`VectorW];
           Assert[MatrixQ[WMassMatrix]];
           extraGaugeCouplings = Cases[SARAH`Gauge, x_ /; FreeQ[x, SARAH`hypercharge] && FreeQ[x, SARAH`left] && FreeQ[x, SARAH`color] :> x[[4]]];
           submatrixList = WMassMatrix[[#, #]] & /@ Flatten[Table[{i, j}, {i, 1, Length[WMassMatrix]}, {j, i + 1, Length[WMassMatrix]}], 1];
           submatrix = Cases[submatrixList, x_ /; And @@ (FreeQ[x, #] & /@ extraGaugeCouplings)];
           If[Length[submatrix] != 1, Print["Error: W mass matrix could not be identified"]; Return[0];];
           mass2Eigenvalues = Eigenvalues[submatrix];
           If[Length[DeleteDuplicates[mass2Eigenvalues]] != 1, Print["Error: Determination of UnmixedWMass2 failed"]; Return[0];];
           Return[mass2Eigenvalues[[1]] /. Parameters`ApplyGUTNormalization[]];
          ];

(*checks whether the Gell-Mann-Nishijima relation is valid*)
GellMannNishijimaRelationHolds[] :=
    Module[{photonMassMatrix, extraGaugeCouplings, submatrixIndices, BW3pos, photonEigenSystem, photonVector},
           If[FreeQ[SARAH`Gauge, SARAH`hypercharge] || FreeQ[SARAH`Gauge, SARAH`left], Print["Error: hypercharge or left gauge group does not exist. Please choose another method for the determination of the Weinberg angle."]; Return[False];];
           photonMassMatrix = SARAH`MassMatrix[SARAH`VectorP];
           Assert[MatrixQ[photonMassMatrix]];
           If[Length[photonMassMatrix] > 4, Print["Error: neutral vector boson mass matrix is too large to be diagonalized"]; Return[False];];
           extraGaugeCouplings = Cases[SARAH`Gauge, x_ /; FreeQ[x, SARAH`hypercharge] && FreeQ[x, SARAH`left] && FreeQ[x, SARAH`color] :> x[[4]]];
           submatrixIndices = Flatten[Table[{i, j}, {i, 1, Length[photonMassMatrix]}, {j, i + 1, Length[photonMassMatrix]}], 1];
           BW3pos = Flatten[Extract[submatrixIndices, Position[photonMassMatrix[[#, #]] & /@ submatrixIndices, x_ /; And @@ (FreeQ[x, #] & /@ extraGaugeCouplings), {1}, Heads -> False]]];
           If[Length[BW3pos] != 2, Print["Error: Photon-Z mass matrix could not be identified"]; Return[False];];
           photonEigenSystem = Eigensystem[photonMassMatrix];
           photonVector = Extract[photonEigenSystem[[2]], Position[photonEigenSystem[[1]], 0]];
           If[!MemberQ[Total[Abs[Part[#, Complement[Range[Length[photonMassMatrix]], BW3pos]] & /@ photonVector], {2}], 0], Print["Error: SM-like photon could not be identified. Please choose another method for the determination of the Weinberg angle."]; Return[False];];
           Return[True];
          ];

(*calculates rho_0 from SU(2)_L representations of the Higgs multipletts as in (16) from 0801.1345 [hep-ph]*)
RhoZero[] :=
    Module[{hyperchargePos, leftPos, vevlist},
           If[!GellMannNishijimaRelationHolds[], Print["Error: the Gell-Mann-Nishijima relation does not hold. Please choose another method for the determination of the Weinberg angle."]; Return[0];];
           hyperchargePos = Position[SARAH`Gauge, x_ /; !FreeQ[x, SARAH`hypercharge], {1}][[1, 1]];
           leftPos = Position[SARAH`Gauge, x_ /; !FreeQ[x, SARAH`left], {1}][[1, 1]];
           vevlist = SARAH`DEFINITION[SARAH`EWSB][SARAH`VEVs];
           (* extract isospin from SU(2)_left representation and its third component from Gell-Mann-Nishijima formula with given hypercharge and electric charge = 0 *)
           vevlist = vevlist /. {fieldname_Symbol, vevinfo_List, comp1_List, __} :> Flatten[{vevinfo Boole[ReleaseHold[SARAH`getElectricCharge[comp1[[1]]]] == 0], (SA`DimensionGG[fieldname, leftPos] - 1) / 2, -SA`ChargeGG[fieldname, hyperchargePos]}];
           If[!FreeQ[vevlist, None], Print["Error: determination of electric charge did not work"]; Return[0];];
           Return[Simplify[Plus @@ ((#[[3]]^2 - #[[4]]^2 + #[[3]]) Abs[#[[1]] #[[2]] Sqrt[2]]^2 & /@ vevlist) / Plus @@ (2 #[[4]]^2 Abs[#[[1]] #[[2]] Sqrt[2]]^2 & /@ vevlist),
                           Element[Alternatives @@ Cases[SARAH`DEFINITION[SARAH`EWSB][SARAH`VEVs][[All, 2, 1]], x_ /; Parameters`IsRealParameter[x]], Reals]]];
          ];

ExpressWeinbergAngleInTermsOfGaugeCouplings[] :=
    Module[{solution},
           Print["Expressing Weinberg angle in terms of model parameters ..."];
           solution = ArcCos[Sqrt[UnmixedWMass2[] / UnmixedZMass2[] / RhoZero[]]];
           Return[Simplify[solution]];
          ];

extPars={SINTHETAW, RHOHATRATIO, GFERMI, MW, MZ, MT, RHO2, DELTARHAT1LOOP, PIZZTMZ};
Do[Format[extPars[[i]],CForm]=Format[ToString[extPars[[i]]],OutputForm],{i,Length[extPars]}];

(*returns coefficients of Higgs-top-top vertices*)
HiggsTopVertices[higgsName_] :=
    Module[{indexRange, indexList, topQuark, higgsVertices, rule},
           If[FreeQ[TreeMasses`GetParticles[], higgsName] || TreeMasses`GetDimensionWithoutGoldstones[higgsName] == 0, Return[{}]];
           indexRange = TreeMasses`GetParticleIndices[higgsName][[All, 2]];
           If[indexRange === {}, indexRange = {1}];
           indexList = Flatten[Table @@ {Table[ToExpression["i" <> ToString[k]], {k, Length[indexRange]}], Sequence @@ Table[{ToExpression["i" <> ToString[k]], 1, indexRange[[k]]}, {k, Length[indexRange]}]}, Length[indexRange] - 1];
           topQuark = Level[TreeMasses`GetUpQuark[{3}], {Boole[ListQ[TreeMasses`GetUpQuark[{3}]]]}][[1]];
           higgsVertices = Vertices`StripGroupStructure[SARAH`Vertex[{bar[topQuark], topQuark, higgsName[#]}] & /@ indexList, SARAH`ctNr /@ Range[4]];
           rule = SARAH`sum[idx_, start_, stop_, expr_] :> Sum[expr, {idx, start, stop}];
           higgsVertices = Cases[higgsVertices, {{__, higgsField_}, {coeffPL_, SARAH`PL}, {coeffPR_, SARAH`PR}}
                                 /; ((coeffPL/I //. rule) * Susyno`LieGroups`conj[coeffPL/I //. rule] === (coeffPR/I //. rule) * Susyno`LieGroups`conj[coeffPR/I //. rule])
                                    && !TreeMasses`IsGoldstone[higgsField] :> {higgsField /. List -> Sequence, coeffPL/I}];
           Return[higgsVertices];
          ];

(*generalize Higgs dependent part of (C.5) and (C.6) in hep-ph/9606211 analogous to (C.9) and (C.10)*)
HiggsContributions2LoopSM[massMatrices_List] :=
    Module[{higgsMatrices, higgsVEVlist, higgsDep},
           If[!ValueQ[SARAH`VEVSM], Print["Error: SM like Higgs vev does not exist."]; Return[0];];
           higgsMatrices = Select[massMatrices, !FreeQ[#, SARAH`HiggsBoson] || !FreeQ[#, SARAH`PseudoScalar] &][[All, 1]];
           If[!(And @@ Parameters`IsRealParameter /@ Parameters`FindAllParameters[higgsMatrices]) || (Parameters`FSPhases /. Parameters`FindAllParametersClassified[higgsMatrices]) =!= {},
              Print["Warning: CP violation in Higgs sector"]; Return[0];];
           higgsVEVlist = Cases[Parameters`GetDependenceSPhenoRules[], RuleDelayed[SARAH`VEVSM, repr_] :> repr];
           If[higgsVEVlist === {}, higgsVEVlist = {SARAH`VEVSM}];
           higgsDep = Abs[#[[2]]]^2 RHO2[FlexibleSUSY`M[#[[1]]]/MT] &;
           Return[Simplify[3 (GFERMI MT higgsVEVlist[[1]] / (8 Pi^2 Sqrt[2]))^2 (Plus @@ (higgsDep /@ HiggsTopVertices[SARAH`HiggsBoson]) - Plus @@ (higgsDep /@ HiggsTopVertices[SARAH`PseudoScalar]))]];
          ];

(*formula according to (C.6) from hep-ph/9606211*)
DeltaRhoHat2LoopSM[massMatrices_List]:=
    Module[{gY, alphaDRbar, expr, result},
           gY = SARAH`hyperchargeCoupling FlexibleSUSY`GUTNormalization[SARAH`hyperchargeCoupling];
           alphaDRbar = gY^2 SARAH`leftCoupling^2 / (4 Pi (gY^2 + SARAH`leftCoupling^2));
           expr = (alphaDRbar SARAH`strongCoupling^2/(16 Pi^3 SINTHETAW^2)(-2.145 MT^2/MW^2 + 1.262 Log[MT/MZ] - 2.24 - 0.85 MZ^2/MT^2) + HiggsContributions2LoopSM[massMatrices]) / (1 + PIZZTMZ / MZ^2);
           result = Parameters`CreateLocalConstRefs[expr] <> "\n";
           result = result <> TreeMasses`ExpressionToString[expr, "deltaRhoHat2LoopSM"];
           Return[result];
          ];

(*formula according to (C.5) from hep-ph/9606211*)
DeltaRHat2LoopSM[massMatrices_List]:=
    Module[{gY, alphaDRbar, expr, result},
           gY = SARAH`hyperchargeCoupling FlexibleSUSY`GUTNormalization[SARAH`hyperchargeCoupling];
           alphaDRbar = gY^2 SARAH`leftCoupling^2 / (4 Pi (gY^2 + SARAH`leftCoupling^2));
           expr = alphaDRbar SARAH`strongCoupling^2/(16 Pi^3 SINTHETAW^2 (1 - SINTHETAW^2))(2.145 MT^2/MZ^2 + 0.575 Log[MT/MZ] - 0.224 - 0.144 MZ^2/MT^2) - HiggsContributions2LoopSM[massMatrices] (1 - DELTARHAT1LOOP) RHOHATRATIO;
           result = Parameters`CreateLocalConstRefs[expr] <> "\n";
           result = result <> TreeMasses`ExpressionToString[expr, "deltaRHat2LoopSM"];
           Return[result];
          ];

(*calculates tree-level value of rhohat parameter from umixed and mixed Z mass as well as RhoZero*)
RhoHatTree[]:=
    Module[{Zmass2unmixed, Zmass2mixed, expr, result},
           Zmass2unmixed = UnmixedZMass2[];
           Zmass2mixed = FindMassZ2[TreeMasses`GetUnmixedParticleMasses[] /. Parameters`ApplyGUTNormalization[]];
           expr = Simplify[RhoZero[] Zmass2unmixed / Zmass2mixed /. SARAH`Weinberg -> ExpressWeinbergAngleInTermsOfGaugeCouplings[], SARAH`hyperchargeCoupling > 0 && SARAH`leftCoupling > 0];
           result = Parameters`CreateLocalConstRefs[expr] <> "\n";
           result = result <> "rhohat_tree = ";
           result = result <> CConversion`RValueToCFormString[expr] <> ";";
           Return[result];
          ];


(*functions for creation of wave-function renormalization, vertex and box corrections:*)

InitGenerationOfDiagrams[eigenstates_:FlexibleSUSY`FSEigenstates] :=
    Module[{},
           SA`CurrentStates = eigenstates;
           SARAH`InitVertexCalculation[eigenstates, False];
           SARAH`ReadVertexList[eigenstates, False, False, True];
           SARAH`MakeCouplingLists;
          ];

ExcludeDiagrams[diagrs_List, excludeif_:(False &)] := Select[diagrs, !Or @@ (excludeif /@ (Cases[#, Rule[Internal[_], x_] :> x, Infinity])) &];

GenerateDiagramsWave[particle_] :=
    Module[{diagrs},
           diagrs = SARAH`InsFields[{{C[particle, SARAH`FieldToInsert[1], SARAH`AntiField[SARAH`FieldToInsert[2]]]},
                                    {SARAH`Internal[1] -> SARAH`FieldToInsert[1], SARAH`Internal[2] -> SARAH`FieldToInsert[2], SARAH`External[1] -> particle}}];
           (*add indices for later summation*)
           diagrs = diagrs /. (Rule[SARAH`Internal[i_], x_] /; TreeMasses`GetDimension[x] > 1) :> Rule[SARAH`Internal[i], x[{ToExpression["SARAH`gI" <> ToString[i]]}]];
           diagrs = diagrs /. (Rule[SARAH`External[i_], x_] /; TreeMasses`GetDimension[x] > 1) :> Rule[SARAH`External[i], x[{ToExpression["SARAH`gO" <> ToString[i]]}]];
           diagrs = ({{C[SARAH`External[1], SARAH`Internal[1], SARAH`AntiField[SARAH`Internal[2]]]} /. #[[2]], #[[2]]}) & /@ diagrs;
           Return[diagrs];
          ];

WaveResult[diagr_List, includeGoldstones_] :=
    Module[{coupl, intparticles, intfermion, intscalar, result, intpartwithindex},
           coupl = (diagr[[1, 1]] /. C[a__] -> SARAH`Cp[a])[SARAH`PL];
           intparticles = {SARAH`Internal[2], SARAH`Internal[1]} /. diagr[[2]];
           intfermion = Select[intparticles, TreeMasses`IsFermion][[1]];
           intscalar = Select[intparticles, TreeMasses`IsScalar][[1]];
           result = -coupl Susyno`LieGroups`conj[coupl] SARAH`B1[0, SARAH`Mass2[intfermion], SARAH`Mass2[intscalar]];
           intpartwithindex = Cases[intparticles /. {Susyno`LieGroups`conj -> Identity, SARAH`bar -> Identity}, _[{_}]];
           Do[result = SARAH`sum[intpartwithindex[[i, 1, 1]], If[includeGoldstones, 1, TreeMasses`GetDimensionStartSkippingGoldstones[intpartwithindex[[i]]]], TreeMasses`GetDimension[intpartwithindex[[i]]], result],
                 {i, Length[intpartwithindex]}];
           Return[result];
          ];

DeltaVBwave[includeGoldstones_:False] :=
    Module[{neutrinodiagrs, electrondiagrs, result},
           neutrinodiagrs = ExcludeDiagrams[GenerateDiagramsWave[SARAH`Neutrino], TreeMasses`IsVector];
           electrondiagrs = ExcludeDiagrams[GenerateDiagramsWave[SARAH`Electron], TreeMasses`IsVector];
           result = {WeinbergAngle`DeltaVB[{WeinbergAngle`wave, {SARAH`gO1}, SARAH`Neutrino}, 1/2 * Plus @@ (WaveResult[#, includeGoldstones] &) /@ neutrinodiagrs],
                     WeinbergAngle`DeltaVB[{WeinbergAngle`wave, {SARAH`gO1}, SARAH`Electron}, 1/2 * Plus @@ (WaveResult[#, includeGoldstones] &) /@ electrondiagrs]};
           Return[result];
          ];

AddIndices[{}] := "";

AddIndices[{ind_}] := "unsigned " <> CConversion`ToValidCSymbolString[ind];

AddIndices[{ind1_, ind2_}] := "unsigned " <> CConversion`ToValidCSymbolString[ind1] <> ", unsigned " <> CConversion`ToValidCSymbolString[ind2];

CreateContributionName[WeinbergAngle`DeltaVB[{type_, {___}}, _]] := "delta_vb_" <> CConversion`ToValidCSymbolString[type];

CreateContributionName[WeinbergAngle`DeltaVB[{type_, {___}, spec_}, _]] := "delta_vb_" <> CConversion`ToValidCSymbolString[type] <> "_" <> CConversion`ToValidCSymbolString[spec];

CreateContributionPrototype[deltaVBcontri_WeinbergAngle`DeltaVB] := CreateContributionName[deltaVBcontri] <> "(" <> AddIndices[deltaVBcontri[[1, 2]]] <> ") const";

(*based on CreateNPointFunction from SelfEnergies.m*)
CreateDeltaVBContribution[deltaVBcontri_WeinbergAngle`DeltaVB, vertexRules_List] :=
    Module[{expr, functionName, type, prototype, decl, body},
           expr = deltaVBcontri[[2]];
           functionName = CreateContributionPrototype[deltaVBcontri];
           type = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
           prototype = type <> " " <> functionName <> ";\n";
           decl = "\n" <> type <> " CLASSNAME::" <> functionName <> "\n{\n";
           body = type <> " result;\n\n" <>
                  CConversion`ExpandSums[Parameters`DecreaseIndexLiterals[Parameters`DecreaseSumIndices[expr], TreeMasses`GetParticles[]] /.
                                         vertexRules /.
                                         a_[List[i__]] :> a[i], "result"] <>
                  "\nreturn result * oneOver16PiSqr;";
           body = TextFormatting`IndentText[TextFormatting`WrapLines[body]];
           decl = decl <> body <> "\n}\n";
           Return[{prototype, decl}];
          ];

PrintDeltaVBContributionName[WeinbergAngle`DeltaVB[{WeinbergAngle`wave, {idx_}, part_}, _]] :=
    "deltaVB wave-function contribution for field " <> CConversion`ToValidCSymbolString[part] <> "[" <> CConversion`ToValidCSymbolString[idx] <> "]";

(*based on CreateNPointFunctions from SelfEnergies.m*)
CreateDeltaVBContributions[deltaVBcontris_List, vertexRules_List] :=
    Module[{relevantVertexRules, prototypes = "", defs = "", vertexFunctionNames = {}, p, d},
           Print["Converting vertex functions ..."];
           relevantVertexRules = Cases[vertexRules, r:(Rule[a_, b_] /; !FreeQ[deltaVBcontris, a])];
           {prototypes, defs, vertexFunctionNames} = SelfEnergies`CreateVertexExpressions[relevantVertexRules];
           Print["Generating C++ code for ..."];
           For[k = 1, k <= Length[deltaVBcontris], k++,
               Print["   ", PrintDeltaVBContributionName[deltaVBcontris[[k]]]];
               {p, d} = CreateDeltaVBContribution[deltaVBcontris[[k]], vertexFunctionNames];
               prototypes = prototypes <> p;
               defs = defs <> d;
              ];
           Return[{prototypes, defs}];
          ];

End[];

EndPackage[];

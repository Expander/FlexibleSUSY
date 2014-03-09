
BeginPackage["SelfEnergies`", {"SARAH`", "TextFormatting`", "CConversion`", "TreeMasses`", "Parameters`", "Vertices`"}];

FSSelfEnergy::usage="self-energy head";
Tadpole::usage="tadpole head";

GetField::usage="Returns field in self-energy or tadpole";

ConvertSarahSelfEnergies::usage="converts SARAH's self-energies to our
own format: SelfEnergies`FSSelfEnergy[particle, expression]";

ConvertSarahTadpoles::usage="converts SARAH's tadpoles to our own
format: SelfEnergies`Tadpole[particle, expression]";

CreateNPointFunctions::usage="creates C/C++ functions for the
given list of self-energies and tadpoles";

CreateSelfEnergyFunctionName::usage="creates self-energy function name
for a given field";

CreateHeavySelfEnergyFunctionName::usage="creates heavy self-energy
function name for a given field";

CreateHeavyRotatedSelfEnergyFunctionName::usage="creates heavy rotated
self-energy function name for a given field";

FillArrayWithOneLoopTadpoles::usage="add one-loop tadpoles to array"

FillArrayWithTwoLoopTadpoles::usage="add two-loop tadpoles to array"

CreateTwoLoopTadpolesMSSM::usage="Creates function prototypes and
definitions for two-loop tadpoles in the MSSM";

CreateTwoLoopSelfEnergiesMSSM::usage="Creates function prototypes and
definitions for two-loop Higgs self-energies in the MSSM";

Begin["`Private`"];

GetExpression[selfEnergy_SelfEnergies`FSSelfEnergy] :=
    selfEnergy[[2]];

GetExpression[selfEnergy_SelfEnergies`FSHeavySelfEnergy] :=
    selfEnergy[[2]];

GetExpression[selfEnergy_SelfEnergies`FSHeavyRotatedSelfEnergy] :=
    selfEnergy[[2]];

GetExpression[tadpole_SelfEnergies`Tadpole] :=
    tadpole[[2]];

GetField[selfEnergy_SelfEnergies`FSSelfEnergy] :=
    selfEnergy[[1]];

GetField[selfEnergy_SelfEnergies`FSHeavySelfEnergy] :=
    selfEnergy[[1]];

GetField[selfEnergy_SelfEnergies`FSHeavyRotatedSelfEnergy] :=
    selfEnergy[[1]];

GetField[tadpole_SelfEnergies`Tadpole] :=
    tadpole[[1]];

GetField[sym_] :=
    Module[{},
           Print["Error: GetField is not called with a SelfEnergies`FSSelfEnergy",
                 "or SelfEnergies`Tadpole head: ", sym];
           Quit[1];
          ];

RemoveSMParticles[SelfEnergies`FSSelfEnergy[p_,expr__], _] :=
    SelfEnergies`FSSelfEnergy[p,expr];

RemoveSMParticles[SelfEnergies`Tadpole[p_,expr__], _] :=
    SelfEnergies`Tadpole[p,expr];

ExprContainsNonOfTheseParticles[expr_, particles_List] :=
    And @@ (FreeQ[expr,#]& /@ particles);

RemoveSMParticles[head_[p_,expr_], removeGoldstones_:True, except_:{}] :=
    Module[{strippedExpr, susyParticles, a, goldstones, g, i},
           susyParticles = Join[TreeMasses`GetSusyParticles[], except];
           strippedExpr = expr /. ReplaceGhosts[];
           strippedExpr = strippedExpr //. {
               SARAH`Cp[a__  /; ExprContainsNonOfTheseParticles[{a},susyParticles]][_] -> 0,
               SARAH`Cp[a__  /; ExprContainsNonOfTheseParticles[{a},susyParticles]] -> 0
                                           };
           (* remove goldstone bosons *)
           If[removeGoldstones,
              goldstones = TreeMasses`GetSMGoldstoneBosons[];
              For[i = 1, i <= Length[goldstones], i++,
                  g = CConversion`GetHead[goldstones[[i]]];
                  strippedExpr = strippedExpr //.
                  SARAH`sum[idx_,_,endIdx_,expression_] /; !FreeQ[expression,g[{idx}]] :> SARAH`sum[idx,TreeMasses`GetDimensionStartSkippingGoldstones[g],endIdx,expression];
                 ];
             ];
           Return[head[p,strippedExpr]];
          ];

ReplaceUnrotatedFields[SelfEnergies`FSSelfEnergy[p_,expr_]] :=
    SelfEnergies`FSSelfEnergy[p,expr];

ReplaceUnrotatedFields[SelfEnergies`FSHeavySelfEnergy[p_,expr_]] :=
    SelfEnergies`FSHeavySelfEnergy[p,expr];

ReplaceUnrotatedFields[SelfEnergies`FSHeavyRotatedSelfEnergy[p_,expr__]] :=
    Module[{result},
           result = expr /. {
               SARAH`Cp[a__][l_] :> ReplaceUnrotatedFields[SARAH`Cp[a][l]],
               SARAH`Cp[a__]     :> ReplaceUnrotatedFields[SARAH`Cp[a]]
                            };
           Return[SelfEnergies`FSHeavyRotatedSelfEnergy[p,result]]
          ];

ConvertSarahTadpoles[DeleteLightFieldContrubtions[tadpoles_,_,_]] :=
    ConvertSarahTadpoles[tadpoles];

ConvertSarahTadpoles[tadpoles_List] :=
    Module[{result = {}, k, massESReplacements},
           massESReplacements = Join[
               Flatten[SARAH`diracSubBack1 /@ SARAH`NameOfStates],
               Flatten[SARAH`diracSubBack2 /@ SARAH`NameOfStates]];
           massESReplacements = Cases[massESReplacements, HoldPattern[Except[0] -> _]];
           result = (SelfEnergies`Tadpole @@ #)& /@ tadpoles /. massESReplacements;
           (* append mass eigenstate indices *)
           For[k = 1, k <= Length[result], k++,
               field = GetField[result[[k]]];
               If[GetDimension[field] > 1,
                  result[[k,1]] = field[SARAH`gO1];
                 ];
              ];
           Return[result /. SARAH`Mass2 -> FlexibleSUSY`M];
          ];

ConvertSarahSelfEnergies[selfEnergies_List] :=
    Module[{result = {}, k, field, fermionSE, left, right, scalar, expr,
            massESReplacements, heavySE},
           massESReplacements = Join[
               Flatten[SARAH`diracSubBack1 /@ SARAH`NameOfStates],
               Flatten[SARAH`diracSubBack2 /@ SARAH`NameOfStates]];
           massESReplacements = Cases[massESReplacements, HoldPattern[Except[0] -> _]];
           result = (SelfEnergies`FSSelfEnergy @@ #)& /@ selfEnergies /. massESReplacements;
           (* append mass eigenstate indices *)
           For[k = 1, k <= Length[result], k++,
               field = GetField[result[[k]]];
               If[GetDimension[field] > 1,
                  result[[k,1]] = field[SARAH`gO1, SARAH`gO2];
                 ];
              ];
           (* filter out all fermionic self-energies *)
           fermionSE = Cases[result, SelfEnergies`FSSelfEnergy[_, List[__]]];
           result = Select[result, (Head[GetExpression[#]] =!= List)&];
           (* decompose fermionic self-energies into L,R,S parts *)
           For[k = 1, k <= Length[fermionSE], k++,
               field = GetField[fermionSE[[k]]];
               expr = GetExpression[fermionSE[[k]]];
               If[Length[expr] != 3,
                  Print["Error: self-energy of ", field, " does not have 3 parts"];
                  Continue[];
                 ];
               AppendTo[result, SelfEnergies`FSSelfEnergy[field[1]       , expr[[1]]]];
               AppendTo[result, SelfEnergies`FSSelfEnergy[field[SARAH`PR], expr[[2]]]];
               AppendTo[result, SelfEnergies`FSSelfEnergy[field[SARAH`PL], expr[[3]]]];
              ];
           (* If the external field has dimension 1, remove it's
              indices.  For example in the Glu self-energy, terms
              appear of the form

                 Cp[Glu[{gO2}], conj[Sd[{gI1}]], Fd[{gI2}]]

              Since Glu has dimension 1, the C variables Glu is a
              double and must therefore not be accessed in the form
              Glu(gO2).
              *)
           For[k = 1, k <= Length[result], k++,
               field = GetHead[GetField[result[[k]]]];
               If[GetDimension[field] == 1,
                  result[[k,2]] = result[[k,2]] /. field[{__}] :> field;
                 ];
              ];
           (* Create W, Z self-energy with only SUSY particles in the loop *)
           heavySE = Cases[result, SelfEnergies`FSSelfEnergy[p:SARAH`VectorZ|SARAH`VectorW, expr__] :>
                           SelfEnergies`FSHeavySelfEnergy[p, expr]];
           result = Join[result, RemoveSMParticles /@ heavySE];
           (* Create Bottom, Tau self-energy with only SUSY
              particles and W and Z bosons in the loop *)
           heavySE = Cases[result, SelfEnergies`FSSelfEnergy[p:SARAH`BottomQuark[__][_]|SARAH`Electron[__][_], expr__] :>
                           SelfEnergies`FSHeavyRotatedSelfEnergy[p, expr]];
           result = Join[result,
                         ReplaceUnrotatedFields /@ (RemoveSMParticles[#,False,{SARAH`VectorZ,SARAH`VectorW}]& /@ heavySE)];
           (* Create Top self-energy with only SUSY
              particles and W, Z and photon bosons in the loop *)
           heavySE = Cases[result, SelfEnergies`FSSelfEnergy[p:SARAH`TopQuark[__][_], expr__] :>
                           SelfEnergies`FSHeavyRotatedSelfEnergy[p, expr]];
           result = Join[result,
                         ReplaceUnrotatedFields /@ (RemoveSMParticles[#,False,{SARAH`VectorZ,SARAH`VectorW,SARAH`VectorP}]& /@ heavySE)];
           Return[result /. SARAH`Mass -> FlexibleSUSY`M];
          ];

GetParticleIndices[Cp[a__]] := Flatten[Cases[{a}, List[__], Infinity]];

GetParticleIndices[Cp[a__][_]] := GetParticleIndices[Cp[a]];

CreateCouplingSymbol[coupling_] :=
    Module[{symbol, indices},
           indices = GetParticleIndices[coupling];
           symbol = ToValidCSymbol[coupling /. a_[List[__]] :> a];
           symbol[Sequence @@ indices]
          ];

(* creates a C++ function that calculates a coupling
 *
 * Return: {prototypes_String, definitions_String, rules_List}
 *
 * prototypes is a string that contains all coupling function
 * prototypes.  definitions is a string that contains all coupling
 * function definitions.
 *
 * rules is a list of replacement rules of the form
 * { Cp[bar[UCha[{gO2}]], VZ, Cha[{gI2}]][PR] :>
 *   CpbarUChaVZChaPR[gO2, gI2],
 *   ...
 * }
 *)
CreateCouplingFunction[coupling_, expr_] :=
    Module[{symbol, prototype = "", definition = "",
            indices = {}, body = "", cFunctionName = "", i,
            type, typeStr, initalValue},
           indices = GetParticleIndices[coupling];
           symbol = CreateCouplingSymbol[coupling];
           cFunctionName = ToValidCSymbolString[GetHead[symbol]];
           cFunctionName = cFunctionName <> "(";
           For[i = 1, i <= Length[indices], i++,
               If[i > 1, cFunctionName = cFunctionName <> ", ";];
               cFunctionName = cFunctionName <> "unsigned ";
               (* variable names must not be integers *)
               If[!IntegerQ[indices[[i]]] && !FreeQ[expr, indices[[i]]],
                  cFunctionName = cFunctionName <> ToValidCSymbolString[indices[[i]]];
                 ];
              ];
           cFunctionName = cFunctionName <> ")";
           If[Parameters`IsRealExpression[expr],
              type = CConversion`realScalarCType;    initalValue = " = 0.0";,
              type = CConversion`complexScalarCType; initalValue = "";];
           type = CConversion`ScalarType[type];
           typeStr = CConversion`CreateCType[type];
           prototype = typeStr <> " " <> cFunctionName <> " const;\n";
           definition = typeStr <> " CLASSNAME::" <> cFunctionName <> " const\n{\n";
           body = Parameters`CreateLocalConstRefsForInputParameters[expr, "LOCALINPUT"] <> "\n" <>
                  typeStr <> " result" <> initalValue <> ";\n\n";
           If[FreeQ[expr,SARAH`sum] && FreeQ[expr,SARAH`ThetaStep],
              body = body <> "result = " <>
                     RValueToCFormString[Simplify[DecreaseIndexLiterals[expr]]] <> ";\n";
              ,
              body = body <> ExpandSums[DecreaseIndexLiterals[DecreaseSumIdices[expr]],
                                        "result", type, initalValue];
             ];
           body = body <> "\nreturn result;\n";
           body = IndentText[WrapLines[body]];
           definition = definition <> body <> "}\n";
           Return[{prototype, definition,
                   RuleDelayed @@ {Vertices`ToCpPattern[coupling], symbol}}];
          ];

GetParticleList[Cp[a__]] := {a};

GetParticleList[Cp[a__][_]] := {a};

IsUnrotated[bar[field_]] := IsUnrotated[field];

IsUnrotated[Susyno`LieGroups`conj[field_]] := IsUnrotated[field];

IsUnrotated[field_[__]] := IsUnrotated[field];

IsUnrotated[field_Symbol] := StringTake[ToString[field],1] == "U";

ToRotatedField[field_Symbol] :=
    Symbol[StringReplace[ToString[field], StartOfString ~~ "U" ~~ rest_ :> rest]];

ToRotatedField[SARAH`bar[field_]] := SARAH`bar[ToRotatedField[field]];

ToRotatedField[Susyno`LieGroups`conj[field_]] := Susyno`LieGroups`conj[ToRotatedField[field]];

ToRotatedField[field_List] := ToRotatedField /@ field;

ToRotatedField[field_[indices__]] := ToRotatedField[field][indices];

GetUnrotatedFields[coupling_] := Select[GetParticleList[coupling], IsUnrotated];

ReplaceUnrotatedFields[SARAH`Cp[p__]] :=
    Cp[Sequence @@ ToRotatedField[{p}]];

ReplaceUnrotatedFields[SARAH`Cp[p__][lorentz_]] :=
    ReplaceUnrotatedFields[Cp[p]][lorentz];

CreateVertexExpressions[vertexRules_List] :=
    Module[{k, prototypes = "", defs = "", rules, coupling, expr,
            p, d, r},
           rules = Table[0, {Length[vertexRules]}];
           For[k = 1, k <= Length[vertexRules], k++,
               coupling = Vertices`ToCp[vertexRules[[k,1]]];
               expr = vertexRules[[k,2]];
               WriteString["stdout", "."];
               If[Mod[k, 50] == 0, WriteString["stdout","\n"]];
               {p,d,r} = CreateCouplingFunction[coupling, expr];
               prototypes = prototypes <> p;
               defs = defs <> d <> "\n";
               rules[[k]] = r;
              ];
           WriteString["stdout","\n"];
           Print["All vertices finished."];
           Return[{prototypes, defs, Flatten[rules]}];
          ];

ReplaceGhosts[states_:FlexibleSUSY`FSEigenstates] :=
    Module[{vectorBosons = {}, ghostStr, ghostSym, ghostCSym, ghosts = {}, k},
           vectorBosons = GetVectorBosons[states];
           For[k = 1, k <= Length[vectorBosons], k++,
               ghostStr = ToString[vectorBosons[[k]]];
               ghostStr = StringReplace[ghostStr, StartOfString ~~ "V" ~~ x_ :> "g" <> x];
               ghostSym = Symbol[ghostStr];
               ghostCSym = Symbol[ghostStr <> "C"];
               AppendTo[ghosts, ghostSym -> vectorBosons[[k]]];
               AppendTo[ghosts, ghostCSym -> vectorBosons[[k]]];
              ];
           Return[ghosts];
          ];

DeclareFieldIndices[field_Symbol] := "";

DeclareFieldIndices[field_[ind1_, ind2_]] :=
    ", unsigned " <> ToValidCSymbolString[ind1] <>
    ", unsigned " <> ToValidCSymbolString[ind2];

DeclareFieldIndices[field_[PL]] := DeclareFieldIndices[field];
DeclareFieldIndices[field_[PR]] := DeclareFieldIndices[field];
DeclareFieldIndices[field_[1]]  := DeclareFieldIndices[field];
DeclareFieldIndices[field_[ind_]] :=
    "unsigned " <> ToValidCSymbolString[ind];

CreateFunctionNamePrefix[field_[idx1_,idx2_]] := CreateFunctionNamePrefix[field];
CreateFunctionNamePrefix[field_[PL]]          := CreateFunctionNamePrefix[field] <> "_PL";
CreateFunctionNamePrefix[field_[PR]]          := CreateFunctionNamePrefix[field] <> "_PR";
CreateFunctionNamePrefix[field_[1]]           := CreateFunctionNamePrefix[field] <> "_1";
CreateFunctionNamePrefix[field_[idx_]]        := CreateFunctionNamePrefix[field];
CreateFunctionNamePrefix[field_]              := ToValidCSymbolString[field];

CreateSelfEnergyFunctionName[field_] :=
    "self_energy_" <> CreateFunctionNamePrefix[field];

CreateHeavySelfEnergyFunctionName[field_] :=
    "self_energy_" <> CreateFunctionNamePrefix[field] <> "_heavy";

CreateHeavyRotatedSelfEnergyFunctionName[field_] :=
    "self_energy_" <> CreateFunctionNamePrefix[field] <> "_heavy_rotated";

CreateTadpoleFunctionName[field_] :=
    "tadpole_" <> CreateFunctionNamePrefix[field];

CreateTwoLoopTadpoleFunctionName[field_] :=
    "tadpole_" <> CreateFunctionNamePrefix[field] <> "_2loop";

CreateTwoLoopSelfEnergyFunctionName[field_] :=
    "self_energy_" <> CreateFunctionNamePrefix[field] <> "_2loop";

CreateFunctionName[selfEnergy_SelfEnergies`FSSelfEnergy] :=
    CreateSelfEnergyFunctionName[GetField[selfEnergy]];

CreateFunctionName[selfEnergy_SelfEnergies`FSHeavySelfEnergy] :=
    CreateHeavySelfEnergyFunctionName[GetField[selfEnergy]];

CreateFunctionName[selfEnergy_SelfEnergies`FSHeavyRotatedSelfEnergy] :=
    CreateHeavyRotatedSelfEnergyFunctionName[GetField[selfEnergy]];

CreateFunctionName[tadpole_SelfEnergies`Tadpole] :=
    CreateTadpoleFunctionName[GetField[tadpole]];

CreateFunctionPrototype[selfEnergy_SelfEnergies`FSSelfEnergy] :=
    CreateFunctionName[selfEnergy] <>
    "(double p " <> DeclareFieldIndices[GetField[selfEnergy]] <> ") const";

CreateFunctionPrototype[selfEnergy_SelfEnergies`FSHeavySelfEnergy] :=
    CreateFunctionName[selfEnergy] <>
    "(double p " <> DeclareFieldIndices[GetField[selfEnergy]] <> ") const";

CreateFunctionPrototype[selfEnergy_SelfEnergies`FSHeavyRotatedSelfEnergy] :=
    CreateFunctionName[selfEnergy] <>
    "(double p " <> DeclareFieldIndices[GetField[selfEnergy]] <> ") const";

CreateFunctionPrototype[tadpole_SelfEnergies`Tadpole] :=
    CreateFunctionName[tadpole] <>
    "(" <> DeclareFieldIndices[GetField[tadpole]] <> ") const";

CreateNPointFunction[nPointFunction_, vertexRules_List] :=
    Module[{decl, expr, field, prototype, body, functionName},
           expr = GetExpression[nPointFunction];
           field = GetField[nPointFunction];
           functionName = CreateFunctionPrototype[nPointFunction];
           type = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
           prototype = type <> " " <> functionName <> ";\n";
           decl = "\n" <> type <> " CLASSNAME::" <> functionName <> "\n{\n";
           body = type <> " result;\n\n" <>
                  ExpandSums[DecreaseIndexLiterals[DecreaseSumIdices[expr], TreeMasses`GetParticles[]] /.
                             vertexRules /.
                             a_[List[i__]] :> a[i] /.
                             ReplaceGhosts[FlexibleSUSY`FSEigenstates] /.
                             C -> 1
                             ,"result"] <>
                  "\nreturn result * oneOver16PiSqr;";
           body = IndentText[WrapLines[body]];
           decl = decl <> body <> "\n}\n";
           Return[{prototype, decl}];
          ];

PrintNPointFunctionName[SelfEnergies`FSHeavySelfEnergy[field_,expr__]] :=
    "heavy " <> PrintNPointFunctionName[SelfEnergies`FSSelfEnergy[field,expr]];

PrintNPointFunctionName[SelfEnergies`FSHeavyRotatedSelfEnergy[field_,expr__]] :=
    "heavy, rotated " <> PrintNPointFunctionName[SelfEnergies`FSSelfEnergy[field,expr]];

PrintNPointFunctionName[SelfEnergies`FSSelfEnergy[field_[idx1_,idx2_][projector:(1|SARAH`PL|SARAH`PR)],__]] :=
    "self-energy Sigma^{" <> RValueToCFormString[field] <> "," <>
    RValueToCFormString[projector] <> "}_{" <>
    RValueToCFormString[idx1] <> "," <> RValueToCFormString[idx1] <> "}";

PrintNPointFunctionName[SelfEnergies`FSSelfEnergy[field_[projector:(1|SARAH`PL|SARAH`PR)],__]] :=
    "self-energy Sigma^{" <> RValueToCFormString[field] <> "," <>
    RValueToCFormString[projector] <> "}";

PrintNPointFunctionName[SelfEnergies`FSSelfEnergy[field_[idx1_,idx2_],__]] :=
    "self-energy Sigma^{" <> RValueToCFormString[field] <> "}_{" <>
    RValueToCFormString[idx1] <> "," <> RValueToCFormString[idx1] <> "}";

PrintNPointFunctionName[SelfEnergies`FSSelfEnergy[field_,__]] :=
    "self-energy Sigma^{" <> RValueToCFormString[field] <> "}";

PrintNPointFunctionName[SelfEnergies`Tadpole[field_[idx_],__]] :=
    "tadpole T^{" <> RValueToCFormString[field] <> "}_{" <> RValueToCFormString[idx] <> "}";

PrintNPointFunctionName[SelfEnergies`Tadpole[field_,__]] :=
    "tadpole T^{" <> RValueToCFormString[field] <> "}";

CreateNPointFunctions[nPointFunctions_List, vertexRules_List] :=
    Module[{prototypes = "", defs = "", vertexFunctionNames = {}, p, d},
           (* create coupling functions for all vertices in the list *)
           Print["Converting vertex functions ..."];
           {prototypes, defs, vertexFunctionNames} = CreateVertexExpressions[vertexRules];
           (* creating n-point functions *)
           Print["Generating C++ code for ..."];
           For[k = 1, k <= Length[nPointFunctions], k++,
               Print["   ", PrintNPointFunctionName[nPointFunctions[[k]]]];
               {p,d} = CreateNPointFunction[nPointFunctions[[k]], vertexFunctionNames];
               prototypes = prototypes <> p;
               defs = defs <> d;
              ];
           Return[{prototypes, defs}];
          ];

FillArrayWithOneLoopTadpoles[vevsAndFields_List, arrayName_String:"tadpole"] :=
    Module[{body = "", v, vev, field, idx, functionName},
           For[v = 1, v <= Length[vevsAndFields], v++,
               vev = vevsAndFields[[v,1]];
               field = vevsAndFields[[v,2]];
               idx = vevsAndFields[[v,3]];
               functionName = CreateTadpoleFunctionName[field];
               body = body <> arrayName <> "[" <> ToString[v-1] <> "] -= " <>
                      "Re(model->" <> functionName <>
                      "(" <> ToString[idx - 1] <> "));\n";
              ];
           Return[IndentText[body]];
          ];

FillArrayWithTwoLoopTadpoles[higgsBoson_, arrayName_String:"tadpole"] :=
    Module[{body, v, field, functionName},
           functionName = CreateTwoLoopTadpoleFunctionName[higgsBoson];
           body = "double two_loop_tadpole[2];\n" <>
                  "model->" <> functionName <>
                  "(two_loop_tadpole);\n";
           For[v = 1, v <= TreeMasses`GetDimension[higgsBoson], v++,
               body = body <> arrayName <> "[" <> ToString[v-1] <> "] -= " <>
                      "two_loop_tadpole[" <> ToString[v-1] <> "];\n";
              ];
           Return[IndentText[IndentText[body]]];
          ];

GetTwoLoopTadpoleCorrections[] :=
    Module[{body,
            g3Str, mtStr, mbStr, mtauStr,
            vev2Str, tanbStr, muStr, m3Str, mA0Str},
           mtStr   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`TopQuark][2]];
           mbStr   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`BottomQuark][2]];
           mtauStr = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Electron][2]];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vev2Str = CConversion`RValueToCFormString[SARAH`VEVSM1^2 + SARAH`VEVSM2^2];
           tanbStr = CConversion`RValueToCFormString[SARAH`VEVSM2 / SARAH`VEVSM1];
           muStr   = CConversion`RValueToCFormString[-Global`Mu];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = CConversion`RValueToCFormString[FlexibleSUSY`M[PseudoScalar][1]];
           body = "\
// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double mstau_1, mstau_2, theta_tau;
double msnu_1, msnu_2, theta_nu;

" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`TopQuark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`BottomQuark, "msb_1", "msb_2", "theta_b"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Electron, "mstau_1", "mstau_2", "theta_tau"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Neutrino, "msnu_1", "msnu_2", "theta_nu"] <>
";

double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
double msnusq = Sqr(msnu_2);
double sxt = Sin(theta_t), cxt = Cos(theta_t);
double sxb = Sin(theta_b), cxb = Cos(theta_b);
double sintau = Sin(theta_tau), costau = Cos(theta_tau);

double gs = " <> g3Str <> ";
double rmtsq = Sqr(" <> mtStr <> ");
double scalesq = Sqr(get_scale());
double vev2 = " <> vev2Str <> ";
double tanb = " <> tanbStr <> ";
double amu = " <> muStr <> ";
double mg = " <> m3Str <> ";
double mAsq = Sqr(" <> mA0Str <> ");
double cotbeta = 1.0 / tanb;
double rmbsq = Sqr(" <> mbStr <> ");
double rmtausq = Sqr(" <> mtauStr <> ");

double s1s = 0., s2s = 0., s1t = 0., s2t = 0.;
double s1b = 0., s2b = 0., s1tau = 0., s2tau = 0.;

ewsb2loop_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq,
           &amu, &tanb, &vev2, &gs, &s1s, &s2s);
ddstad_(&rmtsq, &rmbsq, &mAsq, &mst1sq, &mst2sq, &msb1sq, &msb2sq,
        &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, &s1t,
        &s2t);
ewsb2loop_(&rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq,
           &amu, &cotbeta, &vev2, &gs, &s2b, &s1b);
tausqtad_(&rmtausq, &mAsq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
          &costau, &scalesq, &amu, &tanb, &vev2, &s1tau, &s2tau);

if (!std::isnan(s1s * s1t * s1b * s1tau * s2s * s2t * s2b * s2tau)) {
   result[0] = (- s1s - s1t - s1b - s1tau) * " <> CConversion`ToValidCSymbolString[SARAH`VEVSM1] <> ";
   result[1] = (- s2s - s2t - s2b - s2tau) * " <> CConversion`ToValidCSymbolString[SARAH`VEVSM2] <> ";
} else {
   result[0] = 0.;
   result[1] = 0.;
}
";
           Return[body];
          ];

CreateTwoLoopTadpolesMSSM[higgsBoson_] :=
    Module[{prototype, function, functionName},
           functionName = CreateTwoLoopTadpoleFunctionName[higgsBoson];
           prototype = "void " <> functionName <> "(double result[2]) const;\n";
           body = GetTwoLoopTadpoleCorrections[];
           function = "void CLASSNAME::" <> functionName <>
                      "(double result[2]) const\n{\n" <> IndentText[body] <>
                      "\n}\n";
           Return[{prototype, function}];
          ];

GetTwoLoopSelfEnergyCorrections[] :=
    Module[{body,
            g3Str, mtStr, mbStr, mtauStr,
            vev2Str, vuStr, vdStr, tanbStr, muStr, m3Str, mA0Str},
           mtStr   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`TopQuark][2]];
           mbStr   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`BottomQuark][2]];
           mtauStr = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Electron][2]];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vev2Str = CConversion`RValueToCFormString[SARAH`VEVSM1^2 + SARAH`VEVSM2^2];
           vdStr   = CConversion`RValueToCFormString[SARAH`VEVSM1];
           vuStr   = CConversion`RValueToCFormString[SARAH`VEVSM2];
           tanbStr = CConversion`RValueToCFormString[SARAH`VEVSM2 / SARAH`VEVSM1];
           muStr   = CConversion`RValueToCFormString[-Global`Mu];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = CConversion`RValueToCFormString[FlexibleSUSY`M[PseudoScalar][1]];
           body = "\
// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double mstau_1, mstau_2, theta_tau;
double msnu_1, msnu_2, theta_nu;

" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`TopQuark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`BottomQuark, "msb_1", "msb_2", "theta_b"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Electron, "mstau_1", "mstau_2", "theta_tau"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Neutrino, "msnu_1", "msnu_2", "theta_nu"] <>
";

double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
double msnusq = Sqr(msnu_2);
double sxt = Sin(theta_t), cxt = Cos(theta_t);
double sxb = Sin(theta_b), cxb = Cos(theta_b);
double sintau = Sin(theta_tau), costau = Cos(theta_tau);

double gs = " <> g3Str <> ";
double rmtsq = Sqr(" <> mtStr <> ");
double scalesq = Sqr(get_scale());
double vev2 = " <> vev2Str <> ";
double tanb = " <> tanbStr <> ";
const double tanb2 = Sqr(tanb);
const double sinb = tanb / Sqrt(1. + tanb2);
const double cosb = 1. / Sqrt(1. + tanb2);
double amu = " <> muStr <> ";
double mg = " <> m3Str <> ";
double mAsq = Sqr(" <> mA0Str <> ");
double cotbeta = 1.0 / tanb;
double rmbsq = Sqr(" <> mbStr <> ");
double rmtausq = Sqr(" <> mtauStr <> ");
double fmasq = Abs(mAsq);

double s11s = 0., s22s = 0., s12s = 0.;
double s11b = 0., s12b = 0., s22b = 0.;
double s11tau = 0., s12tau = 0., s22tau = 0.;
double s11w = 0., s22w = 0., s12w = 0.;
int scheme = 0; // chooses DR-bar scheme from slavich et al

// two-loop Higgs corrections: alpha_s alpha_t, alpha_s alpha_b and
// alpha_b^2, alpha_t*2, alpha_b alpha_t
dszhiggs_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq, &amu,
          &tanb, &vev2, &gs, &scheme, &s11s, &s22s, &s12s);
dszhiggs_(&rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq, &amu,
          &cotbeta, &vev2, &gs, &scheme, &s22b, &s11b, &s12b);
ddshiggs_(&rmtsq, &rmbsq, &fmasq, &mst1sq, &mst2sq, &msb1sq, &msb2sq,
          &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, &s11w,
          &s12w, &s22w);
tausqhiggs_(&rmtausq, &fmasq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
            &costau, &scalesq, &amu, &tanb, &vev2, &scheme, &s11tau,
            &s22tau, &s12tau);

result[0] = - s11s - s11w - s11b - s11tau; // 1,1 element
result[1] = - s12s - s12w - s12b - s12tau; // 1,2 element
result[2] = - s22s - s22w - s22b - s22tau; // 2,2 element

// calculate dMA, which is the two loop correction to take the DRbar
// psuedoscalar mass ( = -2m3sq/sin(2beta)) to the pole mass (as in
// Eq. (8) of hep-ph/0305127)

double p2s = 0., p2w = 0., p2b = 0., p2tau = 0.;

dszodd_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq, &amu,
        &tanb, &vev2, &gs, &p2s);
dszodd_(&rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq, &amu,
        &cotbeta, &vev2, &gs, &p2b);
ddsodd_(&rmtsq, &rmbsq, &fmasq, &mst1sq, &mst2sq, &msb1sq, &msb2sq,
        &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, &p2w);
tausqodd_(&rmtausq, &fmasq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
          &costau, &scalesq, &amu, &tanb, &vev2, &p2tau);

const double dMA = p2s + p2w + p2b + p2tau;

// dMA contains two loop tadpoles, which we'll subtract
double tadpole[2];
" <> CreateTwoLoopTadpoleFunctionName[SARAH`HiggsBoson] <> "(tadpole);

result[0] += - dMA * Sqr(sinb) + tadpole[0] / " <> vdStr <> ";
result[1] += + dMA * sinb * cosb;
result[2] += - dMA * Sqr(cosb) + tadpole[1] / " <> vuStr <> ";
";
           Return[body];
          ];

CreateTwoLoopSelfEnergiesMSSM[higgsBoson_] :=
    Module[{prototype, function, functionName},
           functionName = CreateTwoLoopSelfEnergyFunctionName[higgsBoson];
           prototype = "void " <> functionName <> "(double result[3]) const;\n";
           body = GetTwoLoopSelfEnergyCorrections[];
           function = "void CLASSNAME::" <> functionName <>
                      "(double result[3]) const\n{\n" <> IndentText[body] <>
                      "\n}\n";
           Return[{prototype, function}];
          ];

End[];

EndPackage[];


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

CreateTwoLoopTadpolesNMSSM::usage="Creates function prototypes and
definitions for two-loop tadpoles in the NMSSM";

CreateTwoLoopSelfEnergiesSM::usage="Creates function prototypes and
definitions for two-loop Higgs self-energies in the SM";

CreateTwoLoopSelfEnergiesMSSM::usage="Creates function prototypes and
definitions for two-loop Higgs self-energies in the MSSM";

CreateTwoLoopSelfEnergiesNMSSM::usage="Creates function prototypes and
definitions for two-loop Higgs self-energies in the NMSSM";

CreateThreeLoopSelfEnergiesSplit::usage="Creates function prototypes and
definitions for three-loop Higgs self-energies in split-SUSY";

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

ExprContainsParticle[expr_, particle_] :=
    !FreeQ[expr,particle];

RemoveParticle[head_[p_,expr_], particle_] :=
    Module[{strippedExpr, a},
           strippedExpr = expr //. {
               SARAH`Cp[a__  /; ExprContainsParticle[{a},particle]][_] -> 0,
               SARAH`Cp[a__  /; ExprContainsParticle[{a},particle]] -> 0
                                   };
           Return[head[p,strippedExpr]];
          ];

RemoveSMParticles[SelfEnergies`FSSelfEnergy[p_,expr__], _] :=
    SelfEnergies`FSSelfEnergy[p,expr];

RemoveSMParticles[SelfEnergies`Tadpole[p_,expr__], _] :=
    SelfEnergies`Tadpole[p,expr];

ExprContainsNonOfTheseParticles[expr_, particles_List] :=
    And @@ (FreeQ[expr,#]& /@ particles);

RemoveSMParticles[head_[p_,expr_], removeGoldstones_:True, except_:{}] :=
    Module[{strippedExpr, keepParticles, a, goldstones,
            goldstonesWithoutIndex, g, i},
           keepParticles = Join[TreeMasses`GetSusyParticles[], except];
           goldstones = TreeMasses`GetSMGoldstoneBosons[];
           goldstonesWithoutIndex = goldstones /. a_[_] :> a;
           If[removeGoldstones,
              keepParticles = Complement[keepParticles, goldstonesWithoutIndex];,
              keepParticles = DeleteDuplicates[Join[keepParticles, goldstonesWithoutIndex]];
             ];
           strippedExpr = expr /. ReplaceGhosts[];
           strippedExpr = strippedExpr //. {
               SARAH`Cp[a__  /; ExprContainsNonOfTheseParticles[{a},keepParticles]][_] -> 0,
               SARAH`Cp[a__  /; ExprContainsNonOfTheseParticles[{a},keepParticles]] -> 0
                                           };
           (* remove goldstone bosons *)
           If[removeGoldstones,
              For[i = 1, i <= Length[goldstones], i++,
                  g = CConversion`GetHead[goldstones[[i]]];
                  strippedExpr = strippedExpr //.
                  SARAH`sum[idx_,_,endIdx_,expression_] /; !FreeQ[expression,g[{idx}]] :> SARAH`sum[idx,TreeMasses`GetDimensionStartSkippingSMGoldstones[g],endIdx,expression];
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
           Return[result /. SARAH`Mass -> FlexibleSUSY`M];
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
           heavySE = Cases[result, SelfEnergies`FSSelfEnergy[
               p:SARAH`BottomQuark[__][_]|SARAH`BottomQuark[_]|((particle_[__][_]|particle_[_]) /; TreeMasses`IsSMChargedLepton[particle]), expr__] :>
                           SelfEnergies`FSHeavyRotatedSelfEnergy[p, expr]];
           result = Join[result,
                         ReplaceUnrotatedFields /@ (RemoveSMParticles[#,False,{SARAH`VectorZ,SARAH`VectorW,SARAH`HiggsBoson}]& /@ heavySE)];
           (* Create rotated Top self-energy with only SUSY
              particles and W, Z and photon bosons in the loop *)
           heavySE = Cases[result, SelfEnergies`FSSelfEnergy[p:SARAH`TopQuark[__][_]|SARAH`TopQuark[_], expr__] :>
                           SelfEnergies`FSHeavyRotatedSelfEnergy[p, expr]];
           result = Join[result,
                         ReplaceUnrotatedFields /@ (RemoveParticle[#,SARAH`VectorG]& /@ heavySE)];
           (* Create unrotated Top self-energy with only SUSY
              particles and W, Z and photon bosons in the loop *)
           heavySE = Cases[result, SelfEnergies`FSSelfEnergy[p:SARAH`TopQuark[__][_]|SARAH`TopQuark[_], expr__] :>
                           SelfEnergies`FSHeavySelfEnergy[p, expr]];
           result = Join[result, RemoveParticle[#,SARAH`VectorG]& /@ heavySE];
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
              type = CConversion`ScalarType[CConversion`realScalarCType];    initalValue = " = 0.0";,
              type = CConversion`ScalarType[CConversion`complexScalarCType]; initalValue = "";];
           typeStr = CConversion`CreateCType[type];
           prototype = typeStr <> " " <> cFunctionName <> " const;\n";
           definition = typeStr <> " CLASSNAME::" <> cFunctionName <> " const\n{\n";
           body = Parameters`CreateLocalConstRefsForInputParameters[expr, "LOCALINPUT"] <> "\n" <>
                  typeStr <> " result" <> initalValue <> ";\n\n";
           body = body <> TreeMasses`ExpressionToString[expr, "result"];
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

CreateNLoopTadpoleFunctionName[field_, loop_] :=
    CreateTadpoleFunctionName[field] <> "_" <> ToString[loop] <> "loop";

CreateNLoopSelfEnergyFunctionName[field_, loop_] :=
    CreateSelfEnergyFunctionName[field] <> "_" <> ToString[loop] <> "loop";

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
                  ExpandSums[DecreaseIndexLiterals[DecreaseSumIndices[expr], TreeMasses`GetParticles[]] /.
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
    Module[{prototypes = "", defs = "", vertexFunctionNames = {}, p, d,
            relevantVertexRules},
           (* create coupling functions for all vertices in the list *)
           Print["Converting vertex functions ..."];
           (* extract vertex rules needed for the given nPointFunctions *)
           relevantVertexRules = Cases[vertexRules, r:(Rule[a_,b_] /; !FreeQ[nPointFunctions,a]) :> r];
           {prototypes, defs, vertexFunctionNames} = CreateVertexExpressions[relevantVertexRules];
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

FillArrayWithOneLoopTadpoles[higgsAndIdx_List, arrayName_String, sign_String:"-", struct_String:""] :=
    Module[{body = "", v, field, idx, head, functionName},
           For[v = 1, v <= Length[higgsAndIdx], v++,
               field = higgsAndIdx[[v,1]];
               idx = higgsAndIdx[[v,2]];
               head = CConversion`ToValidCSymbolString[higgsAndIdx[[v,3]]];
               functionName = CreateTadpoleFunctionName[field];
               If[TreeMasses`GetDimension[field] == 1,
                  body = body <> arrayName <> "[" <> ToString[v-1] <> "] " <> sign <> "= " <>
                         head <> "(" <> struct <> functionName <> "());\n";
                  ,
                  body = body <> arrayName <> "[" <> ToString[v-1] <> "] " <> sign <> "= " <>
                         head <> "(" <> struct <> functionName <>
                         "(" <> ToString[idx - 1] <> "));\n";
                 ];
              ];
           Return[IndentText[body]];
          ];

FillArrayWithTwoLoopTadpoles[higgsBoson_, arrayName_String, sign_String:"-", struct_String:""] :=
    Module[{body, v, field, functionName, dim, dimStr},
           functionName = CreateNLoopTadpoleFunctionName[higgsBoson,2];
           dim = GetDimension[higgsBoson];
           dimStr = ToString[dim];
           body = "const auto tadpole_2l(" <> struct <> functionName <> "());\n";
           For[v = 1, v <= dim, v++,
               body = body <> arrayName <> "[" <> ToString[v-1] <> "] " <> sign <> "= " <>
                      "tadpole_2l(" <> ToString[v-1] <> ");\n";
              ];
           Return[IndentText[IndentText[body]]];
          ];

AssertFieldDimension[field_, dim_, model_] :=
    Block[{fieldDim},
          fieldDim = GetDimension[field];
          If[fieldDim =!= dim,
             Print["Error: Cannot use ", model, " two-loop routines for ",
                   field, " (multiplet size ", fieldDim, ").  Multiplet size ",
                   dim, " required!"];
             Quit[1];
            ];
         ];

GetTwoLoopTadpoleCorrections[model_String /; model === "MSSM"] :=
    Module[{g3Str, mtStr, mbStr, mtauStr,
            mTop, mBot, mTau,
            vev2Str, tanbStr, muStr, m3Str, mA0Str},
           AssertFieldDimension[SARAH`HiggsBoson, 2, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mBot    = TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]];
           mTau    = TreeMasses`GetMass[TreeMasses`GetDownLepton[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mbStr   = CConversion`RValueToCFormString[mBot];
           mtauStr = CConversion`RValueToCFormString[mTau];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vev2Str = CConversion`RValueToCFormString[SARAH`VEVSM1^2 + SARAH`VEVSM2^2];
           tanbStr = CConversion`RValueToCFormString[SARAH`VEVSM2 / SARAH`VEVSM1];
           muStr   = CConversion`RValueToCFormString[-Parameters`GetEffectiveMu[]];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = TreeMasses`CallPseudoscalarHiggsMassGetterFunction[] <> "(0)";
"\
// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double mstau_1, mstau_2, theta_tau;
double msnu_1, msnu_2, theta_nu;

" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Selectron, "mstau_1", "mstau_2", "theta_tau"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Sneutrino, "msnu_1", "msnu_2", "theta_nu"] <>
";

const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
const double msnusq = Sqr(msnu_2);
const double sxt = Sin(theta_t), cxt = Cos(theta_t);
const double sxb = Sin(theta_b), cxb = Cos(theta_b);
const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
const double gs = " <> g3Str <> ";
const double rmtsq = Sqr(" <> mtStr <> ");
const double scalesq = Sqr(get_scale());
const double vev2 = " <> vev2Str <> ";
const double tanb = " <> tanbStr <> ";
const double amu = Re(" <> muStr <> ");
const double mg = " <> m3Str <> ";
const double mAsq = Sqr(" <> mA0Str <> ");
const double cotbeta = 1.0 / tanb;
const double rmbsq = Sqr(" <> mbStr <> ");
const double rmtausq = Sqr(" <> mtauStr <> ");

" <> GetTadpoleVectorCType[2] <> " tadpole_2l(" <> GetTadpoleVectorCType[2] <> "::Zero());

if (HIGGS_2LOOP_CORRECTION_AT_AS) {
   tadpole_2l += tadpole_higgs_2loop_at_as_mssm(
      rmtsq, mg, mst1sq, mst2sq, sxt, cxt, scalesq,
      amu, tanb, vev2, gs);
}

if (HIGGS_2LOOP_CORRECTION_AT_AT) {
   tadpole_2l += tadpole_higgs_2loop_at_at_mssm(
      rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq, msb2sq,
      sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
}

if (HIGGS_2LOOP_CORRECTION_AB_AS) {
   tadpole_2l += tadpole_higgs_2loop_ab_as_mssm(
      rmbsq, mg, msb1sq, msb2sq, sxb, cxb, scalesq,
      amu, cotbeta, vev2, gs);
}

if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
   tadpole_2l += tadpole_higgs_2loop_atau_atau_mssm(
      rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
      costau, scalesq, amu, tanb, vev2);
}

tadpole_2l(0) *= " <> CConversion`ToValidCSymbolString[SARAH`VEVSM1] <> ";
tadpole_2l(1) *= " <> CConversion`ToValidCSymbolString[SARAH`VEVSM2] <> ";

if (!IsFinite(tadpole_2l))
   tadpole_2l.setZero();

return tadpole_2l;"
          ];

GetTwoLoopTadpoleCorrections[model_String /; model === "NMSSM"] :=
    Module[{mTop, mBot, mTau,
            g3Str, mtStr, mbStr, mtauStr,
            vev2Str, svevStr, tanbStr, muStr, m3Str, mA0Str},
           AssertFieldDimension[SARAH`HiggsBoson, 3, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mBot    = TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]];
           mTau    = TreeMasses`GetMass[TreeMasses`GetDownLepton[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mbStr   = CConversion`RValueToCFormString[mBot];
           mtauStr = CConversion`RValueToCFormString[mTau];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vev2Str = CConversion`RValueToCFormString[SARAH`VEVSM1^2 + SARAH`VEVSM2^2];
           tanbStr = CConversion`RValueToCFormString[SARAH`VEVSM2 / SARAH`VEVSM1];
           muStr   = CConversion`RValueToCFormString[-Parameters`GetEffectiveMu[]];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = CConversion`RValueToCFormString[Parameters`GetEffectiveMASqr[]];
           svevStr = CConversion`RValueToCFormString[Parameters`GetParameterFromDescription["Singlet-VEV"]];
"\
// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double mstau_1, mstau_2, theta_tau;
double msnu_1, msnu_2, theta_nu;

" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Selectron, "mstau_1", "mstau_2", "theta_tau"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Sneutrino, "msnu_1", "msnu_2", "theta_nu"] <>
";

const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
const double msnusq = Sqr(msnu_2);
const double sxt = Sin(theta_t), cxt = Cos(theta_t);
const double sxb = Sin(theta_b), cxb = Cos(theta_b);
const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
const double gs = " <> g3Str <> ";
const double rmtsq = Sqr(" <> mtStr <> ");
const double scalesq = Sqr(get_scale());
const double vev2 = " <> vev2Str <> ";
const double tanb = " <> tanbStr <> ";
const double amu = Re(" <> muStr <> ");
const double mg = " <> m3Str <> ";
const double mAsq = " <> mA0Str <> ";
const double cotbeta = 1.0 / tanb;
const double rmbsq = Sqr(" <> mbStr <> ");
const double rmtausq = Sqr(" <> mtauStr <> ");
const double svevS = " <> svevStr <> " / Sqrt(2.0);

" <> GetTadpoleVectorCType[3] <> " tadpole_2l(" <> GetTadpoleVectorCType[3] <> "::Zero());

if (HIGGS_2LOOP_CORRECTION_AT_AS) {
   tadpole_2l += tadpole_higgs_2loop_at_as_nmssm(
      rmtsq, mg, mst1sq, mst2sq, sxt, cxt, scalesq,
      amu, tanb, vev2, gs, svevS);
}

if (HIGGS_2LOOP_CORRECTION_AT_AT) {
   tadpole_2l.head<2>() += tadpole_higgs_2loop_at_at_mssm(
      rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq, msb2sq,
      sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
}

if (HIGGS_2LOOP_CORRECTION_AB_AS) {
   tadpole_2l += tadpole_higgs_2loop_ab_as_nmssm(
      rmbsq, mg, msb1sq, msb2sq, sxb, cxb, scalesq,
      amu, cotbeta, vev2, gs, svevS);
}

if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
   tadpole_2l.head<2>() += tadpole_higgs_2loop_atau_atau_mssm(
      rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
      costau, scalesq, amu, tanb, vev2);
}

tadpole_2l(0) *= " <> CConversion`ToValidCSymbolString[SARAH`VEVSM1] <> ";
tadpole_2l(1) *= " <> CConversion`ToValidCSymbolString[SARAH`VEVSM2] <> ";
tadpole_2l(2) *= " <>  svevStr <> ";

if (!IsFinite(tadpole_2l))
   tadpole_2l.setZero();

return tadpole_2l;"
          ];

GetTwoLoopTadpoleCorrections[model_] :=
    Module[{},
           Print["Error: two-loop tadpoles for ", model, " not available"];
           ""
          ];

GetTadpoleVectorType[dim_] :=
    CConversion`VectorType[CConversion`realScalarCType, dim];

GetTadpoleVectorCType[dim_] :=
    CConversion`CreateCType[GetTadpoleVectorType[dim]];

CreateTwoLoopTadpoles[higgsBoson_, model_String] :=
    Module[{prototype, function, functionName, dim, dimStr, cType},
           dim = GetDimension[higgsBoson];
           dimStr = ToString[dim];
           functionName = CreateNLoopTadpoleFunctionName[higgsBoson,2];
           cType = GetTadpoleVectorCType[dim];
           prototype = cType <> " " <> functionName <> "() const;\n";
           body = GetTwoLoopTadpoleCorrections[model];
           function = cType <> " CLASSNAME::" <> functionName <> "() const\n{\n" <>
                      IndentText[body] <> "\n}\n";
           Return[{prototype, function}];
          ];

CreateTwoLoopTadpolesMSSM[higgsBoson_] :=
    CreateTwoLoopTadpoles[higgsBoson, "MSSM"];

CreateTwoLoopTadpolesNMSSM[higgsBoson_] :=
    CreateTwoLoopTadpoles[higgsBoson, "NMSSM"];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`HiggsBoson,
                              model_String /; model === "SM", 2] :=
    Module[{mTop, mtStr, yt, ytStr, g3Str},
           AssertFieldDimension[particle, 1, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           yt      = Parameters`GetThirdGeneration[SARAH`UpYukawa];
           ytStr   = CConversion`RValueToCFormString[yt];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
"\
const double mt = " <> mtStr <> ";
const double yt = " <> ytStr <> ";
const double gs = " <> g3Str <> ";
const double scale = get_scale();
double self_energy = 0.;

if (HIGGS_2LOOP_CORRECTION_AT_AT) {
   self_energy += self_energy_higgs_2loop_at_at_sm(scale, mt, yt);
}

if (HIGGS_2LOOP_CORRECTION_AT_AS) {
   self_energy += self_energy_higgs_2loop_at_as_sm(scale, mt, yt, gs);
}

return self_energy;"
          ];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`HiggsBoson,
                              model_String /; model === "Split", 3] :=
    Module[{mTop, mGluino, mtStr, mgStr, yt, ytStr, g3Str},
           AssertFieldDimension[particle, 1, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mGluino = TreeMasses`GetMass[Parameters`GetParticleFromDescription["Gluino"]];
           mgStr   = CConversion`RValueToCFormString[mGluino];
           yt      = Parameters`GetThirdGeneration[SARAH`UpYukawa];
           ytStr   = CConversion`RValueToCFormString[yt];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
"\
const double mt = " <> mtStr <> ";
const double mg = " <> mgStr <> ";
const double yt = " <> ytStr <> ";
const double gs = " <> g3Str <> ";
const double scale = get_scale();
double self_energy = 0.;

if (HIGGS_3LOOP_CORRECTION_AT_AS_AS) {
   self_energy += self_energy_higgs_3loop_gluino_split(scale, mt, yt, gs, mg);
}

return self_energy;"
          ];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`HiggsBoson,
                              model_String /; model === "MSSM", 2] :=
    Module[{g3Str, mtStr, mbStr, mtauStr,
            mTop, mBot, mTau,
            vev2Str, vuStr, vdStr, tanbStr, muStr, m3Str, mA0Str},
           AssertFieldDimension[particle, 2, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mBot    = TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]];
           mTau    = TreeMasses`GetMass[TreeMasses`GetDownLepton[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mbStr   = CConversion`RValueToCFormString[mBot];
           mtauStr = CConversion`RValueToCFormString[mTau];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vev2Str = CConversion`RValueToCFormString[SARAH`VEVSM1^2 + SARAH`VEVSM2^2];
           vdStr   = CConversion`RValueToCFormString[SARAH`VEVSM1];
           vuStr   = CConversion`RValueToCFormString[SARAH`VEVSM2];
           tanbStr = CConversion`RValueToCFormString[SARAH`VEVSM2 / SARAH`VEVSM1];
           muStr   = CConversion`RValueToCFormString[-Parameters`GetEffectiveMu[]];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = TreeMasses`CallPseudoscalarHiggsMassGetterFunction[] <> "(0)";
"\
// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double mstau_1, mstau_2, theta_tau;
double msnu_1, msnu_2, theta_nu;

" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Selectron, "mstau_1", "mstau_2", "theta_tau"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Sneutrino, "msnu_1", "msnu_2", "theta_nu"] <>
";

const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
const double msnusq = Sqr(msnu_2);
const double sxt = Sin(theta_t), cxt = Cos(theta_t);
const double sxb = Sin(theta_b), cxb = Cos(theta_b);
const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
const double gs = " <> g3Str <> ";
const double rmtsq = Sqr(" <> mtStr <> ");
const double scalesq = Sqr(get_scale());
const double vev2 = " <> vev2Str <> ";
const double tanb = " <> tanbStr <> ";
const double amu = Re(" <> muStr <> ");
const double mg = " <> m3Str <> ";
const double mAsq = Sqr(" <> mA0Str <> ");
const double cotbeta = 1.0 / tanb;
const double rmbsq = Sqr(" <> mbStr <> ");
const double rmtausq = Sqr(" <> mtauStr <> ");

" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`HiggsBoson]] <> " self_energy_2l(" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`HiggsBoson]] <> "::Zero());

if (HIGGS_2LOOP_CORRECTION_AT_AS) {
   self_energy_2l += self_energy_higgs_2loop_at_as_mssm(
      rmtsq, mg, mst1sq, mst2sq, sxt, cxt, scalesq, amu,
      tanb, vev2, gs);
}

if (HIGGS_2LOOP_CORRECTION_AB_AS) {
   self_energy_2l += self_energy_higgs_2loop_ab_as_mssm(
      rmbsq, mg, msb1sq, msb2sq, sxb, cxb, scalesq, amu,
      cotbeta, vev2, gs);
}

if (HIGGS_2LOOP_CORRECTION_AT_AT) {
   self_energy_2l += self_energy_higgs_2loop_at_at_mssm(
      rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq, msb2sq,
      sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
}

if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
   self_energy_2l += self_energy_higgs_2loop_atau_atau_mssm(
      rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
      costau, scalesq, amu, tanb, vev2);
}

return self_energy_2l;"
          ];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`PseudoScalar,
                              model_String /; model === "MSSM", 2] :=
    Module[{g3Str, mtStr, mbStr, mtauStr,
            mTop, mBot, mTau,
            vev2Str, vuStr, vdStr, tanbStr, muStr, m3Str, mA0Str},
           AssertFieldDimension[particle, 2, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mBot    = TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]];
           mTau    = TreeMasses`GetMass[TreeMasses`GetDownLepton[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mbStr   = CConversion`RValueToCFormString[mBot];
           mtauStr = CConversion`RValueToCFormString[mTau];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vev2Str = CConversion`RValueToCFormString[SARAH`VEVSM1^2 + SARAH`VEVSM2^2];
           vdStr   = CConversion`RValueToCFormString[SARAH`VEVSM1];
           vuStr   = CConversion`RValueToCFormString[SARAH`VEVSM2];
           tanbStr = CConversion`RValueToCFormString[SARAH`VEVSM2 / SARAH`VEVSM1];
           muStr   = CConversion`RValueToCFormString[-Parameters`GetEffectiveMu[]];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = TreeMasses`CallPseudoscalarHiggsMassGetterFunction[] <> "(0)";
"\
// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double mstau_1, mstau_2, theta_tau;
double msnu_1, msnu_2, theta_nu;

" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Selectron, "mstau_1", "mstau_2", "theta_tau"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Sneutrino, "msnu_1", "msnu_2", "theta_nu"] <>
";

const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
const double msnusq = Sqr(msnu_2);
const double sxt = Sin(theta_t), cxt = Cos(theta_t);
const double sxb = Sin(theta_b), cxb = Cos(theta_b);
const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
const double gs = " <> g3Str <> ";
const double rmtsq = Sqr(" <> mtStr <> ");
const double scalesq = Sqr(get_scale());
const double vev2 = " <> vev2Str <> ";
const double tanb = " <> tanbStr <> ";
const double amu = Re(" <> muStr <> ");
const double mg = " <> m3Str <> ";
const double mAsq = Sqr(" <> mA0Str <> ");
const double cotbeta = 1.0 / tanb;
const double rmbsq = Sqr(" <> mbStr <> ");
const double rmtausq = Sqr(" <> mtauStr <> ");

" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`PseudoScalar]] <> " self_energy_2l(" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`PseudoScalar]] <> "::Zero());

if (HIGGS_2LOOP_CORRECTION_AT_AS) {
   self_energy_2l += self_energy_pseudoscalar_2loop_at_as_mssm(
      rmtsq, mg, mst1sq, mst2sq, sxt, cxt, scalesq, amu,
      tanb, vev2, gs);
}

if (HIGGS_2LOOP_CORRECTION_AB_AS) {
   self_energy_2l += self_energy_pseudoscalar_2loop_ab_as_mssm(
      rmbsq, mg, msb1sq, msb2sq, sxb, cxb, scalesq, amu,
      cotbeta, vev2, gs);
}

if (HIGGS_2LOOP_CORRECTION_AT_AT) {
   self_energy_2l += self_energy_pseudoscalar_2loop_at_at_mssm(
      rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq, msb2sq,
      sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
}

if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
   self_energy_2l += self_energy_pseudoscalar_2loop_atau_atau_mssm(
      rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
      costau, scalesq, amu, tanb, vev2);
}

return self_energy_2l;"
          ];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`HiggsBoson,
                              model_String /; model === "NMSSM", 2] :=
    Module[{g3Str, mtStr, mbStr, mtauStr,
            mTop, mBot, mTau,
            vev2Str, vuStr, vdStr, vsStr, tanbStr, muStr, m3Str, mA0Str,
            lambdaStr},
           AssertFieldDimension[particle, 3, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mBot    = TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]];
           mTau    = TreeMasses`GetMass[TreeMasses`GetDownLepton[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mbStr   = CConversion`RValueToCFormString[mBot];
           mtauStr = CConversion`RValueToCFormString[mTau];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vev2Str = CConversion`RValueToCFormString[SARAH`VEVSM1^2 + SARAH`VEVSM2^2];
           vdStr   = CConversion`RValueToCFormString[SARAH`VEVSM1];
           vuStr   = CConversion`RValueToCFormString[SARAH`VEVSM2];
           tanbStr = CConversion`RValueToCFormString[SARAH`VEVSM2 / SARAH`VEVSM1];
           muStr   = CConversion`RValueToCFormString[-Parameters`GetEffectiveMu[]];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = CConversion`RValueToCFormString[Parameters`GetEffectiveMASqr[]];
           vsStr   = CConversion`RValueToCFormString[Parameters`GetParameterFromDescription["Singlet-VEV"]];
           lambdaStr = CConversion`RValueToCFormString[Parameters`GetParameterFromDescription["Singlet-Higgs-Interaction"]];
"\
// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double mstau_1, mstau_2, theta_tau;
double msnu_1, msnu_2, theta_nu;

" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Selectron, "mstau_1", "mstau_2", "theta_tau"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Sneutrino, "msnu_1", "msnu_2", "theta_nu"] <>
";

const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
const double msnusq = Sqr(msnu_2);
const double sxt = Sin(theta_t), cxt = Cos(theta_t);
const double sxb = Sin(theta_b), cxb = Cos(theta_b);
const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
const double gs = " <> g3Str <> ";
const double as = Sqr(gs) / (4.0 * Pi);
const double rmt = " <> mtStr <> ";
const double rmtsq = Sqr(rmt);
const double scalesq = Sqr(get_scale());
const double vev2 = " <> vev2Str <> ";
const double vev = Sqrt(" <> vev2Str <> ");
const double tanb = " <> tanbStr <> ";
const double amu = Re(" <> muStr <> ");
const double mg = " <> m3Str <> ";
const double mAsq = " <> mA0Str <> ";
const double cotb = 1.0 / tanb;
const double rmb = " <> mbStr <> ";
const double rmbsq = Sqr(rmb);
const double rmtausq = Sqr(" <> mtauStr <> ");
const double lamS = Re(" <> lambdaStr <> ");
static const double root2 = Sqrt(2.0);
const double vevS =  vev / root2;
const double svevS = " <> vsStr <> " / root2;

" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`HiggsBoson]] <> " self_energy_2l(" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`HiggsBoson]] <> "::Zero());

if (HIGGS_2LOOP_CORRECTION_AT_AS) {
   self_energy_2l += self_energy_higgs_2loop_at_as_nmssm(
      rmt, mg, mst1sq, mst2sq, sxt, cxt,
      scalesq, tanb, vevS, lamS, svevS, as, amu);
}

if (HIGGS_2LOOP_CORRECTION_AB_AS) {
   self_energy_2l += self_energy_higgs_2loop_ab_as_nmssm(
      rmb, mg, msb1sq, msb2sq, sxb, cxb,
      scalesq, cotb, vevS, lamS, svevS, as, amu);
}

// Corrections as in MSSM, not corrected for NMSSM,
// should be OK for MSSM states when S state is close to decoupled

if (HIGGS_2LOOP_CORRECTION_AT_AT) {
   self_energy_2l.topLeftCorner<2,2>() += self_energy_higgs_2loop_at_at_mssm(
      rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq,
      msb2sq, sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
}

if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
   self_energy_2l.topLeftCorner<2,2>() += self_energy_higgs_2loop_atau_atau_mssm(
      rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
      costau, scalesq, amu, tanb, vev2);
}

return self_energy_2l;"
          ];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`PseudoScalar,
                              model_String /; model === "NMSSM", 2] :=
    Module[{g3Str, mtStr, mbStr, mtauStr,
            mTop, mBot, mTau,
            vev2Str, vuStr, vdStr, vsStr, tanbStr, muStr, m3Str, mA0Str,
            lambdaStr},
           AssertFieldDimension[particle, 3, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mBot    = TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]];
           mTau    = TreeMasses`GetMass[TreeMasses`GetDownLepton[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mbStr   = CConversion`RValueToCFormString[mBot];
           mtauStr = CConversion`RValueToCFormString[mTau];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vev2Str = CConversion`RValueToCFormString[SARAH`VEVSM1^2 + SARAH`VEVSM2^2];
           vdStr   = CConversion`RValueToCFormString[SARAH`VEVSM1];
           vuStr   = CConversion`RValueToCFormString[SARAH`VEVSM2];
           tanbStr = CConversion`RValueToCFormString[SARAH`VEVSM2 / SARAH`VEVSM1];
           muStr   = CConversion`RValueToCFormString[-Parameters`GetEffectiveMu[]];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = CConversion`RValueToCFormString[Parameters`GetEffectiveMASqr[]];
           vsStr   = CConversion`RValueToCFormString[Parameters`GetParameterFromDescription["Singlet-VEV"]];
           lambdaStr = CConversion`RValueToCFormString[Parameters`GetParameterFromDescription["Singlet-Higgs-Interaction"]];
"\
// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double mstau_1, mstau_2, theta_tau;
double msnu_1, msnu_2, theta_nu;

" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Selectron, "mstau_1", "mstau_2", "theta_tau"] <>
";
" <> TreeMasses`CallThirdGenerationHelperFunctionName[SARAH`Sneutrino, "msnu_1", "msnu_2", "theta_nu"] <>
";

const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
const double msnusq = Sqr(msnu_2);
const double sxt = Sin(theta_t), cxt = Cos(theta_t);
const double sxb = Sin(theta_b), cxb = Cos(theta_b);
const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
const double gs = " <> g3Str <> ";
const double as = Sqr(gs) / (4.0 * Pi);
const double rmt = " <> mtStr <> ";
const double rmtsq = Sqr(rmt);
const double scalesq = Sqr(get_scale());
const double vev2 = " <> vev2Str <> ";
const double vev = Sqrt(" <> vev2Str <> ");
const double tanb = " <> tanbStr <> ";
const double amu = Re(" <> muStr <> ");
const double mg = " <> m3Str <> ";
const double mAsq = " <> mA0Str <> ";
const double cotb = 1.0 / tanb;
const double rmb = " <> mbStr <> ";
const double rmbsq = Sqr(rmb);
const double rmtausq = Sqr(" <> mtauStr <> ");
const double lamS = Re(" <> lambdaStr <> ");
static const double root2 = Sqrt(2.0);
const double vevS =  vev / root2;
const double svevS = " <> vsStr <> " / root2;

" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`PseudoScalar]] <> " self_energy_2l(" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`PseudoScalar]] <> "::Zero());

if (HIGGS_2LOOP_CORRECTION_AT_AS) {
   self_energy_2l += self_energy_pseudoscalar_2loop_at_as_nmssm(
      rmt, mg, mst1sq, mst2sq, sxt, cxt,
      scalesq, tanb, vevS, lamS, svevS, as, amu);
}

if (HIGGS_2LOOP_CORRECTION_AB_AS) {
   self_energy_2l += self_energy_pseudoscalar_2loop_ab_as_nmssm(
      rmb, mg, msb1sq, msb2sq, sxb, cxb,
      scalesq, cotb, vevS, lamS, svevS, as, amu);
}

// Corrections as in MSSM, not corrected for NMSSM,
// should be OK for MSSM states when S state is close to decoupled

if (HIGGS_2LOOP_CORRECTION_AT_AT) {
   self_energy_2l.topLeftCorner<2,2>() += self_energy_pseudoscalar_2loop_at_at_mssm(
      rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq, msb2sq,
      sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
}

if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
   self_energy_2l.topLeftCorner<2,2>() += self_energy_pseudoscalar_2loop_atau_atau_mssm(
      rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
      costau, scalesq, amu, tanb, vev2);
}

return self_energy_2l;"
          ];

GetNLoopSelfEnergyCorrections[particle_, model_, loop_] :=
    Module[{},
           Print["Error: ", loop,"-loop self-energy corrections not available for ",
                 particle, " in the ", model];
           ""
          ];

CreateNLoopSelfEnergy[particle_, model_String, loop_] :=
    Module[{prototype, function, functionName, dim, dimStr, cType},
           dim = Parameters`NumberOfIndependentEntriesOfSymmetricMatrix[GetDimension[particle]];
           dimStr = ToString[dim];
           functionName = CreateNLoopSelfEnergyFunctionName[particle,loop];
           cType = CConversion`CreateCType[TreeMasses`GetMassMatrixType[particle]];
           prototype = cType <> " " <> functionName <> "() const;\n";
           body = GetNLoopSelfEnergyCorrections[particle, model, loop];
           function = cType <> " CLASSNAME::" <> functionName <> "() const\n{\n" <>
                      IndentText[body] <> "\n}\n";
           Return[{prototype, function}];
          ];

CreateNLoopSelfEnergies[particles_List, model_String, loop_] :=
    Module[{prototype = "", function = "", i, p, f},
           For[i = 1, i <= Length[particles], i++,
               {p, f} = CreateNLoopSelfEnergy[particles[[i]], model, loop];
               prototype = prototype <> p;
               function = function <> f <> "\n";
              ];
           Return[{prototype, function}];
          ];

CreateTwoLoopSelfEnergiesSM[particles_List] :=
    CreateNLoopSelfEnergies[particles, "SM", 2];

CreateTwoLoopSelfEnergiesMSSM[particles_List] :=
    CreateNLoopSelfEnergies[particles, "MSSM", 2];

CreateTwoLoopSelfEnergiesNMSSM[particles_List] :=
    CreateNLoopSelfEnergies[particles, "NMSSM", 2];

CreateThreeLoopSelfEnergiesSplit[particles_List] :=
    CreateNLoopSelfEnergies[particles, "Split", 3];

End[];

EndPackage[];

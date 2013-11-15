
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

End[];

EndPackage[];

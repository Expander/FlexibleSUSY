
BeginPackage["SelfEnergies`", {"SARAH`", "TextFormatting`", "CConversion`", "TreeMasses`"}];

FSSelfEnergy::usage="self-energy head";
Tadpole::usage="tadpole head";

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

SetParameterReplacementRules::usage="replacement rules for model
parameters";

SetIndexReplacementRules::usage="set index replacement rules for model
parameters";

FillArrayWithOneLoopTadpoles::usage="add one-loop tadpoles to array"

Begin["Private`"];

parameterReplacementRules = {};

SetParameterReplacementRules[rules_List] :=
    parameterReplacementRules = rules;

SetParameterReplacementRules[rules_Rule] :=
    parameterReplacementRules = {rule};

indexReplacementRules = {};

SetIndexReplacementRules[rules_List] := indexReplacementRules = rules;

(* In this variable we collect the names and number of indices of the
   defined C/C++ functions to avoid double definitions *)
cppFunctionDefinitions = {};

FunctionAlreadyDefined[name_String, numberOfIndices_Integer] :=
    MemberQ[cppFunctionDefinitions, {name, numberOfIndices}];

RegisterFunction[name_String, numberOfIndices_Integer] :=
    AppendTo[cppFunctionDefinitions, {name, numberOfIndices}];

GetExpression[selfEnergy_SelfEnergies`FSSelfEnergy] :=
    selfEnergy[[2]];

GetExpression[selfEnergy_SelfEnergies`FSHeavySelfEnergy] :=
    selfEnergy[[2]];

GetExpression[tadpole_SelfEnergies`Tadpole] :=
    tadpole[[2]];

GetField[selfEnergy_SelfEnergies`FSSelfEnergy] :=
    selfEnergy[[1]];

GetField[selfEnergy_SelfEnergies`FSHeavySelfEnergy] :=
    selfEnergy[[1]];

GetField[tadpole_SelfEnergies`Tadpole] :=
    tadpole[[1]];

GetField[sym_] :=
    Module[{},
           Print["Error: GetField is not called with a SelfEnergies`FSSelfEnergy",
                 "or SelfEnergies`Tadpole head: ", sym];
           Quit[];
          ];

RemoveSMParticles[SelfEnergies`FSSelfEnergy[p_,expr__], _] :=
    SelfEnergies`FSSelfEnergy[p,expr];

RemoveSMParticles[SelfEnergies`Tadpole[p_,expr__], _] :=
    SelfEnergies`Tadpole[p,expr];

ExprContainsNonOfTheseParticles[expr_, particles_List] :=
    And @@ (FreeQ[expr,#]& /@ particles);

RemoveSMParticles[SelfEnergies`FSHeavySelfEnergy[p_,expr_], removeGoldstones_:True] :=
    Module[{strippedExpr, susyParticles, a, goldstones, g, i},
           susyParticles = TreeMasses`GetSusyParticles[];
           strippedExpr = expr /. ReplaceGhosts[];
           strippedExpr = strippedExpr //. {
               SARAH`A0[a__  /; ExprContainsNonOfTheseParticles[{a},susyParticles]] -> 0,
               SARAH`B0[a__  /; ExprContainsNonOfTheseParticles[{a},susyParticles]] -> 0,
               SARAH`B1[a__  /; ExprContainsNonOfTheseParticles[{a},susyParticles]] -> 0,
               SARAH`B00[a__ /; ExprContainsNonOfTheseParticles[{a},susyParticles]] -> 0,
               SARAH`B22[a__ /; ExprContainsNonOfTheseParticles[{a},susyParticles]] -> 0,
               SARAH`F0[a__  /; ExprContainsNonOfTheseParticles[{a},susyParticles]] -> 0,
               SARAH`G0[a__  /; ExprContainsNonOfTheseParticles[{a},susyParticles]] -> 0,
               SARAH`H0[a__  /; ExprContainsNonOfTheseParticles[{a},susyParticles]] -> 0
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
           Return[SelfEnergies`FSHeavySelfEnergy[p,strippedExpr]];
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
           (* Create Bottom self-energy with only SUSY particles in the loop *)
           heavySE = Cases[result, SelfEnergies`FSSelfEnergy[p:SARAH`BottomQuark[__][_], expr__] :>
                           SelfEnergies`FSHeavySelfEnergy[p, expr]];
           result = Join[result, RemoveSMParticles[#,False]& /@ heavySE];
           Return[result /. SARAH`Mass -> FlexibleSUSY`M];
          ];

StripIndices[expr_, indices_List] :=
    Module[{headers = {SARAH`Delta}, h, indexCombinations, removeSymbols = {}},
           indexCombinations = DeleteCases[Subsets[indices],{}];
           For[h = 1, h <= Length[headers], h++,
               AppendTo[removeSymbols, (headers[[h]] @@ #)& /@ indexCombinations];
              ];
           removeSymbols = Flatten[removeSymbols];
           removeSymbols = (Rule[#, 1])& /@ removeSymbols;
           removeSymbols = Join[removeSymbols,
                                { Rule[SARAH`Lam[__], 2],
                                  Rule[SARAH`Sig[__], 2],
                                  Rule[SARAH`fSU2[__], 1],
                                  Rule[SARAH`fSU3[__], 1] }];
           expr /. removeSymbols
          ];

GetLorentzStructure[Cp[__]] := 1;

GetLorentzStructure[Cp[__][a_]] := a;

GetParticleIndices[Cp[a__]] := DeleteDuplicates[Flatten[Cases[{a}, List[__], Infinity]]];

GetParticleIndices[Cp[a__][_]] := GetParticleIndices[Cp[a]];

CreateCouplingSymbol[coupling_] :=
    Module[{symbol, indices},
           indices = GetParticleIndices[coupling];
           symbol = ToValidCSymbol[coupling /. a_[List[__]] :> a];
           symbol[Sequence @@ indices]
          ];

(* creates a C++ function that calculates a coupling *)
CreateCouplingFunction[coupling_, expr_, strippedIndices_] :=
    Module[{symbol, prototype = "", definition = "",
            indices = {}, body = "", cFunctionName = "", i, strippedExpr,
            type, initalValue},
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
           strippedExpr = StripIndices[expr, strippedIndices];
           If[Parameters`IsRealExpression[strippedExpr],
              type = "double";  initalValue = " = 0.0";,
              type = "Complex"; initalValue = "";];
           prototype = type <> " " <> cFunctionName <> " const;\n";
           definition = type <> " CLASSNAME::" <> cFunctionName <> " const\n{\n";
           body = type <> " result" <> initalValue <> ";\n\n";
           If[FreeQ[strippedExpr,SARAH`sum] && FreeQ[strippedExpr,SARAH`ThetaStep],
              body = body <> "result = " <>
                     RValueToCFormString[Simplify[strippedExpr]] <> ";\n";
              ,
              body = body <> ExpandSums[strippedExpr, "result",
                                        type, initalValue];
             ];
           body = body <> "\nreturn result;\n";
           body = IndentText[WrapLines[body]];
           definition = definition <> body <> "}\n";
           Return[{prototype, definition}];
          ];

FindLorentzStructure[list_List, lorentz_Integer] := -I list[[lorentz, 1]];

FindLorentzStructure[list_List, lorentz_] :=
    Module[{result},
           result = Select[list, (!FreeQ[#[[2]],lorentz])&];
           If[Length[result] == 0,
              Print["Error: can't find lorentz structure ", lorentz,
                    " in list ", list];
              Return[0];
             ];
           If[Length[result] > 1,
              Print["Error: lorentz structure ", lorentz,
                    " is not unique in list ", list];
              Return[0];
             ];
           Return[-I result[[1,1]]];
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

ContainsVectorBosons[coupling_] :=
    Select[GetParticleList[coupling], IsVector] =!= {};

CreateMixingMatrixReplacementRulesFor[Cp[fields__]] :=
    CreateMixingMatrixReplacementRulesFor[{fields}];

CreateMixingMatrixReplacementRulesFor[Cp[fields__][_]] :=
    CreateMixingMatrixReplacementRulesFor[{fields}];

CreateMixingMatrixReplacementRulesFor[bar[field_]] :=
    CreateMixingMatrixReplacementRulesFor[field];

CreateMixingMatrixReplacementRulesFor[Susyno`LieGroups`conj[field_]] :=
    CreateMixingMatrixReplacementRulesFor[field];

CreateMixingMatrixReplacementRulesFor[fields_List] :=
    Flatten[CreateMixingMatrixReplacementRulesFor /@ fields];

CreateMixingMatrixReplacementRulesFor[field_[indices__] /; IsUnrotated[field]] :=
    Module[{mixingMatrix, rotatedField, rules = {}, i},
           rotatedField = ToRotatedField[field];
           mixingMatrix = Flatten[{FindMixingMatrixSymbolFor[rotatedField]}];
           (* create replacement rules for the mixing matrix *)
           For[i = 1, i <= Length[mixingMatrix], i++,
               AppendTo[rules,
                        ({mixingMatrix[[i]][#,idx_] :> SARAH`Delta[#,idx],
                          mixingMatrix[[i]][idx_,#] :> SARAH`Delta[idx,#]})& /@ indices
                       ];
              ];
           Return[DeleteDuplicates[Flatten[rules]]];
          ];

CreateMixingMatrixReplacementRulesFor[field_] := {};

(* The IsDiagonal function checks if an expression is diagonal in the
   given list of fields.  One field has be conjugated and the other
   not. *)
IsDiagonal[expr_, {Susyno`LieGroups`conj[field1_], Susyno`LieGroups`conj[field2_]}] :=
    False;

IsDiagonal[expr_, {bar[field1_], bar[field2_]}] :=
    False;

IsDiagonal[expr_, {Susyno`LieGroups`conj[field1_], bar[field2_]}] :=
    False;

IsDiagonal[expr_, {bar[field1_], Susyno`LieGroups`conj[field2_]}] :=
    False;

IsDiagonal[expr_, {field1_, Susyno`LieGroups`conj[field2_]}] :=
    IsDiagonal[expr, {Susyno`LieGroups`conj[field2], field1}];

IsDiagonal[expr_, {field1_, bar[field2_]}] :=
    IsDiagonal[expr, {bar[field2], field1}];

IsDiagonal[expr_, {Susyno`LieGroups`conj[field1_[{indices1_}]], field2_[{indices2_}]} /;
           Xor[IsUnrotated[field1], IsUnrotated[field2]]] :=
    (ToRotatedField[field1] === field2 || ToRotatedField[field2] === field1) &&
    (!FreeQ[expr, SARAH`Delta[indices1, indices2]] ||
     !FreeQ[expr, SARAH`Delta[indices2, indices1]]);

IsDiagonal[expr_, {bar[field1_[{indices1_}]], field2_[{indices2_}]} /;
           Xor[IsUnrotated[field1], IsUnrotated[field2]]] :=
    (ToRotatedField[field1] === field2 || ToRotatedField[field2] === field1) &&
    (!FreeQ[expr, SARAH`Delta[indices1, indices2]] ||
     !FreeQ[expr, SARAH`Delta[indices2, indices1]]);

IsDiagonal[expr_, fields_] := False;

(* Checks if the coupling is diagonal in a rotated and conjugate
   unrotated field. *)
HasDeltaInUnrotatedAndRotatedFields[expr_, Cp[fields__][_]] :=
    HasDeltaInUnrotatedAndRotatedFields[expr, Cp[fields]];

HasDeltaInUnrotatedAndRotatedFields[expr_, Cp[fields__]] :=
    Module[{tuples, matchingPairs = {}},
           tuples = DeleteDuplicates[Sort /@ Tuples[{fields}, 2]];
           Or[Evaluate[Sequence @@ (IsDiagonal[expr,#]& /@ tuples)]]
          ];

(* Sums over mixing matrices of the given fields to rotate the
   coupling back to gauge eigenstates. *)
RotateToGaugeEigenstate[expr_, {}] := expr;

RotateToGaugeEigenstate[expr_, {field_, rest___}] :=
    RotateToGaugeEigenstate[RotateToGaugeEigenstate[expr, field], {rest}];

RotateToGaugeEigenstate[expr_, bar[unrotatedField_[{idx_}]]] :=
    RotateToGaugeEigenstate[expr, Susyno`LieGroups`conj[unrotatedField[{idx}]]];

RotateToGaugeEigenstate[expr_, Susyno`LieGroups`conj[unrotatedField_[{idx_}]]] :=
    Module[{var, dim, mixingMatrix, field},
           field = ToRotatedField[unrotatedField];
           dim = GetDimension[field];
           mixingMatrix = FindMixingMatrixSymbolFor[field];
           If[Head[mixingMatrix] =!= List,
              var = CConversion`MakeUnique["gen"];
              Return[SARAH`sum[var, 1, dim, mixingMatrix[var,idx] (expr /. idx -> var)]];
              ,
              (* Don't know how to deal with SVD here.  So, the vertex
                 {bar[UFe],VP,Fe} will not be rotated correctly.
                 Maybe the PL and PR indicators on the original
                 coupling will help us here? *)
              Return[expr];
             ];
          ];

RotateToGaugeEigenstate[expr_, unrotatedField_[{idx_}]] :=
    Module[{var, dim, mixingMatrix, field},
           field = ToRotatedField[unrotatedField];
           dim = GetDimension[field];
           mixingMatrix = FindMixingMatrixSymbolFor[field];
           If[Head[mixingMatrix] =!= List,
              var = CConversion`MakeUnique["gen"];
              Return[SARAH`sum[var, 1, dim, Susyno`LieGroups`conj[mixingMatrix[var,idx]] (expr /. idx -> var)]];
              ,
              (* Don't know how to deal with SVD here.  So, the vertex
                 {bar[UFe],VP,Fe} will not be rotated correctly.
                 Maybe the PL and PR indicators on the original
                 coupling will help us here? *)
              Return[expr];
             ];
          ];

ReplaceMixingMatrixByIdentityIn[expr_, coupling_] :=
    Module[{unrotatedExpr, mixingMatrixReplacementRules},
           mixingMatrixReplacementRules = CreateMixingMatrixReplacementRulesFor[coupling];
           unrotatedExpr = expr /. mixingMatrixReplacementRules;
           (* It can happen that a coupling (with rotated fields) is
              already diagonal, i.e. doesn't contain mixing matrices
              (for example in VP and VG vertices).  In such cases we
              have to rotate the external fields back into a gauge
              eigenstates by hand. *)
           If[mixingMatrixReplacementRules =!= {} && unrotatedExpr === expr &&
              HasDeltaInUnrotatedAndRotatedFields[expr, coupling] &&
              ContainsVectorBosons[coupling],
              unrotatedExpr = RotateToGaugeEigenstate[expr, GetUnrotatedFields[coupling]];
             ];
           Return[unrotatedExpr];
          ];

ReplaceUnrotatedFields[Cp[p__]] :=
    Cp[Sequence @@ ToRotatedField[{p}]];

ReplaceUnrotatedFields[Cp[p__][lorentz_]] :=
    ReplaceUnrotatedFields[Cp[p]][lorentz];

FindInnerColorIndices[particles_List] :=
    Cases[Flatten[
        Select[particles,
               FreeQ[#, gO1 | gO2 | gO3 | gO4]&] //.
        {symb_[lst_][_] :> lst} //.
        {symb_[lst_] :> lst}], ct1 | ct2 | ct3 | ct4];

(* creates a C++ function that calculates a coupling *)
CreateCouplingFunctions[coupling_, sumOverInternalColors_:False] :=
    Module[{symbol, indices, cFunctionName, prototype = "", definition = "",
            vertex, i, lorentz, expr, rotatedCoupling, particles,
            allIndices = {}, neededIndices = {}, strippedIndices = {},
            innerColorIndices = {}},
           indices = GetParticleIndices[coupling];
           symbol = CreateCouplingSymbol[coupling];
           cFunctionName = ToValidCSymbolString[GetHead[symbol]];
           If[!FunctionAlreadyDefined[cFunctionName, Length[indices]],
              RegisterFunction[cFunctionName, Length[indices]];
              rotatedCoupling = ReplaceUnrotatedFields[coupling];
              vertex = SARAH`Vertex[GetParticleList[rotatedCoupling]];
              If[Head[vertex] =!= List || Length[vertex] < 2,
                 Print["Error: could not find vertex for the coupling ", rotatedCoupling];
                 Return[{0, "", ""}];
                ];
              particles = vertex[[1]];
              allIndices = DeleteDuplicates[Flatten[Cases[Cp @@ particles, List[__], Infinity]]];
              neededIndices = GetParticleIndices[coupling];
              strippedIndices = Complement[allIndices, neededIndices];
              vertex = Drop[vertex, 1]; (* drop particle list *)
              lorentz = GetLorentzStructure[coupling];
              expr = FindLorentzStructure[vertex, lorentz] /. parameterReplacementRules;
              If[sumOverInternalColors,
                 (* SARAH workaround: to resolve the color factor C
                    which appears in the self-energy of Sd and Su

                    -C sum[gI1, 1, 6,
                           A0[Mass2[Sd[{gI1}]]]
                           Cp[USd[{gO1}], conj[USd[{gO1}]], conj[Sd[{gI1}]], Sd[{gI1}]]]

                    we have to sum over the color indces of Sd[{gI1}] and Sd[{gI1}].
                    *)
                 innerColorIndices = FindInnerColorIndices[particles];
                 expr = Sum[expr,
                            {ct1, 1, If[FreeQ[innerColorIndices, ct1], 1, 3]},
                            {ct2, 1, If[FreeQ[innerColorIndices, ct2], 1, 3]},
                            {ct3, 1, If[FreeQ[innerColorIndices, ct3], 1, 3]},
                            {ct4, 1, If[FreeQ[innerColorIndices, ct4], 1, 3]}];
                ];
              (* replace mixing matrices by Delta[] for unrotated fields *)
              expr = TreeMasses`ReplaceDependencies[ReplaceMixingMatrixByIdentityIn[expr, coupling]] /.
                     Parameters`ApplyGUTNormalization[] /.
                     indexReplacementRules;
              {prototype, definition} = CreateCouplingFunction[coupling, expr, strippedIndices];
             ];
           Return[{symbol, prototype, definition}];
          ];

HasUnresolvedColorPrefactor[coupling_, expr_] :=
  !FreeQ[Coefficient[expr //. sum[__, ex_] :> ex, coupling], C];

CreateVertexExpressions[expr_] :=
    Module[{prototypes = "", functions = "", allRules = {}, allDecls = "",
            allProtos = "", i, allCouplings = {}, hasUnresolvedColorFactor,
            symbol, prototype = "", definition = "", coupling},
           allCouplings = DeleteDuplicates[Cases[expr, SARAH`Cp[__] | SARAH`Cp[__][_], Infinity]];
           (* pre-allocate list of rules *)
           allRules = Table[0, {Length[allCouplings]}];
           For[i = 1, i <= Length[allCouplings], i++,
               coupling = allCouplings[[i]];
               hasUnresolvedColorFactor = HasUnresolvedColorPrefactor[coupling, expr];
               {symbol, prototype, definition} = CreateCouplingFunctions[coupling, hasUnresolvedColorFactor];
               allRules[[i]] = Rule[coupling, symbol];
               allDecls = allDecls <> definition <> "\n";
               allProtos = allProtos <> prototype;
              ];
           Return[{allProtos, allDecls, allRules}];
          ];

CreateVertexExpressions[se_SelfEnergies`FSSelfEnergy] :=
    CreateVertexExpressions[GetExpression[se]];

CreateVertexExpressions[se_SelfEnergies`FSHeavySelfEnergy] :=
    CreateVertexExpressions[GetExpression[se]];

CreateVertexExpressions[se_SelfEnergies`Tadpole] :=
    CreateVertexExpressions[GetExpression[se]];

CreateVertexExpressions[nPointFunctions_List] :=
    Module[{k, prototypes = "", decls = "", rules,
            p, d, r},
           rules = Table[0, {Length[nPointFunctions]}];
           For[k = 1, k <= Length[nPointFunctions], k++,
               Print["   Calculating vertices appearing in Gamma[",
                     GetField[nPointFunctions[[k]]], "]"];
               {p,d,r} = CreateVertexExpressions[nPointFunctions[[k]]];
               prototypes = prototypes <> p;
               decls = decls <> d;
               rules[[k]] = r;
              ];
           Return[{prototypes, decls, Flatten[rules]}];
          ];

ReplaceGhosts[states_:SARAH`EWSB] :=
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

CreateTadpoleFunctionName[field_] :=
    "tadpole_" <> CreateFunctionNamePrefix[field];

CreateFunctionName[selfEnergy_SelfEnergies`FSSelfEnergy] :=
    CreateSelfEnergyFunctionName[GetField[selfEnergy]];

CreateFunctionName[selfEnergy_SelfEnergies`FSHeavySelfEnergy] :=
    CreateHeavySelfEnergyFunctionName[GetField[selfEnergy]];

CreateFunctionName[tadpole_SelfEnergies`Tadpole] :=
    CreateTadpoleFunctionName[GetField[tadpole]];

CreateFunctionPrototype[selfEnergy_SelfEnergies`FSSelfEnergy] :=
    CreateFunctionName[selfEnergy] <>
    "(double p " <> DeclareFieldIndices[GetField[selfEnergy]] <> ") const";

CreateFunctionPrototype[selfEnergy_SelfEnergies`FSHeavySelfEnergy] :=
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
           prototype = "Complex " <> functionName <> ";\n";
           decl = "\nComplex CLASSNAME::" <> functionName <> "\n{\n";
           body = "Complex result;\n\n" <>
                  ExpandSums[expr /. vertexRules /.
                             a_[List[i__]] :> a[i] /.
                             ReplaceGhosts[SARAH`EWSB] /.
                             parameterReplacementRules /.
                             C -> 1
                             ,"result"] <>
                  "\nreturn result * oneOver16PiSqr;";
           body = IndentText[WrapLines[body]];
           decl = decl <> body <> "\n}\n";
           Return[{prototype, decl}];
          ];

CreateNPointFunctions[nPointFunctions_List] :=
    Module[{prototypes = "", decls = "", vertexRules = {}, p, d},
           (* create coupling functions for all vertices in the list *)
           Print["Calculating vertices for all n-point functions ..."];
           {prototypes, decls, vertexRules} = CreateVertexExpressions[nPointFunctions];
           (* creating n-point functions *)
           For[k = 1, k <= Length[nPointFunctions], k++,
               Print["Writing n-point function for ", GetField[nPointFunctions[[k]]]];
               {p,d} = CreateNPointFunction[nPointFunctions[[k]], vertexRules];
               prototypes = prototypes <> p;
               decls = decls <> d;
              ];
           Return[{prototypes, decls}];
          ];

FillArrayWithOneLoopTadpoles[tadpoles_List, arrayName_String:"tadpole"] :=
    Module[{result, body, f, i = 0, d, field, functionName},
           result = "if (model->get_ewsb_loop_order() > 0) {\n";
           body = "CALC_DRBAR_PARS(model);\n";
           For[f = 1, f <= Length[tadpoles], f++,
               field = GetField[tadpoles[[f]]];
               functionName = CreateTadpoleFunctionName[field];
               For[d = 1, d <= GetDimension[field], d++; i++,
                   body = body <> arrayName <> "[" <> ToString[i] <> "] -= " <>
                          "Re(model->" <> functionName <>
                          "(" <> ToString[d] <> "));\n";
                  ];
              ];
           Return[result <> IndentText[body] <> "}\n"];
          ];

End[];

EndPackage[];

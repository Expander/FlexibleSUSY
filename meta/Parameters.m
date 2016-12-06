
BeginPackage["Parameters`", {"SARAH`", "CConversion`", "Utils`", "Phases`"}];

{ FSModelParameters, FSInputParameters, FSOutputParameters,
  FSPhysicalOutputParameters, FSPhases, FSDerivedParameters };

FindSymbolDef::usage="";

CreateSetAssignment::usage="";
CreateDisplayAssignment::usage="";
CreateParameterSARAHNames::usage="";
CreateParameterEnums::usage="";
CreateInputParameterEnum::usage="";
CreateInputParameterNames::usage="";
CreateStdVectorNamesAssignment::usage="";

CreateInputParameterArrayGetter::usage="";
CreateInputParameterArraySetter::usage="";

CreateEnumName::usage="Creates enum symbol for given parameter";
DecomposeParameter::usage="decomposes parameter into its real components";

SetParameter::usage="set model parameter";
SetSMParameter::usage="sets a SM input parameter in the QedQcd class";
SetInputParameter::usage="set input parameter to value";
AddInputParameters::usage="add an input parameter";
AddExtraParameters::usage="add an extra parameter";
SetPhases::usage="sets field phases";
GetPhases::usage="returns field phases";
SetPhase::usage="sets a phase to a value";

ApplyGUTNormalization::usage="Returns a list of repacement rules for
gauge couplings, which replace non-normalized gauge couplings
(e.g. gY) by normalized ones (e.g. g1).";

GetGUTNormalization::usage="Returns GUT normalization of the given
coupling";

CreateIndexReplacementRules::usage="";

AddRealParameter::usage="";
SetRealParameters::usage="";

SaveParameterLocally::usage="Save parameters in local variables";

GetType::usage="";
GetPhase::usage="";
HasPhase::usage="";
GetRealTypeFromDimension::usage="";
GetParameterDimensions::usage="";
GetThirdGeneration::usage="returns parameter with third generation index";

IsRealParameter::usage="";
IsComplexParameter::usage="";
IsRealExpression::usage="";
IsMatrix::usage="returns true if parameter is a matrix";
IsSymmetricMatrixParameter::usage="returns true if parameter is a matrix";
IsTensor::usage="returns true if parameter is a matrix";
IsModelParameter::usage="returns True if parameter is a model parameter";
IsInputParameter::usage="returns False if parameter is an input parameter";
IsOutputParameter::usage="returns True if parameter is a defined output parameter";
IsIndex::usage="returns true if given symbol is an index";
IsPhase::usage="returns True if given symbol is a phase";

GetIndices::usage="returns list of indices from a given parameter";

AllModelParametersAreReal::usage="returns True if all model parameters
are real, False otherwise";

SetInputParameters::usage="";
SetModelParameters::usage="";
SetOutputParameters::usage="";
SetExtraParameters::usage="";

GetInputParameters::usage="";
GetInputParametersAndTypes::usage="";
GetModelParameters::usage="";
GetOutputParameters::usage="";
GetExtraParameters::usage="";
GetExtraParametersAndTypes::usage="";
GetModelParametersWithMassDimension::usage="Returns model parameters
with given mass dimension";

GetDependenceSPhenoSymbols::usage="Returns list of symbols for which a
 DependenceSPheno rule is defined";

GetDependenceSPhenoRules::usage="Returns list of replacement rules for
 symbols for which a DependenceSPheno rule is defined";

GetOutputParameterDependencies::usage="Returns list of output
 parameters which appear in the given expression";

GetIntermediateOutputParameterDependencies::usage="Returns list of
 intermediate output parameters which appear in the given expression";

CreateLocalConstRefs::usage="creates local const references to model
parameters / input parameters.";

CreateLocalConstRefsForBetas::usage="";

CreateLocalConstRefsForInputParameters::usage="creates local const
references for all input parameters in the given expression.";

CreateLocalConstRefsForPhysicalParameters::usage="creates local const
references for all physical output parameters in the given
expression.";

FillInputParametersFromTuples::usage="";

DecreaseIndexLiterals::usage="";
IncreaseIndexLiterals::usage="";

DecreaseSumIndices::usage="";

ExpressionToString::usage="Converts an expression to a valid C++
string.";

GetEffectiveMu::usage="";
GetEffectiveMASqr::usage="";

GetParameterFromDescription::usage="Returns model parameter from a
given description string.";
GetParticleFromDescription::usage="Returns particle symbol from a
given description string.";
GetPDGCodesForParticle::usage="Returns the PDG codes for a particle."

NumberOfIndependentEntriesOfSymmetricMatrix::usage="Returns number of
independent parameters of a real symmetric nxn matrix";

ExpandExpressions::usage="";
AppendGenerationIndices::usage="";

StripIndices::usage="removes indices from model parameter";

StripIndicesRules::usage="removes given indices from a symbol";

StripSARAHIndicesRules::usage="removes SARAH-specific indices from a symbol";

FilterOutLinearDependentEqs::usage="returns linear independent equations";

FilterOutIndependentEqs::usage = "returns equations that depend on the
given list of parameters.  I.e. equations, that do not depend on the
given list of parameters are omitted from the output.";

FindAllParameters::usage = "returns list of all parameters contained
in the given expression";

FindAllParametersClassified::usage = "returns list of all parameters
 contained in the given expression, classified by their meaning.";

FindSLHABlock::usage = "returns SLHA input block name for given
 parameter";

Begin["`Private`"];

allInputParameters = {};
allExtraParameters = {};
allModelParameters = {};
allOutputParameters = {};
allPhases = {};

SetInputParameters[pars_List] := allInputParameters = DeleteDuplicates[pars];
AddInputParameters[pars_List] := allInputParameters = DeleteDuplicates[Utils`ForceJoin[allInputParameters, pars]];
SetExtraParameters[pars_List] := allExtraParameters = DeleteDuplicates[pars];
AddExtraParameters[pars_List] := allExtraParameters = DeleteDuplicates[Utils`ForceJoin[allExtraParameters, pars]];
SetModelParameters[pars_List] := allModelParameters = DeleteDuplicates[pars];
SetOutputParameters[pars_List] := allOutputParameters = DeleteDuplicates[pars];
SetPhases[phases_List]        := allPhases = DeleteDuplicates[phases];

GetInputParameters[] := First /@ allInputParameters;
GetInputParametersAndTypes[] := allInputParameters;
GetModelParameters[] := allModelParameters;
GetOutputParameters[] := allOutputParameters;
GetPhases[] := allPhases;
GetExtraParameters[] := First /@ allExtraParameters;
GetExtraParametersAndType[] := allExtraParameters;

additionalRealParameters = {};

AddRealParameter[parameter_List] :=
    additionalRealParameters = DeleteDuplicates[Join[additionalRealParameters, parameter]];

AddRealParameter[parameter_] :=
    additionalRealParameters = DeleteDuplicates[Join[additionalRealParameters, {parameter}]];

SetRealParameters[parameters_List] :=
    additionalRealParameters = parameters;

DebugPrint[msg___] :=
    If[FlexibleSUSY`FSDebugOutput,
       Print["Debug<Parameters>: ", Sequence @@ InputFormOfNonStrings /@ {msg}]];

FindSymbolDef[sym_, opt_:DependenceNum] :=
    Module[{symDef},
           symDef = Cases[SARAH`ParameterDefinitions,
                          {sym, {___, opt -> definition_, ___}} :> definition];
           If[Head[symDef] =!= List || symDef === {},
              Print["Error: Could not find definition of ",
                    sym, " in SARAH`ParameterDefinitions"];
              Return[0];
             ];
           If[Length[symDef] > 1,
              Print["Warning: ", sym, " defined multiple times"];
             ];
           symDef = symDef[[1]];
           Return[symDef];
          ];

(* Returns all parameters within an expression *)
FindAllParameters[expr_] :=
    Module[{symbols, compactExpr, allParameters, allOutPars},
           allOutPars = DeleteDuplicates[Flatten[
               Join[allOutputParameters,
                    allOutputParameters /. FlexibleSUSY`M[{a__}] :> FlexibleSUSY`M[a],
                    allOutputParameters /. FlexibleSUSY`M[{a__}] :> (FlexibleSUSY`M /@ {a})
                   ]]];
           allParameters = DeleteDuplicates[
               Join[allModelParameters, allOutPars,
                    GetInputParameters[], Phases`GetArg /@ allPhases,
                    GetDependenceSPhenoSymbols[]]];
           compactExpr = RemoveProtectedHeads[expr];
           (* find all model parameters with SARAH head *)
           symbols = DeleteDuplicates[Flatten[
               { Cases[compactExpr, SARAH`L[a_][__] /; MemberQ[allParameters,SARAH`L[a]] :> SARAH`L[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`B[a_][__] /; MemberQ[allParameters,SARAH`B[a]] :> SARAH`B[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`T[a_][__] /; MemberQ[allParameters,SARAH`T[a]] :> SARAH`T[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`Q[a_][__] /; MemberQ[allParameters,SARAH`Q[a]] :> SARAH`Q[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`L[a_]     /; MemberQ[allParameters,SARAH`L[a]] :> SARAH`L[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`B[a_]     /; MemberQ[allParameters,SARAH`B[a]] :> SARAH`B[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`T[a_]     /; MemberQ[allParameters,SARAH`T[a]] :> SARAH`T[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`Q[a_]     /; MemberQ[allParameters,SARAH`Q[a]] :> SARAH`Q[a], {0,Infinity}]
               }]];
           (* remove parameters found from compactExpr *)
           compactExpr = compactExpr /. (RuleDelayed[#, CConversion`ToValidCSymbolString[#]]& /@ symbols);
           (* find all model parameters without SARAH head in compactExpr *)
           symbols = Join[symbols,
               { Cases[compactExpr, a_Symbol /; MemberQ[allParameters,a], {0,Infinity}],
                 Cases[compactExpr, a_[__] /; MemberQ[allParameters,a] :> a, {0,Infinity}],
                 Cases[compactExpr, FlexibleSUSY`M[a_]     /; MemberQ[allOutPars,FlexibleSUSY`M[a]], {0,Infinity}],
                 Cases[compactExpr, FlexibleSUSY`M[a_[__]] /; MemberQ[allOutPars,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], {0,Infinity}]
               }];
           DeleteDuplicates[Flatten[symbols]]
          ];

IsScalar[sym_] :=
    Length[SARAH`getDimParameters[sym]] === 1 || Length[SARAH`getDimParameters[sym]] == 0;

IsMatrix[sym_[Susyno`LieGroups`i1, SARAH`i2]] :=
    IsMatrix[sym];

IsMatrix[sym_] :=
    Length[SARAH`getDimParameters[sym]] === 2;

IsSymmetricMatrixParameter[sym_[Susyno`LieGroups`i1, SARAH`i2]] :=
    IsSymmetricMatrixParameter[sym];

IsSymmetricMatrixParameter[sym_] :=
    IsMatrix[sym] && MemberQ[SARAH`ListSoftBreakingScalarMasses, sym];

IsTensor[sym_[Susyno`LieGroups`i1, SARAH`i2, SARAH`i3]] :=
    IsTensor[sym];

IsTensor[sym_] :=
    Length[SARAH`getDimParameters[sym]] > 2;

AllModelParametersAreReal[] := MemberQ[SARAH`RealParameters, All];

sarahIndices = {
    SARAH`gt1, SARAH`gt2, SARAH`gt3, SARAH`gt4,
    Susyno`LieGroups`i1 , SARAH`i2 , SARAH`i3 , SARAH`i4
};

IsIndex[i_?NumberQ] := True;
IsIndex[i_ /; MemberQ[sarahIndices,i]] := True;
IsIndex[_] := False;
IsIndex[indices_List] := And @@ (IsIndex /@ indices);
IsIndex[indices__] := IsIndex[{indices}];

GetIndices[parameter_[indices__] /; And @@ (IsIndex /@ {indices})] := {indices};
GetIndices[parameter_] := {};

IsPhase[parameter_] := MemberQ[Phases`GetArg /@ allPhases, Phases`GetArg[parameter]];

IsModelParameter[parameter_] := MemberQ[allModelParameters, parameter];
IsModelParameter[Re[parameter_]] := IsModelParameter[parameter];
IsModelParameter[Im[parameter_]] := IsModelParameter[parameter];
IsModelParameter[FlexibleSUSY`Temporary[parameter_]] := IsModelParameter[parameter];

IsModelParameter[parameter_[indices__] /; And @@ (IsIndex /@ {indices})] :=
    IsModelParameter[parameter];

IsInputParameter[parameter_] := MemberQ[GetInputParameters[], parameter];

IsOutputParameter[lst_List] := And @@ (IsOutputParameter /@ lst);
IsOutputParameter[sym_]     := MemberQ[GetOutputParameters[],sym];

IsRealParameter[Re[sym_]] := True;
IsRealParameter[Im[sym_]] := True;
IsRealParameter[FlexibleSUSY`M[_]] := True;

IsRealParameter[sym_] :=
    (IsModelParameter[sym] && AllModelParametersAreReal[]) ||
    (IsInputParameter[sym] && CConversion`IsRealType[GetType[sym]]) ||
    MemberQ[Utils`ForceJoin[SARAH`realVar, additionalRealParameters, SARAH`RealParameters], sym];

IsComplexParameter[sym_] :=
    !IsRealParameter[sym];

IsRealExpression[parameter_ /; IsModelParameter[parameter]] :=
    IsRealParameter[parameter];

IsRealExpression[expr_?NumericQ] :=
    Element[expr, Reals];

IsRealExpression[_Complex] := False;

IsRealExpression[_Real] := True;

IsRealExpression[expr_[Susyno`LieGroups`i1, SARAH`i2]] :=
    IsRealExpression[expr];

IsRealExpression[HoldPattern[SARAH`Delta[_,_]]] := True;
IsRealExpression[HoldPattern[SARAH`ThetaStep[_,_]]] := True;
IsRealExpression[Cos[_]] := True;
IsRealExpression[Sin[_]] := True;
IsRealExpression[ArcCos[_]] := True;
IsRealExpression[ArcSin[_]] := True;
IsRealExpression[ArcTan[_]] := True;
IsRealExpression[Power[a_,b_]] := IsRealExpression[a] && IsRealExpression[b];
IsRealExpression[Susyno`LieGroups`conj[expr_]] :=
    IsRealExpression[expr];
IsRealExpression[SARAH`Conj[expr_]] := IsRealExpression[expr];
IsRealExpression[Conjugate[expr_]]  := IsRealExpression[expr];
IsRealExpression[Transpose[expr_]]  := IsRealExpression[expr];
IsRealExpression[SARAH`Tp[expr_]]   := IsRealExpression[expr];
IsRealExpression[SARAH`Adj[expr_]]  := IsRealExpression[expr];
IsRealExpression[bar[expr_]]        := IsRealExpression[expr];

IsRealExpression[expr_Symbol] := IsRealParameter[expr];

IsRealExpression[Times[Conjugate[a_],a_]] := True;
IsRealExpression[Times[a_,Conjugate[a_]]] := True;
IsRealExpression[Times[Susyno`LieGroups`conj[a_],a_]] := True;
IsRealExpression[Times[a_,Susyno`LieGroups`conj[a_]]] := True;
IsRealExpression[Times[SARAH`Conj[a_],a_]]            := True;
IsRealExpression[Times[a_,SARAH`Conj[a_]]]            := True;

IsRealExpression[factors_Times] :=
    And @@ (IsRealExpression[#]& /@ (List @@ factors));

IsRealExpression[terms_Plus] :=
    And @@ (IsRealExpression[#]& /@ (List @@ terms));

IsHermitian[a_] :=
    Or[IsScalar[a] && IsRealParameter[a],
       IsSymmetricMatrixParameter[a] && IsRealParameter[a],
       MemberQ[SARAH`ListSoftBreakingScalarMasses, a]
      ];

(* helper function which calculates the adjoint of an expression *)
FSAdj[Susyno`LieGroups`conj[a_]] := SARAH`Tp[a];
FSAdj[SARAH`Tp[a_]] := a /; IsRealParameter[a];
FSAdj[SARAH`Tp[a_]] := Susyno`LieGroups`conj[a];
FSAdj[a_] := a /; IsHermitian[a];
FSAdj[a_] := SARAH`Adj[a];
FSAdj[a__] := Sequence @@ (FSAdj /@ Reverse[{a}]);

TraceEquality[SARAH`trace[a__], SARAH`trace[b__]] :=
    Or @@ (({a} === #)& /@ NestList[RotateLeft, {b}, Length[{b}] - 1]);

(* if all parameters in the trace are real, the trace is real *)
IsRealExpression[SARAH`trace[a__]] :=
    True /; And @@ (IsRealParameter /@ FindAllParameters[{a}]);

IsRealExpression[SARAH`trace[a__]] :=
    TraceEquality[SARAH`trace[a] /. SARAH`Adj -> FSAdj, SARAH`trace[FSAdj[a]]];

IsRealExpression[terms_SARAH`MatMul] :=
    And @@ (IsRealExpression[#]& /@ (List @@ terms));

IsRealExpression[sum[index_, start_, stop_, expr_]] :=
    IsRealExpression[expr];

IsRealExpression[otherwise_] := False;

HasPhase[particle_] :=
    MemberQ[#[[1]]& /@ SARAH`ParticlePhases, particle];

GetPhase[particle_ /; HasPhase[particle]] :=
    Cases[SARAH`ParticlePhases, {particle, phase_} :> phase][[1]];

GetPhase[_] := Null;

GetTypeFromDimension[sym_, {}] :=
    If[IsRealParameter[sym],
       CConversion`ScalarType[CConversion`realScalarCType],
       CConversion`ScalarType[CConversion`complexScalarCType]
      ];

GetTypeFromDimension[sym_, {0|1}] :=
    GetTypeFromDimension[sym, {}];

GetTypeFromDimension[sym_, {num_?NumberQ}] :=
    Module[{scalarType},
           scalarType = If[IsRealParameter[sym],
                           CConversion`realScalarCType,
                           CConversion`complexScalarCType
                          ];
           CConversion`VectorType[scalarType, num]
          ];

GetTypeFromDimension[sym_, {num1_?NumberQ, num2_?NumberQ}] :=
    If[IsRealParameter[sym],
       CConversion`MatrixType[CConversion`realScalarCType, num1, num2],
       CConversion`MatrixType[CConversion`complexScalarCType, num1, num2]
      ];

GetTypeFromDimension[sym_, {num1_?NumberQ, num2_?NumberQ, num3_?NumberQ}] :=
    If[IsRealParameter[sym],
       CConversion`TensorType[CConversion`realScalarCType, num1, num2, num3],
       CConversion`TensorType[CConversion`complexScalarCType, num1, num2, num3]
      ];

GetTypeFromDimension[sym_, {num1_?NumberQ, num2_?NumberQ, num3_?NumberQ, num4_?NumberQ}] :=
    If[IsRealParameter[sym],
       CConversion`TensorType[CConversion`realScalarCType, num1, num2, num3, num4],
       CConversion`TensorType[CConversion`complexScalarCType, num1, num2, num3, num4]
      ];

GetRealTypeFromDimension[{}] :=
    CConversion`ScalarType[CConversion`realScalarCType];

GetRealTypeFromDimension[{0}] :=
    GetRealTypeFromDimension[{}];

GetRealTypeFromDimension[{1}] :=
    GetRealTypeFromDimension[{}];

GetRealTypeFromDimension[{num_?NumberQ}] :=
    CConversion`VectorType[CConversion`realScalarCType, num];

GetRealTypeFromDimension[{num1_?NumberQ, num2_?NumberQ}] :=
    CConversion`MatrixType[CConversion`realScalarCType, num1, num2];

GetRealTypeFromDimension[{num1_?NumberQ, num2_?NumberQ, num3_?NumberQ}] :=
    CConversion`TensorType[CConversion`realScalarCType, num1, num2, num3];

GetRealTypeFromDimension[{num1_?NumberQ, num2_?NumberQ, num3_?NumberQ, num4_?NumberQ}] :=
    CConversion`TensorType[CConversion`realScalarCType, num1, num2, num3, num4];

GetType[FlexibleSUSY`SCALE] := GetRealTypeFromDimension[{}];

GetType[FlexibleSUSY`M[sym_]] :=
    GetRealTypeFromDimension[{SARAH`getGen[sym, FlexibleSUSY`FSEigenstates]}];

GetType[sym_?IsInputParameter] :=
    Cases[GetInputParametersAndTypes[], {sym, _, type_} :> type][[1]];

GetType[sym_] :=
    GetTypeFromDimension[sym, SARAH`getDimParameters[sym]];

GetType[sym_[indices__] /; And @@ (IsIndex /@ {indices})] :=
    GetTypeFromDimension[sym, SARAH`getDimParameters[sym]];

GetParameterDimensions[sym_] :=
    Module[{dim},
           dim = SARAH`getDimParameters[sym];
           Switch[dim,
                  {}, {1},
                  {0}, {1},
                  _, dim
                 ]
          ];

CreateIndexReplacementRule[{parameter_, CConversion`ScalarType[_]}] := {};

CreateIndexReplacementRule[{parameter_, CConversion`VectorType[_,_] | CConversion`ArrayType[_,_]}] :=
    Module[{i},
           RuleDelayed @@ Rule[parameter[i_], parameter[i-1]]
          ];

CreateIndexReplacementRule[{parameter_, CConversion`MatrixType[_,_,_]}] :=
    Module[{i,j},
           RuleDelayed @@ Rule[parameter[i_,j_], parameter[i-1,j-1]]
          ];

CreateIndexReplacementRule[{parameter_, CConversion`TensorType[_,_,_,_]}] :=
    Module[{i,j,k},
           RuleDelayed @@ Rule[parameter[i_,j_,k_], parameter[i-1,j-1,k-1]]
          ];

CreateIndexReplacementRule[{parameter_, CConversion`TensorType[_,_,_,_,_]}] :=
    Module[{i,j,k,l},
           RuleDelayed @@ Rule[parameter[i_,j_,k_,l_], parameter[i-1,j-1,k-1,l-1]]
          ];

CreateIndexReplacementRule[parameter_] :=
    Module[{i,j,k,l, dim, rule},
           dim = SARAH`getDimParameters[parameter];
           rule = {};
           Switch[Length[dim],
                  1, rule = RuleDelayed @@ Rule[parameter[i_], parameter[i-1]];,
                  2, rule = RuleDelayed @@ Rule[parameter[i_,j_], parameter[i-1,j-1]];,
                  3, rule = RuleDelayed @@ Rule[parameter[i_,j_,k_], parameter[i-1,j-1,k-1]];,
                  4, rule = RuleDelayed @@ Rule[parameter[i_,j_,k_,l_], parameter[i-1,j-1,k-1,l-1]];
                 ];
           rule
          ];

CreateIndexReplacementRules[pars_List] :=
    Flatten[CreateIndexReplacementRule /@ pars];

GetGUTNormalization[coupling_Symbol] :=
    Module[{pos, norm},
           pos = Position[SARAH`Gauge, coupling];
           If[pos =!= {},
              norm = SARAH`GUTren[pos[[1,1]]];
              If[NumericQ[norm],
                 Return[norm];
                ];
             ];
           Return[1];
          ];

ApplyGUTNormalization[] :=
    Module[{i, rules = {}, coupling},
           For[i = 1, i <= Length[SARAH`Gauge], i++,
               If[NumericQ[SARAH`GUTren[i]],
                  coupling = SARAH`Gauge[[i,4]];
                  AppendTo[rules, Rule[coupling, coupling SARAH`GUTren[i]]];
                 ];
              ];
           Return[rules];
          ];

CreateSetAssignment[name_, startIndex_, parameterType_, struct_:"pars"] :=
    Block[{},
          Print["Error: CreateSetAssignment: unknown parameter type: ", ToString[parameterType]];
          Quit[1];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`ScalarType[CConversion`realScalarCType | CConversion`integerScalarCType], struct_:"pars"] :=
    Module[{ass = ""},
           ass = name <> " = " <> struct <> "(" <> ToString[startIndex] <> ");\n";
           Return[{ass, 1}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`ScalarType[CConversion`complexScalarCType], struct_:"pars"] :=
    Module[{ass = "", type},
           type = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
           ass = name <> " = " <> type <> "(" <> struct <> "(" <> ToString[startIndex] <>
                 "), " <> struct <> "(" <> ToString[startIndex + 1] <> "));\n";
           Return[{ass, 2}];
          ];

CreateSetAssignment[name_, startIndex_, (CConversion`VectorType | CConversion`ArrayType)[CConversion`realScalarCType, rows_], struct_:"pars"] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count++,
               ass = ass <> name <> "(" <> ToString[i] <> ") = " <> struct <> "(" <>
                     ToString[startIndex + count] <> ");\n";
              ];
           If[rows != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];];
           Return[{ass, rows}];
          ];

CreateSetAssignment[name_, startIndex_, (CConversion`VectorType | CConversion`ArrayType)[CConversion`complexScalarCType, rows_], struct_:"pars"] :=
    Module[{ass = "", i, count = 0, type},
           type = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]];
           For[i = 0, i < rows, i++; count+=2,
               ass = ass <> name <> "(" <> ToString[i] <> ") = " <>
                     type <> "(" <> struct <>
                     "(" <> ToString[startIndex + count    ] <> "), " <> struct <>
                     "(" <> ToString[startIndex + count + 1] <> "));\n";
              ];
           If[2 * rows != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_], struct_:"pars"] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count++,
                   ass = ass <> name <> "(" <> ToString[i] <> "," <> ToString[j]
                         <> ") = " <> struct <> "(" <> ToString[startIndex + count] <> ");\n";
                  ];
              ];
           If[rows * cols != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`MatrixType[CConversion`complexScalarCType, rows_, cols_], struct_:"pars"] :=
    Module[{ass = "", i, j, count = 0, type},
           type = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]];
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count+=2,
                   ass = ass <> name <> "(" <> ToString[i] <> "," <> ToString[j] <>
                         ") = " <> type <> "(" <> struct <>
                         "(" <> ToString[startIndex + count    ] <> "), " <> struct <>
                         "(" <> ToString[startIndex + count + 1] <> "));\n";
                  ];
              ];
           If[2 * rows * cols != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[2 rows * cols] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`TensorType[CConversion`realScalarCType, dim1_, dim2_, dim3_]] :=
    Module[{ass = "", i, j, k, count = 0},
           For[i = 0, i < dim1, i++,
               For[j = 0, j < dim2, j++,
                   For[k = 0, k < dim3, k++; count++,
                       ass = ass <> name <> "(" <> ToString[i] <> "," <> ToString[j] <>
                             "," <> ToString[k] <>
                             ") = pars(" <> ToString[startIndex + count] <> ");\n";
                      ];
                  ];
              ];
           If[dim1 * dim2 * dim3 != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[dim1 dim2 dim3] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`TensorType[CConversion`complexScalarCType, dim1_, dim2_, dim3_]] :=
    Module[{ass = "", i, j, k, count = 0, type},
           type = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]];
           For[i = 0, i < dim1, i++,
               For[j = 0, j < dim2, j++,
                   For[k = 0, k < dim3, k++; count+=2,
                       ass = ass <> name <> "(" <> ToString[i] <> "," <> ToString[j] <> "," <>
                             ToString[k] <> ") = " <> type <> "(" <>
                             "pars(" <> ToString[startIndex + count    ] <> "), " <>
                             "pars(" <> ToString[startIndex + count + 1] <> "));\n";
                      ];
                  ];
              ];
           If[2 * dim1 * dim2 * dim3 != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[2 dim1 dim2 dim3] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

(* assigns the parameter (str) to the container element *)
DisplayAssignParSet[str_String, startIndex_, currIdx_, struct_String] :=
    struct <> "(" <> ToString[startIndex + currIdx - 1] <> ") = " <> str <> ";\n";

(* creates parameter name and calls SetterFunc[] (= DisplayAssignParSet[] by default) *)
DisplayAssignPar[(h:(Re|Im))[(par_String)[idx__]], startIndex_, currIdx_, struct_String, SetterFunc_:DisplayAssignParSet] :=
    SetterFunc[ToString[h] <> "(" <> par <> "(" <>
         Utils`StringJoinWithSeparator[{idx}, ","] <>
         "))", startIndex, currIdx, struct];

DisplayAssignPar[(par_String)[idx__], startIndex_, currIdx_, struct_String, SetterFunc_:DisplayAssignParSet] :=
    SetterFunc[par <> "(" <>
         Utils`StringJoinWithSeparator[{idx}, ","] <>
         ")", startIndex, currIdx, struct];

DisplayAssignPar[(h:(Re|Im))[par_String], startIndex_, currIdx_, struct_String, SetterFunc_:DisplayAssignParSet] :=
    SetterFunc[ToString[h] <> "(" <> par <> ")", startIndex, currIdx, struct];

DisplayAssignPar[par_, startIndex_, currIdx_, struct_String, SetterFunc_:DisplayAssignParSet] :=
    SetterFunc[CConversion`RValueToCFormString[par], startIndex, currIdx, struct];

CreateDisplayAssignment[name_, startIndex_, type_, struct_:"pars"] :=
    { StringJoin @
        MapIndexed[DisplayAssignPar[#1,startIndex,First[#2],struct]&,
                   DecomposeParameter[name, type]],
      Length[DecomposeParameter[name, type]]
    };

(* assigns the string (str) to the std container element *)
DisplayAssignNameSet[str_String, startIndex_, currIdx_, struct_String] :=
    struct <> "[" <> ToString[startIndex + currIdx - 1] <> "] = \"" <> str <> "\";\n";

(* creates parameter name and calls DisplayAssignNameSet[] *)
DisplayAssignName[par_, startIndex_, currIdx_, struct_String] :=
    DisplayAssignPar[par, startIndex, currIdx, struct, DisplayAssignNameSet];

CreateStdVectorNamesAssignment[name_, startIndex_, type_, struct_:"names"] :=
    { StringJoin @
        MapIndexed[DisplayAssignName[#1,startIndex,First[#2],struct]&,
                   DecomposeParameter[name, type]],
      Length[DecomposeParameter[name, type]]
    };

CreateParameterSARAHNameStr[par_] :=
    "\"" <> CConversion`RValueToCFormString[par] <> "\"";

CreateParameterSARAHNames[name_, type_] :=
    Utils`StringJoinWithSeparator[CreateParameterSARAHNameStr /@ DecomposeParameter[name, type], ", "];

(* decomposes a parameter into its real pieces *)
DecomposeParameter[par_, parameterType_] :=
    Block[{},
          Print["Error: DecomposeParameter: unknown parameter type: ",
                ToString[parameterType]];
          Quit[1];
          ];

DecomposeParameter[name_, CConversion`ScalarType[CConversion`realScalarCType | CConversion`integerScalarCType]] :=
    { name };

DecomposeParameter[name_, CConversion`ScalarType[CConversion`complexScalarCType]] :=
    { Re[name], Im[name] };

DecomposeParameter[name_, (CConversion`VectorType|CConversion`ArrayType)[CConversion`realScalarCType, rows_]] :=
    Array[name, rows, 0];

DecomposeParameter[name_, (CConversion`VectorType|CConversion`ArrayType)[CConversion`complexScalarCType, rows_]] :=
    Flatten[{Re[#], Im[#]}& /@ Array[name, rows, 0]];

DecomposeParameter[name_, CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_]] :=
    Flatten[Array[name, {rows, cols}, 0]];

DecomposeParameter[name_, CConversion`MatrixType[CConversion`complexScalarCType, rows_, cols_]] :=
    Flatten[{Re[#], Im[#]}& /@ Flatten[Array[name, {rows,cols}, 0]]];

DecomposeParameter[name_, CConversion`TensorType[CConversion`realScalarCType, dims__]] :=
    Flatten[Array[name, {dims}, 0]];

DecomposeParameter[name_, CConversion`TensorType[CConversion`complexScalarCType, dims__]] :=
    Flatten[{Re[#], Im[#]}& /@ Flatten[Array[name, {dims}, 0]]];

CreateEnumName[(h:(Re|Im))[par_]] :=
    CreateEnumName[h] <> CreateEnumName[par];

CreateEnumName[par_[idx__]] :=
    CreateEnumName[par] <>
    Utils`StringJoinWithSeparator[{idx}, "_", CConversion`ToValidCSymbolString];

CreateEnumName[par_] :=
    CConversion`ToValidCSymbolString[par];

CreateParameterEnums[name_, type_] :=
    Utils`StringJoinWithSeparator[CreateEnumName /@ DecomposeParameter[name, type], ", "];

CreateInputParameterEnum[inputParameters_List] :=
    Module[{result},
           result = Utils`StringJoinWithSeparator[CreateParameterEnums[#[[1]],#[[3]]]& /@ inputParameters, ", "];
           If[Length[inputParameters] > 0, result = result <> ", ";];
           "enum Input_parameters : unsigned { " <> result <> "NUMBER_OF_INPUT_PARAMETERS };\n"
          ];

CreateInputParameterNames[inputParameters_List] :=
    Module[{result},
           result = Utils`StringJoinWithSeparator[CreateParameterSARAHNames[#[[1]],#[[3]]]& /@ inputParameters, ", "];
           "const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS> input_parameter_names = {" <>
           result <> "};\n"
          ];

SetInputParameter[parameter_, value_, wrapper_String, castToType_:None] :=
    Module[{parameterStr, valueStr},
           If[IsInputParameter[parameter],
              parameterStr = CConversion`ToValidCSymbolString[parameter];
              valueStr = CConversion`RValueToCFormString[value];
              If[wrapper == "",
                 parameterStr <> " = " <> CConversion`CastTo[valueStr,castToType] <> ";\n",
                 wrapper <> "(" <> parameterStr <> ") = " <> CConversion`CastTo[valueStr,castToType] <> ";\n"
                ]
              ,
              ""
             ]
          ];

SetPhase[phase_, value_, class_String] :=
    Module[{phaseStr, valueStr},
           If[IsPhase[phase],
              phaseStr = Phases`CreatePhaseName[phase];
              valueStr = CConversion`RValueToCFormString[value];
              class <> "->set_" <> phaseStr <> "(" <> valueStr <> ");\n",
              ""
             ]
          ];

ConcatIndices[indices_List] :=
    Utils`StringJoinWithSeparator[ToString /@ indices,","];

CreateIndices[parameter_[indices__] /; And @@ (IsIndex /@ {indices})] :=
    "(" <> ConcatIndices[{indices}] <> ")";

CreateIndices[parameter_] := "";

AppendIfNotEmpty[str_String, sym_String] := If[str == "", "", str <> sym];

SetParameter[Re[parameter_], value_String, class_String, castToType_:None] :=
    Module[{parameterStr, newValue},
           parameterStr = CConversion`RValueToCFormString[parameter];
           newValue = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]] <>
               "(" <> value <> ",Im(MODELPARAMETER(" <> parameterStr <> ")))";
           SetParameter[parameter, newValue, class, castToType]
          ];

SetParameter[Im[parameter_], value_String, class_String, castToType_:None] :=
    Module[{parameterStr, newValue},
           parameterStr = CConversion`RValueToCFormString[parameter];
           newValue = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]] <>
               "(Re(MODELPARAMETER(" <> parameterStr <> "))," <> value <> ")";
           SetParameter[parameter, newValue, class, castToType]
          ];

SetParameter[parameter_, value_String, class_String, castToType_:None] :=
    Module[{parameterStr, targetType = castToType},
           If[IsModelParameter[parameter],
              parameterStr = CConversion`ToValidCSymbolString[StripIndices[parameter]];
              (* if the parameter indices, we need to cast to the element type *)
              If[GetIndices[parameter] =!= {} && targetType =!= None,
                 targetType = CConversion`GetScalarElementType[targetType];
                ];
              class <> "set_" <> parameterStr <> "(" <>
              AppendIfNotEmpty[ConcatIndices[GetIndices[parameter]],","] <>
              CConversion`CastTo[value,targetType] <> ");\n"
              ,
              ""
             ]
          ];

SetParameter[parameter_, value_, class_String] :=
    SetParameter[parameter, CConversion`RValueToCFormString[value], class, GetType[parameter]]

SetParameter[Re[parameter_], value_, castToType_:CConversion`ScalarType[CConversion`realScalarCType]] :=
    CConversion`ToValidCSymbolString[StripIndices[parameter]] <> ".real(" <>
    CConversion`CastTo[CConversion`RValueToCFormString[value] <> CreateIndices[parameter],
                       castToType] <> ");\n";

SetParameter[Im[parameter_], value_, castToType_:CConversion`ScalarType[CConversion`realScalarCType]] :=
    CConversion`ToValidCSymbolString[StripIndices[parameter]] <> ".imag(" <>
    CConversion`CastTo[CConversion`RValueToCFormString[value] <> CreateIndices[parameter],
                       castToType] <> ");\n";

SetParameter[parameter_, value_, castToType_:None] :=
    CConversion`ToValidCSymbolString[StripIndices[parameter]] <> CreateIndices[parameter] <> " = " <>
    CConversion`CastTo[CConversion`RValueToCFormString[value],castToType] <> ";\n";

SaveParameterLocally[parameters_List] :=
    StringJoin[SaveParameterLocally /@ parameters];

SaveParameterLocally[parameter_] :=
    Module[{ parStr, parStrSym },
           parStr = CConversion`RValueToCFormString[parameter];
           parStrSym = CConversion`ToValidCSymbolString[parameter];
           "const auto save_" <> parStrSym <> "_raii = make_raii_save(" <> parStr <> ");\n"
          ];

RemoveProtectedHeads[expr_] :=
    expr /. { FlexibleSUSY`LowEnergyConstant[__] -> FlexibleSUSY`LowEnergyConstant[],
              FlexibleSUSY`Pole[__]  -> FlexibleSUSY`Pole[],
              FlexibleSUSY`BETA[__] -> FlexibleSUSY`BETA[]};

DefineLocalConstCopy[parameter_, macro_String, prefix_String:""] :=
    "const auto " <> prefix <> ToValidCSymbolString[parameter] <> " = " <>
    macro <> "(" <> ToValidCSymbolString[parameter] <> ");\n";

PrivateCallLoopMassFunction[FlexibleSUSY`M[particle_Symbol]] :=
    "calculate_" <> ToValidCSymbolString[FlexibleSUSY`M[particle]] <> "_pole();\n";

CalculateLocalPoleMasses[parameter_] :=
    "MODEL->" <> PrivateCallLoopMassFunction[parameter];

FindAllParametersClassified[expr_] :=
    Module[{symbols = DeleteDuplicates[Flatten[FindAllParameters[expr]]],
            inputPars, modelPars, outputPars,
            poleMasses, phases, depNum, allOutPars},
           allOutPars = DeleteDuplicates[Flatten[
               Join[allOutputParameters,
                    allOutputParameters /. FlexibleSUSY`M[{a__}] :> FlexibleSUSY`M[a],
                    allOutputParameters /. FlexibleSUSY`M[{a__}] :> (FlexibleSUSY`M /@ {a})
                   ]]];
           poleMasses = {
               Cases[expr, FlexibleSUSY`Pole[FlexibleSUSY`M[a_]]     /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], {0,Infinity}],
               Cases[expr, FlexibleSUSY`Pole[FlexibleSUSY`M[a_[__]]] /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], {0,Infinity}]
                        };
           poleMasses   = DeleteDuplicates[Flatten[poleMasses]];
           inputPars    = DeleteDuplicates[Select[symbols, (MemberQ[GetInputParameters[],#])&]];
           modelPars    = DeleteDuplicates[Select[symbols, (MemberQ[allModelParameters,#])&]];
           outputPars   = DeleteDuplicates[Select[symbols, (MemberQ[allOutPars,#])&]];
           phases       = DeleteDuplicates[Select[symbols, (MemberQ[Phases`GetArg /@ allPhases,#])&]];
           depNum       = DeleteDuplicates[Select[symbols, (MemberQ[GetDependenceSPhenoSymbols[],#])&]];
           {
               FSModelParameters -> modelPars,
               FSInputParameters -> inputPars,
               FSOutputParameters -> outputPars,
               FSPhysicalOutputParameters -> poleMasses,
               FSPhases -> phases,
               FSDerivedParameters -> depNum
           }
          ];

CreateLocalConstRefs[expr_] :=
    Module[{result = "", pars, inputSymbols, modelPars, outputPars,
            poleMasses, phases, depNum},
           pars = FindAllParametersClassified[expr];
           inputSymbols = FSInputParameters /. pars;
           modelPars    = FSModelParameters /. pars;
           outputPars   = FSOutputParameters /. pars;
           phases       = FSPhases /. pars;
           depNum       = FSDerivedParameters /. pars;
           poleMasses   = FSPhysicalOutputParameters /. pars;
           (result = result <> DefineLocalConstCopy[#,"INPUTPARAMETER"])& /@ inputSymbols;
           (result = result <> DefineLocalConstCopy[#,"MODELPARAMETER"])& /@ modelPars;
           (result = result <> DefineLocalConstCopy[#,"MODELPARAMETER"])& /@ outputPars;
           (result = result <> DefineLocalConstCopy[#,"PHASE"         ])& /@ phases;
           (result = result <> DefineLocalConstCopy[#,"DERIVEDPARAMETER"])& /@ depNum;
           (result = result <> CalculateLocalPoleMasses[#])& /@ poleMasses;
           Return[result];
          ];

CreateLocalConstRefsForPhysicalParameters[expr_] :=
    Module[{result = "", pars, outputPars},
           pars = FindAllParametersClassified[expr];
           outputPars = FSOutputParameters /. pars;
           (result = result <> DefineLocalConstCopy[#,"PHYSICAL"])& /@ outputPars;
           Return[result];
          ];

DefineLocalConstCopyForBeta[{par_, -1}] :=
    DefineLocalConstCopy[par, "BETAPARAMETER", "beta_"];

DefineLocalConstCopyForBeta[{par_, loops_}] :=
    Module[{lstr = ToString[loops], parStr = ToValidCSymbolString[par]},
           "const auto BETA1(" <> lstr <> "," <> parStr <>
           ") = BETAPARAMETER1(" <> lstr <> "," <> parStr <> ");\n"
          ];

CreateLocalConstRefsForBetas[expr_] :=
    Module[{result = "", exprFull, pars},
           (* add index to beta function for current loop level *)
           exprFull = expr /. FlexibleSUSY`BETA[p_] :> FlexibleSUSY`BETA[-1,p];
           pars = DeleteDuplicates @ \
                  Cases[exprFull, FlexibleSUSY`BETA[l_,p_] | FlexibleSUSY`BETA[l_,p_][___] :> {p,l}, {0, Infinity}];
           StringJoin[DefineLocalConstCopyForBeta /@ pars]
          ];

CreateLocalConstRefsForInputParameters[expr_, head_String:"INPUT"] :=
    Module[{result = "", pars, inputPars},
           pars = FindAllParametersClassified[expr];
           inputPars = FSInputParameters /. pars;
           (result = result <> DefineLocalConstCopy[#, head])& /@ inputPars;
           Return[result];
          ];

SetInputParameterTo[Re[par_], value_String] :=
    CConversion`ToValidCSymbolString[par] <> ".real(" <> value <> ");";

SetInputParameterTo[Im[par_], value_String] :=
    CConversion`ToValidCSymbolString[par] <> ".imag(" <> value <> ");";

SetInputParameterTo[par_, value_String] :=
    CConversion`ToValidCSymbolString[par] <> " = " <> value <> ";";

CreateCaseFromTuple[{key_?NumberQ, parameter_}] :=
    "case " <> ToString[key] <> ": input." <>
    SetInputParameterTo[parameter, "value"] <> " break;\n";

CreateCaseFromTuple[expr_] :=
    Block[{},
          Print["Error: not a valid {key,parameter} tuple: ", expr];
          ""
         ];

FillInputParametersFromTuples[minpar_List, blockName_String] :=
    Module[{result = ""},
           (result = result <> CreateCaseFromTuple[#])& /@ minpar;
           result = "switch (key) {\n" <> result <>
                    "default: WARNING(\"Unrecognized entry in block " <>
                    blockName <> ": \" << key); break;\n}\n";
           Return[result];
          ];

IncreaseIndex[ind_Integer, num_Integer] := ind + num;
IncreaseIndex[ind_, _]     := ind;
IncreaseIndices[a_[{ind__}], num_Integer] := a[IncreaseIndex[#,num]& /@ {ind}];
IncreaseIndices[a_[ind__], num_Integer] := a[Sequence @@ (IncreaseIndex[#,num]& /@ {ind})];
IncreaseIndices[a_, _]     := a;
IncreaseIndices[SARAH`Delta[a_, b_], num_Integer] :=
    CConversion`FSKroneckerDelta[IncreaseIndex[a,num], IncreaseIndex[b,num]];

IncreaseIndexLiterals[expr_] :=
    IncreaseIndexLiterals[expr, 1];

IncreaseIndexLiterals[expr_, num_Integer] :=
    IncreaseIndexLiterals[expr, num, Join[GetInputParameters[], allModelParameters,
                                          allOutputParameters]];

IncreaseIndexLiterals[expr_, num_Integer, heads_List] :=
    Module[{indexedSymbols, rules, decrExpr, allHeads},
           allHeads = Join[heads /. FlexibleSUSY`M -> Identity, {SARAH`Delta, SARAH`ThetaStep}];
           indexedSymbols = Cases[{expr}, s_[__] /; MemberQ[allHeads, s], Infinity];
           rules = Rule[#, IncreaseIndices[#,num]] & /@ indexedSymbols;
           decrExpr = expr /. rules;
           Return[decrExpr]
          ];

DecreaseIndexLiterals[expr_] :=
    IncreaseIndexLiterals[expr, -1];

DecreaseIndexLiterals[expr_, heads_List] :=
    IncreaseIndexLiterals[expr, -1, heads];

DecreaseSumIndices[expr_] :=
    expr //. SARAH`sum[idx_, start_, stop_, exp_] :> FlexibleSUSY`SUM[idx, start - 1, stop - 1, exp];

ReplaceThetaStep[expr_] :=
    Expand[expr] //. {
        Times[a__, ThetaStep[i1_, i2_], b__] :> FlexibleSUSY`IF[i1 < i2, a b, 0],
        Times[a__, ThetaStep[i1_, i2_]] :> FlexibleSUSY`IF[i1 < i2, a, 0],
        Times[ThetaStep[i1_, i2_], a__] :> FlexibleSUSY`IF[i1 < i2, a, 0]
    };

(* Converts an expression to a valid C++ string. *)
ExpressionToString[expr_] :=
    CConversion`RValueToCFormString[
        Simplify[Parameters`DecreaseIndexLiterals[Parameters`DecreaseSumIndices[ReplaceThetaStep[expr]]]]];

ExpressionToString[expr_, heads_] :=
    CConversion`RValueToCFormString[
        Simplify[Parameters`DecreaseIndexLiterals[Parameters`DecreaseSumIndices[ReplaceThetaStep[expr]], heads]]];

GetEffectiveMu[] :=
    Module[{},
           If[!ValueQ[FlexibleSUSY`EffectiveMu],
              Print["Error: effective Mu parameter not defined!"];
              Print["   Please set EffectiveMu to the expression of the"];
              Print["   effective Mu parameter."];
              Quit[1];
             ];
           FlexibleSUSY`EffectiveMu
          ];

GetEffectiveMASqr[] :=
    Module[{},
           If[!ValueQ[FlexibleSUSY`EffectiveMASqr],
              Print["Error: effective CP-odd Higgs mass not defined!"];
              Print["   Please set EffectiveMASqr to the expression of the"];
              Print["   effective CP-odd Higgs mass."];
              Quit[1];
             ];
           FlexibleSUSY`EffectiveMASqr
          ];

GetParameterFromDescription[description_String] :=
    Module[{parameter},
           parameter =Cases[SARAH`ParameterDefinitions,
                            {parameter_,
                             {___, SARAH`Description -> description, ___}} :>
                            parameter];
           If[Length[parameter] == 0,
              Print["Error: Parameter with description \"", description,
                    "\" not found."];
              Return[Null];
             ];
           If[Length[parameter] > 1,
              Print["Warning: Parameter with description \"", description,
                    "\" not unique."];
             ];
           parameter[[1]]
          ];

GetParticleFromDescription[description_String, eigenstates_:FlexibleSUSY`FSEigenstates] :=
    Module[{particle},
           particle =Cases[SARAH`ParticleDefinitions[eigenstates],
                            {particle_,
                             {___, SARAH`Description -> description, ___}} :>
                            particle];
           If[Length[particle] == 0,
              DebugPrint["Note: Particle with description \"", description,
                         "\" not found."];
              Return[Null];
             ];
           If[Length[particle] > 1,
              Print["Warning: Particle with description \"", description,
                    "\" not unique."];
             ];
           particle[[1]]
          ];

GetParticleFromDescription[multipletName_String, splitNames_List] :=
    Module[{result},
           result = GetParticleFromDescription[multipletName];
           If[result =!= Null, Return[{result}]];
           DeleteCases[GetParticleFromDescription /@ splitNames, Null]
          ];

GetPDGCodesForParticle[particle_] :=
    Module[{pdgList},
            pdgList = SARAH`getPDGList[particle];
            If[pdgList === None,
               pdgList = {};
              ];
           pdgList
          ];

NumberOfIndependentEntriesOfSymmetricMatrix[n_] := (n^2 + n) / 2;

AppendGenerationIndices[expr_List] :=
    AppendGenerationIndices /@ expr;

AppendGenerationIndices[expr_Symbol] :=
    Switch[SARAH`getDimParameters[expr],
           {}                          , expr,
           {0}                         , expr,
           {1}                         , expr,
           {idx_}                      , expr[SARAH`gt1],
           {idx1_, idx2_}              , expr[SARAH`gt1, SARAH`gt2],
           {idx1_, idx2_, idx3_}       , expr[SARAH`gt1, SARAH`gt2, SARAH`gt3],
           {idx1_, idx2_, idx3_, idx4_}, expr[SARAH`gt1, SARAH`gt2, SARAH`gt3, SARAH`gt4]
          ];

AppendGenerationIndices[expr_] := expr;

(*
 * Expands a list of expressions of the form
 *
 *   { 1 + A[SARAH`gt1] }
 *
 * to
 *
 *   { 1 + A[1], 1 + A[2], 1 + A[3] }
 *
 * where the indices SARAH`gt1 and SARAH`gt2 are assumed to run from 1
 * to their maximum value.
 *)
ExpandExpressions[eqs_List] :=
    Join[Flatten[ExpandExpressions /@ eqs]];

GetIdxRange[{idx_,par_}] :=
    Module[{dim = GetParameterDimensions[par]},
           Switch[dim,
                  {}  , {idx,1,1},
                  {0} , {idx,1,1},
                  {n_?NumberQ}, {idx,1,dim[[1]]},
                  {__}, {idx,1,1}
                 ]
          ];

ExpandExpressions[eq_] :=
    Module[{par, indexSymbols, indexRanges, indices = {SARAH`gt1, SARAH`gt2}},
           indexSymbols = DeleteDuplicates[
               Join @@ (Cases[eq, par_[idx_ /; !FreeQ[idx,#]] /;
                                  MemberQ[Join[GetModelParameters[], GetOutputParameters[]], par] :> {#,par},
                              {0,Infinity}]& /@ indices)];
           indexRanges  = DeleteDuplicates[GetIdxRange /@ indexSymbols];
           If[indexRanges === {},
              eq,
              Table[eq, Evaluate[Sequence @@ indexRanges]]
             ]
          ];

(* removes indices from model Parameter, taking SARAH's L, B, T, Q
   into account.  *)
StripIndices[par_[idx___] /; And @@ (NumberQ /@ {idx})] := par;

StripIndices[par_[idx___] /; MemberQ[Join[allModelParameters,allOutputParameters],par]] := par;

StripIndices[par_] := par;

StripIndicesRules[indices_List, numberOfIndices_] :=
    Flatten[{ RuleDelayed[a_[Sequence @@ #], a]& /@ Permutations[indices, {numberOfIndices}] }];

StripSARAHIndicesRules[numberOfIndices_] :=
    StripIndicesRules[sarahIndices, numberOfIndices];

ExtractParametersFromSARAHBetaLists[beta_List] :=
    StripIndices[#[[1]]]& /@ beta;

GetModelParametersWithMassDimension[dim_?IntegerQ] :=
    Module[{dimPars},
           Switch[dim,
                  0, dimPars = Join[SARAH`BetaGauge, SARAH`BetaLijkl, SARAH`BetaYijk, SARAH`BetaQijkl];,
                  1, dimPars = Join[SARAH`BetaMuij, SARAH`BetaTijk, SARAH`BetaMi, SARAH`BetaDGi, SARAH`BetaVEV];,
                  2, dimPars = Join[SARAH`BetaLi, SARAH`BetaBij, SARAH`Betam2ij];,
                  3, dimPars = Join[SARAH`BetaLSi];,
                  _,
                  Print["Error: GetModelParametersWithMassDimension: ", dim,
                        " is not a valid dimension"];
                  Return[{}];
                 ];
           ExtractParametersFromSARAHBetaLists[dimPars]
          ];

AreLinearDependent[{eq1_, eq2_}, parameters_List] :=
    Module[{frac = Simplify[eq1/eq2 /. FlexibleSUSY`tadpole[_] -> 0],
            pars},
           (* ignore parameter heads Re[], Im[], Abs[], Phase[] *)
           pars = parameters /. { Re[p_] :> p, Im[p_] :> p,
                                  Abs[p_] :> p, FlexibleSUSY`Phase[p_] :> p };
           And @@ (FreeQ[frac,#]& /@ pars)
          ];

FilterOutLinearDependentEqs[{}, _List] := {};

FilterOutLinearDependentEqs[{eq_}, _List] := {eq};

FilterOutLinearDependentEqs[{eq_, rest__}, parameters_List] :=
    If[Or @@ (AreLinearDependent[#,parameters]& /@ ({eq,#}& /@ {rest})),
       (* leave out eq and check rest *)
       FilterOutLinearDependentEqs[{rest}, parameters],
       (* keep eq and check rest *)
       {eq, Sequence @@ FilterOutLinearDependentEqs[{rest}, parameters]}
      ];

FilterOutIndependentEqs[eqs_List, pars_List] :=
    DeleteDuplicates @ Flatten @ Join[FilterOutIndependentEqs[eqs,#]& /@ pars];

FilterOutIndependentEqs[eqs_List, p_] :=
    Select[eqs, (!FreeQ[#,p])&];

GetThirdGeneration[par_] :=
    Which[IsScalar[par], par,
          IsMatrix[par], par[2,2],
          True, Print["Warning: GetThirdGeneration[",par,"]: unknown type"]; par
         ];

GetSARAHParameters[] :=
    (#[[1]])& /@ SARAH`SARAHparameters;

GetAllDependenceSPhenoSymbols[] :=
    DeleteDuplicates @ Flatten @
    Cases[SARAH`ParameterDefinitions,
          {parameter_, {___, SARAH`DependenceSPheno -> value:Except[None], ___}} :> parameter];

GetAllDependenceSPhenoRules[] :=
    Cases[SARAH`ParameterDefinitions,
          {parameter_, {___, SARAH`DependenceSPheno -> value:Except[None], ___}} :> RuleDelayed[parameter, value]];

GetDependenceSPhenoSymbols[] :=
    Module[{sarahPars = GetSARAHParameters[]},
           Select[GetAllDependenceSPhenoSymbols[], MemberQ[sarahPars,#]&]
          ];

GetDependenceSPhenoRules[] :=
    Module[{sarahPars = GetSARAHParameters[]},
           Select[GetAllDependenceSPhenoRules[], MemberQ[sarahPars,#[[1]]]&]
          ];

GetAllOutputParameterDependencies[expr_] :=
    Complement[Select[Join[GetSARAHParameters[],
                           GetDependenceSPhenoSymbols[]],
                      (!FreeQ[expr,#])&],
               GetModelParameters[]];

GetAllOutputParameterDependenciesReplaced[expr_] :=
    DeleteCases[GetAllOutputParameterDependencies[expr /. GetDependenceSPhenoRules[]], _?NumericQ];

GetOutputParameterDependencies[expr_] :=
    Select[GetOutputParameters[],
           (!FreeQ[GetAllOutputParameterDependenciesReplaced[expr],#])&];

GetExponent[a_^b_] := -I b;
GetExponent[a_]    := a;

GetIntermediateOutputParameterDependencies[expr_] :=
    Complement[
        GetAllOutputParameterDependenciesReplaced[expr],
        Join[GetOutputParameters[], GetInputParameters[], GetExponent /@ GetPhases[]]
    ];

CreateInputParameterArrayGetter[{}] :=
    "return Eigen::ArrayXd();\n";

CreateInputParameterArrayGetter[inputParameters_List] :=
    Module[{get = "", paramCount = 0, name = "", par,
            type, i, assignment = "", nAssignments = 0},
           For[i = 1, i <= Length[inputParameters], i++,
               par  = inputParameters[[i,1]];
               type = inputParameters[[i,2]];
               name = CConversion`ToValidCSymbolString[par];
               {assignment, nAssignments} = Parameters`CreateDisplayAssignment[name, paramCount, type];
               get = get <> assignment;
               paramCount += nAssignments;
              ];
           get = "Eigen::ArrayXd pars(" <> ToString[paramCount] <> ");\n\n" <>
                 get <> "\n" <>
                 "return pars;";
           Return[get];
          ];

CreateInputParameterArraySetter[inputParameters_List] :=
    Module[{set = "", paramCount = 0, name = "", par,
            type, i, assignment = "", nAssignments = 0},
           For[i = 1, i <= Length[inputParameters], i++,
               par  = inputParameters[[i,1]];
               type = inputParameters[[i,2]];
               name = CConversion`ToValidCSymbolString[par];
               {assignment, nAssignments} = Parameters`CreateSetAssignment[name, paramCount, type];
               set = set <> assignment;
               paramCount += nAssignments;
              ];
           Return[set];
          ];

FindSLHABlock[blockList_List, par_] :=
    Module[{foundBlocks},
           foundBlocks = Cases[blockList, {par, block_, ___} :> block];
           If[foundBlocks === {},
              Print["Error: FindSLHABlock: no input block defined for parameter ", par];
              Quit[1];
             ];
           foundBlocks[[1]]
          ];

SetSMParameter[FlexibleSUSY`AlphaEMInvInput    , value_String, struct_String] := struct <> ".setAlphaEmInput(1./(" <> value <> "))";
SetSMParameter[FlexibleSUSY`GFermiInput        , value_String, struct_String] := struct <> ".setFermiConstant(" <> value <> ")";
SetSMParameter[FlexibleSUSY`AlphaSInput        , value_String, struct_String] := struct <> ".setAlphaSInput(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MZPoleInput        , value_String, struct_String] := struct <> ".setPoleMZ(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MBottomMbottomInput, value_String, struct_String] := struct <> ".setMbMb(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MTopPoleInput      , value_String, struct_String] := struct <> ".setPoleMt(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MTauPoleInput      , value_String, struct_String] := struct <> ".setPoleMtau(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MNeutrino3PoleInput, value_String, struct_String] := struct <> ".setNeutrinoPoleMass(3, " <> value <> ")";
SetSMParameter[FlexibleSUSY`MWPoleInput        , value_String, struct_String] := struct <> ".setPoleMW(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MElectronPoleInput , value_String, struct_String] := struct <> ".setPoleMel(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MNeutrino1PoleInput, value_String, struct_String] := struct <> ".setNeutrinoPoleMass(1, " <> value <> ")";
SetSMParameter[FlexibleSUSY`MMuonPoleInput     , value_String, struct_String] := struct <> ".setPoleMmuon(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MNeutrino2PoleInput, value_String, struct_String] := struct <> ".setNeutrinoPoleMass(2, " <> value <> ")";
SetSMParameter[FlexibleSUSY`MDown2GeVInput     , value_String, struct_String] := struct <> ".setMass(softsusy::mDown, " <> value <> ")";
SetSMParameter[FlexibleSUSY`MUp2GeVInput       , value_String, struct_String] := struct <> ".setMass(softsusy::mUp, " <> value <> ")";
SetSMParameter[FlexibleSUSY`MStrange2GeVInput  , value_String, struct_String] := struct <> ".setMass(softsusy::mStrange, " <> value <> ")";
SetSMParameter[FlexibleSUSY`MCharmMCharm       , value_String, struct_String] := struct <> ".setMass(softsusy::mCharm, " <> value <> ")";

End[];

EndPackage[];

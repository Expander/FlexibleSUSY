
BeginPackage["Parameters`", {"SARAH`", "CConversion`"}];

FindSymbolDef::usage="";

CreateSetAssignment::usage="";
CreateDisplayAssignment::usage="";
CreateParameterNamesStr::usage="";
CreateParameterEnums::usage="";

SetParameter::usage="set model parameter";

ApplyGUTNormalization::usage="Returns a list of repacement rules for
gauge couplings, which replace non-normalized gauge couplings
(e.g. gY) by normalized ones (e.g. g1).";

GetGUTNormalization::usage="Returns GUT normalization of the given
coupling";

CreateIndexReplacementRules::usage="";

AddRealParameter::usage="";

SaveParameterLocally::usage="Save parameters in local variables";
RestoreParameter::usage="Restore parameters from local variables";

GetType::usage="";
GetPhase::usage="";
HasPhase::usage="";

IsRealParameter::usage="";
IsComplexParameter::usage="";
IsRealExpression::usage="";
IsMatrix::usage="returns true if parameter is a matrix";
IsSymmetricMatrixParameter::usage="returns true if parameter is a matrix";

SetInputParameters::usage="";
SetModelParameters::usage="";
SetOutputParameters::usage="";

GetInputParameters::usage="";
GetModelParameters::usage="";
GetOutputParameters::usage="";

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

DecreaseSumIdices::usage="";

Begin["`Private`"];

allInputParameters = {};
allModelParameters = {};
allOutputParameters = {};

SetInputParameters[pars_List] := allInputParameters = pars;
SetModelParameters[pars_List] := allModelParameters = pars;
SetOutputParameters[pars_List] := allOutputParameters = pars;

GetInputParameters[] := allInputParameters;
GetModelParameters[] := allModelParameters;
GetOutputParameters[] := allOutputParameters;

additionalRealParameters = {};

AddRealParameter[parameter_List] :=
    additionalRealParameters = DeleteDuplicates[Join[additionalRealParameters, parameter]];

AddRealParameter[parameter_] :=
    additionalRealParameters = DeleteDuplicates[Join[additionalRealParameters, {parameter}]];

FindSymbolDef[sym_] :=
    Module[{symDef},
           symDef = Cases[SARAH`ParameterDefinitions,
                          {sym, {___, DependenceNum -> definition_, ___}} :> definition];
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
    Module[{symbols, compactExpr},
           compactExpr = RemoveProtectedHeads[expr];
           symbols = { Cases[compactExpr, _Symbol, Infinity],
                       Cases[compactExpr, a_[__] /; MemberQ[allModelParameters,a] :> a, Infinity],
                       Cases[compactExpr, a_[__] /; MemberQ[allOutputParameters,a] :> a, Infinity],
                       Cases[compactExpr, FlexibleSUSY`M[a_]     /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]], Infinity],
                       Cases[compactExpr, FlexibleSUSY`M[a_[__]] /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], Infinity]
                     };
           Return[DeleteDuplicates[Flatten[symbols]]];
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

IsRealParameter[sym_] :=
    MemberQ[Join[SARAH`realVar, additionalRealParameters], sym];

IsComplexParameter[sym_] :=
    !IsRealParameter[sym];

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
    If[True || IsRealParameter[sym],
       CConversion`ScalarType[CConversion`realScalarCType],
       CConversion`ScalarType[CConversion`complexScalarCType]
      ];

GetTypeFromDimension[sym_, {1}] :=
    GetTypeFromDimension[sym, {}];

GetTypeFromDimension[sym_, {num_?NumberQ}] :=
    Module[{scalarType},
           scalarType = If[True || IsRealParameter[sym],
                           CConversion`realScalarCType,
                           CConversion`complexScalarCType
                          ];
           CConversion`VectorType[scalarType, num]
          ];

GetTypeFromDimension[sym_, {num1_?NumberQ, num2_?NumberQ}] :=
    If[True || IsRealParameter[sym],
       CConversion`MatrixType[CConversion`realScalarCType, num1, num2],
       CConversion`MatrixType[CConversion`complexScalarCType, num1, num2]
      ];

GetType[FlexibleSUSY`M[sym_]] :=
    GetTypeFromDimension[sym, {SARAH`getGen[sym, FlexibleSUSY`FSEigenstates]}];

GetType[sym_] :=
    GetTypeFromDimension[sym, SARAH`getDimParameters[sym]];

CreateIndexReplacementRules[pars_List] :=
    Module[{indexReplacementRules, p, i,j,k,l, dim, rule, parameter},
           indexReplacementRules = {};
           For[p = 1, p <= Length[pars], p++,
               parameter = pars[[p]];
               dim = SARAH`getDimParameters[parameter];
               rule = {};
               Switch[Length[dim],
                      1, rule = RuleDelayed @@ Rule[parameter[i_], parameter[i-1]];,
                      2, rule = RuleDelayed @@ Rule[parameter[i_,j_], parameter[i-1,j-1]];,
                      3, rule = RuleDelayed @@ Rule[parameter[i_,j_,k_], parameter[i-1,j-1,k-1]];,
                      4, rule = RuleDelayed @@ Rule[parameter[i_,j_,k_,l_], parameter[i-1,j-1,k-1,l-1]];
                     ];
               AppendTo[indexReplacementRules, rule];
              ];
           Return[Flatten[indexReplacementRules]]
          ];

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

CreateSetAssignment[name_, startIndex_, parameterType_] :=
    Block[{},
          Print["Error: CreateSetAssignment: unknown parameter type: ", ToString[parameterType]];
          Quit[1];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`ScalarType[CConversion`realScalarCType]] :=
    Module[{ass = ""},
           ass = name <> " = v(" <> ToString[startIndex] <> ");\n";
           Return[{ass, 1}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`ScalarType[CConversion`complexScalarCType]] :=
    Module[{ass = "", type},
           type = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
           ass = name <> " = " <> type <> "(v(" <> ToString[startIndex] <>
                 ", v(" <> ToString[startIndex + 1] <> "));\n";
           Return[{ass, 2}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`VectorType[CConversion`realScalarCType, rows_]] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count++,
               ass = ass <> name <> "(" <> ToString[i] <> ") = v(" <>
                     ToString[startIndex + count] <> ");\n";
              ];
           If[rows != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];];
           Return[{ass, rows}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_]] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count++,
                   ass = ass <> name <> "(" <> ToString[i] <> "," <> ToString[j]
                         <> ") = v(" <> ToString[startIndex + count] <> ");\n";
                  ];
              ];
           If[rows * cols != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];];
           Return[{ass, rows * cols}];
          ];

CreateDisplayAssignment[name_, startIndex_, parameterType_] :=
    Block[{},
          Print["Error: CreateDisplayAssignment: unknown parameter type: ",
                ToString[parameterType]];
          Quit[1];
          ];

CreateDisplayAssignment[name_, startIndex_, CConversion`ScalarType[CConversion`realScalarCType]] :=
    Module[{ass = ""},
           ass = "pars(" <> ToString[startIndex] <> ") = "
                 <> name <> ";\n";
           Return[{ass, 1}];
          ];

CreateDisplayAssignment[name_, startIndex_, CConversion`ScalarType[CConversion`complexScalarCType]] :=
    Module[{ass = ""},
           ass = "pars(" <> ToString[startIndex] <> ") = Re(" <> name <> ");\n" <>
                 "pars(" <> ToString[startIndex + 1] <> ") = Im(" <> name <> ");\n";
           Return[{ass, 2}];
          ];

CreateDisplayAssignment[name_, startIndex_, CConversion`VectorType[CConversion`realScalarCType, rows_]] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count++,
               ass = ass <> "pars(" <> ToString[startIndex + count] <> ") = "
                     <> name <> "(" <> ToString[i] <> ");\n";
              ];
           If[rows != count,
              Print["Error: CreateDisplayAssignment: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];];
           Return[{ass, rows}];
          ];

CreateDisplayAssignment[name_, startIndex_, CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_]] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count++,
                   ass = ass <> "pars(" <> ToString[startIndex + count] <> ") = "
                          <> name <> "(" <> ToString[i] <> "," <> ToString[j]
                          <> ");\n";
                  ];
              ];
           If[rows * cols != count,
              Print["Error: CreateDisplayAssignment: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];];
           Return[{ass, rows * cols}];
          ];

CreateParameterNamesStr[name_, parameterType_] :=
    Block[{},
          Print["Error: CreateParameterNamesStr: unknown parameter type: ",
                ToString[parameterType]];
          Quit[1];
          ];

CreateParameterNamesStr[name_, CConversion`ScalarType[CConversion`realScalarCType]] :=
    "\"" <> CConversion`ToValidCSymbolString[name] <> "\"";

CreateParameterNamesStr[name_, CConversion`ScalarType[CConversion`complexScalarCType]] :=
    "\"Re(" <> CConversion`ToValidCSymbolString[name] <>
    "), Im(" <> CConversion`ToValidCSymbolString[name] <> ")\"";

CreateParameterNamesStr[name_, CConversion`VectorType[CConversion`realScalarCType, rows_]] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count++,
               If[ass != "", ass = ass <> ", ";];
               ass = ass <> "\"" <> CConversion`ToValidCSymbolString[name] <>
                     "(" <> ToString[i] <> ")\"";
              ];
           If[rows != count,
              Print["Error: CreateParameterNamesStr: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];
             ];
           Return[ass];
          ];

CreateParameterNamesStr[name_, CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_]] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count++,
                   If[ass != "", ass = ass <> ", ";];
                   ass = ass <> "\"" <> CConversion`ToValidCSymbolString[name] <>
                         "(" <> ToString[i] <> "," <> ToString[j] <> ")\"";
                  ];
              ];
           If[rows * cols != count,
              Print["Error: CreateParameterNamesStr: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];
             ];
           Return[ass];
          ];

CreateParameterEnums[name_, parameterType_] :=
    Block[{},
          Print["Error: CreateParameterEnums: unknown parameter type: ",
                ToString[parameterType]];
          Quit[1];
          ];

CreateParameterEnums[name_, CConversion`ScalarType[CConversion`realScalarCType]] :=
    CConversion`ToValidCSymbolString[name];

CreateParameterEnums[name_, CConversion`ScalarType[CConversion`complexScalarCType]] :=
    CConversion`ToValidCSymbolString[Re[name]] <> ", " <>
    CConversion`ToValidCSymbolString[Im[name]];

CreateParameterEnums[name_, CConversion`VectorType[CConversion`realScalarCType, rows_]] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count++,
               If[ass != "", ass = ass <> ", ";];
               ass = ass <> CConversion`ToValidCSymbolString[name[i]];
              ];
           If[rows != count,
              Print["Error: CreateParameterEnums: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];
             ];
           Return[ass];
          ];

CreateParameterEnums[name_, CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_]] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count++,
                   If[ass != "", ass = ass <> ", ";];
                   ass = ass <> CConversion`ToValidCSymbolString[name[i,j]];
                  ];
              ];
           If[rows * cols != count,
              Print["Error: CreateParameterEnums: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];
             ];
           Return[ass];
          ];

CheckParameter[parameter_] :=
    If[!MemberQ[allModelParameters, parameter] &&
       !MemberQ[allInputParameters, parameter],
       Print["Warning: Trying to set unknown parameter ", parameter];
      ];

SetParameter[parameter_, value_String, class_String] :=
    Module[{parameterStr},
           CheckParameter[parameter];
           parameterStr = CConversion`ToValidCSymbolString[parameter];
           class <> "->set_" <> parameterStr <> "(" <> value <> ");\n"
          ];

SetParameter[parameter_[idx_Integer], value_String, class_String] :=
    Module[{parameterStr},
           CheckParameter[parameter];
           parameterStr = CConversion`ToValidCSymbolString[parameter];
           class <> "->set_" <> parameterStr <> "(" <> ToString[idx] <> ", " <>
           value <> ");\n"
          ];

SetParameter[parameter_[idx1_Integer, idx2_Integer], value_String, class_String] :=
    Module[{parameterStr},
           CheckParameter[parameter];
           parameterStr = CConversion`ToValidCSymbolString[parameter];
           class <> "->set_" <> parameterStr <> "(" <> ToString[idx1] <> ", " <>
           ToString[idx2] <> ", " <> value <> ");\n"
          ];

SetParameter[parameter_, value_, class_String] :=
    SetParameter[parameter, CConversion`RValueToCFormString[value], class];

SaveParameterLocally[parameters_List, prefix_String, caller_String] :=
    Module[{i, result = ""},
           For[i = 1, i <= Length[parameters], i++,
               result = result <> SaveParameterLocally[parameters[[i]], prefix, caller];
              ];
           Return[result];
          ];

SaveParameterLocally[parameter_, prefix_String, caller_String] :=
    Module[{ parStr },
           parStr = CConversion`ToValidCSymbolString[parameter];
           "const auto " <> prefix <> parStr <> " = " <>
           If[caller != "", caller <> "(" <> parStr <> ")", parStr] <> ";\n"
          ];

RestoreParameter[parameters_List, prefix_String, modelPtr_String] :=
    Module[{i, result = ""},
           For[i = 1, i <= Length[parameters], i++,
               result = result <> RestoreParameter[parameters[[i]], prefix, modelPtr];
              ];
           Return[result];
          ];

RestoreParameter[parameter_, prefix_String, modelPtr_String] :=
    Module[{ parStr },
           parStr = CConversion`ToValidCSymbolString[parameter];
           If[modelPtr != "",
              SetParameter[parameter, prefix <> parStr, modelPtr],
              CConversion`ToValidCSymbolString[parameter] <> " = " <> prefix <> parStr <> ";\n"]
          ];

RemoveProtectedHeads[expr_] :=
    expr /. { SARAH`SM[__] -> SARAH`SM[],
              FlexibleSUSY`Pole[__]  -> FlexibleSUSY`Pole[] };

DefineLocalConstCopy[parameter_, macro_String, prefix_String:""] :=
    "const auto " <> prefix <> ToValidCSymbolString[parameter] <> " = " <>
    macro <> "(" <> ToValidCSymbolString[parameter] <> ");\n";

PrivateCallLoopMassFunction[FlexibleSUSY`M[particle_Symbol]] :=
    "calculate_" <> ToValidCSymbolString[FlexibleSUSY`M[particle]] <> "_pole_1loop();\n";

CalculateLocalPoleMasses[parameter_] :=
    "MODEL->" <> PrivateCallLoopMassFunction[parameter];

CreateLocalConstRefs[expr_] :=
    Module[{result = "", symbols, inputSymbols, modelPars, outputPars,
            poleMasses},
           symbols = FindAllParameters[expr];
           poleMasses = {
               Cases[expr, FlexibleSUSY`Pole[FlexibleSUSY`M[a_]]     /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], Infinity],
               Cases[expr, FlexibleSUSY`Pole[FlexibleSUSY`M[a_[__]]] /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], Infinity]
                        };
           symbols = DeleteDuplicates[Flatten[symbols]];
           poleMasses = DeleteDuplicates[Flatten[poleMasses]];
           inputSymbols = DeleteDuplicates[Select[symbols, (MemberQ[allInputParameters,#])&]];
           modelPars    = DeleteDuplicates[Select[symbols, (MemberQ[allModelParameters,#])&]];
           outputPars   = DeleteDuplicates[Select[symbols, (MemberQ[allOutputParameters,#])&]];
           (result = result <> DefineLocalConstCopy[#,"INPUTPARAMETER"])& /@ inputSymbols;
           (result = result <> DefineLocalConstCopy[#,"MODELPARAMETER"])& /@ modelPars;
           (result = result <> DefineLocalConstCopy[#,"MODELPARAMETER"])& /@ outputPars;
           (result = result <> CalculateLocalPoleMasses[#])& /@ poleMasses;
           Return[result];
          ];

CreateLocalConstRefsForPhysicalParameters[expr_] :=
    Module[{result = "", symbols, outputPars, compactExpr},
           compactExpr = RemoveProtectedHeads[expr];
           symbols = { Cases[compactExpr, _Symbol, Infinity],
                       Cases[compactExpr, a_[__] /; MemberQ[allOutputParameters,a] :> a, Infinity],
                       Cases[compactExpr, FlexibleSUSY`M[a_]     /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]], Infinity],
                       Cases[compactExpr, FlexibleSUSY`M[a_[__]] /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], Infinity]
                     };
           symbols = DeleteDuplicates[Flatten[symbols]];
           outputPars = DeleteDuplicates[Select[symbols, (MemberQ[allOutputParameters,#])&]];
           (result = result <> DefineLocalConstCopy[#,"PHYSICAL"])& /@ outputPars;
           Return[result];
          ];

CreateLocalConstRefsForBetas[expr_] :=
    Module[{result = "", symbols, modelPars, compactExpr},
           compactExpr = RemoveProtectedHeads[expr];
           symbols = { Cases[compactExpr, _Symbol, Infinity],
                       Cases[compactExpr, a_[__] /; MemberQ[allModelParameters,a] :> a, Infinity] };
           symbols = DeleteDuplicates[Flatten[symbols]];
           modelPars = DeleteDuplicates[Select[symbols, (MemberQ[allModelParameters,#])&]];
           (result = result <> DefineLocalConstCopy[#, "BETAPARAMETER", "beta_"])& /@ modelPars;
           Return[result];
          ];

CreateLocalConstRefsForInputParameters[expr_, head_String:"INPUT"] :=
    Module[{result = "", symbols, inputPars, compactExpr},
           compactExpr = RemoveProtectedHeads[expr];
           symbols = Cases[compactExpr, _Symbol, Infinity];
           symbols = DeleteDuplicates[Flatten[symbols]];
           inputPars = DeleteDuplicates[Select[symbols, (MemberQ[allInputParameters,#])&]];
           (result = result <> DefineLocalConstCopy[#, head])& /@ inputPars;
           Return[result];
          ];

CreateCaseFromTuple[{key_?NumberQ, parameter_}] :=
    "case " <> ToString[key] <> ": input." <>
    CConversion`ToValidCSymbolString[parameter] <> " = value; break;\n";

CreateCaseFromTuple[expr_] :=
    Block[{},
          Print["Error: not a valid {key,parameter} tuple: ", expr];
          ""
         ];

FillInputParametersFromTuples[minpar_List] :=
    Module[{result = ""},
           (result = result <> CreateCaseFromTuple[#])& /@ minpar;
           result = "switch (key) {\n" <> result <>
                    "default: WARNING(\"Unrecognized key: \" << key); break;\n}\n";
           Return[result];
          ];

DecreaseIndex[ind_Integer] := ind - 1;
DecreaseIndex[ind_]        := ind;
DecreaseIndices[a_[{ind__}]] := a[DecreaseIndex /@ {ind}];
DecreaseIndices[a_[ind__]] := a[Sequence @@ (DecreaseIndex /@ {ind})];
DecreaseIndices[a_]        := a;
DecreaseIndices[SARAH`Delta[a_, b_]] :=
    CConversion`FSKroneckerDelta[DecreaseIndex[a], DecreaseIndex[b]];

DecreaseIndexLiterals[expr_] :=
    DecreaseIndexLiterals[expr, Join[allInputParameters, allModelParameters,
                                     allOutputParameters]];

DecreaseIndexLiterals[expr_, heads_List] :=
    Module[{indexedSymbols, rules, decrExpr, allHeads},
           allHeads = Join[heads, {SARAH`Delta, SARAH`ThetaStep}];
           indexedSymbols = Cases[{expr}, s_[__] /; MemberQ[allHeads, s], Infinity];
           rules = Rule[#, DecreaseIndices[#]] & /@ indexedSymbols;
           decrExpr = expr /. rules;
           Return[decrExpr]
          ];

DecreaseSumIdices[expr_] :=
    expr //. SARAH`sum[idx_, start_, stop_, exp_] :> CConversion`IndexSum[idx, start - 1, stop - 1, exp];

End[];

EndPackage[];


BeginPackage["Constraint`", {"CConversion`", "BetaFunction`", "Parameters`"}];

ApplyConstraints::usage="";
CalculateScale::usage="";
DefineInputParameters::usage="";
InitializeInputParameters::usage="";

SetBetaFunctions::usage=""

Begin["Private`"];

allBetaFunctions = {};

SetBetaFunctions[pars_List] := allBetaFunctions = pars;

GetParameter[parameter_Symbol, macro_String, namePrefix_:""] :=
    Module[{parameterStr, cSymbol, decl},
           parameterStr = ToValidCSymbolString[parameter];
           cSymbol = namePrefix <> ToValidCSymbolString[parameter];
           decl = "const double " <> cSymbol <> " = " <>
                  macro <> "(" <> parameterStr <> ");\n";
           Return[{cSymbol, decl}];
          ];

GetParameter[parameter_[idx_], macro_String, namePrefix_:""] :=
    Module[{parameterStr, cSymbol, idxStr, decl},
           parameterStr = ToValidCSymbolString[parameter];
           cSymbol = namePrefix <> ToValidCSymbolString[parameter[idx]];
           idxStr = ToValidCSymbolString[idx];
           decl = "const double " <> cSymbol <> " = " <>
                  macro <> "2(" <> parameterStr <> "," <> idxStr <> ");\n";
           Return[{cSymbol, decl}];
          ];

GetParameter[parameter_[idx1_,idx2_], macro_String, namePrefix_:""] :=
    Module[{parameterStr, cSymbol, idx1Str, idx2Str, decl},
           parameterStr = ToValidCSymbolString[parameter];
           cSymbol = namePrefix <> ToValidCSymbolString[parameter[idx1,idx2]];
           idx1Str = ToValidCSymbolString[idx1];
           idx2Str = ToValidCSymbolString[idx2];
           decl = "const double " <> cSymbol <> " = " <>
                  macro <> "3(" <> parameterStr <> "," <> idx1Str <>
                  "," <> idx2Str <> ");\n";
           Return[{cSymbol, decl}];
          ];

ApplyConstraint[{parameter_, value_}, modelName_String] :=
    Parameters`SetParameter[parameter, value, modelName];

CreateStartPoint[parameters_List, name_String] :=
    Module[{dim, dimStr, startPoint = "", i},
           dim = Length[parameters];
           dimStr = ToString[dim];
           For[i = 1, i <= dim, i++,
               startPoint = startPoint <> If[i==1," ",", "] <> "MODELPARAMETER(" <>
                            CConversion`ToValidCSymbolString[parameters[[i]]] <> ")";
              ];
           startPoint = "const double " <> name <> "[" <> dimStr <> "] = {" <>
                        startPoint <> " };\n";
           Return[startPoint];
          ];

ApplyConstraint[FlexibleSUSY`FSMinimize[parameters_List, function_], modelName_String] :=
    Module[{callMinimizer, dim, dimStr, startPoint},
           dim = Length[parameters];
           dimStr = ToString[dim];
           startPoint = CreateStartPoint[parameters, "start_point"];
           callMinimizer = startPoint <>
                           "// Minimizer<" <> dimStr <>
                           "> minimizer(func, " <> modelName <> ", 100, 1.0e-2);\n" <>
                           "// const int error = minimizer.minimize(start_point);\n";
           Return[callMinimizer];
          ];

ApplyConstraints[settings_List] :=
    Module[{result, noMacros},
           noMacros = DeleteCases[settings, FlexibleSUSY`FSMinimize[__] | FlexibleSUSY`FSFindRoot[__]];
           result = Parameters`CreateLocalConstRefs[(#[[2]])& /@ noMacros];
           result = result <> "\n";
           (result = result <> ApplyConstraint[#, "MODEL"])& /@ settings;
           Return[result];
          ];

CalculateScale[Null, _] := "";

CalculateScale[False, _] :=
    "ERROR(\"scale condition is allways false!\");\n";

CalculateScale[True, _] :=
    "WARNING(\"scale condition is allways true!\");\n";

CalculateScale[expr_, scaleName_String] :=
    Module[{result},
           result = Parameters`CreateLocalConstRefs[expr];
           result = result <> "\n";
           result = result <> CalculateScaleFromExpr[expr, scaleName];
           Return[result];
          ];

CalculateScale[expr_Equal, scaleName_String] :=
    Module[{result},
           result = Parameters`CreateLocalConstRefs[expr];
           result = result <> Parameters`CreateLocalConstRefsForBetas[expr];
           result = result <> "\n";
           result = result <> CalculateScaleFromExpr[expr, scaleName];
           Return[result];
          ];

GetIndexedParameter[BetaFunction`BetaFunction[name_, CConversion`ScalarType[_], _]] := name;

GetIndexedParameter[BetaFunction`BetaFunction[name_[__], CConversion`MatrixType[_,rows_,cols_], _]] :=
    Flatten[Table[name[i,j], {i,1,rows}, {j,1,cols}]];

GetListOfIndexedParameters[] :=
    Flatten[GetIndexedParameter[#]& /@ allBetaFunctions];

(* Approximately find scale where equation
     expr1 == expr2
   is fulfilled, using the beta function for each parameter p
     p(MX) = p(mu) + BETA[p] log(MX/mu)
*)
CalculateScaleFromExprSymb[Equal[expr1_, expr2_]] :=
    Module[{result, F1, F2, betaFunctions, solution, scale, parameters, parSeq},
           parameters = GetListOfIndexedParameters[];
           parSeq = Sequence @@ parameters;
           F1[parSeq] := expr1;
           F2[parSeq] := expr2;
           betaFunctions = Global`BETA /@ parameters;
           solution = Solve[Log[scale/Global`currentScale] *
                            (betaFunctions . D[F1[parSeq] - F2[parSeq],
                                               {parameters}])
                            == F2[parSeq] - F1[parSeq], scale];
           If[solution === {{}},
              Print["Error: no solution found for ", expr1 == expr2];
              result = Null;,
              result = FullSimplify[solution[[1,1,2]]];
             ];
           Return[result];
          ];

CalculateScaleFromExprSymb[False] := Null;

CalculateScaleFromExprSymb[True] := Global`currentScale;

CalculateScaleFromExpr[Equal[expr1_, expr2_], scaleName_String] :=
    Module[{result, solution},
           solution = CalculateScaleFromExprSymb[Equal[expr1, expr2]];
           If[solution === Null,
              result = "ERROR(\"no solution found for the equation " <>
                        ToString[expr1] <> " == " <> ToString[expr2] <> "\");\n";
              ,
              result = scaleName <> " = " <> RValueToCFormString[solution] <> ";\n";
             ];
           Return[result];
          ];

CalculateScaleFromExpr[expr_, scaleName_String] :=
    scaleName <> " = " <> RValueToCFormString[expr] <> ";\n";

DefineParameter[parameter_Symbol] :=
    "double " <> ToValidCSymbolString[parameter] <> ";\n";

DefineParameter[FlexibleSUSY`Phase[phase_]] :=
    "Complex " <> ToValidCSymbolString[FlexibleSUSY`Phase[phase]] <> ";\n";

DefineParameter[FlexibleSUSY`Sign[phase_]] :=
    "int " <> ToValidCSymbolString[FlexibleSUSY`Sign[phase]] <> ";\n";

DefineInputParameters[inputParameters_List] :=
    Module[{result = ""},
           (result = result <> DefineParameter[#])& /@ inputParameters;
           Return[result];
          ];

InitializeInputParameter[{parameter_, value_?NumericQ}] :=
    ToValidCSymbolString[parameter] <> "(" <> RValueToCFormString[value] <> ")";

InitializeInputParameter[FlexibleSUSY`Phase[phase_]] :=
    ToValidCSymbolString[FlexibleSUSY`Phase[phase]] <> "(1.,.0)";

InitializeInputParameter[FlexibleSUSY`Sign[phase_]] :=
    ToValidCSymbolString[FlexibleSUSY`Sign[phase]] <> "(1)";

InitializeInputParameter[pars__] :=
    Module[{},
           Print["Error: Default values for parameters must be given in the",
                 " form {parameter, value} where value is a number."];
           Return[""];
          ];

InitializeInputParameters[defaultValues_List] :=
    Module[{result = "", i},
           For[i = 1, i <= Length[defaultValues], i++,
               If[i == 1,
                  result = ": ";,
                  result = result <> ", ";
                 ];
               result = result <> InitializeInputParameter[defaultValues[[i]]];
              ];
           Return[result];
          ];

End[];

EndPackage[];


BeginPackage["Constraint`", {"CConversion`", "BetaFunction`"}];

ApplyConstraints::usage="";
CalculateScale::usage="";
DefineInputParameters::usage="";
InitializeInputParameters::usage="";

SetInputParameters::usage="";
SetModelParameters::usage="";
SetOutputParameters::usage="";
SetBetaFunctions::usage=""

CreateLocalConstRefs::usage="creates local const references to model
parameters / input parameters.";

SetParameter::usage="set model parameter";

Begin["Private`"];

allInputParameters = {};
allModelParameters = {};
allOutputParameters = {};
allBetaFunctions = {};

SetInputParameters[pars_List] := allInputParameters = pars;
SetModelParameters[pars_List] := allModelParameters = pars;
SetOutputParameters[pars_List] := allOutputParameters = pars;
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

SetParameter[parameter_, value_String, class_String] :=
    Module[{parameterStr},
           parameterStr = ToValidCSymbolString[parameter];
           class <> "->set_" <> parameterStr <> "(" <> value <> ");\n"
          ];

SetParameter[parameter_, value_, class_String] :=
    SetParameter[parameter, RValueToCFormString[value], class];

ApplyConstraints[settings_List] :=
    Module[{result},
           result = CreateLocalConstRefs[(#[[2]])& /@ settings];
           result = result <> "\n";
           (result = result <> SetParameter[#[[1]], #[[2]], "model"])& /@ settings;
           Return[result];
          ];

RemoveProtectedHeads[expr_] :=
    expr /. SARAH`SM[__] -> SARAH`SM[];

DefineLocalConstCopy[parameter_, macro_String, prefix_String:""] :=
    "const auto " <> prefix <> ToValidCSymbolString[parameter] <> " = " <>
    macro <> "(" <> ToValidCSymbolString[parameter] <> ");\n";

CreateLocalConstRefs[expr_] :=
    Module[{result = "", symbols, inputSymbols, modelPars, outputPars,
            compactExpr},
           compactExpr = RemoveProtectedHeads[expr];
           symbols = { Cases[compactExpr, _Symbol, Infinity],
                       Cases[compactExpr, a_[__] /; MemberQ[allModelParameters,a] :> a, Infinity],
                       Cases[compactExpr, a_[__] /; MemberQ[allOutputParameters,a] :> a, Infinity],
                       Cases[compactExpr, FlexibleSUSY`M[a_]     /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]], Infinity],
                       Cases[compactExpr, FlexibleSUSY`M[a_[__]] /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], Infinity]
                     };
           symbols = DeleteDuplicates[Flatten[symbols]];
           inputSymbols = DeleteDuplicates[Select[symbols, (MemberQ[allInputParameters,#])&]];
           modelPars    = DeleteDuplicates[Select[symbols, (MemberQ[allModelParameters,#])&]];
           outputPars   = DeleteDuplicates[Select[symbols, (MemberQ[allOutputParameters,#])&]];
           (result = result <> DefineLocalConstCopy[#,"INPUTPARAMETER"])& /@ inputSymbols;
           (result = result <> DefineLocalConstCopy[#,"MODELPARAMETER"])& /@ modelPars;
           (result = result <> DefineLocalConstCopy[#,"MODELPARAMETER"])& /@ outputPars;
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

CalculateScale[Null] := "";

CalculateScale[False] :=
    "ERROR(\"GUT scale condition is allways false!\");\n";

CalculateScale[True] :=
    "WARNING(\"GUT scale condition is allways true!\");\n";

CalculateScale[expr_] :=
    Module[{result},
           result = CreateLocalConstRefs[expr];
           result = result <> "\n";
           result = result <> CalculateScaleFromExpr[expr];
           Return[result];
          ];

CalculateScale[expr_Equal] :=
    Module[{result},
           result = CreateLocalConstRefs[expr];
           result = result <> CreateLocalConstRefsForBetas[expr];
           result = result <> "\n";
           result = result <> CalculateScaleFromExpr[expr];
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

CalculateScaleFromExpr[Equal[expr1_, expr2_]] :=
    Module[{result, solution},
           solution = CalculateScaleFromExprSymb[Equal[expr1, expr2]];
           If[solution === Null,
              result = "ERROR(\"no solution found for the equation " <>
                        ToString[expr1] <> " == " <> ToString[expr2] <> "\");\n";
              ,
              result = "scale = " <> RValueToCFormString[solution] <> ";\n";
             ];
           Return[result];
          ];

CalculateScaleFromExpr[expr_] :=
    "scale = " <> RValueToCFormString[expr] <> ";\n";

DefineParameter[parameter_Symbol] :=
    "double " <> ToValidCSymbolString[parameter] <> ";\n";

DefineInputParameters[inputParameters_List] :=
    Module[{result = ""},
           (result = result <> DefineParameter[#])& /@ inputParameters;
           Return[result];
          ];

InitializeInputParameter[{parameter_, value_?NumberQ}] :=
    ToValidCSymbolString[parameter] <> "(" <> RValueToCFormString[value] <> ")";

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

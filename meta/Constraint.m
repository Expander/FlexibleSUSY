
BeginPackage["Constraint`", {"CConversion`", "BetaFunction`", "Parameters`", "TextFormatting`", "TreeMasses`", "Utils`"}];

ApplyConstraints::usage="";
CalculateScale::usage="";
DefineInputParameters::usage="";
InitializeInputParameters::usage="";
InitialGuessAtLowScaleGaugeCouplings::usage="";
IsFixed::usage="returns true if given parameter is fixed in given constraint";

FindFixedParametersFromConstraint::usage="Returns a list of all
parameters which are fixed by the given constraint";

CheckConstraint::usage="Checks a given constraint for syntax errors.";

SanityCheck::usage="Checks general features of the constraints.";

SetBetaFunctions::usage=""

RestrictScale::usage="";

CheckPerturbativityForParameters::usage="";

GetSMMatchingScale::usage="returns SM matching scale from low-energy data set";

SetTemporarily::usage="set temporary variables";

ResetTemporarily::usage="set temporary variables";

Begin["`Private`"];

allBetaFunctions = {};

SetBetaFunctions[pars_List] := allBetaFunctions = pars;

ApplyConstraint[{parameter_, value_}, modelName_String] :=
    Which[Parameters`IsModelParameter[parameter],
          Parameters`SetParameter[parameter, value, modelName <> "->"],
          Parameters`IsInputParameter[parameter],
          Parameters`SetInputParameter[parameter, value, "INPUTPARAMETER"],
          Parameters`IsPhase[parameter],
          Parameters`SetPhase[parameter, value, modelName],
          True,
          Print["Error: ", parameter, " is neither a model nor an input parameter!"];
          ""
         ];

ApplyConstraint[{parameter_ | parameter_[__] /; parameter === SARAH`UpYukawa,
                 value_ /; (!FreeQ[value, Global`upQuarksDRbar] || value === Automatic)},
                modelName_String] :=
    "calculate_" <> CConversion`ToValidCSymbolString[parameter] <> "_DRbar();\n";

ApplyConstraint[{parameter_ | parameter_[__] /; parameter === SARAH`DownYukawa,
                 value_ /; (!FreeQ[value, Global`downQuarksDRbar] || value === Automatic)},
                modelName_String] :=
    "calculate_" <> CConversion`ToValidCSymbolString[parameter] <> "_DRbar();\n";

ApplyConstraint[{parameter_ | parameter_[__] /; parameter === SARAH`ElectronYukawa,
                 value_ /; (!FreeQ[value, Global`downLeptonsDRbar] || value === Automatic)},
                modelName_String] :=
    "calculate_" <> CConversion`ToValidCSymbolString[parameter] <> "_DRbar();\n";

ApplyConstraint[{parameter_,
                 value_ /; !FreeQ[value, Global`neutrinoDRbar]}, modelName_String] :=
    "calculate_MNeutrino_DRbar();\n" <>
    Parameters`SetParameter[parameter, value, modelName <> "->"];

ApplyConstraint[{parameter_ /; !MemberQ[{SARAH`UpYukawa, SARAH`DownYukawa, SARAH`ElectronYukawa}, parameter], value_ /; value === Automatic}, modelName_String] :=
    Block[{},
          Print["Error: cannot determine ", parameter, " automatically!"];
          Quit[1];
         ];

CreateStartPoint[parameters_List, name_String] :=
    Module[{dim, dimStr, startPoint = "", i},
           dim = Length[parameters];
           dimStr = ToString[dim];
           For[i = 1, i <= dim, i++,
               startPoint = startPoint <> If[i==1," ",", "] <> "MODELPARAMETER(" <>
                            CConversion`ToValidCSymbolString[parameters[[i]]] <> ")";
              ];
           "Eigen::VectorXd " <> name <> "(" <> dimStr <> ");\n" <>
           name <> " << " <> startPoint <> " ;\n"
          ];

SetModelParametersFromVector[model_String, vector_String, parameters_List] :=
    Module[{result = "", i, gslElement},
           For[i = 1, i <= Length[parameters], i++,
               gslElement = vector <> "(" <> ToString[i-1] <> ")";
               result = result <> Parameters`SetParameter[parameters[[i]],gslElement,model];
              ];
           Return[result];
          ];

SetVectorFromExpressions[vector_String, expressions_List] :=
    Module[{result = "", i, gslElement},
           For[i = 1, i <= Length[expressions], i++,
               gslElement = vector <> "(" <> ToString[i-1] <> ") = " <>
                            CConversion`RValueToCFormString[expressions[[i]]] <> ";\n";
               result = result <> gslElement;
              ];
           Return[result];
          ];

CreateMinimizationFunctionWrapper[functionName_String, dim_Integer, parameters_List, function_] :=
    Module[{type, stype},
           type  = CConversion`CreateCType[CConversion`MatrixType[CConversion`realScalarCType, dim, 1]];
           stype = CConversion`CreateCType[CConversion`ScalarType[CConversion`realScalarCType]];
"auto " <> functionName <> " = [this](const "<> type <> "& x) -> " <> stype <> " {
" <> TextFormatting`IndentText[SetModelParametersFromVector["MODEL->","x",parameters]] <> "
   MODEL->calculate_DRbar_masses();
" <> TextFormatting`IndentText[Parameters`CreateLocalConstRefs[function]] <> "
   return " <> CConversion`RValueToCFormString[function] <> ";
};
"
          ];

localFunctionWrapper = 0;

CreateSolverName[] := "solver_" <> ToString[localFunctionWrapper++];

ApplyConstraint[FlexibleSUSY`FSMinimize[parameters_List, function_], modelName_String] :=
    Module[{callMinimizer, dim, dimStr, startPoint, functionWrapper, functionName},
           dim = Length[parameters];
           dimStr = ToString[dim];
           startPoint = CreateStartPoint[parameters, "start_point"];
           functionName = CreateSolverName[];
           functionWrapper = CreateMinimizationFunctionWrapper[functionName,dim,parameters,function];
           callMinimizer = functionWrapper <> "\n" <> startPoint <>
                           "Minimizer<" <> dimStr <>
                           "> minimizer(" <> functionName <> ", 100, 1.0e-2);\n" <>
                           "const int status = minimizer.minimize(start_point);\n" <>
                           "VERBOSE_MSG(\"\\tminimizer status: \" << gsl_strerror(status));\n";
           Return[callMinimizer];
          ];

CreateRootFinderFunctionWrapper[functionName_String, dim_Integer, parameters_List, function_List] :=
    Module[{type},
           type = CConversion`CreateCType[CConversion`MatrixType[CConversion`realScalarCType, dim, 1]];
"auto " <> functionName <> " = [this](const "<> type <> "& x) -> " <> type <> " {
" <> TextFormatting`IndentText[SetModelParametersFromVector["MODEL->","x",parameters]] <> "
   MODEL->calculate_DRbar_masses();
" <> TextFormatting`IndentText[Parameters`CreateLocalConstRefs[function]] <> "
   "<> type <> " f;
" <> TextFormatting`IndentText[SetVectorFromExpressions["f",function]] <> "
   return f;
};
"
          ];

ApplyConstraint[FlexibleSUSY`FSFindRoot[parameters_List, function_List], modelName_String] :=
    Module[{callRootFinder, dim, dimStr, startPoint, functionWrapper, functionName},
           dim = Length[parameters];
           dimStr = ToString[dim];
           startPoint = CreateStartPoint[parameters, "start_point"];
           functionName = CreateSolverName[];
           functionWrapper = CreateRootFinderFunctionWrapper[functionName,dim,parameters,function];
           callRootFinder = functionWrapper <> "\n" <> startPoint <>
                           "Root_finder<" <> dimStr <>
                           "> root_finder(" <> functionName <> ", 100, 1.0e-2);\n" <>
                           "const int status = root_finder.find_root(start_point);\n" <>
                           "VERBOSE_MSG(\"\\troot finder status: \" << gsl_strerror(status));\n";
           Return[callRootFinder];
          ];

ApplyConstraint[FlexibleSUSY`FSSolveEWSBFor[___], modelName_String] :=
    modelName <> "->solve_ewsb();\n";

ApplyConstraint[Null, _] :=
    Block[{},
          Print["Error: Null is not a valid constraint setting!"];
          Print["   Maybe there is a trailing comma in the constraint list?"];
          Quit[1];
         ];

ApplyConstraint[p_, _] :=
    Block[{},
          Print["Error: This is not a valid constraint setting: ", p];
          Quit[1];
         ];

ApplyConstraints[settings_List] :=
    Module[{result, noMacros, noTemp},
           noMacros = DeleteCases[
               settings,
               FlexibleSUSY`FSMinimize[__] | \
               FlexibleSUSY`FSFindRoot[__] | \
               FlexibleSUSY`FSSolveEWSBFor[__] | \
               {FlexibleSUSY`Temporary[_], _}
           ];
           noTemp = DeleteCases[
               settings,
               {FlexibleSUSY`Temporary[_], _}
           ];
           result = Parameters`CreateLocalConstRefs[(#[[2]])& /@ noMacros];
           result = result <> AddBetas[noMacros];
           result = result <> "\n";
           (result = result <> ApplyConstraint[#, "MODEL"])& /@ noTemp;
           Return[result];
          ];

ContainsBetas[expr_] := !FreeQ[expr, FlexibleSUSY`BETA];

(* -1 = current beta-function loop order *)
MaxBetaLoopOrder[expr_] :=
    Sort @ Cases[expr /. FlexibleSUSY`BETA[p_] :> FlexibleSUSY`BETA[-1,p],
                 FlexibleSUSY`BETA[l_, _] | FlexibleSUSY`BETA[l_, _][___] :> l, {0, Infinity}];

CallCalcBeta[-1] :=
    "const " <> FlexibleSUSY`FSModelName <> "_soft_parameters beta_functions(MODEL->calc_beta());\n";

CallCalcBeta[l_?IntegerQ] :=
    Module[{lstr = ToString[l]},
           "const " <> FlexibleSUSY`FSModelName <>
           "_soft_parameters beta_functions_" <> lstr <>
           "L(MODEL->calc_beta(" <> lstr <> "));\n"
          ];

AddBetas[expr_?ContainsBetas] :=
    StringJoin[CallCalcBeta /@ MaxBetaLoopOrder[expr]] <>
    Parameters`CreateLocalConstRefsForBetas[expr];

AddBetas[_] := "";

FindFixedParametersFromSetting[{parameter_, value_}] := Parameters`StripIndices[parameter];
FindFixedParametersFromSetting[FlexibleSUSY`FSMinimize[parameters_List, value_]] := parameters;
FindFixedParametersFromSetting[FlexibleSUSY`FSFindRoot[parameters_List, value_]] := parameters;
FindFixedParametersFromSetting[FlexibleSUSY`FSSolveEWSBFor[parameters_List]] := parameters;

FindFixedParametersFromConstraint[settings_List] :=
    DeleteDuplicates[Flatten[FindFixedParametersFromSetting /@ settings]];

CheckSetting[patt:(FlexibleSUSY`FSMinimize|FlexibleSUSY`FSFindRoot)[parameters_, value_],
             constraintName_String] :=
    Module[{modelParameters, unknownParameters},
           modelParameters = Parameters`GetModelParameters[];
           If[Head[parameters] =!= List,
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   First parameter must be a list!"];
              Return[False];
             ];
           If[parameters === {},
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   List of parameters is empty!"];
              Return[False];
             ];
           If[MatchQ[patt, FlexibleSUSY`FSFindRoot[___]],
              If[!MatchQ[patt, FlexibleSUSY`FSFindRoot[_List, _List]],
                 Print["Error: Syntax error in constraint ", constraintName, ": ", InputForm[patt]];
                 Print["   Correct syntax: FSFindRoot[{a,b}, {f[a],f[b]}]"];
                 Return[False];
                 ,
                 If[Length[parameters] =!= Length[value],
                    Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
                    Print["   Argument lists must have the same length!"];
                    Return[False];
                   ];
                ];
             ];
           If[MatchQ[patt, FlexibleSUSY`FSMinimize[_List, _List]],
              Print["Error: Syntax error in constraint ", constraintName, ": ", InputForm[patt]];
              Print["   Correct syntax: FSMinimize[{a,b}, f[a,b]]"];
              Return[False];
             ];
           unknownParameters = Complement[parameters, modelParameters];
           If[unknownParameters =!= {},
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   Unknown parameters: ", unknownParameters];
              Return[False];
             ];
           True
          ];

CheckSetting[patt:FlexibleSUSY`FSSolveEWSBFor[parameters_List],
             constraintName_String] :=
    Module[{unknownParameters = Complement[parameters, Parameters`GetModelParameters[]]},
           If[unknownParameters =!= {},
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt],
                    "   Unknown parameters: ", unknownParameters];
              Return[False];
             ];
           True
          ];

CheckSetting[patt:{parameter_[idx_Integer], value_}, constraintName_String] :=
    Module[{modelParameters, dim},
           modelParameters = Parameters`GetModelParameters[];
           If[!CheckSetting[{parameter, value}],
              Return[False];
             ];
           dim = SARAH`getDimParameters[parameter];
           If[Length[dim] =!= 1,
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   ", parameter, " has ", Length[dim],
                    " dimensions, but one index is accessed!"];
              Return[False];
             ];
           dim = dim[[1]];
           If[1 > idx || idx > dim,
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["    ", parameter, " index out of range!",
                    " Allowed range: 1 ... ", dim];
              Return[False];
             ];
           True
          ];

CheckSetting[patt:{parameter_[idx1_Integer, idx2_Integer], value_}, constraintName_String] :=
    Module[{modelParameters, dim},
           modelParameters = Parameters`GetModelParameters[];
           If[!CheckSetting[{parameter, value}],
              Return[False];
             ];
           dim = SARAH`getDimParameters[parameter];
           If[Length[dim] =!= 2,
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   ", parameter, " has ", Length[dim],
                    " dimensions, but two indices are accessed!"];
              Return[False];
             ];
           If[1 > idx1 || idx1 > dim[[1]],
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   ", parameter, " first index out of range! ",
                    "Allowed range: 1 ... ", dim[[1]]];
              Return[False];
             ];
           If[1 > idx2 || idx2 > dim[[2]],
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   ", parameter, " second index out of range! ",
                    "Allowed range: 1 ... ", dim[[2]]];
              Return[False];
             ];
           True
          ];

CheckSetting[patt:{parameter_, value_}, constraintName_String] :=
    Module[{outputParameters},
           outputParameters = Parameters`GetOutputParameters[];
           If[MemberQ[outputParameters, parameter],
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   ", parameter, " is a output parameter!"];
              Return[False];
             ];
           True
          ];

CheckSetting[patt_, constraintName_String] :=
    Module[{},
           Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
           Print["   This is not a valid constraint setting!"];
           If[patt === Null,
              Print["   Maybe there is a trailing comma in the constraint list?"];
             ];
           False
          ];

CheckConstraint[settings_List, constraintName_String] :=
    CheckSetting[#,constraintName]& /@ settings;

SanityCheck[settings_List, constraintName_String:""] :=
    Module[{setParameters, y,
            yukawas = {SARAH`UpYukawa, SARAH`DownYukawa, SARAH`ElectronYukawa}},
           setParameters = #[[1]]& /@ settings;
           (* check for unset Yukawa couplings *)
           For[y = 1, y <= Length[yukawas], y++,
               If[(ValueQ /@ yukawas)[[y]] &&
                  FreeQ[setParameters, yukawas[[y]]],
                  Print["Warning: Yukawa coupling ", yukawas[[y]],
                        " not set",
                        If[constraintName != "", " in the " <> constraintName, ""],
                        "."];
                 ];
              ];
          ];

CalculateScale[Null, _] := "";

CalculateScale[False, _] :=
    "ERROR(\"scale condition is allways false!\");\n";

CalculateScale[True, _] :=
    "WARNING(\"scale condition is allways true!\");\n";

GetSMMatchingScale[FlexibleSUSY`LowEnergyConstant[FlexibleSUSY`MT], qedqcd_String] :=
    qedqcd <> ".displayPoleMt()";

GetSMMatchingScale[FlexibleSUSY`LowEnergyConstant[FlexibleSUSY`MZ], qedqcd_String] :=
    qedqcd <> ".displayPoleMZ()";

CalculateScale[s:FlexibleSUSY`LowEnergyConstant[_], scaleName_String] :=
    scaleName <> " = " <> GetSMMatchingScale[s, "qedqcd"] <> ";\n";

CalculateScale[expr_, scaleName_String] :=
    Module[{result},
           result = Parameters`CreateLocalConstRefs[expr];
           result = result <> "\n";
           result = result <> CalculateScaleFromExpr[expr, scaleName];
           Return[result];
          ];

CreateBetasForParsIn[expr_] :=
    Module[{pars},
           pars = Parameters`FSModelParameters /. Parameters`FindAllParametersClassified[expr];
           FlexibleSUSY`BETA /@ pars
          ];

CalculateScale[expr_Equal, scaleName_String] :=
    Module[{result},
           result = Parameters`CreateLocalConstRefs[expr];
           result = result <> Parameters`CreateLocalConstRefsForBetas[CreateBetasForParsIn[expr]];
           result = result <> "\n";
           result = result <> CalculateScaleFromExpr[expr, scaleName];
           Return[result];
          ];

GetIndexedParameter[BetaFunction`BetaFunction[name_, CConversion`ScalarType[_], _]] := name;

GetIndexedParameter[BetaFunction`BetaFunction[name_[__], CConversion`MatrixType[_,rows_,cols_], _]] :=
    Flatten[Table[name[i,j], {i,1,rows}, {j,1,cols}]];

GetListOfIndexedParameters[] :=
    Flatten[GetIndexedParameter[#]& /@ allBetaFunctions];

StripConditionalExpression[ConditionalExpression[expr_, ___]] :=
    expr;

StripConditionalExpression[expr_] :=
    expr;

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
           betaFunctions = FlexibleSUSY`BETA /@ parameters;
           solution = Solve[Log[scale/Global`currentScale] *
                            (betaFunctions . D[F1[parSeq] - F2[parSeq],
                                               {parameters}])
                            == F2[parSeq] - F1[parSeq], scale];
           If[solution === {{}},
              Print["Error: no solution found for ", expr1 == expr2];
              result = Null;,
              result = FullSimplify[StripConditionalExpression[solution[[1,1,2]]]];
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
    scaleName <> " = " <> CConversion`RValueToCFormString[Parameters`DecreaseIndexLiterals[expr, Parameters`GetOutputParameters[]]] <> ";\n";

DefineParameter[{parameter_, type_}] :=
    CConversion`CreateCType[type] <> " " <> CConversion`ToValidCSymbolString[parameter] <> ";\n";

DefineInputParameters[inputParameters_List] :=
    Module[{result = ""},
           (result = result <> DefineParameter[#])& /@ inputParameters;
           Return[result];
          ];

InitializeInputParameter[{FlexibleSUSY`Phase[phase_], _}] :=
    ToValidCSymbolString[FlexibleSUSY`Phase[phase]] <> "(1.,.0)";

InitializeInputParameter[{Sign[phase_], _}] :=
    ToValidCSymbolString[Sign[phase]] <> "(1)";

InitializeInputParameter[{parameter_, type_}] :=
    CConversion`CreateDefaultConstructor[CConversion`ToValidCSymbolString[parameter],type];

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

InitialGuessAtLowScaleGaugeCouplings[] :=
    Module[{result = ""},
           If[ValueQ[SARAH`hyperchargeCoupling],
              result = result <> Parameters`SetParameter[SARAH`hyperchargeCoupling,
                                                         "Sqrt(4. * Pi * alpha_sm(0))", "MODEL->"];
             ];
           If[ValueQ[SARAH`leftCoupling],
              result = result <> Parameters`SetParameter[SARAH`leftCoupling,
                                                         "Sqrt(4. * Pi * alpha_sm(1))", "MODEL->"];
             ];
           If[ValueQ[SARAH`strongCoupling],
              result = result <> Parameters`SetParameter[SARAH`strongCoupling,
                                                         "Sqrt(4. * Pi * alpha_sm(2))", "MODEL->"];
             ];
           result
          ];

IsFixedIn[par_, {p_, _}] :=
    Parameters`StripIndices[par] === Parameters`StripIndices[p];

IsFixedIn[par_, FlexibleSUSY`FSMinimize[parameters_List, _]] :=
    MemberQ[Parameters`StripIndices /@ parameters, Parameters`StripIndices[par]];

IsFixedIn[par_, FlexibleSUSY`FSFindRoot[parameters_List, _]] :=
    MemberQ[Parameters`StripIndices /@ parameters, Parameters`StripIndices[par]];

IsFixedIn[par_, FlexibleSUSY`FSSolveEWSBFor[parameters___]] :=
    MemberQ[Parameters`StripIndices /@ Flatten[{parameters}], Parameters`StripIndices[par]];

IsFixedIn[par_, p___] :=
    Block[{},
          Print["Error: This is not a valid constraint setting: ", p];
          Quit[1];
         ];

IsFixed[par_, constraint_List] :=
    Or @@ (IsFixedIn[par, #]& /@ constraint);

RestrictScale[{minimumScale_, maximumScale_}, scaleName_String:"scale"] :=
    Module[{result = "", value},
           If[NumericQ[minimumScale],
              value = CConversion`RValueToCFormString[minimumScale];
              result = result <>
                       "\
if (" <> scaleName <> " < " <> value <> ") {
#ifdef ENABLE_VERBOSE
" <> IndentText["WARNING(\"" <> scaleName <> " < " <> value <> "\");"] <> "
#endif
" <> IndentText[scaleName <> " = " <> value] <> ";
}
";
             ];
           If[NumericQ[maximumScale],
              value = CConversion`RValueToCFormString[maximumScale];
              result = result <>
                       "\
if (" <> scaleName <> " > " <> value <> ") {
#ifdef ENABLE_VERBOSE
" <> IndentText["WARNING(\"" <> scaleName <> " > " <> value <> "\");"] <> "
#endif
" <> IndentText[scaleName <> " = " <> value] <> ";
}
";
             ];
           Return[result];
          ];

PerturbativityCheckSnippet[par_, thresh_, model_String, problem_String] :=
    Module[{parStr, parEnum, threshStr},
           parStr = CConversion`RValueToCFormString[par];
           parEnum = FlexibleSUSY`FSModelName <> "_info::" <> Parameters`CreateEnumName[par];
           threshStr = CConversion`RValueToCFormString[thresh];
           "\
if (MaxAbsValue(" <> parStr <> ") > " <> threshStr <> ") {
   " <> problem <> " = true;
   " <> model <> "->get_problems().flag_non_perturbative_parameter(" <> parEnum <> ", MaxAbsValue(" <> parStr <> "), " <> model <> "->get_scale(), " <> threshStr <> ");
} else {
   " <> model <> "->get_problems().unflag_non_perturbative_parameter(" <> parEnum <> ");
}
"
          ];

CheckPerturbativityForParameter[par_, thresh_, model_String:"model", problem_String:"problem"] :=
    Module[{snippet, parStr, threshStr},
           Utils`StringJoinWithSeparator[
               PerturbativityCheckSnippet[#, thresh, model, problem]& /@ Parameters`DecomposeParameter[par, GetType[par]], "\n"]
          ];

CheckPerturbativityForParameters[pars_List, thresh_] :=
    Parameters`CreateLocalConstRefs[pars] <> "\n" <>
    StringJoin[CheckPerturbativityForParameter[#,thresh]& /@ pars];

SaveValue[par_[idx__], prefix_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par],
            parStrIdx = CConversion`ToValidCSymbolString[par[idx]]},
           "const auto " <> prefix <> parStrIdx <> " = MODELPARAMETER(" <> parStr <> ")" <>
           "(" <> Utils`StringJoinWithSeparator[ToString /@ {idx},","] <> ");"
          ];

SaveValue[par_, prefix_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "const auto " <> prefix <> parStr <> " = MODELPARAMETER(" <> parStr <> ");"
          ];

RestoreValue[par_, prefix_String] :=
    Parameters`SetParameter[
        par,
        prefix <> CConversion`ToValidCSymbolString[par],
        "MODEL->"];

SetTemporarily[settings_List] :=
    Module[{tempSettings = Cases[settings, {FlexibleSUSY`Temporary[p_], v_} :> {p,v}],
            set, savedVals},
           If[tempSettings === {}, Return[""];];
           set = ApplyConstraints[tempSettings];
           savedVals = Utils`StringJoinWithSeparator[SaveValue[#[[1]], "old_"]& /@ tempSettings, "\n"];
           "// temporary parameter re-definitons\n" <>
           savedVals <> "\n{\n" <> IndentText[set] <> "}"
          ];

ResetTemporarily[settings_List] :=
    Module[{tempSettings = Cases[settings, {FlexibleSUSY`Temporary[p_], v_} :> {p,v}]},
           If[tempSettings === {}, Return[""];];
           "// reset temporary parameter re-definitons\n" <>
           StringJoin[RestoreValue[#[[1]], "old_"]& /@ tempSettings]
          ];

End[];

EndPackage[];

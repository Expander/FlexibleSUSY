
BeginPackage["BetaFunction`", {"SARAH`", "TextFormatting`", "CConversion`", "Parameters`", "Traces`"}];

BetaFunction[];

ConvertSarahRGEs::usage="converts SARAH's beta functions";
CreateBetaFunction::usage="";
GetAllBetaFunctions::usage="";
CreateBetaFunction::usage="";
CountNumberOfParameters::usage="";
CreateDisplayFunction::usage="";
CreateSetFunction::usage="";
CreateSetters::usage="";
CreateGetters::usage="";
CreateParameterDefinitions::usage="";
CreateParameterDefaultInitialization::usage="";
CreateCCtorInitialization::usage="";
CreateCCtorParameterList::usage="";
CreateParameterList::usage="";
ClearParameters::usage="";

CreateParameterEnum::usage="";
CreateParameterNames::usage="";

GetName::usage="returns parameter name from beta function";

CreateSingleBetaFunctionDecl::usage="";
CreateSingleBetaFunctionDefs::usage="";

Begin["`Private`"];

GetName[BetaFunction[name_, type_, beta_List]] := name;

GetType[BetaFunction[name_, type_, beta_List]] := type;

GetBeta1Loop[BetaFunction[name_, type_, beta_List]] := beta[[1]];

GetBeta2Loop[BetaFunction[name_, type_, beta_List]] := beta[[2]];

GetAllBetaFunctions[BetaFunction[name_, type_, beta_List]] := beta;

GuessType[sym_[Susyno`LieGroups`i1, SARAH`i2]] :=
    Parameters`GetType[sym];

GuessType[sym_[Susyno`LieGroups`i1]] :=
    Parameters`GetType[sym];

GuessType[sym_] :=
    Parameters`GetType[sym];

CreateSingleBetaFunctionDecl[betaFun_List] :=
    Module[{result = ""},
           (result = result <> CConversion`CreateCType[GetType[#]] <>
                     " calc_beta_" <> CConversion`ToValidCSymbolString[GetName[#]] <>
                     "_one_loop(const TRACE_STRUCT_TYPE&) const;\n" <>
                     CConversion`CreateCType[GetType[#]] <>
                     " calc_beta_" <> CConversion`ToValidCSymbolString[GetName[#]] <>
                     "_two_loop(const TRACE_STRUCT_TYPE&) const;\n";)& /@ betaFun;
           Return[result];
          ];

CreateSingleBetaFunctionDefs[betaFun_List, templateFile_String, sarahTraces_List] :=
    Module[{b, para, type, paraStr, typeStr, files = {},
            inputFile, outputFile,
            localDeclOneLoop, localDeclTwoLoop, betaOneLoop, betaTwoLoop},
           For[b = 1, b <= Length[betaFun], b++,
               para = GetName[betaFun[[b]]];
               type = GetType[betaFun[[b]]];
               paraStr = CConversion`ToValidCSymbolString[para];
               typeStr = CConversion`CreateCType[type];
               inputFile  = FileNameJoin[{Global`$flexiblesusyTemplateDir, templateFile}];
               outputFile = FileNameJoin[{Global`$flexiblesusyOutputDir,
                                          FlexibleSUSY`FSModelName <> "_" <>
                                          StringReplace[templateFile,
                                                        {".cpp.in" -> paraStr <> ".cpp"}]}];
               {localDeclOneLoop, betaOneLoop} = CreateBetaFunction[betaFun[[b]], 1, sarahTraces];
               {localDeclTwoLoop, betaTwoLoop} = CreateBetaFunction[betaFun[[b]], 2, sarahTraces];
               WriteOut`ReplaceInFiles[{{inputFile, outputFile}},
                     { "@ModelName@"     -> FlexibleSUSY`FSModelName,
                       "@parameterType@" -> typeStr,
                       "@parameterName@" -> paraStr,
                       "@localDeclOneLoop@" -> WrapLines[IndentText[localDeclOneLoop]],
                       "@localDeclTwoLoop@" -> WrapLines[IndentText[localDeclTwoLoop]],
                       "@betaOneLoop@"   -> WrapLines[IndentText[betaOneLoop]],
                       "@betaTwoLoop@"   -> WrapLines[IndentText[betaTwoLoop]],
                       "@DateAndTime@"   -> DateString[]
                     } ];
               AppendTo[files, outputFile];
              ];
           Return[files];
          ];

(*
 * Create one-loop and two-loop beta function assignments and local definitions.
 *)
CreateBetaFunction[betaFunction_BetaFunction, loopOrder_Integer, sarahTraces_List] :=
     Module[{beta, betaName, name, betaStr, dataType, unitMatrix,
             type = ErrorType, localDecl, traceRules, expr, loopFactor},
            Switch[loopOrder,
                   1, loopFactor = CConversion`oneOver16PiSqr;,
                   2, loopFactor = CConversion`twoLoop;,
                   _, loopFactor = CConversion`oneOver16PiSqr^loopOrder];
            type = GetType[betaFunction];
            expr = GetAllBetaFunctions[betaFunction];
            If[loopOrder > Length[expr],
               Print["Warning: there is no beta function of ", loopOrder, " loop order."];
               Return[{"", RValueToCFormString[CConversion`CreateZero[type]]}];
              ];
            expr       = expr[[loopOrder]];
            dataType   = CConversion`CreateCType[type];
            unitMatrix = CreateUnitMatrix[type];
            (* convert beta function expressions to C form *)
            name       = ToValidCSymbolString[GetName[betaFunction]];
            betaName   = "beta_" <> name;
            beta       = (loopFactor * expr) /.
                            { Kronecker[Susyno`LieGroups`i1,SARAH`i2] -> unitMatrix,
                              a_[Susyno`LieGroups`i1] :> a,
                              a_[Susyno`LieGroups`i1,SARAH`i2] :> a };
            {localDecl, traceRules} = Traces`CreateDoubleTraceAbbrs[{expr}];
            localDecl  = Traces`CreateLocalCopiesOfTraces[{expr}, "TRACE_STRUCT"];
            beta = beta /. traceRules;
            (* replace SARAH traces in expr *)
            traceRules = Rule[#,ToValidCSymbol[#]]& /@ (Traces`FindSARAHTraces[expr, sarahTraces]);
            beta = beta /. traceRules;
            (* declare SARAH traces locally *)
            localDecl  = localDecl <> Traces`CreateLocalCopiesOfSARAHTraces[expr, sarahTraces, "TRACE_STRUCT"];
            If[beta === 0,
               beta = CConversion`CreateZero[type];
              ];
            betaStr    = RValueToCFormString[beta];
            betaStr    = betaName <> " = " <> betaStr <> ";\n";
            localDecl  = Parameters`CreateLocalConstRefsForInputParameters[expr] <>
                         localDecl;
            Return[{localDecl, betaStr}];
           ];

CreateBetaFunctionCall[betaFunction_BetaFunction] :=
     Module[{beta1L, beta2L = "", betaName = "", name = "",
             oneLoopBetaStr, dataType, localDecl = "",
             twoLoopBeta, twoLoopBetaStr, type = ErrorType},
            name          = ToValidCSymbolString[GetName[betaFunction]];
            dataType      = CConversion`CreateCType[GetType[betaFunction]];
            betaName      = "beta_" <> name;
            oneLoopBetaStr = "calc_beta_" <> name <> "_one_loop(TRACE_STRUCT)";
            beta1L        = dataType <> " " <> betaName <> "(" <> oneLoopBetaStr <> ");\n";
           If[Length[GetAllBetaFunctions[betaFunction]] > 1,
              twoLoopBetaStr = "calc_beta_" <> name <> "_two_loop(TRACE_STRUCT)";
              beta2L = beta2L <> betaName <> " += " <> twoLoopBetaStr <> ";\n";
             ];
            Return[{localDecl, beta1L, beta2L}];
          ];

CreateBetaFunction[betaFunctions_List, sarahTraces_List] :=
    Module[{def = "",
            localDecl = "", beta1L = "", beta2L = "", allDecl = "", allBeta = "",
            allBeta1L = "", allBeta2L = "", i, inputParsDecl},
           For[i = 1, i <= Length[betaFunctions], i++,
               {localDecl, beta1L, beta2L} = CreateBetaFunctionCall[betaFunctions[[i]]];
               allDecl = allDecl <> localDecl;
               allBeta1L = allBeta1L <> beta1L;
               allBeta2L = allBeta2L <> "   " <> beta2L;
              ];
           allBeta = allDecl <> "\n" <> allBeta1L <> "\nif (get_loops() > 1) {\n" <>
                     allBeta2L <> "\n}\n";
           Return[allBeta];
          ];

(* Converts SARAH beta functions to our own format.
 *
 * SARAH format:
 *   { name, one-loop beta, two-loop beta }
 *
 * Our format:
 *   BetaFunction[name, type, {one-loop beta, two-loop beta}]
 *
 * @param betaFunctions list of SARAH-like formated beta functions
 *)
ConvertSarahRGEs[betaFunctions_List] :=
    Module[{lst = {}, beta, i, k, name, type, expr},
           For[i = 1, i <= Length[betaFunctions], i++,
               beta = betaFunctions[[i]];
               (* extract all beta functions and guess type *)
               For[k = 1, k <= Length[beta], k++,
                   If[Length[beta[[k]] < 2], Continue[];];
                   (* beta[[k,1]] == name, beta[[k,2]] == 1 loop beta function *)
                   name = beta[[k,1]];
                   type = GuessType[name];
                   expr = Drop[beta[[k]], 1];
                   (* protect tensor products *)
                   expr = CConversion`ProtectTensorProducts[#, name]& /@ expr;
                   (* simplify expressions *)
                   expr = TimeConstrained[Simplify[#], FlexibleSUSY`FSSimplifyBetaFunctionsTimeConstraint, #]& /@ expr;
                   AppendTo[lst, BetaFunction[name, type, expr]];
                  ];
              ];
           Return[lst];
          ];

(* count number of parameters in beta functions list *)
CountNumberOfParameters[CConversion`ScalarType[CConversion`realScalarCType]] := 1;
CountNumberOfParameters[CConversion`ScalarType[CConversion`complexScalarCType]] := 2;
CountNumberOfParameters[CConversion`ArrayType[CConversion`realScalarCType, entries_]] := entries;
CountNumberOfParameters[CConversion`ArrayType[CConversion`complexScalarCType, entries_]] := 2 * entries;
CountNumberOfParameters[CConversion`VectorType[CConversion`realScalarCType, entries_]] := entries;
CountNumberOfParameters[CConversion`VectorType[CConversion`complexScalarCType, entries_]] := 2 * entries;
CountNumberOfParameters[CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_]] := rows * cols;
CountNumberOfParameters[CConversion`MatrixType[CConversion`complexScalarCType, rows_, cols_]] := rows * cols;

CountNumberOfParameters[betaFunctions_List] :=
    Module[{num = 0},
           (num += CountNumberOfParameters[GetType[#]])& /@ betaFunctions;
           Return[num];
          ];

(* creating set function *)
CreateSetFunction[betaFunctions_List, parameterNumberOffset_:0] :=
    Module[{set = "", paramCount = parameterNumberOffset, name = "", beta = {},
            type = ErrorType, i, numberOfParameters, assignment = "",
            nAssignments = 0},
           For[i = 1, i <= Length[betaFunctions], i++,
               beta = GetAllBetaFunctions[betaFunctions[[i]]];
               type = GetType[betaFunctions[[i]]];
               name = ToValidCSymbolString[GetName[betaFunctions[[i]]]];
               {assignment, nAssignments} = Parameters`CreateSetAssignment[name, paramCount, type];
               set = set <> assignment;
               paramCount += nAssignments;
              ];
           (* sanity check *)
           numberOfParameters = CountNumberOfParameters[betaFunctions] + parameterNumberOffset;
           If[paramCount != numberOfParameters,
              Print["Error: CreateSetFunction: number of parameters does not match: ", paramCount,
                    " != ", numberOfParameters]; Quit[1];];
           Return[set];
          ];

(* creating display function *)
CreateDisplayFunction[betaFunctions_List, parameterNumberOffset_:0] :=
    Module[{display = "", paramCount = parameterNumberOffset, name = "",
            beta = {}, type = ErrorType, i, numberOfParameters, assignment = "",
            nAssignments = 0},
           numberOfParameters = CountNumberOfParameters[betaFunctions] + parameterNumberOffset;
           For[i = 1, i <= Length[betaFunctions], i++,
               beta = GetAllBetaFunctions[betaFunctions[[i]]];
               type = GetType[betaFunctions[[i]]];
               name = ToValidCSymbolString[GetName[betaFunctions[[i]]]];
               {assignment, nAssignments} = Parameters`CreateDisplayAssignment[name, paramCount, type];
               display = display <> assignment;
               paramCount += nAssignments;
              ];
           (* sanity check *)
           If[paramCount != numberOfParameters,
              Print["Error: CreateDisplayFunction: number of parameters does not match: ", paramCount,
                    " != ", numberOfParameters]; Quit[1];];
           Return[display];
          ];

CreateParameterNames[betaFunctions_List] :=
    Module[{i, par, type, name, result = ""},
           For[i = 1, i <= Length[betaFunctions], i++,
               par = GetName[betaFunctions[[i]]];
               type = GetType[betaFunctions[[i]]];
               name = Parameters`CreateParameterNamesStr[par, type];
               If[i > 1, result = result <> ", ";];
               result = result <> name;
              ];
           result = "const char* parameter_names[NUMBER_OF_PARAMETERS] = {" <>
                    result <> "};\n";
           Return[result];
          ];

CreateParameterEnum[betaFunctions_List] :=
    Module[{i, par, type, name, result = ""},
           For[i = 1, i <= Length[betaFunctions], i++,
               par = GetName[betaFunctions[[i]]];
               type = GetType[betaFunctions[[i]]];
               name = Parameters`CreateParameterEnums[par, type];
               If[i > 1, result = result <> ", ";];
               result = result <> name;
              ];
           (* append enum state for the number of betaFunctions *)
           If[Length[betaFunctions] > 0, result = result <> ", ";];
           result = result <> "NUMBER_OF_PARAMETERS";
           result = "enum Parameters : unsigned {" <>
                    result <> "};\n";
           Return[result];
          ];

(* create setters *)
CreateElementSetter[name_String, CConversion`ScalarType[_]] := "";

CreateElementSetter[name_String, type_] :=
    CConversion`CreateInlineElementSetter[name, type];

CreateSetters[betaFunction_BetaFunction] :=
    Module[{setter = "", name = "", type},
           name = ToValidCSymbolString[GetName[betaFunction]];
           type = GetType[betaFunction];
           setter = setter <> CConversion`CreateInlineSetter[name, type];
           setter = setter <> CreateElementSetter[name, type];
           Return[setter];
          ];

CreateSetters[betaFunctions_List] :=
    Module[{setter = ""},
           (setter = setter <> CreateSetters[#])& /@ betaFunctions;
           Return[setter];
          ];

(* create getters *)
CreateElementGetter[name_String, CConversion`ScalarType[_]] := "";

CreateElementGetter[name_String, type_] :=
    CConversion`CreateInlineElementGetter[name, type];

CreateGetters[betaFunction_BetaFunction] :=
    Module[{getter = "", name = "", type},
           name = ToValidCSymbolString[GetName[betaFunction]];
           type = GetType[betaFunction];
           getter = getter <> CConversion`CreateInlineGetter[name, type];
           getter = getter <> CreateElementGetter[name, type];
           Return[getter];
          ];

CreateGetters[betaFunctions_List] :=
    Module[{getter = ""},
           (getter = getter <> CreateGetters[#])& /@ betaFunctions;
           Return[getter];
          ];

(* create parameter definition in C++ class *)
CreateParameterDefinitions[betaFunction_BetaFunction] :=
    Module[{def = "", name = "", dataType = ""},
           dataType = CConversion`CreateCType[GetType[betaFunction]];
           name = ToValidCSymbolString[GetName[betaFunction]];
           def  = def <> dataType <> " " <> name <> ";\n";
           Return[def];
          ];

CreateParameterDefinitions[betaFunctions_List] :=
    Module[{def = ""},
           (def = def <> CreateParameterDefinitions[#])& /@ betaFunctions;
           Return[def];
          ];

(* create parameter default initialization list *)
CreateParameterDefaultInitialization[betaFunction_BetaFunction] :=
    Module[{def = "", name = "", dataType = ""},
           dataType = CConversion`CreateCType[GetType[betaFunction]];
           name = ToValidCSymbolString[GetName[betaFunction]];
           def  = def <> ", " <> CreateDefaultConstructor[name, GetType[betaFunction]];
           Return[def];
          ];

CreateParameterDefaultInitialization[betaFunctions_List] :=
    Module[{def = ""},
           (def = def <> CreateParameterDefaultInitialization[#])& /@ betaFunctions;
           Return[def];
          ];

(* create copy constructor initialization list *)
CreateCCtorInitialization[betaFunction_BetaFunction] :=
    Module[{def = "", name = "", dataType = ""},
           dataType = CConversion`CreateCType[GetType[betaFunction]];
           name = ToValidCSymbolString[GetName[betaFunction]];
           def  = def <> ", " <> name <> "(" <> name <> "_)";
           Return[def];
          ];

CreateCCtorInitialization[betaFunctions_List] :=
    Module[{def = ""},
           (def = def <> CreateCCtorInitialization[#])& /@ betaFunctions;
           Return[def];
          ];

(* create copy constructor initialization list *)
CreateCCtorParameterList[betaFunction_BetaFunction] :=
    Module[{def = "", name = "", dataType = ""},
           dataType = CreateGetterReturnType[GetType[betaFunction]];
           name = ToValidCSymbolString[GetName[betaFunction]];
           def = def <> ", " <> dataType <> " " <> name <> "_";
           Return[def];
          ];

CreateCCtorParameterList[betaFunctions_List] :=
    Module[{def = "", i},
           For[i = 1, i <= Length[betaFunctions], i++,
               def = def <> CreateCCtorParameterList[betaFunctions[[i]]];
              ];
           Return[def];
          ];

(* create parameter list *)
CreateParameterList[betaFunction_BetaFunction, prefix_] :=
    Module[{def = "", name = "", dataType = ""},
           dataType = CConversion`CreateCType[GetType[betaFunction]];
           name = ToValidCSymbolString[GetName[betaFunction]];
           If[def != "", def = def <> ", "];
           def = def <> prefix <> name;
           Return[def];
          ];

CreateParameterList[betaFunctions_List, prefix_:""] :=
    Module[{def = "", i},
           For[i = 1, i <= Length[betaFunctions], i++,
               If[def != "", def = def <> ", "];
               def = def <> CreateParameterList[betaFunctions[[i]], prefix];
              ];
           Return[def];
          ];

ClearParameter[betaFunction_BetaFunction] :=
    Module[{def, type},
           def = CConversion`SetToDefault[ToValidCSymbolString[GetName[betaFunction]],
                                          GetType[betaFunction]];
           Return[def];
          ];

ClearParameters[betaFunctions_List] :=
    Module[{def = ""},
           (def = def <> ClearParameter[#])& /@ betaFunctions;
           Return[def];
          ];

End[];

EndPackage[];

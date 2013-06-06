
BeginPackage["BetaFunction`", {"SARAH`", "TextFormatting`", "CConversion`", "Parameters`"}];

BetaFunction[];

ConvertSarahRGEs::usage="converts SARAH's beta functions";
CreateBetaFunction::usage="";
GetAllBetaFunctions::usage="";
CreateBetaFunction::usage="";
CountNumberOfParameters::usage="";
CreateDisplayFunction::usage="";
CreateSetFunction::usage="";
ConvertParameterNames::usage="";
CreateSetters::usage="";
CreateGetters::usage="";
CreateParameterDefinitions::usage="";
CreateParameterDefaultInitialization::usage="";
CreateCCtorInitialization::usage="";
CreateCCtorParameterList::usage="";
CreateParameterList::usage="";

Begin["Private`"];

GetName[BetaFunction[name_, type_, beta_List]] := name;

GetType[BetaFunction[name_, type_, beta_List]] := type;

GetBeta1Loop[BetaFunction[name_, type_, beta_List]] := beta[[1]];

GetBeta2Loop[BetaFunction[name_, type_, beta_List]] := beta[[2]];

GetAllBetaFunctions[BetaFunction[name_, type_, beta_List]] := beta;

GuessType[sym_[i1,i2]] := CConversion`MatrixType["DoubleMatrix", Sequence @@ SARAH`getDimParameters[sym]];

GuessType[sym_Symbol] := CConversion`ScalarType["double"];

(*
 * Create one-loop and two-loop beta function assignments and local definitions.
 *)
CreateBetaFunction[betaFunction_BetaFunction] :=
     Module[{beta1L = "", beta2L = "", betaName = "", name = "",
             oneLoopBeta, localDecl = "", dataType, unitMatrix,
             twoLoopBeta = "", type = ErrorType},
            type = GetType[betaFunction];
            dataType = GetCParameterType[type];
            unitMatrix = CreateUnitMatrix[type];
           (* convert beta function expressions to C form *)
           name          = ToValidCSymbolString[GetName[betaFunction]];
           betaName     := "beta_" <> name;
           oneLoopBeta  := RValueToCFormString[(CConversion`oneOver16PiSqr * GetBeta1Loop[betaFunction])
                                         /. { Kronecker[i1,i2] -> unitMatrix }
                                         /. { a_[i1,i2] :> a }];
           beta1L        = beta1L <> betaName <> " = " <> oneLoopBeta <> ";\n";
           If[Length[GetAllBetaFunctions[betaFunction]] > 1,
              twoLoopBeta  := RValueToCFormString[(CConversion`twoLoop * GetBeta2Loop[betaFunction])
                                            /. { Kronecker[i1,i2] -> unitMatrix }
                                            /. { a_[i1,i2] :> a }];
              beta2L        = beta2L <> betaName <> " = " <> betaName <>
                              " + " <> twoLoopBeta <> ";\n";
              ];
           localDecl     = localDecl <> dataType <> " "
                           <> CreateDefaultConstructor[betaName, type] <> ";\n";
           Return[{localDecl, beta1L, beta2L}];
          ];

CreateBetaFunction[betaFunctions_List, additionalDecl_String] :=
    Module[{def = "",
            localDecl = "", beta1L = "", beta2L = "", allDecl = "", allBeta = "",
            allBeta1L = "", allBeta2L = "", i},
           allDecl = "const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);\n" <>
                     "const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;\n";
           allDecl = allDecl <> "\n" <> additionalDecl <> "\n";
           For[i = 1, i <= Length[betaFunctions], i++,
               {localDecl, beta1L, beta2L} = CreateBetaFunction[betaFunctions[[i]]];
               allDecl = allDecl <> localDecl;
               allBeta1L = allBeta1L <> beta1L;
               allBeta2L = allBeta2L <> "   " <> beta2L;
              ];
           allBeta = allDecl <> "\n" <> allBeta1L <> "\nif (displayLoops() > 1) {\n" <>
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
    Module[{lst = {}, beta, i, k, name, type},
           For[i = 1, i <= Length[betaFunctions], i++,
               beta = betaFunctions[[i]];
               (* extract all beta functions and guess type *)
               For[k = 1, k <= Length[beta], k++,
                   If[Length[beta[[k]] < 2], Continue[];];
                   (* beta[[k,1]] == name, beta[[k,2]] == 1 loop beta function *)
                   name = beta[[k,1]];
                   type = GuessType[name];
                   AppendTo[lst, BetaFunction[name, type, Drop[beta[[k]], 1]]];
                  ];
              ];
           Return[lst];
          ];

(* count number of parameters in beta functions list *)
CountNumberOfParameters[CConversion`ScalarType[type_]] := 1;

CountNumberOfParameters[CConversion`VectorType[type_, entries_]] := entries;

CountNumberOfParameters[CConversion`MatrixType[type_, rows_, cols_]] := rows * cols;

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
               (* start counting with 1 *)
               {assignment, nAssignments} = Parameters`CreateSetAssignment[name, paramCount + 1, type];
               set = set <> assignment;
               paramCount += nAssignments;
              ];
           (* sanity check *)
           numberOfParameters = CountNumberOfParameters[betaFunctions] + parameterNumberOffset;
           If[paramCount != numberOfParameters,
              Print["Error: CreateSetFunction: number of parameters does not match: ", paramCount,
                    " != ", numberOfParameters]; Quit[];];
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
               (* start counting with 1 *)
               {assignment, nAssignments} = Parameters`CreateDisplayAssignment[name, paramCount + 1, type];
               display = display <> assignment;
               paramCount += nAssignments;
              ];
           (* sanity check *)
           If[paramCount != numberOfParameters,
              Print["Error: CreateDisplayFunction: number of parameters does not match: ", paramCount,
                    " != ", numberOfParameters]; Quit[];];
           Return[display];
          ];

ConvertParameterNames[betaFunctions_List] :=
    Module[{names = {}, rules = {}},
           names = (GetName[#])& /@ betaFunctions;
           names = Flatten[names /. a_[i1,i2] :> a];
           rules = (Rule[#, ToValidCSymbol[#]])& /@ names;
           rules = Cases[rules, HoldPattern[Rule[a_,b_]] /; a=!=b, 1];
           Return[rules];
           ];

(* create setters *)
CreateSetters[betaFunction_BetaFunction] :=
    Module[{setter = "", name = ""},
           name = ToValidCSymbolString[GetName[betaFunction]];
           setter = setter <> CConversion`CreateInlineSetter[name, GetType[betaFunction]];
           Return[setter];
          ];

CreateSetters[betaFunctions_List] :=
    Module[{setter = ""},
           (setter = setter <> CreateSetters[#])& /@ betaFunctions;
           Return[setter];
          ];

(* create getters *)
CreateGetters[betaFunction_BetaFunction] :=
    Module[{getter = "", name = ""},
           name = ToValidCSymbolString[GetName[betaFunction]];
           getter = getter <> CConversion`CreateInlineGetter[name, GetType[betaFunction]];
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
           dataType = GetCParameterType[GetType[betaFunction]];
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
           dataType = GetCParameterType[GetType[betaFunction]];
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
           dataType = GetCParameterType[GetType[betaFunction]];
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
           dataType = GetCParameterType[GetType[betaFunction]];
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
           dataType = GetCParameterType[GetType[betaFunction]];
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

End[];

EndPackage[];


BeginPackage["EWSB`", {"SARAH`", "TextFormatting`", "CConversion`", "Parameters`"}];

CheckEWSBEquations::usage="";

CreateEWSBEqPrototype::usage="creates C function prototype for a
given EWSB equation";

CreateEWSBEqFunction::usage="creates C function definition for a
given EWSB equation";

FillArrayWithEWSBEqs::usage="fills array of doubles with the values
of the EWSB equations";

FillInitialGuessArray::usage="fills a C array with initial values for the
EWSB eqs. solver";

SolveTreeLevelEwsb::usage="Solves tree-level EWSB equations";

SolveTreeLevelEwsbVia::usage="Solves tree-level EWSB equations for the
given list of parameters.  Retuns an empty string if no unique
solution can be found";

Begin["`Private`"];

freePhases = {};

AppearsInEquationOnlyAs[parameter_, equation_, function_] :=
    FreeQ[equation /. function[parameter] :> Unique[ToValidCSymbolString[parameter]], parameter];

AppearsOnlySquaredInEquation[parameter_, equation_] :=
    AppearsInEquationOnlyAs[parameter, equation, Power[#,2]&];

AppearsOnlyAbsSquaredInEquation[parameter_, equation_] :=
    AppearsInEquationOnlyAs[parameter, equation, (Susyno`LieGroups`conj[#] #)&] ||
    AppearsInEquationOnlyAs[parameter, equation, (# Susyno`LieGroups`conj[#])&];

AppearsNotInEquation[parameter_, equation_] :=
    FreeQ[equation, parameter];

CheckInEquations[parameter_, statement_, equations_List] :=
    And @@ (statement[parameter,#]& /@ equations);

CheckEWSBEquations[ewsbEqs_List, outputParameters_List] :=
    Module[{i, par, uniquePar, uniqueEqs, rules, newPhases = {}},
           For[i = 1, i <= Length[outputParameters], i++,
               par = outputParameters[[i]];
               uniquePar = Unique[CConversion`ToValidCSymbol[par]];
               rules = Join[{ par -> uniquePar},
                            Rule[#,Unique[CConversion`ToValidCSymbol[#]]]& /@ Select[outputParameters, (# =!= par)&]
                           ];
               uniqueEqs = ewsbEqs /. rules;
               If[CheckInEquations[uniquePar, AppearsNotInEquation, uniqueEqs],
                  Print["Error: ", par, " does not appear in EWSB equations!"];
                  Continue[];
                 ];
               If[Parameters`IsRealParameter[par],
                  If[CheckInEquations[uniquePar, AppearsOnlySquaredInEquation, uniqueEqs] ||
                     CheckInEquations[uniquePar, AppearsOnlyAbsSquaredInEquation, uniqueEqs],
                     Print["   Note: ", par, " appears only squared in EWSB equations."];
                     AppendTo[newPhases, FlexibleSUSY`Sign[par]];
                    ];
                  ,
                  If[CheckInEquations[uniquePar, AppearsOnlyAbsSquaredInEquation, uniqueEqs],
                     Print["   Note: ", par, " is complex and appears only ",
                           "absolute squared in EWSB equations."];
                     AppendTo[newPhases, FlexibleSUSY`Phase[par]];
                     ,
                     Print["   Note: ", par, " is complex and appears in EWSB equations."];
                     AppendTo[newPhases, FlexibleSUSY`Phase[par]];
                    ];
                 ];
              ];
           If[Length[newPhases] > 0,
              Print["   Introducing new free parameters: ", newPhases];
             ];
           freePhases = newPhases;
           Return[newPhases];
          ];

CreateEWSBEqPrototype[vev_Symbol] :=
    Module[{result = ""},
           result = "double get_ewsb_eq_" <> ToValidCSymbolString[vev] <>
                    "() const;\n";
           Return[result];
          ];

CreateEWSBEqFunction[vev_Symbol, equation_] :=
    Module[{result = "", body = "", variableName = ""},
           variableName = "ewsb_eq_" <> ToValidCSymbolString[vev];
           result = "double CLASSNAME::get_" <> variableName <>
                    "() const\n{\n";
           body = Parameters`CreateLocalConstRefsForInputParameters[equation, "LOCALINPUT"] <>
                  "\n" <> "double " <> variableName <> " = " <>
                  RValueToCFormString[equation] <> ";\n";
           body = body <> "\nreturn " <> variableName <> ";\n";
           body = IndentText[WrapLines[body]];
           Return[result <> body <> "}\n\n"];
          ];

FindFreePhase[parameter_] :=
    Module[{phases},
           phases = Cases[freePhases, FlexibleSUSY`Phase[parameter] | FlexibleSUSY`Sign[parameter]];
           If[phases === {}, Null, phases[[1]]]
          ];

SetParameterWithPhase[parameter_, gslIntputVector_String, index_Integer] :=
    Module[{result, parameterStr, freePhase, gslInput},
           parameterStr = ToValidCSymbolString[parameter];
           freePhase = FindFreePhase[parameter];
           gslInput = "gsl_vector_get(" <> gslIntputVector <> ", " <> ToString[index] <> ")";
           If[freePhase =!= Null,
              result = "model->set_" <> parameterStr <> "(" <>
                       "INPUT(" <> ToValidCSymbolString[freePhase] <> ") * " <>
                       "Abs(" <> gslInput <> "));\n";
              ,
              result = "model->set_" <> parameterStr <> "(" <> gslInput <> ");\n";
             ];
           Return[result];
          ];

FillArrayWithEWSBEqs[tadpoleEquations_List, parametersFixedByEWSB_List,
                      gslIntputVector_String:"x", gslOutputVector_String:"tadpole"] :=
    Module[{i, result = "", vev, par},
           If[Length[tadpoleEquations] =!= Length[parametersFixedByEWSB],
              Print["Error: number of EWSB equations (",Length[tadpoleEquations],
                    ") is not equal to the number of fixed parameters (",
                    Length[parametersFixedByEWSB],")"];
              Return[""];
             ];
           For[i = 1, i <= Length[parametersFixedByEWSB], i++,
               par = parametersFixedByEWSB[[i]];
               result = result <> SetParameterWithPhase[par, gslIntputVector, i-1];
              ];
           result = result <> "\n";
           For[i = 1, i <= Length[tadpoleEquations], i++,
               vev = tadpoleEquations[[i,1]];
               result = result <> gslOutputVector <> "[" <> ToString[i-1] <>
                        "] = " <> "model->get_ewsb_eq_" <>
                        ToValidCSymbolString[vev] <> "();\n";
              ];
           Return[result];
          ];

FillInitialGuessArray[parametersFixedByEWSB_List, arrayName_String:"x_init"] :=
    Module[{i, result = ""},
           For[i = 1, i <= Length[parametersFixedByEWSB], i++,
               result = result <> arrayName <> "[" <> ToString[i-1] <> "] = " <>
                        ToValidCSymbolString[parametersFixedByEWSB[[i]]] <>
                        ";\n";
              ];
           Return[result];
          ];

SimplifyEwsbEqs[equations_List, parametersFixedByEWSB_List] :=
    Module[{realParameters, simplificationRules},
           realParameters = Select[parametersFixedByEWSB, Parameters`IsRealParameter[#]&];
           simplificationRules = Flatten[{Rule[SARAH`Conj[#],#], Rule[Susyno`LieGroups`conj[#],#]}& /@ realParameters];
           equations /. simplificationRules
          ];

ExpressionsDifferBySignAtMost[{expr1_, expr2_}] :=
    (Expand[expr1 - expr2] === 0) || (Expand[expr1 + expr2] === 0);

ExpressionsDifferBySignAtMost[exprs_List] :=
    Module[{i, k},
           For[i = 1, i <= Length[exprs], i++,
               For[k = i+1, k <= Length[exprs], k++,
                   If[!ExpressionsDifferBySignAtMost[{exprs[[i]], exprs[[k]]}],
                      Return[False]
                     ];
                  ];
              ];
           Return[True];
          ];

ExpressionsAreEqual[{expr1_, expr2_}] :=
    Expand[expr1 - expr2] === 0;

ExpressionsAreEqual[exprs_List] :=
    Module[{i, k},
           For[i = 1, i <= Length[exprs], i++,
               For[k = i+1, k <= Length[exprs], k++,
                   If[!ExpressionsAreEqual[{exprs[[i]], exprs[[k]]}],
                      Return[False]
                     ];
                  ];
              ];
           Return[True];
          ];

CanReduceSolution[solution_List, signs_List] :=
    Module[{allParameters, signedParameters, signedSolutions, i, par,
            signCheck = True, unsignedCheck = True,
            unsignedParameters, unsignedSolutions, reducedSolution},
           If[Length[solution] <= 1, Return[True]];
           allParameters = DeleteDuplicates[Flatten[solution /. Rule[p_,_] :> p]];
           signedParameters = signs /. FlexibleSUSY`Sign -> Identity;
           unsignedParameters = Complement[allParameters, signedParameters];
           (* Check 1: Do the solutions for the signed parameters
              differ by a sign only? *)
           For[i = 1, i <= Length[signedParameters], i++,
               par = signedParameters[[i]];
               signedSolutions = Cases[Flatten[solution], Rule[par,expr_] :> expr];
               If[!ExpressionsDifferBySignAtMost[signedSolutions],
                  signCheck = False;
                 ];
              ];
           (* Check 2: Are the solutions for the other parameters
              equal? *)
           For[i = 1, i <= Length[unsignedParameters], i++,
               par = unsignedParameters[[i]];
               unsignedSolutions = Cases[Flatten[solution], Rule[par,expr_] :> expr];
               If[!ExpressionsAreEqual[unsignedSolutions],
                  unsignedCheck = False;
                 ];
              ];
           Return[signCheck && unsignedCheck];
          ];

ReduceSolution[solution_List, signs_List] :=
    Module[{signedParameters, reducedSolution},
           signedParameters = signs /. FlexibleSUSY`Sign -> Identity;
           reducedSolution = solution[[1]] /.
              Rule[p_, expr_] /; MemberQ[signedParameters,p] :>
              Rule[p, Global`LOCALINPUT[CConversion`ToValidCSymbol[FlexibleSUSY`Sign[p]]] Sqrt[Simplify[expr^2]]];
           Return[reducedSolution];
          ];

SolveTreeLevelEwsb[equations_List, parametersFixedByEWSB_List] :=
    Module[{result = "", simplifiedEqs, parameters,
            solution, signs, reducedSolution, i, par, expr, parStr},
           simplifiedEqs = SimplifyEwsbEqs[equations, parametersFixedByEWSB];
           simplifiedEqs = (#[[2]] == 0)& /@ simplifiedEqs;
           parameters = parametersFixedByEWSB;
           solution = Solve[simplifiedEqs, parameters];
           signs = Cases[FindFreePhase /@ parametersFixedByEWSB, FlexibleSUSY`Sign[_]];
           (* Try to reduce the solution *)
           If[CanReduceSolution[solution, signs],
              reducedSolution = ReduceSolution[solution, signs];
              (* create local const refs to input parameters appearing
                 in the solution *)
              result = Parameters`CreateLocalConstRefsForInputParameters[reducedSolution, "LOCALINPUT"] <> "\n";
              For[i = 1, i <= Length[reducedSolution], i++,
                  par  = reducedSolution[[i,1]];
                  expr = reducedSolution[[i,2]];
                  parStr = "new_" <> CConversion`ToValidCSymbolString[par];
                  result = result <>
                           "const double " <> parStr <> " = " <>
                           CConversion`RValueToCFormString[expr] <> ";\n";
                 ];
              result = result <> "\n";
              For[i = 1, i <= Length[reducedSolution], i++,
                  par  = reducedSolution[[i,1]];
                  parStr = CConversion`ToValidCSymbolString[par];
                  result = result <>
                           "if (std::isfinite(new_" <> parStr <> "))\n" <>
                               IndentText[parStr <> " = new_" <> parStr <> ";"] <> "\n" <>
                           "else\n" <>
                               IndentText["error = 1;"] <> "\n";
                 ];
              ,
              result = "error = solve_ewsb_iteratively(0);\n";
             ];
           Return[result];
          ];

SolveTreeLevelEwsbVia[equations_List, parameters_List] :=
    Module[{result = "", simplifiedEqs, solution, i, par, expr, parStr},
           simplifiedEqs = (#[[2]] == 0)& /@ equations;
           solution = Solve[simplifiedEqs, parameters];
           If[solution === {} || Length[solution] > 1,
              Print["Error: can't solve the EWSB equations for the parameters ",
                    parameters, " uniquely"];
              Print["Here are the EWSB equations we have: ", simplifiedEqs];
              Print["Here is the solution we get: ", solution];
              Return[result];
             ];
           solution = solution[[1]]; (* select first solution *)
           (* create local const refs to input parameters appearing
              in the solution *)
           result = Parameters`CreateLocalConstRefsForInputParameters[solution, "LOCALINPUT"] <> "\n";
           For[i = 1, i <= Length[solution], i++,
               par  = solution[[i,1]];
               expr = solution[[i,2]];
               parStr = "new_" <> CConversion`ToValidCSymbolString[par];
               result = result <>
               "const double " <> parStr <> " = " <>
               CConversion`RValueToCFormString[expr] <> ";\n";
              ];
           result = result <> "\n";
           For[i = 1, i <= Length[solution], i++,
               par  = solution[[i,1]];
               parStr = CConversion`ToValidCSymbolString[par];
               result = result <>
               "if (std::isfinite(new_" <> parStr <> "))\n" <>
               IndentText[parStr <> " = new_" <> parStr <> ";"] <> "\n" <>
               "else\n" <>
               IndentText["error = 1;"] <> "\n";
              ];
           Return[result];
          ];

End[];

EndPackage[];

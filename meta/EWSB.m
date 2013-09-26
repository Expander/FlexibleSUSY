
BeginPackage["EWSB`", {"SARAH`", "TextFormatting`", "CConversion`", "Parameters`"}];

FindSolutionAndFreePhases::usage="Finds solution to the EWSB and free
phases / signs."

CreateEWSBEqPrototype::usage="creates C function prototype for a
given EWSB equation";

CreateEWSBEqFunction::usage="creates C function definition for a
given EWSB equation";

FillArrayWithEWSBEqs::usage="fills array of doubles with the values
of the EWSB equations";

FillInitialGuessArray::usage="fills a C array with initial values for the
EWSB eqs. solver";

CreateTreeLevelEwsbSolver::usage="Converts tree-level EWSB solutions
to C form";

SolveTreeLevelEwsbVia::usage="Solves tree-level EWSB equations for the
given list of parameters.  Retuns an empty string if no unique
solution can be found";

Begin["`Private`"];

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

FindFreePhase[parameter_, freePhases_] :=
    Module[{phases},
           phases = Cases[freePhases, FlexibleSUSY`Phase[parameter] | FlexibleSUSY`Sign[parameter]];
           If[phases === {}, Null, phases[[1]]]
          ];

SetParameterWithPhase[parameter_, gslIntputVector_String, index_Integer, freePhases_List] :=
    Module[{result, parameterStr, freePhase, gslInput},
           parameterStr = ToValidCSymbolString[parameter];
           freePhase = FindFreePhase[parameter, freePhases];
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

FillArrayWithEWSBEqs[vevs_List, parametersFixedByEWSB_List, freePhases_List,
                     gslIntputVector_String:"x", gslOutputVector_String:"tadpole"] :=
    Module[{i, result = "", vev, par},
           If[Length[vevs] =!= Length[parametersFixedByEWSB],
              Print["Error: number of EWSB equations (",Length[vevs],
                    ") is not equal to the number of fixed parameters (",
                    Length[parametersFixedByEWSB],")"];
              Return[""];
             ];
           For[i = 1, i <= Length[parametersFixedByEWSB], i++,
               par = parametersFixedByEWSB[[i]];
               result = result <> SetParameterWithPhase[par, gslIntputVector, i-1, freePhases];
              ];
           result = result <> "\n";
           For[i = 1, i <= Length[vevs], i++,
               vev = vevs[[i]];
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

FindIndependentSubset[equations_List, parameters_List] :=
    Module[{equationSubsets, numberOfEquations, parameterSubsets,
            numberOfParameters, e, p, result = {}, isFreeOf},
           numberOfEquations = Length[equations];
           numberOfParameters = Length[parameters];
           equationSubsets = Subsets[equations, {1, numberOfEquations - 1}];
           parameterSubsets = Subsets[parameters, {1, numberOfParameters - 1}];
           For[e = 1, e <= Length[equationSubsets], e++,
               For[p = 1, p <= Length[parameterSubsets], p++,
                   isFreeOf =
                   And @@ (FreeQ[equationSubsets[[e]], #] & /@
                           parameterSubsets[[p]]);
                   If[isFreeOf,
                      AppendTo[
                          result, {equationSubsets[[e]],
                                   Complement[parameters, parameterSubsets[[p]]]}]];
                  ];
              ];

           result = Select[result, (Length[#[[1]]] == Length[#[[2]]]) &];
           Return[result];
          ];

FindMinimumByteCount[{}] := Null;

FindMinimumByteCount[lst_List] :=
    First[Sort[lst, ByteCount[#1] < ByteCount[#2]]];

EliminateOneParameter[{}, {}] := {};

EliminateOneParameter[{eq_}, {p_}] :=
    TimeConstrained[{Solve[eq, p]}, FlexibleSUSY`FSSolveEWSBTimeConstraint, {}];

EliminateOneParameter[{eq1_, eq2_}, {p1_, p2_}] :=
    Module[{reduction = {{}, {}}, rest = {}, solution},
           If[FreeQ[{eq1, eq2}, p1] || FreeQ[{eq1, eq2}, p2],
              Return[{}];
             ];
           reduction[[1]] =
           TimeConstrained[Solve[Eliminate[{eq1, eq2}, p1], p2],
                           FlexibleSUSY`FSSolveEWSBTimeConstraint, {}];
           reduction[[2]] =
           TimeConstrained[Solve[Eliminate[{eq1, eq2}, p2], p1],
                           FlexibleSUSY`FSSolveEWSBTimeConstraint, {}];
           If[reduction[[1]] === {} || reduction[[2]] === {} ||

              reduction[[1]] === {{}} || reduction[[2]] === {{}},
              Return[{}];
             ];
           If[ByteCount[reduction[[1]]] <= ByteCount[reduction[[2]]],
              solution = TimeConstrained[Solve[eq1, p1],
                                         FlexibleSUSY`FSSolveEWSBTimeConstraint, {}];
              If[solution =!= {} && solution =!= {{}}, AppendTo[rest, solution]];
              solution = TimeConstrained[Solve[eq2, p1],
                                         FlexibleSUSY`FSSolveEWSBTimeConstraint, {}];
              If[solution =!= {} && solution =!= {{}}, AppendTo[rest, solution]];
              If[rest === {}, Return[{}];];
              rest = FindMinimumByteCount[rest];
              Return[{reduction[[1]], rest}];
              ,
              solution = TimeConstrained[Solve[eq1, p2],
                                         FlexibleSUSY`FSSolveEWSBTimeConstraint, {}];
              If[solution =!= {{}}, AppendTo[rest, solution]];
              solution = TimeConstrained[Solve[eq2, p2],
                                         FlexibleSUSY`FSSolveEWSBTimeConstraint, {}];
              If[solution =!= {{}}, AppendTo[rest, solution]];
              If[rest === {}, Return[{}];];
              rest = FindMinimumByteCount[rest];
              Return[{reduction[[2]], rest}];
             ];
          ];

EliminateOneParameter[equations_List, parameters_List] :=
    Module[{independentSubset, reducedEqs, reducedPars, reducedSolution,
            complementEq, complementPar, complementSolution,
            largestIndependentSubset, s},
           independentSubset = FindIndependentSubset[equations, parameters];
           If[independentSubset === {},
              Print["EWSB equations are not reducible"];
              Return[{}];
             ];
           If[Length[independentSubset] == 1,
              largestIndependentSubset = independentSubset[[1]];,
              (* search for the largest independent equation subset *)
              largestIndependentSubset = independentSubset[[1]];
              For[s = 1, s <= Length[independentSubset], s++,
                  If[Length[independentSubset[[s,2]]] > Length[largestIndependentSubset[[2]]],
                     largestIndependentSubset = independentSubset[[s]];
                    ];
                 ];
             ];
           reducedEqs = largestIndependentSubset[[1]];
           reducedPars = largestIndependentSubset[[2]];
           reducedSolution = EliminateOneParameter[reducedEqs, reducedPars];
           If[reducedSolution === {},
              Print["Warning: could not solve reduced EWSB eqs. subset"];
              Return[{}];
             ];
           complementEq = Complement[equations, reducedEqs];
           complementPar = Complement[parameters, reducedPars];
           complementSolution = EliminateOneParameter[complementEq, complementPar];
           Join[reducedSolution, complementSolution]
          ];

MakeParameterUnique[SARAH`L[par_]] := Rule[SARAH`L[par], CConversion`ToValidCSymbol[SARAH`L[par]]];
MakeParameterUnique[SARAH`B[par_]] := Rule[SARAH`B[par], CConversion`ToValidCSymbol[SARAH`B[par]]];
MakeParameterUnique[SARAH`T[par_]] := Rule[SARAH`T[par], CConversion`ToValidCSymbol[SARAH`T[par]]];
MakeParameterUnique[par_]          :=
    { MakeParameterUnique[SARAH`L[par]],
      MakeParameterUnique[SARAH`B[par]],
      MakeParameterUnique[SARAH`T[par]] };

MakeParametersUnique[parameters_List] :=
    Flatten[MakeParameterUnique /@ parameters];

FindSolution[equations_List, parametersFixedByEWSB_List] :=
    Module[{simplifiedEqs, uniqueParameters, solution},
           simplifiedEqs = SimplifyEwsbEqs[equations, parametersFixedByEWSB];
           simplifiedEqs = (# == 0)& /@ simplifiedEqs;
           (* replace non-symbol parameters by unique symbols *)
           uniqueParameters = MakeParametersUnique[parametersFixedByEWSB];
           solution = EliminateOneParameter[
                          simplifiedEqs /. uniqueParameters,
                          parametersFixedByEWSB /. uniqueParameters];
           (* substitute back unique parameters *)
           uniqueParameters = Reverse /@ uniqueParameters;
           solution /. uniqueParameters
          ];

StripSign[Times[int_?NumericQ,expr_]] := Abs[int] expr;

StripSign[expr_] := expr;

ReduceSolution[{}] := {{},{}};

ReduceSolution[{{}}] := {{},{}};

ReduceSolution[solution_List] :=
    Module[{reducedSolution = {}, s, flattenedSolution, freePhases = {}},
           For[s = 1, s <= Length[solution], s++,
               flattenedSolution = Flatten[solution[[s]]];
               Switch[Length[flattenedSolution],
                      0,
                      Print["Warning: no solution found for the EWSB eqs."];,
                      1,
                      AppendTo[reducedSolution, flattenedSolution];,
                      2,
                      If[PossibleZeroQ[
                          flattenedSolution[[1, 2]] + flattenedSolution[[2, 2]]],
                         flattenedSolution[[1]] = flattenedSolution[[1]] /.
                         Rule[p_, expr_] :>
                         Rule[p, FlexibleSUSY`Sign[flattenedSolution[[1,1]]] StripSign[expr]];
                         AppendTo[reducedSolution, {flattenedSolution[[1]]}];
                         AppendTo[freePhases, FlexibleSUSY`Sign[flattenedSolution[[1,1]]]];
                         ,
                         Print["Warning: cannot reduce solution for ", flattenedSolution[[1,1]]];
                         Print["   because the two solutions are not related by a global sign."];
                        ];,
                      _,
                      Print["Warning: cannot reduce solution for ", flattenedSolution];
                      Print["   because there are more than two solutions"];
                     ];
              ];
           If[Length[reducedSolution] != Length[solution],
              Print["Warning: analytic reduction of EWSB solutions failed."];
              Return[{{},{}}];
             ];
           Return[{reducedSolution, freePhases}];
          ];

FindSolutionAndFreePhases[equations_List, parametersFixedByEWSB_List, outputFile_String:""] :=
    Module[{solution, reducedSolution, freePhases},
           solution = FindSolution[equations, parametersFixedByEWSB];
           {reducedSolution, freePhases} = ReduceSolution[solution];
           If[reducedSolution === {} && outputFile != "",
              Put[solution, outputFile];
             ];
           Return[{Flatten[reducedSolution], freePhases}];
          ];

CreateTreeLevelEwsbSolver[solution_List] :=
    Module[{result = "", body = "",
            i, par, expr, parStr, reducedSolution},
           reducedSolution = solution;
           If[reducedSolution =!= {},
              (* create local const refs to input parameters appearing
                 in the solution *)
              reducedSolution = reducedSolution /. {
                  FlexibleSUSY`Sign[p_]  :> Global`LOCALINPUT[CConversion`ToValidCSymbol[FlexibleSUSY`Sign[p]]],
                  FlexibleSUSY`Phase[p_] :> Global`LOCALINPUT[CConversion`ToValidCSymbol[FlexibleSUSY`Phase[p]]]
                                                   };
              result = Parameters`CreateLocalConstRefsForInputParameters[reducedSolution, "LOCALINPUT"] <> "\n";
              (* save old parameters *)
              For[i = 1, i <= Length[reducedSolution], i++,
                  par  = reducedSolution[[i,1]];
                  expr = reducedSolution[[i,2]];
                  parStr = CConversion`ToValidCSymbolString[par];
                  result = result <>
                           "const double old_" <> parStr <> " = " <> parStr <> ";\n";
                 ];
              result = result <> "\n";
              (* write solution *)
              For[i = 1, i <= Length[reducedSolution], i++,
                  par  = reducedSolution[[i,1]];
                  expr = reducedSolution[[i,2]];
                  parStr = CConversion`ToValidCSymbolString[par];
                  result = result <> parStr <> " = " <>
                           CConversion`RValueToCFormString[expr] <> ";\n";
                 ];
              result = result <> "\n";
              (* check for errors *)
              result = result <> "const bool is_finite = ";
              For[i = 1, i <= Length[reducedSolution], i++,
                  par    = reducedSolution[[i,1]];
                  parStr = CConversion`ToValidCSymbolString[par];
                  result = result <> "std::isfinite(" <> parStr <> ")";
                  If[i != Length[reducedSolution],
                     result = result <> " && ";
                    ];
                 ];
              result = result <> ";\n\n";
              body = "";
              For[i = 1, i <= Length[reducedSolution], i++,
                  par    = reducedSolution[[i,1]];
                  parStr = CConversion`ToValidCSymbolString[par];
                  body = body <> parStr <> " = old_" <> parStr <> ";\n";
                 ];
              body = body <> "error = 1;\n";
              result = result <>
                       "if (!is_finite) {\n" <>
                       IndentText[body] <>
                       "}";
              ,
              result = "error = solve_ewsb_iteratively(0);\n";
             ];
           Return[result];
          ];

SolveTreeLevelEwsbVia[equations_List, parameters_List] :=
    Module[{result = "", simplifiedEqs, solution, i, par, expr, parStr},
           simplifiedEqs = (# == 0)& /@ equations;
           solution = TimeConstrained[Solve[simplifiedEqs, parameters],
                                      FlexibleSUSY`FSSolveEWSBTimeConstraint, {}];
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


BeginPackage["EWSB`", {"SARAH`", "TextFormatting`", "CConversion`", "Parameters`", "TreeMasses`", "WriteOut`", "Utils`"}];

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

CreateEWSBRootFinders::usage="Creates comma separated list of GSL root
finders";

SetEWSBSolution::usage="sets the model parameters to the solution
provided by the solver";

FillArrayWithParameters::usage="fill an array with parameters";

DivideTadpolesByVEV::usage="Divides an array of tadpoles by their
corresponding VEV";

CreateEwsbSolverWithTadpoles::usage="Solve EWSB eqs. including
tadpoles (one step, no iteration)";

GetEWSBParametersFromGSLVector::usage="Create local copies of EWSB
output parameters from GSL vector";

SetEWSBParametersFromLocalCopies::usage="Set model parameters from
local copies";

SetEWSBParametersFromGSLVector::usage="Set model parameters from GSL
vector.";

CreateEWSBParametersInitializationList::usage="Creates initialization
list with EWSB output parameters";

Begin["`Private`"];

InputFormOfNonStrings[a_String] := a;
InputFormOfNonStrings[a_] := InputForm[a];

DebugPrint[msg___] :=
    If[FlexibleSUSY`FSDebugOutput,
       Print["Debug<EWSB>: ", Sequence @@ InputFormOfNonStrings /@ {msg}]];

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

CreateEWSBEqPrototype[higgs_] :=
    Module[{result = "", i, ctype},
           ctype = CConversion`CreateCType[CConversion`ScalarType[CConversion`realScalarCType]];
           For[i = 1, i <= TreeMasses`GetDimension[higgs], i++,
               result = result <> ctype <> " get_ewsb_eq_" <>
                        ToValidCSymbolString[higgs] <>
                        "_" <> ToString[i] <> "() const;\n";
              ];
           Return[result];
          ];

(* Creates a get_ewsb_eq_<higgs>() function for each Higgs multiplet
   entry.

   This function assumes, that the ordering of the EWSB
   eqs. corresponds to the Higgs multiplet entries.  If there are less
   equations than Higgs multiplet entries, the remaining equations
   will be set to 0.  If there are more EWSB equations than Higgs
   bosons, the remaining equations will be ignored.
 *)
CreateEWSBEqFunction[higgs_, equation_List] :=
    Module[{result = "", body = "", dim, i, eq, dimEq, type, ctype},
           dim = TreeMasses`GetDimension[higgs];
           dimEq = Length[equation];
           type = CConversion`ScalarType[CConversion`realScalarCType];
           ctype = CConversion`CreateCType[type];
           If[dim < dimEq,
              Print["Warning: number of Higgs bosons (", dim,
                    ") < number of EWSB eqs. (",dimEq,")."
                    "The EWSB eqs. ", dim, "...", dimEq,
                    " will be ignored"];
             ];
           If[dim > dimEq,
              Print["Warning: number of physical Higgs bosons (", dim,
                    ") > number of EWSB eqs. (",dimEq,").",
                    "The EWSB eqs. for the fields ", higgs, "(n), n >= ",
                    dimEq, ", will be set to zero."];
             ];
           For[i = 1, i <= dim, i++,
               result = result <>
                        ctype <> " CLASSNAME::get_ewsb_eq_" <>
                        ToValidCSymbolString[higgs] <>
                        "_" <> ToString[i] <> "() const\n{\n";
               If[i <= dimEq,
                  eq = equation[[i]];,
                  eq = 0;
                 ];
               body = Parameters`CreateLocalConstRefsForInputParameters[eq, "LOCALINPUT"] <>
                      "\n" <> ctype <> " result = " <>
                      CastTo[RValueToCFormString[eq], type] <> ";\n";
               body = body <> "\nreturn result;\n";
               body = IndentText[WrapLines[body]];
               result = result <> body <> "}\n\n"
              ];
           Return[result];
          ];

FindFreePhase[parameter_, freePhases_] :=
    Module[{phases},
           phases = Cases[freePhases, FlexibleSUSY`Phase[parameter] | FlexibleSUSY`Sign[parameter]];
           If[phases === {}, Null, phases[[1]]]
          ];

GetValueWithPhase[parameter_, gslIntputVector_String, index_Integer, freePhases_List] :=
    Module[{result, parameterStr, freePhase, gslInput},
           parameterStr = ToValidCSymbolString[parameter];
           freePhase = FindFreePhase[parameter, freePhases];
           gslInput = "gsl_vector_get(" <> gslIntputVector <> ", " <> ToString[index] <> ")";
           If[freePhase =!= Null,
              result = "INPUT(" <> ToValidCSymbolString[freePhase] <> ") * " <> "Abs(" <> gslInput <> ")";
              ,
              result = gslInput;
             ];
           Return[result];
          ];

SetParameterWithPhase[parameter_, gslIntputVector_String, index_Integer, freePhases_List] :=
    Module[{value},
           value = GetValueWithPhase[parameter, gslIntputVector, index, freePhases];
           Parameters`SetParameter[parameter, value, "model"]
          ];

FillArrayWithEWSBEqs[higgs_, gslOutputVector_String] :=
    Module[{i, result = "", par, dim},
           dim = TreeMasses`GetDimension[higgs];
           For[i = 1, i <= dim, i++,
               result = result <> gslOutputVector <> "[" <> ToString[i-1] <>
                        "] = " <> "get_ewsb_eq_" <>
                        ToValidCSymbolString[higgs] <> "_" <>
                        ToString[i] <> "();\n";
              ];
           Return[result];
          ];

InitialGuessFor[par_] :=
    If[Parameters`IsRealParameter[par], par, Abs[par]];

FillInitialGuessArray[parametersFixedByEWSB_List, arrayName_String:"x_init"] :=
    Module[{i, result = ""},
           For[i = 1, i <= Length[parametersFixedByEWSB], i++,
               result = result <> arrayName <> "[" <> ToString[i-1] <> "] = " <>
                        CConversion`RValueToCFormString[InitialGuessFor[parametersFixedByEWSB[[i]]]] <>
                        ";\n";
              ];
           Return[result];
          ];

MakeParameterUnique[Re[par_]]      := Rule[Re[par]     , CConversion`ToValidCSymbol[Re[par]]];
MakeParameterUnique[Im[par_]]      := Rule[Im[par]     , CConversion`ToValidCSymbol[Im[par]]];
MakeParameterUnique[SARAH`L[par_]] := Rule[SARAH`L[par], CConversion`ToValidCSymbol[SARAH`L[par]]];
MakeParameterUnique[SARAH`B[par_]] := Rule[SARAH`B[par], CConversion`ToValidCSymbol[SARAH`B[par]]];
MakeParameterUnique[SARAH`T[par_]] := Rule[SARAH`T[par], CConversion`ToValidCSymbol[SARAH`T[par]]];
MakeParameterUnique[SARAH`Q[par_]] := Rule[SARAH`Q[par], CConversion`ToValidCSymbol[SARAH`Q[par]]];
MakeParameterUnique[par_]          :=
    { MakeParameterUnique[SARAH`L[par]],
      MakeParameterUnique[SARAH`B[par]],
      MakeParameterUnique[SARAH`T[par]],
      MakeParameterUnique[SARAH`Q[par]] };

MakeParametersUnique[parameters_List] :=
    Flatten[MakeParameterUnique /@ parameters];

ComplexParameterReplacementRules[eqs_List, pars_List] :=
    Join[ComplexParameterReplacementRules[eqs,#]& /@ pars];

(* returns replacement rules which, if appied to eqs, lead to
   equations that are free of FlexibleSUSY`Phase[par] *)
ComplexParameterReplacementRules[eqs_List, par_] :=
    Module[{rules, replacedEqs},
           rules = {Rule[par,Abs[par] FlexibleSUSY`Phase[par]],
                    Rule[SARAH`Conj[par],Abs[par]/FlexibleSUSY`Phase[par]],
                    Rule[Susyno`LieGroups`conj[par],Abs[par]/FlexibleSUSY`Phase[par]]};
           replacedEqs = Simplify[eqs /. rules];
           If[FreeQ[replacedEqs, FlexibleSUSY`Phase[par]],
              rules,
              {}
             ]
          ];

SplitRealAndImagParts[eqs_List, pars_List] :=
    eqs /. DeleteDuplicates[Cases[pars, Re[p_] | Im[p_] :> Rule[p,Re[p]+I Im[p]]]];

SimplifyEwsbEqs[equations_List, parametersFixedByEWSB_List] :=
    Module[{realParameters, complexParameters, simplificationRules,
            renamedEqs, splitEqs},
           DebugPrint["Splitting Re[] and Im[] within EWSB eqs. ..."];
           splitEqs = SplitRealAndImagParts[equations, parametersFixedByEWSB];
           realParameters = Select[parametersFixedByEWSB, Parameters`IsRealParameter[#]&];
           DebugPrint["real parameters: ", realParameters];
           complexParameters = Complement[parametersFixedByEWSB, realParameters];
           DebugPrint["complex parameters: ", complexParameters];
           (* make parameters unique *)
           uniqueParameters = MakeParametersUnique[parametersFixedByEWSB];
           DebugPrint["Making parameters unique via: ", uniqueParameters];
           realParameters = realParameters /. uniqueParameters;
           complexParameters = complexParameters /. uniqueParameters;
           renamedEqs = splitEqs /. uniqueParameters;
           DebugPrint["EWSB eqs. with unique parameters: ", renamedEqs];
           simplificationRules =
               Flatten[Join[{Rule[SARAH`Conj[#],#],
                             Rule[Susyno`LieGroups`conj[#],#]}& /@ realParameters,
                            ComplexParameterReplacementRules[renamedEqs, complexParameters]
                           ]
                      ];
           DebugPrint["Simplification rules: ", simplificationRules];
           (* substitute back *)
           uniqueParameters = Reverse /@ uniqueParameters;
           renamedEqs /. simplificationRules /. uniqueParameters
          ];

FindIndependentSubset[equations_List, {}] := {};

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

FindMinimumByteCount[{}] := {};

FindMinimumByteCount[lst_List] :=
    First[Sort[lst, ByteCount[#1] < ByteCount[#2]]];

IsNoSolution[expr_] :=
    Head[expr] === Solve || Flatten[expr] === {};

TimeConstrainedSolve[eq_, par_] :=
    TimeConstrained[Solve[eq, par], FlexibleSUSY`FSSolveEWSBTimeConstraint, {}];

EliminateOneParameter[{}, {}] := {};

EliminateOneParameter[{eq_}, {p_}] :=
    Block[{},
          DebugPrint["eliminating ", p, ": ", eq];
          {TimeConstrainedSolve[eq, p]}
         ];

(* solves the two equations for `par' and returns the simpler solution *)
SolveRest[eq1_, eq2_, par_] :=
    Module[{rest = {}, solution},
           DebugPrint["solving rest for ", par, ": ", eq1];
           solution = TimeConstrainedSolve[eq1, par];
           If[IsNoSolution[solution],
              DebugPrint["Failed"];,
              DebugPrint["Solution found: ", solution];
              AppendTo[rest, solution];
             ];
           DebugPrint["solving rest for ", par, ": ", eq2];
           solution = TimeConstrainedSolve[eq2, par];
           If[IsNoSolution[solution],
              DebugPrint["Failed"];,
              DebugPrint["Solution found: ", solution];
              AppendTo[rest, solution];
             ];
           FindMinimumByteCount[rest]
          ];

EliminateOneParameter[{eq1_, eq2_}, {p1_, p2_}] :=
    Module[{reduction = {{}, {}}, solution, rest},
           DebugPrint["Tying to elimiate one of the parameters ",
                      {p1,p2}, " from the eqs.: ",
                      {eq1,eq2}];
           If[FreeQ[{eq1, eq2}, p1],
              Print["Error: EWSB output parameter ", p1, " does not appear in the EWSB eqs."];
              Return[{}];
             ];
           If[FreeQ[{eq1, eq2}, p2],
              Print["Error: EWSB output parameter ", p2, " does not appear in the EWSB eqs."];
              Return[{}];
             ];
           (* special case: no elimination needed *)
           If[!FreeQ[{eq1},p1] && FreeQ[{eq1},p2] &&
              !FreeQ[{eq2},p2] && FreeQ[{eq2},p1],
              DebugPrint["The two equations are independent of each other."];
              DebugPrint["Step 1: solving for ", p1, ": ", eq1];
              reduction[[1]] = TimeConstrainedSolve[{eq1}, p1];
              If[IsNoSolution[reduction[[1]]],
                 DebugPrint["Failed"];
                 Return[{}];,
                 DebugPrint["Solution: ", reduction[[1]]];
                ];
              DebugPrint["Step 2: solving for ", p2, ": ", eq2];
              reduction[[2]] = TimeConstrainedSolve[{eq2}, p2];
              If[IsNoSolution[reduction[[2]]],
                 DebugPrint["Failed"];
                 Return[{}];,
                 DebugPrint["Solution: ", reduction[[2]]];
                ];
              DebugPrint["Full solution: ", reduction];
              Return[reduction];
             ];
           If[!FreeQ[{eq1},p2] && FreeQ[{eq1},p1] &&
              !FreeQ[{eq2},p1] && FreeQ[{eq2},p2],
              DebugPrint["The two equations are independent of each other."];
              DebugPrint["Step 1: solving for ", p2, ": ", eq1];
              reduction[[1]] = TimeConstrainedSolve[{eq1}, p2];
              If[IsNoSolution[reduction[[1]]],
                 DebugPrint["Failed"];
                 Return[{}];,
                 DebugPrint["Solution: ", reduction[[1]]];
                ];
              DebugPrint["Step 2: solving for ", p1, ": ", eq2];
              reduction[[2]] = TimeConstrainedSolve[{eq2}, p1];
              If[IsNoSolution[reduction[[2]]],
                 DebugPrint["Failed"];
                 Return[{}];,
                 DebugPrint["Solution: ", reduction[[2]]];
                ];
              DebugPrint["Full solution: ", reduction];
              Return[reduction];
             ];
           DebugPrint["eliminate ", p1, " and solve for ", p2, ": ", {eq1, eq2}];
           reduction[[1]] =
           TimeConstrainedSolve[Eliminate[{eq1, eq2}, p1], p2];
           DebugPrint["eliminate ", p2, " and solve for ", p1, ": ", {eq1, eq2}];
           reduction[[2]] =
           TimeConstrainedSolve[Eliminate[{eq1, eq2}, p2], p1];
           If[IsNoSolution[reduction[[1]]] || IsNoSolution[reduction[[2]]],
              DebugPrint["Failed"];
              Return[{}];
             ];
           If[ByteCount[reduction[[1]]] <= ByteCount[reduction[[2]]],
              DebugPrint["continue with solution 1, because it is simpler:",
                         reduction[[1]]];
              rest = SolveRest[eq1, eq2, p1];
              If[IsNoSolution[rest],
                 DebugPrint["could not solve rest for solution 1"];
                 Return[{}];
                ];
              Return[{reduction[[1]], rest}];
              ,
              DebugPrint["continue with solution 2, because it is simpler:",
                         reduction[[2]]];
              rest = SolveRest[eq1, eq2, p2];
              If[IsNoSolution[rest],
                 DebugPrint["could not solve rest for solution 2"];
                 Return[{}];
                ];
              Return[{reduction[[2]], rest}];
             ];
          ];

EliminateOneParameter[equations_List, parameters_List] :=
    Module[{independentSubset, reducedEqs, reducedPars, reducedSolution,
            complementEq, complementPar, complementSolution,
            largestIndependentSubset, s},
           independentSubset = FindIndependentSubset[equations, parameters];
           If[independentSubset === {},
              DebugPrint["EWSB equations are not reducible"];
              Return[{}];
             ];
           (* search for the largest independent equation subset *)
           If[Length[independentSubset] == 1,
              largestIndependentSubset = independentSubset[[1]];,
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
              DebugPrint["Could not solve reduced EWSB eqs. subset"];
              Return[{}];
             ];
           complementEq = Complement[equations, reducedEqs];
           complementPar = Complement[parameters, reducedPars];
           complementSolution = EliminateOneParameter[complementEq, complementPar];
           If[complementSolution === {},
              DebugPrint["Could not solve remaining EWSB eqs. subset"];
              Return[{}];
             ];
           DebugPrint["Solution = ", Join[reducedSolution, complementSolution]];
           Join[reducedSolution, complementSolution]
          ];

ToMathematicaSolutionFormat[{}] := {};

ToMathematicaSolutionFormat[sol_List] :=
    Tuples[Flatten /@ sol];

FindSolution[equations_List, parametersFixedByEWSB_List] :=
    Module[{simplifiedEqs, makeParsUnique, solution,
            uniquePars, uniqueEqs},
           DebugPrint["Simplifying the EWSB eqs. ..."];
           simplifiedEqs = SimplifyEwsbEqs[equations, parametersFixedByEWSB];
           simplifiedEqs = (# == 0)& /@ simplifiedEqs;
           DebugPrint["Simplified EWSB eqs.: ", simplifiedEqs];
           (* replace non-symbol parameters by unique symbols *)
           makeParsUnique = MakeParametersUnique[parametersFixedByEWSB];
           uniquePars = parametersFixedByEWSB /. makeParsUnique;
           uniqueEqs = simplifiedEqs /. makeParsUnique;
           DebugPrint["Eliminating the parameters ", uniquePars];
           solution = ToMathematicaSolutionFormat @ EliminateOneParameter[uniqueEqs, uniquePars];
           If[solution === {},
              solution = TimeConstrainedSolve[uniqueEqs, uniquePars];
             ];
           (* substitute back unique parameters *)
           makeParsUnique = Reverse /@ makeParsUnique;
           solution /. makeParsUnique
          ];

StripSign[Times[int_?NumericQ,expr_]] := Abs[int] expr;

StripSign[expr_] := expr;

ReduceSolution[{}] := {{},{}};

ReduceSolution[{{}}] := {{},{}};

SignOrPhase[par_] :=
    If[Parameters`IsRealParameter[par],
       FlexibleSUSY`Sign[par],
       FlexibleSUSY`Phase[par]];

ReduceTwoSolutions[sol1_, sol2_] :=
    Module[{par, signOrPhase, reducedSolution},
           par = sol1[[1]];
           DebugPrint["Reducing solutions for ", par, "..."];
           signOrPhase = SignOrPhase[par];
           If[PossibleZeroQ[sol1[[2]] - sol2[[2]]],
              DebugPrint["The two solutions for ", par, " are identical"];
              Return[{sol1}];
             ];
           If[!PossibleZeroQ[sol1[[2]] + sol2[[2]]],
              Print["Warning: cannot reduce solution for ", par];
              Print["   because the two solutions are not related by a global sign."];
              Return[{}];
             ];
           DebugPrint["The two solutions for ", par,
                      " are related by a global sign/phase: ",
                      sol1, ", ", sol2];
           reducedSolution = sol1 /.
               Rule[p_, expr_] :>
               Rule[p, signOrPhase StripSign[expr]];
           DebugPrint["=> the reduced solution is: ",
                      {reducedSolution, signOrPhase}];
           {reducedSolution, signOrPhase}
          ];

ReduceSolution[{sol_}] := {{sol},{}};

ReduceSolution[{sol1_, sol2_}] :=
    Module[{reducedSolution = {}, freePhases = {}, s, red},
           DebugPrint["Reducing the two solutions: ", {sol1,sol2}];
           For[s = 1, s <= Length[sol1], s++,
               red = ReduceTwoSolutions[sol1[[s]], sol2[[s]]];
               Switch[Length[red],
                      1, AppendTo[reducedSolution, red[[1]]];,
                      2, AppendTo[reducedSolution, red[[1]]];
                         AppendTo[freePhases, red[[2]]];
                     ];
              ];
           If[Length[reducedSolution] != Length[sol1],
              Print["Warning: analytic reduction of EWSB solutions failed."];
              Return[{{},{}}];
             ];
           Return[{reducedSolution, freePhases}];
          ];

ReduceSolution[solution_List] :=
    Module[{},
           Print["Error: cannot reduce the solution ", solution];
           Print["   because there are more than two solutions"];
           {{},{}}
          ];

FindSolutionAndFreePhases[equations_List, parametersFixedByEWSB_List, outputFile_String:""] :=
    Module[{solution, reducedSolution, freePhases},
           solution = FindSolution[equations, parametersFixedByEWSB];
           {reducedSolution, freePhases} = ReduceSolution[solution];
           DebugPrint["The full, reduced solution to the EWSB eqs. is:",
                      reducedSolution];
           If[reducedSolution === {} && outputFile != "",
              Put[solution, outputFile];
             ];
           Return[{Flatten[reducedSolution], freePhases}];
          ];

CreateTreeLevelEwsbSolver[solution_List] :=
    Module[{result = "", body = "",
            i, par, expr, parStr, oldParStr, reducedSolution,
            type},
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
                  type = CConversion`CreateCType[Parameters`GetType[par]];
                  parStr = CConversion`RValueToCFormString[par];
                  oldParStr = "old_" <> CConversion`ToValidCSymbolString[par];
                  result = result <>
                           "const " <> type <>" " <> oldParStr <> " = " <> parStr <> ";\n";
                 ];
              result = result <> "\n";
              (* write solution *)
              For[i = 1, i <= Length[reducedSolution], i++,
                  par  = reducedSolution[[i,1]];
                  expr = reducedSolution[[i,2]];
                  result = result <> Parameters`SetParameter[par, expr];
                 ];
              result = result <> "\n";
              (* check for errors *)
              result = result <> "const bool is_finite = ";
              For[i = 1, i <= Length[reducedSolution], i++,
                  par    = reducedSolution[[i,1]];
                  parStr = CConversion`RValueToCFormString[par];
                  result = result <> "IsFinite(" <> parStr <> ")";
                  If[i != Length[reducedSolution],
                     result = result <> " && ";
                    ];
                 ];
              result = result <> ";\n\n";
              body = "";
              For[i = 1, i <= Length[reducedSolution], i++,
                  par    = reducedSolution[[i,1]];
                  oldParStr = "old_" <> CConversion`ToValidCSymbolString[par];
                  body = body <> Parameters`SetParameter[par, oldParStr, None];
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

SolveTreeLevelEwsbVia[equations_List, {}] :=
    Module[{},
           Print["Error: SolveTreeLevelEwsbVia: list of output parameters is empty"];
           Quit[1];
          ];

SolveTreeLevelEwsbVia[equations_List, parameters_List] :=
    Module[{result = "", simplifiedEqs, solution, i, par, expr, parStr, type, ctype},
           If[Length[equations] =!= Length[parameters],
              Print["Warning: SolveTreeLevelEwsbVia: trying to solve ",
                    Length[equations], " equations for ", Length[parameters],
                    " parameters ", parameters];
             ];
           simplifiedEqs = (# == 0)& /@ equations;
           simplifiedEqs = Parameters`FilterOutIndependentEqs[simplifiedEqs, parameters];
           solution = TimeConstrainedSolve[simplifiedEqs, parameters];
           If[solution === {} || Length[solution] > 1,
              Print["Error: can't solve the EWSB equations for the parameters ",
                    parameters, " uniquely"];
              Print["Here are the EWSB equations we have: ", InputForm[simplifiedEqs]];
              Print["Here is the solution we get: ", InputForm[solution]];
              Return[result];
             ];
           solution = solution[[1]]; (* select first solution *)
           (* create local const refs to input parameters appearing
              in the solution *)
           result = Parameters`CreateLocalConstRefsForInputParameters[solution, "LOCALINPUT"] <> "\n";
           For[i = 1, i <= Length[solution], i++,
               par  = solution[[i,1]];
               expr = solution[[i,2]];
               type = Parameters`GetType[par];
               ctype = CConversion`CreateCType[type];
               parStr = "new_" <> CConversion`ToValidCSymbolString[par];
               result = result <>
               "const " <> ctype <> " " <> parStr <> " = " <>
               CConversion`CastTo[CConversion`RValueToCFormString[expr],type] <> ";\n";
              ];
           result = result <> "\n";
           For[i = 1, i <= Length[solution], i++,
               par  = solution[[i,1]];
               parStr = CConversion`ToValidCSymbolString[par];
               result = result <>
               "if (IsFinite(new_" <> parStr <> "))\n" <>
               IndentText[CConversion`RValueToCFormString[par] <>
                          " = new_" <> parStr <> ";"] <> "\n" <>
               "else\n" <>
               IndentText["error = 1;"] <> "\n";
               If[i < Length[solution],
                  result = result <> "\n";
                 ];
              ];
           Return[result];
          ];

CreateNewEWSBRootFinder[] :=
    "new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, ";

CreateEWSBRootFinder[rootFinder_ /; rootFinder === FlexibleSUSY`FPIRelative] :=
    "new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_relative>(CLASSNAME::ewsb_step, &params, number_of_ewsb_iterations, ewsb_iteration_precision)";

CreateEWSBRootFinder[rootFinder_ /; rootFinder === FlexibleSUSY`FPIAbsolute] :=
    "new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_absolute>(CLASSNAME::ewsb_step, &params, number_of_ewsb_iterations, ewsb_iteration_precision)";

CreateEWSBRootFinder[rootFinder_ /; rootFinder === FlexibleSUSY`FPITadpole] :=
    "new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_tadpole>(CLASSNAME::ewsb_step, &params, number_of_ewsb_iterations, fixed_point_iterator::Convergence_tester_tadpole(ewsb_iteration_precision, CLASSNAME::tadpole_equations, &params))";

CreateEWSBRootFinder[rootFinder_ /; rootFinder === FlexibleSUSY`GSLHybrid] :=
    CreateNewEWSBRootFinder[] <> "gsl_multiroot_fsolver_hybrid)";

CreateEWSBRootFinder[rootFinder_ /; rootFinder === FlexibleSUSY`GSLHybridS] :=
    CreateNewEWSBRootFinder[] <> "gsl_multiroot_fsolver_hybrids)";

CreateEWSBRootFinder[rootFinder_ /; rootFinder === FlexibleSUSY`GSLBroyden] :=
    CreateNewEWSBRootFinder[] <> "gsl_multiroot_fsolver_broyden)";

CreateEWSBRootFinder[rootFinder_ /; rootFinder === FlexibleSUSY`GSLNewton] :=
    CreateNewEWSBRootFinder[] <> "gsl_multiroot_fsolver_dnewton)";

CreateEWSBRootFinders[{}] :=
    Block[{},
          Print["Error: List of EWSB root finders must not be empty!"];
          Quit[1];
         ];

CreateEWSBRootFinders[rootFinders_List] :=
    Utils`StringJoinWithSeparator[CreateEWSBRootFinder /@ rootFinders, ",\n"];

WrapPhase[phase_ /; phase === Null, str_String] :=
    str;

WrapPhase[phase_, str_String] :=
    "LOCALINPUT(" <> CConversion`ToValidCSymbolString[phase] <> ")*Abs(" <> str <> ")";

SetEWSBSolution[par_, idx_, phase_, func_String] :=
    Parameters`SetParameter[par, WrapPhase[phase, func <> "(" <> ToString[idx-1] <> ")"], None];

SetEWSBSolution[parametersFixedByEWSB_List, freePhases_List, func_String] :=
    Module[{result = "", i, phase},
           For[i = 1, i <= Length[parametersFixedByEWSB], i++,
               phase = FindFreePhase[parametersFixedByEWSB[[i]], freePhases];
               result = result <> SetEWSBSolution[parametersFixedByEWSB[[i]], i, phase, func];
              ];
           result
          ];

ConvertToReal[par_] :=
    If[Parameters`IsRealParameter[par],
       CConversion`ToValidCSymbolString[par],
       "Abs(" <> CConversion`ToValidCSymbolString[par] <> ")"
      ];

FillArrayEntryWithParameter[arrayName_String, par_, idx_] :=
    arrayName <> "[" <> ToString[idx-1] <> "] = " <> ConvertToReal[par] <> ";\n";

FillArrayWithParameters[arrayName_String, parameters_List] :=
    Module[{result = "", i},
           For[i = 1, i <= Length[parameters], i++,
               result = result <> FillArrayEntryWithParameter[arrayName, parameters[[i]], i];
              ];
           result
          ];

DivideArrayEntryByParameter[arrayName_String, par_, idx_] :=
    arrayName <> "[" <> ToString[idx-1] <> "] /= " <> CConversion`ToValidCSymbolString[par] <> ";\n";

DivideTadpolesByVEV[arrayName_String, vevToTadpoleAssociation_List] :=
    Module[{result = "", i, vevs},
           vevs = #[[3]]& /@ vevToTadpoleAssociation;
           For[i = 1, i <= Length[vevs], i++,
               result = result <> DivideArrayEntryByParameter[arrayName, vevs[[i]], i];
              ];
           result
          ];

CreateEwsbSolverWithTadpoles[solution_List, softHiggsMassToTadpoleAssociation_List] :=
    Module[{result = "", i, par, expr, parStr, reducedSolution, rules, type},
           reducedSolution = solution /.
               FlexibleSUSY`tadpole[p_] :> CConversion`ReleaseHoldAt[HoldForm[FlexibleSUSY`tadpole[[p-1]]], {1,2}];
           If[reducedSolution =!= {},
              (* create local const refs to input parameters appearing
                 in the solution *)
              reducedSolution = reducedSolution /. {
                  FlexibleSUSY`Sign[p_]  :> Global`LOCALINPUT[CConversion`ToValidCSymbol[FlexibleSUSY`Sign[p]]],
                  FlexibleSUSY`Phase[p_] :> Global`LOCALINPUT[CConversion`ToValidCSymbol[FlexibleSUSY`Phase[p]]]
                                                   };
              result = Parameters`CreateLocalConstRefsForInputParameters[reducedSolution, "LOCALINPUT"];
              (* define variables for new parameters *)
              For[i = 1, i <= Length[reducedSolution], i++,
                  par  = reducedSolution[[i,1]];
                  expr = reducedSolution[[i,2]];
                  type = CConversion`CreateCType[Parameters`GetType[par]];
                  parStr = CConversion`ToValidCSymbolString[par];
                  result = result <> type <> " " <> parStr <> ";\n";
                 ];
              result = result <> "\n";
              (* write solution *)
              For[i = 1, i <= Length[reducedSolution], i++,
                  par  = reducedSolution[[i,1]];
                  expr = reducedSolution[[i,2]];
                  type = Parameters`GetType[par];
                  parStr = CConversion`ToValidCSymbolString[par];
                  result = result <> parStr <> " = " <>
                           CConversion`CastTo[CConversion`RValueToCFormString[expr],type] <> ";\n";
                 ];
              result = result <> "\n";
              (* check for errors *)
              result = result <> "const bool is_finite = ";
              For[i = 1, i <= Length[reducedSolution], i++,
                  par    = reducedSolution[[i,1]];
                  parStr = CConversion`ToValidCSymbolString[par];
                  result = result <> "IsFinite(" <> parStr <> ")";
                  If[i != Length[reducedSolution],
                     result = result <> " && ";
                    ];
                 ];
              result = result <> ";\n";
              ,
              result = "const bool is_finite = false;\n";
             ];
           Return[result];
          ];

GetEWSBParametersFromGSLVector[parametersFixedByEWSB_List, freePhases_List,
                               gslIntputVector_String:"x"] :=
    Module[{i, result = "", par, parStr, type},
           For[i = 1, i <= Length[parametersFixedByEWSB], i++,
               par = parametersFixedByEWSB[[i]];
               type = CConversion`CreateCType[Parameters`GetType[par]];
               parStr = CConversion`ToValidCSymbolString[par];
               result = result <>
                        "const " <> type <> " " <> parStr <> " = " <>
                        GetValueWithPhase[par, gslIntputVector, i-1, freePhases] <> ";\n";
              ];
           Return[result];
          ];

SetEWSBParametersFromLocalCopies[parameters_List, struct_String] :=
    Module[{result = ""},
           (result = result <> Parameters`SetParameter[#, CConversion`ToValidCSymbolString[#], struct])& /@ parameters;
           result
          ];

SetEWSBParametersFromGSLVector[parametersFixedByEWSB_List, freePhases_List,
                               gslIntputVector_String] :=
    Module[{i, result = "", par},
           For[i = 1, i <= Length[parametersFixedByEWSB], i++,
               par = parametersFixedByEWSB[[i]];
               result = result <> SetParameterWithPhase[par, gslIntputVector, i-1, freePhases];
              ];
           Return[result];
          ];

CreateEWSBParametersInitializationList[parameters_List] :=
    Module[{result = "{}"},
           If[Length[parameters] > 0,
              result = Utils`StringJoinWithSeparator[
                  ConvertToReal /@ parameters,
                  ", "
              ];
              result = "{ " <> result <> " }";
             ];
           result
          ];

End[];

EndPackage[];

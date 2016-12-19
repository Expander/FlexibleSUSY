
BeginPackage["SemiAnalytic`", {"SARAH`", "CConversion`", "Constraint`", "Parameters`",
                               "TextFormatting`", "Utils`", "WriteOut`"}];

SemiAnalyticSolution::usage="Head of a semi-analytic solution.
A semi-analytic solution to the RGEs has the structure
SemiAnalyticSolution[parameter, {basis}]";

CheckSemiAnalyticBoundaryConditions::usage="";
SelectSemiAnalyticConstraint::usage="";

SetSemiAnalyticParameters::usage="";
GetSemiAnalyticParameters::usage="";

IsAllowedSemiAnalyticParameter::usage="";
IsSemiAnalyticParameter::usage="";

GetSemiAnalyticSolutions::usage="Constructs the semi-analytic
solutions implied by the given list of boundary conditions.";

CreateSemiAnalyticSolutionsDefinitions::usage="";
CreateSemiAnalyticSolutionsInitialization::usage="";

Begin["`Private`"];

allSemiAnalyticParameters = {};

GetName[SemiAnalyticSolution[name_, basis_List]] := name;

GetBasis[SemiAnalyticSolution[name_, basis_List]] := basis;

IsDimensionOne[par_] :=
    Module[{dimOnePars},
           dimOnePars = { SARAH`BetaTijk };
           If[SARAH`SupersymmetricModel,
              dimOnePars = Append[dimOnePars, SARAH`BetaMi];,
              dimOnePars = Append[dimOnePars, SARAH`BetaMuij];
             ];
           dimOnePars = (Parameters`StripIndices[#[[1]]])& /@ (Join @@ dimOnePars);
           MemberQ[dimOnePars, Parameters`StripIndices[par]]
          ];

IsDiracGauginoMass[par_] :=
    Module[{diracMasses = {}},
           If[SARAH`SupersymmetricModel,
              diracMasses = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaDGi;
             ];
           MemberQ[diracMasses, Parameters`StripIndices[par]]
          ];

IsScalarMass[par_] :=
    Module[{scalarMasses},
           If[SARAH`SupersymmetricModel,
              scalarMasses = Parameters`StripIndices[#[[1]]]& /@ SARAH`Betam2ij;,
              scalarMasses = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaBij;
             ];
           MemberQ[scalarMasses, Parameters`StripIndices[par]]
          ];

IsSoftBilinear[par_] :=
    Module[{softBilinears = {}},
           If[SARAH`SupersymmetricModel,
              softBilinears = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaBij;
             ];
           MemberQ[softBilinears, Parameters`StripIndices[par]]
          ];

IsSoftLinear[par_] :=
    Module[{softLinears = {}},
           If[SARAH`SupersymmetricModel,
              softLinears = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaLSi;
             ];
           MemberQ[softLinears, Parameters`StripIndices[par]]
          ];

IsAllowedSemiAnalyticParameter[par_] :=
    Or[IsDimensionOne[par],
       IsDiracGauginoMass[par],
       IsScalarMass[par],
       IsSoftBilinear[par],
       IsSoftLinear[par]];

IsSemiAnalyticParameter[par_] := MemberQ[allSemiAnalyticParameters, Parameters`StripIndices[par]];

SetSemiAnalyticParameters[parameters_List] :=
    Module[{},
           allSemiAnalyticParameters = DeleteDuplicates[Select[parameters, (IsAllowedSemiAnalyticParameter[#])&]];
           allSemiAnalyticParameters = (Parameters`StripIndices[#])& /@ allSemiAnalyticParameters;
          ];

GetSemiAnalyticParameters[] := allSemiAnalyticParameters;

CheckSemiAnalyticBoundaryConditions[constraints_List] :=
    Module[{i, sortedPars, fixedPars},
           sortedPars = Sort[allSemiAnalyticParameters];
           For[i = 1, i <= Length[constraints], i++,
               fixedPars = Sort[Intersection[sortedPars,
                                             Constraint`FindFixedParametersFromConstraint[constraints[[i]]]]];
               If[fixedPars =!= {} && fixedPars =!= sortedPars,
                  Print["Error: all semi-analytic parameters must be set at the same scale."];
                  Return[False];
                 ];
              ];
           True
          ];

SelectSemiAnalyticConstraint[constraints_List] :=
    Module[{i, sortedPars, fixedPars, result = {}},
           sortedPars = Sort[allSemiAnalyticParameters];
           For[i = 1, i <= Length[constraints], i++,
               fixedPars = Sort[Intersection[sortedPars,
                                             Constraint`FindFixedParametersFromConstraint[constraints[[i]]]]];
               If[fixedPars =!= {} && fixedPars === sortedPars,
                  result = constraints[[i]];
                 ];
              ];
           result
          ];

SelectParametersWithMassDimension[parameters_List, dim_?IntegerQ] :=
    Module[{allParameters},
           allParameters = Parameters`GetModelParametersWithMassDimension[dim];
           Select[parameters, MemberQ[allParameters, #]&]
          ];

GetParameterExpansionRules[par_] :=
    Module[{dims, numDims, i, indices, rules = {}},
           dims = Parameters`GetParameterDimensions[par];
           If[dims =!= {1},
              numDims = Length[dims];
              indices = Table["i" <> ToString[i], {i, 1, numDims}];
              (* do not expand when indices are explicitly given *)
              rules = { RuleDelayed @@ Rule[par[Sequence @@ (Pattern[Evaluate[Symbol[#]], Blank[]]& /@ indices)],
                                            par[Sequence @@ (Symbol[#] & /@ indices)]],
                        RuleDelayed @@ Rule[par, Table[par[Sequence @@ (Symbol[#] & /@ indices)],
                                                       Evaluate[Sequence @@ (Inner[Function[{idx, bound}, {Symbol[idx], bound}],
                                                                                   indices, dims, List])]]]
                       };
             ];
           rules
          ];

GetMacroExpansionRules[] :=
    { RuleDelayed[CConversion`UNITMATRIX[n_Integer], IdentityMatrix[n]],
      RuleDelayed[CConversion`UNITMATRIXCOMPLEX[n_Integer], IdentityMatrix[n]],
      RuleDelayed[CConversion`ZEROARRAY[n_Integer], Table[0, {n}]],
      RuleDelayed[CConversion`ZEROARRAYCOMPLEX[n_Integer], Table[0, {n}]],
      RuleDelayed[CConversion`ZEROVECTOR[n_Integer], Table[0, {n}]],
      RuleDelayed[CConversion`ZEROVECTORCOMPLEX[n_Integer], Table[0, {n}]],
      RuleDelayed[CConversion`ZEROMATRIX[m_Integer, n_Integer], Table[0, {m}, {n}]],
      RuleDelayed[CConversion`ZEROMATRIXCOMPLEX[m_Integer, n_Integer], Table[0, {m}, {n}]],
      RuleDelayed[CConversion`ZEROTENSOR3[l_Integer, m_Integer, n_Integer], Table[0, {l}, {m}, {n}]],
      RuleDelayed[CConversion`ZEROTENSOR3COMPLEX[l_Integer, m_Integer, n_Integer], Table[0, {l}, {m}, {n}]],
      RuleDelayed[CConversion`ZEROTENSOR4[l_Integer, m_Integer, n_Integer, o_Integer], Table[0, {l}, {m}, {n}, {o}]],
      RuleDelayed[CConversion`ZEROTENSOR4COMPLEX[l_Integer, m_Integer, n_Integer, o_Integer], Table[0, {l}, {m}, {n}, {o}]],
      RuleDelayed[CConversion`ZEROTENSOR[dimensions__ /; And @@ (IntegerQ /@ {dimensions})],
                  Table[0, Evaluate[Sequence @@ Map[List, List[dimensions]]]]],
      RuleDelayed[CConversion`ZEROTENSORCOMPLEX[dimensions__ /; And @@ (IntegerQ /@ {dimensions})],
                  Table[0, Evaluate[Sequence @@ Map[List, List[dimensions]]]]]
    };

GetBoundaryValueFromSetting[par_, value_] :=
    Module[{dims, numDims, i, indices,
            parameters, parameterExpansions, lhs,
            rhs, result},
           dims = Parameters`GetParameterDimensions[par];
           parameters = Parameters`FindAllParameters[value];
           parameterExpansions = Flatten[GetParameterExpansionRules /@ parameters, 1];
           rhs = value /. Join[GetMacroExpansionRules[], parameterExpansions];
           lhs = par /. GetParameterExpansionRules[par];
           (* construct rules *)
           If[dims =!= {1},
              (* check dimensions are consistent *)
              If[Dimensions[rhs] =!= dims,
                 Print["Error: dimensions do not match in setting: " {par, value}];
                 Quit[1];
                ];
              numDims = Length[dims];
              indices = Table["i" <> ToString[i], {i, 1, numDims}];
              result = Flatten[Table[Rule[Part[lhs,Sequence @@ (Symbol[#]& /@ indices)],
                                          Part[rhs, Sequence @@ (Symbol[#]& /@ indices)]],
                                     Evaluate[Sequence @@ (Inner[Function[{idx, bound}, {Symbol[idx], bound}],
                                                                 indices, dims, List])]]];,
              (* @todo distinguish length 1 vectors from scalars? *)
              result = {Rule[lhs, rhs]};
             ];
           result
          ];

(* Returns a list of rules of the form {parameter -> value, ...}
   determined from the provided constraint settings.

   The list of settings are assumed to apply at a single scale.
   The value in each rule is the final result after applying all
   of the settings in the list, in order.
*)
GetBoundaryValueSubstitutions[settings_List] :=
    Module[{i, noTemp, test, fixedPars, boundaryValues, rules, finalRules = {}},
           noTemp = DeleteCases[settings, {FlexibleSUSY`Temporary[_], _}];
           test = Function[{first, second}, first[[1]] === second[[1]]];
           For[i = 1, i <= Length[noTemp], i++,
               Switch[noTemp[[i]],
                      {_, _},
                      rules = GetBoundaryValueFromSetting[noTemp[[i,1]], noTemp[[i,2]]];
                      finalRules = Last[(finalRules = Utils`AppendOrReplaceInList[finalRules, #, test])& /@ rules];,
                      FlexibleSUSY`FSMinimize[__],
                      fixedPars = Constraint`FindFixedParametersFromConstraint[{noTemp[[i]]}];
                      boundaryValues = Symbol[CConversion`ToValidCSymbolString[#] <> "MinSol"]& /@ fixedPars;
                      finalRules = Last[(finalRules = Utils`AppendOrReplaceInList[finalRules, Rule[#[[1]], #[[2]]], test])&
                                        /@ Utils`Zip[fixedPars, boundaryValues]];,
                      FlexibleSUSY`FSFindRoot[__],
                      fixedPars = Constraint`FindFixedParametersFromConstraint[{noTemp[[i]]}];
                      boundaryValues = Symbol[CConversion`ToValidCSymbolString[#] <> "RootSol"]& /@ fixedPars;
                      finalRules = Last[(finalRules = Utils`AppendOrReplaceInList[finalRules, Rule[#[[1]], #[[2]]], test])&
                                        /@ Utils`Zip[fixedPars, boundaryValues]];,
                      FlexibleSUSY`FSSolveEWSBFor[__],
                      fixedPars = Constraint`FindFixedParametersFromConstraint[{noTemp[[i]]}];
                      boundaryValues = Symbol[CConversion`ToValidCSymbolString[#] <> "EWSBSol"]& /@ fixedPars;
                      finalRules = Last[(finalRules = Utils`AppendOrReplaceInList[finalRules, Rule[#[[1]], #[[2]]], test])&
                                        /@ Utils`Zip[fixedPars, boundaryValues]];
                     ];
              ];
           finalRules
          ];

SumOverElements[par_] :=
    Module[{i, dims, numDims, indices, result},
           dims = Parameters`GetParameterDimensions[par];
           If[dims === {1},
              result = par;,
              numDims = Length[dims];
              indices = Table["i" <> ToString[i], {i, 1, numDims}];
              result = Sum[par[Evaluate[Sequence @@ (Symbol[#]& /@ indices)]],
                           Evaluate[Sequence @@ (Inner[Function[{idx, bound}, {Symbol[idx], bound}],
                                                       indices, dims, List])]];
             ];
           result
          ];

(* @todo handle cases like TYu = Yu *)
GetCoefficientAndTerm[monomial_NumericQ, dimZeroPars_List] := {monomial, Null};

GetCoefficientAndTerm[monomial_Times, dimZeroPars_List] :=
    Module[{i, operands, coeff = {}, term = {}},
           operands = List @@ monomial;
           For[i = 1, i <= Length[operands], i++,
               If[NumericQ[operands[[i]]] || DependsAtMostOn[operands[[i]], dimZeroPars],
                  coeff = Append[coeff, operands[[i]]];,
                  term = Append[term, operands[[i]]];
                 ];
              ];
           Which[term === {} && coeff =!= {},
                 coeff = Times @@ coeff;
                 term = Null;,
                 term =!= {} && coeff === {},
                 coeff = 1;
                 term = Times @@ term;,
                 True,
                 coeff = Times @@ coeff;
                 term = Times @@ term;
                ];
           {coeff, term}
          ];

GetCoefficientAndTerm[monomial_, dimZeroPars_List] :=
    Module[{coeff, term},
           If[DependsAtMostOn[monomial, dimZeroPars],
              coeff = monomial;
              term = Null;,
              coeff = 1;
              term = monomial;
             ];
           {coeff, term}
          ];

GetSolutionBasis[pars_List, subs_List] :=
    Module[{i, dimZeroPars, dimensions, expr,
            monomials = {}, coeff, term, basis = {}},
           If[pars =!= {},
              dimZeroPars = Parameters`GetModelParametersWithMassDimension[0];
              dimensions = Parameters`GetParameterDimensions[#]& /@ pars;
              expr = Plus @@ (SumOverElements[#]& /@ pars);
              expr = Expand[expr /. subs];
              (* get independent basis elements from result *)
              If[Head[expr] === Plus,
                 monomials = List @@ expr;,
                 monomials = {expr};
                ];
              For[i = 1, i <= Length[monomials], i++,
                  {coeff, term} = GetCoefficientAndTerm[monomials[[i]], dimZeroPars];
                  If[term =!= Null,
                     basis = Append[basis, term];,
                     Print["Error: cannot interpret settings for parameters ", InputForm[pars]];
                     Print["   because a boundary value cannot be detemined."];
                     Quit[1];
                    ];
                 ];
              basis = DeleteDuplicates[basis];
             ];
           basis
          ];

ListAllTermsOfForm[par_, termTypes_List] :=
    Module[{i, sets = {}, possibleTerms, terms = {}},
           For[i = 1, i <= Length[termTypes], i++,
               sets = termTypes[[i]] /. {f_, n_Integer?Positive} :> Sequence @@ ConstantArray[f, n];
               sets = With[{set = #}, {#[[1]], #[[2]]}& /@ set]& /@ sets;
               possibleTerms = Tuples[sets];
               possibleTerms = Flatten[With[{tuple = #}, Tuples[#[[2]]& /@ tuple]]& /@ possibleTerms, 1];
               possibleTerms = DeleteDuplicates[Apply[Times, possibleTerms, 1]];
               terms = Join[terms, possibleTerms];
              ];
           DeleteDuplicates[terms]
          ];

GetLinearSystemSolutions[pars_List, subs_List, nonHomogeneousTerms_List:{}] :=
    Module[{basis, solns = {}},
           basis = GetSolutionBasis[pars, subs];
           solns = SemiAnalyticSolution[#, basis]& /@ pars;
           If[nonHomogeneousTerms =!= {},
              solns = SemiAnalyticSolution[#[[1]], Join[#[[2]], ListAllTermsOfForm[#[[1]], nonHomogeneousTerms]]]& /@ solns;
             ];
           solns
          ];

GetBoundaryConditionsFor[pars_List, subs_List] := Select[subs, MemberQ[pars, Parameters`StripIndices[#[[1]]]]&];

GetScalarMassesSolutions[scalarMasses_List, boundaryConditionSubs_List, dimOneSolns_List] :=
    Module[{extraTerms},
           extraTerms = {{{dimOneSolns, 2}}};
           GetLinearSystemSolutions[scalarMasses, boundaryConditionSubs, extraTerms]
          ];

GetSoftBilinearsSolutions[softBilinears_List, boundaryConditionSubs_List, dimOneSolns_List] :=
    Module[{susyBilinears, extraTerms},
           susyBilinears = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaMuij;
           susyBilinears = Flatten[(# /. GetParameterExpansionRules[#])& /@ susyBilinears];
           (* dimensionful SUSY parameters have trivial semi-analytic solutions *)
           susyBilinears = SemiAnalyticSolution[#, {#}]& /@ susyBilinears;
           extraTerms = {{{dimOneSolns, 1}, {susyBilinears, 1}}};
           GetLinearSystemSolutions[softBilinears, boundaryConditionSubs, extraTerms]
          ];

GetSoftLinearsSolutions[softLinears_List, boundaryConditionSubs_List, dimOneSolns_List,
                        diracSolns_List, scalarMassesSolns_List, softBilinearsSolns_List] :=
    Module[{susyBilinears, susyLinears, extraTerms},
           susyBilinears = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaMuij;
           susyBilinears = Flatten[(# /. GetParameterExpansionRules[#])& /@ susyBilinears];
           susyBilinears = SemiAnalyticSolution[#, {#}]& /@ susyBilinears;
           susyLinears = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaLi;
           susyLinears = Flatten[(# /. GetParameterExpansionRules[#])& /@ susyLinears];
           susyLinears = SemiAnalyticSolution[#, {#}]& /@ susyLinears;
           extraTerms = {{{susyLinears, 1}, {dimOneSolns, 1}},
                         {{susyBilinears, 2}, {dimOneSolns, 1}},
                         {{susyBilinears, 1}, {dimOneSolns, 2}},
                         {{susyBilinears, 1}, {scalarMassesSolns, 1}},
                         {{susyBilinears, 1}, {softBilinears, 1}},
                         {{softBilinears, 1}, {dimOneSolns, 1}}
                        };
           If[diracSolns =!= {},
              extraTerms = Join[extraTerms,
                                {{{diracSolns, 1}, {scalarMassesSolns, 1}},
                                 {{diracSolns, 2}, {dimOneSolns, 1}},
                                 {{diracSolns, 2}, {susyBilinears, 1}}
                                }
                               ];
             ];
           GetLinearSystemSolutions[softLinears, boundaryConditionSubs, extraTerms]
          ];

GetSemiAnalyticSolutions[settings_List] :=
    Module[{boundaryValSubs, dimOnePars, dimOneBCs, dimOneSolns = {},
            diracPars, diracBCs, diracSolns = {}, scalarMasses, scalarMassesBCs,
            scalarMassesSolns = {}, softBilinears, softBilinearsBCs,
            softBilinearsSolns = {}, softLinears, softLinearBCs,
            softLinearSolns = {}, result = {}},

           boundaryValSubs = GetBoundaryValueSubstitutions[settings];

           dimOnePars = SelectParametersWithMassDimension[allSemiAnalyticParameters, 1];
           If[SARAH`SupersymmetricModel,
              (* do Dirac gaugino masses separately *)
              diracPars = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaDGi;
              diracPars = Select[dimOnePars, MemberQ[diracPars, #]&];
              diracBCs = GetBoundaryConditionsFor[diracPars, boundaryValSubs];
              dimOnePars = Complement[dimOnePars, diracPars];
             ];
           dimOneBCs = GetBoundaryConditionsFor[dimOnePars, boundaryValSubs];
           dimOneSolns = GetLinearSystemSolutions[dimOnePars, dimOneBCs];
           result = Join[result, dimOneSolns];

           If[diracPars =!= {},
              diracSolns = GetLinearSystemSolutions[diracPars, diracBCs];
              result = Join[result, diracSolns];
             ];

           scalarMasses = SelectParametersWithMassDimension[allSemiAnalyticParameters, 2];
           If[SARAH`SupersymmetricModel,
              softBilinears = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaBij;
              softBilinears = Select[scalarMasses, MemberQ[softBilinears, #]&];
              softBilinearsBCs = GetBoundaryConditionsFor[softBilinears, boundaryValSubs];
              scalarMasses = Complement[scalarMasses, softBilinears];
             ];
           scalarMassesBCs = GetBoundaryConditionsFor[scalarMasses, boundaryValSubs];
           scalarMassesSolns
               = GetScalarMassesSolutions[scalarMasses, scalarMassesBCs, dimOneSolns];
           result = Join[result, scalarMassesSolns];

           If[SARAH`SupersymmetricModel,
              softBilinearsSolns
                  = GetSoftBilinearsSolutions[softBilinears, softBilinearsBCs, dimOneSolns];
              result = Join[result, softBilinearsSolns];
             ];

           If[SARAH`SupersymmetricModel,
              softLinears = SelectParametersWithMassDimension[allSemiAnalyticParameters, 3];
              softLinearsBCs = GetBoundaryConditionsFor[softLinears, boundaryValSubs];
              softLinearsSolns
                  = GetSoftLinearsSolutions[softLinears, softLinearsBCs, dimOneSolns,
                                            diracSolns, scalarMassesSolns, softBilinearsSolns];
              result = Join[result, softLinearsSolns];
             ];

           result
          ];

CreateCoefficientNames[solution_SemiAnalyticSolution] :=
    Module[{par, basis, basisSize, i},
           par = CConversion`ToValidCSymbolString[GetName[solution]];
           basis = GetBasis[solution];
           basisSize = Length[basis];
	   Table[par <> "_coeff_" <> ToString[i], {i, 1, basisSize}]
          ];

CreateSemiAnalyticSolutionsDefinitions[solution_SemiAnalyticSolution] :=
    Module[{par, type, coeffs, defs = ""},
           par = GetName[solution];
           type = CConversion`CreateCType[Parameters`GetType[par]];
           coeffs = CreateCoefficientNames[solution];
           Utils`StringJoinWithSeparator[(type <> " " <> #)& /@ coeffs, ";\n"] <> ";\n"
          ];

CreateSemiAnalyticSolutionsDefinitions[solutions_List] :=
    Module[{def = ""},
           (def = def <> CreateSemiAnalyticSolutionsDefinitions[#])& /@ solutions;
           Return[def];
          ];

CreateSemiAnalyticSolutionDefaultInitialization[solution_SemiAnalyticSolution] :=
    Module[{par, type = "", coeffs, defaults},
           par = GetName[solution];
           type = Parameters`GetType[par];
	   coeffs = CreateCoefficientNames[solution];
           defaults = CConversion`CreateDefaultConstructor[#, type]& /@ coeffs;
	   "," <> Utils`StringJoinWithSeparator[defaults, ","]
          ];

CreateSemiAnalyticSolutionsInitialization[solutions_List] :=
    Module[{def = ""},
           (def = def <> CreateSemiAnalyticSolutionDefaultInitialization[#])& /@ solutions;
           Return[def];
          ];

DependsAtMostOn[num_?NumericQ, pars_List] := True;
DependsAtMostOn[True, pars_List] := True;
DependsAtMostOn[False, pars_List] := True;
DependsAtMostOn[sym_Symbol, pars_List] := MemberQ[pars, sym];
DependsAtMostOn[sym_[indices__] /; And @@ (Parameters`IsIndex /@ {indices}), pars_List] := DependsAtMostOn[sym, pars];
DependsAtMostOn[expr_Times, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_Plus, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_Power, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[Factorial[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Conjugate[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[SARAH`Conj[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Susyno`LieGroups`conj[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ConjugateTranspose[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Transpose[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[SARAH`Tp[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[SARAH`Adj[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Re[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Im[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Abs[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Arg[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Sign[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Round[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Floor[expr_, mult_], pars_List] := DependsAtMostOn[expr, pars] && DependsAtMostOn[mult, pars];
DependsAtMostOn[Floor[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Ceiling[expr_, mult_], pars_List] := DependsAtMostOn[expr, pars] && DependsAtMostOn[mult, pars];
DependsAtMostOn[Ceiling[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[expr_Min, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_Max, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_UnitStep, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_UnitBox, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_UnitTriangle, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_Piecewise, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[HoldPattern[(Equal|Unequal|Greater|Less|GreaterEqual|LessEqual)[expr__]], pars_List] := DependsAtMostOn[List[expr], pars];
DependsAtMostOn[HoldPattern[(And|Or|Xor|Nand|Nor)[expr__]], pars_List] := DependsAtMostOn[List[expr], pars];
DependsAtMostOn[Not[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[expr_Which, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_If, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_Switch, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[Log[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Sin[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Cos[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Tan[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Csc[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Sec[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Cot[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcSin[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcCos[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcTan[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcCsc[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcSec[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcCot[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Sinh[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Cosh[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Tanh[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Csch[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Sech[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Coth[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcSinh[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcCosh[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcTanh[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcCsch[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcSech[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcCoth[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[IsClose[first_, second_], pars_List] := DependsAtMostOn[first, pars] && DependsAtMostOn[second, pars];
DependsAtMostOn[IsClose[first_, second_, tol_], pars_List] :=
    DependsAtMostOn[first, pars] && DependsAtMostOn[second, pars] && DependsAtMostOn[tol, pars];
DependsAtMostOn[IsCloseRel[first_, second_], pars_List] := DependsAtMostOn[first, pars] && DependsAtMostOn[second, pars];
DependsAtMostOn[IsCloseRel[first_, second_, tol_], pars_List] :=
    DependsAtMostOn[first, pars] && DependsAtMostOn[second, pars] && DependsAtMostOn[tol, pars];
DependsAtMostOn[expr_List, pars_List] := And @@ (DependsAtMostOn[#, pars]& /@ expr);
DependsAtMostOn[expr_, pars_List] := False;

End[];

EndPackage[];

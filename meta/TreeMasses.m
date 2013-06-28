
BeginPackage["TreeMasses`", {"SARAH`", "TextFormatting`", "CConversion`", "Parameters`"}];

FSMassMatrix::usage="Head of a mass matrix";

ConvertSarahMassMatrices::usage="creates list of mass matrices using
SARAH's MassMatrix[] function";

CreateMassGetter::usage="creates a C function for
the mass getter";

CreateMixingMatrixGetter::usage="creates a getter for the mixing
matrix";

CreateMassCalculationPrototype::usage="creates a C function prototype
from a mass matrix";

CreateMassCalculationFunction::usage="creates a C function that
calculates the mass eigenstates by diagonalizing the mass matrix";

CallMassCalculationFunction::usage="creates a C function call of a
mass matrix calcualtion function";

CreatePhysicalMassDefinition::usage="creates definition of physical
mass.";

CreatePhysicalMassInitialization::usage="creates default
initialization of the given mass eigenstate";

CreateMixingMatrixDefinition::usage="creates definition of mixing
matrix";

CreateMixingMatrixInitialization::usage="creates default
initialization of the given mixing matrix";

ClearOutputParameters::usage="clears masses and mixing matrices";

GetParticles::usage="returns list of particles";

GetSusyParticles::usage="returns list of susy particles";

GetSMParticles::usage="returns list of Standard Model particles";

GetDimension::usage="returns the size of the particle multiplet";

GetDimensionStartSkippingGoldstones::usage="return first index,
skipping goldstone bosons";

FindMixingMatrixSymbolFor::usage="returns the mixing matrix symbol for
a given field";

SetUnrotatedParticles::usage="set list of unrotated particles in SARAH
format (result of ListUnmixed[EWSB] or content of
UnrotatedParticles.m)";

SetModelParameters::usage="set model parameters";

GetMassEigenstate::usage="get mass eigenstates symbol from mass
matrix";

GetMixingMatrixSymbol::usage="get mixing matrix symbol from mass
matrix";

ReplaceDependencies::usage="returs expression with dependencies
(ThetaW etc.) replaced by the user-defined expressions (";

FindColorGaugeGroup::usage="returns triplet of color gauge coupling,
group and SARAH name";

FindColorGaugeCoupling::usage="returns symbol of color gauge coupling";
FindLeftGaugeCoupling::usage="returns symbol of weak (left) gauge coupling";
FindHyperchargeGaugeCoupling::usage="returns symbol of hypercharge gauge coupling";

CreateDependenceNumPrototypes::usage="";
CreateDependenceNumFunctions::usage="";

IsScalar::usage="";
IsFermion::usage="";
IsVector::usage="";
IsGhost::usage="";
IsAuxiliary::usage="";
IsVEV::usage="";
IsMajoranaFermion::usage="";
IsDiracFermion::usage="";
IsComplexScalar::usage="";
IsRealScalar::usage="";
IsMassless::usage="";

Begin["Private`"];

allModelParameters = {};

SetModelParameters[pars_List] := allModelParameters = pars;

unrotatedParticles = {};

SetUnrotatedParticles[list_List] :=
    unrotatedParticles = ({#[[1]], #[[4]]})& /@ list;

GetVectorBosons[states_:SARAH`EWSB] :=
    #[[1]]& /@ Cases[SARAH`Particles[states], {_,_,_,V,__}] /.
       SARAH`diracSubBack1[SARAH`ALL] /.
       SARAH`diracSubBack2[SARAH`ALL];

(* Create list of mass eigenstate particles *)
GetParticles[states_:SARAH`EWSB] :=
    Module[{particles = {}},
           particles = (GetHead[#[[1,1]]])& /@ SARAH`Masses[states];
           particles = particles /.
                       SARAH`diracSubBack1[SARAH`ALL] /.
                       SARAH`diracSubBack2[SARAH`ALL];
           particles = Join[particles, GetVectorBosons[states]];
           particles = Select[particles, (!IsOfType[#,SARAH`NoField])&];
           Return[DeleteDuplicates[particles]];
          ];

GetSusyParticles[states_:SARAH`EWSB] :=
    Select[GetParticles[states], (!SARAH`SMQ[#])&];

GetSMParticles[states_:SARAH`EWSB] :=
    Select[GetParticles[states], (SARAH`SMQ[#])&];

IsOfType[sym_Symbol, type_Symbol, states_:SARAH`EWSB] :=
    SARAH`getType[sym, False, states] === type;

IsOfType[sym_[__], type_Symbol, states_:SARAH`EWSB] :=
    IsOfType[sym, type, states];

IsScalar[sym_Symbol] := IsOfType[sym, S];

IsFermion[sym_Symbol] := IsOfType[sym, F];

IsVector[sym_Symbol] := IsOfType[sym, V];

IsGhost[sym_Symbol] := IsOfType[sym, G];

IsAuxiliary[sym_Symbol] := IsOfType[sym, A];

IsVEV[sym_Symbol] := IsOfType[sym, VEV];

IsMajoranaFermion[sym_Symbol] :=
    And[IsFermion[sym], MemberQ[SARAH`MajoranaPart, sym]];

IsDiracFermion[sym_Symbol] :=
    And[IsFermion[sym], !MemberQ[SARAH`MajoranaPart, sym]];

IsComplexScalar[sym_Symbol] :=
    And[IsScalar[sym], Parameters`IsComplexParameter[sym]];

IsRealScalar[sym_Symbol] :=
    And[IsScalar[sym], Parameters`IsRealParameter[sym]];

IsMassless[sym_Symbol, states_:SARAH`EWSB] :=
    MemberQ[SARAH`Massless[states], sym];

GetDimension[sym_[__], states_:SARAH`EWSB] := GetDimension[sym, states];

GetDimension[sym_Symbol, states_:SARAH`EWSB] :=
    SARAH`getGen[sym, states];

GetDimensionStartSkippingGoldstones[sym_[__], states_:SARAH`EWSB] :=
    GetDimensionStartSkippingGoldstones[sym, states];

GetDimensionStartSkippingGoldstones[sym_, states_:SARAH`EWSB] :=
    Module[{goldstones, max = 1},
           goldstones = Transpose[SARAH`GoldstoneGhost][[2]];
           If[FreeQ[goldstones, sym],
              Return[1];,
              If[GetDimension[sym] === 1,
                 Return[2];,
                 While[!FreeQ[goldstones, sym[{max}]],
                       max++];
                 Return[max];
                ];
             ];
          ];

FindMixingMatrixSymbolFor[particle_Symbol] :=
    Module[{k, i, l, mixingMatrixSymbol = Null, mixingList = {}, mixingScheme = {},
            currentName},
           For[k = 1, k <= Length[NameOfStates], k++,
               If[Head[DEFINITION[NameOfStates[[k]]][MatterSector]] === List,
                  mixingList = DEFINITION[NameOfStates[[k]]][MatterSector];
                  For[i = 1, i <= Length[mixingList], i++,
                      If[Length[mixingList[[i]]] != 2, Continue[]];
                      mixingScheme = mixingList[[i,2]];
                      If[Length[mixingScheme] != 2, Continue[]];
                      If[Head[mixingScheme[[1]]] === Symbol,
                         currentName = mixingScheme[[1]] /.
                            SARAH`diracSubBack1[NameOfStates[[k]]] /.
                            SARAH`diracSubBack2[NameOfStates[[k]]];
                         If[currentName === particle,
                            mixingMatrixSymbol = mixingScheme[[2]];
                            Return[mixingMatrixSymbol];
                           ];
                         ,
                         currentName = mixingScheme[[1,1]] /.
                            SARAH`diracSubBack1[NameOfStates[[k]]] /.
                            SARAH`diracSubBack2[NameOfStates[[k]]];
                         If[currentName === particle,
                            mixingMatrixSymbol = {mixingScheme[[1,2]], mixingScheme[[2,2]]};
                            Return[mixingMatrixSymbol];
                           ];
                        ];
                     ];
                 ];
               If[Head[DEFINITION[NameOfStates[[k]]][GaugeSector]] === List,
                  mixingList = DEFINITION[NameOfStates[[k]]][GaugeSector];
                  For[i = 1, i <= Length[mixingList], i++,
                      mixingScheme = mixingList[[i]];
                      If[Length[mixingScheme] != 3, Continue[]];
                      If[MemberQ[mixingScheme[[2]], particle],
                         mixingMatrixSymbol = mixingScheme[[3]];
                         Return[mixingMatrixSymbol];
                        ];
                     ];
                 ];
              ];
           Return[mixingMatrixSymbol];
          ];

IsUnmixed[particle_Symbol] :=
    MemberQ[(#[[1]])& /@ unrotatedParticles, particle];

GetMassOfUnmixedParticle[particle_Symbol] :=
    Cases[unrotatedParticles, {particle, _}][[1,2]];

ConvertSarahMassMatrices[] :=
    Module[{particles = {}, result = {}, eigenstateName, massMatrix,
            gaugeDefs = {}, gaugeMassES = {}, multiplet = {},
            k, multipletName, rules = {}, mixingMatrixSymbol, dim},
           particles = GetParticles[];
           For[k = 1, k <= Length[particles], k++,
               eigenstateName = particles[[k]];
               If[IsUnmixed[eigenstateName],
                  massMatrix = { ReplaceDependencies[GetMassOfUnmixedParticle[eigenstateName]] };
                  AppendTo[result, TreeMasses`FSMassMatrix[massMatrix, eigenstateName, Null]];
                  ,
                  massMatrix = ReplaceDependencies[SARAH`MassMatrix[eigenstateName]];
                  mixingMatrixSymbol = FindMixingMatrixSymbolFor[eigenstateName];
                  If[Head[massMatrix] === List,
                     AppendTo[result, TreeMasses`FSMassMatrix[massMatrix, eigenstateName, mixingMatrixSymbol]];
                    ];
                 ];
              ];
           Return[result];
          ];

GetMixingMatrixSymbol[massMatrix_TreeMasses`FSMassMatrix] := massMatrix[[3]];

GetMassEigenstate[massMatrix_TreeMasses`FSMassMatrix] := massMatrix[[2]];

GetMassMatrix[massMatrix_TreeMasses`FSMassMatrix] := massMatrix[[1]];

GetMixingMatrixType[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{type, eigenstate, mixingMatrixSymbol, dim},
           eigenstate = GetMassEigenstate[massMatrix];
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           dim = Length[GetMassMatrix[massMatrix]];
           Which[Parameters`IsRealParameter[mixingMatrixSymbol], type = "DoubleMatrix";,
                 IsFermion[eigenstate],                     type = "ComplexMatrix";,
                 True,                                      type = "DoubleMatrix";
                ];
           Return[CConversion`MatrixType[type, dim, dim]];
          ];

CreateMassGetter[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{massESSymbol, returnType, dim, massESSymbolStr},
           massESSymbol = GetMassEigenstate[massMatrix];
           massESSymbolStr = ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]];
           dim = GetDimension[massESSymbol];
           If[dim == 1,
              returnType = CConversion`ScalarType["double"];,
              returnType = CConversion`VectorType["DoubleVector", dim];
             ];
           CConversion`CreateInlineGetter[massESSymbolStr, returnType]
          ];

CreateMixingMatrixGetter[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{mixingMatrixSymbol, returnType},
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           returnType = GetMixingMatrixType[massMatrix];
           CreateMixingMatrixGetter[mixingMatrixSymbol, returnType]
          ];

CreateMixingMatrixGetter[mixingMatrixSymbol_List, returnType_] :=
    Module[{result = ""},
           (result = result <> CreateMixingMatrixGetter[#,returnType])& /@ mixingMatrixSymbol;
           Return[result];
          ];

CreateMixingMatrixGetter[Null, returnType_] := "";

CreateMixingMatrixGetter[mixingMatrixSymbol_Symbol, returnType_] :=
    CConversion`CreateInlineGetter[ToValidCSymbolString[mixingMatrixSymbol], returnType];

CreateMassCalculationPrototype[TreeMasses`FSMassMatrix[_, massESSymbol_, Null]] :=
    Module[{result, ev = ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]]},
           result = "void calculate_" <> ev <> "();\n";
           Return[result];
          ];

CreateMassCalculationPrototype[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result = "", massESSymbol},
           massESSymbol = GetMassEigenstate[massMatrix];
           result = CreateMassMatrixGetterPrototype[massMatrix] <>
                    "void calculate_" <> ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]] <>
                    "();\n";
           Return[result];
          ];

CallMassCalculationFunction[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result = "", k, massESSymbol},
           massESSymbol = GetMassEigenstate[massMatrix];
           result = "calculate_" <> ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]]
                    <> "();\n";
           Return[result];
          ];

(* These rules are needed to check if a matrix is hermitian *)
SARAH`sum /: Susyno`LieGroups`conj[SARAH`sum[ind_,a_,b_,expr_]] := SARAH`sum[ind,a,b,Susyno`LieGroups`conj[expr]];

Susyno`LieGroups`conj[m_[a_,b_]] := m[b,a] /; MemberQ[SARAH`ListSoftBreakingScalarMasses, m];

IsSymmetric[matrix_List] := IsHermitian[matrix, Identity];

IsHermitian[matrix_List, op_:Susyno`LieGroups`conj] :=
    Module[{rows, cols, i, k},
           rows = Length[matrix];
           For[i = 1, i <= rows, i++,
               cols = Length[matrix[[i]]];
               If[rows =!= cols, Return[False];];
               For[k = 1, k < i, k++,
                   If[matrix[[i,k]] =!= op[matrix[[k,i]]], Return[False];];
                  ];
              ];
           Return[True];
          ];

MatrixToCFormString[matrix_List, symbol_String, matrixType_String:"DoubleMatrix"] :=
    Module[{dim, result = "", i, k, isSymmetric = IsSymmetric[matrix],
            isHermitian = IsHermitian[matrix]},
           dim = Length[matrix];
           result = matrixType <> " " <> symbol <> "(" <> ToString[dim] <>
                    "," <> ToString[dim] <> ");\n";
           For[i = 1, i <= dim, i++,
               For[k = 1, k <= dim, k++,
                   result = result <> symbol <> "(" <> ToString[i] <>
                            "," <> ToString[k] <> ") = ";
                   Which[isSymmetric && i > k,
                         result = result <> symbol <> "(" <> ToString[k] <>
                                  "," <> ToString[i] <> ");\n"
                         ,
                         isHermitian && i > k,
                         result = result <> "Conj(" <> symbol <> "(" <> ToString[k] <>
                                  "," <> ToString[i] <> "));\n"
                         ,
                         True,
                         result = result <>
                                  RValueToCFormString[matrix[[i,k]]] <> ";\n";
                     ];
                  ];
              ];
           Return[result];
          ];

CreateMassMatrixGetterFunction[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result, body, ev, matrixSymbol, matrix, massESSymbol},
           massESSymbol = GetMassEigenstate[massMatrix];
           ev = ToValidCSymbolString[GetHead[massESSymbol]];
           matrixSymbol = "mass_matrix_" <> ev;
           matrix = GetMassMatrix[massMatrix];
           body = MatrixToCFormString[matrix, matrixSymbol] <> "\n";
           result = "DoubleMatrix CLASSNAME::get_" <> matrixSymbol <> "() const\n{\n" <>
                    IndentText[body] <>
                    "return " <> matrixSymbol <> ";\n}\n";
           Return[result];
          ];

CreateMassMatrixGetterPrototype[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result, ev, matrixSymbol, massESSymbol},
           massESSymbol = GetMassEigenstate[massMatrix];
           ev = ToValidCSymbolString[GetHead[massESSymbol]];
           matrixSymbol = "mass_matrix_" <> ev;
           result = "DoubleMatrix get_" <> matrixSymbol <> "() const;\n";
           Return[result];
          ];

CreateDiagonalizationFunction[matrix_List, eigenVector_, mixingMatrixSymbol_] :=
    Module[{dim, body = "", result, U = "", V = "", dimStr = "", ev = "", k},
           dim = Length[matrix];
           dimStr = ToString[dim];
           ev = ToValidCSymbolString[GetHead[eigenVector]];
           matrixSymbol = "mass_matrix_" <> ev;
           ev = ToValidCSymbolString[FlexibleSUSY`M[GetHead[eigenVector]]];
           result = "void CLASSNAME::calculate_" <> ev <> "()\n{\n";
           body = "DoubleMatrix " <> matrixSymbol <> "(get_" <> matrixSymbol <> "());\n";
           If[Head[mixingMatrixSymbol] === List && Length[mixingMatrixSymbol] == 2,
              (* use SVD *)
              U = ToValidCSymbolString[mixingMatrixSymbol[[1]]];
              V = ToValidCSymbolString[mixingMatrixSymbol[[2]]];
              If[dim == 2,
                 body = body <> "Diagonalize2by2(" <>
                        matrixSymbol <> ", " <> U <> ", " <> V <> ", " <> ev <> ");\n";
                 ,
                 body = body <> "Diagonalize(" <>
                        matrixSymbol <> ", " <> U <> ", " <> V <> ", " <> ev <> ");\n";
                ];
              ,
              (* use conventional diagonalization *)
              U = ToValidCSymbolString[mixingMatrixSymbol];
              If[dim == 2,
                 body = body <> "Diagonalize2by2(" <> matrixSymbol <> ", " <>
                        U <> ", " <> ev <> ");\n";
                 ,
                 body = body <> "Diagonalize(" <> matrixSymbol <> ", " <>
                        U <> ", " <> ev <> ");\n";
                ];
             ];
           If[IsScalar[eigenVector] || IsVector[eigenVector],
              (* check for tachyons *)
              body = body <> "\nint min_element;\n" <>
                     "if (" <> ev <> ".min(min_element) < 0.)\n" <>
                     IndentText["throw TachyonError(this, \"" <> ev <> "\", min_element);\n\n"];
              body = body <> ev <> " = " <> ev <> ".apply(ZeroSqrt);\n";
             ];
           Return[result <> IndentText[body] <> "}\n"];
          ];

CreateMassCalculationFunction[TreeMasses`FSMassMatrix[mass_, massESSymbol_, Null]] :=
    Module[{result, ev = ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]], body,
            trans = Identity},
           result = "void CLASSNAME::calculate_" <> ev <> "()\n{\n";
           If[IsVector[massESSymbol] || IsScalar[massESSymbol],
              trans = Sqrt;
             ];
           body = ev <> " = " <>
                  RValueToCFormString[trans[mass[[1]]]] <> ";\n";
           body = IndentText[body];
           Return[result <> body <> "}\n\n"];
          ];

CreateMassCalculationFunction[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result = "", massESSymbol, mixingMatrixSymbol, matrix},
           massESSymbol = GetMassEigenstate[massMatrix];
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           matrix = GetMassMatrix[massMatrix];
           result = result <>
                    CreateMassMatrixGetterFunction[massMatrix] <> "\n" <>
                    CreateDiagonalizationFunction[matrix, massESSymbol,
                                                  mixingMatrixSymbol]
                    <> "\n";
           Return[result];
          ];

CreatePhysicalMassDefinition[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result = "", massESSymbol, returnType = "DoubleVector"},
           massESSymbol = GetMassEigenstate[massMatrix];
           If[GetDimension[massESSymbol] == 1,
              returnType = "double";
             ];
           result = returnType <> " " <>
                    ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]] <> ";\n";
           Return[result];
          ];

CreatePhysicalMassInitialization[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result = "", massESSymbol, dim},
           massESSymbol = GetMassEigenstate[massMatrix];
           dim = GetDimension[massESSymbol];
           result = ", " <> ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]];
           If[dim == 1,
              result = result <> "(0)";,
              result = result <> "(" <> ToString[dim] <> ")";
             ];
           Return[result];
          ];

DefineMatrix[Null, _] := "";

DefineMatrix[matrix_Symbol, type_String] :=
    type <> " " <> ToValidCSymbolString[matrix] <> ";\n";

DefineMatrix[matrix_List, type_String] :=
    Module[{result = ""},
           (result = result <> DefineMatrix[#,type])& /@ matrix;
           Return[result];
          ];

CreateMixingMatrixDefinition[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result, mixingMatrixSymbol, matrixType},
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           matrixType = GetMixingMatrixType[massMatrix][[1]];
           result = DefineMatrix[mixingMatrixSymbol, matrixType];
           Return[result];
          ];

ClearOutputParameters[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result, massESSymbol, mixingMatrixSymbol, matrixType, dim, dimStr, i},
           massESSymbol = GetMassEigenstate[massMatrix];
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           dim = GetDimension[massESSymbol];
           dimStr = ToString[dim];
           If[dim == 1,
              result = ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]] <> " = 0.0;\n";
              ,
              result = ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]] <> " = DoubleVector(" <> dimStr <> ");\n";
             ];
           If[mixingMatrixSymbol =!= Null,
              matrixType = GetCParameterType[GetMixingMatrixType[massMatrix]];
              If[Head[mixingMatrixSymbol] === List,
                 For[i = 1, i <= Length[mixingMatrixSymbol], i++,
                     result = result <> ToValidCSymbolString[mixingMatrixSymbol[[i]]] <>
                              " = " <> matrixType <> "(" <> dimStr <> "," <> dimStr <> ");\n";
                    ];
                 ,
                 result = result <> ToValidCSymbolString[mixingMatrixSymbol] <>
                          " = " <> matrixType <> "(" <> dimStr <> "," <> dimStr <> ");\n";
                ];
             ];
           Return[result];
          ];

InitializeMatrix[Null, _] := "";

InitializeMatrix[matrix_Symbol, dim_Integer] :=
    Module[{},
           Return[", " <> ToValidCSymbolString[matrix] <> "(" <> ToString[dim] <> "," <> ToString[dim] <> ")"];
          ];

InitializeMatrix[matrix_List, dim_Integer] :=
    Module[{result = ""},
           (result = result <> InitializeMatrix[#, dim])& /@ matrix;
           Return[result];
          ];

CreateMixingMatrixInitialization[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result, mixingMatrixSymbol, dim},
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           dim = Length[GetMassMatrix[massMatrix]];
           result = InitializeMatrix[mixingMatrixSymbol, dim];
           Return[result];
          ];

FindColorGaugeGroup[] :=
    Module[{coupling, gaugeGroup, result},
           coupling = FindColorGaugeCoupling[];
           gaugeGroup = Cases[SARAH`Gauge, {_, group_, name_, coupling, ___}];
           If[gaugeGroup === {},
              Print["Error: could not find color gauge group"];
              result = Null;
              ,
              result = {coupling, gaugeGroup[[1,3]], gaugeGroup[[1,2]]};
             ];
           result
          ];

FindLeftGaugeGroup[] :=
    Module[{coupling, gaugeGroup, result},
           coupling = FindHyperchargeGaugeCoupling[];
           gaugeGroup = Cases[SARAH`Gauge, {_, group_, name_, coupling, ___}];
           If[gaugeGroup === {},
              Print["Error: could not weak gauge group"];
              result = Null;
              ,
              result = {coupling, gaugeGroup[[1,3]], gaugeGroup[[1,2]]};
             ];
           result
          ];

FindHyperchargeGaugeGroup[] :=
    Module[{coupling, gaugeGroup, result},
           coupling = FindColorGaugeCoupling[];
           gaugeGroup = Cases[SARAH`Gauge, {_, group_, name_, coupling, ___}];
           If[gaugeGroup === {},
              Print["Error: could not find Hypercharge gauge group"];
              result = Null;
              ,
              result = {coupling, gaugeGroup[[1,3]], gaugeGroup[[1,2]]};
             ];
           result
          ];

FindColorGaugeCoupling[] := SARAH`strongCoupling;

FindLeftGaugeCoupling[] := SARAH`leftCoupling;

FindHyperchargeGaugeCoupling[] := SARAH`hyperchargeCoupling;

dependenceNumsUpToDate = False;
dependenceNumRulesUpToDate = False;
dependenceNums = {}; (* replacement rules for all DependenceNum *)
dependenceNumRules = {}; (* replacement rules for all DependenceNum *)

FindDependenceNums[] :=
    Module[{hyperchargeCoupling, leftCoupling},
           If[!dependenceNumsUpToDate,
              hyperchargeCoupling = FindHyperchargeGaugeCoupling[];
              leftCoupling = FindLeftGaugeCoupling[];
              dependenceNums = Join[
                  { Rule[SARAH`Weinberg,
                         ArcSin[hyperchargeCoupling / Sqrt[hyperchargeCoupling^2 + leftCoupling^2]] /.
                         Parameters`ApplyGUTNormalization[]] },
                  Cases[SARAH`ParameterDefinitions,
                        {parameter_ /; !MemberQ[allModelParameters, parameter] &&
                         parameter =!= SARAH`Weinberg && parameter =!= SARAH`electricCharge,
                         {___, SARAH`DependenceNum -> value:Except[None], ___}} :>
                        Rule[parameter, value /. Parameters`ApplyGUTNormalization[]]]
                                   ];
              dependenceNumsUpToDate = True;
             ];
           dependenceNums
          ];

FindDependenceNumRules[] :=
    Module[{hyperchargeCoupling, leftCoupling},
           If[!dependenceNumRulesUpToDate,
              hyperchargeCoupling = FindHyperchargeGaugeCoupling[];
              leftCoupling = FindLeftGaugeCoupling[];
              dependenceNumRules = Join[
                  { SARAH`Weinberg -> SARAH`Weinberg[] },
                  Cases[SARAH`ParameterDefinitions,
                        {parameter_ /; !MemberQ[allModelParameters, parameter] &&
                         parameter =!= SARAH`Weinberg && parameter =!= SARAH`electricCharge,
                         {___, SARAH`DependenceNum -> value:Except[None], ___}} :>
                        Rule[parameter, parameter[]]]
                                   ];
              dependenceNumRulesUpToDate = True;
             ];
           dependenceNumRules
          ];

CreateDependenceNumPrototype[Rule[parameter_, _]] :=
    "double " <> ToValidCSymbolString[parameter] <> "() const;\n";

CreateDependenceNumPrototypes[] :=
    Module[{dependenceNums, result = ""},
           dependenceNums = FindDependenceNums[];
           (result = result <> CreateDependenceNumPrototype[#])& /@ dependenceNums;
           Return[result];
          ];

CreateDependenceNumFunction[Rule[parameter_, value_]] :=
    Module[{result, body, parStr},
           parStr = ToValidCSymbolString[parameter];
           body = "return " <> RValueToCFormString[value] <> ";\n";
           result = "double CLASSNAME::" <> parStr <> "() const\n{\n" <>
                    IndentText[body] <> "}\n\n";
           Return[result];
          ];

CreateDependenceNumFunctions[] :=
    Module[{dependenceNums, result = ""},
           dependenceNums = FindDependenceNums[];
           (result = result <> CreateDependenceNumFunction[#])& /@ dependenceNums;
           Return[result];
          ];

ReplaceDependencies[expr_] :=
    expr /. FindDependenceNumRules[];

End[];

EndPackage[];

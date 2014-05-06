
BeginPackage["TreeMasses`", {"SARAH`", "TextFormatting`", "CConversion`", "Parameters`", "WeinbergAngle`"}];

FSMassMatrix::usage="Head of a mass matrix";

ConvertSarahMassMatrices::usage="creates list of mass matrices using
SARAH's MassMatrix[] function";

CreateMassGetter::usage="creates a C function for
the mass getter";

CreateParticleLaTeXNames::usage="creates a list of the particle's
LaTeX names";

CreateParticleNames::usage="creates a list of the particle's names";

CreateParticleEnum::usage="creates an enum of the particles";

CreateParticleMultiplicity::usage="creates array of the particle
multiplicities";

FillSpectrumVector::usage="";

CreateMixingMatrixGetter::usage="creates a getter for the mixing
matrix";

CreateMassCalculationPrototype::usage="creates a C function prototype
from a mass matrix";

CreateMassCalculationFunction::usage="creates a C function that
calculates the mass eigenstates by diagonalizing the mass matrix";

CallMassCalculationFunctions::usage="creates C function calls of all
mass matrix calcualtion functions";

CreatePhysicalMassDefinition::usage="creates definition of physical
mass.";

CreatePhysicalMassInitialization::usage="creates default
initialization of the given mass eigenstate";

CreateMixingMatrixDefinition::usage="creates definition of mixing
matrix";

CreateMixingMatrixInitialization::usage="creates default
initialization of the given mixing matrix";

ClearOutputParameters::usage="clears masses and mixing matrices";

CopyDRBarMassesToPoleMasses::usage="copies DRbar mass to pole mass";

GetParticles::usage="returns list of particles";

GetSusyParticles::usage="returns list of susy particles";

GetSMParticles::usage="returns list of Standard Model particles";

GetGoldstoneBosons::usage="returns list of all goldstone bosons";

GetSMGoldstoneBosons::usage="returns list of all Standard Model
goldstone bosons";

GetVectorBosons::usage="returns list of all vector bosons";

GetDimension::usage="returns the size of the particle multiplet";

GetDimensionStartSkippingGoldstones::usage="return first index,
skipping goldstone bosons";

FindMixingMatrixSymbolFor::usage="returns the mixing matrix symbol for
a given field";

SetUnrotatedParticles::usage="set list of unrotated particles in SARAH
format (result of ListUnmixed[EWSB] or content of
UnrotatedParticles.m)";

GetMassEigenstate::usage="get mass eigenstates symbol from mass
matrix";

GetMixingMatrixSymbol::usage="get mixing matrix symbol from mass
matrix";

GetMassOfUnmixedParticle::usage="returns mass of unmixed particle";

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
IsGolstone::usage="";
IsAuxiliary::usage="";
IsVEV::usage="";
IsMajoranaFermion::usage="";
IsDiracFermion::usage="";
IsComplexScalar::usage="";
IsRealScalar::usage="";
IsMassless::usage="";
IsUnmixed::usage="";

StripGenerators::usage="removes all generators Lam, Sig, fSU2, fSU3
and removes Delta with the given indices";

CreateThirdGenerationHelpers::usage="";
CallThirdGenerationHelperFunctionName::usage="";

Begin["`Private`"];

unrotatedParticles = {};

SetUnrotatedParticles[list_List] :=
    unrotatedParticles = ({#[[1]], #[[4]]})& /@ list;

GetVectorBosons[states_:FlexibleSUSY`FSEigenstates] :=
    #[[1]]& /@ Cases[SARAH`Particles[states], {_,_,_,V,__}] /.
       SARAH`diracSubBack1[SARAH`ALL] /.
       SARAH`diracSubBack2[SARAH`ALL];

(* Create list of mass eigenstate particles *)
GetParticles[states_:FlexibleSUSY`FSEigenstates] :=
    Module[{particles = {}},
           particles = (GetHead[#[[1,1]]])& /@ SARAH`Masses[states];
           particles = particles /.
                       SARAH`diracSubBack1[SARAH`ALL] /.
                       SARAH`diracSubBack2[SARAH`ALL];
           particles = Join[particles, GetVectorBosons[states]];
           particles = Select[particles, (!IsOfType[#,SARAH`NoField])&];
           Return[DeleteDuplicates[particles]];
          ];

GetSusyParticles[states_:FlexibleSUSY`FSEigenstates] :=
    Select[GetParticles[states], (!SARAH`SMQ[#] && !IsGhost[#])&];

GetSMParticles[states_:FlexibleSUSY`FSEigenstates] :=
    Select[GetParticles[states], (SARAH`SMQ[#])&];

IsOfType[sym_Symbol, type_Symbol, states_:FlexibleSUSY`FSEigenstates] :=
    SARAH`getType[sym, False, states] === type;

IsOfType[sym_[__], type_Symbol, states_:FlexibleSUSY`FSEigenstates] :=
    IsOfType[sym, type, states];

IsScalar[sym_Symbol] := IsOfType[sym, S];

IsFermion[sym_Symbol] := IsOfType[sym, F];

IsVector[sym_Symbol] := IsOfType[sym, V];

IsGhost[sym_Symbol] := IsOfType[sym, G];

IsGolstone[sym_] := MemberQ[GetGoldstoneBosons[] /. a_[{idx__}] :> a[idx], sym];

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

IsMassless[sym_Symbol, states_:FlexibleSUSY`FSEigenstates] :=
    MemberQ[SARAH`Massless[states], sym];

(* Returns list of pairs {p,v}, where p is the given golstone
   boson and v is the corresponding vector boson.

   Example (MSSM):
     GetCorrespondingVectorBosons[Ah]      ->  {{Ah[1], VZ}}
     GetCorrespondingVectorBosons[Ah[1]]   ->  {{Ah[1], VZ}}
     GetCorrespondingVectorBosons[Ah[{1}]] ->  {{Ah[1], VZ}}
*)
GetCorrespondingVectorBosons[goldstone_[idx_Integer]] :=
    GetCorrespondingVectorBosons[goldstone[{idx}]];

GetCorrespondingVectorBosons[goldstone_[idx_Symbol]] :=
    GetCorrespondingVectorBosons[goldstone[{idx}]];

GetCorrespondingVectorBosons[goldstone_] :=
    Module[{vector, idx, sym, association},
           association = Cases[SARAH`GoldstoneGhost, {vector_, goldstone | goldstone[{idx_}]}];
           Reverse /@ association /. sym_[{idx_}] :> sym[idx]
          ];

GetGoldstoneBosons[] :=
    Transpose[SARAH`GoldstoneGhost][[2]];

GetSMGoldstoneBosons[] :=
    Cases[SARAH`GoldstoneGhost, {vector_?SARAH`SMQ, goldstone_} :> goldstone];

GetDimension[sym_[__], states_:FlexibleSUSY`FSEigenstates] := GetDimension[sym, states];

GetDimension[sym_Symbol, states_:FlexibleSUSY`FSEigenstates] :=
    SARAH`getGen[sym, states];

GetDimensionStartSkippingGoldstones[sym_[__]] :=
    GetDimensionStartSkippingGoldstones[sym];

GetDimensionStartSkippingGoldstones[sym_] :=
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

(* Removes generators and Delta with the given indices.
 * Especially the following replacements are done:
 *
 * SARAH`Lam[__] -> 2    corresponds to SU(3) generator T[__] -> 1
 * SARAH`Sig[__] -> 2    corresponds to SU(2) generator T[__] -> 1
 * SARAH`fSU2[__] -> 1   SU(2) structure function
 * SARAH`fSU3[__] -> 1   SU(3) structure function
 * SARAH`Delta[a,b] -> 1 where a and b are from the `indices' list
 *)
StripGenerators[expr_, indices_List] :=
    Module[{headers = {SARAH`Delta}, h, indexCombinations, removeSymbols = {}},
           indexCombinations = DeleteCases[Subsets[indices],{}];
           For[h = 1, h <= Length[headers], h++,
               AppendTo[removeSymbols, (headers[[h]] @@ #)& /@ indexCombinations];
              ];
           removeSymbols = (Rule[#, 1])& /@ Flatten[removeSymbols];
           removeSymbols = Join[removeSymbols,
                                { Rule[SARAH`Lam[__], 2],
                                  Rule[SARAH`Sig[__], 2],
                                  Rule[SARAH`fSU2[__], 1],
                                  Rule[SARAH`fSU3[__], 1] }];
           expr /. removeSymbols
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
           (* the ghost masses will later be replaced explicitely by the
              corresponding vector boson masses *)
           particles = Select[GetParticles[], (!IsGhost[#])&];
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
           Which[Parameters`IsRealParameter[mixingMatrixSymbol],
                 type = CConversion`realScalarCType;,
                 IsFermion[eigenstate],
                 type = CConversion`complexScalarCType;,
                 True,
                 type = CConversion`realScalarCType;
                ];
           Return[CConversion`MatrixType[type, dim, dim]];
          ];

CreateMassGetter[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{massESSymbol, returnType, dim, dimStr, massESSymbolStr},
           massESSymbol = GetMassEigenstate[massMatrix];
           massESSymbolStr = ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]];
           dim = GetDimension[massESSymbol];
           dimStr = ToString[dim];
           If[dim == 1,
              returnType = CConversion`ScalarType[CConversion`realScalarCType];,
              returnType = CConversion`ArrayType[CConversion`realScalarCType, dim];
             ];
           CConversion`CreateInlineGetter[massESSymbolStr, returnType]
          ];

CreateParticleEnum[particles_List] :=
    Module[{i, par, name, result = ""},
           For[i = 1, i <= Length[particles], i++,
               par = particles[[i]];
               name = CConversion`ToValidCSymbolString[par];
               If[i > 1, result = result <> ", ";];
               result = result <> name;
              ];
           (* append enum state for the number of particles *)
           If[Length[particles] > 0, result = result <> ", ";];
           result = result <> "NUMBER_OF_PARTICLES";
           result = "enum Particles : unsigned {" <>
                    result <> "};\n";
           Return[result];
          ];

CreateParticleNames[particles_List] :=
    Module[{i, par, name, result = ""},
           For[i = 1, i <= Length[particles], i++,
               par = particles[[i]];
               name = CConversion`ToValidCSymbolString[par];
               If[i > 1, result = result <> ", ";];
               result = result <> "\"" <> name <> "\"";
              ];
           result = "const char* particle_names[NUMBER_OF_PARTICLES] = {" <>
                    result <> "};\n";
           Return[result];
          ];

CreateParticleMultiplicity[particles_List] :=
    Module[{i, par, mult, result = ""},
           For[i = 1, i <= Length[particles], i++,
               par = particles[[i]];
               mult = CConversion`ToValidCSymbolString[GetDimension[par]];
               If[i > 1, result = result <> ", ";];
               result = result <> mult;
              ];
           result = "const unsigned particle_multiplicities[NUMBER_OF_PARTICLES] = {" <>
                    result <> "};\n";
           Return[result];
          ];

CreateParticleLaTeXNames[particles_List] :=
    Module[{i, par, latexName, result = ""},
           For[i = 1, i <= Length[particles], i++,
               par = particles[[i]];
               latexName = StringReplace[SARAH`getLaTeXField[par], "\\" -> "\\\\"];
               If[i > 1, result = result <> ", ";];
               result = result <> "\"" <> latexName <> "\"";
              ];
           result = "const char* particle_latex_names[NUMBER_OF_PARTICLES] = {" <>
                    IndentText[result] <> "};\n";
           Return[result];
          ];

FillSpectrumVector[particles_List] :=
    Module[{par, parStr, massStr, latexName, result = ""},
           For[i = 1, i <= Length[particles], i++,
               par = particles[[i]];
               parStr = ToValidCSymbolString[par];
               massStr = ToValidCSymbolString[FlexibleSUSY`M[par]];
               latexName = StringReplace[SARAH`getLaTeXField[par], "\\" -> "\\\\"];
               result = result <> "spectrum.push_back(TParticle(\"" <> parStr <>
                        "\", \"" <> latexName <> "\", to_valarray(PHYSICAL(" <>
                        massStr <> "))));\n";
              ];
           Return[result];
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

CreateFSMassMatrixForUnmixedParticle[TreeMasses`FSMassMatrix[expr_, massESSymbol_, Null]] :=
    Module[{matrix, dim},
           dim = GetDimension[massESSymbol];
           If[dim == 1,
              Print["Warning: trying to create a mass matrix from the 1-plet ", massESSymbol];
             ];
           matrix = Table[expr /. List -> Identity,
                          {SARAH`gt1, 1, dim}, {SARAH`gt2, 1, dim}];
           TreeMasses`FSMassMatrix[matrix, massESSymbol, Null]
          ];

CreateMassCalculationPrototype[m:TreeMasses`FSMassMatrix[expr_, massESSymbol_, Null]] :=
    Module[{result, ev = ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]],
            massMatrix},
           result = "void calculate_" <> ev <> "();\n";
           If[!FreeQ[expr, SARAH`gt1] && !FreeQ[expr, SARAH`gt2],
              massMatrix = CreateFSMassMatrixForUnmixedParticle[m];
              result = CreateMassMatrixGetterPrototype[massMatrix] <>
                       result;
             ];
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

CallMassCalculationFunctions[massMatrices_List] :=
    Module[{result = "", k, sortedMassMatrices, matrix, PredVectorsFirst},
           (* Predicate function which returns false if m2 is a vector
              boson and m1 is not.  True otherwise. *)
           PredVectorsFirst[m1_TreeMasses`FSMassMatrix, m2_TreeMasses`FSMassMatrix] :=
               Module[{es1, es2},
                      es1 = GetMassEigenstate[m1];
                      es2 = GetMassEigenstate[m2];
                      IsVector[es1] || !IsVector[es2]
                     ];
           (* Sort mass matrices such that vector boson masses get
              calculated first.  This is necessary because the later
              calculated masses might depend on some SM mixing angles,
              as ThetaW. *)
           sortedMassMatrices = Sort[massMatrices, PredVectorsFirst];
           For[k = 1, k <= Length[sortedMassMatrices], k++,
               matrix = sortedMassMatrices[[k]];
               result = result <> CallMassCalculationFunction[matrix];
              ];
           Return[result];
          ];

CallMassCalculationFunction[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result = "", k, massESSymbol},
           massESSymbol = GetMassEigenstate[massMatrix];
           result = "calculate_" <> ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]]
                    <> "();\n";
           Return[result];
          ];

IsSymmetric[matrix_List] := IsHermitian[matrix, Identity];

IsHermitian[matrix_List, op_:Susyno`LieGroups`conj] :=
    Module[{rows, cols, i, k, difference},
           rows = Length[matrix];
           For[i = 1, i <= rows, i++,
               cols = Length[matrix[[i]]];
               If[rows =!= cols, Return[False];];
               For[k = 1, k <= i, k++,
                   difference = matrix[[i,k]] - op[matrix[[k,i]]] //. {
                       Susyno`LieGroups`conj[SARAH`sum[ind_,a_,b_,expr_]] :>
                           SARAH`sum[ind,a,b,Susyno`LieGroups`conj[expr]],
                       Susyno`LieGroups`conj[m_[a_,b_]] :>
                           m[b,a] /; MemberQ[SARAH`ListSoftBreakingScalarMasses, m]
                   };
                   If[!PossibleZeroQ[difference], Return[False];];
                  ];
              ];
           Return[True];
          ];

MatrixToCFormString[matrix_List, symbol_String, matrixElementType_:CConversion`realScalarCType] :=
    Module[{dim, result = "", i, k, isSymmetric = IsSymmetric[matrix],
            isHermitian = IsHermitian[matrix], matrixType, dimStr},
           dim = Length[matrix];
           dimStr = ToString[dim];
           matrixType = CreateCType[CConversion`MatrixType[matrixElementType, dim, dim]];
           result = matrixType <> " " <> symbol <> ";\n"; (* not initialized *)
           For[i = 1, i <= dim, i++,
               For[k = 1, k <= dim, k++,
                   result = result <> symbol <> "(" <> ToString[i-1] <>
                            "," <> ToString[k-1] <> ") = ";
                   Which[isSymmetric && i > k,
                         result = result <> symbol <> "(" <> ToString[k-1] <>
                                  "," <> ToString[i-1] <> ");\n"
                         ,
                         isHermitian && i > k,
                         result = result <> "Conj(" <> symbol <> "(" <> ToString[k-1] <>
                                  "," <> ToString[i-1] <> "));\n"
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
    Module[{result, body, ev, matrixSymbol, matrix, massESSymbol,
            inputParsDecl, matrixType, dim, dimStr},
           massESSymbol = GetMassEigenstate[massMatrix];
           ev = ToValidCSymbolString[GetHead[massESSymbol]];
           matrixSymbol = "mass_matrix_" <> ev;
           matrix = GetMassMatrix[massMatrix];
           dim = Length[matrix];
           dimStr = ToString[dim];
           matrixType = CreateCType[CConversion`MatrixType[CConversion`realScalarCType, dim, dim]];
           inputParsDecl = Parameters`CreateLocalConstRefsForInputParameters[matrix, "LOCALINPUT"];
           body = inputParsDecl <> "\n" <> MatrixToCFormString[matrix, matrixSymbol] <> "\n";
           result = matrixType <> " CLASSNAME::get_" <> matrixSymbol <> "() const\n{\n" <>
                    IndentText[body] <>
                    "return " <> matrixSymbol <> ";\n}\n";
           Return[result];
          ];

CreateMassMatrixGetterPrototype[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result, ev, matrixSymbol, matrix, massESSymbol, matrixType,
            dim, dimStr},
           massESSymbol = GetMassEigenstate[massMatrix];
           ev = ToValidCSymbolString[GetHead[massESSymbol]];
           matrix = GetMassMatrix[massMatrix];
           dim = Length[matrix];
           dimStr = ToString[dim];
           matrixType = CreateCType[CConversion`MatrixType[CConversion`realScalarCType, dim, dim]];
           matrixSymbol = "mass_matrix_" <> ev;
           result = matrixType <> " get_" <> matrixSymbol <> "() const;\n";
           Return[result];
          ];

CreateDiagonalizationFunction[matrix_List, eigenVector_, mixingMatrixSymbol_] :=
    Module[{dim, body = "", result, U = "", V = "", dimStr = "", ev, particle, k},
           dim = Length[matrix];
           dimStr = ToString[dim];
           particle = ToValidCSymbolString[GetHead[eigenVector]];
           matrixSymbol = "mass_matrix_" <> particle;
           ev = ToValidCSymbolString[FlexibleSUSY`M[GetHead[eigenVector]]];
           result = "void CLASSNAME::calculate_" <> ev <> "()\n{\n";
           body = "const auto " <> matrixSymbol <> "(get_" <> matrixSymbol <> "());\n";
           If[Head[mixingMatrixSymbol] === List && Length[mixingMatrixSymbol] == 2,
              (* use SVD *)
              U = ToValidCSymbolString[mixingMatrixSymbol[[1]]];
              V = ToValidCSymbolString[mixingMatrixSymbol[[2]]];
              body = body <> "fs_svd(" <>
                     matrixSymbol <> ", " <> ev <> ", " <> U <> ", " <> V <> ");\n";
              ,
              (* use conventional diagonalization *)
              U = ToValidCSymbolString[mixingMatrixSymbol];
              If[IsSymmetric[matrix] && IsFermion[GetHead[eigenVector]],
                 body = body <> "fs_diagonalize_symmetric(" <> matrixSymbol <> ", " <>
                        ev <> ", " <> U <> ");\n";
                 ,
                 body = body <> "fs_diagonalize_hermitian(" <> matrixSymbol <> ", " <>
                        ev <> ", " <> U <> ");\n";
                ];
             ];
           If[IsScalar[eigenVector] || IsVector[eigenVector],
              (* check for tachyons *)
              body = body <> "\n" <>
                     "if (" <> ev <> ".minCoeff() < 0.)\n" <>
                     IndentText["problems.flag_tachyon(" <> particle <> ");"] <> "\n" <>
                     "else\n" <>
                     IndentText["problems.unflag_tachyon(" <> particle <> ");"] <> "\n\n";
              body = body <> ev <> " = AbsSqrt(" <> ev <> ");\n";
             ];
           (* Set the goldstone boson masses equal to the
              corresponding vector boson masses *)
           Return[result <> IndentText[body] <> "}\n"];
          ];

CreateMassCalculationFunction[m:TreeMasses`FSMassMatrix[mass_, massESSymbol_, Null]] :=
    Module[{result, ev = ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]], body,
            inputParsDecl, expr, particle, dim, dimStr, phase, massMatrix},
           result = "void CLASSNAME::calculate_" <> ev <> "()\n{\n";
           (* Remove color SU(3) generators, structure functions and
              Kronecker delta with color indices.
              Note: ct1 ... ct4 are reserved SU(3) color indices of
              the fundamental representation of SU(3) in SARAH.
           *)
           expr = StripGenerators[mass[[1]],
                                  {SARAH`ct1, SARAH`ct2, SARAH`ct3, SARAH`ct4}];
           dim = GetDimension[massESSymbol];
           dimStr = ToString[dim];
           inputParsDecl = Parameters`CreateLocalConstRefsForInputParameters[expr, "LOCALINPUT"];
           If[dim == 1,
              body = inputParsDecl <> "\n" <> ev <> " = " <>
                     RValueToCFormString[expr] <> ";\n";,
              If[FreeQ[expr, SARAH`gt1] && FreeQ[expr, SARAH`gt2],
                 body = inputParsDecl <> "\n" <> ev <>
                        ".setConstant(" <> RValueToCFormString[expr] <> ");\n";,
                 body = inputParsDecl <> "\n" <>
                        "for (int gt1 = 1; gt1 <= " <> dimStr <> "; gt1++) {\n" <>
                        IndentText[ev <> "(gt1) = " <> RValueToCFormString[expr /. SARAH`gt2 -> SARAH`gt1] <> ";"] <>
                        "\n}\n";
                ];
             ];
           phase = Parameters`GetPhase[massESSymbol];
           If[IsFermion[massESSymbol] && phase =!= Null &&
              !IsMassless[massESSymbol],
              particle = ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]];
              body = body <> "\n" <> "if (" <> ev <> " < 0.) {\n" <>
                     IndentText[particle <> " *= -1;\n" <>
                                CConversion`ToValidCSymbolString[phase] <> " = " <>
                                CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]] <>
                                "(0., 1.);"] <> "\n}\n";
             ];
           If[(IsVector[massESSymbol] || IsScalar[massESSymbol]) &&
              !IsMassless[massESSymbol],
              (* check for tachyons *)
              particle = ToValidCSymbolString[massESSymbol];
              If[dim == 1,
                 body = body <> "\n" <> "if (" <> ev <> " < 0.)\n";,
                 body = body <> "\n" <> "if (" <> ev <> ".minCoeff() < 0.)\n";
                ];
              body = body <>
                     IndentText["problems.flag_tachyon(" <> particle <> ");"] <> "\n" <>
                     "else\n" <>
                     IndentText["problems.unflag_tachyon(" <> particle <> ");"] <> "\n\n";
              body = body <> ev <> " = AbsSqrt(" <> ev <> ");\n";
             ];
           body = IndentText[body];
           result = result <> body <> "}\n\n";
           If[!FreeQ[mass, SARAH`gt1] && !FreeQ[mass, SARAH`gt2],
              massMatrix = CreateFSMassMatrixForUnmixedParticle[m];
              result = CreateMassMatrixGetterFunction[massMatrix] <>
                       "\n" <> result;
             ];
           Return[result];
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
    Module[{result = "", massESSymbol, dim, dimStr, returnType},
           massESSymbol = GetMassEigenstate[massMatrix];
           dim = GetDimension[massESSymbol];
           dimStr = ToString[dim];
           If[dim == 1,
              returnType = CConversion`ScalarType[CConversion`realScalarCType];,
              returnType = CConversion`ArrayType[CConversion`realScalarCType, dim];
             ];
           result = CreateCType[returnType] <> " " <>
                    ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]] <> ";\n";
           Return[result];
          ];

CreatePhysicalMassInitialization[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result = "", massESSymbol, dim, matrixType},
           massESSymbol = GetMassEigenstate[massMatrix];
           dim = GetDimension[massESSymbol];
           matrixType = CreateCType[CConversion`ArrayType[CConversion`realScalarCType, dim]];
           result = ", " <> ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]];
           If[dim == 1,
              result = result <> "(0)";,
              result = result <> "(" <> matrixType <> "::Zero())";
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
           matrixType = CreateCType[GetMixingMatrixType[massMatrix]];
           result = DefineMatrix[mixingMatrixSymbol, matrixType];
           Return[result];
          ];

ClearOutputParameters[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result, massESSymbol, mixingMatrixSymbol, matrixType,
            dim, i, massESType},
           massESSymbol = GetMassEigenstate[massMatrix];
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           dim = GetDimension[massESSymbol];
           massESType = CreateCType[CConversion`ArrayType[CConversion`realScalarCType, dim]];
           If[dim == 1,
              result = ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]] <> " = 0.0;\n";
              ,
              result = ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]] <> " = " <> massESType <> "::Zero();\n";
             ];
           If[mixingMatrixSymbol =!= Null,
              matrixType = CreateCType[GetMixingMatrixType[massMatrix]];
              If[Head[mixingMatrixSymbol] === List,
                 For[i = 1, i <= Length[mixingMatrixSymbol], i++,
                     result = result <> ToValidCSymbolString[mixingMatrixSymbol[[i]]] <>
                              " = " <> matrixType <> "::Zero();\n";
                    ];
                 ,
                 result = result <> ToValidCSymbolString[mixingMatrixSymbol] <>
                          " = " <> matrixType <> "::Zero();\n";
                ];
             ];
           Return[result];
          ];

CopyDRBarMassesToPoleMasses[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result, massESSymbol, mixingMatrixSymbol, dim, dimStr,
            i, massStr, mixStr},
           massESSymbol = GetMassEigenstate[massMatrix];
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           dim = GetDimension[massESSymbol];
           dimStr = ToString[dim];
           massStr = ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]];
           (* copy mass *)
           result = "PHYSICAL(" <> massStr <> ") = " <> massStr <> ";\n";
           If[mixingMatrixSymbol =!= Null,
              If[Head[mixingMatrixSymbol] === List,
                 For[i = 1, i <= Length[mixingMatrixSymbol], i++,
                     mixStr = ToValidCSymbolString[mixingMatrixSymbol[[i]]];
                     result = result <> "PHYSICAL(" <> mixStr <> ") = " <> mixStr <> ";\n";
                    ];
                 ,
                 mixStr = ToValidCSymbolString[mixingMatrixSymbol];
                 result = result <> "PHYSICAL(" <> mixStr <> ") = " <> mixStr <> ";\n";
                ];
             ];
           Return[result];
          ];

InitializeMatrix[Null, _] := "";

InitializeMatrix[matrix_Symbol, type_String] :=
    ", " <> ToValidCSymbolString[matrix] <> "(" <> type <> "::Zero())";

InitializeMatrix[matrix_List, type_String] :=
    Module[{result = ""},
           (result = result <> InitializeMatrix[#, type])& /@ matrix;
           Return[result];
          ];

CreateMixingMatrixInitialization[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result, mixingMatrixSymbol, matrixType},
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           matrixType = CreateCType[GetMixingMatrixType[massMatrix]];
           result = InitializeMatrix[mixingMatrixSymbol, matrixType];
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

FindDependenceNums[massMatrices_List] :=
    Module[{hyperchargeCoupling, leftCoupling},
           If[!dependenceNumsUpToDate,
              hyperchargeCoupling = FindHyperchargeGaugeCoupling[];
              leftCoupling = FindLeftGaugeCoupling[];
              (* @todo derive Weinberg angle in terms of fundamental model
                 parameters from SARAH's expressions.  The definition below
                 might not be true in a general model. *)
              dependenceNums = Join[
                  { Rule[SARAH`Weinberg,
                         WeinbergAngle`ExpressWeinbergAngleInTermsOfGaugeCouplings[massMatrices]] },
                  Cases[SARAH`ParameterDefinitions,
                        {parameter_ /; !MemberQ[Parameters`GetModelParameters[], parameter] &&
                         parameter =!= SARAH`Weinberg &&
                         parameter =!= SARAH`electricCharge,
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
                        {parameter_ /; !MemberQ[Parameters`GetModelParameters[], parameter] &&
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

CreateDependenceNumPrototypes[massMatrices_List] :=
    Module[{dependenceNums, result = ""},
           dependenceNums = FindDependenceNums[massMatrices];
           (result = result <> CreateDependenceNumPrototype[#])& /@ dependenceNums;
           Return[result];
          ];

CreateDependenceNumFunction[Rule[parameter_, value_]] :=
    Module[{result, body, parStr},
           parStr = ToValidCSymbolString[parameter];
           body = "return " <> RValueToCFormString[Simplify[value]] <> ";\n";
           result = "double CLASSNAME::" <> parStr <> "() const\n{\n" <>
                    IndentText[body] <> "}\n\n";
           Return[result];
          ];

CreateDependenceNumFunctions[massMatrices_List] :=
    Module[{dependenceNums, result = ""},
           dependenceNums = FindDependenceNums[massMatrices];
           (result = result <> CreateDependenceNumFunction[#])& /@ dependenceNums;
           Return[result];
          ];

ReplaceDependencies[expr_] :=
    expr /. FindDependenceNumRules[];

CallThirdGenerationHelperFunctionName[fermion_, msf1_String, msf2_String, theta_] :=
    "calculate_" <>
    CConversion`ToValidCSymbolString[FlexibleSUSY`M[fermion]] <>
    "_3rd_generation(" <> msf1 <> ", " <> msf2 <> ", " <> theta <> ")";

CreateThirdGenerationHelperPrototype[fermion_] :=
    "void calculate_" <>
    CConversion`ToValidCSymbolString[FlexibleSUSY`M[fermion]] <>
    "_3rd_generation(double&, double&, double&) const;\n";

CreateThirdGenerationHelperFunction[fermion_ /; fermion === SARAH`TopQuark] :=
    "void CLASSNAME::calculate_" <>
    CConversion`ToValidCSymbolString[FlexibleSUSY`M[fermion]] <>
    "_3rd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = " <> CConversion`RValueToCFormString[SARAH`SoftSquark[2,2]] <> ";
   sf_data.mr2 = " <> CConversion`RValueToCFormString[SARAH`SoftUp[2,2]] <> ";
   sf_data.yf  = " <> CConversion`RValueToCFormString[SARAH`UpYukawa[2,2]] <> ";
   sf_data.vd  = " <> CConversion`RValueToCFormString[SARAH`VEVSM1] <> ";
   sf_data.vu  = " <> CConversion`RValueToCFormString[SARAH`VEVSM2] <> ";
   sf_data.gY  = " <> CConversion`RValueToCFormString[SARAH`hyperchargeCoupling /.
                                                      Parameters`ApplyGUTNormalization[]] <> ";
   sf_data.g2  = " <> CConversion`RValueToCFormString[SARAH`leftCoupling] <> ";
   sf_data.Tyf = " <> CConversion`RValueToCFormString[SARAH`TrilinearUp[2,2]] <> ";
   sf_data.mu  = " <> CConversion`RValueToCFormString[Parameters`GetEffectiveMu[]] <> ";
   sf_data.T3  = sfermions::Isospin[sfermions::up];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::up];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::up];

   Eigen::Array<double,2,1> msf;

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf);
   msf1  = msf(0);
   msf2  = msf(1);
}
";

CreateThirdGenerationHelperFunction[fermion_ /; fermion === SARAH`BottomQuark] :=
    "void CLASSNAME::calculate_" <>
    CConversion`ToValidCSymbolString[FlexibleSUSY`M[fermion]] <>
    "_3rd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = " <> CConversion`RValueToCFormString[SARAH`SoftSquark[2,2]] <> ";
   sf_data.mr2 = " <> CConversion`RValueToCFormString[SARAH`SoftDown[2,2]] <> ";
   sf_data.yf  = " <> CConversion`RValueToCFormString[SARAH`DownYukawa[2,2]] <> ";
   sf_data.vd  = " <> CConversion`RValueToCFormString[SARAH`VEVSM1] <> ";
   sf_data.vu  = " <> CConversion`RValueToCFormString[SARAH`VEVSM2] <> ";
   sf_data.gY  = " <> CConversion`RValueToCFormString[SARAH`hyperchargeCoupling /.
                                                      Parameters`ApplyGUTNormalization[]] <> ";
   sf_data.g2  = " <> CConversion`RValueToCFormString[SARAH`leftCoupling] <> ";
   sf_data.Tyf = " <> CConversion`RValueToCFormString[SARAH`TrilinearDown[2,2]] <> ";
   sf_data.mu  = " <> CConversion`RValueToCFormString[Parameters`GetEffectiveMu[]] <> ";
   sf_data.T3  = sfermions::Isospin[sfermions::down];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::down];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::down];

   Eigen::Array<double,2,1> msf;

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf);
   msf1  = msf(0);
   msf2  = msf(1);
}
";

CreateThirdGenerationHelperFunction[fermion_ /; fermion === SARAH`Electron] :=
    "void CLASSNAME::calculate_" <>
    CConversion`ToValidCSymbolString[FlexibleSUSY`M[fermion]] <>
    "_3rd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = " <> CConversion`RValueToCFormString[SARAH`SoftLeftLepton[2,2]] <> ";
   sf_data.mr2 = " <> CConversion`RValueToCFormString[SARAH`SoftRightLepton[2,2]] <> ";
   sf_data.yf  = " <> CConversion`RValueToCFormString[SARAH`ElectronYukawa[2,2]] <> ";
   sf_data.vd  = " <> CConversion`RValueToCFormString[SARAH`VEVSM1] <> ";
   sf_data.vu  = " <> CConversion`RValueToCFormString[SARAH`VEVSM2] <> ";
   sf_data.gY  = " <> CConversion`RValueToCFormString[SARAH`hyperchargeCoupling /.
                                                      Parameters`ApplyGUTNormalization[]] <> ";
   sf_data.g2  = " <> CConversion`RValueToCFormString[SARAH`leftCoupling] <> ";
   sf_data.Tyf = " <> CConversion`RValueToCFormString[SARAH`TrilinearLepton[2,2]] <> ";
   sf_data.mu  = " <> CConversion`RValueToCFormString[Parameters`GetEffectiveMu[]] <> ";
   sf_data.T3  = sfermions::Isospin[sfermions::electron];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::electron];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::electron];

   Eigen::Array<double,2,1> msf;

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf);
   msf1  = msf(0);
   msf2  = msf(1);
}
";

CreateThirdGenerationHelperFunction[fermion_ /; fermion === SARAH`Neutrino] :=
    "void CLASSNAME::calculate_" <>
    CConversion`ToValidCSymbolString[FlexibleSUSY`M[fermion]] <>
    "_3rd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = " <> CConversion`RValueToCFormString[SARAH`SoftLeftLepton[2,2]] <> ";
   sf_data.mr2 = 0.;
   sf_data.yf  = 0.;
   sf_data.vd  = " <> CConversion`RValueToCFormString[SARAH`VEVSM1] <> ";
   sf_data.vu  = " <> CConversion`RValueToCFormString[SARAH`VEVSM2] <> ";
   sf_data.gY  = " <> CConversion`RValueToCFormString[SARAH`hyperchargeCoupling /.
                                                      Parameters`ApplyGUTNormalization[]] <> ";
   sf_data.g2  = " <> CConversion`RValueToCFormString[SARAH`leftCoupling] <> ";
   sf_data.Tyf = 0.;
   sf_data.mu  = " <> CConversion`RValueToCFormString[Parameters`GetEffectiveMu[]] <> ";
   sf_data.T3  = sfermions::Isospin[sfermions::neutrino];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::neutrino];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::neutrino];

   Eigen::Array<double,2,1> msf;

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf);
   msf1  = msf(0);
   msf2  = msf(1);
}
";

CreateThirdGenerationHelperFunction[fermion_] :=
    Module[{},
           Print["Error: ", fermion, " does not seem to be a",
                 " 3rd generation SM fermion."];
           "void CLASSNAME::calculate_" <>
           CConversion`ToValidCSymbolString[FlexibleSUSY`M[fermion]] <>
           "_3rd_generation(double& msf1, double& msf2, double& theta) const {}"
          ];

CreateThirdGenerationHelpers[] :=
    Module[{prototypes, functions},
           functions = CreateThirdGenerationHelperFunction[SARAH`TopQuark] <> "\n" <>
                       CreateThirdGenerationHelperFunction[SARAH`BottomQuark] <> "\n" <>
                       CreateThirdGenerationHelperFunction[SARAH`Neutrino] <> "\n" <>
                       CreateThirdGenerationHelperFunction[SARAH`Electron];
           prototypes = CreateThirdGenerationHelperPrototype[SARAH`TopQuark] <>
                        CreateThirdGenerationHelperPrototype[SARAH`BottomQuark] <>
                        CreateThirdGenerationHelperPrototype[SARAH`Neutrino] <>
                        CreateThirdGenerationHelperPrototype[SARAH`Electron];
           {prototypes, functions}
          ];

End[];

EndPackage[];

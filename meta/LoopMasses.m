
BeginPackage["LoopMasses`", {"SARAH`", "TextFormatting`",
                             "CConversion`", "TreeMasses`",
                             "SelfEnergies`", "TwoLoop`", "Parameters`"}];

CreateOneLoopPoleMassFunctions::usage="";
CreateOneLoopPoleMassPrototypes::usage="";
CallAllPoleMassFunctions::usage="";
CreateRunningDRbarMassPrototypes::usage="";
CreateRunningDRbarMassFunctions::usage="";

GetLoopCorrectedParticles::usage="Returns list of all particles that
get loop corrected masses.  These are all particles, except for
ghosts.";

CreateLSPFunctions::usage="";

Begin["`Private`"];

GetLoopCorrectedParticles[states_] :=
    Module[{particles},
           particles = GetParticles[states];
           Select[particles, (!IsGhost[#])&]
          ];

FillTadpoleMatrix[{}, _] := "";

FillTadpoleMatrix[tadpoles_List, matrixName_:"tadpoles"] :=
    Module[{result, dim, dimStr, particle, i, particleIndex, vev,
            tadpoleMatrixType},
           particle = tadpoles[[1,1]];
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           tadpoleMatrixType = CreateCType[CConversion`MatrixType[CConversion`realScalarCType, dim, dim]];
           If[dim > 1,
              result = tadpoleMatrixType <> " " <> matrixName <> "(" <> dimStr <> "," <> dimStr <> ");\n";
              For[i = 1, i <= Length[tadpoles], i++,
                  particleIndex = ToString[tadpoles[[i,2]]];
                  vev           = ToValidCSymbolString[tadpoles[[i,3]]];
                  result = result <>
                           matrixName <> "(" <> particleIndex <> "," <> particleIndex <> ") = Re(" <>
                           CreateTadpoleFunctionName[particle] <> "(" <> particleIndex <> ")) / " <>
                           vev <> ";\n";
                 ];
              ,
              particleIndex = ToString[tadpoles[[1,2]]];
              vev           = ToValidCSymbolString[tadpoles[[1,3]]];
              result = "const double " <> matrixName <> " = " <>
                       CreateTadpoleFunctionName[particle] <> "(" <> particleIndex <> ");\n";
             ];
           Return[result];
          ];

Do1DimScalar[particleName_String, massName_String, selfEnergyFunction_String,
             momentum_String, tadpole_String:""] :=
    "const double p = " <> momentum <> ";\n" <>
    "const double self_energy = Re(" <> selfEnergyFunction <> "(p));\n" <>
    "const double mass_sqr = Sqr(" <> massName <> ") - self_energy" <>
    If[tadpole == "", "", " + " <> tadpole] <> ";\n\n" <>
    "if (mass_sqr < 0.)\n" <>
    IndentText["problems.flag_tachyon(" <> particleName <> ");"] <> "\n\n" <>
    "PHYSICAL(" <> massName <> ") = ZeroSqrt(mass_sqr);\n";

Do1DimFermion[massName_String, selfEnergyFunctionS_String,
              selfEnergyFunctionPL_String, selfEnergyFunctionPR_String, momentum_String] :=
    "const double p = " <> momentum <> ";\n" <>
    "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p));\n" <>
    "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p));\n" <>
    "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p));\n" <>
    "PHYSICAL(" <> massName <> ") = " <> massName <>
    " - self_energy_1 - " <> massName <> " * (self_energy_PL + self_energy_PR);\n";

Do1DimVector[particleName_String, massName_String, selfEnergyFunction_String,
             momentum_String] :=
    "const double p = " <> momentum <> ";\n" <>
    "const double self_energy = Re(" <> selfEnergyFunction <> "(p));\n" <>
    "const double mass_sqr = Sqr(" <> massName <> ") - self_energy;\n\n" <>
    "if (mass_sqr < 0.)\n" <>
    IndentText["problems.flag_tachyon(" <> particleName <> ");"] <> "\n\n" <>
    "PHYSICAL(" <> massName <> ") = ZeroSqrt(mass_sqr);\n";


(* ********** fast diagonalization routines ********** *)

DoFastDiagonalization[particle_Symbol /; IsScalar[particle], tadpoles_List] :=
    Module[{result, dim, dimStr, massName, particleName, mixingMatrix, selfEnergyFunction,
            tadpoleMatrix, U, V, massMatrixStr, selfEnergyIsSymmetric, selfEnergyMatrixType},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           mixingMatrix = FindMixingMatrixSymbolFor[particle];
           massMatrixStr = "get_mass_matrix_" <> ToValidCSymbolString[particle];
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle];
           selfEnergyMatrixType = CreateCType[CConversion`MatrixType[CConversion`realScalarCType, dim, dim]];
           tadpoleMatrix = FillTadpoleMatrix[tadpoles, "tadpoles"];
           If[dim > 1,
              selfEnergyIsSymmetric = Length[Flatten[{mixingMatrix}]] === 1;
              result = tadpoleMatrix <>
                       selfEnergyMatrixType <> " self_energy;\n" <>
                       "for (unsigned i1 = 0; i1 < " <> dimStr <>"; ++i1) {\n" <>
                       IndentText["for (unsigned i2 = " <> If[selfEnergyIsSymmetric,"i1","0"] <>
                                  "; i2 < " <> dimStr <>"; ++i2) {\n" <>
                                  IndentText["const double p = AbsSqrt(" <> massName <> "(i1) * " <> 
                                             massName <> "(i2));\n" <>
                                             "self_energy(i1,i2) = Re(" <>
                                             selfEnergyFunction <> "(p,i1,i2));\n"] <>
                                  "}\n"
                                 ] <>
                       "}\n" <>
                       If[selfEnergyIsSymmetric, "Symmetrize(self_energy);\n", ""] <>
                       "const " <> selfEnergyMatrixType <>" M_1loop(" <> massMatrixStr <>
                       "() - self_energy" <>
                       If[tadpoleMatrix == "", "", " + tadpoles"] <> ");\n";
              If[Head[mixingMatrix] === List,
                 U = ToValidCSymbolString[mixingMatrix[[1]]];
                 V = ToValidCSymbolString[mixingMatrix[[2]]];
                 result = result <>
                          "fs_svd(M_1loop, PHYSICAL(" <> massName <> "), " <>
                          "PHYSICAL(" <> U <> "), " <>
                          "PHYSICAL(" <> V <> "));\n";
                 ,
                 U = ToValidCSymbolString[mixingMatrix];
                 result = result <>
                       "fs_diagonalize_hermitian(M_1loop, PHYSICAL(" <> massName <> "), " <>
                       "PHYSICAL(" <> U <> "));\n";
                ];
              result = result <>
                       "\n" <>
                       "if (" <> massName <> ".minCoeff() < 0.)\n" <>
                       IndentText["problems.flag_tachyon(" <> particleName <> ");"] <> "\n\n" <>
                       "PHYSICAL(" <> massName <> ") = ZeroSqrt(PHYSICAL(" <>
                       massName <> "));\n";
              ,
              result = Do1DimScalar[particleName, massName, selfEnergyFunction, massName,
                                    If[tadpoleMatrix == "", "", "tadpoles"]];
             ];
           Return[result];
          ];

DoFastDiagonalization[particle_Symbol /; IsFermion[particle], _] :=
    Module[{result, dim, dimStr, massName, mixingMatrix, U, V,
            selfEnergyFunctionS, selfEnergyFunctionPL, selfEnergyFunctionPR,
            massMatrixStr, selfEnergyMatrixType},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           massMatrixStr = "get_mass_matrix_" <> ToValidCSymbolString[particle];
           If[IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0,
              If[dim == 1,
                 Return["PHYSICAL(" <> massName <> ") = 0.;\n"];,
                 Return["PHYSICAL(" <> massName <> ").setConstant(0.);\n"];
                ];
             ];
           mixingMatrix = FindMixingMatrixSymbolFor[particle];
           selfEnergyMatrixType = CreateCType[CConversion`MatrixType[CConversion`realScalarCType, dim, dim]];
           selfEnergyFunctionS  = SelfEnergies`CreateSelfEnergyFunctionName[particle[1]];
           selfEnergyFunctionPL = SelfEnergies`CreateSelfEnergyFunctionName[particle[PL]];
           selfEnergyFunctionPR = SelfEnergies`CreateSelfEnergyFunctionName[particle[PR]];
           If[dim > 1,
              result = selfEnergyMatrixType <> " self_energy_1;\n" <>
                       selfEnergyMatrixType <> " self_energy_PL;\n" <>
                       selfEnergyMatrixType <> " self_energy_PR;\n" <>
                       "for (unsigned i1 = 0; i1 < " <> dimStr <>"; ++i1) {\n" <>
                       IndentText["for (unsigned i2 = 0; i2 < " <> dimStr <>"; ++i2) {\n" <>
                                  IndentText["const double p = AbsSqrt(" <> massName <> "(i1) * " <> 
                                             massName <> "(i2));\n" <>
                                             "self_energy_1(i1,i2)  = Re(" <>
                                             selfEnergyFunctionS <> "(p,i1,i2));\n" <>
                                             "self_energy_PL(i1,i2) = Re(" <>
                                             selfEnergyFunctionPL <> "(p,i1,i2));\n" <>
                                             "self_energy_PR(i1,i2) = Re(" <>
                                             selfEnergyFunctionPR <> "(p,i1,i2));\n"
                                            ] <>
                                  "}\n"
                                 ] <>
                       "}\n" <>
                       "const " <> selfEnergyMatrixType <> " M_tree(" <> massMatrixStr <> "());\n" <>
                       "const " <> selfEnergyMatrixType <> " delta_M(- self_energy_PR * M_tree " <>
                       "- M_tree * self_energy_PL - self_energy_1);\n";
              If[IsMajoranaFermion[particle],
                 result = result <>
                          "const " <> selfEnergyMatrixType <>
                          " M_1loop(M_tree + 0.5 * (delta_M + delta_M.transpose()));\n";
                 ,
                 result = result <>
                          "const " <> selfEnergyMatrixType <>
                          " M_1loop(M_tree + delta_M);\n";
                ];
              If[Head[mixingMatrix] === List,
                 (* two mixing matrixs => SVD *)
                 U = ToValidCSymbolString[mixingMatrix[[1]]];
                 V = ToValidCSymbolString[mixingMatrix[[2]]];
                 result = result <>
                          "fs_svd(M_1loop, " <>
                          "PHYSICAL(" <> massName <> "), " <>
                          "PHYSICAL(" <> U <> "), " <>
                          "PHYSICAL(" <> V <> "));\n";
                 ,
                 U = ToValidCSymbolString[mixingMatrix];
                 result = result <>
                          "fs_diagonalize_symmetric(M_1loop, " <>
                          "PHYSICAL(" <> massName <> "), " <>
                          "PHYSICAL(" <> U <> "));\n";
                ];
              ,
              (* for a dimension 1 fermion it plays not role if it's a
                 Majorana ferimion or not *)
              result = Do1DimFermion[massName, selfEnergyFunctionS,
                                     selfEnergyFunctionPL, selfEnergyFunctionPR,
                                     massName];
             ];
           Return[result];
          ];

DoFastDiagonalization[particle_Symbol /; IsVector[particle], _] :=
    Module[{result, dim, dimStr, massName, particleName, mixingMatrix, selfEnergyFunction},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           mixingMatrix = ToValidCSymbolString[FindMixingMatrixSymbolFor[particle]];
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle];
           If[IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0,
              If[dim == 1,
                 Return["PHYSICAL(" <> massName <> ") = 0.;\n"];,
                 Return["PHYSICAL(" <> massName <> ").setConstant(0.);\n"];
                ];
             ];
           If[dim > 1,
              result = "WARNING(\"diagonalization of " <> ToString[particle] <> " not implemented\");\n";
              ,
              result = Do1DimVector[particleName, massName, selfEnergyFunction, massName];
             ];
           Return[result];
          ];

DoFastDiagonalization[particle_Symbol, _] :=
    "ERROR(\"fast diagonalization of " <> ToString[particle] <> " not implemented\");\n";

(* ********** medium diagonalization routines ********** *)

DoMediumDiagonalization[particle_Symbol /; IsScalar[particle], inputMomentum_, tadpole_List] :=
    Module[{result, dim, dimStr, massName, particleName, mixingMatrix, selfEnergyFunction,
            momentum = inputMomentum, U, V, Utemp, Vtemp, tadpoleMatrix, diagSnippet,
            massMatrixStr, selfEnergyIsSymmetric,
            selfEnergyMatrixType, eigenArrayType,
            addTwoLoopHiggsContributions = "", calcTwoLoopHiggsContributions = "",
            numberOfIndependentMatrixEntries, numberOfIndependentMatrixEntriesStr, n, l, k},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[inputMomentum == "", momentum = massName];
           mixingMatrix = FindMixingMatrixSymbolFor[particle];
           mixingMatrixType = CreateCType[CConversion`MatrixType[CConversion`realScalarCType, dim, dim]];
           selfEnergyMatrixType = mixingMatrixType;
           eigenArrayType = CreateCType[CConversion`ArrayType[CConversion`realScalarCType, dim]];
           (* create diagonalisation code snippet *)
           If[Head[mixingMatrix] === List,
              U = ToValidCSymbolString[mixingMatrix[[1]]];
              V = ToValidCSymbolString[mixingMatrix[[2]]];
              Utemp = "mix_" <> U;
              Vtemp = "mix_" <> V;
              diagSnippet = mixingMatrixType <> " " <> Utemp <> ", " <> Vtemp <> ";\n" <>
                            "fs_svd(M_1loop, eigen_values, " <> Utemp <> ", " <> Vtemp <> ");\n\n" <>
                            "if (eigen_values(es) < 0.)\n" <>
                            IndentText["problems.flag_tachyon(" <> particleName <> ");"] <> "\n\n" <>
                            "PHYSICAL(" <> massName <> "(es)) = ZeroSqrt(eigen_values(es));\n" <>
                            "if (es == 0) {\n" <>
                            IndentText["PHYSICAL(" <> U <> ") = " <> Utemp <> ";\n" <>
                                       "PHYSICAL(" <> V <> ") = " <> Vtemp <> ";\n"] <>
                            "}\n";
              ,
              U = ToValidCSymbolString[mixingMatrix];
              Utemp = "mix_" <> U;
              diagSnippet = mixingMatrixType <> " " <> Utemp <> ";\n" <>
                            "fs_diagonalize_hermitian(M_1loop, eigen_values, " <> Utemp <> ");\n\n" <>
                            "if (eigen_values(es) < 0.)\n" <>
                            IndentText["problems.flag_tachyon(" <> particleName <> ");"] <> "\n\n" <>
                            "PHYSICAL(" <> massName <> "(es)) = ZeroSqrt(eigen_values(es));\n";
              If[mixingMatrix =!= Null,
                 diagSnippet = diagSnippet <>
                               "if (es == 0)\n" <>
                               IndentText["PHYSICAL(" <> U <> ") = " <> Utemp <> ";\n"];
                ];
             ];
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle];
           tadpoleMatrix = FillTadpoleMatrix[tadpole, "tadpoles"];
           (* fill self-energy and do diagonalisation *)
           If[dim > 1,
              If[SARAH`UseHiggs2LoopMSSM === True ||
                 FlexibleSUSY`UseHiggs2LoopNMSSM === True,
                 If[MemberQ[{SARAH`HiggsBoson, SARAH`PseudoScalar}, particle],
                    numberOfIndependentMatrixEntries = Parameters`NumberOfIndependentEntriesOfSymmetricMatrix[dim];
                    numberOfIndependentMatrixEntriesStr = ToString[numberOfIndependentMatrixEntries];
                    addTwoLoopHiggsContributions = "";
                    For[k = 0; n = 0, k < dim, k++,
                        For[l = k, l < dim, l++; n++,
                            addTwoLoopHiggsContributions = addTwoLoopHiggsContributions <>
                               "self_energy(" <> ToString[k] <> ", " <>
                               ToString[l] <> ") += two_loop[" <> ToString[n] <> "];\n";
                           ];
                       ];
                    addTwoLoopHiggsContributions = "
if (pole_mass_loop_order > 1) {
" <> IndentText[addTwoLoopHiggsContributions] <> "\
}
";
                    calcTwoLoopHiggsContributions = "
// two-loop Higgs self-energy contributions
double two_loop[" <> numberOfIndependentMatrixEntriesStr <> "] = { 0. };
if (pole_mass_loop_order > 1)
" <> IndentText["\
self_energy_" <> CConversion`ToValidCSymbolString[particle] <> "_2loop(two_loop);
"] <> "\
";
                   ];
                ];
              selfEnergyIsSymmetric = Length[Flatten[{mixingMatrix}]] === 1;
              massMatrixStr = "get_mass_matrix_" <> ToValidCSymbolString[particle];
              result = tadpoleMatrix <>
                       selfEnergyMatrixType <> " self_energy;\n" <>
                       "const " <> selfEnergyMatrixType <> " M_tree(" <> massMatrixStr <> "());\n" <>
                       calcTwoLoopHiggsContributions <> "\n" <>
                       "for (unsigned es = 0; es < " <> dimStr <> "; ++es) {\n" <>
                       IndentText["const double p = Abs(" <> momentum <> "(es));\n" <>
                                  "for (unsigned i1 = 0; i1 < " <> dimStr <> "; ++i1) {\n" <>
                                  IndentText["for (unsigned i2 = " <> If[selfEnergyIsSymmetric,"i1","0"] <>
                                             "; i2 < " <> dimStr <> "; ++i2) {\n" <>
                                             IndentText["self_energy(i1,i2) = Re(" <>
                                                        selfEnergyFunction <> "(p,i1,i2));\n"
                                                       ] <>
                                             "}\n"
                                            ] <>
                                  "}\n" <>
                                  addTwoLoopHiggsContributions <> "\n" <>
                                  If[selfEnergyIsSymmetric, "Symmetrize(self_energy);\n", ""] <>
                                  "const " <> selfEnergyMatrixType <> " M_1loop(M_tree - self_energy" <>
                                  If[tadpoleMatrix == "", "", " + tadpoles"] <> ");\n" <>
                                  eigenArrayType <> " eigen_values;\n" <>
                                  diagSnippet
                                 ] <>
                       "}\n";
              ,
              result = tadpoleMatrix <>
                       Do1DimScalar[particleName, massName, selfEnergyFunction, momentum,
                                    If[tadpoleMatrix == "", "", "tadpoles"]];
             ];
           Return[result];
          ];

DoMediumDiagonalization[particle_Symbol /; IsFermion[particle], inputMomentum_, _] :=
    Module[{result, dim, dimStr, massName, mixingMatrix, U, V,
            selfEnergyFunctionS, selfEnergyFunctionPL, selfEnergyFunctionPR,
            momentum = inputMomentum, massMatrixStr, selfEnergyMatrixType,
            eigenArrayType, mixingMatrixType},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[inputMomentum == "", momentum = massName];
           If[IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0,
              If[dim == 1,
                 Return["PHYSICAL(" <> massName <> ") = 0.;\n"];,
                 Return["PHYSICAL(" <> massName <> ").setConstant(0.);\n"];
                ];
             ];
           mixingMatrix = FindMixingMatrixSymbolFor[particle];
           selfEnergyMatrixType = CreateCType[CConversion`MatrixType[CConversion`realScalarCType, dim, dim]];
           eigenArrayType = CreateCType[CConversion`ArrayType[CConversion`realScalarCType, dim]];
           selfEnergyFunctionS  = SelfEnergies`CreateSelfEnergyFunctionName[particle[1]];
           selfEnergyFunctionPL = SelfEnergies`CreateSelfEnergyFunctionName[particle[PL]];
           selfEnergyFunctionPR = SelfEnergies`CreateSelfEnergyFunctionName[particle[PR]];
           If[dim > 1,
              massMatrixStr = "get_mass_matrix_" <> ToValidCSymbolString[particle];
              result = selfEnergyMatrixType <> " self_energy_1;\n" <>
                       selfEnergyMatrixType <> " self_energy_PL;\n" <>
                       selfEnergyMatrixType <> " self_energy_PR;\n" <>
                       "const " <> selfEnergyMatrixType <> " M_tree(" <> massMatrixStr <> "());\n" <>
                       "for (unsigned es = 0; es < " <> dimStr <> "; ++es) {\n" <>
                       IndentText["const double p = Abs(" <> momentum <> "(es));\n" <>
                                  "for (unsigned i1 = 0; i1 < " <> dimStr <>"; ++i1) {\n" <>
                                  IndentText["for (unsigned i2 = 0; i2 < " <> dimStr <>"; ++i2) {\n" <>
                                             IndentText["self_energy_1(i1,i2)  = Re(" <>
                                                        selfEnergyFunctionS <> "(p,i1,i2));\n" <>
                                                        "self_energy_PL(i1,i2) = Re(" <>
                                                        selfEnergyFunctionPL <> "(p,i1,i2));\n" <>
                                                        "self_energy_PR(i1,i2) = Re(" <>
                                                        selfEnergyFunctionPR <> "(p,i1,i2));\n"
                                                       ] <>
                                             "}\n"
                                            ] <>
                                  "}\n" <>
                                  "const " <> selfEnergyMatrixType <> " delta_M(- self_energy_PR * M_tree " <>
                                  "- M_tree * self_energy_PL - self_energy_1);\n"
                                 ];
              If[IsMajoranaFermion[particle],
                 result = result <>
                          IndentText["const " <> selfEnergyMatrixType <> " M_1loop(M_tree + 0.5 * (delta_M + delta_M.transpose()));\n"];
                 ,
                 result = result <>
                          IndentText["const " <> selfEnergyMatrixType <> " M_1loop(M_tree + delta_M);\n"];
                ];
              result = result <>
                       IndentText[eigenArrayType <> " eigen_values;\n"];
              If[Head[mixingMatrix] === List,
                 (* two mixing matrixs => SVD *)
                 U = ToValidCSymbolString[mixingMatrix[[1]]];
                 V = ToValidCSymbolString[mixingMatrix[[2]]];
                 result = result <>
                          IndentText["decltype(" <> U <> ") mix_" <> U <> ";\n" <>
                                     "decltype(" <> V <> ") mix_" <> V <> ";\n"];
                 result = result <>
                          IndentText["fs_svd(M_1loop, eigen_values, " <>
                                     "mix_" <> U <> ", " <> "mix_" <> V <> ");\n"];
                 result = result <>
                          IndentText["if (es == 0) {\n" <>
                                     IndentText["PHYSICAL(" <> U <> ") = mix_" <> U <> ";\n" <>
                                                "PHYSICAL(" <> V <> ") = mix_" <> V <> ";\n"] <>
                                     "}\n"
                                    ];
                 ,
                 U = ToValidCSymbolString[mixingMatrix];
                 If[mixingMatrix =!= Null,
                    result = result <>
                             IndentText["decltype(" <> U <> ") mix_" <> U <> ";\n" <>
                                        "fs_diagonalize_symmetric(M_1loop, eigen_values, mix_" <> U <> ");\n" <>
                                        "if (es == 0)\n" <>
                                        IndentText["PHYSICAL(" <> U <> ") = mix_" <> U <> ";\n"]
                                       ];
                    ,
                    mixingMatrixType = CreateCType[CConversion`MatrixType[CConversion`complexScalarCType, dim, dim]];
                    result = result <>
                             IndentText[mixingMatrixType <> " mix_" <> U <> ";\n" <>
                                        "fs_diagonalize_symmetric(M_1loop, eigen_values, mix_" <> U <> ");\n"];
                   ];
                ];
              result = result <>
                       IndentText["PHYSICAL(" <> massName <>
                                  "(es)) = Abs(eigen_values(es));\n"];
              result = result <> "}\n";
              ,
              (* for a dimension 1 fermion it plays not role if it's a
                 Majorana fermion or not *)
              result = Do1DimFermion[massName, selfEnergyFunctionS,
                                     selfEnergyFunctionPL, selfEnergyFunctionPR,
                                     momentum];
             ];
           Return[result];
          ];

DoMediumDiagonalization[particle_Symbol /; IsVector[particle], inputMomentum_, _] :=
    Module[{result, dim, dimStr, massName, particleName, mixingMatrix, selfEnergyFunction,
            momentum = inputMomentum},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[inputMomentum == "", momentum = massName];
           mixingMatrix = ToValidCSymbolString[FindMixingMatrixSymbolFor[particle]];
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle];
           If[IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0,
              If[dim == 1,
                 Return["PHYSICAL(" <> massName <> ") = 0.;\n"];,
                 Return["PHYSICAL(" <> massName <> ").setConstant(0.);\n"];
                ];
             ];
           If[dim > 1,
              result = "WARNING(\"diagonalization of " <> ToString[particle] <> " not implemented\");\n";
              ,
              result = Do1DimVector[particleName, massName, selfEnergyFunction, momentum];
             ];
           Return[result];
          ];

DoMediumDiagonalization[particle_Symbol, inputMomentum_:"", _] :=
    "ERROR(\"medium diagonalization of " <> ToString[particle] <> " not implemented\");\n";

(* ********** high precision diagonalization routines ********** *)

DoSlowDiagonalization[particle_Symbol, tadpole_] :=
    Module[{result, dim, dimStr, massName, inputMomenta, outputMomenta, body},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           inputMomenta = "old_" <> massName;
           outputMomenta = "new_" <> massName;
           body = DoMediumDiagonalization[particle, inputMomenta, tadpole] <> "\n" <>
                  outputMomenta <> " = PHYSICAL(" <> massName <> ");\n" <>
                  "diff = MaxRelDiff(" <> outputMomenta <> ", " <> inputMomenta <> ");\n" <>
                  inputMomenta <> " = " <> outputMomenta <> ";\n" <>
                  "iteration++;\n";
           result = "unsigned iteration = 0;\n" <>
                    "double diff = 0.0;\n" <>
                    "decltype(" <> massName <> ") " <>
                    inputMomenta  <> "(" <> massName <> "), " <>
                    outputMomenta <> "(" <> massName <> ");\n\n" <>
                    "do {\n" <>
                    IndentText[body] <>
                    "} while (diff > precision\n" <>
                    "         && iteration < number_of_mass_iterations);\n";
           Return[result];
          ];

DoDiagonalization[particle_Symbol, FlexibleSUSY`LowPrecision, tadpole_] :=
    "// diagonalization with low precision\n" <> DoFastDiagonalization[particle, tadpole];

DoDiagonalization[particle_Symbol, FlexibleSUSY`MediumPrecision, tadpole_] :=
    "// diagonalization with medium precision\n" <> DoMediumDiagonalization[particle, "", tadpole];

DoDiagonalization[particle_Symbol, FlexibleSUSY`HighPrecision, tadpole_] :=
    "// diagonalization with high precision\n" <> DoSlowDiagonalization[particle, tadpole];

CreateLoopMassFunctionName[particle_Symbol] :=
    "calculate_" <> ToValidCSymbolString[FlexibleSUSY`M[particle]] <> "_pole";

CallPoleMassFunction[particle_Symbol] :=
    CreateLoopMassFunctionName[particle] <> "();\n";

CreateLoopMassPrototype[particle_Symbol] :=
    "void " <> CreateLoopMassFunctionName[particle] <> "();\n";

CreateLoopMassFunction[particle_Symbol, precision_Symbol, tadpole_] :=
    Module[{result, body = ""},
           If[!IsFermion[particle] &&
              !(IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0),
              body = "if (problems.is_tachyon(" <> ToValidCSymbolString[particle] <> "))\n" <>
                     IndentText["return;"] <> "\n";
             ];
           body = body <> DoDiagonalization[particle, precision, tadpole];
           result = "void CLASSNAME::" <> CreateLoopMassFunctionName[particle] <>
                    "()\n{\n" <> IndentText[body] <> "}\n\n";
           Return[result];
          ];

CreateOneLoopPoleMassFunctions[precision_List, oneLoopTadpoles_List, vevs_List] :=
    Module[{result = "", particle, prec, i = 1, f, d, tadpole, fieldsAndVevs = {}},
           (* create list that associates fields at vevs *)
           For[f = 1, f <= Length[oneLoopTadpoles], f++,
               particle = GetField[oneLoopTadpoles[[f]]];
               For[d = 1, d <= GetDimension[particle], d++; i++,
                   AppendTo[fieldsAndVevs, {particle, d, vevs[[i]]}];
                  ];
              ];
           For[i = 1, i <= Length[precision], i++,
               particle = precision[[i,1]];
               prec     = precision[[i,2]];
               tadpole  = Cases[fieldsAndVevs, {particle[_], __}];
               result   = result <> CreateLoopMassFunction[particle, prec, tadpole];
              ];
           Return[result];
          ];

CreateOneLoopPoleMassPrototypes[states_:FlexibleSUSY`FSEigenstates] :=
    Module[{particles, result = ""},
           particles = GetLoopCorrectedParticles[states];
           (result = result <> CreateLoopMassPrototype[#])& /@ particles;
           Return[result];
          ];

CallThreadedPoleMassFunction[particle_Symbol] :=
    Module[{massStr},
           massStr = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           "std::thread thread_" <> massStr <> "(Thread(this, &CLASSNAME::" <>
           CreateLoopMassFunctionName[particle] <> "));\n"
          ];

JoinLoopMassFunctionThread[particle_Symbol] :=
    "thread_" <> ToValidCSymbolString[FlexibleSUSY`M[particle]] <> ".join();\n";

CallAllPoleMassFunctions[states_, enablePoleMassThreads_] :=
    Module[{particles, susyParticles, smParticles, callSusy = "",
            callSM = "", result, joinSmThreads = "", joinSusyThreads = ""},
           particles = GetLoopCorrectedParticles[states];
           susyParticles = Select[particles, (!SARAH`SMQ[#])&];
           smParticles = Complement[particles, susyParticles];
           If[enablePoleMassThreads =!= True,
              (callSusy = callSusy <> CallPoleMassFunction[#])& /@ susyParticles;
              (callSM   = callSM   <> CallPoleMassFunction[#])& /@ smParticles;
              result = callSusy <> "\n" <>
                       "if (calculate_sm_pole_masses) {\n" <>
                       IndentText[callSM] <>
                       "}\n";
              ,
              (callSusy = callSusy <> CallThreadedPoleMassFunction[#])& /@ susyParticles;
              (callSM   = callSM   <> CallThreadedPoleMassFunction[#])& /@ smParticles;
              (joinSmThreads   = joinSmThreads   <> JoinLoopMassFunctionThread[#])& /@ smParticles;
              (joinSusyThreads = joinSusyThreads <> JoinLoopMassFunctionThread[#])& /@ susyParticles;
              result = callSusy <> "\n" <>
                       "if (calculate_sm_pole_masses) {\n" <>
                       IndentText[callSM] <>
                       IndentText[joinSmThreads] <>
                       "}\n\n" <>
                       joinSusyThreads <> "\n";
             ];
           Return[result];
          ];

GetRunningOneLoopDRbarParticles[] :=
    {SARAH`TopQuark, SARAH`BottomQuark, SARAH`Electron, SARAH`Neutrino,
     SARAH`VectorP, SARAH`VectorZ, SARAH`VectorW};

CreateRunningDRbarMassPrototype[particle_ /; IsFermion[particle]] :=
    "double calculate_" <> ToValidCSymbolString[FlexibleSUSY`M[particle]] <>
    "_DRbar(double, int) const;\n";

CreateRunningDRbarMassPrototype[particle_] :=
    "double calculate_" <> ToValidCSymbolString[FlexibleSUSY`M[particle]] <>
    "_DRbar(double) const;\n";

CreateRunningDRbarMassPrototypes[] :=
    Module[{result = "", particles},
           particles = GetRunningOneLoopDRbarParticles[];
           (result = result <> CreateRunningDRbarMassPrototype[#])& /@ particles;
           Return[result];
          ];

CreateRunningDRbarMassFunction[particle_ /; particle === SARAH`BottomQuark] :=
    Module[{result, body, selfEnergyFunctionS, selfEnergyFunctionPL,
            selfEnergyFunctionPR, name, alphaS, drbarConversion, gPrime,
            dimParticle, thirdGenMass},
           dimParticle = TreeMasses`GetDimension[particle];
           thirdGenMass = TreeMasses`GetThirdGenerationMass[particle];
           selfEnergyFunctionS  = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[1]];
           selfEnergyFunctionPL = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[PL]];
           selfEnergyFunctionPR = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[PR]];
           name = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[IsMassless[particle],
              result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double, int) const\n{\n";
              body = "return 0.0;\n";
              ,
              alphaS = SARAH`strongCoupling^2/(4 Pi);
              gPrime = SARAH`hyperchargeCoupling /. Parameters`ApplyGUTNormalization[];
              (* convert MSbar to DRbar mass hep-ph/0207126 *)
              drbarConversion = 1 - alphaS / (3 Pi) - 23 / 72 alphaS^2 / Pi^2 + 3 SARAH`leftCoupling^2 / (128 Pi^2) + 13 gPrime^2 / (1152 Pi^2);
              If[dimParticle == 1,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_sm_msbar, int) const\n{\n";
                 body = "const double p = m_sm_msbar;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p));\n";
                 ,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_sm_msbar, int idx) const\n{\n";
                 body = "const double p = m_sm_msbar;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p, idx, idx));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p, idx, idx));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p, idx, idx));\n";
                ];
              body = body <>
              "const double m_tree = " <> RValueToCFormString[thirdGenMass] <> ";\n" <>
              "const double m_sm_drbar = m_sm_msbar * (" <> RValueToCFormString[drbarConversion] <> ");\n\n" <>
              "const double m_susy_drbar = m_sm_drbar / (1.0 - self_energy_1/m_tree " <>
              "- self_energy_PL - self_energy_PR);\n\n" <>
              "return m_susy_drbar;\n";
             ];
           Return[result <> IndentText[body] <> "}\n\n"];
          ];

CreateRunningDRbarMassFunction[particle_ /; particle === SARAH`Electron] :=
    Module[{result, body, selfEnergyFunctionS, selfEnergyFunctionPL,
            selfEnergyFunctionPR, name, drbarConversion, gPrime,
            dimParticle},
           dimParticle = TreeMasses`GetDimension[particle];
           selfEnergyFunctionS  = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[1]];
           selfEnergyFunctionPL = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[PL]];
           selfEnergyFunctionPR = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[PR]];
           name = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[IsMassless[particle],
              result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double, int) const\n{\n";
              body = "return 0.0;\n";
              ,
              (* convert MSbar to DRbar mass *)
              gPrime = SARAH`hyperchargeCoupling /. Parameters`ApplyGUTNormalization[];
              drbarConversion = 1 - 3 (gPrime^2 - SARAH`leftCoupling^2) / (128 Pi^2);
              If[dimParticle == 1,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_sm_msbar, int) const\n{\n";
                 body = "const double p = m_sm_msbar;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p));\n";
                 ,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_sm_msbar, int idx) const\n{\n";
                 body = "const double p = m_sm_msbar;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p, idx, idx));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p, idx, idx));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p, idx, idx));\n";
                ];
              body = body <>
              "const double m_sm_drbar = m_sm_msbar * (" <> RValueToCFormString[drbarConversion] <> ");\n\n" <>
              "const double m_susy_drbar = m_sm_drbar + self_energy_1 " <>
              "+ m_sm_drbar * (self_energy_PL + self_energy_PR);\n\n" <>
              "return m_susy_drbar;\n";
             ];
           Return[result <> IndentText[body] <> "}\n\n"];
          ];

CreateRunningDRbarMassFunction[particle_ /; particle === SARAH`TopQuark] :=
    Module[{result, body, selfEnergyFunctionS, selfEnergyFunctionPL,
            selfEnergyFunctionPR, name, qcdOneLoop, qcdTwoLoop,
            dimParticle, thirdGenMass},
           dimParticle = TreeMasses`GetDimension[particle];
           thirdGenMass = TreeMasses`GetThirdGenerationMass[particle];
           selfEnergyFunctionS  = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[1]];
           selfEnergyFunctionPL = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[PL]];
           selfEnergyFunctionPR = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[PR]];
           name = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[IsMassless[particle],
              result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double, int) const\n{\n";
              body = "return 0.0;\n";
              ,
              qcdOneLoop = - TwoLoop`GetDeltaMOverMQCDOneLoop[particle, Global`currentScale];
              qcdTwoLoop = N[Expand[qcdOneLoop^2 - TwoLoop`GetDeltaMOverMQCDTwoLoop[particle, Global`currentScale]]];
              If[dimParticle == 1,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_pole, int) const\n{\n";
                 body = "const double p = m_pole;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p));\n\n";
                 ,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_pole, int idx) const\n{\n";
                 body = "const double p = m_pole;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p, idx, idx));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p, idx, idx));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p, idx, idx));\n\n";
                ];
              body = body <>
              "const double currentScale = get_scale();\n" <>
              "const double qcd_1l = " <> CConversion`RValueToCFormString[qcdOneLoop /. FlexibleSUSY`M[particle] -> thirdGenMass] <> ";\n" <>
              "const double qcd_2l = " <> CConversion`RValueToCFormString[qcdTwoLoop /. FlexibleSUSY`M[particle] -> thirdGenMass] <> ";\n\n" <>
              "const double m_susy_drbar = m_pole + self_energy_1 " <>
              "+ m_pole * (self_energy_PL + self_energy_PR + qcd_1l + qcd_2l);\n\n" <>
              "return m_susy_drbar;\n";
             ];
           Return[result <> IndentText[body] <> "}\n\n"];
          ];

CreateRunningDRbarMassFunction[particle_ /; IsFermion[particle]] :=
    Module[{result, body, selfEnergyFunctionS, selfEnergyFunctionPL,
            selfEnergyFunctionPR, name,
            twoLoopCorrection, twoLoopCorrectionDecl = "", addTwoLoopCorrection = False,
            dimParticle},
           dimParticle = TreeMasses`GetDimension[particle];
           selfEnergyFunctionS  = SelfEnergies`CreateSelfEnergyFunctionName[particle[1]];
           selfEnergyFunctionPL = SelfEnergies`CreateSelfEnergyFunctionName[particle[PL]];
           selfEnergyFunctionPR = SelfEnergies`CreateSelfEnergyFunctionName[particle[PR]];
           name = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           twoLoopCorrection = 0; (* disable corrections until checked against Softsusy *)
           (* twoLoopCorrection = TwoLoop`GetDeltaMQCD[particle, Global`displayMu[]] /. *)
           (*                     FlexibleSUSY`M[p_] :> FlexibleSUSY`M[p[Global`idx]]; *)
           addTwoLoopCorrection = twoLoopCorrection =!= 0;
           If[addTwoLoopCorrection,
              twoLoopCorrectionDecl = "const double two_loop = " <> RValueToCFormString[twoLoopCorrection] <> ";\n";
             ];
           If[IsMassless[particle],
              result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double, int) const\n{\n";
              body = "return 0.0;\n";
              ,
              If[dimParticle == 1,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_pole, int) const\n{\n";
                 body = "const double p = m_pole;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p));\n";
                 ,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_pole, int idx) const\n{\n";
                 body = "const double p = m_pole;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p, idx, idx));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p, idx, idx));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p, idx, idx));\n";
                ];
              If[addTwoLoopCorrection,
                 body = body <> twoLoopCorrectionDecl;
                ];
              body = body <> "\n" <>
                     "const double m_drbar = m_pole + self_energy_1 + m_pole * (self_energy_PL + self_energy_PR)";
              If[addTwoLoopCorrection,
                 body = body <> " - two_loop";
                ];
              body = body <> ";\n\n";
              body = body <> "return m_drbar;\n";
             ];
           Return[result <> IndentText[body] <> "}\n\n"];
          ];

CreateRunningDRbarMassFunction[particle_] :=
    Module[{result, body, selfEnergyFunction, name},
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle];
           name = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[IsMassless[particle],
              result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double) const \n{\n";
              body = "return 0.0;\n";
              ,
              result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_pole) const\n{\n";
              body = "const double p = m_pole;\n" <>
              "const double self_energy = Re(" <> selfEnergyFunction <> "(p));\n" <>
              "const double mass_sqr = Sqr(m_pole) + self_energy;\n\n" <>
              "return ZeroSqrt(mass_sqr);\n";
             ];
           Return[result <> IndentText[body] <> "}\n\n"];
          ];

CreateRunningDRbarMassFunctions[] :=
    Module[{result = "", particles},
           particles = GetRunningOneLoopDRbarParticles[];
           (result = result <> CreateRunningDRbarMassFunction[#])& /@ particles;
           Return[result];
          ];

GetLightestMassEigenstate[FlexibleSUSY`M[mass_]] :=
    GetLightestMassEigenstate[mass];

GetLightestMassEigenstate[mass_] :=
    If[GetDimension[mass] == 1,
       FlexibleSUSY`M[mass],
       FlexibleSUSY`M[mass][GetDimensionStartSkippingGoldstones[mass] - 1]];

CreateLSPFunctions[{}] := {"", ""};

CreateLSPFunctions[masses_List] :=
    Module[{prototype, function, mass, info, particleType, m},
           info = FlexibleSUSY`FSModelName <> "_info";
           particleType = info <> "::Particles";
           body = "double lsp_mass = std::numeric_limits<double>::max();
double tmp_mass;
particle_type = " <> info <> "::NUMBER_OF_PARTICLES;

";
           For[m = 1, m <= Length[masses], m++,
               mass = masses[[m]];
               body = body <> "\
tmp_mass = Abs(PHYSICAL(" <>
CConversion`RValueToCFormString[GetLightestMassEigenstate[mass]] <>
"));
if (tmp_mass < lsp_mass) {
" <> IndentText["\
lsp_mass = tmp_mass;
particle_type = " <> info <> "::" <>
CConversion`ToValidCSymbolString[mass /. FlexibleSUSY`M -> Identity] <>
";"] <>
"
}

";
              ];
           body = body <> "return lsp_mass;\n";
           prototype = "double get_lsp(" <> particleType <> "&);\n";
           function = "double CLASSNAME::get_lsp(" <> particleType <>
                      "& particle_type)\n{\n" <> IndentText[body] <> "}\n";
           Return[{prototype, function}];
          ];

CreateLSPFunctions[_] := {"", ""};

End[];

EndPackage[];

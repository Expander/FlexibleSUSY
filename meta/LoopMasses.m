
BeginPackage["LoopMasses`", {"SARAH`", "TextFormatting`",
                             "CConversion`", "TreeMasses`",
                             "SelfEnergies`"}];

CreateLoopMassFunctions::usage="";
CreateLoopMassPrototypes::usage="";
CallAllLoopMassFunctions::usage="";
CreateRunningDRbarMassPrototypes::usage="";
CreateRunningDRbarMassFunctions::usage="";

GetLoopCorrectedParticles::usage="Returns list of all particles that
get loop corrected masses.  These are all particles, except for
ghosts.";

Begin["Private`"];

GetLoopCorrectedParticles[states_] :=
    Module[{particles},
           particles = GetParticles[states];
           Select[particles, (!IsGhost[#])&]
          ];

FillTadpoleMatrix[{}, _] := "";

FillTadpoleMatrix[tadpoles_List, matrixName_:"tadpoles"] :=
    Module[{result, dim, dimStr, particle, i, particleIndex, vev},
           particle = tadpoles[[1,1]];
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           If[dim > 1,
              result = "DoubleMatrix " <> matrixName <> "(" <> dimStr <> "," <> dimStr <> ");\n";
              For[i = 1, i <= Length[tadpoles], i++,
                  particleIndex = ToString[tadpoles[[i,2]]];
                  vev           = ToValidCSymbolString[tadpoles[[i,3]]];
                  result = result <>
                           matrixName <> "(" <> particleIndex <> "," <> particleIndex <> ") = " <>
                           CreateTadpoleFunctionName[particle] <> "(" <> particleIndex <> ").imag() / " <>
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

Do1DimScalar[particleName_String, selfEnergyFunction_String, momentum_String, tadpole_String:""] :=
    "const double p = " <> momentum <> ";\n" <>
    "const Complex self_energy = " <> selfEnergyFunction <> "(p)\n" <>
    "PHYSICAL(" <> particleName <> ") = ZeroSqrt(" <> particleName <>
    " - self_energy.real()" <> If[tadpole == "", "", " + " <> tadpole] <> ");\n";

Do1DimFermion[particleName_String, selfEnergyFunctionS_String,
              selfEnergyFunctionPL_String, selfEnergyFunctionPR_String, momentum_String] :=
    "const double p = " <> momentum <> ";\n" <>
    "const double self_energy_1  = " <> selfEnergyFunctionS  <> "(p).real();\n" <>
    "const double self_energy_PL = " <> selfEnergyFunctionPL <> "(p).real();\n" <>
    "const double self_energy_PR = " <> selfEnergyFunctionPR <> "(p).real();\n" <>
    "PHYSICAL(" <> particleName <> ") = " <> particleName <>
    " - self_energy_1 - " <> particleName <> " * (self_energy_PL + self_energy_PR);\n";

Do1DimVector[particleName_String, selfEnergyFunction_String, momentum_String] :=
    "const double p = " <> momentum <> ";\n" <>
    "const double self_energy = " <> selfEnergyFunction <> "(p).real();\n" <>
    "PHYSICAL(" <> particleName <> ") = ZeroSqrt(Sqr(" <> particleName <>
    ") - self_energy);\n";


(* ********** fast diagonalization routines ********** *)

DoFastDiagonalization[particle_Symbol /; IsScalar[particle], tadpoles_List] :=
    Module[{result, dim, dimStr, particleName, mixingMatrix, selfEnergyFunction,
            tadpoleMatrix, U, V},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           mixingMatrix = FindMixingMatrixSymbolFor[particle];
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle];
           tadpoleMatrix = FillTadpoleMatrix[tadpoles, "tadpoles"];
           If[dim > 1,
              result = tadpoleMatrix <>
                       "ComplexMatrix self_energy(" <> dimStr <> "," <> dimStr <> ");\n" <>
                       "for (unsigned i1 = 1; i1 <= " <> dimStr <>"; ++i1) {\n" <>
                       IndentText["for (unsigned i2 = 1; i2 <= " <> dimStr <>"; ++i2) {\n" <>
                                  IndentText["const double p = ZeroSqrt(" <> particleName <> "(i1) * " <> 
                                             particleName <> "(i2));\n" <>
                                             "self_energy(i1,i2) = " <>
                                             selfEnergyFunction <> "(p,i1,i2);\n"] <>
                                  "}\n"
                                 ] <>
                       "}\n" <>
                       "const DoubleMatrix M_1loop(get_mass_matrix_" <> particleName <>
                       "() - self_energy.real()" <>
                       If[tadpoleMatrix == "", "", " + tadpoles"] <> ");\n";
              If[Head[mixingMatrix] === List,
                 U = ToValidCSymbolString[mixingMatrix[[1]]];
                 V = ToValidCSymbolString[mixingMatrix[[2]]];
                 result = result <>
                          "DiagonaliseSVD(M_1loop, PHYSICAL(" <> U <> "), " <>
                          "PHYSICAL(" <> V <> "), " <>
                          "PHYSICAL(" <> particleName <> "));\n";
                 ,
                 U = ToValidCSymbolString[mixingMatrix];
                 result = result <>
                          "DiagonalizeUnsorted(M_1loop, PHYSICAL(" <> U <> "), " <>
                          "PHYSICAL(" <> particleName <> "));\n";
                ];
              result = result <>
                       "PHYSICAL(" <> particleName <> ") = PHYSICAL(" <>
                       particleName <> ").apply(ZeroSqrt);\n";
              ,
              result = Do1DimScalar[particleName, selfEnergyFunction, particleName, "tadpoles"];
             ];
           Return[result];
          ];

DoFastDiagonalization[particle_Symbol /; IsFermion[particle], _] :=
    Module[{result, dim, dimStr, particleName, mixingMatrix, U, V,
            selfEnergyFunctionS, selfEnergyFunctionPL, selfEnergyFunctionPR},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           If[IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0,
              Return["PHYSICAL(" <> particleName <> ") = 0;\n"];
             ];
           mixingMatrix = FindMixingMatrixSymbolFor[particle];
           selfEnergyFunctionS  = SelfEnergies`CreateSelfEnergyFunctionName[particle[1]];
           selfEnergyFunctionPL = SelfEnergies`CreateSelfEnergyFunctionName[particle[PL]];
           selfEnergyFunctionPR = SelfEnergies`CreateSelfEnergyFunctionName[particle[PR]];
           If[dim > 1,
              result = "ComplexMatrix self_energy_1(" <> dimStr <> "," <> dimStr <> ");\n" <>
                       "ComplexMatrix self_energy_PL(" <> dimStr <> "," <> dimStr <> ");\n" <>
                       "ComplexMatrix self_energy_PR(" <> dimStr <> "," <> dimStr <> ");\n" <>
                       "for (unsigned i1 = 1; i1 <= " <> dimStr <>"; ++i1) {\n" <>
                       IndentText["for (unsigned i2 = 1; i2 <= " <> dimStr <>"; ++i2) {\n" <>
                                  IndentText["const double p = ZeroSqrt(" <> particleName <> "(i1) * " <> 
                                             particleName <> "(i2));\n" <>
                                             "self_energy_1(i1,i2) = " <>
                                             selfEnergyFunctionS <> "(p,i1,i2);\n" <>
                                             "self_energy_PL(i1,i2) = " <>
                                             selfEnergyFunctionPL <> "(p,i1,i2);\n" <>
                                             "self_energy_PR(i1,i2) = " <>
                                             selfEnergyFunctionPR <> "(p,i1,i2);\n"
                                            ] <>
                                  "}\n"
                                 ] <>
                       "}\n" <>
                       "const DoubleMatrix M_tree(get_mass_matrix_" <> particleName <> "());\n" <>
                       "const DoubleMatrix delta_M(- self_energy_PR.real() * M_tree " <>
                       "- M_tree * self_energy_PL.real() - self_energy_1.real());\n";
              If[IsMajoranaFermion[particle],
                 result = result <>
                          "const DoubleMatrix M_1loop(M_tree + 0.5 * (delta_M + delta_M.transpose()));\n";
                 ,
                 result = result <>
                          "const DoubleMatrix M_1loop(M_tree + delta_M);\n";
                ];
              If[Head[mixingMatrix] === List,
                 (* two mixing matrixs => SVD *)
                 U = ToValidCSymbolString[mixingMatrix[[1]]];
                 V = ToValidCSymbolString[mixingMatrix[[2]]];
                 If[dim == 2,
                    result = result <>
                             "Diagonalize2by2(M_1loop, " <>
                             "PHYSICAL(" <> U <> "), " <>
                             "PHYSICAL(" <> V <> "), " <>
                             "PHYSICAL(" <> particleName <> "));\n";
                    ,
                    result = result <>
                             "Diagonalize(M_1loop, " <>
                             "PHYSICAL(" <> U <> "), " <>
                             "PHYSICAL(" <> V <> "), " <>
                             "PHYSICAL(" <> particleName <> "));\n";
                   ];
                 ,
                 result = result <>
                          "Diagonalize(M_1loop, " <>
                          "PHYSICAL(" <> ToValidCSymbolString[mixingMatrix] <> "), " <>
                          "PHYSICAL(" <> particleName <> "));\n";
                ];
              ,
              (* for a dimension 1 fermion it plays not role if it's a
                 Majorana ferimion or not *)
              result = Do1DimFermion[particleName, selfEnergyFunctionS,
                                     selfEnergyFunctionPL, selfEnergyFunctionPR,
                                     particleName];
             ];
           Return[result];
          ];

DoFastDiagonalization[particle_Symbol /; IsVector[particle], _] :=
    Module[{result, dim, dimStr, particleName, mixingMatrix, selfEnergyFunction},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           mixingMatrix = ToValidCSymbolString[FindMixingMatrixSymbolFor[particle]];
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle];
           If[IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0,
              Return["PHYSICAL(" <> particleName <> ") = 0;\n"];
             ];
           If[dim > 1,
              result = "WARNING(\"diagonalization of " <> ToString[particle] <> " not implemented\");\n";
              ,
              result = Do1DimVector[particleName, selfEnergyFunction, particleName];
             ];
           Return[result];
          ];

DoFastDiagonalization[particle_Symbol, _] :=
    "ERROR(\"fast diagonalization of " <> ToString[particle] <> " not implemented\");\n";

(* ********** medium diagonalization routines ********** *)

DoMediumDiagonalization[particle_Symbol /; IsScalar[particle], inputMomentum_, tadpole_List] :=
    Module[{result, dim, dimStr, particleName, mixingMatrix, selfEnergyFunction,
            momentum = inputMomentum, U, V, Utemp, Vtemp, tadpoleMatrix, diagSnippet},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           If[inputMomentum == "", momentum = particleName];
           mixingMatrix = FindMixingMatrixSymbolFor[particle];
           (* create diagonalisation code snippet *)
           If[Head[mixingMatrix] === List,
              U = ToValidCSymbolString[mixingMatrix[[1]]];
              V = ToValidCSymbolString[mixingMatrix[[2]]];
              Utemp = "mix_" <> U;
              Vtemp = "mix_" <> V;
              diagSnippet = "DoubleMatrix " <> Utemp <> "(" <> dimStr <> "," <> dimStr <> "), " <>
                            Vtemp <> "(" <> dimStr <> "," <> dimStr <> ");\n" <>
                            "Diagonalise(M_1loop, " <> Utemp <> ", " <> Vtemp <> ", eigen_values);\n" <>
                            "PHYSICAL(" <> particleName <> "(es)) = ZeroSqrt(eigen_values(es));\n" <>
                            "if (es == 1) {\n" <>
                            IndentText["PHYSICAL(" <> U <> ") = " <> Utemp <> ";\n" <>
                                       "PHYSICAL(" <> V <> ") = " <> Vtemp <> ";\n"] <>
                            "}\n";
              ,
              U = ToValidCSymbolString[mixingMatrix];
              Utemp = "mix_" <> U;
              diagSnippet = "DoubleMatrix " <> Utemp <> "(" <> dimStr <> "," <> dimStr <> ");\n" <>
                            "DiagonalizeUnsorted(M_1loop, " <> Utemp <> ", eigen_values);\n" <>
                            "PHYSICAL(" <> particleName <> "(es)) = ZeroSqrt(eigen_values(es));\n" <>
                            "if (es == 1)\n" <>
                            IndentText["PHYSICAL(" <> U <> ") = " <> Utemp <> ";\n"];
             ];
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle];
           tadpoleMatrix = FillTadpoleMatrix[tadpole, "tadpoles"];
           (* fill self-energy and do diagonalisation *)
           If[dim > 1,
              result = tadpoleMatrix <>
                       "ComplexMatrix self_energy(" <> dimStr <> "," <> dimStr <> ");\n" <>
                       "const DoubleMatrix M_tree(get_mass_matrix_" <> particleName <> "());\n" <>
                       "for (unsigned es = 1; es <= " <> dimStr <> "; ++es) {\n" <>
                       IndentText["const double p = std::fabs(" <> momentum <> "(es));\n" <>
                                  "for (unsigned i1 = 1; i1 <= " <> dimStr <> "; ++i1) {\n" <>
                                  IndentText["for (unsigned i2 = 1; i2 <= " <> dimStr <> "; ++i2) {\n" <>
                                             IndentText["self_energy(i1,i2) = " <>
                                                        selfEnergyFunction <> "(p,i1,i2);\n"
                                                       ] <>
                                             "}\n"
                                            ] <>
                                  "}\n" <>
                                  "const DoubleMatrix M_1loop(M_tree - self_energy.real()" <>
                                  If[tadpoleMatrix == "", "", " + tadpoles"] <> ");\n" <>
                                  "DoubleVector eigen_values(" <> dimStr <> ");\n" <>
                                  diagSnippet
                                 ] <>
                       "}\n";
              ,
              result = tadpoleMatrix <>
                       Do1DimScalar[particleName, selfEnergyFunction, momentum, "tadpoles"];
             ];
           Return[result];
          ];

DoMediumDiagonalization[particle_Symbol /; IsFermion[particle], inputMomentum_, _] :=
    Module[{result, dim, dimStr, particleName, mixingMatrix, U, V,
            selfEnergyFunctionS, selfEnergyFunctionPL, selfEnergyFunctionPR,
            momentum = inputMomentum},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           If[inputMomentum == "", momentum = particleName];
           If[IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0,
              Return["PHYSICAL(" <> particleName <> ") = 0;\n"];
             ];
           mixingMatrix = FindMixingMatrixSymbolFor[particle];
           selfEnergyFunctionS  = SelfEnergies`CreateSelfEnergyFunctionName[particle[1]];
           selfEnergyFunctionPL = SelfEnergies`CreateSelfEnergyFunctionName[particle[PL]];
           selfEnergyFunctionPR = SelfEnergies`CreateSelfEnergyFunctionName[particle[PR]];
           If[dim > 1,
              result = "ComplexMatrix self_energy_1(" <> dimStr <> "," <> dimStr <> ");\n" <>
                       "ComplexMatrix self_energy_PL(" <> dimStr <> "," <> dimStr <> ");\n" <>
                       "ComplexMatrix self_energy_PR(" <> dimStr <> "," <> dimStr <> ");\n" <>
                       "const DoubleMatrix M_tree(get_mass_matrix_" <> particleName <> "());\n" <>
                       "for (unsigned es = 1; es <= " <> dimStr <> "; ++es) {\n" <>
                       IndentText["const double p = std::fabs(" <> momentum <> "(es));\n" <>
                                  "for (unsigned i1 = 1; i1 <= " <> dimStr <>"; ++i1) {\n" <>
                                  IndentText["for (unsigned i2 = 1; i2 <= " <> dimStr <>"; ++i2) {\n" <>
                                             IndentText["self_energy_1(i1,i2) = " <>
                                                        selfEnergyFunctionS <> "(p,i1,i2);\n" <>
                                                        "self_energy_PL(i1,i2) = " <>
                                                        selfEnergyFunctionPL <> "(p,i1,i2);\n" <>
                                                        "self_energy_PR(i1,i2) = " <>
                                                        selfEnergyFunctionPR <> "(p,i1,i2);\n"
                                                       ] <>
                                             "}\n"
                                            ] <>
                                  "}\n" <>
                                  "const DoubleMatrix delta_M(- self_energy_PR.real() * M_tree " <>
                                  "- M_tree * self_energy_PL.real() - self_energy_1.real());\n"
                                 ];
              If[IsMajoranaFermion[particle],
                 result = result <>
                          IndentText["const DoubleMatrix M_1loop(M_tree + 0.5 * (delta_M + delta_M.transpose()));\n"];
                 ,
                 result = result <>
                          IndentText["const DoubleMatrix M_1loop(M_tree + delta_M);\n"];
                ];
              result = result <>
                       IndentText["DoubleVector eigen_values(" <> dimStr <> ");\n"];
              If[Head[mixingMatrix] === List,
                 (* two mixing matrixs => SVD *)
                 U = ToValidCSymbolString[mixingMatrix[[1]]];
                 V = ToValidCSymbolString[mixingMatrix[[2]]];
                 result = result <>
                          IndentText["decltype(" <> U <> ") mix_" <> U <> "(" <> dimStr <> "," <> dimStr <> ");\n" <>
                                     "decltype(" <> V <> ") mix_" <> V <> "(" <> dimStr <> "," <> dimStr <> ");\n"];
                 If[dim == 2,
                    result = result <>
                             IndentText["Diagonalize2by2(M_1loop, " <>
                                        "mix_" <> U <> ", " <> "mix_" <> V <> ", eigen_values);\n"];
                    ,
                    result = result <>
                             IndentText["Diagonalize(M_1loop, " <>
                                        "mix_" <> U <> ", " <> "mix_" <> V <> ", eigen_values);\n"];
                   ];
                 result = result <>
                          IndentText["if (es == 1) {\n" <>
                                     IndentText["PHYSICAL(" <> U <> ") = mix_" <> U <> ";\n" <>
                                                "PHYSICAL(" <> V <> ") = mix_" <> V <> ";\n"] <>
                                     "}\n"
                                    ];
                 ,
                 U = ToValidCSymbolString[mixingMatrix];
                 result = result <>
                          IndentText["decltype(" <> U <> ") mix_" <> U <> "(" <> dimStr <> "," <> dimStr <> ");\n" <>
                                     "Diagonalize(M_1loop, mix_" <> U <> ", eigen_values);\n" <>
                                     "if (es == 1)\n" <>
                                     IndentText["PHYSICAL(" <> U <> ") = mix_" <> U <> ";\n"]
                                    ];
                ];
              result = result <>
                       IndentText["PHYSICAL(" <> particleName <>
                                  "(es)) = std::fabs(eigen_values(es));\n"];
              result = result <> "}\n";
              ,
              (* for a dimension 1 fermion it plays not role if it's a
                 Majorana fermion or not *)
              result = Do1DimFermion[particleName, selfEnergyFunctionS,
                                     selfEnergyFunctionPL, selfEnergyFunctionPR,
                                     momentum];
             ];
           Return[result];
          ];

DoMediumDiagonalization[particle_Symbol /; IsVector[particle], inputMomentum_, _] :=
    Module[{result, dim, dimStr, particleName, mixingMatrix, selfEnergyFunction,
            momentum = inputMomentum},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           If[inputMomentum == "", momentum = particleName];
           mixingMatrix = ToValidCSymbolString[FindMixingMatrixSymbolFor[particle]];
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle];
           If[IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0,
              Return["PHYSICAL(" <> particleName <> ") = 0;\n"];
             ];
           If[dim > 1,
              result = "WARNING(\"diagonalization of " <> ToString[particle] <> " not implemented\");\n";
              ,
              result = Do1DimVector[particleName, selfEnergyFunction, momentum];
             ];
           Return[result];
          ];

DoMediumDiagonalization[particle_Symbol, inputMomentum_:"", _] :=
    "ERROR(\"medium diagonalization of " <> ToString[particle] <> " not implemented\");\n";

(* ********** high precision diagonalization routines ********** *)

DoSlowDiagonalization[particle_Symbol, tadpole_] :=
    Module[{result, dim, dimStr, particleName, inputMomenta, outputMomenta, body},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           inputMomenta = "old_" <> particleName;
           outputMomenta = "new_" <> particleName;
           body = DoMediumDiagonalization[particle, inputMomenta, tadpole] <> "\n" <>
                  outputMomenta <> " = PHYSICAL(" <> particleName <> ");\n" <>
                  "diff = MaxRelDiff(" <> outputMomenta <> ", " <> inputMomenta <> ");\n" <>
                  inputMomenta <> " = " <> outputMomenta <> ";\n" <>
                  "iteration++;\n";
           result = "unsigned iteration = 0;\n" <>
                    "double diff = 0.0;\n" <>
                    "decltype(" <> particleName <> ") " <>
                    inputMomenta  <> "(" <> particleName <> "), " <>
                    outputMomenta <> "(" <> particleName <> ");\n\n" <>
                    "do {\n" <>
                    IndentText[body] <>
                    "} while (diff > mass_iteration_precision\n" <>
                    "         && iteration < number_of_mass_iterations);\n";
           Return[result];
          ];

DoDiagonalization[particle_Symbol, FlexibleSUSY`LowPrecision, tadpole_] :=
    "// diagonalization with low precision\n" <> DoFastDiagonalization[particle, tadpole];

DoDiagonalization[particle_Symbol, FlexibleSUSY`MediumPrecision, tadpole_] :=
    "// diagonalization with medium precision\n" <> DoMediumDiagonalization[particle, "", tadpole];

DoDiagonalization[particle_Symbol, FlexibleSUSY`HighPrecision, tadpole_] :=
    "// diagonalization with high precision\n" <> DoSlowDiagonalization[particle, tadpole];

CreateLoopMassFunction[particle_Symbol, precision_Symbol, tadpole_] :=
    Module[{result, body = ""},
           body = DoDiagonalization[particle, precision, tadpole];
           result = "void CLASSNAME::calculate_" <> ToString[particle] <>
                    "_onshell_1loop()\n{\n" <> IndentText[body] <> "}\n\n";
           Return[result];
          ];

CreateLoopMassFunctions[precision_List, oneLoopTadpoles_List, vevs_List] :=
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

CreateLoopMassPrototype[particle_Symbol] :=
    "void calculate_" <> ToString[particle] <> "_onshell_1loop();\n";

CreateLoopMassPrototypes[states_:SARAH`EWSB] :=
    Module[{particles, result = ""},
           particles = GetLoopCorrectedParticles[states];
           (result = result <> CreateLoopMassPrototype[#])& /@ particles;
           Return[result];
          ];

CallLoopMassFunction[particle_Symbol] :=
    "calculate_" <> ToString[particle] <> "_onshell_1loop();\n";

CallAllLoopMassFunctions[states_:SARAH`EWSB] :=
    Module[{particles, result = ""},
           particles = GetLoopCorrectedParticles[states];
           (result = result <> CallLoopMassFunction[#])& /@ particles;
           Return[result];
          ];

GetRunningOneLoopDRbarParticles[] :=
    {SARAH`TopQuark, SARAH`BottomQuark, SARAH`Electron, SARAH`Neutrino,
     SARAH`VectorP, SARAH`VectorZ, SARAH`VectorW};

CreateRunningDRbarMassPrototype[particle_ /; IsFermion[particle]] :=
    "double calculate_" <> ToValidCSymbolString[particle] <> "_DRbar_1loop(double, int) const;\n";

CreateRunningDRbarMassPrototype[particle_] :=
    "double calculate_" <> ToValidCSymbolString[particle] <> "_DRbar_1loop(double) const;\n";

CreateRunningDRbarMassPrototypes[] :=
    Module[{result = "", particles},
           particles = GetRunningOneLoopDRbarParticles[];
           (result = result <> CreateRunningDRbarMassPrototype[#])& /@ particles;
           Return[result];
          ];

CreateRunningDRbarMassFunction[particle_ /; IsFermion[particle]] :=
    Module[{result, body, selfEnergyFunctionS, selfEnergyFunctionPL,
            selfEnergyFunctionPR, name},
           selfEnergyFunctionS  = SelfEnergies`CreateSelfEnergyFunctionName[particle[1]];
           selfEnergyFunctionPL = SelfEnergies`CreateSelfEnergyFunctionName[particle[PL]];
           selfEnergyFunctionPR = SelfEnergies`CreateSelfEnergyFunctionName[particle[PR]];
           name = ToValidCSymbolString[particle];
           If[IsMassless[particle],
              result = "double CLASSNAME::calculate_" <> name <> "_DRbar_1loop(double, int) const\n{\n";
              body = "return 0.0;\n";
              ,
              result = "double CLASSNAME::calculate_" <> name <> "_DRbar_1loop(double m_onshell, int index) const\n{\n";
              body = "const double p = m_onshell;\n" <>
              "const double self_energy_1  = " <> selfEnergyFunctionS  <> "(p, index, index).real();\n" <>
              "const double self_energy_PL = " <> selfEnergyFunctionPL <> "(p, index, index).real();\n" <>
              "const double self_energy_PR = " <> selfEnergyFunctionPR <> "(p, index, index).real();\n" <>
              "return m_onshell + self_energy_1 + m_onshell * (self_energy_PL + self_energy_PR);\n";
             ];
           Return[result <> IndentText[body] <> "}\n\n"];
          ];

CreateRunningDRbarMassFunction[particle_] :=
    Module[{result, body, selfEnergyFunction, name},
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle];
           name = ToValidCSymbolString[particle];
           If[IsMassless[particle],
              result = "double CLASSNAME::calculate_" <> name <> "_DRbar_1loop(double) const \n{\n";
              body = "return 0.0;\n";
              ,
              result = "double CLASSNAME::calculate_" <> name <> "_DRbar_1loop(double m_onshell) const\n{\n";
              body = "const double p = m_onshell;\n" <>
              "const double self_energy = " <> selfEnergyFunction <> "(p).real();\n" <>
              "return ZeroSqrt(Sqr(m_onshell) + self_energy);\n";
             ];
           Return[result <> IndentText[body] <> "}\n\n"];
          ];

CreateRunningDRbarMassFunctions[] :=
    Module[{result = "", particles},
           particles = GetRunningOneLoopDRbarParticles[];
           (result = result <> CreateRunningDRbarMassFunction[#])& /@ particles;
           Return[result];
          ];

End[];

EndPackage[];


BeginPackage["WriteOut`", {"SARAH`", "TextFormatting`", "CConversion`",
                           "Parameters`", "TreeMasses`", "LatticeUtils`"}];

ReplaceInFiles::usage="Replaces tokens in files.";
PrintParameters::usage="Creates parameter printout statements";
WriteSLHAExtparBlock::usage="";
WriteSLHAMassBlock::usage="";
WriteSLHAMixingMatricesBlocks::usage="";
WriteSLHAModelParametersBlocks::usage="";
WriteSLHAMinparBlock::usage="";
ReadLesHouchesInputParameters::usage="";
ReadLesHouchesOutputParameters::usage="";
ReadLesHouchesPhysicalParameters::usage="";
ConvertMixingsToSLHAConvention::usage="";
GetDRbarBlockNames::usage="";
GetNumberOfDRbarBlocks::usage="";
StringJoinWithSeparator::usage="Joins a list of strings with a given separator string";
ParseCmdLineOptions::usage="";
PrintCmdLineOptions::usage="";

Begin["`Private`"];

StringJoinWithSeparator[list_List, separator_String, transformer_:ToString] :=
    Module[{result = "", i},
           For[i = 1, i <= Length[list], i++,
               If[i > 1, result = result <> separator;];
               result = result <> transformer[list[[i]]];
              ];
           Return[result];
          ];

(*
 * @brief Replaces tokens in files.
 *
 * @param files list of two-element lists.  The first entry is the
 * input file and the second entry is the output file.
 *
 * Example:
 *    files = {{"input.hpp", "output.hpp"},
 *             {"input.cpp", "output.cpp"}}
 *
 * @param replacementList list of string replacement rules
 *
 * Example:
 *    replacementList = { "@token@" -> "1+2", "@bar@" -> "2+3" }
 *)
ReplaceInFiles[files_List, replacementList_List] :=
    Module[{cppFileName, cppTemplateFileName, cppFile, modifiedCppFile, f},
          For[f = 1, f <= Length[files], f++,
              cppFileName         = files[[f,1]];
              cppTemplateFileName = files[[f,2]];
              cppFile             = Import[cppFileName, "String"];
              modifiedCppFile     = StringReplace[cppFile, replacementList];
              Print["   Writing file ", cppTemplateFileName];
              Export[cppTemplateFileName, modifiedCppFile, "String"];
             ];
          ];

TransposeIfVector[parameter_, CConversion`ArrayType[__]] :=
    SARAH`Tp[parameter];

TransposeIfVector[parameter_, CConversion`VectorType[__]] :=
    SARAH`Tp[parameter];

TransposeIfVector[parameter_, _] := parameter;

PrintParameter[Null, streamName_String] := "";

PrintParameter[parameter_, streamName_String] :=
    Module[{parameterName, parameterNameWithoutIndices, expr, type},
           parameterNameWithoutIndices = parameter /.
                                         a_[Susyno`LieGroups`i1,SARAH`i2] :> a;
           parameterName = CConversion`ToValidCSymbolString[parameterNameWithoutIndices];
           type = Parameters`GetType[parameterNameWithoutIndices];
           expr = TransposeIfVector[parameterNameWithoutIndices, type];
           Return[streamName <> " << \"" <> parameterName <> " = \" << " <>
                  CConversion`RValueToCFormString[expr] <> " << '\\n';\n"];
          ];

PrintParameters[parameters_List, streamName_String] :=
    Module[{result = ""},
           (result = result <> PrintParameter[#,streamName])& /@ parameters;
           Return[result];
          ];

WriteSLHAMass[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result = "", eigenstateName, eigenstateNameStr, massNameStr,
            pdgList, pdg, dim, i},
           eigenstateName = TreeMasses`GetMassEigenstate[massMatrix];
           dim = TreeMasses`GetDimension[eigenstateName];
           pdgList = SARAH`getPDGList[eigenstateName];
           If[Length[pdgList] != dim,
              Print["Error: length of PDG number list != dimension of particle ", eigenstateName];
              Print["       PDG number list = ", pdgList];
              Print["       dimension of particle ", eigenstateName, " = ", dim];
             ];
           If[Length[pdgList] < dim,
              Return[""];
             ];
           If[dim == 1,
              pdg = Abs[pdgList[[1]]];
              If[pdg != 0,
                 eigenstateNameStr = CConversion`RValueToCFormString[eigenstateName];
                 massNameStr = CConversion`RValueToCFormString[FlexibleSUSY`M[eigenstateName]];
                 result = "<< FORMAT_MASS(" <> ToString[pdg] <>
                          ", LOCALPHYSICAL(" <> massNameStr <> "), \"" <> eigenstateNameStr <> "\")\n";
                ];
              ,
              For[i = 1, i <= dim, i++,
                  pdg = Abs[pdgList[[i]]];
                  If[pdg != 0,
                     eigenstateNameStr = CConversion`RValueToCFormString[eigenstateName] <> "_" <> ToString[i];
                     massNameStr = CConversion`RValueToCFormString[FlexibleSUSY`M[eigenstateName[i-1]]];
                     result = result <> "<< FORMAT_MASS(" <> ToString[pdg] <>
                              ", LOCALPHYSICAL(" <> massNameStr <> "), \"" <> eigenstateNameStr <> "\")\n";
                    ];
                 ];
             ];
           Return[result];
          ];

WriteSLHAMassBlock[massMatrices_List] :=
    Module[{result, allMasses, smMasses, susyMasses,
            smMassesStr = "", susyMassesStr = ""},
           allMasses = FlexibleSUSY`M[TreeMasses`GetMassEigenstate[#]]& /@ massMatrices;
           smMasses = Select[massMatrices, (SARAH`SMQ[TreeMasses`GetMassEigenstate[#]])&];
           susyMasses = Complement[massMatrices, smMasses];
           (smMassesStr = smMassesStr <> WriteSLHAMass[#])& /@ smMasses;
           (susyMassesStr = susyMassesStr <> WriteSLHAMass[#])& /@ susyMasses;
           susyMassesStr = "mass << \"Block MASS\\n\"\n" <>
                           TextFormatting`IndentText[susyMassesStr] <> ";\n\n";
           smMassesStr = "if (write_sm_masses) {\n" <>
                         TextFormatting`IndentText["mass\n" <>
                             TextFormatting`IndentText[smMassesStr] <> ";"] <>
                         "\n}\n\n";
           result = "std::ostringstream mass;\n\n" <>
                    susyMassesStr <> smMassesStr <>
                    "slha_io.set_block(mass);\n";
           Return[result];
          ];

WriteParameterTuple[{key_?NumberQ, parameter_}, streamName_String] :=
    Module[{parameterStr},
           parameterStr = CConversion`ToValidCSymbolString[parameter];
           streamName <> " << FORMAT_ELEMENT(" <> ToString[key] <> ", input." <>
           parameterStr <> ", \"" <> parameterStr <> "\");\n"
          ];

WriteParameterTuple[expr_, _] :=
    Block[{},
          Print["Error: not a valid {key,parameter} tuple: ", expr];
          ""
         ];

WriteSLHAExtparBlock[{}] := "";

WriteSLHAExtparBlock[extpar_List] :=
    Module[{result, body = ""},
           (body = body <> WriteParameterTuple[#, "extpar"])& /@ extpar;
           result = "std::ostringstream extpar;\n\n" <>
                    "extpar << \"Block EXTPAR\\n\";\n" <>
                    body <>
                    "slha_io.set_block(extpar);\n";
           Return[result];
          ];

WriteSLHAMinparBlock[{}] := "";

WriteSLHAMinparBlock[minpar_List] :=
    Module[{result, body = ""},
           (body = body <> WriteParameterTuple[#, "minpar"])& /@ minpar;
           result = "std::ostringstream minpar;\n\n" <>
                    "minpar << \"Block MINPAR\\n\";\n" <>
                    body <>
                    "slha_io.set_block(minpar);\n";
           Return[result];
          ];

GetSLHAMixinMatrices[] :=
    DeleteCases[Select[FlexibleSUSY`FSLesHouchesList,
                       MemberQ[Parameters`GetOutputParameters[],#[[1]]]&],
                {_,None}];

GetSLHAModelParameters[] :=
    DeleteCases[Select[FlexibleSUSY`FSLesHouchesList,
                       MemberQ[Parameters`GetModelParameters[],#[[1]]]&],
                {_,None}];

WriteSLHAMatrix[{mixingMatrix_, lesHouchesName_}, head_String] :=
    WriteSLHAMatrix[{mixingMatrix, lesHouchesName}, head, ""];

WriteSLHAMatrix[{mixingMatrix_, lesHouchesName_}, head_String, scale_String] :=
    Module[{str, lhs, wrapper, dim, diagMatrix},
           If[SARAH`getDimParameters[mixingMatrix] === {} ||
              SARAH`getDimParameters[mixingMatrix] === {1},
              Print["Warning: You are trying to create a SLHA matrix block for"];
              Print["   ", mixingMatrix, ", which is not a matrix!"];
              Print["   Please specify a Les Houches index in the SARAH model file."];
             ];
           str = CConversion`ToValidCSymbolString[mixingMatrix];
           lhs = ToString[lesHouchesName];
           wrapper = If[head == "", str, head <> "(" <> str <> ")"];
           If[mixingMatrix === SARAH`UpYukawa ||
              mixingMatrix === SARAH`DownYukawa ||
              mixingMatrix === SARAH`ElectronYukawa
              ,
              (* In SLHA-2 the diagonalized Yukawa couplings must be
                 output.  See SLHA-2 standard, Sec. 4.1.3, p.15 *)
              dim = SARAH`getDimParameters[mixingMatrix][[1]];
              diagMatrix = str <> "_diag";
              "{\n" <> IndentText[
                  "// diagonalize the " <> lhs <> " Yukawa coupling for SLHA-2 compatible output\n" <>
                  CreateCType[CConversion`ArrayType[CConversion`realScalarCType, dim]] <> " " <>
                  diagMatrix <> ";\n" <>
                  CreateCType[CConversion`MatrixType[CConversion`complexScalarCType, dim, dim]] <> " " <>
                  " u, v;\n" <>
                  "fs_svd(" <> wrapper <> ", " <> diagMatrix <> ", u, v);\n" <>
                  "slha_io.set_block(\"" <> lhs <> "\", " <>
                  CreateCType[CConversion`MatrixType[CConversion`realScalarCType, dim, dim]] <>
                  "(" <> diagMatrix <> ".matrix().asDiagonal()), \"" <> lhs <> "\", model.get_scale());"
              ] <> "\n}\n"
              ,
              "slha_io.set_block(\"" <> lhs <> "\", " <> wrapper <> ", \"" <> str <>
              "\"" <> If[scale != "", ", " <> scale, ""] <> ");\n"
             ]
          ];

WriteSLHAMixingMatricesBlocks[] :=
    Module[{result, mixingMatrices, smMix, susyMix, smMixStr = "", susyMixStr = ""},
           mixingMatrices = GetSLHAMixinMatrices[];
           smMix = Flatten[TreeMasses`FindMixingMatrixSymbolFor /@ SARAH`SMParticles];
           smMix = Select[mixingMatrices, MemberQ[smMix,#[[1]]]&];
           susyMix = Complement[mixingMatrices, smMix];
           (smMixStr = smMixStr <> WriteSLHAMatrix[#,"LOCALPHYSICAL"])& /@ smMix;
           (susyMixStr = susyMixStr <> WriteSLHAMatrix[#,"LOCALPHYSICAL"])& /@ susyMix;
           result = susyMixStr <> "\n" <>
                    "if (write_sm_mixing_matrics) {\n" <>
                    TextFormatting`IndentText[smMixStr] <> "}\n";
           Return[result];
          ];

LesHouchesNameToFront[{parameter_, {lh_,idx_}}] :=
    {lh, {parameter, idx}};

LesHouchesNameToFront[{parameter_, lh_}] :=
    {lh, parameter};

SortBlocks[modelParameters_List] :=
    Module[{reformed, allBlocks, collected},
           reformed = LesHouchesNameToFront /@ modelParameters;
           allBlocks = DeleteDuplicates[Transpose[reformed][[1]]];
           collected = {#, Cases[reformed, {#, a_} :> a]}& /@ allBlocks;
           Return[collected];
          ];

WriteSLHABlock[{blockName_, tuples_List}] :=
    Module[{result = "", blockNameStr, t, pdg, par, parStr, parVal, scale},
           blockNameStr = ToString[blockName];
           scale = "model.get_scale()";
           result = "std::ostringstream block;\n" <>
                    "block << \"Block " <> blockNameStr <> " Q= \" << FORMAT_NUMBER(" <>
                    scale <> ") << '\\n'\n";
           For[t = 1, t <= Length[tuples], t++,
               If[Head[tuples[[t]]] =!= List || Length[tuples[[t]]] < 2,
                  Print["WriteSLHABlock: Error: tuple ", tuples[[t]],
                        " is not a list of lenght 2."];
                  Print["  Tuples list: ", tuples];
                  Continue[];
                 ];
               par = tuples[[t,1]];
               parStr = CConversion`ToValidCSymbolString[par];
               parVal = "MODELPARAMETER(" <> parStr <> ")";
               (* print unnormalized hypercharge gauge coupling *)
               If[par === SARAH`hyperchargeCoupling,
                  parVal = "(" <> parVal <> " * " <>
                           CConversion`RValueToCFormString[
                               Parameters`GetGUTNormalization[par]] <> ")";
                  parStr = "gY";
                 ];
               pdg = ToString[tuples[[t,2]]];
               result = result <> "      << FORMAT_ELEMENT(" <> pdg <> ", " <>
                        parVal <> ", \"" <> parStr <> "\")" <>
                        If[t == Length[tuples], ";", ""] <> "\n";
              ];
           result = result <> "slha_io.set_block(block);\n";
           result = "{\n" <> TextFormatting`IndentText[result] <> "}\n";
           Return[result];
          ];

WriteSLHABlock[{blockName_, parameter_}] :=
    WriteSLHAMatrix[{parameter, blockName}, "MODELPARAMETER", "model.get_scale()"];

WriteSLHABlock[{blockName_, {parameter_ /; Head[parameter] =!= List}}] :=
    WriteSLHABlock[{blockName, parameter}];

WriteSLHAModelParametersBlocks[] :=
    Module[{result = "", modelParameters, blocks},
           modelParameters = GetSLHAModelParameters[];
           blocks = SortBlocks[modelParameters];
           (result = result <> WriteSLHABlock[#])& /@ blocks;
           Return[result];
          ];

ReadSLHAInputBlock[{parameter_, {blockName_Symbol, pdg_?NumberQ}}] :=
    Module[{result, blockNameStr, parmStr, pdgStr},
           blockNameStr = ToString[blockName] <> "IN";
           parmStr = CConversion`ToValidCSymbolString[parameter];
           pdgStr = ToString[pdg];
           result = "input." <> parmStr <>
                    " = slha_io.read_entry(\"" <> blockNameStr <> "\", " <>
                    pdgStr <> ");\n";
           Return[result];
          ];

ReadSLHAInputBlock[{parameter_, blockName_Symbol}] :=
    Module[{paramStr, blockNameStr},
           paramStr = CConversion`ToValidCSymbolString[parameter];
           blockNameStr = ToString[blockName] <> "IN";
           "slha_io.read_block(\"" <> blockNameStr <> "\", input." <>
           paramStr <> ");\n"
          ];

ReadLesHouchesInputParameters[lesHouchesInputParameters_List] :=
    Module[{result = "", modelParameters, names, rules},
           names = (#[[1]])& /@ lesHouchesInputParameters;
           rules = Rule[#[[1]], #[[2]]]& /@ lesHouchesInputParameters;
           (* get block names of all les Houches input parameters (names) *)
           modelParameters = Select[GetSLHAModelParameters[], MemberQ[names,#[[1]]]&];
           modelParameters = {#[[1]] /. rules, #[[2]]}& /@ modelParameters;
           (result = result <> ReadSLHAInputBlock[#])& /@ modelParameters;
           Return[result];
          ];

ReadSLHAOutputBlock[{parameter_, {blockName_Symbol, pdg_?NumberQ}}] :=
    Module[{result, blockNameStr, parmStr, pdgStr, gutNorm = ""},
           blockNameStr = ToString[blockName];
           parmStr = CConversion`ToValidCSymbolString[parameter];
           pdgStr = ToString[pdg];
           If[parameter === SARAH`hyperchargeCoupling,
              gutNorm = " * " <> CConversion`RValueToCFormString[
                  1/Parameters`GetGUTNormalization[parameter]];
             ];
           result = "model.set_" <> parmStr <>
                    "(slha_io.read_entry(\"" <> blockNameStr <> "\", " <>
                    pdgStr <> ")" <> gutNorm <> ");\n";
           Return[result];
          ];

ReadSLHAOutputBlock[{parameter_, blockName_Symbol}] :=
    Module[{paramStr, blockNameStr},
           paramStr = CConversion`ToValidCSymbolString[parameter];
           blockNameStr = ToString[blockName];
           "{\n" <> IndentText[
               "DEFINE_PARAMETER(" <> paramStr <> ");\n" <>
               "slha_io.read_block(\"" <> blockNameStr <> "\", " <>
               paramStr <> ");\n" <>
               "model.set_" <> paramStr <> "(" <> paramStr <> ");"] <> "\n" <>
           "}\n"
          ];

ReadSLHAPhysicalMixingMatrixBlock[{parameter_, blockName_Symbol}] :=
    Module[{paramStr, blockNameStr},
           paramStr = CConversion`ToValidCSymbolString[parameter];
           blockNameStr = ToString[blockName];
           "{\n" <> IndentText[
               "DEFINE_PARAMETER(" <> paramStr <> ");\n" <>
               "slha_io.read_block(\"" <> blockNameStr <> "\", " <>
               paramStr <> ");\n" <>
               "PHYSICAL(" <> paramStr <> ") = " <> paramStr <> ";"] <> "\n" <>
           "}\n"
          ];

ReadSLHAPhysicalMass[particle_] :=
    Module[{result = "", mass, massStr, dim, pdgList, pdg, pdgStr, i},
           mass = FlexibleSUSY`M[particle];
           massStr = CConversion`ToValidCSymbolString[mass];
           dim = TreeMasses`GetDimension[particle];
           pdgList = SARAH`getPDGList[particle];
           If[Head[pdgList] =!= List || Length[pdgList] < dim,
              Return[""];
             ];
           If[dim == 1,
              pdg = Abs[pdgList[[1]]];
              pdgStr = ToString[pdg];
              If[pdg != 0,
                 result = "PHYSICAL(" <> massStr <>
                          ") = slha_io.read_entry(\"MASS\", " <> pdgStr <> ");\n";
                ];
              ,
              For[i = 1, i <= dim, i++,
                  pdg = Abs[pdgList[[i]]];
                  pdgStr = ToString[pdg];
                  If[pdg != 0,
                     result = result <>
                              "PHYSICAL(" <> massStr <> ")(" <> ToString[i-1] <>
                              ") = slha_io.read_entry(\"MASS\", " <> pdgStr <> ");\n";
                    ];
                 ];
             ];
           Return[result];
          ];

ReadSLHAPhysicalMassBlock[] :=
    Module[{result = "", particles},
           particles = TreeMasses`GetParticles[];
           (result = result <> ReadSLHAPhysicalMass[#])& /@ particles;
           Return[result];
          ];

ReadLesHouchesOutputParameters[] :=
    Module[{result = "", modelParameters},
           modelParameters = GetSLHAModelParameters[];
           (result = result <> ReadSLHAOutputBlock[#])& /@ modelParameters;
           Return[result];
          ];

ReadLesHouchesPhysicalParameters[] :=
    Module[{result = "", physicalParameters},
           physicalParameters = GetSLHAMixinMatrices[];
           (result = result <> ReadSLHAPhysicalMixingMatrixBlock[#])& /@ physicalParameters;
           result = result <> "\n" <> ReadSLHAPhysicalMassBlock[];
           Return[result];
          ];

GetDRbarBlocks[] :=
    Module[{modelParameters},
           modelParameters = GetSLHAModelParameters[];
           DeleteDuplicates[Cases[modelParameters, {_, blockName_Symbol} | {_, {blockName_Symbol, _?NumberQ}} :> blockName]]
          ];

GetDRbarBlockNames[] :=
    Module[{blocks, transformer},
           blocks = GetDRbarBlocks[];
           transformer = ("\"" <> ToString[#] <> "\"")&;
           "{ " <> StringJoinWithSeparator[blocks, ", ", transformer] <> " }"
          ];

GetNumberOfDRbarBlocks[] := Length[GetDRbarBlocks[]];

ConvertMixingsToSLHAConvention[massMatrices_List] :=
    Module[{result = "", i,
            eigenstateName, mixingMatrixSym,
            eigenstateNameStr, mixingMatrixSymStr},
           For[i = 1, i <= Length[massMatrices], i++,
               eigenstateName = TreeMasses`GetMassEigenstate[massMatrices[[i]]];
               mixingMatrixSym = TreeMasses`GetMixingMatrixSymbol[massMatrices[[i]]];
               If[LatticeUtils`MajoranaMassMatrixQ[massMatrices[[i]]],
                  eigenstateNameStr  = CConversion`ToValidCSymbolString[FlexibleSUSY`M[eigenstateName]];
                  mixingMatrixSymStr = CConversion`ToValidCSymbolString[mixingMatrixSym];
                  result = result <>
                           "SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(" <>
                           eigenstateNameStr <> "), LOCALPHYSICAL(" <>
                           mixingMatrixSymStr <> "));\n";
                 ];
              ];
           Return[result];
          ];

ParseCmdLineOption[parameter_Symbol] :=
    Module[{parameterStr},
           parameterStr = ToValidCSymbolString[parameter];
           "\
if(Command_line_options::get_parameter_value(option, \"--" <> parameterStr <> "=\", input." <> parameterStr <>"))
   continue;

"
          ];

ParseCmdLineOption[FlexibleSUSY`Sign[phase_]] :=
    Module[{parameterStr},
           parameterStr = ToValidCSymbolString[FlexibleSUSY`Sign[phase]];
           "\
if(Command_line_options::get_parameter_value(option, \"--" <> parameterStr <> "=\", input." <> parameterStr <>"))
   continue;

"
          ];

ParseCmdLineOption[_] := "";

ParseCmdLineOptions[inputParameters_List] :=
    StringJoin[ParseCmdLineOption /@ inputParameters];

PrintCmdLineOption[parameter_Symbol] :=
    "\"  --" <> ToValidCSymbolString[parameter] <> "=<value>\\n\"\n";

PrintCmdLineOption[FlexibleSUSY`Sign[phase_]] :=
    "\"  --" <> ToValidCSymbolString[FlexibleSUSY`Sign[phase]] <> "=<value>\\n\"\n";

PrintCmdLineOption[_] := "";

PrintCmdLineOptions[inputParameters_List] :=
    StringJoin[PrintCmdLineOption /@ inputParameters];

End[];

EndPackage[];

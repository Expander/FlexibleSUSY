
BeginPackage["WriteOut`", {"SARAH`", "TextFormatting`", "CConversion`",
                           "Parameters`", "TreeMasses`", "LatticeUtils`",
                           "Utils`"}];

ReplaceInFiles::usage="Replaces tokens in files.";
PrintParameters::usage="Creates parameter printout statements";
PrintInputParameters::usage="Creates input parameter printout statements";
WriteSLHAExtparBlock::usage="";
WriteSLHAMassBlock::usage="";
WriteSLHAMixingMatricesBlocks::usage="";
WriteSLHAModelParametersBlocks::usage="";
WriteSLHAMinparBlock::usage="";
WriteExtraSLHAOutputBlock::usage="";
ReadLesHouchesInputParameters::usage="";
ReadLesHouchesOutputParameters::usage="";
ReadLesHouchesPhysicalParameters::usage="";
ConvertMixingsToSLHAConvention::usage="";
GetDRbarBlockNames::usage="";
GetNumberOfDRbarBlocks::usage="";
ParseCmdLineOptions::usage="";
PrintCmdLineOptions::usage="";
GetGaugeCouplingNormalizationsDecls::usage="";
GetGaugeCouplingNormalizationsDefs::usage="";

CreateSLHAYukawaDefinition::usage="";
CreateSLHAYukawaGetters::usage="";
ConvertYukawaCouplingsToSLHA::usage="";
CreateSLHAFermionMixingMatricesDef::usage="";
CreateSLHAFermionMixingMatricesGetters::usage=""
CreateSLHATrilinearCouplingDefinition::usage="";
CreateSLHATrilinearCouplingGetters::usage="";
ConvertTrilinearCouplingsToSLHA::usage="";
CreateSLHASoftSquaredMassesDefinition::usage="";
CreateSLHASoftSquaredMassesGetters::usage="";
ConvertSoftSquaredMassesToSLHA::usage="";

CalculateCKMMatrix::usage="";
CalculatePMNSMatrix::usage="";

Begin["`Private`"];

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

TransposeIfVector[p:FlexibleSUSY`Sign[parameter_], _] :=
    CConversion`ToValidCSymbolString[p];

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

PrintInputParameter[Null, _] := "";

PrintInputParameter[{parameter_, type_}, streamName_String] :=
    Module[{parameterStr, expr},
           parameterStr = CConversion`ToValidCSymbolString[parameter];
           expr = TransposeIfVector[parameter, type];
           Return[streamName <> " << \"" <> parameterStr <> " = \" << " <>
                  "INPUT(" <> CConversion`RValueToCFormString[expr] <> ") << \", \";\n"];
          ];

PrintInputParameter[parameter_, streamName_String] :=
    PrintInputParameter[{parameter, Parameters`GetType[parameter]}, streamName];

PrintInputParameters[parameters_List, streamName_String] :=
    Module[{result = ""},
           (result = result <> PrintInputParameter[#,streamName])& /@ parameters;
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
                     eigenstateNameStr = CConversion`RValueToCFormString[eigenstateName] <> "(" <> ToString[i] <> ")";
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
           (* filter out MW, because MW should always appear in the output *)
           smMasses = Select[smMasses, (TreeMasses`GetMassEigenstate[#] =!= SARAH`VectorW)&];
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

GetSLHAInputParameters[] :=
    DeleteCases[Select[FlexibleSUSY`FSLesHouchesList,
                       MemberQ[Parameters`GetInputParameters[],#[[1]]]&],
                {_,None}];

WriteSLHAMatrix[{mixingMatrix_, lesHouchesName_}, head_String] :=
    WriteSLHAMatrix[{mixingMatrix, lesHouchesName}, head, ""];

WriteSLHAMatrix[{mixingMatrix_, lesHouchesName_}, head_String, scale_String] :=
    Module[{str, strSLHA, lhs, wrapper},
           If[SARAH`getDimParameters[mixingMatrix] === {} ||
              SARAH`getDimParameters[mixingMatrix] === {1},
              Print["Warning: You are trying to create a SLHA matrix block for"];
              Print["   ", mixingMatrix, ", which is not a matrix!"];
              Print["   Please specify a Les Houches index in the SARAH model file."];
             ];
           str = CConversion`ToValidCSymbolString[mixingMatrix];
           (* use SLHA compliant yukawas, trilinears, soft-squared masses *)
           strSLHA = If[mixingMatrix === SARAH`UpYukawa ||
                        mixingMatrix === SARAH`DownYukawa ||
                        mixingMatrix === SARAH`ElectronYukawa ||
                        mixingMatrix === SARAH`TrilinearUp ||
                        mixingMatrix === SARAH`TrilinearDown ||
                        mixingMatrix === SARAH`TrilinearLepton ||
                        mixingMatrix === SARAH`SoftSquark ||
                        mixingMatrix === SARAH`SoftUp ||
                        mixingMatrix === SARAH`SoftDown ||
                        mixingMatrix === SARAH`SoftLeftLepton ||
                        mixingMatrix === SARAH`SoftRightLepton,
                        str <> "_slha",
                        str
                       ];
           lhs = ToString[lesHouchesName];
           wrapper = If[head == "", strSLHA, head <> "(" <> strSLHA <> ")"];
           (* convert SLHA yukawa vectors to matrix *)
           wrapper = If[mixingMatrix === SARAH`UpYukawa ||
                        mixingMatrix === SARAH`DownYukawa ||
                        mixingMatrix === SARAH`ElectronYukawa,
                        "ToMatrix(" <> wrapper <> ")",
                        wrapper
                       ];
           "slha_io.set_block(\"" <> lhs <> "\", " <> wrapper <> ", \"" <> str <>
           "\"" <> If[scale != "", ", " <> scale, ""] <> ");\n"
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

CreateRulesForProtectedHead[expr_, protectedHead_Symbol] :=
    Cases[expr, protectedHead[p__] :> Rule[protectedHead[p],Symbol["x$" <> ToString[Hash[p]]]], {0, Infinity}];

CreateRulesForProtectedHead[expr_, protectedHeads_List] :=
    Flatten @ Join[CreateRulesForProtectedHead[expr,#]& /@ protectedHeads];

WrapPreprocessorMacroAround[expr_String, ___] := expr;

WrapPreprocessorMacroAround[expr_, symbols_, macroSymbol_,
                             protectedHeads_List:{FlexibleSUSY`Pole, SARAH`SM}] :=
    Module[{replacements, protectionRules, exprWithoutProtectedSymbols},
           replacements = Join[
               RuleDelayed[#     , macroSymbol[#]   ]& /@ symbols,
               RuleDelayed[#[i__], macroSymbol[#][i]]& /@ symbols,
               {RuleDelayed[FlexibleSUSY`M[p_[i__]], macroSymbol[FlexibleSUSY`M[p]][i]]}
           ];
           protectionRules = CreateRulesForProtectedHead[expr, protectedHeads];
           exprWithoutProtectedSymbols = expr /. protectionRules;
           (* substitute back protected symbols *)
           exprWithoutProtectedSymbols /. replacements /. (Reverse /@ protectionRules)
          ];

SetAttributes[WriteSLHABlockEntry, HoldFirst];

WriteSLHABlockEntry[{Hold[par_], idx___}, comment_String:""] :=
    Module[{parStr, commentStr},
           parStr = ToString[Unevaluated[par]];
           commentStr = If[comment == "", parStr, comment];
           parStr = StringReplace[parStr,
                                  {RegularExpression["\\bSUSYScale\\b"] -> "SCALES(SUSYScale)",
                                   RegularExpression["\\bHighScale\\b"] -> "SCALES(HighScale)",
                                   RegularExpression["\\bLowScale\\b"]  -> "SCALES(LowScale)"}
                                 ];
           WriteSLHABlockEntry[{parStr, idx}, commentStr]
          ];

ClearAttributes[WriteSLHABlockEntry, HoldFirst];

WriteSLHABlockEntry[{par_, idx1_?NumberQ, idx2_?NumberQ}, comment_String:""] :=
    Module[{parStr, parVal, idx1Str, idx2Str, commentStr},
           parStr = CConversion`RValueToCFormString[par];
           parVal = CConversion`RValueToCFormString[
               WrapPreprocessorMacroAround[par, Join[Parameters`GetModelParameters[],
                                                     Parameters`GetOutputParameters[]],
                                           Global`MODELPARAMETER]];
           idx1Str = ToString[idx1];
           idx2Str = ToString[idx2];
           commentStr = If[comment == "", parStr, comment];
           (* result *)
           "      << FORMAT_MIXING_MATRIX(" <> idx1Str <> ", " <> idx2Str <>
           ", (" <> parVal <> "), \"" <> commentStr <> "\")" <> "\n"
          ];

WriteSLHABlockEntry[{par_, pdg_?NumberQ}, comment_String:""] :=
    Module[{parStr, parVal, pdgStr, commentStr},
           parStr = CConversion`RValueToCFormString[par];
           parVal = CConversion`RValueToCFormString[
               WrapPreprocessorMacroAround[par, Join[Parameters`GetModelParameters[],
                                                     Parameters`GetOutputParameters[]],
                                           Global`MODELPARAMETER]];
           (* print unnormalized hypercharge gauge coupling *)
           If[par === SARAH`hyperchargeCoupling,
              parVal = parVal <> " * " <>
                           CConversion`RValueToCFormString[
                               Parameters`GetGUTNormalization[par]];
              parStr = "gY";
             ];
           pdgStr = ToString[pdg];
           commentStr = If[comment == "", parStr, comment];
           (* result *)
           "      << FORMAT_ELEMENT(" <> pdgStr <> ", (" <> parVal <>
           "), \"" <> commentStr <> "\")" <> "\n"
          ];

WriteSLHABlockEntry[{par_}, comment_String:""] :=
    Module[{parStr, parVal, commentStr},
           parStr = CConversion`RValueToCFormString[par];
           parVal = CConversion`RValueToCFormString[
               WrapPreprocessorMacroAround[par, Join[Parameters`GetModelParameters[],
                                                     Parameters`GetOutputParameters[]],
                                           Global`MODELPARAMETER]];
           commentStr = If[comment == "", parStr, comment];
           (* result *)
           "      << FORMAT_NUMBER((" <> parVal <> "), \"" <> commentStr <> "\")\n"
          ];

WriteSLHABlockEntry[tuple___] :=
    Block[{},
          Print["WriteSLHABlockEntry: Error: malformed entry ", tuple];
          ""
         ];

WriteSLHABlock[{blockName_, tuples_List}] :=
    Module[{result = "", blockNameStr, scale},
           blockNameStr = ToString[blockName];
           scale = "model.get_scale()";
           result = "std::ostringstream block;\n" <>
                    "block << \"Block " <> blockNameStr <> " Q= \" << FORMAT_SCALE(" <>
                    scale <> ") << '\\n'\n";
           (result = result <> WriteSLHABlockEntry[#])& /@ tuples;
           result = result <> ";\n" <> "slha_io.set_block(block);\n";
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

WriteExtraSLHAOutputBlock[outputBlocks_List] :=
    Module[{result = "", reformed},
           ReformeBlocks[{block_, tuples_List}] := {block, ReformeBlocks /@ tuples};
           ReformeBlocks[{expr_}]               := {expr};
           ReformeBlocks[{idx_, expr_}]         := {expr, idx};
           ReformeBlocks[{idx1_, idx2_, expr_}] := {expr, idx1, idx2};
           reformed = ReformeBlocks /@ outputBlocks;
           (result = result <> WriteSLHABlock[#])& /@ reformed;
           Return[result];
          ];

ReadSLHAInputBlock[{parameter_, {blockName_, pdg_?NumberQ}}] :=
    Module[{result, blockNameStr, parmStr, pdgStr},
           blockNameStr = ToString[blockName];
           parmStr = CConversion`ToValidCSymbolString[parameter];
           pdgStr = ToString[pdg];
           result = "input." <> parmStr <>
                    " = slha_io.read_entry(\"" <> blockNameStr <> "\", " <>
                    pdgStr <> ");\n";
           Return[result];
          ];

ReadSLHAInputBlock[{parameter_, blockName_}] :=
    Module[{paramStr, blockNameStr},
           paramStr = CConversion`ToValidCSymbolString[parameter];
           blockNameStr = ToString[blockName];
           "slha_io.read_block(\"" <> blockNameStr <> "\", input." <>
           paramStr <> ");\n"
          ];

CreateInputBlockName[{blockName_, pdg_?NumberQ}] :=
    {ToString[blockName] <> "IN", pdg};

CreateInputBlockName[blockName_] :=
    ToString[blockName] <> "IN";

ReadLesHouchesInputParameters[lesHouchesInputParameters_List] :=
    Module[{result = "", parameters, names, rules},
           names = (#[[1]])& /@ lesHouchesInputParameters;
           rules = Cases[lesHouchesInputParameters, {p_, block_, _} /; MemberQ[Parameters`GetModelParameters[],p] :> Rule[p,block]];
           (* get block names of all les Houches input parameters (names) *)
           parameters = Select[Join[GetSLHAModelParameters[],GetSLHAInputParameters[]], MemberQ[names,#[[1]]]&];
           parameters = {#[[1]] /. rules,
                         If[MemberQ[Parameters`GetModelParameters[],#[[1]]], CreateInputBlockName[#[[2]]], #[[2]]]}& /@ parameters;
           (result = result <> ReadSLHAInputBlock[#])& /@ parameters;
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
           "{ " <> Utils`StringJoinWithSeparator[blocks, ", ", transformer] <> " }"
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

(* SLHA CKM conversion *)

GetYukawas[] :=
    Select[{SARAH`UpYukawa, SARAH`DownYukawa, SARAH`ElectronYukawa},
           MemberQ[Parameters`GetModelParameters[],#]&];

GetFermionMixingMatrices[] :=
    Select[{SARAH`DownMatrixL, SARAH`UpMatrixL,
            SARAH`DownMatrixR, SARAH`UpMatrixR,
            SARAH`ElectronMatrixL, SARAH`ElectronMatrixR,
            SARAH`NeutrinoMM},
           MemberQ[Parameters`GetOutputParameters[],#]&];

GetMixingMatricesFor[yuk_] :=
    Switch[yuk,
           SARAH`UpYukawa       , {SARAH`UpMatrixL      , SARAH`UpMatrixR      },
           SARAH`DownYukawa     , {SARAH`DownMatrixL    , SARAH`DownMatrixR    },
           SARAH`ElectronYukawa , {SARAH`ElectronMatrixL, SARAH`ElectronMatrixR},
           SARAH`TrilinearUp    , {SARAH`UpMatrixL      , SARAH`UpMatrixR      },
           SARAH`TrilinearDown  , {SARAH`DownMatrixL    , SARAH`DownMatrixR    },
           SARAH`TrilinearLepton, {SARAH`ElectronMatrixL, SARAH`ElectronMatrixR},
           SARAH`SoftSquark     , {SARAH`DownMatrixL    , SARAH`DownMatrixR    },
           SARAH`SoftUp         , {SARAH`UpMatrixL      , SARAH`UpMatrixR      },
           SARAH`SoftDown       , {SARAH`DownMatrixL    , SARAH`DownMatrixR    },
           SARAH`SoftLeftLepton , {SARAH`ElectronMatrixL, SARAH`ElectronMatrixR},
           SARAH`SoftRightLepton, {SARAH`ElectronMatrixL, SARAH`ElectronMatrixR}
          ];

IsLeftHanded[c_] :=
    Switch[c,
           SARAH`SoftSquark     , True,
           SARAH`SoftUp         , False,
           SARAH`SoftDown       , False,
           SARAH`SoftLeftLepton , True,
           SARAH`SoftRightLepton, False
          ];

GetTrilinearCouplings[] :=
    Select[{SARAH`TrilinearUp, SARAH`TrilinearDown, SARAH`TrilinearLepton},
           MemberQ[Parameters`GetModelParameters[],#]&];

GetSoftSquaredMasses[] :=
    Select[{SARAH`SoftSquark, SARAH`SoftUp, SARAH`SoftDown,
            SARAH`SoftLeftLepton, SARAH`SoftRightLepton},
           MemberQ[Parameters`GetModelParameters[],#]&];

(* SLHA Yukawa couplings *)

CreateSLHAYukawaName[yuk_] :=
    CConversion`ToValidCSymbolString[yuk] <> "_slha";

GetSLHAYukawaType[yuk_] :=
    CConversion`ArrayType[CConversion`realScalarCType,
                          SARAH`getDimParameters[yuk][[1]]];

CreateSLHAFermionMixingMatrixName[m_] :=
    CConversion`ToValidCSymbolString[m] <> "_slha";

GetSLHAFermionMixingMatrixType[m_] :=
    CConversion`MatrixType[CConversion`complexScalarCType,
                           SARAH`getDimParameters[m][[1]],
                           SARAH`getDimParameters[m][[2]]];

CreateSLHATrilinearCouplingName[c_] :=
    CConversion`ToValidCSymbolString[c] <> "_slha";

GetSLHATrilinearCouplingType[c_] :=
    Parameters`GetType[c];

CreateSLHASoftSquaredMassName[c_] :=
    CConversion`ToValidCSymbolString[c] <> "_slha";

GetSLHASoftSquaredMassType[c_] :=
    Parameters`GetType[c];

CreateSLHAYukawaDefinition[] :=
    Module[{result = "", yuks},
           yuks = GetYukawas[];
           Block[{},
               result = result <>
                        CConversion`CreateCType[GetSLHAYukawaType[#]] <>
                        " " <> CreateSLHAYukawaName[#] <> ";\n";
           ]& /@ yuks;
           result
          ];

CreateSLHAYukawaGetters[] :=
    Module[{result = "", yuks},
           yuks = GetYukawas[];
           Block[{},
              result = result <>
                       CConversion`CreateInlineGetter[
                           CreateSLHAYukawaName[#], GetSLHAYukawaType[#]
                       ] <>
                       CConversion`CreateInlineElementGetter[
                           CreateSLHAYukawaName[#], GetSLHAYukawaType[#]
                       ];
           ]& /@ yuks;
           result
          ];

ConvertYukawaCouplingsToSLHA[] :=
    Module[{result = ""},
           yuks = GetYukawas[];
           Module[{dim, vL, vR},
                  dim = SARAH`getDimParameters[#][[1]];
                  {vL, vR} = GetMixingMatricesFor[#];
                  result = result <>
                           "fs_svd(" <> CConversion`ToValidCSymbolString[#] <> ", " <>
                                   CreateSLHAYukawaName[#] <> ", " <>
                                   CreateSLHAFermionMixingMatrixName[vR] <> ", " <>
                                   CreateSLHAFermionMixingMatrixName[vL] <> ");\n";
           ]& /@ yuks;
           result
          ];

(* SLHA fermion mixing matrices *)

CreateSLHAFermionMixingMatricesDef[] :=
    Module[{result = "", yuks},
           yuks = GetFermionMixingMatrices[];
           Block[{},
                 result = result <>
                          CConversion`CreateCType[GetSLHAFermionMixingMatrixType[#]] <>
                          " " <> CreateSLHAFermionMixingMatrixName[#] <> ";\n";
           ]& /@ yuks;
           result
          ];

CreateSLHATrilinearCouplingDefinition[] :=
    Module[{result = "", tril},
           tril = GetTrilinearCouplings[];
           Block[{},
               result = result <>
                        CConversion`CreateCType[GetSLHATrilinearCouplingType[#]] <>
                        " " <> CreateSLHATrilinearCouplingName[#] <> ";\n";
           ]& /@ tril;
           result
          ];

CreateSLHATrilinearCouplingGetters[] :=
    Module[{result = "", tril},
           tril = GetTrilinearCouplings[];
           Block[{},
              result = result <>
                       CConversion`CreateInlineGetter[
                           CreateSLHATrilinearCouplingName[#], GetSLHATrilinearCouplingType[#]
                       ] <>
                       CConversion`CreateInlineElementGetter[
                           CreateSLHATrilinearCouplingName[#], GetSLHATrilinearCouplingType[#]
                       ];
           ]& /@ tril;
           result
          ];

ConvertTrilinearCouplingsToSLHA[] :=
    Module[{result = "", tril},
           tril = GetTrilinearCouplings[];
           Module[{vL, vR},
                  {vL, vR} = GetMixingMatricesFor[#];
                  result = result <>
                           CreateSLHATrilinearCouplingName[#] <> " = (" <>
                           CreateSLHAFermionMixingMatrixName[vR] <> ".conjugate() * " <>
                           CConversion`ToValidCSymbolString[#] <> " * " <>
                           CreateSLHAFermionMixingMatrixName[vL] <> ".adjoint()" <>
                           ").real();\n";
           ]& /@ tril;
           result
          ];

CreateSLHAFermionMixingMatricesGetters[] :=
    Module[{result = "", mix},
           mix = GetFermionMixingMatrices[];
           Block[{},
              result = result <>
                       CConversion`CreateInlineGetter[
                           CreateSLHAFermionMixingMatrixName[#], GetSLHAFermionMixingMatrixType[#]
                       ] <>
                       CConversion`CreateInlineElementGetter[
                           CreateSLHAFermionMixingMatrixName[#], GetSLHAFermionMixingMatrixType[#]
                       ];
           ]& /@ mix;
           result
          ];

CreateSLHASoftSquaredMassesDefinition[] :=
    Module[{result = "", massSq},
           massSq = GetSoftSquaredMasses[];
           Block[{},
               result = result <>
                        CConversion`CreateCType[GetSLHASoftSquaredMassType[#]] <>
                        " " <> CreateSLHASoftSquaredMassName[#] <> ";\n";
           ]& /@ massSq;
           result
          ];

CreateSLHASoftSquaredMassesGetters[] :=
    Module[{result = "", massSq},
           massSq = GetSoftSquaredMasses[];
           Block[{},
              result = result <>
                       CConversion`CreateInlineGetter[
                           CreateSLHASoftSquaredMassName[#], GetSLHASoftSquaredMassType[#]
                       ] <>
                       CConversion`CreateInlineElementGetter[
                           CreateSLHASoftSquaredMassName[#], GetSLHASoftSquaredMassType[#]
                       ];
           ]& /@ massSq;
           result
          ];

ConvertSoftSquaredMassesToSLHA[] :=
    Module[{result = "", massSq},
           massSq = GetSoftSquaredMasses[];
           Module[{vL, vR},
                  {vL, vR} = GetMixingMatricesFor[#];
                  If[IsLeftHanded[#],
                     result = result <>
                              CreateSLHASoftSquaredMassName[#] <> " = (" <>
                              CreateSLHAFermionMixingMatrixName[vL] <> " * " <>
                              CConversion`ToValidCSymbolString[#] <> " * " <>
                              CreateSLHAFermionMixingMatrixName[vL] <> ".adjoint()" <>
                              ").real();\n";
                     ,
                     result = result <>
                              CreateSLHASoftSquaredMassName[#] <> " = (" <>
                              CreateSLHAFermionMixingMatrixName[vR] <> ".conjugate() * " <>
                              CConversion`ToValidCSymbolString[#] <> " * " <>
                              CreateSLHAFermionMixingMatrixName[vR] <> ".transpose()" <>
                              ").real();\n";
                    ];
           ]& /@ massSq;
           result
          ];

CalculateCKMMatrix[] :=
    Module[{result = ""},
           If[MemberQ[Parameters`GetOutputParameters[], SARAH`DownMatrixL] &&
              MemberQ[Parameters`GetOutputParameters[], SARAH`UpMatrixL]
              ,
              result = result <> "ckm = " <>
              CreateSLHAFermionMixingMatrixName[SARAH`UpMatrixL  ] <> " * " <>
              CreateSLHAFermionMixingMatrixName[SARAH`DownMatrixL] <> ".adjoint();\n";
             ];
           (* convert CKM matrix to PDG convention *)
           If[MemberQ[Parameters`GetOutputParameters[], SARAH`DownMatrixL] &&
              MemberQ[Parameters`GetOutputParameters[], SARAH`UpMatrixL  ] &&
              MemberQ[Parameters`GetOutputParameters[], SARAH`DownMatrixR] &&
              MemberQ[Parameters`GetOutputParameters[], SARAH`UpMatrixR  ]
              ,
              result = result <> "CKM_parameters::to_pdg_convention(ckm, " <>
              CreateSLHAFermionMixingMatrixName[SARAH`UpMatrixL  ] <> ", " <>
              CreateSLHAFermionMixingMatrixName[SARAH`DownMatrixL] <> ", " <>
              CreateSLHAFermionMixingMatrixName[SARAH`UpMatrixR  ] <> ", " <>
              CreateSLHAFermionMixingMatrixName[SARAH`DownMatrixR] <> ");\n";
             ];
           result
          ];

CalculatePMNSMatrix[] :=
    Module[{result = ""},
           If[MemberQ[Parameters`GetOutputParameters[], SARAH`ElectronMatrixL] &&
              MemberQ[Parameters`GetOutputParameters[], SARAH`NeutrinoMM] &&
              SARAH`getDimParameters[SARAH`ElectronMatrixL] === SARAH`getDimParameters[SARAH`NeutrinoMM]
              ,
              result = "pmns = " <>
              CreateSLHAFermionMixingMatrixName[SARAH`ElectronMatrixL] <> " * " <>
              CreateSLHAFermionMixingMatrixName[SARAH`NeutrinoMM] <> ".adjoint();\n";
              ,
              result = "pmns << 1, 0, 0, 0, 1, 0, 0, 0, 1;\n";
             ];
           result
          ];

GetGaugeCouplingNormalizationsDecls[gauge_List] :=
    StringJoin[
        CConversion`CreateConstExternDecl[
            "normalization_" <> CConversion`ToValidCSymbolString[#[[4]]],
            CConversion`ScalarType[realScalarCType]
        ]& /@ gauge
    ];

GetGaugeCouplingNormalizationsDefs[gauge_List] :=
    StringJoin[
        CConversion`CreateConstDef[
            "normalization_" <> CConversion`ToValidCSymbolString[#[[4]]],
            CConversion`ScalarType[realScalarCType],
            Parameters`GetGUTNormalization[#[[4]]]
        ]& /@ gauge
    ];

End[];

EndPackage[];

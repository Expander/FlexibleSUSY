
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
CreateSLHAMassBlockStream::usage="creates ostringstream with masses";
ReadLesHouchesInputParameters::usage="";
ReadLesHouchesOutputParameters::usage="";
ReadLesHouchesPhysicalParameters::usage="";
ConvertMixingsToSLHAConvention::usage="";
ConvertMixingsToHKConvention::usage="";
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

CreateInputBlockName::usage="Creates an SLHA input block name for a
 given SLHA block name";

CreateFormattedSLHABlocks::usage = "";

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

TransposeIfVector[p:FlexibleSUSY`Phase[parameter_], _] :=
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

PrintInputParameters[parameters_List, streamName_String] :=
    Module[{result = ""},
           (result = result <> PrintInputParameter[#,streamName])& /@ parameters;
           Return[result];
          ];

WriteSLHAMass[p:TreeMasses`FSMassMatrix[_,massESSymbols_List,_]] :=
    Module[{massMatrices},
           massMatrices = DeleteDuplicates[TreeMasses`FSMassMatrix[0, #, Null]& /@ massESSymbols];
           StringJoin[WriteSLHAMass /@ massMatrices]
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

CreateSLHAMassBlockStream[massMatrices_List, blockName_String:"MASS", streamName_String:"mass"] :=
    Module[{smMasses, susyMasses,
            smMassesStr = "", susyMassesStr = ""},
           smMasses = Select[massMatrices, (IsSMParticle[TreeMasses`GetMassEigenstate[#]])&];
           (* filter out MW, because MW should always appear in the output *)
           smMasses = Select[smMasses, (TreeMasses`GetMassEigenstate[#] =!= SARAH`VectorW)&];
           susyMasses = Complement[massMatrices, smMasses];
           (smMassesStr = smMassesStr <> WriteSLHAMass[#])& /@ smMasses;
           (susyMassesStr = susyMassesStr <> WriteSLHAMass[#])& /@ susyMasses;
           susyMassesStr = streamName <> " << \"Block " <> blockName <> "\\n\"\n" <>
                           TextFormatting`IndentText[susyMassesStr] <> ";\n\n";
           smMassesStr = "if (write_sm_masses) {\n" <>
                         TextFormatting`IndentText[streamName <> "\n" <>
                             TextFormatting`IndentText[smMassesStr] <> ";"] <>
                         "\n}\n\n";
           "std::ostringstream " <> streamName <> ";\n\n" <> susyMassesStr <> smMassesStr
          ];

WriteSLHAMassBlock[massMatrices_List, blockName_String:"MASS", streamName_String:"mass"] :=
    CreateSLHAMassBlockStream[massMatrices, blockName, streamName] <>
    "slha_io.set_block(" <> streamName <> ");\n";

ConvertToRealInputParameter[parameter_, struct_String] :=
    struct <> CConversion`ToValidCSymbolString[parameter];

ConvertToRealInputParameter[FlexibleSUSY`Phase[parameter_], struct_String] :=
    "Re(" <> struct <> CConversion`ToValidCSymbolString[FlexibleSUSY`Phase[parameter]] <> ")";

WriteParameterTuple[{key_?NumberQ, parameter_}, streamName_String] :=
    Module[{parameterStr},
           parameterStr = CConversion`ToValidCSymbolString[parameter];
           streamName <> " << FORMAT_ELEMENT(" <> ToString[key] <> ", " <>
           ConvertToRealInputParameter[parameter,"input."] <>
           ", \"" <> parameterStr <> "\");\n"
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

WriteSLHAMatrix[{mixingMatrix_, lesHouchesName_}, head_String, scale_String, setter_String:"set_block"] :=
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
           "slha_io." <> setter <> "(\"" <> lhs <> "\", " <> wrapper <> ", \"" <> str <>
           "\"" <> If[scale != "", ", " <> scale, ""] <> ");\n"
          ];

WriteSLHAMixingMatricesBlocks[] :=
    Module[{result, mixingMatrices, smMix, susyMix, majoranaMix,
            smMixStr = "", susyMixStr = "", majoranaMixStr = ""},
           mixingMatrices = GetSLHAMixinMatrices[];
           smMix = Flatten[TreeMasses`FindMixingMatrixSymbolFor /@ SARAH`SMParticles];
           smMix = Select[mixingMatrices, MemberQ[smMix,#[[1]]]&];
           majoranaMix = Flatten[TreeMasses`FindMixingMatrixSymbolFor /@
                                 Select[TreeMasses`GetParticles[], TreeMasses`IsMajoranaFermion]];
           majoranaMix = Select[mixingMatrices, MemberQ[majoranaMix,#[[1]]]&];
           majoranaMix = {#[[1]], Symbol["IM" <> ToString[#[[2]]]]}& /@ majoranaMix;
           susyMix = Complement[mixingMatrices, smMix];
           (smMixStr = smMixStr <> WriteSLHAMatrix[#,"LOCALPHYSICAL",""])& /@ smMix;
           (susyMixStr = susyMixStr <> WriteSLHAMatrix[#,"LOCALPHYSICAL",""])& /@ susyMix;
           (majoranaMixStr = majoranaMixStr <> WriteSLHAMatrix[#,"LOCALPHYSICAL","","set_block_imag"])& /@ majoranaMix;
           result = susyMixStr <> "\n" <>
                    "if (write_sm_mixing_matrics) {\n" <>
                    TextFormatting`IndentText[smMixStr] <> "}\n\n" <>
                    "if (print_imaginary_parts_of_majorana_mixings) {\n" <>
                    TextFormatting`IndentText[majoranaMixStr] <> "}\n";
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

WriteEffectiveCouplingsSLHABlockEntry[particle_, vectorBoson_] :=
    Module[{i, dim, dimWithoutGoldstones, start, particlePDG, vectorPDG,
            struct, comment, value, result = ""},
           vectorPDG = Parameters`GetPDGCodesForParticle[vectorBoson][[1]];
           particlePDG = Parameters`GetPDGCodesForParticle[particle];
           dim = TreeMasses`GetDimension[particle];
           dimWithoutGoldstones = TreeMasses`GetDimensionWithoutGoldstones[particle];
           If[Length[particlePDG] != dim,
              Print["Warning: length of PDG number list != dimension of particle ", particle];
              Print["       PDG number list = ", particlePDG];
              Print["       dimension of particle ", particle, " = ", dim];
             ];
           If[Length[particlePDG] < dim,
              Return[""];
             ];
           start = TreeMasses`GetDimensionStartSkippingGoldstones[particle];
           Which[particle === SARAH`HiggsBoson && vectorBoson === SARAH`VectorP,
                 struct = "OBSERVABLES.eff_cp_higgs_photon_photon";
                 comment = "Abs(effective H-Photon-Photon coupling)";,
                 particle === SARAH`HiggsBoson && vectorBoson === SARAH`VectorG,
                 struct = "OBSERVABLES.eff_cp_higgs_gluon_gluon";
                 comment = "Abs(effective H-Gluon-Gluon coupling)";,
                 particle === SARAH`PseudoScalar && vectorBoson === SARAH`VectorP,
                 struct = "OBSERVABLES.eff_cp_pseudoscalar_photon_photon";
                 comment = "Abs(effective A-Photon-Photon coupling)";,
                 particle === SARAH`PseudoScalar && vectorBoson === SARAH`VectorG,
                 struct = "OBSERVABLES.eff_cp_pseudoscalar_gluon_gluon";
                 comment = "Abs(effective A-Gluon-Gluon coupling)";,
                 True,
                 Print["Error: unsupported effective coupling ",
                       particle, "-", vectorBoson, "-", vectorBoson,
                       "requested!"];
                 Quit[1]
                ];
           If[dimWithoutGoldstones == 1 || start == dim,
              value = "Abs(" <> struct <> ")";
              result = result
                        <> WriteSLHABlockEntry[{value, particlePDG[[start]], vectorPDG, vectorPDG},
                                                comment];,
              For[i = start, i <= Length[particlePDG], i++,
                  value = "Abs(" <> struct <> "(" <> ToString[i-start] <> "))";
                  result = result
                           <> WriteSLHABlockEntry[{value, particlePDG[[i]], vectorPDG, vectorPDG},
                                                  comment];
                 ];
             ];
           result
          ];

WriteSLHABlockEntry[{par_ /; MemberQ[FlexibleSUSYObservable`FSObservables,par], idx___}, comment_String:""] :=
    Module[{i, dim, scalarPDG, vectorPDG, result = ""},
           Switch[par,
                  FlexibleSUSYObservable`aMuonGM2Calc,
                      result = WriteSLHABlockEntry[{"OBSERVABLES.a_muon_gm2calc", idx}, "Delta(g-2)_muon/2 GM2Calc"],
                  FlexibleSUSYObservable`aMuonGM2CalcUncertainty,
                      result = WriteSLHABlockEntry[{"OBSERVABLES.a_muon_gm2calc_uncertainty", idx}, "Delta(g-2)_muon/2 GM2Calc uncertainty"],
                  FlexibleSUSYObservable`CpHiggsPhotonPhoton,
                      result = WriteEffectiveCouplingsSLHABlockEntry[SARAH`HiggsBoson, SARAH`VectorP],
                  FlexibleSUSYObservable`CpHiggsGluonGluon,
                      result = WriteEffectiveCouplingsSLHABlockEntry[SARAH`HiggsBoson, SARAH`VectorG],
                  FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton,
                      result = WriteEffectiveCouplingsSLHABlockEntry[SARAH`PseudoScalar, SARAH`VectorP],
                  FlexibleSUSYObservable`CpPseudoScalarGluonGluon,
                      result = WriteEffectiveCouplingsSLHABlockEntry[SARAH`PseudoScalar, SARAH`VectorG],
                  _,
                     result = WriteSLHABlockEntry[{"", idx}, ""]
                 ];
           result
          ];

WriteSLHABlockEntry[{par_, idx1_?NumberQ, idx2_?NumberQ, idx3_?NumberQ}, comment_String:""] :=
    Module[{parStr, parVal, idx1Str, idx2Str, idx3Str, commentStr},
           parStr = CConversion`RValueToCFormString[Parameters`IncreaseIndexLiterals[par]];
           parVal = CConversion`RValueToCFormString[
               WrapPreprocessorMacroAround[par, Join[Parameters`GetModelParameters[],
                                                     Parameters`GetOutputParameters[]],
                                           Global`MODELPARAMETER]];
           idx1Str = ToString[idx1];
           idx2Str = ToString[idx2];
           idx3Str = ToString[idx3];
           commentStr = If[comment == "", parStr, comment];
           (* result *)
           "      << FORMAT_RANK_THREE_TENSOR(" <> idx1Str <> ", " <> idx2Str <> ", "
           <> idx3Str <> ", (" <> parVal <> "), \"" <> commentStr <> "\")" <> "\n"
          ];

WriteSLHABlockEntry[{par_, idx1_?NumberQ, idx2_?NumberQ}, comment_String:""] :=
    Module[{parStr, parVal, idx1Str, idx2Str, commentStr},
           parStr = CConversion`RValueToCFormString[Parameters`IncreaseIndexLiterals[par]];
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
           parStr = CConversion`RValueToCFormString[Parameters`IncreaseIndexLiterals[par]];
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
           parStr = CConversion`RValueToCFormString[Parameters`IncreaseIndexLiterals[par]];
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

WriteSLHABlock[{blockName_, tuples_List}, scale_String:"model.get_scale()"] :=
    Module[{result = "", blockNameStr},
           blockNameStr = ToString[blockName];
           result = "std::ostringstream block;\n" <>
                    "block << \"Block " <> blockNameStr <>
                    If[scale != "",
                       " Q= \" << FORMAT_SCALE(" <> scale <> ")",
                       "\""
                      ] <>
                    " << '\\n'\n";
           (result = result <> WriteSLHABlockEntry[#])& /@ tuples;
           result = result <> ";\n" <> "slha_io.set_block(block);\n";
           result = "{\n" <> TextFormatting`IndentText[result] <> "}\n";
           Return[result];
          ];

WriteSLHABlock[{blockName_, parameter_}, scale_String:"model.get_scale()"] :=
    WriteSLHAMatrix[{parameter, blockName}, "MODELPARAMETER", scale];

WriteSLHABlock[{blockName_, {parameter_ /; Head[parameter] =!= List}}, scale_String:"model.get_scale()"] :=
    WriteSLHABlock[{blockName, parameter}, scale];

WriteSLHAModelParametersBlocks[] :=
    Module[{result = "", modelParameters, blocks},
           modelParameters = GetSLHAModelParameters[];
           blocks = SortBlocks[modelParameters];
           (result = result <> WriteSLHABlock[#])& /@ blocks;
           Return[result];
          ];

GetExtraSLHAOutputBlockScale[scale_ /; scale === FlexibleSUSY`NoScale] := "";

GetExtraSLHAOutputBlockScale[scale_ /; scale === FlexibleSUSY`CurrentScale] := "model.get_scale()";

GetExtraSLHAOutputBlockScale[scale_?NumericQ] := ToString[scale];

GetExtraSLHAOutputBlockScale[scale_] :=
    Module[{result},
           scaleStr = CConversion`RValueToCFormString[
               WrapPreprocessorMacroAround[Parameters`DecreaseIndexLiterals[scale],
                                           Join[Parameters`GetModelParameters[],
                                                Parameters`GetOutputParameters[]],
                                           Global`MODELPARAMETER]];
           StringReplace[scaleStr, "CurrentScale" -> "model.get_scale()"]
          ];

WriteExtraSLHAOutputBlock[outputBlocks_List] :=
    Module[{result = "", reformed},
           ReformeBlocks[{block_, tuples_List}] := {{block, ReformeBlocks /@ tuples}, "model.get_scale()"};
           ReformeBlocks[{block_, scale_, tuples_List}] := {{block, ReformeBlocks /@ tuples},
                                                            GetExtraSLHAOutputBlockScale[scale]};
           ReformeBlocks[{expr_}]               := {expr};
           ReformeBlocks[{idx_, expr_}]         := {expr, idx};
           ReformeBlocks[{idx1_, idx2_, expr_}] := {expr, idx1, idx2};
           reformed = ReformeBlocks /@ outputBlocks;
           (result = result <> WriteSLHABlock[#[[1]], #[[2]]])& /@ reformed;
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

ReadLesHouchesInputParameters[slhaInputParameters_List] :=
    Module[{result = ""},
           (result = result <> ReadSLHAInputBlock[#])& /@ slhaInputParameters;
           Return[result];
          ];

ReadSLHAOutputBlock[{parameter_, {blockName_String, pdg_?NumberQ}}] :=
    Module[{result, parmStr, pdgStr, gutNorm = ""},
           parmStr = CConversion`ToValidCSymbolString[parameter];
           pdgStr = ToString[pdg];
           If[parameter === SARAH`hyperchargeCoupling,
              gutNorm = " * " <> CConversion`RValueToCFormString[
                  1/Parameters`GetGUTNormalization[parameter]];
             ];
           result = "model.set_" <> parmStr <>
                    "(slha_io.read_entry(\"" <> blockName <> "\", " <>
                    pdgStr <> ")" <> gutNorm <> ");\n";
           Return[result];
          ];

ReadSLHAOutputBlock[{parameter_, {blockName_Symbol, pdg_?NumberQ}}] :=
    ReadSLHAOutputBlock[{parameter, {ToString[blockName], pdg}}];

ReadSLHAOutputBlock[{parameter_, {blockName_, pdg_?NumberQ}}] :=
    Block[{},
          Print["Warning: SLHA block name is not a symbol: ", blockName];
          Print["   I'm using: ", CConversion`RValueToCFormString[blockName]];
          ReadSLHAOutputBlock[{parameter, {CConversion`RValueToCFormString[blockName], pdg}}]
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

ReadSLHAPhysicalMixingMatrixBlock[{parameter_, blockName_Symbol}, struct_String:"PHYSICAL", defMacro_String:"DEFINE_PHYSICAL_PARAMETER"] :=
    Module[{paramStr, blockNameStr},
           paramStr = CConversion`ToValidCSymbolString[parameter];
           blockNameStr = ToString[blockName];
           "{\n" <> IndentText[
               defMacro <> "(" <> paramStr <> ");\n" <>
               "slha_io.read_block(\"" <> blockNameStr <> "\", " <>
               paramStr <> ");\n" <>
               struct <> "(" <> paramStr <> ") = " <> paramStr <> ";"] <> "\n" <>
           "}\n"
          ];

ReadSLHAPhysicalMass[particle_,struct_String:"PHYSICAL"] :=
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
                 result = struct <> "(" <> massStr <>
                          ") = slha_io.read_entry(\"MASS\", " <> pdgStr <> ");\n";
                ];
              ,
              For[i = 1, i <= dim, i++,
                  pdg = Abs[pdgList[[i]]];
                  pdgStr = ToString[pdg];
                  If[pdg != 0,
                     result = result <>
                              struct <> "(" <> massStr <> ")(" <> ToString[i-1] <>
                              ") = slha_io.read_entry(\"MASS\", " <> pdgStr <> ");\n";
                    ];
                 ];
             ];
           Return[result];
          ];

ReadSLHAPhysicalMassBlock[struct_String:"PHYSICAL"] :=
    Module[{result = "", particles},
           particles = TreeMasses`GetParticles[];
           (result = result <> ReadSLHAPhysicalMass[#,struct])& /@ particles;
           Return[result];
          ];

ReadLesHouchesOutputParameters[] :=
    Module[{result = "", modelParameters},
           modelParameters = GetSLHAModelParameters[];
           (result = result <> ReadSLHAOutputBlock[#])& /@ modelParameters;
           Return[result];
          ];

ReadLesHouchesPhysicalParameters[struct_String:"PHYSICAL", defMacro_String:"DEFINE_PHYSICAL_PARAMETER"] :=
    Module[{result = "", physicalParameters},
           physicalParameters = GetSLHAMixinMatrices[];
           (result = result <> ReadSLHAPhysicalMixingMatrixBlock[#,struct,defMacro])& /@ physicalParameters;
           result = result <> "\n" <> ReadSLHAPhysicalMassBlock[struct];
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
    ConvertMixingsToConvention[massMatrices, "slha"];

ConvertMixingsToHKConvention[massMatrices_List] :=
    ConvertMixingsToConvention[massMatrices, "hk"];

ConvertMixingsToConvention[massMatrices_List, convention_String] :=
    Module[{result = "", i,
            eigenstateName, mixingMatrixSym,
            eigenstateNameStr, mixingMatrixSymStr},
           For[i = 1, i <= Length[massMatrices], i++,
               eigenstateName = TreeMasses`GetMassEigenstate[massMatrices[[i]]];
               mixingMatrixSym = TreeMasses`GetMixingMatrixSymbol[massMatrices[[i]]];
               If[IsMajoranaFermion[eigenstateName] && mixingMatrixSym =!= Null,
                  eigenstateNameStr  = CConversion`ToValidCSymbolString[FlexibleSUSY`M[eigenstateName]];
                  mixingMatrixSymStr = CConversion`ToValidCSymbolString[mixingMatrixSym];
                  result = result <>
                           "SLHA_io::convert_symmetric_fermion_mixings_to_" <> convention <> "(LOCALPHYSICAL(" <>
                           eigenstateNameStr <> "), LOCALPHYSICAL(" <>
                           mixingMatrixSymStr <> "));\n";
                 ];
              ];
           Return[result];
          ];

ParseCmdLineOption[{parameter_, CConversion`ScalarType[CConversion`realScalarCType | CConversion`integerScalarCType]}] :=
    Module[{parameterStr},
           parameterStr = CConversion`ToValidCSymbolString[parameter];
           "\
if(Command_line_options::get_parameter_value(option, \"--" <> parameterStr <> "=\", input." <> parameterStr <>"))
   continue;

"
          ];

ParseCmdLineOption[_] := "";

ParseCmdLineOptions[inputParameters_List] :=
    StringJoin[ParseCmdLineOption /@ inputParameters];

PrintCmdLineOption[{parameter_, CConversion`ScalarType[CConversion`realScalarCType | CConversion`integerScalarCType]}] :=
    "\"  --" <> ToValidCSymbolString[parameter] <> "=<value>\\n\"\n";

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
    CConversion`ToRealType[Parameters`GetType[c]];

CreateSLHASoftSquaredMassName[c_] :=
    CConversion`ToValidCSymbolString[c] <> "_slha";

GetSLHASoftSquaredMassType[c_] :=
    CConversion`ToRealType[Parameters`GetType[c]];

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

ParametersHaveSameDimension[pars_List] :=
    CConversion`HaveSameDimension[Parameters`GetType /@ pars];

ConvertYukawaCouplingsToSLHA[] :=
    Module[{result = ""},
           yuks = GetYukawas[];
           Module[{dim, vL, vR},
                  dim = SARAH`getDimParameters[#][[1]];
                  {vL, vR} = GetMixingMatricesFor[#];
                  If[Parameters`IsOutputParameter[{vL, vR}] &&
                     ParametersHaveSameDimension[{vL, vR, #}],
                     result = result <>
                              "fs_svd(" <> CConversion`ToValidCSymbolString[#] <> ", " <>
                                      CreateSLHAYukawaName[#] <> ", " <>
                                      CreateSLHAFermionMixingMatrixName[vR] <> ", " <>
                                      CreateSLHAFermionMixingMatrixName[vL] <> ");\n";
                     ,
                     Print["Warning: Cannot convert Yukawa coupling ", #,
                           " to SLHA, because ", {vL,vR}, " are not defined",
                           " or have incompatible dimension."];
                     result = result <>
                              CreateSLHAYukawaName[#] <> " = " <>
                              CConversion`ToValidCSymbolString[#] <> ".diagonal().real();\n";
                    ];
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
                  If[Parameters`IsOutputParameter[{vL, vR}] &&
                     ParametersHaveSameDimension[{vL, vR, #}],
                     result = result <>
                              CreateSLHATrilinearCouplingName[#] <> " = (" <>
                              CreateSLHAFermionMixingMatrixName[vR] <> ".conjugate() * " <>
                              CConversion`ToValidCSymbolString[#] <> " * " <>
                              CreateSLHAFermionMixingMatrixName[vL] <> ".adjoint()" <>
                              ").real();\n";
                     ,
                     Print["Warning: Cannot convert Trilinear coupling ", #,
                           " to SLHA, because ", {vL,vR}, " are not defined",
                           " or have incompatible dimension."];
                     result = result <>
                              CreateSLHATrilinearCouplingName[#] <> " = " <>
                              CConversion`ToValidCSymbolString[#] <> ".real();\n";
                    ];
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
                  If[Parameters`IsOutputParameter[{vL, vR}] &&
                     ParametersHaveSameDimension[{vL, vR, #}],
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
                     ,
                     Print["Warning: Cannot convert soft squared mass ", #,
                           " to SLHA, because ", {vL,vR}, " are not defined",
                           " or have incompatible dimension."];
                     result = result <>
                              CreateSLHASoftSquaredMassName[#] <> " = " <>
                              CConversion`ToValidCSymbolString[#] <> ".real();\n";
                    ];
           ]& /@ massSq;
           result
          ];

CalculateCKMMatrix[] :=
    Module[{result = ""},
           If[MemberQ[Parameters`GetOutputParameters[], SARAH`DownMatrixL] &&
              MemberQ[Parameters`GetOutputParameters[], SARAH`UpMatrixL] &&
              SARAH`getDimParameters[SARAH`DownMatrixL] === {3,3} &&
              SARAH`getDimParameters[SARAH`UpMatrixL] === {3,3}
              ,
              result = result <> "ckm = " <>
              CreateSLHAFermionMixingMatrixName[SARAH`UpMatrixL  ] <> " * " <>
              CreateSLHAFermionMixingMatrixName[SARAH`DownMatrixL] <> ".adjoint();\n";
             ];
           (* convert CKM matrix to PDG convention *)
           If[MemberQ[Parameters`GetOutputParameters[], SARAH`DownMatrixL] &&
              MemberQ[Parameters`GetOutputParameters[], SARAH`UpMatrixL  ] &&
              MemberQ[Parameters`GetOutputParameters[], SARAH`DownMatrixR] &&
              MemberQ[Parameters`GetOutputParameters[], SARAH`UpMatrixR  ] &&
              SARAH`getDimParameters[SARAH`DownMatrixL] === {3,3} &&
              SARAH`getDimParameters[SARAH`UpMatrixL] === {3,3} &&
              SARAH`getDimParameters[SARAH`DownMatrixR] === {3,3} &&
              SARAH`getDimParameters[SARAH`UpMatrixR] === {3,3}
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
              SARAH`getDimParameters[SARAH`ElectronMatrixL] === {3,3} &&
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

CreateFormattedSLHABlockEntry[{par_, CConversion`ScalarType[_]}] :=
    "   0   # " <> CConversion`RValueToCFormString[par] <> "\n";

CreateFormattedSLHABlockEntry[{par_, (CConversion`ArrayType | CConversion`VectorType)[_,n_]}] :=
    Module[{i, result = ""},
           For[i = 1, i <= n, i++,
               result = result <>
                        "   " <> ToString[i] <> "   0   # " <>
                        CConversion`RValueToCFormString[par[i]] <> "\n";
              ];
           result
          ];

CreateFormattedSLHABlockEntry[{par_, CConversion`MatrixType[_,m_,n_]}] :=
    Module[{i, k, result = ""},
           For[i = 1, i <= m, i++,
               For[k = 1, k <= n, k++,
                   result = result <>
                            "   " <> ToString[i] <> "   " <> ToString[k] <>
                            "   0   # " <>
                            CConversion`RValueToCFormString[par[i,k]] <> "\n";
                  ];
              ];
           result
          ];

CreateFormattedSLHABlockEntry[{par_, CConversion`TensorType[_,m_,n_,o_]}] :=
    Module[{i, k, l, result = ""},
           For[i = 1, i <= m, i++,
               For[k = 1, k <= n, k++,
                   For[l = 1, l <= o, l++,
                       result = result <>
                                "   " <> ToString[i] <> "   " <> ToString[k] <>
                                "   " <> ToString[l] <> "   0   # " <>
                                CConversion`RValueToCFormString[par[i,k,l]] <> "\n";
                      ];
                  ];
              ];
           result
          ];

CreateFormattedSLHABlockEntry[{par_, CConversion`TensorType[_,__]}] := "\n";

CreateFormattedSLHABlockEntry[{par_, _, idx_}] :=
    "   " <> ToString[idx] <> "   0   # " <> CConversion`RValueToCFormString[par] <> "\n";

CreateFormattedSLHABlock[{block_, parameters_List}] :=
    Module[{head = "Block " <> ToString[block] <> "\n", body},
           body = StringJoin[CreateFormattedSLHABlockEntry /@ parameters];
           head <> body
          ];

FindParametersInBlock[inputParameters_List, block_] :=
    Join[Cases[inputParameters, {p_, block, t_} :> {p, t}],
         Cases[inputParameters, {p_, {block, idx_}, t_} :> {p, t, idx}]];

CreateFormattedSLHABlocks[inputPars_List] :=
    Module[{blocks, sortForBlocks},
           blocks = DeleteDuplicates @ Cases[inputPars, {_, {block_, _} | block_, ___} :> block];
           sortForBlocks = {#, FindParametersInBlock[inputPars, #]}& /@ blocks;
           StringJoin[CreateFormattedSLHABlock /@ sortForBlocks]
          ];

End[];

EndPackage[];


BeginPackage["WriteOut`", {"SARAH`", "TextFormatting`", "CConversion`", "Parameters`", "TreeMasses`"}];

ReplaceInFiles::usage="Replaces tokens in files.";
PrintParameters::usage="Creates parameter printout statements";
WriteSLHAMassBlock::usage="";
WriteSLHAMixingMatricesBlocks::usage="";
WriteSLHAModelParametersBlocks::usage="";
ReadUnfixedParameters::usage="";

Begin["Private`"];

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

PrintParameter[Null, streamName_String] := "";

PrintParameter[parameter_, streamName_String] :=
    Module[{parameterName},
           parameterName = CConversion`ToValidCSymbolString[parameter /. a_[i1,i2] :> a];
           Return[streamName <> " << \"" <> parameterName <> " = \" << " <>
                  parameterName <> " << '\\n';\n"];
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
                          ", " <> massNameStr <> ", \"" <> eigenstateNameStr <> "\")\n";
                ];
              ,
              For[i = 1, i <= dim, i++,
                  pdg = Abs[pdgList[[i]]];
                  If[pdg != 0,
                     eigenstateNameStr = CConversion`RValueToCFormString[eigenstateName] <> "_" <> ToString[i];
                     massNameStr = CConversion`RValueToCFormString[FlexibleSUSY`M[eigenstateName[i]]];
                     result = result <> "<< FORMAT_MASS(" <> ToString[pdg] <>
                              ", " <> massNameStr <> ", \"" <> eigenstateNameStr <> "\")\n";
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
           smMassesStr = "if (model.do_calculate_sm_pole_masses()) {\n" <>
                         TextFormatting`IndentText["mass\n" <>
                             TextFormatting`IndentText[smMassesStr] <> ";"] <>
                         "\n}\n\n";
           result = Parameters`CreateLocalConstRefsForPhysicalParameters[allMasses] <> "\n" <>
                    "std::ostringstream mass;\n\n" <>
                    susyMassesStr <> smMassesStr <>
                    "slha_io.set_block(mass);\n";
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
    Module[{str, lhs, wrapper},
           If[SARAH`getDimParameters[mixingMatrix] === {} ||
              SARAH`getDimParameters[mixingMatrix] === {1},
              Print["Warning: You are trying to create a SLHA matrix block for"];
              Print["   ", mixingMatrix, ", which is not a matrix!"];
              Print["   Please specify a Les Houches index in the SARAH model file."];
             ];
           str = CConversion`ToValidCSymbolString[mixingMatrix];
           lhs = ToString[lesHouchesName];
           wrapper = If[head == "", str, head <> "(" <> str <> ")"];
           "slha_io.set_block(\"" <> lhs <> "\", " <> wrapper <> ", \"" <> str <>
           "\"" <> If[scale != "", ", " <> scale, ""] <> ");\n"
          ];

WriteSLHAMixingMatricesBlocks[] :=
    Module[{result, mixingMatrices, smMix, susyMix, smMixStr = "", susyMixStr = ""},
           mixingMatrices = GetSLHAMixinMatrices[];
           smMix = Flatten[TreeMasses`FindMixingMatrixSymbolFor /@ SARAH`SMParticles];
           smMix = Select[mixingMatrices, MemberQ[smMix,#[[1]]]&];
           susyMix = Complement[mixingMatrices, smMix];
           (smMixStr = smMixStr <> WriteSLHAMatrix[#,"PHYSICAL"])& /@ smMix;
           (susyMixStr = susyMixStr <> WriteSLHAMatrix[#,"PHYSICAL"])& /@ susyMix;
           result = susyMixStr <> "\n" <>
                    "if (model.do_calculate_sm_pole_masses()) {\n" <>
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
           collected = collected /. {a_Symbol} :> a
          ];

WriteSLHABlock[{blockName_, tuples_List}] :=
    Module[{result = "", blockNameStr, t, pdg, parmStr, scale},
           blockNameStr = ToString[blockName];
           scale = "model.get_scale()";
           result = "std::ostringstream block;\n" <>
                    "block << \"Block " <> blockNameStr <> " Q= \" << FORMAT_NUMBER(" <>
                    scale <> ") << '\\n'\n";
           For[t = 1, t <= Length[tuples], t++,
               parmStr = CConversion`ToValidCSymbolString[tuples[[t,1]]];
               pdg = ToString[tuples[[t,2]]];
               result = result <> "      << FORMAT_ELEMENT(" <> pdg <> ", " <>
                        "MODELPARAMETER(" <> parmStr <> "), \"" <> parmStr <> "\")" <>
                        If[t == Length[tuples], ";", ""] <> "\n";
              ];
           result = result <> "slha_io.set_block(block);\n";
           result = "{\n" <> TextFormatting`IndentText[result] <> "}\n";
           Return[result];
          ];

WriteSLHABlock[{blockName_, parameter_}] :=
    WriteSLHAMatrix[{parameter, blockName}, "MODELPARAMETER", "model.get_scale()"];

WriteSLHAModelParametersBlocks[] :=
    Module[{result = "", modelParameters, blocks},
           modelParameters = GetSLHAModelParameters[];
           blocks = SortBlocks[modelParameters];
           (result = result <> WriteSLHABlock[#])& /@ blocks;
           Return[result];
          ];

ReadSLHABlock[{parameter_, {blockName_Symbol, pdg_?NumberQ}}] :=
    Module[{result, blockNameStr, parmStr, pdgStr},
           blockNameStr = ToString[blockName] <> "IN";
           parmStr = CConversion`ToValidCSymbolString[parameter];
           pdgStr = ToString[pdg];
           result = "input." <> parmStr <>
                    " = slha_io.read_entry(\"" <> blockNameStr <> "\", " <>
                    pdgStr <> ");\n";
           Return[result];
          ];

ReadSLHABlock[{parameter_, blockName_Symbol}] :=
    Module[{paramStr, blockNameStr},
           paramStr = CConversion`ToValidCSymbolString[parameter];
           blockNameStr = ToString[blockName] <> "IN";
           "slha_io.read_block(\"" <> blockNameStr <> "\", input." <>
           paramStr <> ");\n"
          ];

ReadUnfixedParameters[unfixedParameters_List] :=
    Module[{result = "", modelParameters, unfixedNames, rules},
           unfixedNames = (#[[1]])& /@ unfixedParameters;
           rules = Rule[#[[1]], #[[2]]]& /@ unfixedParameters;
           (* get block names of all unfixed parameters (unfixedNames) *)
           modelParameters = Select[GetSLHAModelParameters[], MemberQ[unfixedNames,#[[1]]]&];
           modelParameters = {#[[1]] /. rules, #[[2]]}& /@ modelParameters;
           (result = result <> ReadSLHABlock[#])& /@ modelParameters;
           Return[result];
          ];

End[];

EndPackage[];

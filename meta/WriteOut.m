
BeginPackage["WriteOut`", {"SARAH`", "TextFormatting`", "CConversion`", "Parameters`", "TreeMasses`"}];

ReplaceInFiles::usage="Replaces tokens in files.";
PrintParameters::usage="Creates parameter printout statements";
WriteSLHAExtparBlock::usage="";
WriteSLHAMassBlock::usage="";
WriteSLHAMixingMatricesBlocks::usage="";
WriteSLHAModelParametersBlocks::usage="";
WriteSLHAMinparBlock::usage="";
ReadLesHouchesInputParameters::usage="";
StringJoinWithSeparator::usage="Joins a list of strings with a given separator string";

Begin["`Private`"];

StringJoinWithSeparator[strings_List, separator_String] :=
    Module[{result = "", i},
           For[i = 1, i <= Length[strings], i++,
               If[i > 1, result = result <> separator;];
               result = result <> ToString[strings[[i]]];
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
                          ", " <> massNameStr <> ", \"" <> eigenstateNameStr <> "\")\n";
                ];
              ,
              For[i = 1, i <= dim, i++,
                  pdg = Abs[pdgList[[i]]];
                  If[pdg != 0,
                     eigenstateNameStr = CConversion`RValueToCFormString[eigenstateName] <> "_" <> ToString[i];
                     massNameStr = CConversion`RValueToCFormString[FlexibleSUSY`M[eigenstateName[i-1]]];
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

ReadLesHouchesInputParameters[lesHouchesInputParameters_List] :=
    Module[{result = "", modelParameters, names, rules},
           names = (#[[1]])& /@ lesHouchesInputParameters;
           rules = Rule[#[[1]], #[[2]]]& /@ lesHouchesInputParameters;
           (* get block names of all les Houches input parameters (names) *)
           modelParameters = Select[GetSLHAModelParameters[], MemberQ[names,#[[1]]]&];
           modelParameters = {#[[1]] /. rules, #[[2]]}& /@ modelParameters;
           (result = result <> ReadSLHABlock[#])& /@ modelParameters;
           Return[result];
          ];

End[];

EndPackage[];

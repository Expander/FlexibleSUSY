
BeginPackage["WriteOut`", {"SARAH`", "TextFormatting`", "CConversion`", "Parameters`", "TreeMasses`"}];

PrintParameters::usage="Creates parameter printout statements";
WriteSLHAMassBlock::usage="";
WriteSLHAMixingMatricesBlocks::usage="";
WriteSLHAModelParametersBlocks::usage="";

Begin["Private`"];

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
              Print["Error: length of PDG list != dimension of particle"];
              Print["       PDG list = ", pdgList];
              Quit[1];
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
    Module[{result = "", allMasses},
           allMasses = FlexibleSUSY`M[TreeMasses`GetMassEigenstate[#]]& /@ massMatrices;
           (result = result <> WriteSLHAMass[#])& /@ massMatrices;
           result = Parameters`CreateLocalConstRefsForPhysicalParameters[allMasses] <> "\n" <>
                    "std::ostringstream mass;\n\n" <>
                    "mass << \"Block MASS\\n\"\n" <>
                    TextFormatting`IndentText[result] <> ";\n\n" <>
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
    Module[{str, lhs, wrapper},
           str = CConversion`ToValidCSymbolString[mixingMatrix];
           lhs = ToString[lesHouchesName];
           wrapper = If[head == "", str, head <> "(" <> str <> ")"];
           "set_block(\"" <> lhs <> "\", " <> wrapper <> ", \"" <> str <> "\");\n"
          ];

WriteSLHAMatrix[{mixingMatrix_, lesHouchesName_}, head_String, scale_String] :=
    Module[{str, lhs, wrapper},
           str = CConversion`ToValidCSymbolString[mixingMatrix];
           lhs = ToString[lesHouchesName];
           wrapper = If[head == "", str, head <> "(" <> str <> ")"];
           "set_block(\"" <> lhs <> "\", " <> wrapper <> ", \"" <> str <>
           "\", " <> scale <> ");\n"
          ];

WriteSLHAMixingMatricesBlocks[] :=
    Module[{result = "", mixingMatrices},
           mixingMatrices = GetSLHAMixinMatrices[];
           (result = result <> WriteSLHAMatrix[#,"PHYSICAL"])& /@ mixingMatrices;
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
           collected = collected /. {a_} :> a
          ];

WriteSLHABlock[{blockName_, tuples_List}] :=
    Module[{result = "", blockNameStr, t, pdg, parmStr},
           blockNameStr = ToString[blockName];
           result = "std::ostringstream block;\n" <>
                    "block << \"Block " <> blockNameStr <> "\\n\"\n";
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

End[];

EndPackage[];


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
                          "," <> massNameStr <> ",\"" <> eigenstateNameStr <> "\")\n";
                ];
              ,
              For[i = 1, i <= dim, i++,
                  pdg = Abs[pdgList[[i]]];
                  If[pdg != 0,
                     eigenstateNameStr = CConversion`RValueToCFormString[eigenstateName] <> "_" <> ToString[i];
                     massNameStr = CConversion`RValueToCFormString[FlexibleSUSY`M[eigenstateName[i]]];
                     result = result <> "<< FORMAT_MASS(" <> ToString[pdg] <>
                              "," <> massNameStr <> ",\"" <> eigenstateNameStr <> "\")\n";
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
    Select[FlexibleSUSY`FSLesHouchesList,
           MemberQ[Parameters`GetModelParameters[],#[[1]]]&];

WriteSLHAMixingMatrix[{mixingMatrix_, lesHouchesName_}] :=
    Module[{str, lhs},
           str = CConversion`ToValidCSymbolString[mixingMatrix];
           lhs = ToString[lesHouchesName];
           "set_block(\"" <> lhs <> "\", PHYSICAL(" <> str <> "), \"" <> str <> "\");\n"
          ];

WriteSLHAMixingMatricesBlocks[] :=
    Module[{result = "", mixingMatrices},
           mixingMatrices = GetSLHAMixinMatrices[];
           (result = result <> WriteSLHAMixingMatrix[#])& /@ mixingMatrices;
           result = Parameters`CreateLocalConstRefsForPhysicalParameters[(#[[1]])& /@ mixingMatrices] <>
                    "\n" <> result;
           Return[result];
          ];

WriteSLHAModelParametersBlocks[] :=
    Module[{result = "", modelParameters},
           modelParameters = GetSLHAModelParameters[];
           Print["modelParameters = ", modelParameters];
           Return[result];
          ];

End[];

EndPackage[];

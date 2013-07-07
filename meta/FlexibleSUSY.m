
BeginPackage["FlexibleSUSY`", {"SARAH`", "AnomalousDimension`", "BetaFunction`", "TextFormatting`", "CConversion`", "TreeMasses`", "EWSB`", "Traces`", "SelfEnergies`", "Phases`", "LoopMasses`", "WriteOut`", "Constraint`", "ThresholdCorrections`", "ConvergenceTester`"}];

MakeFlexibleSUSY::usage="";

LowPrecision::usage="";
MediumPrecision::usage="";
HighPrecision::usage="";
softSusyCompatibleRGEs::usage="";

Begin["Private`"];

allParameters = {};

GetParameters[] := allParameters;

allIndexReplacementRules = {};

GetIndexReplacementRules[] := allIndexReplacementRules;

allBetaFunctions = {};

GetBetaFunctions[] := allBetaFunctions;

allOutputParameters = {};

GetOutputParameters[] := allOutputParameters;

CheckModelFileSettings[] :=
    Module[{},
           If[Head[Global`InitialGuess] =!= List,
              Global`InitialGuess = {};
             ];
           If[!NameQ["Global`BoundaryHighScaleFirstGuess"],
              Print["Warning: Global`BoundaryHighScaleFirstGuess should be",
                    " set in the model file!"];
              Global`BoundaryHighScaleFirstGuess = 1.0 10^14;
             ];
           If[!NameQ[Global`BoundaryLowScale],
              Print["Warning: Global`BoundaryLowScale should be",
                    " set in the model file!"];
              Global`BoundaryLowScale = Global`MZ;
             ];
           If[Head[Global`DefaultParameterPoint] =!= List,
              Global`DefaultParameterPoint = {};
             ];
           If[Head[Global`InputParameters] =!= List,
              Global`InputParameters = {};
             ];
          ];

(*
   @brief Replaces tokens in files.

   @param files list of two-element lists.  The first entry is the
   input file and the second entry is the output file.
   Example:
      files = {{"input.hpp", "output.hpp"},
               {"input.cpp", "output.cpp"}}

   @param replacementList list of string replacement rules
   Example:
      replacementList = { "@token@" -> "1+2", "@bar@" -> "2+3" }
 *)
ReplaceInFiles[files_List, replacementList_List] :=
    Module[{cppFileName, cppTemplateFileName, cppFile, modifiedCppFile, f},
          For[f = 1, f <= Length[files], f++,
              cppFileName         = files[[f,1]];
              cppTemplateFileName = files[[f,2]];
              cppFile             = Import[cppFileName, "String"];
              modifiedCppFile     = StringReplace[cppFile, replacementList];
              Print["   writing file ", cppTemplateFileName];
              Export[cppTemplateFileName, modifiedCppFile, "String"];
             ];
          ];

GeneralReplacementRules[] :=
    { "@VectorZ@"     -> ToValidCSymbolString[SARAH`VectorZ],
      "@VectorP@"     -> ToValidCSymbolString[SARAH`VectorP],
      "@VectorW@"     -> ToValidCSymbolString[SARAH`VectorW],
      "@VectorG@"     -> ToValidCSymbolString[SARAH`VectorG],
      "@TopQuark@"    -> ToValidCSymbolString[SARAH`TopQuark],
      "@BottomQuark@" -> ToValidCSymbolString[SARAH`BottomQuark],
      "@Electron@"    -> ToValidCSymbolString[SARAH`Electron],
      "@Neutrino@"    -> ToValidCSymbolString[SARAH`Neutrino],
      "@UpYukawa@"       -> ToValidCSymbolString[SARAH`UpYukawa],
      "@DownYukawa@"     -> ToValidCSymbolString[SARAH`DownYukawa],
      "@ElectronYukawa@" -> ToValidCSymbolString[SARAH`ElectronYukawa],
      "@hyperchargeCoupling@" -> ToValidCSymbolString[SARAH`hyperchargeCoupling],
      "@leftCoupling@"        -> ToValidCSymbolString[SARAH`leftCoupling],
      "@strongCoupling@"      -> ToValidCSymbolString[SARAH`strongCoupling],
      "@ModelName@"           -> Model`Name,
      "@DateAndTime@"         -> DateString[]
    }


WriteRGEClass[betaFun_List, anomDim_List, files_List,
              additionalDecl_:"", numberOfBaseClassParameters_:0] :=
   Module[{setter, getter, parameterDef, set,
           display, parameterDefaultInit,
           cCtorParameterList, parameterCopyInit, betaParameterList,
           anomDimPrototypes, anomDimFunctions, printParameters, parameters,
           clearParameters, parameterNames},
          (* extract list of parameters from the beta functions *)
          parameters = GetName[#]& /@ betaFun;
          (* count number of parameters *)
          numberOfParameters = BetaFunction`CountNumberOfParameters[betaFun] + numberOfBaseClassParameters;
          (* create C++ functions and parameter declarations *)
          beta                 = BetaFunction`CreateBetaFunction[betaFun, additionalDecl];
          setter               = BetaFunction`CreateSetters[betaFun];
          getter               = BetaFunction`CreateGetters[betaFun];
          parameterDef         = BetaFunction`CreateParameterDefinitions[betaFun];
          set                  = BetaFunction`CreateSetFunction[betaFun, numberOfBaseClassParameters];
          display              = BetaFunction`CreateDisplayFunction[betaFun, numberOfBaseClassParameters];
          parameterDefaultInit = BetaFunction`CreateParameterDefaultInitialization[betaFun];
          cCtorParameterList   = BetaFunction`CreateCCtorParameterList[betaFun];
          parameterCopyInit    = BetaFunction`CreateCCtorInitialization[betaFun];
          betaParameterList    = BetaFunction`CreateParameterList[betaFun, "beta_"];
          clearParameters      = BetaFunction`ClearParameters[betaFun];
          parameterNames       = BetaFunction`CreateParameterNamesFunction[betaFun, numberOfBaseClassParameters];
          anomDimPrototypes    = AnomalousDimension`CreateAnomDimPrototypes[anomDim];
          anomDimFunctions     = AnomalousDimension`CreateAnomDimFunctions[anomDim];
          printParameters      = WriteOut`PrintParameters[parameters, "ostr"];
          ReplaceInFiles[files,
                 { "@beta@"                 -> IndentText[WrapLines[beta]],
                   "@clearParameters@"      -> IndentText[WrapLines[clearParameters]],
                   "@parameterDefaultInit@" -> WrapLines[parameterDefaultInit],
                   "@display@"              -> IndentText[display],
                   "@set@"                  -> IndentText[set],
                   "@cCtorParameterList@"   -> WrapLines[cCtorParameterList],
                   "@parameterCopyInit@"    -> WrapLines[parameterCopyInit],
                   "@betaParameterList@"    -> betaParameterList,
                   "@parameterDef@"         -> IndentText[parameterDef],
                   "@cCtorParameterList@"   -> WrapLines[cCtorParameterList],
                   "@setter@"               -> IndentText[setter],
                   "@getter@"               -> IndentText[getter],
                   "@anomDimPrototypes@"    -> IndentText[anomDimPrototypes],
                   "@anomDimFunctions@"     -> WrapLines[anomDimFunctions],
                   "@numberOfParameters@"   -> RValueToCFormString[numberOfParameters],
                   "@printParameters@"      -> IndentText[printParameters],
                   "@parameterNames@"       -> IndentText[parameterNames],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteInputParameterClass[inputParameters_List, freePhases_List,
                         defaultValues_List, files_List] :=
   Module[{defineInputParameters, defaultInputParametersInit},
          defineInputParameters = Constraint`DefineInputParameters[Join[inputParameters,freePhases]];
          defaultInputParametersInit = Constraint`InitializeInputParameters[Join[defaultValues, freePhases]];
          ReplaceInFiles[files,
                         { "@defineInputParameters@" -> IndentText[defineInputParameters],
                           "@defaultInputParametersInit@" -> WrapLines[defaultInputParametersInit],
                           Sequence @@ GeneralReplacementRules[]
                         } ];
          ];

WriteConstraintClass[condition_, settings_List, scaleFirstGuess_,
                     inputParameters_List, files_List] :=
   Module[{applyConstraint = "", calculateScale, scaleGuess,
           setDRbarGaugeCouplings, setDRbarYukawaCouplings},
          Constraint`SetInputParameters[inputParameters];
          Constraint`SetModelParameters[GetParameters[]];
          Constraint`SetOutputParameters[GetOutputParameters[]];
          Constraint`SetBetaFunctions[GetBetaFunctions[]];
          applyConstraint = Constraint`ApplyConstraints[settings];
          calculateScale  = Constraint`CalculateScale[condition];
          scaleGuess      = Constraint`CalculateScale[scaleFirstGuess];
          setDRbarGaugeCouplings  = ThresholdCorrections`SetDRbarGaugeCouplings[];
          setDRbarYukawaCouplings = ThresholdCorrections`SetDRbarYukawaCouplings[];
          ReplaceInFiles[files,
                 { "@applyConstraint@"      -> IndentText[WrapLines[applyConstraint]],
                   "@calculateScale@"       -> IndentText[WrapLines[calculateScale]],
                   "@scaleGuess@"           -> IndentText[WrapLines[scaleGuess]],
                   "@setDRbarGaugeCouplings@"  -> IndentText[WrapLines[setDRbarGaugeCouplings]],
                   "@setDRbarYukawaCouplings@" -> IndentText[WrapLines[setDRbarYukawaCouplings]],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteInitialGuesserClass[settings_List, files_List] :=
   Module[{applyConstraint, setDRbarYukawaCouplings},
          applyConstraint = Constraint`ApplyConstraints[settings /. { Global`MZDRbar -> Global`MZ,
                                                                      Global`MWDRbar -> Global`MW }];
          setDRbarYukawaCouplings = ThresholdCorrections`SetDRbarYukawaCouplings[];
          ReplaceInFiles[files,
                 { "@applyConstraint@"      -> IndentText[WrapLines[applyConstraint]],
                   "@setDRbarYukawaCouplings@" -> IndentText[WrapLines[setDRbarYukawaCouplings]],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteConvergenceTesterClass[particles_List, files_List] :=
   Module[{compareFunction},
          compareFunction = ConvergenceTester`CreateCompareFunction[particles];
          ReplaceInFiles[files,
                 { "@compareFunction@"      -> IndentText[WrapLines[compareFunction]],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteModelClass[massMatrices_List, tadpoleEquations_List,
                parametersFixedByEWSB_List, nPointFunctions_List, phases_List,
                files_List, diagonalizationPrecision_List] :=
    Module[{massGetters = "", k,
            mixingMatrixGetters = "",
            tadpoleEqPrototypes = "", tadpoleEqFunctions = "",
            numberOfEWSBEquations = Length[tadpoleEquations], calculateTreeLevelTadpoles = "",
            initialGuess = "", physicalMassesDef = "", mixingMatricesDef = "",
            physicalMassesInit = "", physicalMassesInitNoLeadingComma = "", mixingMatricesInit = "",
            massCalculationPrototypes = "", massCalculationFunctions = "",
            calculateAllMasses = "", calculateOneLoopTadpoles = "",
            selfEnergyPrototypes = "", selfEnergyFunctions = "",
            phasesDefinition = "", phasesGetterSetters = "",
            phasesInit = "", vevs,
            loopMassesPrototypes = "", loopMassesFunctions = "",
            runningDRbarMassesPrototypes = "", runningDRbarMassesFunctions = "",
            callAllLoopMassFunctions = "", printMasses = "", printMixingMatrices = "",
            masses, mixingMatrices, oneLoopTadpoles,
            dependenceNumPrototypes, dependenceNumFunctions,
            clearOutputParameters = ""
           },
           vevs = #[[1]]& /@ tadpoleEquations; (* list of VEVs *)
           For[k = 1, k <= Length[massMatrices], k++,
               massGetters = massGetters <> TreeMasses`CreateMassGetter[massMatrices[[k]]];
               mixingMatrixGetters = mixingMatrixGetters <> TreeMasses`CreateMixingMatrixGetter[massMatrices[[k]]];
               physicalMassesDef    = physicalMassesDef <> TreeMasses`CreatePhysicalMassDefinition[massMatrices[[k]]];
               mixingMatricesDef    = mixingMatricesDef <> TreeMasses`CreateMixingMatrixDefinition[massMatrices[[k]]];
               physicalMassesInit   = physicalMassesInit <> TreeMasses`CreatePhysicalMassInitialization[massMatrices[[k]]];
               physicalMassesInitNoLeadingComma = StringTrim[physicalMassesInit, StartOfString ~~ ","];
               mixingMatricesInit   = mixingMatricesInit <> TreeMasses`CreateMixingMatrixInitialization[massMatrices[[k]]];
               clearOutputParameters = clearOutputParameters <> TreeMasses`ClearOutputParameters[massMatrices[[k]]];
               massCalculationPrototypes = massCalculationPrototypes <> TreeMasses`CreateMassCalculationPrototype[massMatrices[[k]]];
               massCalculationFunctions  = massCalculationFunctions  <> TreeMasses`CreateMassCalculationFunction[massMatrices[[k]]];
               calculateAllMasses        = calculateAllMasses <> TreeMasses`CallMassCalculationFunction[massMatrices[[k]]];
              ];
           For[k = 1, k <= Length[tadpoleEquations], k++,
               tadpoleEqPrototypes = tadpoleEqPrototypes <> EWSB`CreateEWSBEqPrototype[tadpoleEquations[[k,1]]];
               tadpoleEqFunctions  = tadpoleEqFunctions  <> EWSB`CreateEWSBEqFunction[tadpoleEquations[[k,1]], tadpoleEquations[[k,2]]];
              ];
           If[Length[parametersFixedByEWSB] != numberOfEWSBEquations,
              Print["Error: There are ", numberOfEWSBEquations, " EWSB ",
                    "equations, but you want to fix ", Length[parametersFixedByEWSB],
                    " parameters: ", parametersFixedByEWSB];
             ];
           calculateTreeLevelTadpoles = EWSB`FillArrayWithEWSBEqs[tadpoleEquations, parametersFixedByEWSB];
           oneLoopTadpoles = Cases[nPointFunctions, SelfEnergies`Tadpole[___]];
           calculateOneLoopTadpoles   = SelfEnergies`FillArrayWithOneLoopTadpoles[oneLoopTadpoles];
           initialGuess = EWSB`FillInitialGuessArray[parametersFixedByEWSB];
           {selfEnergyPrototypes, selfEnergyFunctions} = SelfEnergies`CreateNPointFunctions[nPointFunctions];
           phasesDefinition = Phases`CreatePhasesDefinition[phases];
           phasesGetterSetters          = Phases`CreatePhasesGetterSetters[phases];
           phasesInit                   = Phases`CreatePhasesInitialization[phases];
           loopMassesPrototypes         = LoopMasses`CreateLoopMassPrototypes[];
           (* If you want to add tadpoles, call the following routine like this:
              CreateLoopMassFunctions[diagonalizationPrecision, oneLoopTadpoles, vevs];
              *)
           loopMassesFunctions          = LoopMasses`CreateLoopMassFunctions[diagonalizationPrecision, {}, {}];
           runningDRbarMassesPrototypes = LoopMasses`CreateRunningDRbarMassPrototypes[];
           runningDRbarMassesFunctions  = LoopMasses`CreateRunningDRbarMassFunctions[];
           callAllLoopMassFunctions     = LoopMasses`CallAllLoopMassFunctions[];
           masses                       = FlexibleSUSY`M[TreeMasses`GetMassEigenstate[#]]& /@ massMatrices;
           printMasses                  = WriteOut`PrintParameters[masses, "ostr"];
           mixingMatrices               = Flatten[TreeMasses`GetMixingMatrixSymbol[#]& /@ massMatrices];
           printMixingMatrices          = WriteOut`PrintParameters[mixingMatrices, "ostr"];
           dependenceNumPrototypes      = TreeMasses`CreateDependenceNumPrototypes[];
           dependenceNumFunctions       = TreeMasses`CreateDependenceNumFunctions[];
           ReplaceInFiles[files,
                          { "@massGetters@" -> IndentText[massGetters],
                            "@mixingMatrixGetters@" -> IndentText[mixingMatrixGetters],
                            "@tadpoleEqPrototypes@"  -> IndentText[tadpoleEqPrototypes],
                            "@tadpoleEqFunctions@"   -> tadpoleEqFunctions,
                            "@numberOfEWSBEquations@"-> ToString[numberOfEWSBEquations],
                            "@calculateTreeLevelTadpoles@" -> IndentText[calculateTreeLevelTadpoles],
                            "@calculateOneLoopTadpoles@"   -> IndentText[calculateOneLoopTadpoles],
                            "@clearOutputParameters@"  -> IndentText[clearOutputParameters],
                            "@initialGuess@"           -> IndentText[initialGuess],
                            "@physicalMassesDef@"      -> IndentText[physicalMassesDef],
                            "@mixingMatricesDef@"      -> IndentText[mixingMatricesDef],
                            "@physicalMassesInit@"     -> IndentText[WrapLines[physicalMassesInit]],
                            "@physicalMassesInitNoLeadingComma@" -> IndentText[WrapLines[physicalMassesInitNoLeadingComma]],
                            "@mixingMatricesInit@"     -> IndentText[WrapLines[mixingMatricesInit]],
                            "@massCalculationPrototypes@" -> IndentText[massCalculationPrototypes],
                            "@massCalculationFunctions@"  -> WrapLines[massCalculationFunctions],
                            "@calculateAllMasses@"        -> IndentText[calculateAllMasses],
                            "@selfEnergyPrototypes@"      -> IndentText[selfEnergyPrototypes],
                            "@selfEnergyFunctions@"       -> selfEnergyFunctions,
                            "@phasesDefinition@"          -> IndentText[phasesDefinition],
                            "@phasesGetterSetters@"          -> IndentText[phasesGetterSetters],
                            "@phasesInit@"                   -> IndentText[WrapLines[phasesInit]],
                            "@loopMassesPrototypes@"         -> IndentText[WrapLines[loopMassesPrototypes]],
                            "@loopMassesFunctions@"          -> WrapLines[loopMassesFunctions],
                            "@runningDRbarMassesPrototypes@" -> IndentText[runningDRbarMassesPrototypes],
                            "@runningDRbarMassesFunctions@"  -> WrapLines[runningDRbarMassesFunctions],
                            "@callAllLoopMassFunctions@"     -> IndentText[callAllLoopMassFunctions],
                            "@printMasses@"                  -> IndentText[printMasses],
                            "@printMixingMatrices@"          -> IndentText[printMixingMatrices],
                            "@dependenceNumPrototypes@"      -> IndentText[dependenceNumPrototypes],
                            "@dependenceNumFunctions@"       -> WrapLines[dependenceNumFunctions],
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteUserExample[files_List] :=
    Module[{},
           ReplaceInFiles[files,
                          { Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

FileExists[fileName_String] := FileExistsQ[fileName];

FileExists[path_String, fileName_String] :=
    Module[{fileExists},
           fileExists = FileExists[FileNameJoin[{path, fileName}]];
           If[!fileExists,
              Print["File not found: ", fileName, " in directory ", path];
             ];
           fileExists
          ];

FilesExist[path_String, fileNames_List] :=
    Module[{filesExist},
           filesExist = FileExists[path,#]& /@ fileNames;
           And @@ filesExist
          ];

RGEsHaveBeenCalculated[outputDir_String] :=
    Module[{rgeDir, fileNames},
           rgeDir = FileNameJoin[{outputDir, "RGEs"}];
           fileNames = { "BetaYijk.m", "BetaGauge.m", "BetaWijkl.m",
                         "BetaMuij.m", "BetaLi.m", "BetaQijkl.m",
                         "BetaTijk.m", "BetaBij.m", "BetaLSi.m",
                         "Betam2ij.m", "BetaMi.m", "BetaVEV.m" };
           FilesExist[rgeDir, fileNames]
          ];

SearchSelfEnergies[outputDir_String, eigenstates_] :=
    Module[{seDir, fileName = "SelfEnergy.m"},
           seDir = FileNameJoin[{outputDir, ToString[eigenstates], "One-Loop"}];
           If[FileExists[seDir, fileName],
              Return[FileNameJoin[{seDir,fileName}]];
             ];
           Return[""];
          ];

SearchUnrotatedParticles[outputDir_String, eigenstates_] :=
    Module[{dir, fileName = "UnrotatedParticles.m"},
           dir = FileNameJoin[{outputDir, ToString[eigenstates], "One-Loop"}];
           If[FileExists[dir, fileName],
              Return[FileNameJoin[{dir,fileName}]];
             ];
           Return[""];
          ];

SearchTadpoles[outputDir_String, eigenstates_] :=
    Module[{tadpoleDir, fileName = "Tadpoles1Loop.m"},
           tadpoleDir = FileNameJoin[{outputDir, ToString[eigenstates], "One-Loop"}];
           If[FileExists[tadpoleDir, fileName],
              Return[FileNameJoin[{tadpoleDir,fileName}]];
             ];
           Return[""];
          ];

PrepareRGEs[] :=
    Module[{rgesHaveBeenCalculated},
           rgesHaveBeenCalculated = RGEsHaveBeenCalculated[$sarahCurrentOutputMainDir];
           If[rgesHaveBeenCalculated,
              Print["RGEs have already been calculated, reading them from file ..."];,
              Print["RGEs have not been calculated yet, calculating them ..."];
             ];
           SARAH`CalcRGEs[ReadLists -> rgesHaveBeenCalculated,
                          TwoLoop -> True,
                          NoMatrixMultiplication -> False];
          ];

PrepareSelfEnergies[eigenstates_] :=
    Module[{selfEnergies = {}, selfEnergiesFile},
           selfEnergiesFile = SearchSelfEnergies[$sarahCurrentOutputMainDir, eigenstates];
           If[selfEnergiesFile != "",
              Print["Self-energies have already been calculated, reading them from file ", selfEnergiesFile, " ..."];
              selfEnergies = Get[selfEnergiesFile];
              ,
              Print["Self-energies have not been calculated yet, calculating them ..."];
              SARAH`CalcLoopCorrections[eigenstates];
              selfEnergiesFile = SearchSelfEnergies[$sarahCurrentOutputMainDir, eigenstates];
              selfEnergies = Get[selfEnergiesFile];
             ];
           Print["Converting self-energies ..."];
           ConvertSarahSelfEnergies[selfEnergies]
          ];

PrepareTadpoles[eigenstates_] :=
    Module[{tadpoles = {}, tadpolesFile},
           tadpolesFile = SearchTadpoles[$sarahCurrentOutputMainDir, eigenstates];
           If[tadpolesFile != "",
              Print["Tadpoles have already been calculated, reading them from file ", tadpolesFile, " ..."];
              tadpoles = Get[tadpolesFile];
              ,
              Print["Tadpoles have not been calculated yet, calculating them ..."];
              SARAH`CalcLoopCorrections[eigenstates];
              tadpolesFile = SearchTadpoles[$sarahCurrentOutputMainDir, eigenstates];
              tadpoles = Get[tadpolesFile];
             ];
           Print["Converting tadpoles ..."];
           ConvertSarahTadpoles[tadpoles]
          ];

PrepareUnrotatedParticles[eigenstates_] :=
    Module[{nonMixedParticles = {}, nonMixedParticlesFile},
           nonMixedParticlesFile = SearchUnrotatedParticles[$sarahCurrentOutputMainDir, eigenstates];
           If[nonMixedParticlesFile != "",
              Print["Unrotated particles have already been calculated, reading them from file ",
                    nonMixedParticlesFile, " ..."];
              nonMixedParticles = Get[nonMixedParticlesFile];
              ,
              Print["Unrotated particles have not been calculated yet, calculating them ..."];
              SARAH`CalcLoopCorrections[eigenstates];
              nonMixedParticlesFile = SearchUnrotatedParticles[$sarahCurrentOutputMainDir, eigenstates];
              nonMixedParticles = Get[nonMixedParticlesFile];
             ];
           TreeMasses`SetUnrotatedParticles[nonMixedParticles];
          ];

ReadDiagonalizationPrecisions[defaultPrecision_Symbol, highPrecisionList_List,
                              mediumPrecisionList_List, lowPrecisionList_List, eigenstates_] :=
    Module[{particles, particle, i, precisionList = {}},
           particles = LoopMasses`GetLoopCorrectedParticles[eigenstates];
           For[i = 1, i <= Length[particles], i++,
               particle = particles[[i]];
               Which[MemberQ[highPrecisionList  , particle], AppendTo[precisionList, {particle, HighPrecision}],
                     MemberQ[mediumPrecisionList, particle], AppendTo[precisionList, {particle, MediumPrecision}],
                     MemberQ[lowPrecisionList   , particle], AppendTo[precisionList, {particle, LowPrecision}],
                     True, AppendTo[precisionList, {particle, defaultPrecision}]
                    ];
              ];
           Return[precisionList];
          ];

LoadModelFile[file_String] :=
    Module[{},
           If[FileExists[file],
              Print["Loading model file ", file];
              Get[file];
              CheckModelFileSettings[];
              ,
              Print["Error: model file not found: ", file];
              Quit[];
             ];
          ];

Options[MakeFlexibleSUSY] :=
    {
        Eigenstates -> SARAH`EWSB,
        InputFile -> "FlexibleSUSY.m",
        softSusyCompatibleRGEs -> True,
        defaultDiagonalizationPrecision -> MediumPrecision,
        highPrecision -> {},
        mediumPrecision -> {},
        lowPrecision -> {}
    };

MakeFlexibleSUSY[OptionsPattern[]] :=
    Module[{nPointFunctions, eigenstates = OptionValue[Eigenstates],
            susyBetaFunctions, susyBreakingBetaFunctions,
            susyParameterReplacementRules, susyBreakingParameterReplacementRules,
            susyTraceDecl, susyTraceRules,
            nonSusyTraceDecl, nonSusyTraceRules,
            numberOfSusyParameters, anomDim,
            ewsbEquations, massMatrices, phases, vevs,
            diagonalizationPrecision, allParticles, freePhases},
           (* check if SARAH`Start[] was called *)
           If[!ValueQ[Model`Name],
              Print["Error: Model`Name is not defined.  Did you call SARAH`Start[\"Model\"]?"];
              Quit[];
             ];
           (* load model file *)
           LoadModelFile[OptionValue[InputFile]];
           (* get RGEs *)
           PrepareRGEs[];
           nPointFunctions = Join[PrepareSelfEnergies[eigenstates], PrepareTadpoles[eigenstates]];
           PrepareUnrotatedParticles[eigenstates];
           (* adapt SARAH`Conj to our needs *)
           (* Clear[Conj]; *)
           SARAH`Conj[(B_)[b__]] = .;
           SARAH`Conj /: SARAH`Conj[SARAH`Conj[x_]] := x;
           RXi[_] = 1;
           SARAH`Xi = 1;
           SARAH`Xip = 1;

           If[OptionValue[softSusyCompatibleRGEs] === True && Model`Name == "MSSM",
              Print["Generating SoftSusy compatible beta function ..."];
              Print["Note: SoftSusy is missing g^4 two-loop contributions to the VEVs, I'm disabling them."];
              SARAH`BetaVEV = SARAH`BetaVEV /.
              { SARAH`g1^4 -> 0,
                Susyno`LieGroups`g2^4 -> 0,
                SARAH`g3^4 -> 0,
                SARAH`g1^2 Susyno`LieGroups`g2^2 -> 0 };
             ];

           (* pick beta functions of supersymmetric parameters *)
           susyBetaFunctions = { SARAH`BetaWijkl,
                                 SARAH`BetaYijk ,
                                 SARAH`BetaMuij ,
                                 SARAH`BetaLi   ,
                                 SARAH`BetaGauge,
                                 SARAH`BetaVEV  };

           (* pick beta functions of non-supersymmetric parameters *)
           susyBreakingBetaFunctions = { SARAH`BetaQijkl,
                                         SARAH`BetaTijk ,
                                         SARAH`BetaBij  ,
                                         SARAH`BetaSLi  ,
                                         SARAH`Betam2ij ,
                                         SARAH`BetaMi   };

           susyBetaFunctions = ConvertSarahRGEs[susyBetaFunctions];
           susyBetaFunctions = Select[susyBetaFunctions, (GetAllBetaFunctions[#]!={})&];
           Parameters`AddRealParameter[(GetName /@ susyBetaFunctions) /. a_[i1,i2] :> a];
           susyParameterReplacementRules = BetaFunction`ConvertParameterNames[susyBetaFunctions];
           susyBetaFunctions = susyBetaFunctions /. susyParameterReplacementRules;

           {susyTraceDecl, susyTraceRules} = CreateDoubleTraceAbbrs[FindMultipleTraces[GetAllBetaFunctions[#]& /@ susyBetaFunctions]];
           susyBetaFunctions = susyBetaFunctions /. susyTraceRules;

           numberOfSusyParameters = BetaFunction`CountNumberOfParameters[susyBetaFunctions];
           anomDim = ConvertSarahAnomDim[SARAH`Gij];
           anomDim = anomDim /. BetaFunction`ConvertParameterNames[anomDim] /. susyParameterReplacementRules;

           Print["Creating class for susy parameters ..."];
           WriteRGEClass[susyBetaFunctions, anomDim,
                         {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "susy_parameters.hpp.in"}],
                           FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_susy_parameters.hpp"}]},
                          {FileNameJoin[{Global`$flexiblesusyTemplateDir, "susy_parameters.cpp.in"}],
                           FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_susy_parameters.cpp"}]}},
                         susyTraceDecl];

           susyBreakingBetaFunctions = ConvertSarahRGEs[susyBreakingBetaFunctions];
           susyBreakingBetaFunctions = Select[susyBreakingBetaFunctions, (GetAllBetaFunctions[#]!={})&];
           Parameters`AddRealParameter[(GetName /@ susyBreakingBetaFunctions) /. a_[i1,i2] :> a];
           susyBreakingParameterReplacementRules = Flatten[{susyParameterReplacementRules, BetaFunction`ConvertParameterNames[susyBreakingBetaFunctions]}];
           susyBreakingBetaFunctions = susyBreakingBetaFunctions /. susyBreakingParameterReplacementRules;

           {nonSusyTraceDecl, nonSusyTraceRules} = CreateDoubleTraceAbbrs[FindMultipleTraces[GetAllBetaFunctions[#]& /@ susyBreakingBetaFunctions]];

           allBetaFunctions = Join[susyBetaFunctions, susyBreakingBetaFunctions];

           {traceDecl, traceRules} = CreateTraceAbbr[SARAH`TraceAbbr /. susyBreakingParameterReplacementRules];
           traceRules = DeleteDuplicates[Flatten[{traceRules, susyTraceRules, nonSusyTraceRules}]];

           susyBreakingBetaFunctions = susyBreakingBetaFunctions /. traceRules;

           (* store all model parameters *)
           allParameters = Join[GetName /@ susyBetaFunctions,
                                GetName /@ susyBreakingBetaFunctions] /. a_[i1,i2] :> a;
           allIndexReplacementRules = Parameters`CreateIndexReplacementRules[allParameters];

           TreeMasses`SetModelParameters[allParameters];

           Print["Creating class for soft parameters ..."];
           WriteRGEClass[susyBreakingBetaFunctions, {},
                         {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "soft_parameters.hpp.in"}],
                           FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_soft_parameters.hpp"}]},
                          {FileNameJoin[{Global`$flexiblesusyTemplateDir, "soft_parameters.cpp.in"}],
                           FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_soft_parameters.cpp"}]}},
                         traceDecl <> "\n" <> nonSusyTraceDecl, numberOfSusyParameters];

           Print["Checking EWSB equations ..."];
           freePhases = EWSB`CheckEWSBEquations[SARAH`TadpoleEquations[eigenstates], ParametersToSolveTadpoles];

           Print["Creating class for input parameters ..."];
           WriteInputParameterClass[Global`InputParameters, freePhases,
                                    Global`DefaultParameterPoint,
                                    {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "input_parameters.hpp.in"}],
                                      FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_input_parameters.hpp"}]}}
                                   ];

           massMatrices = ConvertSarahMassMatrices[] /.
                          susyBreakingParameterReplacementRules /.
                          Parameters`ApplyGUTNormalization[] /.
                          { SARAH`sum[j_, start_, end_, expr_] :> (Sum[expr, {j,start,end}]) } /.
                          allIndexReplacementRules;

           allParticles = FlexibleSUSY`M[GetMassEigenstate[#]]& /@ massMatrices;
           allOutputParameters = DeleteCases[DeleteDuplicates[
               Join[allParticles,
                    Flatten[GetMixingMatrixSymbol[#]& /@ massMatrices]]], Null];

           Print["Creating class for convergence tester ..."];
           WriteConvergenceTesterClass[allParticles,
               {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "convergence_tester.hpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_convergence_tester.hpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "convergence_tester.cpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_convergence_tester.cpp"}]}}
                                      ];

           Print["Creating class for high-scale constraint ..."];
           WriteConstraintClass[SARAH`ConditionGUTscale /. susyBreakingParameterReplacementRules,
                                SARAH`BoundaryHighScale /. susyBreakingParameterReplacementRules,
                                Global`BoundaryHighScaleFirstGuess /. susyBreakingParameterReplacementRules,
                                Global`InputParameters,
                                {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "high_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_high_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "high_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_high_scale_constraint.cpp"}]}}
                               ];

           Print["Creating class for susy-scale constraint ..."];
           WriteConstraintClass[SARAH`RenormalizationScale /. susyBreakingParameterReplacementRules,
                                SARAH`BoundarySUSYScale /. susyBreakingParameterReplacementRules,
                                SARAH`RenormalizationScaleFirstGuess /. susyBreakingParameterReplacementRules,
                                Global`InputParameters,
                                {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "susy_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_susy_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "susy_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_susy_scale_constraint.cpp"}]}}
                               ];

           Print["Creating class for low-scale constraint ..."];
           WriteConstraintClass[Global`BoundaryLowScale /. susyBreakingParameterReplacementRules,
                                SARAH`BoundaryLowScaleInput /. susyBreakingParameterReplacementRules,
                                Global`MZ,
                                Global`InputParameters,
                                {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "low_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_low_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "low_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_low_scale_constraint.cpp"}]}}
                               ];

           Print["Creating class for initial guesser ..."];
           WriteInitialGuesserClass[Global`InitialGuess /. susyBreakingParameterReplacementRules,
                                    {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "initial_guesser.hpp.in"}],
                                      FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_initial_guesser.hpp"}]},
                                     {FileNameJoin[{Global`$flexiblesusyTemplateDir, "initial_guesser.cpp.in"}],
                                      FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_initial_guesser.cpp"}]}}
                                   ];

           vevs = #[[1]]& /@ SARAH`BetaVEV;
           ewsbEquations = SARAH`TadpoleEquations[eigenstates] /.
                           susyBreakingParameterReplacementRules /.
                           Parameters`ApplyGUTNormalization[];
           ewsbEquations = MapThread[List, {vevs, ewsbEquations}];

           SelfEnergies`SetParameterReplacementRules[susyBreakingParameterReplacementRules];
           SelfEnergies`SetIndexReplacementRules[allIndexReplacementRules];

           phases = ConvertSarahPhases[SARAH`ParticlePhases];

           (* determin diagonalization precision for each particle *)
           diagonalizationPrecision = ReadDiagonalizationPrecisions[
               OptionValue[defaultDiagonalizationPrecision],
               Flatten[{OptionValue[highPrecision]}],
               Flatten[{OptionValue[mediumPrecision]}],
               Flatten[{OptionValue[lowPrecision]}],
               eigenstates];

           Print["Creating class for model ..."];
           WriteModelClass[massMatrices, ewsbEquations,
                           ParametersToSolveTadpoles,
                           nPointFunctions, phases,
                           {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "model.hpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_model.hpp"}]},
                            {FileNameJoin[{Global`$flexiblesusyTemplateDir, "model.cpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_model.cpp"}]},
                            {FileNameJoin[{Global`$flexiblesusyTemplateDir, "physical.hpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_physical.hpp"}]},
                            {FileNameJoin[{Global`$flexiblesusyTemplateDir, "physical.cpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, Model`Name <> "_physical.cpp"}]}},
                           diagonalizationPrecision];

           Print["Creating user example spectrum generator program ..."];
           WriteUserExample[{{FileNameJoin[{Global`$flexiblesusyTemplateDir, "run.cpp.in"}],
                              FileNameJoin[{Global`$flexiblesusyOutputDir, "run_" <> Model`Name <> ".cpp"}]}}];
          ];

End[];

EndPackage[];

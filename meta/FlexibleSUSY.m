
BeginPackage["FlexibleSUSY`", {"SARAH`", "AnomalousDimension`", "BetaFunction`", "TextFormatting`", "CConversion`", "TreeMasses`", "EWSB`", "Traces`", "SelfEnergies`", "Phases`", "LoopMasses`", "WriteOut`", "Constraint`", "ThresholdCorrections`", "ConvergenceTester`"}];

MakeFlexibleSUSY::usage="";

LowPrecision::usage="";
MediumPrecision::usage="";
HighPrecision::usage="";
softSusyCompatibleRGEs::usage="";
GUTNormalization::usage="Returns GUT normalization of a given coupling";

FSModelName;
InputParameters;
DefaultParameterPoint;
ParametersToSolveTadpoles;
SUSYScale;
SUSYScaleFirstGuess;
SUSYScaleInput;
HighScale;
HighScaleFirstGuess;
HighScaleInput;
LowScale;
LowScaleFirstGuess;
LowScaleInput;
InitialGuessAtLowScale;
InitialGuessAtHighScale;
MZ;

Begin["Private`"];

allParameters = {};

allIndexReplacementRules = {};

GetIndexReplacementRules[] := allIndexReplacementRules;

allBetaFunctions = {};

GetBetaFunctions[] := allBetaFunctions;

allOutputParameters = {};

numberOfModelParameters = 0;

CheckModelFileSettings[] :=
    Module[{},
           (* FlexibleSUSY model name *)
           If[!NameQ["FlexibleSUSY`FSModelName"] || Head[FlexibleSUSY`FSModelName] =!= String,
              Print["Warning: FlexibleSUSY`FSModelName not defined!",
                    " I'm using Model`Name from SARAH: ", Model`Name];
              FlexibleSUSY`FSModelName = Model`Name;
             ];
           If[Head[FlexibleSUSY`InitialGuessAtLowScale] =!= List,
              FlexibleSUSY`InitialGuessAtLowScale = {};
             ];
           If[Head[FlexibleSUSY`InitialGuessAtHighScale] =!= List,
              FlexibleSUSY`InitialGuessAtHighScale = {};
             ];
           (* HighScale *)
           If[!NameQ["FlexibleSUSY`HighScale"],
              Print["Warning: FlexibleSUSY`HighScale should be",
                    " set in the model file!"];
              FlexibleSUSY`HighScale := SARAH`hyperchargeCoupling == SARAH`leftCoupling;
             ];
           If[!NameQ["FlexibleSUSY`HighScaleFirstGuess"],
              Print["Warning: FlexibleSUSY`HighScaleFirstGuess should be",
                    " set in the model file!"];
              FlexibleSUSY`HighScaleFirstGuess = 1.0 10^14;
             ];
           If[Head[FlexibleSUSY`HighScaleInput] =!= List,
              FlexibleSUSY`HighScaleInput = {};
             ];
           (* LowScale *)
           If[!NameQ["FlexibleSUSY`LowScale"],
              Print["Warning: FlexibleSUSY`LowScale should be",
                    " set in the model file!"];
              FlexibleSUSY`LowScale := FlexibleSUSY`MZ;
             ];
           If[!NameQ["FlexibleSUSY`LowScaleFirstGuess"],
              Print["Warning: FlexibleSUSY`LowScaleFirstGuess should be",
                    " set in the model file!"];
              FlexibleSUSY`LowScaleFirstGuess = FlexibleSUSY`MZ;
             ];
           If[Head[FlexibleSUSY`LowScaleInput] =!= List,
              FlexibleSUSY`LowScaleInput = {};
             ];
           (* SusyScale *)
           If[!NameQ["FlexibleSUSY`SusyScale"],
              Print["Warning: FlexibleSUSY`SusyScale should be",
                    " set in the model file!"];
              FlexibleSUSY`SusyScale := 1000;
             ];
           If[!NameQ["FlexibleSUSY`SusyScaleFirstGuess"],
              Print["Warning: FlexibleSUSY`SusyScaleFirstGuess should be",
                    " set in the model file!"];
              FlexibleSUSY`SusyScaleFirstGuess = 1000;
             ];
           If[Head[FlexibleSUSY`SusyScaleInput] =!= List,
              FlexibleSUSY`SusyScaleInput = {};
             ];

           If[Head[FlexibleSUSY`DefaultParameterPoint] =!= List,
              FlexibleSUSY`DefaultParameterPoint = {};
             ];
           If[Head[FlexibleSUSY`InputParameters] =!= List,
              FlexibleSUSY`InputParameters = {};
             ];
          ];

GUTNormalization[coupling_] :=
    Parameters`GetGUTNormalization[coupling];

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
              Print["   Writing file ", cppTemplateFileName];
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
      "@hyperchargeCouplingGutNormalization@"  -> RValueToCFormString[Parameters`GetGUTNormalization[SARAH`hyperchargeCoupling]],
      "@leftCouplingGutNormalization@"  -> RValueToCFormString[Parameters`GetGUTNormalization[SARAH`leftCoupling]],
      "@hyperchargeCouplingInverseGutNormalization@" -> RValueToCFormString[1/Parameters`GetGUTNormalization[SARAH`hyperchargeCoupling]],
      "@leftCouplingInverseGutNormalization@" -> RValueToCFormString[1/Parameters`GetGUTNormalization[SARAH`leftCoupling]],
      "@ModelName@"           -> FlexibleSUSY`FSModelName,
      "@numberOfModelParameters@" -> ToString[numberOfModelParameters],
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

WriteConstraintClass[condition_, settings_List, scaleFirstGuess_, files_List] :=
   Module[{applyConstraint = "", calculateScale, scaleGuess,
           setDRbarYukawaCouplings,
           calculateDeltaAlphaEm, calculateDeltaAlphaS,
           saveEwsbOutputParameters, restoreEwsbOutputParameters},
          Constraint`SetBetaFunctions[GetBetaFunctions[]];
          applyConstraint = Constraint`ApplyConstraints[settings];
          calculateScale  = Constraint`CalculateScale[condition, "scale"];
          scaleGuess      = Constraint`CalculateScale[scaleFirstGuess, "initial_scale_guess"];
          calculateDeltaAlphaEm   = ThresholdCorrections`CalculateDeltaAlphaEm[];
          calculateDeltaAlphaS    = ThresholdCorrections`CalculateDeltaAlphaS[];
          setDRbarYukawaCouplings = ThresholdCorrections`SetDRbarYukawaCouplings[];
          saveEwsbOutputParameters    = Parameters`SaveParameterLocally[ParametersToSolveTadpoles, "old_", "MODELPARAMETER"];
          restoreEwsbOutputParameters = Parameters`RestoreParameter[ParametersToSolveTadpoles, "old_", "model"];
          ReplaceInFiles[files,
                 { "@applyConstraint@"      -> IndentText[WrapLines[applyConstraint]],
                   "@calculateScale@"       -> IndentText[WrapLines[calculateScale]],
                   "@scaleGuess@"           -> IndentText[WrapLines[scaleGuess]],
                   "@calculateDeltaAlphaEm@" -> IndentText[WrapLines[calculateDeltaAlphaEm]],
                   "@calculateDeltaAlphaS@"  -> IndentText[WrapLines[calculateDeltaAlphaS]],
                   "@setDRbarYukawaCouplings@" -> IndentText[WrapLines[setDRbarYukawaCouplings]],
                   "@saveEwsbOutputParameters@"    -> IndentText[saveEwsbOutputParameters],
                   "@restoreEwsbOutputParameters@" -> IndentText[restoreEwsbOutputParameters],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteInitialGuesserClass[lowScaleGuess_List, highScaleGuess_List, files_List] :=
   Module[{initialGuessAtLowScale, initialGuessAtHighScale, setDRbarYukawaCouplings},
          initialGuessAtLowScale  = Constraint`ApplyConstraints[lowScaleGuess];
          initialGuessAtHighScale = Constraint`ApplyConstraints[highScaleGuess];
          setDRbarYukawaCouplings = ThresholdCorrections`SetDRbarYukawaCouplings[];
          ReplaceInFiles[files,
                 { "@initialGuessAtLowScale@"  -> IndentText[WrapLines[initialGuessAtLowScale]],
                   "@initialGuessAtHighScale@" -> IndentText[WrapLines[initialGuessAtHighScale]],
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

WriteModelClass[massMatrices_List, ewsbEquations_List,
                parametersFixedByEWSB_List, nPointFunctions_List, phases_List,
                files_List, diagonalizationPrecision_List] :=
    Module[{massGetters = "", k,
            mixingMatrixGetters = "",
            tadpoleEqPrototypes = "", tadpoleEqFunctions = "",
            numberOfEWSBEquations = Length[ewsbEquations], calculateTreeLevelTadpoles = "",
            ewsbInitialGuess = "", physicalMassesDef = "", mixingMatricesDef = "",
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
            clearOutputParameters = "", solveEwsbTreeLevel = "",
            saveEwsbOutputParameters, restoreEwsbOutputParameters,
            softScalarMasses, softHiggsMasses,
            saveSoftHiggsMasses, restoreSoftHiggsMasses,
            solveTreeLevelEWSBviaSoftHiggsMasses
           },
           vevs = #[[1]]& /@ ewsbEquations; (* list of VEVs *)
           For[k = 1, k <= Length[massMatrices], k++,
               massGetters          = massGetters <> TreeMasses`CreateMassGetter[massMatrices[[k]]];
               mixingMatrixGetters  = mixingMatrixGetters <> TreeMasses`CreateMixingMatrixGetter[massMatrices[[k]]];
               physicalMassesDef    = physicalMassesDef <> TreeMasses`CreatePhysicalMassDefinition[massMatrices[[k]]];
               mixingMatricesDef    = mixingMatricesDef <> TreeMasses`CreateMixingMatrixDefinition[massMatrices[[k]]];
               physicalMassesInit   = physicalMassesInit <> TreeMasses`CreatePhysicalMassInitialization[massMatrices[[k]]];
               physicalMassesInitNoLeadingComma = StringTrim[physicalMassesInit, StartOfString ~~ ","];
               mixingMatricesInit   = mixingMatricesInit <> TreeMasses`CreateMixingMatrixInitialization[massMatrices[[k]]];
               clearOutputParameters = clearOutputParameters <> TreeMasses`ClearOutputParameters[massMatrices[[k]]];
               massCalculationPrototypes = massCalculationPrototypes <> TreeMasses`CreateMassCalculationPrototype[massMatrices[[k]]];
               massCalculationFunctions  = massCalculationFunctions  <> TreeMasses`CreateMassCalculationFunction[massMatrices[[k]]];
              ];
           calculateAllMasses = TreeMasses`CallMassCalculationFunctions[massMatrices];
           For[k = 1, k <= Length[ewsbEquations], k++,
               tadpoleEqPrototypes = tadpoleEqPrototypes <> EWSB`CreateEWSBEqPrototype[ewsbEquations[[k,1]]];
               tadpoleEqFunctions  = tadpoleEqFunctions  <> EWSB`CreateEWSBEqFunction[ewsbEquations[[k,1]], ewsbEquations[[k,2]]];
              ];
           If[Length[parametersFixedByEWSB] != numberOfEWSBEquations,
              Print["Error: There are ", numberOfEWSBEquations, " EWSB ",
                    "equations, but you want to fix ", Length[parametersFixedByEWSB],
                    " parameters: ", parametersFixedByEWSB];
             ];
           oneLoopTadpoles              = Cases[nPointFunctions, SelfEnergies`Tadpole[___]];
           calculateOneLoopTadpoles     = SelfEnergies`FillArrayWithOneLoopTadpoles[oneLoopTadpoles];
           calculateTreeLevelTadpoles   = EWSB`FillArrayWithEWSBEqs[ewsbEquations, parametersFixedByEWSB];
           ewsbInitialGuess             = EWSB`FillInitialGuessArray[parametersFixedByEWSB];
           solveEwsbTreeLevel           = EWSB`SolveTreeLevelEwsb[ewsbEquations, parametersFixedByEWSB];
           {selfEnergyPrototypes, selfEnergyFunctions} = SelfEnergies`CreateNPointFunctions[nPointFunctions];
           phasesDefinition             = Phases`CreatePhasesDefinition[phases];
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
           saveEwsbOutputParameters     = Parameters`SaveParameterLocally[ParametersToSolveTadpoles, "one_loop_", ""];
           restoreEwsbOutputParameters  = Parameters`RestoreParameter[ParametersToSolveTadpoles, "one_loop_", ""];
           softScalarMasses             = DeleteDuplicates[SARAH`ListSoftBreakingScalarMasses];
           softHiggsMasses              = Select[softScalarMasses, (!FreeQ[ewsbEquations, #])&];
           saveSoftHiggsMasses          = Parameters`SaveParameterLocally[softHiggsMasses, "old_", ""];
           restoreSoftHiggsMasses       = Parameters`RestoreParameter[softHiggsMasses, "old_", ""];
           solveTreeLevelEWSBviaSoftHiggsMasses = EWSB`SolveTreeLevelEwsbVia[ewsbEquations, softHiggsMasses];
           ReplaceInFiles[files,
                          { "@massGetters@"          -> IndentText[massGetters],
                            "@mixingMatrixGetters@"  -> IndentText[mixingMatrixGetters],
                            "@tadpoleEqPrototypes@"  -> IndentText[tadpoleEqPrototypes],
                            "@tadpoleEqFunctions@"   -> tadpoleEqFunctions,
                            "@numberOfEWSBEquations@"-> ToString[numberOfEWSBEquations],
                            "@calculateTreeLevelTadpoles@" -> IndentText[calculateTreeLevelTadpoles],
                            "@calculateOneLoopTadpoles@"   -> IndentText[calculateOneLoopTadpoles],
                            "@clearOutputParameters@"  -> IndentText[clearOutputParameters],
                            "@ewsbInitialGuess@"       -> IndentText[ewsbInitialGuess],
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
                            "@solveEwsbTreeLevel@"           -> IndentText[WrapLines[solveEwsbTreeLevel]],
                            "@saveEwsbOutputParameters@"     -> IndentText[saveEwsbOutputParameters],
                            "@restoreEwsbOutputParameters@"  -> IndentText[restoreEwsbOutputParameters],
                            "@saveSoftHiggsMasses@"          -> IndentText[saveSoftHiggsMasses],
                            "@restoreSoftHiggsMasses@"       -> IndentText[restoreSoftHiggsMasses],
                            "@solveTreeLevelEWSBviaSoftHiggsMasses@" -> IndentText[WrapLines[solveTreeLevelEWSBviaSoftHiggsMasses]],
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteUserExample[files_List] :=
    Module[{},
           ReplaceInFiles[files,
                          { Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WritePlotScripts[files_List] :=
    Module[{},
           ReplaceInFiles[files,
                          { Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteUtilitiesClass[massMatrices_List, files_List] :=
    Module[{k, particles, susyParticles, smParticles,
            fillSpectrumVectorWithSusyParticles = "",
            fillSpectrumVectorWithSMParticles = "",
            particleLaTeXNames = "",
            particleNames = "", particleEnum = "", particleMultiplicity = ""},
           particles = GetMassEigenstate /@ massMatrices;
           susyParticles = Select[particles, (!SARAH`SMQ[#])&];
           smParticles   = Complement[particles, susyParticles];
           particleEnum       = TreeMasses`CreateParticleEnum[particles];
           particleMultiplicity = TreeMasses`CreateParticleMultiplicity[particles];
           particleNames      = TreeMasses`CreateParticleNames[particles];
           particleLaTeXNames = TreeMasses`CreateParticleLaTeXNames[particles];
           fillSpectrumVectorWithSusyParticles = TreeMasses`FillSpectrumVector[susyParticles];
           fillSpectrumVectorWithSMParticles   = TreeMasses`FillSpectrumVector[smParticles];
           ReplaceInFiles[files,
                          { "@fillSpectrumVectorWithSusyParticles@" -> IndentText[fillSpectrumVectorWithSusyParticles],
                            "@fillSpectrumVectorWithSMParticles@"   -> IndentText[IndentText[fillSpectrumVectorWithSMParticles]],
                            "@particleEnum@"       -> IndentText[WrapLines[particleEnum]],
                            "@particleMultiplicity@" -> IndentText[WrapLines[particleMultiplicity]],
                            "@particleNames@"      -> IndentText[WrapLines[particleNames]],
                            "@particleLaTeXNames@" -> IndentText[particleLaTeXNames],
                            Sequence @@ GeneralReplacementRules[]
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

           Parameters`SetInputParameters[FlexibleSUSY`InputParameters];

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
           susyBetaFunctions = Select[susyBetaFunctions, (BetaFunction`GetAllBetaFunctions[#]!={})&];
           Parameters`AddRealParameter[(GetName /@ susyBetaFunctions) /. a_[i1,i2] :> a];
           susyParameterReplacementRules = BetaFunction`ConvertParameterNames[susyBetaFunctions];
           susyBetaFunctions = susyBetaFunctions /. susyParameterReplacementRules;

           numberOfSusyParameters = BetaFunction`CountNumberOfParameters[susyBetaFunctions];
           anomDim = ConvertSarahAnomDim[SARAH`Gij];
           anomDim = anomDim /. BetaFunction`ConvertParameterNames[anomDim] /. susyParameterReplacementRules;

           Print["Creating class for susy parameters ..."];
           WriteRGEClass[susyBetaFunctions, anomDim,
                         {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "susy_parameters.hpp.in"}],
                           FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_susy_parameters.hpp"}]},
                          {FileNameJoin[{Global`$flexiblesusyTemplateDir, "susy_parameters.cpp.in"}],
                           FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_susy_parameters.cpp"}]}}];

           susyBreakingBetaFunctions = ConvertSarahRGEs[susyBreakingBetaFunctions];
           susyBreakingBetaFunctions = Select[susyBreakingBetaFunctions, (BetaFunction`GetAllBetaFunctions[#]!={})&];
           Parameters`AddRealParameter[(GetName /@ susyBreakingBetaFunctions) /. a_[i1,i2] :> a];
           susyBreakingParameterReplacementRules = Flatten[{susyParameterReplacementRules, BetaFunction`ConvertParameterNames[susyBreakingBetaFunctions]}];
           susyBreakingBetaFunctions = susyBreakingBetaFunctions /. susyBreakingParameterReplacementRules;

           allBetaFunctions = Join[susyBetaFunctions, susyBreakingBetaFunctions];

           {traceDecl, traceRules} = CreateTraceAbbr[SARAH`TraceAbbr /. susyBreakingParameterReplacementRules];
           susyBreakingBetaFunctions = susyBreakingBetaFunctions /. traceRules;

           (* store all model parameters *)
           allParameters = Join[GetName /@ susyBetaFunctions,
                                GetName /@ susyBreakingBetaFunctions] /. a_[i1,i2] :> a;
           allIndexReplacementRules = Parameters`CreateIndexReplacementRules[allParameters];
           Parameters`SetModelParameters[allParameters];

           numberOfSusyBreakingParameters = BetaFunction`CountNumberOfParameters[susyBreakingBetaFunctions];
           numberOfModelParameters = numberOfSusyParameters + numberOfSusyBreakingParameters;

           Print["Creating class for soft parameters ..."];
           WriteRGEClass[susyBreakingBetaFunctions, {},
                         {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "soft_parameters.hpp.in"}],
                           FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_soft_parameters.hpp"}]},
                          {FileNameJoin[{Global`$flexiblesusyTemplateDir, "soft_parameters.cpp.in"}],
                           FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_soft_parameters.cpp"}]}},
                         traceDecl, numberOfSusyParameters];

           Print["Checking EWSB equations ..."];
           freePhases = EWSB`CheckEWSBEquations[SARAH`TadpoleEquations[eigenstates], ParametersToSolveTadpoles];

           Print["Creating plot scripts ..."];
           WritePlotScripts[{{FileNameJoin[{Global`$flexiblesusyTemplateDir, "plot_spectrum.gnuplot.in"}],
                              FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_plot_spectrum.gnuplot"}]},
                             {FileNameJoin[{Global`$flexiblesusyTemplateDir, "plot_rge_running.gnuplot.in"}],
                              FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_plot_rge_running.gnuplot"}]}}
                           ];

           Print["Creating class for input parameters ..."];
           WriteInputParameterClass[FlexibleSUSY`InputParameters, freePhases,
                                    FlexibleSUSY`DefaultParameterPoint,
                                    {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "input_parameters.hpp.in"}],
                                      FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_input_parameters.hpp"}]}}
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

          Parameters`SetOutputParameters[allOutputParameters];

           Print["Creating class for convergence tester ..."];
           WriteConvergenceTesterClass[allParticles,
               {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "convergence_tester.hpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_convergence_tester.hpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "convergence_tester.cpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_convergence_tester.cpp"}]}}
                                      ];

           Print["Creating utilities class ..."];
           WriteUtilitiesClass[massMatrices,
               {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "utilities.hpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_utilities.hpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "utilities.cpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_utilities.cpp"}]}}
                              ];

           Print["Creating class for high-scale constraint ..."];
           WriteConstraintClass[FlexibleSUSY`HighScale /. susyBreakingParameterReplacementRules,
                                FlexibleSUSY`HighScaleInput /. susyBreakingParameterReplacementRules,
                                FlexibleSUSY`HighScaleFirstGuess /. susyBreakingParameterReplacementRules,
                                {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "high_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_high_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "high_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_high_scale_constraint.cpp"}]}}
                               ];

           Print["Creating class for susy-scale constraint ..."];
           WriteConstraintClass[FlexibleSUSY`SUSYScale /. susyBreakingParameterReplacementRules,
                                FlexibleSUSY`SUSYScaleInput /. susyBreakingParameterReplacementRules,
                                FlexibleSUSY`SUSYScaleFirstGuess /. susyBreakingParameterReplacementRules,
                                {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "susy_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_susy_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "susy_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_susy_scale_constraint.cpp"}]}}
                               ];

           Print["Creating class for low-scale constraint ..."];
           WriteConstraintClass[FlexibleSUSY`LowScale /. susyBreakingParameterReplacementRules,
                                FlexibleSUSY`LowScaleInput /. susyBreakingParameterReplacementRules,
                                FlexibleSUSY`LowScaleFirstGuess,
                                {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "low_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_low_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "low_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_low_scale_constraint.cpp"}]}}
                               ];

           Print["Creating class for initial guesser ..."];
           WriteInitialGuesserClass[FlexibleSUSY`InitialGuessAtLowScale /. susyBreakingParameterReplacementRules,
                                    FlexibleSUSY`InitialGuessAtHighScale /. susyBreakingParameterReplacementRules,
                                    {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "initial_guesser.hpp.in"}],
                                      FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_initial_guesser.hpp"}]},
                                     {FileNameJoin[{Global`$flexiblesusyTemplateDir, "initial_guesser.cpp.in"}],
                                      FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_initial_guesser.cpp"}]}}
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
                             FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_model.hpp"}]},
                            {FileNameJoin[{Global`$flexiblesusyTemplateDir, "model.cpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_model.cpp"}]},
                            {FileNameJoin[{Global`$flexiblesusyTemplateDir, "physical.hpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_physical.hpp"}]},
                            {FileNameJoin[{Global`$flexiblesusyTemplateDir, "physical.cpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_physical.cpp"}]}},
                           diagonalizationPrecision];

           Print["Creating user example spectrum generator program ..."];
           WriteUserExample[{{FileNameJoin[{Global`$flexiblesusyTemplateDir, "run.cpp.in"}],
                              FileNameJoin[{Global`$flexiblesusyOutputDir, "run_" <> FlexibleSUSY`FSModelName <> ".cpp"}]}}];
          ];

End[];

EndPackage[];

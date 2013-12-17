
BeginPackage["FlexibleSUSY`", {"SARAH`", "AnomalousDimension`", "BetaFunction`", "TextFormatting`", "CConversion`", "TreeMasses`", "EWSB`", "Traces`", "SelfEnergies`", "Vertices`", "Phases`", "LoopMasses`", "WriteOut`", "Constraint`", "ThresholdCorrections`", "ConvergenceTester`"}];

FS`Version = StringTrim[Import[FileNameJoin[{Global`$flexiblesusyConfigDir,"version"}], "String"]];
FS`Authors = {"P. Athron", "J. Park", "D. Stöckinger", "A. Voigt"};
FS`Years   = {2013};

Print["*****************************************************"];
Print["FlexibleSUSY ", FS`Version];
Print["by " <> WriteOut`StringJoinWithSeparator[FS`Authors, ", "] <> ", " <>
      WriteOut`StringJoinWithSeparator[FS`Years, ", "]];
Print["*****************************************************"];
Print[""];

MakeFlexibleSUSY::usage="";

LowPrecision::usage="";
MediumPrecision::usage="";
HighPrecision::usage="";
GUTNormalization::usage="Returns GUT normalization of a given coupling";

FSModelName;
FSLesHouchesList;
FSUnfixedParameters;
InputParameters;
DefaultParameterPoint;
EWSBOutputParameters;
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
OnlyLowEnergyFlexibleSUSY;
TreeLevelEWSBSolution;
Pole;
FSMinimize;
FSFindRoot;
MZ;

FSEigenstates;
FSSolveEWSBTimeConstraint = 120;

Begin["`Private`"];

allParameters = {};

allIndexReplacementRules = {};

GetIndexReplacementRules[] := allIndexReplacementRules;

allBetaFunctions = {};

GetBetaFunctions[] := allBetaFunctions;

allOutputParameters = {};

numberOfModelParameters = 0;

PrintHeadline[text_] :=
    Block[{},
          Print[""];
          Print["---------------------------------"];
          Print[text];
          Print["---------------------------------"];
         ];

DecomposeVersionString[version_String] :=
    ToExpression /@ StringSplit[version, "."];

ToVersionString[{major_Integer, minor_Integer, patch_Integer}] :=
    ToString[major] <> "." <> ToString[minor] <> "." <> ToString[patch];

CheckSARAHVersion[] :=
    Module[{minimRequired, minimRequiredVersionFile, sarahVersion},
           Print["Checking SARAH version ..."];
           minimRequiredVersionFile = FileNameJoin[{Global`$flexiblesusyConfigDir,
                                                    "required_sarah_version.m"}];
           (* reading minimum required SARAH version from config file *)
           minimRequired = Get[minimRequiredVersionFile];
           If[minimRequired === $Failed,
              Print["Error: Cannot read required SARAH version from file ",
                    minimRequiredVersionFile];
              Print["   Did you run configure?"];
              Quit[1];
             ];
           sarahVersion = DecomposeVersionString[SA`Version];
           If[sarahVersion[[1]] < minimRequired[[1]] ||
              (sarahVersion[[1]] == minimRequired[[1]] &&
               sarahVersion[[2]] < minimRequired[[2]]) ||
              (sarahVersion[[1]] == minimRequired[[1]] &&
               sarahVersion[[2]] == minimRequired[[2]] &&
               sarahVersion[[3]] < minimRequired[[3]]),
              Print["Error: SARAH version ", SA`Version, " no longer supported!"];
              Print["Please use version ", ToVersionString[minimRequired],
                    " or higher"];
              Quit[1];
             ];
          ];

CheckModelFileSettings[] :=
    Module[{},
           (* FlexibleSUSY model name *)
           If[!ValueQ[FlexibleSUSY`FSModelName] || Head[FlexibleSUSY`FSModelName] =!= String,
              Print["Warning: FlexibleSUSY`FSModelName not defined!",
                    " I'm using Model`Name from SARAH: ", Model`Name];
              FlexibleSUSY`FSModelName = Model`Name;
             ];
           (* Set OnlyLowEnergyFlexibleSUSY to False by default *)
           If[!ValueQ[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY] ||
              (FlexibleSUSY`OnlyLowEnergyFlexibleSUSY =!= True &&
               FlexibleSUSY`OnlyLowEnergyFlexibleSUSY =!= False),
              FlexibleSUSY`OnlyLowEnergyFlexibleSUSY = False;
             ];
           If[Head[FlexibleSUSY`InitialGuessAtLowScale] =!= List,
              FlexibleSUSY`InitialGuessAtLowScale = {};
             ];
           If[Head[FlexibleSUSY`InitialGuessAtHighScale] =!= List,
              FlexibleSUSY`InitialGuessAtHighScale = {};
             ];
           (* HighScale *)
           If[!ValueQ[FlexibleSUSY`HighScale],
              Print["Warning: FlexibleSUSY`HighScale should be",
                    " set in the model file!"];
              FlexibleSUSY`HighScale := SARAH`hyperchargeCoupling == SARAH`leftCoupling;
             ];
           If[!ValueQ[FlexibleSUSY`HighScaleFirstGuess],
              Print["Warning: FlexibleSUSY`HighScaleFirstGuess should be",
                    " set in the model file!"];
              FlexibleSUSY`HighScaleFirstGuess = 1.0 10^14;
             ];
           If[Head[FlexibleSUSY`HighScaleInput] =!= List,
              FlexibleSUSY`HighScaleInput = {};
             ];
           (* LowScale *)
           If[!ValueQ[FlexibleSUSY`LowScale],
              Print["Warning: FlexibleSUSY`LowScale should be",
                    " set in the model file!"];
              FlexibleSUSY`LowScale := SM[MZ];
             ];
           If[!ValueQ[FlexibleSUSY`LowScaleFirstGuess],
              Print["Warning: FlexibleSUSY`LowScaleFirstGuess should be",
                    " set in the model file!"];
              FlexibleSUSY`LowScaleFirstGuess = SM[MZ];
             ];
           If[Head[FlexibleSUSY`LowScaleInput] =!= List,
              FlexibleSUSY`LowScaleInput = {};
             ];
           (* SUSYScale *)
           If[!ValueQ[FlexibleSUSY`SUSYScale],
              Print["Warning: FlexibleSUSY`SUSYScale should be",
                    " set in the model file!"];
              FlexibleSUSY`SUSYScale := 1000;
             ];
           If[!ValueQ[FlexibleSUSY`SUSYScaleFirstGuess],
              Print["Warning: FlexibleSUSY`SUSYScaleFirstGuess should be",
                    " set in the model file!"];
              FlexibleSUSY`SUSYScaleFirstGuess = 1000;
             ];
           If[Head[FlexibleSUSY`SUSYScaleInput] =!= List,
              FlexibleSUSY`SUSYScaleInput = {};
             ];

           If[Head[FlexibleSUSY`DefaultParameterPoint] =!= List,
              FlexibleSUSY`DefaultParameterPoint = {};
             ];
           If[Head[SARAH`MINPAR] =!= List,
              SARAH`MINPAR = {};
             ];
           If[Head[SARAH`EXTPAR] =!= List,
              SARAH`EXTPAR = {};
             ];
           If[Head[FlexibleSUSY`TreeLevelEWSBSolution] =!= List,
              FlexibleSUSY`TreeLevelEWSBSolution = {};
             ];
          ];

ReplaceIndicesInUserInput[] :=
    Block[{},
          FlexibleSUSY`InitialGuessAtLowScale  = FlexibleSUSY`InitialGuessAtLowScale  /. allIndexReplacementRules;
          FlexibleSUSY`InitialGuessAtHighScale = FlexibleSUSY`InitialGuessAtHighScale /. allIndexReplacementRules;
          FlexibleSUSY`HighScale               = FlexibleSUSY`HighScale               /. allIndexReplacementRules;
          FlexibleSUSY`HighScaleFirstGuess     = FlexibleSUSY`HighScaleFirstGuess     /. allIndexReplacementRules;
          FlexibleSUSY`HighScaleInput          = FlexibleSUSY`HighScaleInput          /. allIndexReplacementRules;
          FlexibleSUSY`LowScale                = FlexibleSUSY`LowScale                /. allIndexReplacementRules;
          FlexibleSUSY`LowScaleFirstGuess      = FlexibleSUSY`LowScaleFirstGuess      /. allIndexReplacementRules;
          FlexibleSUSY`LowScaleInput           = FlexibleSUSY`LowScaleInput           /. allIndexReplacementRules;
          FlexibleSUSY`SUSYScale               = FlexibleSUSY`SUSYScale               /. allIndexReplacementRules;
          FlexibleSUSY`SUSYScaleFirstGuess     = FlexibleSUSY`SUSYScaleFirstGuess     /. allIndexReplacementRules;
          FlexibleSUSY`SUSYScaleInput          = FlexibleSUSY`SUSYScaleInput          /. allIndexReplacementRules;
         ];

GUTNormalization[coupling_] :=
    Parameters`GetGUTNormalization[coupling];

GeneralReplacementRules[] :=
    { "@VectorZ@"     -> ToValidCSymbolString[SARAH`VectorZ],
      "@VectorP@"     -> ToValidCSymbolString[SARAH`VectorP],
      "@VectorW@"     -> ToValidCSymbolString[SARAH`VectorW],
      "@VectorG@"     -> ToValidCSymbolString[SARAH`VectorG],
      "@TopQuark@"    -> ToValidCSymbolString[SARAH`TopQuark],
      "@BottomQuark@" -> ToValidCSymbolString[SARAH`BottomQuark],
      "@Electron@"    -> ToValidCSymbolString[SARAH`Electron],
      "@Neutrino@"    -> ToValidCSymbolString[SARAH`Neutrino],
      "@HiggsBoson@"  -> ToValidCSymbolString[SARAH`HiggsBoson],
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
      "@InputParameter_" ~~ num_ ~~ "@" /; IntegerQ[ToExpression[num]] :> CConversion`ToValidCSymbolString[FlexibleSUSY`InputParameters[[ToExpression[num]]]],
      "@DateAndTime@"         -> DateString[]
    }


WriteRGEClass[betaFun_List, anomDim_List, files_List,
              additionalDecl_:"", numberOfBaseClassParameters_:0] :=
   Module[{beta, setter, getter, parameterDef, set,
           display, parameterDefaultInit,
           cCtorParameterList, parameterCopyInit, betaParameterList,
           anomDimPrototypes, anomDimFunctions, printParameters, parameters,
           numberOfParameters, clearParameters},
          (* extract list of parameters from the beta functions *)
          parameters = BetaFunction`GetName[#]& /@ betaFun;
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
          anomDimPrototypes    = AnomalousDimension`CreateAnomDimPrototypes[anomDim];
          anomDimFunctions     = AnomalousDimension`CreateAnomDimFunctions[anomDim];
          printParameters      = WriteOut`PrintParameters[parameters, "ostr"];
          WriteOut`ReplaceInFiles[files,
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
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteInputParameterClass[inputParameters_List, freePhases_List,
                         unfixedParameters_List,
                         defaultValues_List, files_List] :=
   Module[{defineInputParameters, defaultInputParametersInit},
          defineInputParameters = Constraint`DefineInputParameters[Join[inputParameters,freePhases,unfixedParameters]];
          defaultInputParametersInit = Constraint`InitializeInputParameters[Join[defaultValues,freePhases,unfixedParameters]];
          WriteOut`ReplaceInFiles[files,
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
          saveEwsbOutputParameters    = Parameters`SaveParameterLocally[FlexibleSUSY`EWSBOutputParameters, "old_", "MODELPARAMETER"];
          restoreEwsbOutputParameters = Parameters`RestoreParameter[FlexibleSUSY`EWSBOutputParameters, "old_", "model"];
          WriteOut`ReplaceInFiles[files,
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
          WriteOut`ReplaceInFiles[files,
                 { "@initialGuessAtLowScale@"  -> IndentText[WrapLines[initialGuessAtLowScale]],
                   "@initialGuessAtHighScale@" -> IndentText[WrapLines[initialGuessAtHighScale]],
                   "@setDRbarYukawaCouplings@" -> IndentText[WrapLines[setDRbarYukawaCouplings]],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteConvergenceTesterClass[particles_List, files_List] :=
   Module[{compareFunction},
          compareFunction = ConvergenceTester`CreateCompareFunction[particles];
          WriteOut`ReplaceInFiles[files,
                 { "@compareFunction@"      -> IndentText[WrapLines[compareFunction]],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

(* Returns a list of three-component list where the information is stored
   which vev corresponds to which CP even mass eigenstate.

   Example: MRSSM
   It[] := vevs = {vd,vu,vT,vS};
   It[] := CreateVEVsToFieldsAssociation[vevs]
   Out[] = {{vd, hh, 1}, {vu, hh, 2}, {vT, hh, 4}, {vS, hh, 3}}
 *)
CreateVEVsToFieldsAssociation[vevs_List] :=
    Module[{association = {}, v, phi, higgs},
           For[v = 1, v <= Length[vevs], v++,
               (* find CP even gauge-eigenstate Higgs for the vev *)
               phi = Cases[SARAH`DEFINITION[FlexibleSUSY`FSEigenstates][SARAH`VEVs],
                           {_, {vevs[[v]], _}, {__}, {p_,_}} :> p
                          ];
               If[Head[phi] =!= List || Length[phi] != 1,
                  Print["Error: could not find CP even Higgs field for vev ", vevs[[v]]];
                  Quit[1];
                 ];
               phi = phi[[1]];
               (* find position of phi in the CP even mass eigenstate vector *)
               higgs = Cases[SARAH`DEFINITION[FlexibleSUSY`FSEigenstates][SARAH`MatterSector],
                             {ps__ /; MemberQ[ps, phi], {h_,_}} :> {h, Position[ps, phi][[1,1]]}
                            ];
               If[Head[higgs] =!= List || Length[higgs] != 1,
                  Print["Error: could not find CP even Higgs field ", phi,
                        " in MatterSector definitions "];
                  Quit[1];
                 ];
               higgs = higgs[[1]];
               AppendTo[association, {vevs[[v]], higgs[[1]], higgs[[2]]}];
              ];
           Return[association];
          ];


WriteModelClass[massMatrices_List, vevs_List, ewsbEquations_List,
                parametersFixedByEWSB_List, ewsbSolution_List, freePhases_List,
                nPointFunctions_List, vertexRules_List, phases_List,
                enablePoleMassThreads_,
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
            phasesInit = "",
            loopMassesPrototypes = "", loopMassesFunctions = "",
            runningDRbarMassesPrototypes = "", runningDRbarMassesFunctions = "",
            callAllLoopMassFunctions = "", printMasses = "", printMixingMatrices = "",
            masses, mixingMatrices, oneLoopTadpoles,
            dependenceNumPrototypes, dependenceNumFunctions,
            clearOutputParameters = "", solveEwsbTreeLevel = "",
            saveEwsbOutputParameters, restoreEwsbOutputParameters,
            softScalarMasses, softHiggsMasses,
            saveSoftHiggsMasses, restoreSoftHiggsMasses,
            solveTreeLevelEWSBviaSoftHiggsMasses,
            copyDRbarMassesToPoleMasses = ""
           },
           If[Length[vevs] != Length[ewsbEquations],
              Print["Error: number of vevs != number of EWSB equations"];
              Print["   vevs = ", vevs];
              Print["   EWSB equations = ", ewsbEquations];
              Quit[1];
             ];
           For[k = 1, k <= Length[massMatrices], k++,
               massGetters          = massGetters <> TreeMasses`CreateMassGetter[massMatrices[[k]]];
               mixingMatrixGetters  = mixingMatrixGetters <> TreeMasses`CreateMixingMatrixGetter[massMatrices[[k]]];
               physicalMassesDef    = physicalMassesDef <> TreeMasses`CreatePhysicalMassDefinition[massMatrices[[k]]];
               mixingMatricesDef    = mixingMatricesDef <> TreeMasses`CreateMixingMatrixDefinition[massMatrices[[k]]];
               physicalMassesInit   = physicalMassesInit <> TreeMasses`CreatePhysicalMassInitialization[massMatrices[[k]]];
               physicalMassesInitNoLeadingComma = StringTrim[physicalMassesInit, StartOfString ~~ ","];
               mixingMatricesInit   = mixingMatricesInit <> TreeMasses`CreateMixingMatrixInitialization[massMatrices[[k]]];
               clearOutputParameters = clearOutputParameters <> TreeMasses`ClearOutputParameters[massMatrices[[k]]];
               copyDRbarMassesToPoleMasses = copyDRbarMassesToPoleMasses <> TreeMasses`CopyDRBarMassesToPoleMasses[massMatrices[[k]]];
               massCalculationPrototypes = massCalculationPrototypes <> TreeMasses`CreateMassCalculationPrototype[massMatrices[[k]]];
               massCalculationFunctions  = massCalculationFunctions  <> TreeMasses`CreateMassCalculationFunction[massMatrices[[k]]];
              ];
           calculateAllMasses = TreeMasses`CallMassCalculationFunctions[massMatrices];
           For[k = 1, k <= Length[ewsbEquations], k++,
               tadpoleEqPrototypes = tadpoleEqPrototypes <> EWSB`CreateEWSBEqPrototype[vevs[[k]]];
               tadpoleEqFunctions  = tadpoleEqFunctions  <> EWSB`CreateEWSBEqFunction[vevs[[k]], ewsbEquations[[k]]];
              ];
           If[Length[parametersFixedByEWSB] != numberOfEWSBEquations,
              Print["Error: There are ", numberOfEWSBEquations, " EWSB ",
                    "equations, but you want to fix ", Length[parametersFixedByEWSB],
                    " parameters: ", parametersFixedByEWSB];
             ];
           oneLoopTadpoles              = Cases[nPointFunctions, SelfEnergies`Tadpole[___]];
           calculateOneLoopTadpoles     = SelfEnergies`FillArrayWithOneLoopTadpoles[CreateVEVsToFieldsAssociation[vevs]];
           calculateTreeLevelTadpoles   = EWSB`FillArrayWithEWSBEqs[vevs, parametersFixedByEWSB, freePhases];
           ewsbInitialGuess             = EWSB`FillInitialGuessArray[parametersFixedByEWSB];
           solveEwsbTreeLevel           = EWSB`CreateTreeLevelEwsbSolver[ewsbSolution];
           {selfEnergyPrototypes, selfEnergyFunctions} = SelfEnergies`CreateNPointFunctions[nPointFunctions, vertexRules];
           phasesDefinition             = Phases`CreatePhasesDefinition[phases];
           phasesGetterSetters          = Phases`CreatePhasesGetterSetters[phases];
           phasesInit                   = Phases`CreatePhasesInitialization[phases];
           loopMassesPrototypes         = LoopMasses`CreateOneLoopPoleMassPrototypes[];
           (* If you want to add tadpoles, call the following routine like this:
              CreateOneLoopPoleMassFunctions[diagonalizationPrecision, oneLoopTadpoles, vevs];
              *)
           loopMassesFunctions          = LoopMasses`CreateOneLoopPoleMassFunctions[diagonalizationPrecision, {}, {}];
           runningDRbarMassesPrototypes = LoopMasses`CreateRunningDRbarMassPrototypes[];
           runningDRbarMassesFunctions  = LoopMasses`CreateRunningDRbarMassFunctions[];
           callAllLoopMassFunctions     = LoopMasses`CallAllOneLoopPoleMassFunctions[FlexibleSUSY`FSEigenstates, enablePoleMassThreads];
           masses                       = FlexibleSUSY`M[TreeMasses`GetMassEigenstate[#]]& /@ massMatrices;
           printMasses                  = WriteOut`PrintParameters[masses, "ostr"];
           mixingMatrices               = Flatten[TreeMasses`GetMixingMatrixSymbol[#]& /@ massMatrices];
           printMixingMatrices          = WriteOut`PrintParameters[mixingMatrices, "ostr"];
           dependenceNumPrototypes      = TreeMasses`CreateDependenceNumPrototypes[];
           dependenceNumFunctions       = TreeMasses`CreateDependenceNumFunctions[];
           saveEwsbOutputParameters     = Parameters`SaveParameterLocally[FlexibleSUSY`EWSBOutputParameters, "one_loop_", ""];
           restoreEwsbOutputParameters  = Parameters`RestoreParameter[FlexibleSUSY`EWSBOutputParameters, "one_loop_", ""];
           If[Head[SARAH`ListSoftBreakingScalarMasses] === List,
              softScalarMasses          = DeleteDuplicates[SARAH`ListSoftBreakingScalarMasses];,
              Print["Error: no soft breaking scalar masses found!"];
              softScalarMasses          = {};
             ];
           softHiggsMasses              = Select[softScalarMasses, (!FreeQ[ewsbEquations, #])&];
           saveSoftHiggsMasses          = Parameters`SaveParameterLocally[softHiggsMasses, "old_", ""];
           restoreSoftHiggsMasses       = Parameters`RestoreParameter[softHiggsMasses, "old_", ""];
           solveTreeLevelEWSBviaSoftHiggsMasses = EWSB`SolveTreeLevelEwsbVia[ewsbEquations, softHiggsMasses];
           WriteOut`ReplaceInFiles[files,
                          { "@massGetters@"          -> IndentText[massGetters],
                            "@mixingMatrixGetters@"  -> IndentText[mixingMatrixGetters],
                            "@tadpoleEqPrototypes@"  -> IndentText[tadpoleEqPrototypes],
                            "@tadpoleEqFunctions@"   -> tadpoleEqFunctions,
                            "@numberOfEWSBEquations@"-> ToString[numberOfEWSBEquations],
                            "@calculateTreeLevelTadpoles@" -> IndentText[calculateTreeLevelTadpoles],
                            "@calculateOneLoopTadpoles@"   -> IndentText[calculateOneLoopTadpoles],
                            "@clearOutputParameters@"  -> IndentText[clearOutputParameters],
                            "@copyDRbarMassesToPoleMasses@" -> IndentText[copyDRbarMassesToPoleMasses],
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
           WriteOut`ReplaceInFiles[files,
                          { Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WritePlotScripts[files_List] :=
    Module[{},
           WriteOut`ReplaceInFiles[files,
                          { Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteUtilitiesClass[massMatrices_List, betaFun_List, minpar_List, extpar_List,
                    unfixedParameters_List, files_List] :=
    Module[{k, particles, susyParticles, smParticles,
            fillSpectrumVectorWithSusyParticles = "",
            fillSpectrumVectorWithSMParticles = "",
            particleLaTeXNames = "",
            particleNames = "", particleEnum = "", particleMultiplicity = "",
            parameterNames = "", parameterEnum = "", numberOfParameters = 0,
            fillInputParametersFromMINPAR = "", fillInputParametersFromEXTPAR = "",
            writeSLHAMassBlock = "", writeSLHAMixingMatricesBlocks = "",
            writeSLHAModelParametersBlocks = "", writeSLHAMinparBlock = "",
            writeSLHAExtparBlock = "", readUnfixedParameters},
           particles = GetMassEigenstate /@ massMatrices;
           susyParticles = Select[particles, (!SARAH`SMQ[#])&];
           smParticles   = Complement[particles, susyParticles];
           particleEnum       = TreeMasses`CreateParticleEnum[particles];
           particleMultiplicity = TreeMasses`CreateParticleMultiplicity[particles];
           particleNames      = TreeMasses`CreateParticleNames[particles];
           particleLaTeXNames = TreeMasses`CreateParticleLaTeXNames[particles];
           fillSpectrumVectorWithSusyParticles = TreeMasses`FillSpectrumVector[susyParticles];
           fillSpectrumVectorWithSMParticles   = TreeMasses`FillSpectrumVector[smParticles];
           numberOfParameters = BetaFunction`CountNumberOfParameters[betaFun];
           parameterEnum      = BetaFunction`CreateParameterEnum[betaFun];
           parameterNames     = BetaFunction`CreateParameterNames[betaFun];
           fillInputParametersFromMINPAR = Parameters`FillInputParametersFromTuples[minpar];
           fillInputParametersFromEXTPAR = Parameters`FillInputParametersFromTuples[extpar];
           readUnfixedParameters         = WriteOut`ReadUnfixedParameters[unfixedParameters];
           writeSLHAMassBlock = WriteOut`WriteSLHAMassBlock[massMatrices];
           writeSLHAMixingMatricesBlocks  = WriteOut`WriteSLHAMixingMatricesBlocks[];
           writeSLHAModelParametersBlocks = WriteOut`WriteSLHAModelParametersBlocks[];
           writeSLHAMinparBlock = WriteOut`WriteSLHAMinparBlock[minpar];
           writeSLHAExtparBlock = WriteOut`WriteSLHAExtparBlock[extpar];
           WriteOut`ReplaceInFiles[files,
                          { "@fillSpectrumVectorWithSusyParticles@" -> IndentText[fillSpectrumVectorWithSusyParticles],
                            "@fillSpectrumVectorWithSMParticles@"   -> IndentText[IndentText[fillSpectrumVectorWithSMParticles]],
                            "@particleEnum@"       -> IndentText[WrapLines[particleEnum]],
                            "@particleMultiplicity@" -> IndentText[WrapLines[particleMultiplicity]],
                            "@particleNames@"      -> IndentText[WrapLines[particleNames]],
                            "@particleLaTeXNames@" -> IndentText[WrapLines[particleLaTeXNames]],
                            "@parameterEnum@"     -> IndentText[WrapLines[parameterEnum]],
                            "@parameterNames@"     -> IndentText[WrapLines[parameterNames]],
                            "@fillInputParametersFromMINPAR@" -> IndentText[fillInputParametersFromMINPAR],
                            "@fillInputParametersFromEXTPAR@" -> IndentText[fillInputParametersFromEXTPAR],
                            "@readUnfixedParameters@"         -> IndentText[readUnfixedParameters],
                            "@writeSLHAMassBlock@" -> IndentText[writeSLHAMassBlock],
                            "@writeSLHAMixingMatricesBlocks@"  -> IndentText[writeSLHAMixingMatricesBlocks],
                            "@writeSLHAModelParametersBlocks@" -> IndentText[writeSLHAModelParametersBlocks],
                            "@writeSLHAMinparBlock@"           -> IndentText[writeSLHAMinparBlock],
                            "@writeSLHAExtparBlock@"           -> IndentText[writeSLHAExtparBlock],
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
    And @@ (FileExists[path,#]& /@ fileNames);

FilesExist[fileNames_List] :=
    And @@ (FileExists /@ fileNames);

LatestModificationTimeInSeconds[file_String] :=
    If[FileExists[file],
       AbsoluteTime[FileDate[file, "Modification"]], 0];

LatestModificationTimeInSeconds[files_List] :=
    Max[LatestModificationTimeInSeconds /@ files];

SARAHModelFileModificationTimeInSeconds[] :=
    Module[{files},
           files = Join[{SARAH`ModelFile},
                        FileNameJoin[{$sarahCurrentModelDir, #}]& /@ {"parameters.m", "particles.m"}
                       ];
           Return[LatestModificationTimeInSeconds[files]];
          ];

GetRGEFileNames[outputDir_String] :=
    Module[{rgeDir, fileNames},
           rgeDir = FileNameJoin[{outputDir, "RGEs"}];
           fileNames = { "BetaYijk.m", "BetaGauge.m", "BetaWijkl.m",
                         "BetaMuij.m", "BetaLi.m", "BetaQijkl.m",
                         "BetaTijk.m", "BetaBij.m", "BetaLSi.m",
                         "Betam2ij.m", "BetaMi.m", "BetaVEV.m" };
           If[SARAH`AddDiracGauginos === True,
              AppendTo[fileNames, "BetaDGi.m"];
             ];
           FileNameJoin[{rgeDir, #}]& /@ fileNames
          ];

RGEFilesExist[outputDir_String] :=
    FilesExist[GetRGEFileNames[outputDir]];

RGEsModificationTimeInSeconds[outputDir_String] :=
    LatestModificationTimeInSeconds[GetRGEFileNames[outputDir]];

GetSelfEnergyFileNames[outputDir_String, eigenstates_] :=
    FileNameJoin[{outputDir, ToString[eigenstates],
                  "One-Loop", "SelfEnergy.m"}];

SelfEnergyFilesExist[outputDir_String, eigenstates_] :=
    FileExists[GetSelfEnergyFileNames[outputDir, eigenstates]];

SelfEnergyFilesModificationTimeInSeconds[outputDir_String, eigenstates_] :=
    LatestModificationTimeInSeconds[GetSelfEnergyFileNames[outputDir, eigenstates]];

NeedToCalculateSelfEnergies[eigenstates_] :=
    Module[{seFilesExist, seFilesTimeStamp, sarahModelFileTimeStamp,
            needToCalculateSEs},
           seFilesExist = SelfEnergyFilesExist[$sarahCurrentOutputMainDir, eigenstates];
           seFilesTimeStamp = SelfEnergyFilesModificationTimeInSeconds[$sarahCurrentOutputMainDir, eigenstates];
           sarahModelFileTimeStamp = SARAHModelFileModificationTimeInSeconds[];
           needToCalculateSEs = Or[!seFilesExist,
                                    seFilesExist && (sarahModelFileTimeStamp > seFilesTimeStamp)];
           If[!seFilesExist,
              Print["Self-energies have not been calculated yet, calculating them ..."];
             ];
           If[seFilesExist && (sarahModelFileTimeStamp > seFilesTimeStamp),
              Print["SARAH model files are newer than self-energy files, recalculating them ..."];
             ];
           Return[needToCalculateSEs];
          ];

GetTadpoleFileName[outputDir_String, eigenstates_] :=
    FileNameJoin[{outputDir, ToString[eigenstates],
                  "One-Loop", "Tadpoles1Loop.m"}];

TadpoleFileExists[outputDir_String, eigenstates_] :=
    FileExists[GetTadpoleFileName[outputDir, eigenstates]];

TadpoleFilesModificationTimeInSeconds[outputDir_String, eigenstates_] :=
    LatestModificationTimeInSeconds[GetTadpoleFileName[outputDir, eigenstates]];

NeedToCalculateTadpoles[eigenstates_] :=
    Module[{tadpoleFilesExist, tadpoleFilesTimeStamp, sarahModelFileTimeStamp,
            needToCalculateTadpoles},
           tadpoleFilesExist = TadpoleFileExists[$sarahCurrentOutputMainDir, eigenstates];
           tadpoleFilesTimeStamp = TadpoleFilesModificationTimeInSeconds[$sarahCurrentOutputMainDir, eigenstates];
           sarahModelFileTimeStamp = SARAHModelFileModificationTimeInSeconds[];
           needToCalculateTadpoles = Or[!tadpoleFilesExist,
                                        tadpoleFilesExist && (sarahModelFileTimeStamp > tadpoleFilesTimeStamp)];
           If[!tadpoleFilesExist,
              Print["Tadpoles have not been calculated yet, calculating them ..."];
             ];
           If[tadpoleFilesExist && (sarahModelFileTimeStamp > tadpoleFilesTimeStamp),
              Print["SARAH model files are newer than tadpoles file, recalculating them ..."];
             ];
           Return[needToCalculateTadpoles];
          ];

GetUnrotatedParticlesFileName[outputDir_String, eigenstates_] :=
    FileNameJoin[{outputDir, ToString[eigenstates],
                  "One-Loop", "UnrotatedParticles.m"}];

UnrotatedParticlesFilesExist[outputDir_String, eigenstates_] :=
    FileExists[GetUnrotatedParticlesFileName[outputDir, eigenstates]];

UnrotatedParticlesFilesModificationTimeInSeconds[outputDir_String, eigenstates_] :=
    LatestModificationTimeInSeconds[GetUnrotatedParticlesFileName[outputDir, eigenstates]];

NeedToCalculateUnrotatedParticles[eigenstates_] :=
    Module[{unrotatedParticlesFilesExist, unrotatedParticlesFilesTimeStamp, sarahModelFileTimeStamp,
            needToCalculateSEs},
           unrotatedParticlesFilesExist = UnrotatedParticlesFilesExist[$sarahCurrentOutputMainDir, eigenstates];
           unrotatedParticlesFilesTimeStamp = UnrotatedParticlesFilesModificationTimeInSeconds[$sarahCurrentOutputMainDir, eigenstates];
           sarahModelFileTimeStamp = SARAHModelFileModificationTimeInSeconds[];
           needToCalculateSEs = Or[!unrotatedParticlesFilesExist,
                                    unrotatedParticlesFilesExist && (sarahModelFileTimeStamp > unrotatedParticlesFilesTimeStamp)];
           If[!unrotatedParticlesFilesExist,
              Print["Unrotated particles have not been calculated yet, calculating them ..."];
             ];
           If[unrotatedParticlesFilesExist && (sarahModelFileTimeStamp > unrotatedParticlesFilesTimeStamp),
              Print["SARAH model files are newer than unrotated particles file, recalculating them ..."];
             ];
           Return[needToCalculateSEs];
          ];

SearchSelfEnergies[outputDir_String, eigenstates_] :=
    Module[{fileName},
           fileName = GetSelfEnergyFileNames[outputDir, eigenstates];
           If[FileExists[fileName], fileName, ""]
          ];

SearchUnrotatedParticles[outputDir_String, eigenstates_] :=
    Module[{fileName},
           fileName = GetUnrotatedParticlesFileName[outputDir, eigenstates];
           If[FileExists[fileName], fileName, ""]
          ];

SearchTadpoles[outputDir_String, eigenstates_] :=
    Module[{fileName},
           fileName = GetTadpoleFileName[outputDir, eigenstates];
           If[FileExists[fileName], fileName, ""]
          ];

NeedToCalculateRGEs[] :=
    Module[{rgeFilesExist, rgeFilesTimeStamp, sarahModelFileTimeStamp,
            needToCalculateRGEs},
           rgeFilesExist = RGEFilesExist[$sarahCurrentOutputMainDir];
           rgeFilesTimeStamp = RGEsModificationTimeInSeconds[$sarahCurrentOutputMainDir];
           sarahModelFileTimeStamp = SARAHModelFileModificationTimeInSeconds[];
           needToCalculateRGEs = Or[!rgeFilesExist,
                                    rgeFilesExist && (sarahModelFileTimeStamp > rgeFilesTimeStamp)];
           If[!rgeFilesExist,
              Print["RGEs have not been calculated yet, calculating them ..."];
             ];
           If[rgeFilesExist && (sarahModelFileTimeStamp > rgeFilesTimeStamp),
              Print["SARAH model files are newer than RGE files, recalculating them ..."];
             ];
           If[!needToCalculateRGEs,
              Print["Reading RGEs from files."];
             ];
           Return[needToCalculateRGEs];
          ];

FSPrepareRGEs[] :=
    Module[{needToCalculateRGEs, betas},
           needToCalculateRGEs = NeedToCalculateRGEs[];
           SARAH`CalcRGEs[ReadLists -> !needToCalculateRGEs,
                          TwoLoop -> True,
                          NoMatrixMultiplication -> False];
           (* check if the beta functions were calculated correctly *)
           betas = { SARAH`BetaWijkl, SARAH`BetaYijk, SARAH`BetaMuij,
                     SARAH`BetaLi, SARAH`BetaGauge, SARAH`BetaVEV,
                     SARAH`BetaQijkl, SARAH`BetaTijk, SARAH`BetaBij,
                     SARAH`BetaLSi, SARAH`Betam2ij, SARAH`BetaMi,
                     SARAH`BetaDGi };
           If[Head[#] === Symbol && !ValueQ[#], Set[#,{}]]& /@ betas;
           If[!ValueQ[SARAH`Gij] || Head[SARAH`Gij] =!= List,
              SARAH`Gij = {};
             ];
          ];

FSCheckLoopCorrections[eigenstates_] :=
    Module[{needToCalculateLoopCorrections},
           needToCalculateLoopCorrections = Or[
               NeedToCalculateSelfEnergies[eigenstates],
               NeedToCalculateTadpoles[eigenstates],
               NeedToCalculateUnrotatedParticles[eigenstates]
                                              ];
           If[needToCalculateLoopCorrections,
              SARAH`CalcLoopCorrections[eigenstates];
             ];
          ];

PrepareSelfEnergies[eigenstates_] :=
    Module[{selfEnergies = {}, selfEnergiesFile},
           selfEnergiesFile = SearchSelfEnergies[$sarahCurrentOutputMainDir, eigenstates];
           If[selfEnergiesFile == "",
              Print["Error: self-energy files not found: ", selfEnergiesFile];
              Quit[1];
             ];
           Print["Reading self-energies from file ", selfEnergiesFile, " ..."];
           selfEnergies = Get[selfEnergiesFile];
           Print["Converting self-energies ..."];
           ConvertSarahSelfEnergies[selfEnergies]
          ];

PrepareTadpoles[eigenstates_] :=
    Module[{tadpoles = {}, tadpolesFile},
           tadpolesFile = SearchTadpoles[$sarahCurrentOutputMainDir, eigenstates];
           If[tadpolesFile == "",
              Print["Error: tadpole file not found: ", tadpolesFile];
              Quit[1];
             ];
           Print["Reading tadpoles from file ", tadpolesFile, " ..."];
           tadpoles = Get[tadpolesFile];
           Print["Converting tadpoles ..."];
           ConvertSarahTadpoles[tadpoles]
          ];

PrepareUnrotatedParticles[eigenstates_] :=
    Module[{nonMixedParticles = {}, nonMixedParticlesFile},
           nonMixedParticlesFile = SearchUnrotatedParticles[$sarahCurrentOutputMainDir, eigenstates];
           If[nonMixedParticlesFile == "",
              Print["Error: file with unrotated fields not found: ", tadpolesFile];
              Quit[1];
             ];
           Print["Reading unrotated particles from file ", nonMixedParticlesFile, " ..."];
           nonMixedParticles = Get[nonMixedParticlesFile];
           TreeMasses`SetUnrotatedParticles[nonMixedParticles];
          ];

ReadDiagonalizationPrecisions[defaultPrecision_Symbol, highPrecisionList_List,
                              mediumPrecisionList_List, lowPrecisionList_List, eigenstates_] :=
    Module[{particles, particle, i, precisionList = {}},
           If[!MemberQ[{LowPrecision, MediumPrecision, HighPrecision}, defaultPrecision],
              Print["Error: ", defaultPrecision, " is not a valid",
                    " diagonalization precision!"];
              Print["   Available are: LowPrecision, MediumPrecision, HighPrecision"];
              Quit[1];
             ];
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
              Quit[1];
             ];
          ];

FindUnfixedParameters[fixed_List] :=
    Module[{fixedParameters},
           fixedParameters = DeleteDuplicates[Flatten[Join[fixed,
                                          { SARAH`hyperchargeCoupling, SARAH`leftCoupling,
                                            SARAH`strongCoupling, SARAH`UpYukawa, SARAH`DownYukawa,
                                            SARAH`ElectronYukawa }]]];
           Complement[allParameters, fixedParameters]
          ];

Options[MakeFlexibleSUSY] :=
    {
        Eigenstates -> SARAH`EWSB,
        InputFile -> "FlexibleSUSY.m",
        DefaultDiagonalizationPrecision -> MediumPrecision,
        HighDiagonalizationPrecision -> {},
        MediumDiagonalizationPrecision -> {},
        LowDiagonalizationPrecision -> {},
        EnablePoleMassThreads -> True,
        SolveEWSBTimeConstraint -> 120 (* in seconds *)
    };

MakeFlexibleSUSY[OptionsPattern[]] :=
    Module[{nPointFunctions, runInputFile, initialGuesserInputFile,
            susyBetaFunctions, susyBreakingBetaFunctions,
            numberOfSusyParameters, anomDim,
            ewsbEquations, massMatrices, phases, vevs,
            diagonalizationPrecision, allParticles, freePhases, ewsbSolution,
            fixedParameters, treeLevelEwsbOutputFile,
	    vertexRules,
	    Lat$massMatrices},
           (* check if SARAH`Start[] was called *)
           If[!ValueQ[Model`Name],
              Print["Error: Model`Name is not defined.  Did you call SARAH`Start[\"Model\"]?"];
              Quit[1];
             ];
           CheckSARAHVersion[];
           FSEigenstates = OptionValue[Eigenstates];
           FSSolveEWSBTimeConstraint = OptionValue[SolveEWSBTimeConstraint];
           (* load model file *)
           LoadModelFile[OptionValue[InputFile]];
           Print["FlexibleSUSY model file loaded"];
           Print["  Model: ", FlexibleSUSY`FSModelName];
           Print["  Model file: ", OptionValue[InputFile]];
           Print["  Model output directory: ", Global`$flexiblesusyOutputDir];

           PrintHeadline["Reading SARAH output files"];
           (* get RGEs *)
           FSPrepareRGEs[];
           FSCheckLoopCorrections[FSEigenstates];
           nPointFunctions = Join[PrepareSelfEnergies[FSEigenstates], PrepareTadpoles[FSEigenstates]];
           PrepareUnrotatedParticles[FSEigenstates];
           (* adapt SARAH`Conj to our needs *)
           (* Clear[Conj]; *)
           SARAH`Conj[(B_)[b__]] = .;
           SARAH`Conj /: SARAH`Conj[SARAH`Conj[x_]] := x;
           RXi[_] = 1;
           SARAH`Xi = 1;
           SARAH`Xip = 1;
           SARAH`rMS = 0;

           FlexibleSUSY`InputParameters = Join[(#[[2]])& /@ SARAH`MINPAR, (#[[2]])& /@ SARAH`EXTPAR];
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
                                         SARAH`BetaLSi  ,
                                         SARAH`Betam2ij ,
                                         SARAH`BetaMi   ,
                                         SARAH`BetaDGi  };

           susyBetaFunctions = BetaFunction`ConvertSarahRGEs[susyBetaFunctions];
           susyBetaFunctions = Select[susyBetaFunctions, (BetaFunction`GetAllBetaFunctions[#]!={})&];
           Parameters`AddRealParameter[(BetaFunction`GetName /@ susyBetaFunctions) /.
                                       a_[Susyno`LieGroups`i1] :> a /.
                                       a_[Susyno`LieGroups`i1,SARAH`i2] :> a];

           numberOfSusyParameters = BetaFunction`CountNumberOfParameters[susyBetaFunctions];
           anomDim = AnomalousDimension`ConvertSarahAnomDim[SARAH`Gij];

           PrintHeadline["Creating model parameter classes"];
           Print["Creating class for susy parameters ..."];
           WriteRGEClass[susyBetaFunctions, anomDim,
                         {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_susy_parameters.hpp.in"}],
                           FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_susy_parameters.hpp"}]},
                          {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_susy_parameters.cpp.in"}],
                           FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_susy_parameters.cpp"}]}}];

           susyBreakingBetaFunctions = ConvertSarahRGEs[susyBreakingBetaFunctions];
           susyBreakingBetaFunctions = Select[susyBreakingBetaFunctions, (BetaFunction`GetAllBetaFunctions[#]!={})&];
           Parameters`AddRealParameter[(BetaFunction`GetName /@ susyBreakingBetaFunctions) /.
                                       a_[Susyno`LieGroups`i1] :> a /.
                                       a_[Susyno`LieGroups`i1,SARAH`i2] :> a];

           allBetaFunctions = Join[susyBetaFunctions, susyBreakingBetaFunctions];

           {traceDecl, traceRules} = CreateTraceAbbr[SARAH`TraceAbbr];
           susyBreakingBetaFunctions = susyBreakingBetaFunctions /. traceRules;

           (* store all model parameters *)
           allParameters = Join[BetaFunction`GetName /@ susyBetaFunctions,
                                BetaFunction`GetName /@ susyBreakingBetaFunctions] /.
                               a_[Susyno`LieGroups`i1] :> a /.
                               a_[Susyno`LieGroups`i1,SARAH`i2] :> a;
           allIndexReplacementRules = Parameters`CreateIndexReplacementRules[allParameters];
           Parameters`SetModelParameters[allParameters];
           FlexibleSUSY`FSLesHouchesList = SA`LHList;

           (* search for unfixed parameters *)
           fixedParameters = Join[FlexibleSUSY`EWSBOutputParameters,
                                  Constraint`FindFixedParametersFromConstraint[FlexibleSUSY`LowScaleInput],
                                  Constraint`FindFixedParametersFromConstraint[FlexibleSUSY`SUSYScaleInput],
                                  Constraint`FindFixedParametersFromConstraint[FlexibleSUSY`HighScaleInput]
                                 ];
           FlexibleSUSY`FSUnfixedParameters = FindUnfixedParameters[fixedParameters];
           If[FlexibleSUSY`FSUnfixedParameters =!= {} &&
              FlexibleSUSY`OnlyLowEnergyFlexibleSUSY =!= True,
              Print["Warning: the following parameters are not fixed by any constraint:"];
              Print["  ", FlexibleSUSY`FSUnfixedParameters];
             ];
           (* adding the types and their input names to the parameters *)
           FlexibleSUSY`FSUnfixedParameters = Select[Join[{BetaFunction`GetName[#], Symbol[ToValidCSymbolString[BetaFunction`GetName[#]] <> "Input"], #[[2]]}& /@ susyBetaFunctions,
                                                          {BetaFunction`GetName[#], Symbol[ToValidCSymbolString[BetaFunction`GetName[#]] <> "Input"], #[[2]]}& /@ susyBreakingBetaFunctions] /.
                                                     a_[Susyno`LieGroups`i1] :> a /.
                                                     a_[Susyno`LieGroups`i1,SARAH`i2] :> a,
                                                     MemberQ[FlexibleSUSY`FSUnfixedParameters,#[[1]]]&];
           (* add the unfixed parameters to the susy scale constraint *)
           If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY === True,
              FlexibleSUSY`SUSYScaleInput = Join[FlexibleSUSY`SUSYScaleInput,
                                                 {#[[1]],#[[2]]}& /@ FlexibleSUSY`FSUnfixedParameters];
              Parameters`SetInputParameters[Join[FlexibleSUSY`InputParameters,
                                                 (#[[2]])& /@ FlexibleSUSY`FSUnfixedParameters]];
             ];

           (* replace all indices in the user-defined model file variables *)
           ReplaceIndicesInUserInput[];

           numberOfSusyBreakingParameters = BetaFunction`CountNumberOfParameters[susyBreakingBetaFunctions];
           numberOfModelParameters = numberOfSusyParameters + numberOfSusyBreakingParameters;

           Print["Creating class for soft parameters ..."];
           WriteRGEClass[susyBreakingBetaFunctions, {},
                         {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_soft_parameters.hpp.in"}],
                           FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_soft_parameters.hpp"}]},
                          {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_soft_parameters.cpp.in"}],
                           FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_soft_parameters.cpp"}]}},
                         traceDecl, numberOfSusyParameters];

           vevs = #[[1]]& /@ SARAH`BetaVEV;
           ewsbEquations = SARAH`TadpoleEquations[FSEigenstates] /.
                           Parameters`ApplyGUTNormalization[] /.
                           SARAH`sum[idx_, start_, stop_, expr_] :> Sum[expr, {idx,start,stop}];
           If[Head[ewsbEquations] =!= List,
              Print["Error: Could not find EWSB equations for eigenstates ",
                    FSEigenstates];
              Quit[1];
             ];
           If[Length[vevs] =!= Length[ewsbEquations],
              Print["Error: There are ", Length[ewsbEquations],
                    " EWSB equations but ", Length[vevs], " VEVs"];
              ewsbEquations = {};
              vevs = {};
              FlexibleSUSY`EWSBOutputParameters = {};
              Quit[1];
             ];

           If[FlexibleSUSY`TreeLevelEWSBSolution === {},
              (* trying to find an analytic solution for the EWSB eqs. *)
              treeLevelEwsbOutputFile = FileNameJoin[{Global`$flexiblesusyOutputDir,
                                                      FlexibleSUSY`FSModelName <> "_tree_level_EWSB_solution.m"}];
              Print["Solving EWSB equations ..."];
              {ewsbSolution, freePhases} = EWSB`FindSolutionAndFreePhases[ewsbEquations,
                                                                          FlexibleSUSY`EWSBOutputParameters,
                                                                          treeLevelEwsbOutputFile];
              If[ewsbSolution === {},
                 Print["Warning: could not find an analytic solution to the EWSB eqs."];
                 Print["   An iterative algorithm will be used.  You can try to set"];
                 Print["   the solution by hand in the model file like this:"];
                 Print[""];
                 Print["   TreeLevelEWSBSolution = {"];
                 For[i = 1, i <= Length[FlexibleSUSY`EWSBOutputParameters], i++,
                     Print["      { ", FlexibleSUSY`EWSBOutputParameters[[i]], ", ... }" <>
                           If[i != Length[FlexibleSUSY`EWSBOutputParameters], ",", ""]];
                    ];
                 Print["   };\n"];
                 Print["   The tree-level EWSB solution was written to the file:"];
                 Print["      ", treeLevelEwsbOutputFile];
                ];
              ,
              If[Length[FlexibleSUSY`TreeLevelEWSBSolution] != Length[ewsbEquations],
                 Print["Error: not enough EWSB solutions given!"];
                 Quit[1];
                ];
              If[Sort[#[[1]]& /@ FlexibleSUSY`TreeLevelEWSBSolution] =!= Sort[FlexibleSUSY`EWSBOutputParameters],
                 Print["Error: Parameters given in TreeLevelEWSBSolution, do not match"];
                 Print["   the Parameters given in FlexibleSUSY`EWSBOutputParameters!"];
                 Quit[1];
                ];
              Print["Using user-defined EWSB eqs. solution"];
              freePhases = {};
              ewsbSolution = FlexibleSUSY`TreeLevelEWSBSolution;
             ];
           If[freePhases =!= {},
              Print["Note: adding free phases: ", freePhases];
             ];

           Print["Creating class for input parameters ..."];
           WriteInputParameterClass[FlexibleSUSY`InputParameters, Complement[freePhases, FlexibleSUSY`InputParameters],
                                    If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY =!= True, {},
                                       {#[[2]], #[[3]]}& /@ FlexibleSUSY`FSUnfixedParameters],
                                    FlexibleSUSY`DefaultParameterPoint,
                                    {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "input_parameters.hpp.in"}],
                                      FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_input_parameters.hpp"}]}}
                                   ];

	   On[Assert];

           Lat$massMatrices = ConvertSarahMassMatrices[] /.
                          Parameters`ApplyGUTNormalization[] /.
                          { SARAH`sum[j_, start_, end_, expr_] :> (Sum[expr, {j,start,end}]) };
           massMatrices = Lat$massMatrices /. allIndexReplacementRules;
	   Lat$massMatrices = LatticeUtils`FixDiagonalization[Lat$massMatrices];

           allParticles = FlexibleSUSY`M[GetMassEigenstate[#]]& /@ massMatrices;
           allOutputParameters = DeleteCases[DeleteDuplicates[
               Join[allParticles,
                    Flatten[GetMixingMatrixSymbol[#]& /@ massMatrices]]], Null];

           Parameters`SetOutputParameters[allOutputParameters];

           PrintHeadline["Creating utilities"];
           Print["Creating class for convergence tester ..."];
           WriteConvergenceTesterClass[allParticles,
               {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "convergence_tester.hpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_convergence_tester.hpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_convergence_tester.hpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_convergence_tester.hpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_convergence_tester.cpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_convergence_tester.cpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "lattice_convergence_tester.hpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_lattice_convergence_tester.hpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "lattice_convergence_tester.cpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_lattice_convergence_tester.cpp"}]}
               }
                                      ];

           Print["Creating utilities class ..."];
           WriteUtilitiesClass[massMatrices, Join[susyBetaFunctions, susyBreakingBetaFunctions],
                               MINPAR, EXTPAR,
                               If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY =!= True, {},
                                  FlexibleSUSY`FSUnfixedParameters],
               {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "info.hpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_info.hpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "info.cpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_info.cpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "utilities.hpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_utilities.hpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "utilities.cpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_utilities.cpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "slha_io.hpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_slha_io.hpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "slha_io.cpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_slha_io.cpp"}]}
               }
                              ];

           Print["Creating plot scripts ..."];
           WritePlotScripts[{{FileNameJoin[{Global`$flexiblesusyTemplateDir, "plot_spectrum.gnuplot.in"}],
                              FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_plot_spectrum.gnuplot"}]},
                             {FileNameJoin[{Global`$flexiblesusyTemplateDir, "plot_rge_running.gnuplot.in"}],
                              FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_plot_rge_running.gnuplot"}]}}
                           ];

           PrintHeadline["Creating constraints"];
           Print["Creating class for high-scale constraint ..."];
           WriteConstraintClass[FlexibleSUSY`HighScale,
                                FlexibleSUSY`HighScaleInput,
                                FlexibleSUSY`HighScaleFirstGuess,
                                {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "high_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_high_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_high_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_high_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_high_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_high_scale_constraint.cpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "lattice_high_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_lattice_high_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "lattice_high_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_lattice_high_scale_constraint.cpp"}]}
                                }
                               ];

           Print["Creating class for susy-scale constraint ..."];
           WriteConstraintClass[FlexibleSUSY`SUSYScale,
                                FlexibleSUSY`SUSYScaleInput,
                                FlexibleSUSY`SUSYScaleFirstGuess,
                                {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "susy_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_susy_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_susy_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_susy_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_susy_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_susy_scale_constraint.cpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "lattice_susy_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_lattice_susy_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "lattice_susy_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_lattice_susy_scale_constraint.cpp"}]}
                                }
                               ];

           Print["Creating class for low-scale constraint ..."];
           WriteConstraintClass[FlexibleSUSY`LowScale,
                                FlexibleSUSY`LowScaleInput,
                                FlexibleSUSY`LowScaleFirstGuess,
                                {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "low_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_low_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_low_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_low_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_low_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_low_scale_constraint.cpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "lattice_low_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_lattice_low_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "lattice_low_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_lattice_low_scale_constraint.cpp"}]}
                                }
                               ];

           Print["Creating class for initial guesser ..."];
           initialGuesserInputFile = "initial_guesser";
           If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY,
              initialGuesserInputFile = "initial_guesser_low_scale_model";
             ];
           WriteInitialGuesserClass[FlexibleSUSY`InitialGuessAtLowScale,
                                    FlexibleSUSY`InitialGuessAtHighScale,
                                    {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "initial_guesser.hpp.in"}],
                                      FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_initial_guesser.hpp"}]},
                                     {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_" <> initialGuesserInputFile <> ".hpp.in"}],
                                      FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_initial_guesser.hpp"}]},
                                     {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_" <> initialGuesserInputFile <> ".cpp.in"}],
                                      FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_initial_guesser.cpp"}]},
                                     {FileNameJoin[{Global`$flexiblesusyTemplateDir, "lattice_" <> initialGuesserInputFile <> ".hpp.in"}],
                                      FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_lattice_initial_guesser.hpp"}]},
                                     {FileNameJoin[{Global`$flexiblesusyTemplateDir, "lattice_" <> initialGuesserInputFile <> ".cpp.in"}],
                                      FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_lattice_initial_guesser.cpp"}]}
                                    }
                                   ];

           phases = ConvertSarahPhases[SARAH`ParticlePhases];

           (* determin diagonalization precision for each particle *)
           diagonalizationPrecision = ReadDiagonalizationPrecisions[
               OptionValue[DefaultDiagonalizationPrecision],
               Flatten[{OptionValue[HighDiagonalizationPrecision]}],
               Flatten[{OptionValue[MediumDiagonalizationPrecision]}],
               Flatten[{OptionValue[LowDiagonalizationPrecision]}],
               FSEigenstates];

	   vertexRules = Vertices`VertexRules[nPointFunctions, Lat$massMatrices];

           PrintHeadline["Creating model"];
           Print["Creating class for model ..."];
           WriteModelClass[massMatrices, vevs, ewsbEquations,
                           FlexibleSUSY`EWSBOutputParameters, ewsbSolution, freePhases,
                           nPointFunctions, vertexRules, phases, OptionValue[EnablePoleMassThreads],
                           {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "model.hpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_model.hpp"}]},
                            {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_model.hpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_model.hpp"}]},
                            {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_model.cpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_model.cpp"}]},
                            {FileNameJoin[{Global`$flexiblesusyTemplateDir, "lattice_model.hpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_lattice_model.hpp"}]},
                            {FileNameJoin[{Global`$flexiblesusyTemplateDir, "lattice_model.cpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_lattice_model.cpp"}]},
                            {FileNameJoin[{Global`$flexiblesusyTemplateDir, "physical.hpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_physical.hpp"}]},
                            {FileNameJoin[{Global`$flexiblesusyTemplateDir, "physical.cpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_physical.cpp"}]}
                           },
                           diagonalizationPrecision];

           Print["Creating user example spectrum generator program ..."];
           spectrumGeneratorInputFile = "spectrum_generator.hpp.in";
           If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY,
              spectrumGeneratorInputFile = "low_scale_spectrum_generator.hpp.in";];
           WriteUserExample[{{FileNameJoin[{Global`$flexiblesusyTemplateDir, spectrumGeneratorInputFile}],
                              FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_spectrum_generator.hpp"}]},
                             {FileNameJoin[{Global`$flexiblesusyTemplateDir, "run.cpp.in"}],
                              FileNameJoin[{Global`$flexiblesusyOutputDir, "run_" <> FlexibleSUSY`FSModelName <> ".cpp"}]},
                             {FileNameJoin[{Global`$flexiblesusyTemplateDir, "scan.cpp.in"}],
                              FileNameJoin[{Global`$flexiblesusyOutputDir, "scan_" <> FlexibleSUSY`FSModelName <> ".cpp"}]}
                            }];

           PrintHeadline["FlexibleSUSY has finished"];
          ];

End[];

EndPackage[];

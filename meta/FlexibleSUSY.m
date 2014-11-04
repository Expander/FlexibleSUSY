
BeginPackage["FlexibleSUSY`", {"SARAH`", "AnomalousDimension`", "BetaFunction`", "TextFormatting`", "CConversion`", "TreeMasses`", "EWSB`", "Traces`", "SelfEnergies`", "Vertices`", "Phases`", "LoopMasses`", "WriteOut`", "Constraint`", "ThresholdCorrections`", "ConvergenceTester`"}];

FS`Version = StringTrim[Import[FileNameJoin[{Global`$flexiblesusyConfigDir,"version"}], "String"]];
FS`Authors = {"P. Athron", "Jae-hyeon Park", "D. Stöckinger", "A. Voigt"};
FS`Years   = {2013, 2014};
FS`References = Get[FileNameJoin[{Global`$flexiblesusyConfigDir,"references"}]];

Print["*****************************************************************"];
Print["FlexibleSUSY ", FS`Version];
Print["by " <> WriteOut`StringJoinWithSeparator[FS`Authors, ", "] <> ", " <>
      WriteOut`StringJoinWithSeparator[FS`Years, ", "]];
Print[""];
Print["References:"];
Print["  " <> #]& /@ FS`References;
Print["*****************************************************************"];
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
EWSBOutputParameters = {};
SUSYScale;
SUSYScaleFirstGuess;
SUSYScaleInput = {};
SUSYScaleMinimum;
SUSYScaleMaximum;
HighScale;
HighScaleFirstGuess;
HighScaleInput = {};
HighScaleMinimum;
HighScaleMaximum;
LowScale;
LowScaleFirstGuess;
LowScaleInput = {};
LowScaleMinimum;
LowScaleMaximum;
InitialGuessAtLowScale = {};
InitialGuessAtHighScale = {};
OnlyLowEnergyFlexibleSUSY = False;
TreeLevelEWSBSolution = {};
Pole;
FSMinimize;
FSFindRoot;
MZ;
MZDRbar;
MWDRbar;
EDRbar;
UseHiggs2LoopNMSSM;
EffectiveMu;
PotentialLSPParticles = {};
ExtraSLHAOutputBlocks = {};

(* precision of pole mass calculation *)
DefaultPoleMassPrecision = MediumPrecision;
HighPoleMassPrecision    = {SARAH`HiggsBoson, SARAH`PseudoScalar, SARAH`ChargedHiggs};
MediumPoleMassPrecision  = {};
LowPoleMassPrecision     = {};

FSEigenstates = SARAH`EWSB;
FSSolveEWSBTimeConstraint = 120;
FSSimplifyBetaFunctionsTimeConstraint = 120;
FSSolveWeinbergAngleTimeConstraint = 120;

(* list of soft breaking Higgs masses for solving EWSB eqs. *)
FSSoftHiggsMasses = {};

(* EWSB solvers *)
GSLHybrid;   (* hybrid method *)
GSLHybridS;  (* hybrid method with dynamic step size *)
GSLBroyden;  (* Broyden method *)
GSLNewton;   (* Newton method *)
FSEWSBSolvers = { GSLHybrid, GSLHybridS, GSLBroyden };

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
              If[!FlexibleSUSY`OnlyLowEnergyFlexibleSUSY,
                 Print["Warning: FlexibleSUSY`HighScale should be",
                       " set in the model file!"];
                ];
              FlexibleSUSY`HighScale := SARAH`hyperchargeCoupling == SARAH`leftCoupling;
             ];
           If[!ValueQ[FlexibleSUSY`HighScaleFirstGuess],
              If[!FlexibleSUSY`OnlyLowEnergyFlexibleSUSY,
                 Print["Warning: FlexibleSUSY`HighScaleFirstGuess should be",
                       " set in the model file!"];
                ];
              FlexibleSUSY`HighScaleFirstGuess = 2.0 10^16;
             ];
           If[Head[FlexibleSUSY`HighScaleInput] =!= List,
              FlexibleSUSY`HighScaleInput = {};
             ];
           (* LowScale *)
           If[!ValueQ[FlexibleSUSY`LowScale],
              Print["Warning: FlexibleSUSY`LowScale should be",
                    " set in the model file!"];
              FlexibleSUSY`LowScale := SM[MZ];
              ,
              If[FlexibleSUSY`LowScale =!= SM[MZ],
                 Print["Error: The low-scale was set differently from MZ!"];
                 Print["   LowScale = ", FlexibleSUSY`LowScale];
                 Print["   This is currently not supported."];
                 Print["   Please set: LowScale = ", SM[MZ], ";"];
                 Quit[1];
                ];
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

           If[Head[SARAH`MINPAR] =!= List,
              SARAH`MINPAR = {};
             ];
           If[Head[SARAH`EXTPAR] =!= List,
              SARAH`EXTPAR = {};
             ];
           If[Head[FlexibleSUSY`TreeLevelEWSBSolution] =!= List,
              FlexibleSUSY`TreeLevelEWSBSolution = {};
             ];
           If[Head[FlexibleSUSY`ExtraSLHAOutputBlocks] =!= List,
              FlexibleSUSY`ExtraSLHAOutputBlocks = {};
             ];
           If[Head[FlexibleSUSY`EWSBOutputParameters] =!= List,
              Print["Error: EWSBOutputParameters has to be set to a list",
                    " of model parameters chosen to be output of the EWSB eqs."];
              Quit[1];
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
      "@HiggsBoson_" ~~ num_ ~~ "@" /; IntegerQ[ToExpression[num]] :> ToValidCSymbolString[SARAH`HiggsBoson] <> If[TreeMasses`GetDimension[SARAH`HiggsBoson] > 1, "(" <> num <> ")", ""],
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
              templateFile_String, makefileModuleTemplates_List,
              additionalTraces_List:{}, numberOfBaseClassParameters_:0] :=
   Module[{beta, setter, getter, parameterDef, set,
           display, parameterDefaultInit,
           cCtorParameterList, parameterCopyInit, betaParameterList,
           anomDimPrototypes, anomDimFunctions, printParameters, parameters,
           numberOfParameters, clearParameters,
           singleBetaFunctionsDecls, singleBetaFunctionsDefsFiles,
           traceDefs, calcTraces, sarahTraces},
          (* extract list of parameters from the beta functions *)
          parameters = BetaFunction`GetName[#]& /@ betaFun;
          (* count number of parameters *)
          numberOfParameters = BetaFunction`CountNumberOfParameters[betaFun] + numberOfBaseClassParameters;
          (* create C++ functions and parameter declarations *)
          sarahTraces          = Traces`ConvertSARAHTraces[additionalTraces];
          beta                 = BetaFunction`CreateBetaFunction[betaFun, sarahTraces];
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
          singleBetaFunctionsDecls = BetaFunction`CreateSingleBetaFunctionDecl[betaFun];
          traceDefs            = Traces`CreateTraceDefs[betaFun];
          traceDefs            = traceDefs <> Traces`CreateSARAHTraceDefs[sarahTraces];
          calcTraces           = Traces`CreateTraceCalculation[betaFun, "TRACE_STRUCT"];
          calcTraces           = calcTraces <> "\n" <>
                                 Traces`CreateSARAHTraceCalculation[sarahTraces, "TRACE_STRUCT"];
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
                   "@singleBetaFunctionsDecls@" -> IndentText[singleBetaFunctionsDecls],
                   "@traceDefs@"            -> IndentText[IndentText[traceDefs]],
                   "@calcTraces@"           -> IndentText[WrapLines[calcTraces]],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          singleBetaFunctionsDefsFiles = BetaFunction`CreateSingleBetaFunctionDefs[betaFun, templateFile, sarahTraces];
          Print["Creating makefile module for the two-scale method ..."];
          WriteMakefileModule[singleBetaFunctionsDefsFiles,
                              makefileModuleTemplates];
         ];

DefaultValueOf[FlexibleSUSY`Sign[_]] := 1;
DefaultValueOf[_] := 0;

WriteInputParameterClass[inputParameters_List, freePhases_List,
                         lesHouchesInputParameters_List,
                         files_List] :=
   Module[{defineInputParameters, defaultInputParametersInit},
          defaultValues = {#, DefaultValueOf[#]}& /@ inputParameters;
          defineInputParameters = Constraint`DefineInputParameters[Join[inputParameters,freePhases,lesHouchesInputParameters]];
          defaultInputParametersInit = Constraint`InitializeInputParameters[Join[defaultValues,freePhases,lesHouchesInputParameters]];
          WriteOut`ReplaceInFiles[files,
                         { "@defineInputParameters@" -> IndentText[defineInputParameters],
                           "@defaultInputParametersInit@" -> WrapLines[defaultInputParametersInit],
                           Sequence @@ GeneralReplacementRules[]
                         } ];
          ];

WriteConstraintClass[condition_, settings_List, scaleFirstGuess_,
                     {minimumScale_, maximumScale_}, files_List] :=
   Module[{applyConstraint = "", calculateScale, scaleGuess,
           restrictScale,
           setDRbarYukawaCouplings,
           calculateDeltaAlphaEm, calculateDeltaAlphaS,
           calculateGaugeCouplings,
           saveEwsbOutputParameters, restoreEwsbOutputParameters},
          Constraint`SetBetaFunctions[GetBetaFunctions[]];
          applyConstraint = Constraint`ApplyConstraints[settings];
          calculateScale  = Constraint`CalculateScale[condition, "scale"];
          scaleGuess      = Constraint`CalculateScale[scaleFirstGuess, "initial_scale_guess"];
          restrictScale   = Constraint`RestrictScale[{minimumScale, maximumScale}];
          calculateDeltaAlphaEm   = ThresholdCorrections`CalculateDeltaAlphaEm[];
          calculateDeltaAlphaS    = ThresholdCorrections`CalculateDeltaAlphaS[];
          calculateGaugeCouplings = ThresholdCorrections`CalculateGaugeCouplings[];
          setDRbarYukawaCouplings = {
              ThresholdCorrections`SetDRbarYukawaCouplingTop[settings],
              ThresholdCorrections`SetDRbarYukawaCouplingBottom[settings],
              ThresholdCorrections`SetDRbarYukawaCouplingElectron[settings]
          };
          saveEwsbOutputParameters    = Parameters`SaveParameterLocally[FlexibleSUSY`EWSBOutputParameters, "old_", "MODELPARAMETER"];
          restoreEwsbOutputParameters = Parameters`RestoreParameter[FlexibleSUSY`EWSBOutputParameters, "old_", "model"];
          WriteOut`ReplaceInFiles[files,
                 { "@applyConstraint@"      -> IndentText[WrapLines[applyConstraint]],
                   "@calculateScale@"       -> IndentText[WrapLines[calculateScale]],
                   "@scaleGuess@"           -> IndentText[WrapLines[scaleGuess]],
                   "@restrictScale@"        -> IndentText[WrapLines[restrictScale]],
                   "@calculateGaugeCouplings@" -> IndentText[WrapLines[calculateGaugeCouplings]],
                   "@calculateDeltaAlphaEm@" -> IndentText[WrapLines[calculateDeltaAlphaEm]],
                   "@calculateDeltaAlphaS@"  -> IndentText[WrapLines[calculateDeltaAlphaS]],
                   "@setDRbarUpQuarkYukawaCouplings@"   -> IndentText[WrapLines[setDRbarYukawaCouplings[[1]]]],
                   "@setDRbarDownQuarkYukawaCouplings@" -> IndentText[WrapLines[setDRbarYukawaCouplings[[2]]]],
                   "@setDRbarElectronYukawaCouplings@"  -> IndentText[WrapLines[setDRbarYukawaCouplings[[3]]]],
                   "@saveEwsbOutputParameters@"    -> IndentText[saveEwsbOutputParameters],
                   "@restoreEwsbOutputParameters@" -> IndentText[restoreEwsbOutputParameters],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteInitialGuesserClass[lowScaleGuess_List, highScaleGuess_List, files_List] :=
   Module[{initialGuessAtLowScale, initialGuessAtHighScale, setDRbarYukawaCouplings,
           allSettings},
          initialGuessAtLowScale  = Constraint`ApplyConstraints[lowScaleGuess];
          initialGuessAtHighScale = Constraint`ApplyConstraints[highScaleGuess];
          allSettings             = Join[lowScaleGuess, highScaleGuess];
          setDRbarYukawaCouplings = {
              ThresholdCorrections`SetDRbarYukawaCouplingTop[allSettings],
              ThresholdCorrections`SetDRbarYukawaCouplingBottom[allSettings],
              ThresholdCorrections`SetDRbarYukawaCouplingElectron[allSettings]
          };
          WriteOut`ReplaceInFiles[files,
                 { "@initialGuessAtLowScale@"  -> IndentText[WrapLines[initialGuessAtLowScale]],
                   "@initialGuessAtHighScale@" -> IndentText[WrapLines[initialGuessAtHighScale]],
                   "@setDRbarUpQuarkYukawaCouplings@"   -> IndentText[WrapLines[setDRbarYukawaCouplings[[1]]]],
                   "@setDRbarDownQuarkYukawaCouplings@" -> IndentText[WrapLines[setDRbarYukawaCouplings[[2]]]],
                   "@setDRbarElectronYukawaCouplings@"  -> IndentText[WrapLines[setDRbarYukawaCouplings[[3]]]],
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

FindVEV[gauge_] :=
    Module[{result, vev},
           vev = Cases[SARAH`DEFINITION[FlexibleSUSY`FSEigenstates][SARAH`VEVs],
                       {_,{v_,_},{gauge,_},{p_,_},___} | {_,{v_,_},{s_,_},{gauge,_},___} :> v];
           If[vev === {},
              Print["Error: could not find VEV for gauge eigenstate ", gauge];
              Quit[1];
             ];
           vev[[1]]
          ];

GetDimOfVEV[vev_] :=
    Switch[SARAH`getDimParameters[vev],
           {}                         , 1,
           {0}                        , 1,
           {1}                        , 1,
           {idx_}                     , SARAH`getDimParameters[vev][[1]]
          ];

ExpandGaugeIndices[gauge_, number_] :=
    Table[gauge[i], {i,1,number}];

ExpandGaugeIndices[gauge_List] :=
    Flatten[ExpandGaugeIndices[#,GetDimOfVEV[FindVEV[#]]]& /@ gauge];

(* Returns a list of two-component list where the information is stored
   which Higgs corresponds to which EWSB eq.

   Example: MRSSM
   It[] := CreateHiggsToEWSBEqAssociation[]
   Out[] = {{hh, 1}, {hh, 2}, {hh, 4}, {hh, 3}}
 *)
CreateHiggsToEWSBEqAssociation[] :=
    Module[{association = {}, v, phi, sigma, higgs, numberOfVEVs, numberOfHiggses, vevs,
            vev, dimVev},
           vevs = Cases[SARAH`DEFINITION[FlexibleSUSY`FSEigenstates][SARAH`VEVs],
                        {_,{v_,_},{s_,_},{p_,_},___} :> {v,s,p}];
           (* list of phi fields, ordered according to VEVs / EWSB eqs. *)
           phi = Transpose[vevs][[3]];
           phi = ExpandGaugeIndices[phi];
           (* list of phi fields, ordered according to Higgs mixing *)
           higgsGaugeES = Cases[SARAH`DEFINITION[FlexibleSUSY`FSEigenstates][SARAH`MatterSector],
                                {gauge_List, {SARAH`HiggsBoson, _}} :> gauge][[1]];
           higgsGaugeES = ExpandGaugeIndices[higgsGaugeES];
           (* find positions of phi in higgsGaugeES *)
           {SARAH`HiggsBoson,#}& /@ (Flatten[Position[higgsGaugeES, #]& /@ phi])
          ];

WriteModelClass[massMatrices_List, ewsbEquations_List,
                parametersFixedByEWSB_List, ewsbSolution_List, freePhases_List,
                nPointFunctions_List, vertexRules_List, phases_List,
                files_List, diagonalizationPrecision_List] :=
    Module[{massGetters = "", k,
            mixingMatrixGetters = "",
            tadpoleEqPrototypes = "", tadpoleEqFunctions = "",
            numberOfEWSBEquations = Length[ewsbEquations], calculateTreeLevelTadpoles = "",
            ewsbInitialGuess = "", physicalMassesDef = "", mixingMatricesDef = "",
            physicalMassesInit = "", physicalMassesInitNoLeadingComma = "", mixingMatricesInit = "",
            massCalculationPrototypes = "", massCalculationFunctions = "",
            calculateAllMasses = "",
            calculateOneLoopTadpoles = "", calculateTwoLoopTadpoles = "",
            selfEnergyPrototypes = "", selfEnergyFunctions = "",
            twoLoopTadpolePrototypes = "", twoLoopTadpoleFunctions = "",
            twoLoopSelfEnergyPrototypes = "", twoLoopSelfEnergyFunctions = "",
            thirdGenerationHelperPrototypes = "", thirdGenerationHelperFunctions = "",
            phasesDefinition = "", phasesGetterSetters = "",
            phasesInit = "",
            loopMassesPrototypes = "", loopMassesFunctions = "",
            runningDRbarMassesPrototypes = "", runningDRbarMassesFunctions = "",
            callAllLoopMassFunctions = "",
            callAllLoopMassFunctionsInThreads = "",
            printMasses = "", printMixingMatrices = "",
            masses, mixingMatrices, oneLoopTadpoles,
            dependenceNumPrototypes, dependenceNumFunctions,
            clearOutputParameters = "", solveEwsbTreeLevel = "",
            clearPhases = "",
            saveEwsbOutputParameters, restoreEwsbOutputParameters,
            softScalarMasses, softHiggsMasses,
            saveSoftHiggsMasses, restoreSoftHiggsMasses,
            solveTreeLevelEWSBviaSoftHiggsMasses,
            copyDRbarMassesToPoleMasses = "",
            reorderDRbarMasses = "", reorderPoleMasses = "",
            higgsToEWSBEqAssociation,
            twoLoopHiggsHeaders = "",
            lspGetters = "", lspFunctions = "",
            gslEWSBRootFinders = "",
            enablePoleMassThreads = True
           },
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
           clearPhases = Parameters`ClearPhases[phases];
           calculateAllMasses = TreeMasses`CallMassCalculationFunctions[massMatrices];
           tadpoleEqPrototypes = EWSB`CreateEWSBEqPrototype[SARAH`HiggsBoson];
           tadpoleEqFunctions  = EWSB`CreateEWSBEqFunction[SARAH`HiggsBoson, ewsbEquations];
           If[Length[parametersFixedByEWSB] != numberOfEWSBEquations,
              Print["Error: There are ", numberOfEWSBEquations, " EWSB ",
                    "equations, but you want to fix ", Length[parametersFixedByEWSB],
                    " parameters: ", parametersFixedByEWSB];
             ];
           higgsToEWSBEqAssociation     = CreateHiggsToEWSBEqAssociation[];
           oneLoopTadpoles              = Cases[nPointFunctions, SelfEnergies`Tadpole[___]];
           calculateOneLoopTadpoles     = SelfEnergies`FillArrayWithOneLoopTadpoles[higgsToEWSBEqAssociation];
           If[SARAH`UseHiggs2LoopMSSM === True ||
              FlexibleSUSY`UseHiggs2LoopNMSSM === True,
              calculateTwoLoopTadpoles  = SelfEnergies`FillArrayWithTwoLoopTadpoles[SARAH`HiggsBoson];
              {thirdGenerationHelperPrototypes, thirdGenerationHelperFunctions} = TreeMasses`CreateThirdGenerationHelpers[];
             ];
           If[SARAH`UseHiggs2LoopMSSM === True,
              {twoLoopTadpolePrototypes, twoLoopTadpoleFunctions} = SelfEnergies`CreateTwoLoopTadpolesMSSM[SARAH`HiggsBoson];
              {twoLoopSelfEnergyPrototypes, twoLoopSelfEnergyFunctions} = SelfEnergies`CreateTwoLoopSelfEnergiesMSSM[{SARAH`HiggsBoson, SARAH`PseudoScalar}];
              twoLoopHiggsHeaders = "#include \"sfermions.hpp\"\n#include \"mssm_twoloophiggs.h\"\n";
             ];
           If[FlexibleSUSY`UseHiggs2LoopNMSSM === True,
              {twoLoopTadpolePrototypes, twoLoopTadpoleFunctions} = SelfEnergies`CreateTwoLoopTadpolesNMSSM[SARAH`HiggsBoson];
              {twoLoopSelfEnergyPrototypes, twoLoopSelfEnergyFunctions} = SelfEnergies`CreateTwoLoopSelfEnergiesNMSSM[{SARAH`HiggsBoson, SARAH`PseudoScalar}];
              twoLoopHiggsHeaders = "#include \"sfermions.hpp\"\n#include \"nmssm_twoloophiggs.h\"\n";
             ];
           calculateTreeLevelTadpoles   = EWSB`FillArrayWithEWSBEqs[SARAH`HiggsBoson, parametersFixedByEWSB, freePhases];
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
           enablePoleMassThreads = False;
           callAllLoopMassFunctions     = LoopMasses`CallAllPoleMassFunctions[FlexibleSUSY`FSEigenstates, enablePoleMassThreads];
           enablePoleMassThreads = True;
           callAllLoopMassFunctionsInThreads = LoopMasses`CallAllPoleMassFunctions[FlexibleSUSY`FSEigenstates, enablePoleMassThreads];
           masses                       = FlexibleSUSY`M[TreeMasses`GetMassEigenstate[#]]& /@ massMatrices;
           {lspGetters, lspFunctions}   = LoopMasses`CreateLSPFunctions[FlexibleSUSY`PotentialLSPParticles];
           printMasses                  = WriteOut`PrintParameters[masses, "ostr"];
           mixingMatrices               = Flatten[TreeMasses`GetMixingMatrixSymbol[#]& /@ massMatrices];
           printMixingMatrices          = WriteOut`PrintParameters[mixingMatrices, "ostr"];
           dependenceNumPrototypes      = TreeMasses`CreateDependenceNumPrototypes[massMatrices];
           dependenceNumFunctions       = TreeMasses`CreateDependenceNumFunctions[massMatrices];
           saveEwsbOutputParameters     = Parameters`SaveParameterLocally[FlexibleSUSY`EWSBOutputParameters, "one_loop_", ""];
           restoreEwsbOutputParameters  = Parameters`RestoreParameter[FlexibleSUSY`EWSBOutputParameters, "one_loop_", ""];
           If[Head[SARAH`ListSoftBreakingScalarMasses] === List,
              softScalarMasses          = DeleteDuplicates[SARAH`ListSoftBreakingScalarMasses];,
              softScalarMasses          = {};
             ];
           (* find soft Higgs masses that appear in tree-level EWSB eqs. *)
           If[Head[FlexibleSUSY`FSSoftHiggsMasses] =!= List ||
              FlexibleSUSY`FSSoftHiggsMasses === {},
              softHiggsMasses = Select[softScalarMasses, (!FreeQ[ewsbEquations, #])&];
              ,
              softHiggsMasses = FlexibleSUSY`FSSoftHiggsMasses;
             ];
           softHiggsMasses              = Parameters`DecreaseIndexLiterals[Parameters`ExpandExpressions[Parameters`AppendGenerationIndices[softHiggsMasses]]];
           saveSoftHiggsMasses          = Parameters`SaveParameterLocally[softHiggsMasses, "old_", ""];
           restoreSoftHiggsMasses       = Parameters`RestoreParameter[softHiggsMasses, "old_", ""];
           If[Head[softHiggsMasses] === List && Length[softHiggsMasses] > 0,
              solveTreeLevelEWSBviaSoftHiggsMasses = EWSB`SolveTreeLevelEwsbVia[ewsbEquations, softHiggsMasses];,
              solveTreeLevelEWSBviaSoftHiggsMasses = "";
             ];
           gslEWSBRootFinders           = EWSB`CreateEWSBRootFinders[FlexibleSUSY`FSEWSBSolvers];
           reorderDRbarMasses           = TreeMasses`ReorderGoldstoneBosons[""];
           reorderPoleMasses            = TreeMasses`ReorderGoldstoneBosons["PHYSICAL"];
           WriteOut`ReplaceInFiles[files,
                          { "@lspGetters@"           -> IndentText[lspGetters],
                            "@lspFunctions@"         -> lspFunctions,
                            "@massGetters@"          -> IndentText[massGetters],
                            "@mixingMatrixGetters@"  -> IndentText[mixingMatrixGetters],
                            "@tadpoleEqPrototypes@"  -> IndentText[tadpoleEqPrototypes],
                            "@tadpoleEqFunctions@"   -> tadpoleEqFunctions,
                            "@numberOfEWSBEquations@"-> ToString[numberOfEWSBEquations],
                            "@calculateTreeLevelTadpoles@" -> IndentText[calculateTreeLevelTadpoles],
                            "@calculateOneLoopTadpoles@"   -> IndentText[calculateOneLoopTadpoles],
                            "@calculateTwoLoopTadpoles@"   -> IndentText[calculateTwoLoopTadpoles],
                            "@clearOutputParameters@"  -> IndentText[clearOutputParameters],
                            "@clearPhases@"            -> IndentText[clearPhases],
                            "@copyDRbarMassesToPoleMasses@" -> IndentText[copyDRbarMassesToPoleMasses],
                            "@reorderDRbarMasses@"     -> IndentText[reorderDRbarMasses],
                            "@reorderPoleMasses@"      -> IndentText[reorderPoleMasses],
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
                            "@twoLoopTadpolePrototypes@"  -> IndentText[twoLoopTadpolePrototypes],
                            "@twoLoopTadpoleFunctions@"   -> twoLoopTadpoleFunctions,
                            "@twoLoopSelfEnergyPrototypes@" -> IndentText[twoLoopSelfEnergyPrototypes],
                            "@twoLoopSelfEnergyFunctions@"  -> twoLoopSelfEnergyFunctions,
                            "@twoLoopHiggsHeaders@"       -> twoLoopHiggsHeaders,
                            "@thirdGenerationHelperPrototypes@" -> IndentText[thirdGenerationHelperPrototypes],
                            "@thirdGenerationHelperFunctions@"  -> thirdGenerationHelperFunctions,
                            "@phasesDefinition@"          -> IndentText[phasesDefinition],
                            "@phasesGetterSetters@"          -> IndentText[phasesGetterSetters],
                            "@phasesInit@"                   -> IndentText[WrapLines[phasesInit]],
                            "@loopMassesPrototypes@"         -> IndentText[WrapLines[loopMassesPrototypes]],
                            "@loopMassesFunctions@"          -> WrapLines[loopMassesFunctions],
                            "@runningDRbarMassesPrototypes@" -> IndentText[runningDRbarMassesPrototypes],
                            "@runningDRbarMassesFunctions@"  -> WrapLines[runningDRbarMassesFunctions],
                            "@callAllLoopMassFunctions@"     -> IndentText[callAllLoopMassFunctions],
                            "@callAllLoopMassFunctionsInThreads@" -> IndentText[callAllLoopMassFunctionsInThreads],
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
                            "@gslEWSBRootFinders@"           -> IndentText[IndentText[gslEWSBRootFinders]],
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteUserExample[inputParameters_List, freePhases_List,
                 lesHouchesInputParameters_List, files_List] :=
    Module[{parseCmdLineOptions, printCommandLineOptions},
           allParameters = Join[inputParameters,freePhases,lesHouchesInputParameters];
           parseCmdLineOptions = WriteOut`ParseCmdLineOptions[allParameters];
           printCommandLineOptions = WriteOut`PrintCmdLineOptions[allParameters];
           WriteOut`ReplaceInFiles[files,
                          { "@parseCmdLineOptions@" -> IndentText[IndentText[parseCmdLineOptions]],
                            "@printCommandLineOptions@" -> IndentText[IndentText[printCommandLineOptions]],
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WritePlotScripts[files_List] :=
    Module[{},
           WriteOut`ReplaceInFiles[files,
                          { Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteMakefileModule[rgeFile_List, files_List] :=
    Module[{concatenatedFileList},
           concatenatedFileList = "\t" <> WriteOut`StringJoinWithSeparator[rgeFile, " \\\n\t"];
           WriteOut`ReplaceInFiles[files,
                          { "@generatedBetaFunctionModules@" -> concatenatedFileList,
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteUtilitiesClass[massMatrices_List, betaFun_List, minpar_List, extpar_List,
                    lesHouchesInputParameters_List, extraSLHAOutputBlocks_List,
                    files_List] :=
    Module[{k, particles, susyParticles, smParticles,
            fillSpectrumVectorWithSusyParticles = "",
            fillSpectrumVectorWithSMParticles = "",
            particleLaTeXNames = "",
            particleNames = "", particleEnum = "", particleMultiplicity = "",
            parameterNames = "", parameterEnum = "", numberOfParameters = 0,
            isLowEnergyModel = "false",
            fillInputParametersFromMINPAR = "", fillInputParametersFromEXTPAR = "",
            writeSLHAMassBlock = "", writeSLHAMixingMatricesBlocks = "",
            writeSLHAModelParametersBlocks = "", writeSLHAMinparBlock = "",
            writeSLHAExtparBlock = "", readLesHouchesInputParameters,
            writeExtraSLHAOutputBlock = "",
            readLesHouchesOutputParameters, readLesHouchesPhyicalParameters,
            convertMixingsToSLHAConvention = "",
            numberOfDRbarBlocks, drBarBlockNames},
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
           isLowEnergyModel = If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY === True, "true", "false"];
           fillInputParametersFromMINPAR = Parameters`FillInputParametersFromTuples[minpar];
           fillInputParametersFromEXTPAR = Parameters`FillInputParametersFromTuples[extpar];
           readLesHouchesInputParameters = WriteOut`ReadLesHouchesInputParameters[lesHouchesInputParameters];
           readLesHouchesOutputParameters = WriteOut`ReadLesHouchesOutputParameters[];
           readLesHouchesPhyicalParameters = WriteOut`ReadLesHouchesPhysicalParameters[];
           writeSLHAMassBlock = WriteOut`WriteSLHAMassBlock[massMatrices];
           writeSLHAMixingMatricesBlocks  = WriteOut`WriteSLHAMixingMatricesBlocks[];
           writeSLHAModelParametersBlocks = WriteOut`WriteSLHAModelParametersBlocks[];
           writeSLHAMinparBlock = WriteOut`WriteSLHAMinparBlock[minpar];
           writeSLHAExtparBlock = WriteOut`WriteSLHAExtparBlock[extpar];
           writeExtraSLHAOutputBlock = WriteOut`WriteExtraSLHAOutputBlock[extraSLHAOutputBlocks];
           numberOfDRbarBlocks  = WriteOut`GetNumberOfDRbarBlocks[];
           drBarBlockNames      = WriteOut`GetDRbarBlockNames[];
           convertMixingsToSLHAConvention = WriteOut`ConvertMixingsToSLHAConvention[massMatrices];
           WriteOut`ReplaceInFiles[files,
                          { "@fillSpectrumVectorWithSusyParticles@" -> IndentText[fillSpectrumVectorWithSusyParticles],
                            "@fillSpectrumVectorWithSMParticles@"   -> IndentText[IndentText[fillSpectrumVectorWithSMParticles]],
                            "@particleEnum@"       -> IndentText[WrapLines[particleEnum]],
                            "@particleMultiplicity@" -> IndentText[WrapLines[particleMultiplicity]],
                            "@particleNames@"      -> IndentText[WrapLines[particleNames]],
                            "@particleLaTeXNames@" -> IndentText[WrapLines[particleLaTeXNames]],
                            "@parameterEnum@"     -> IndentText[WrapLines[parameterEnum]],
                            "@parameterNames@"     -> IndentText[WrapLines[parameterNames]],
                            "@isLowEnergyModel@"   -> isLowEnergyModel,
                            "@fillInputParametersFromMINPAR@" -> IndentText[fillInputParametersFromMINPAR],
                            "@fillInputParametersFromEXTPAR@" -> IndentText[fillInputParametersFromEXTPAR],
                            "@readLesHouchesInputParameters@" -> IndentText[readLesHouchesInputParameters],
                            "@readLesHouchesOutputParameters@" -> IndentText[readLesHouchesOutputParameters],
                            "@readLesHouchesPhyicalParameters@" -> IndentText[readLesHouchesPhyicalParameters],
                            "@writeSLHAMassBlock@" -> IndentText[writeSLHAMassBlock],
                            "@writeSLHAMixingMatricesBlocks@"  -> IndentText[writeSLHAMixingMatricesBlocks],
                            "@writeSLHAModelParametersBlocks@" -> IndentText[writeSLHAModelParametersBlocks],
                            "@writeSLHAMinparBlock@"           -> IndentText[writeSLHAMinparBlock],
                            "@writeSLHAExtparBlock@"           -> IndentText[writeSLHAExtparBlock],
                            "@writeExtraSLHAOutputBlock@"      -> IndentText[writeExtraSLHAOutputBlock],
                            "@convertMixingsToSLHAConvention@" -> IndentText[convertMixingsToSLHAConvention],
                            "@numberOfDRbarBlocks@"            -> ToString[numberOfDRbarBlocks],
                            "@drBarBlockNames@"                -> WrapLines[drBarBlockNames],
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
           If[SARAH`SupersymmetricModel === False,
              AppendTo[fileNames, "BetaLijkl.m"];
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
    NeedToUpdateTarget[
	"self-energy",
	GetSelfEnergyFileNames[$sarahCurrentOutputMainDir, eigenstates]];

GetTadpoleFileName[outputDir_String, eigenstates_] :=
    FileNameJoin[{outputDir, ToString[eigenstates],
                  "One-Loop", "Tadpoles1Loop.m"}];

TadpoleFileExists[outputDir_String, eigenstates_] :=
    FileExists[GetTadpoleFileName[outputDir, eigenstates]];

TadpoleFilesModificationTimeInSeconds[outputDir_String, eigenstates_] :=
    LatestModificationTimeInSeconds[GetTadpoleFileName[outputDir, eigenstates]];

NeedToCalculateTadpoles[eigenstates_] :=
    NeedToUpdateTarget[
	"tadpole",
	GetTadpoleFileName[$sarahCurrentOutputMainDir, eigenstates]];

GetUnrotatedParticlesFileName[outputDir_String, eigenstates_] :=
    FileNameJoin[{outputDir, ToString[eigenstates],
                  "One-Loop", "UnrotatedParticles.m"}];

UnrotatedParticlesFilesExist[outputDir_String, eigenstates_] :=
    FileExists[GetUnrotatedParticlesFileName[outputDir, eigenstates]];

UnrotatedParticlesFilesModificationTimeInSeconds[outputDir_String, eigenstates_] :=
    LatestModificationTimeInSeconds[GetUnrotatedParticlesFileName[outputDir, eigenstates]];

NeedToCalculateUnrotatedParticles[eigenstates_] :=
    NeedToUpdateTarget[
	"unrotated particle",
	GetUnrotatedParticlesFileName[$sarahCurrentOutputMainDir,eigenstates]];

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
    NeedToUpdateTarget["RGE", GetRGEFileNames[$sarahCurrentOutputMainDir]];

GetVertexRuleFileName[outputDir_String, eigenstates_] :=
    FileNameJoin[{outputDir, ToString[eigenstates], "Vertices",
		  "FSVertexRules.m"}];

NeedToCalculateVertices[eigenstates_] :=
    NeedToUpdateTarget[
	"vertex",
	GetVertexRuleFileName[$sarahCurrentOutputMainDir, eigenstates]];

NeedToUpdateTarget[name_String, targets_List] := Module[{
	targetsExist = FilesExist[targets],
	targetTimeStamp = LatestModificationTimeInSeconds[targets],
	sarahModelFileTimeStamp = SARAHModelFileModificationTimeInSeconds[],
	files = If[Length[targets] === 1, "file", "files"],
	them = If[Length[targets] === 1, "it", "them"]
    },
    If[targetsExist,
       If[sarahModelFileTimeStamp > targetTimeStamp,
	  Print["SARAH model files are newer than ", name,
		" ", files, ", updating ", them, " ..."];
	  True,
	  Print["Found up-to-date ", name, " ", files, "."];
	  False
       ],
       Print[name, " ", files, " not found, producing ", them, " ..."];
       True
    ]
];

NeedToUpdateTarget[name_String, target_] :=
    NeedToUpdateTarget[name, {target}];

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
                     SARAH`BetaDGi, SARAH`BetaLijkl };
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

ReadPoleMassPrecisions[defaultPrecision_Symbol, highPrecisionList_List,
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
                                            SARAH`strongCoupling }]]];
           Complement[allParameters, fixedParameters]
          ];

Options[MakeFlexibleSUSY] :=
    {
        InputFile -> "FlexibleSUSY.m"
    };

MakeFlexibleSUSY[OptionsPattern[]] :=
    Module[{nPointFunctions, runInputFile, initialGuesserInputFile,
            susyBetaFunctions, susyBreakingBetaFunctions,
            numberOfSusyParameters, anomDim,
            ewsbEquations, massMatrices, phases,
            diagonalizationPrecision, allParticles, freePhases, ewsbSolution,
            fixedParameters, treeLevelEwsbOutputFile,
            lesHouchesInputParameters,
            extraSLHAOutputBlocks,
	    vertexRules, vertexRuleFileName,
	    Lat$massMatrices},
           (* check if SARAH`Start[] was called *)
           If[!ValueQ[Model`Name],
              Print["Error: Model`Name is not defined.  Did you call SARAH`Start[\"Model\"]?"];
              Quit[1];
             ];
           CheckSARAHVersion[];
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
           nPointFunctions = EnforceCpColorStructures @ StripInvalidFieldIndices @
	      Join[PrepareSelfEnergies[FSEigenstates], PrepareTadpoles[FSEigenstates]];
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
           susyBetaFunctions = { SARAH`BetaLijkl,
                                 SARAH`BetaWijkl,
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
                           FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_susy_parameters.cpp"}]}},
                         "two_scale_susy_beta_.cpp.in",
                         {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale.mk.in"}],
                           FileNameJoin[{Global`$flexiblesusyOutputDir, "two_scale_susy.mk"}]}}
                        ];

           susyBreakingBetaFunctions = ConvertSarahRGEs[susyBreakingBetaFunctions];
           susyBreakingBetaFunctions = Select[susyBreakingBetaFunctions, (BetaFunction`GetAllBetaFunctions[#]!={})&];
           Parameters`AddRealParameter[(BetaFunction`GetName /@ susyBreakingBetaFunctions) /.
                                       a_[Susyno`LieGroups`i1] :> a /.
                                       a_[Susyno`LieGroups`i1,SARAH`i2] :> a];

           allBetaFunctions = Join[susyBetaFunctions, susyBreakingBetaFunctions];

           (* store all model parameters *)
           allParameters = Join[BetaFunction`GetName /@ susyBetaFunctions,
                                BetaFunction`GetName /@ susyBreakingBetaFunctions] /.
                               a_[Susyno`LieGroups`i1] :> a /.
                               a_[Susyno`LieGroups`i1,SARAH`i2] :> a;
           allIndexReplacementRules = Parameters`CreateIndexReplacementRules[allParameters];
           Parameters`SetModelParameters[allParameters];
           FlexibleSUSY`FSLesHouchesList = SA`LHList;

           (* search for unfixed parameters *)
           Constraint`CheckConstraint[FlexibleSUSY`LowScaleInput, "LowScaleInput"];
           Constraint`CheckConstraint[FlexibleSUSY`SUSYScaleInput, "SUSYScaleInput"];
           Constraint`CheckConstraint[FlexibleSUSY`HighScaleInput, "HighScaleInput"];
           Constraint`CheckConstraint[FlexibleSUSY`InitialGuessAtLowScale, "InitialGuessAtLowScale"];
           Constraint`CheckConstraint[FlexibleSUSY`InitialGuessAtHighScale, "InitialGuessAtHighScale"];
           Constraint`SanityCheck[Join[FlexibleSUSY`InitialGuessAtLowScale,
                                       FlexibleSUSY`InitialGuessAtHighScale],
                                  "initial guess"
                                 ];
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

           lesHouchesInputParameters = DeleteDuplicates[Flatten[Cases[Join[FlexibleSUSY`LowScaleInput,
                                                                           FlexibleSUSY`SUSYScaleInput,
                                                                           FlexibleSUSY`HighScaleInput],
                                                                      SARAH`LHInput[p_] :> p,
                                                                      Infinity
                                                                     ]
                                                               ]
                                                       ];

           lesHouchesInputParameters = Select[Join[{BetaFunction`GetName[#], Symbol[ToValidCSymbolString[BetaFunction`GetName[#]] <> "Input"], #[[2]]}& /@ susyBetaFunctions,
                                                   {BetaFunction`GetName[#], Symbol[ToValidCSymbolString[BetaFunction`GetName[#]] <> "Input"], #[[2]]}& /@ susyBreakingBetaFunctions] /.
                                              a_[Susyno`LieGroups`i1] :> a /.
                                              a_[Susyno`LieGroups`i1,SARAH`i2] :> a,
                                              MemberQ[lesHouchesInputParameters,#[[1]]]&];

           Parameters`SetInputParameters[Join[Parameters`GetInputParameters[],
                                              (#[[2]])& /@ lesHouchesInputParameters]];

           (* replace LHInput[p] by pInput in the constraints *)
           FlexibleSUSY`LowScaleInput = FlexibleSUSY`LowScaleInput /.
               (Rule[SARAH`LHInput[#[[1]]], #[[2]]]& /@ lesHouchesInputParameters);
           FlexibleSUSY`SUSYScaleInput = FlexibleSUSY`SUSYScaleInput /.
               (Rule[SARAH`LHInput[#[[1]]], #[[2]]]& /@ lesHouchesInputParameters);
           FlexibleSUSY`HighScaleInput = FlexibleSUSY`HighScaleInput /.
               (Rule[SARAH`LHInput[#[[1]]], #[[2]]]& /@ lesHouchesInputParameters);

           If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY === True,
              lesHouchesInputParameters = Join[FlexibleSUSY`FSUnfixedParameters,
                                               lesHouchesInputParameters];
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
                         "two_scale_soft_beta_.cpp.in",
                         {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale.mk.in"}],
                           FileNameJoin[{Global`$flexiblesusyOutputDir, "two_scale_soft.mk"}]}},
                         If[Head[SARAH`TraceAbbr] === List, SARAH`TraceAbbr, {}],
                         numberOfSusyParameters];

           ewsbEquations = SARAH`TadpoleEquations[FSEigenstates] /.
                           Parameters`ApplyGUTNormalization[] /.
                           allIndexReplacementRules /.
                           SARAH`sum[idx_, start_, stop_, expr_] :> Sum[expr, {idx,start,stop}];
           If[Head[ewsbEquations] =!= List,
              Print["Error: Could not find EWSB equations for eigenstates ",
                    FSEigenstates];
              Quit[1];
             ];

           ewsbEquations = Parameters`ExpandExpressions[ewsbEquations];
           FlexibleSUSY`EWSBOutputParameters = Parameters`DecreaseIndexLiterals[FlexibleSUSY`EWSBOutputParameters];

           If[FlexibleSUSY`TreeLevelEWSBSolution === {},
              (* trying to find an analytic solution for the EWSB eqs. *)
              treeLevelEwsbOutputFile = FileNameJoin[{Global`$flexiblesusyOutputDir,
                                                      FlexibleSUSY`FSModelName <> "_tree_level_EWSB_solution.m"}];
              Print["Solving EWSB equations for ", FlexibleSUSY`EWSBOutputParameters," ..."];
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
                                    {#[[2]], #[[3]]}& /@ lesHouchesInputParameters,
                                    {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "input_parameters.hpp.in"}],
                                      FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_input_parameters.hpp"}]}}
                                   ];

	   On[Assert];

           Lat$massMatrices = ConvertSarahMassMatrices[] /.
                          Parameters`ApplyGUTNormalization[] //.
                          { SARAH`sum[j_, start_, end_, expr_] :> (Sum[expr, {j,start,end}]) };
           massMatrices = Lat$massMatrices /. allIndexReplacementRules;
	   Lat$massMatrices = LatticeUtils`FixDiagonalization[Lat$massMatrices];

           allParticles = FlexibleSUSY`M[GetMassEigenstate[#]]& /@ massMatrices;
           allOutputParameters = DeleteCases[DeleteDuplicates[
               Join[allParticles,
                    Flatten[GetMixingMatrixSymbol[#]& /@ massMatrices]]], Null];

           Parameters`SetOutputParameters[allOutputParameters];

           extraSLHAOutputBlocks = Parameters`DecreaseIndexLiterals[
               FlexibleSUSY`ExtraSLHAOutputBlocks,
               Join[Parameters`GetOutputParameters[], Parameters`GetModelParameters[]]
           ];

           PrintHeadline["Creating utilities"];
           Print["Creating class for convergence tester ..."];
           WriteConvergenceTesterClass[allParticles,
               {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "convergence_tester.hpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_convergence_tester.hpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_convergence_tester.hpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_convergence_tester.hpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_convergence_tester.cpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_convergence_tester.cpp"}]}
               }
                                      ];

           Print["Creating utilities class ..."];
           WriteUtilitiesClass[massMatrices, Join[susyBetaFunctions, susyBreakingBetaFunctions],
                               MINPAR, EXTPAR,
                               lesHouchesInputParameters, extraSLHAOutputBlocks,
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
                             {FileNameJoin[{Global`$flexiblesusyTemplateDir, "plot_rgflow.gnuplot.in"}],
                              FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_plot_rgflow.gnuplot"}]}}
                           ];

           PrintHeadline["Creating constraints"];
           Print["Creating class for high-scale constraint ..."];
           WriteConstraintClass[FlexibleSUSY`HighScale,
                                FlexibleSUSY`HighScaleInput,
                                FlexibleSUSY`HighScaleFirstGuess,
                                {FlexibleSUSY`HighScaleMinimum, FlexibleSUSY`HighScaleMaximum},
                                {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "high_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_high_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_high_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_high_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_high_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_high_scale_constraint.cpp"}]}
                                }
                               ];

           Print["Creating class for susy-scale constraint ..."];
           WriteConstraintClass[FlexibleSUSY`SUSYScale,
                                FlexibleSUSY`SUSYScaleInput,
                                FlexibleSUSY`SUSYScaleFirstGuess,
                                {FlexibleSUSY`SUSYScaleMinimum, FlexibleSUSY`SUSYScaleMaximum},
                                {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "susy_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_susy_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_susy_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_susy_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_susy_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_susy_scale_constraint.cpp"}]}
                                }
                               ];

           Print["Creating class for low-scale constraint ..."];
           WriteConstraintClass[FlexibleSUSY`LowScale,
                                FlexibleSUSY`LowScaleInput,
                                FlexibleSUSY`LowScaleFirstGuess,
                                {FlexibleSUSY`LowScaleMinimum, FlexibleSUSY`LowScaleMaximum},
                                {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "low_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_low_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_low_scale_constraint.hpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_low_scale_constraint.hpp"}]},
                                 {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_low_scale_constraint.cpp.in"}],
                                  FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_low_scale_constraint.cpp"}]}
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
                                      FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_initial_guesser.cpp"}]}
                                    }
                                   ];

           phases = ConvertSarahPhases[SARAH`ParticlePhases];

           (* determin diagonalization precision for each particle *)
           diagonalizationPrecision = ReadPoleMassPrecisions[
               DefaultPoleMassPrecision,
               Flatten[{HighPoleMassPrecision}],
               Flatten[{MediumPoleMassPrecision}],
               Flatten[{LowPoleMassPrecision}],
               FSEigenstates];

	   vertexRuleFileName =
	      GetVertexRuleFileName[$sarahCurrentOutputMainDir, FSEigenstates];
	   If[NeedToCalculateVertices[FSEigenstates],
	      Put[vertexRules =
		      Vertices`VertexRules[nPointFunctions, Lat$massMatrices],
		  vertexRuleFileName],
	      vertexRules = Get[vertexRuleFileName]];

           PrintHeadline["Creating model"];
           Print["Creating class for model ..."];
           WriteModelClass[massMatrices, ewsbEquations,
                           FlexibleSUSY`EWSBOutputParameters, ewsbSolution, freePhases,
                           nPointFunctions, vertexRules, phases,
                           {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "model.hpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_model.hpp"}]},
                            {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_model.hpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_model.hpp"}]},
                            {FileNameJoin[{Global`$flexiblesusyTemplateDir, "two_scale_model.cpp.in"}],
                             FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_model.cpp"}]},
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
           WriteUserExample[FlexibleSUSY`InputParameters,
                            Complement[freePhases, FlexibleSUSY`InputParameters],
                            {#[[2]], #[[3]]}& /@ lesHouchesInputParameters,
                            {{FileNameJoin[{Global`$flexiblesusyTemplateDir, spectrumGeneratorInputFile}],
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

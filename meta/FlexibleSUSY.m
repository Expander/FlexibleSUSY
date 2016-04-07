BeginPackage["FlexibleSUSY`", {"SARAH`", "AnomalousDimension`", "BetaFunction`", "TextFormatting`", "CConversion`", "TreeMasses`", "EWSB`", "Traces`", "SelfEnergies`", "Vertices`", "Phases`", "LoopMasses`", "WriteOut`", "Constraint`", "ThresholdCorrections`", "ConvergenceTester`", "Utils`", "ThreeLoopSM`", "ThreeLoopMSSM`", "Observables`"}];

$flexiblesusyMetaDir     = DirectoryName[FindFile[$Input]];
$flexiblesusyConfigDir   = FileNameJoin[{ParentDirectory[$flexiblesusyMetaDir], "config"}];
$flexiblesusyTemplateDir = FileNameJoin[{ParentDirectory[$flexiblesusyMetaDir], "templates"}];

FS`Version = StringTrim[FSImportString[FileNameJoin[{$flexiblesusyConfigDir,"version"}]]];
FS`GitCommit = StringTrim[FSImportString[FileNameJoin[{$flexiblesusyConfigDir,"git_commit"}]]];
FS`Authors = {"P. Athron", "J.-h. Park", "D. Stöckinger", "A. Voigt"};
FS`Years   = "2013-2015";
FS`References = Get[FileNameJoin[{$flexiblesusyConfigDir,"references"}]];

Print["==================================================================="];
Print[Style["FlexibleSUSY " <> FS`Version, Larger, Bold, Blue]];
Print["  by " <> Utils`StringJoinWithSeparator[FS`Authors, ", "] <> ", " <>
      FS`Years];
Print[""];
Print[Style["References:", Blue]];
Print["  " <> #]& /@ FS`References;
Print[""];
Print[Style["Download and Documentation:", Blue]];
Print["  https://flexiblesusy.hepforge.org"];
Print["==================================================================="];
Print[""];

Print["meta code directory: ", $flexiblesusyMetaDir];
Print["config directory   : ", $flexiblesusyConfigDir];
Print["templates directory: ", $flexiblesusyTemplateDir];

MakeFlexibleSUSY::usage="Creates a spectrum generator given a
 FlexibleSUSY model file (FlexibleSUSY.m).

Example:

  MakeFlexibleSUSY[
      InputFile -> \"models/<model>/FlexibleSUSY.m\",
      OutputDirectory -> \"models/<model>/\",
      DebugOutput -> False];

Options:

  InputFile: The name of the model file.

  OutputDirectory: The output directory for the generated code.

  DebugOutput (True|False): Enable/Disable debug output while running
    the Mathematica meta code.
";

LowPrecision::usage="";
MediumPrecision::usage="";
HighPrecision::usage="";
GUTNormalization::usage="Returns GUT normalization of a given coupling";

FSModelName;
FSOutputDir = ""; (* directory for generated code *)
FSLesHouchesList;
FSUnfixedParameters;
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
InitialGuessAtSUSYScale = {};
InitialGuessAtHighScale = {};
OnlyLowEnergyFlexibleSUSY = False;
SMTower = False;
SUSYScaleUserMatching={};
AutomaticInputAtMSUSY = True; (* input unfixed parameters at MSUSY *)
TreeLevelEWSBSolution = {};
Pole;
LowEnergyConstant;
FSMinimize;
FSFindRoot;
FSSolveEWSBFor;
MZ;
MT;
MZDRbar;
MWDRbar;
EDRbar;
ThetaWDRbar;
SCALE;
THRESHOLD;
UseHiggs2LoopNMSSM;
EffectiveMu;
EffectiveMASqr;
UseSM3LoopRGEs = False;
UseMSSM3LoopRGEs = False;
UseHiggs2LoopSM;
UseHiggs3LoopSplit;
PotentialLSPParticles = {};
ExtraSLHAOutputBlocks = {};
FSExtraInputParameters = {};

(* renormalization schemes *)
DRbar;
MSbar;
FSRenormalizationScheme = DRbar;

(* all model parameters are real by default *)
SARAH`RealParameters = { All };

(* precision of pole mass calculation *)
DefaultPoleMassPrecision = MediumPrecision;
HighPoleMassPrecision    = {SARAH`HiggsBoson, SARAH`PseudoScalar, SARAH`ChargedHiggs};
MediumPoleMassPrecision  = {};
LowPoleMassPrecision     = {};

FSEigenstates = SARAH`EWSB;
FSSolveEWSBTimeConstraint = 120;
FSSimplifyBetaFunctionsTimeConstraint = 120;
FSSolveWeinbergAngleTimeConstraint = 120;
FSCheckPerturbativityOfDimensionlessParameters = True;
FSPerturbativityThreshold = N[Sqrt[4 Pi]];
FSMaximumExpressionSize = 100;

(* list of soft breaking Higgs masses for solving EWSB eqs. *)
FSSoftHiggsMasses = {};

(* list of masses and parameters to check for convergence

   Example:

   FSConvergenceCheck = {
      M[hh], g3, Yu, Yd[3,3], Ye, B[\[Mu]]
   };
*)
FSConvergenceCheck = Automatic;

(* EWSB solvers *)
GSLHybrid;   (* hybrid method *)
GSLHybridS;  (* hybrid method with dynamic step size *)
GSLBroyden;  (* Broyden method *)
GSLNewton;   (* Newton method *)
FPIRelative; (* Fixed point iteration, convergence crit. relative step size *)
FPIAbsolute; (* Fixed point iteration, convergence crit. absolute step size *)
FPITadpole;  (* Fixed point iteration, convergence crit. relative step size + tadpoles *)
FSEWSBSolvers = { FPIRelative, GSLHybridS, GSLBroyden };

(* input value for the calculation of the weak mixing angle *)
FSFermiConstant;
FSMassW;

{FSHiggs, FSHyperchargeCoupling,
 FSLeftCoupling, FSStrongCoupling, FSVEVSM1, FSVEVSM2, FSNeutralino,
 FSChargino, FSNeutralinoMM, FSCharginoMinusMM, FSCharginoPlusMM,
 FSHiggsMM, FSSelectronL, FSSelectronNeutrinoL, FSSmuonL,
 FSSmuonNeutrinoL, FSVectorW, FSVectorZ, FSElectronYukawa};

FSWeakMixingAngleOptions = {
    FlexibleSUSY`FSWeakMixingAngleInput -> FSFermiConstant, (* or FSMassW *)
    FlexibleSUSY`FSHiggs                -> SARAH`HiggsBoson,
    FlexibleSUSY`FSHyperchargeCoupling  -> SARAH`hyperchargeCoupling,
    FlexibleSUSY`FSLeftCoupling         -> SARAH`leftCoupling,
    FlexibleSUSY`FSStrongCoupling       -> SARAH`strongCoupling,
    FlexibleSUSY`FSVEVSM1               -> SARAH`VEVSM1,
    FlexibleSUSY`FSVEVSM2               -> SARAH`VEVSM2,
    FlexibleSUSY`FSNeutralino           :> Parameters`GetParticleFromDescription["Neutralinos"],
    FlexibleSUSY`FSChargino             :> Parameters`GetParticleFromDescription["Charginos"],
    FlexibleSUSY`FSNeutralinoMM         -> SARAH`NeutralinoMM,
    FlexibleSUSY`FSCharginoMinusMM      -> SARAH`CharginoMinusMM,
    FlexibleSUSY`FSCharginoPlusMM       -> SARAH`CharginoPlusMM,
    FlexibleSUSY`FSHiggsMM              -> SARAH`HiggsMixingMatrix,
    FlexibleSUSY`FSSelectronL           :> Sum[Susyno`LieGroups`conj[SARAH`SleptonMM[Susyno`LieGroups`i,1]] SARAH`SleptonMM[Susyno`LieGroups`i,1] FlexibleSUSY`M[SARAH`Selectron[Susyno`LieGroups`i]], {Susyno`LieGroups`i,1,TreeMasses`GetDimension[SARAH`Selectron]}],
    FlexibleSUSY`FSSelectronNeutrinoL   :> Sum[Susyno`LieGroups`conj[SARAH`SneutrinoMM[Susyno`LieGroups`i,1]] SARAH`SneutrinoMM[Susyno`LieGroups`i,1] FlexibleSUSY`M[SARAH`Sneutrino[Susyno`LieGroups`i]], {Susyno`LieGroups`i,1,TreeMasses`GetDimension[SARAH`Sneutrino]}],
    FlexibleSUSY`FSSmuonL               :> Sum[Susyno`LieGroups`conj[SARAH`SleptonMM[Susyno`LieGroups`i,2]] SARAH`SleptonMM[Susyno`LieGroups`i,2] FlexibleSUSY`M[SARAH`Selectron[Susyno`LieGroups`i]], {Susyno`LieGroups`i,1,TreeMasses`GetDimension[SARAH`Selectron]}],
    FlexibleSUSY`FSSmuonNeutrinoL       :> Sum[Susyno`LieGroups`conj[SARAH`SneutrinoMM[Susyno`LieGroups`i,2]] SARAH`SneutrinoMM[Susyno`LieGroups`i,2] FlexibleSUSY`M[SARAH`Sneutrino[Susyno`LieGroups`i]], {Susyno`LieGroups`i,1,TreeMasses`GetDimension[SARAH`Sneutrino]}],
    FlexibleSUSY`FSVectorW              -> SARAH`VectorW,
    FlexibleSUSY`FSVectorZ              -> SARAH`VectorZ,
    FlexibleSUSY`FSElectronYukawa       -> SARAH`ElectronYukawa
};

ReadPoleMassPrecisions::ImpreciseHiggs="Warning: Calculating the Higgs pole mass M[`1`] with `2` will lead to an inaccurate result!  Please select MediumPrecision or HighPrecision (recommended) for `1`.";

tadpole::usage="symbolic expression for a tadpole contribution in the
EWSB eqs.  The index corresponds to the ordering of the tadpole
equations in SARAH`TadpoleEquations[] .";

FSDebugOutput = False;

Begin["`Private`"];

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
           minimRequiredVersionFile = FileNameJoin[{$flexiblesusyConfigDir,
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

CheckFermiConstantInputRequirements[requiredSymbols_List, printout_:True] :=
    Module[{resolvedSymbols, symbols, areDefined, availPars},
           resolvedSymbols = requiredSymbols /. FlexibleSUSY`FSWeakMixingAngleOptions;
           resolvedSymbols = resolvedSymbols /. {
               a_[idx__] :> a /; And @@ (NumberQ /@ {idx})
           };
           symbols = DeleteDuplicates[Cases[resolvedSymbols, _Symbol, {0,Infinity}]];
           availPars = Join[TreeMasses`GetParticles[],
                            Parameters`GetInputParameters[],
                            Parameters`GetModelParameters[],
                            Parameters`GetOutputParameters[]];
           areDefined = MemberQ[availPars, #]& /@ symbols;
           If[printout,
              Print["Unknown symbol: ", #]& /@
              Cases[Utils`Zip[areDefined, symbols], {False, p_} :> p];
             ];
           And @@ areDefined
          ];

CheckFermiConstantInputRequirementsForSUSYModel[] :=
    CheckFermiConstantInputRequirements[
        {FSHiggs, FSHyperchargeCoupling,
         FSLeftCoupling, FSStrongCoupling, FSVEVSM1, FSVEVSM2,
         FSNeutralino, FSChargino, FSNeutralinoMM, FSCharginoMinusMM,
         FSCharginoPlusMM, FSHiggsMM, FSSelectronL, FSSelectronNeutrinoL,
         FSSmuonL, FSSmuonNeutrinoL, FSVectorW, FSVectorZ,
         FSElectronYukawa}
    ];

CheckFermiConstantInputRequirementsForNonSUSYModel[] :=
    CheckFermiConstantInputRequirements[
        {FSHiggs, FSHyperchargeCoupling,
         FSLeftCoupling, FSStrongCoupling, FSVectorW, FSVectorZ,
         FSElectronYukawa}
    ];

CheckWeakMixingAngleInputRequirements[input_] :=
    Switch[input,
           FlexibleSUSY`FSFermiConstant,
               Switch[SARAH`SupersymmetricModel,
                      True,
                          If[CheckFermiConstantInputRequirementsForSUSYModel[],
                             input
                             ,
                             Print["Error: cannot use ", input, " because model"
                                   " requirements are not fulfilled"];
                             Print["   Using default input: ", FlexibleSUSY`FSMassW];
                             FlexibleSUSY`FSMassW
                          ],
                      False,
                          If[CheckFermiConstantInputRequirementsForNonSUSYModel[],
                             input
                             ,
                             Print["Error: cannot use ", input, " because model"
                                   " requirements are not fulfilled"];
                             Print["   Using default input: ", FlexibleSUSY`FSMassW];
                             FlexibleSUSY`FSMassW
                          ],
                      _,
                          Print["Error: model type: ", SARAH`SupersymmetricModel];
                          Print["   Using default input: ", FlexibleSUSY`FSMassW];
                          FlexibleSUSY`FSMassW
               ],
           FlexibleSUSY`FSMassW,
               input,
           _,
               Print["Error: unknown input ", input];
               Print["   Using default input: ", FlexibleSUSY`FSMassW];
               FlexibleSUSY`FSMassW
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
           If[Head[FlexibleSUSY`InitialGuessAtSUSYScale] =!= List,
              FlexibleSUSY`InitialGuessAtSUSYScale = {};
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
              FlexibleSUSY`LowScale := LowEnergyConstant[MZ];
             ];
           If[!ValueQ[FlexibleSUSY`LowScaleFirstGuess],
              Print["Warning: FlexibleSUSY`LowScaleFirstGuess should be",
                    " set in the model file!"];
              FlexibleSUSY`LowScaleFirstGuess = LowEnergyConstant[MZ];
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
           If[Head[FlexibleSUSY`FSExtraInputParameters] =!= List,
              Print["Error: FSExtraInputParameters has to be set to a list!"];
              Quit[1];
              ,
              If[!(And @@ (MatchQ[#,{_,_,_}]& /@ FlexibleSUSY`FSExtraInputParameters)),
                 Print["Error: FSExtraInputParameters must be of the form",
                       " {{A, AInput, {3,3}}, ... }"];
                ];
             ];
          ];

ReplaceIndicesInUserInput[rules_] :=
    Block[{},
          FlexibleSUSY`InitialGuessAtLowScale  = FlexibleSUSY`InitialGuessAtLowScale  /. rules;
          FlexibleSUSY`InitialGuessAtSUSYScale = FlexibleSUSY`InitialGuessAtSUSYScale /. rules;
          FlexibleSUSY`InitialGuessAtHighScale = FlexibleSUSY`InitialGuessAtHighScale /. rules;
          FlexibleSUSY`HighScale               = FlexibleSUSY`HighScale               /. rules;
          FlexibleSUSY`HighScaleFirstGuess     = FlexibleSUSY`HighScaleFirstGuess     /. rules;
          FlexibleSUSY`HighScaleInput          = FlexibleSUSY`HighScaleInput          /. rules;
          FlexibleSUSY`LowScale                = FlexibleSUSY`LowScale                /. rules;
          FlexibleSUSY`LowScaleFirstGuess      = FlexibleSUSY`LowScaleFirstGuess      /. rules;
          FlexibleSUSY`LowScaleInput           = FlexibleSUSY`LowScaleInput           /. rules;
          FlexibleSUSY`SUSYScale               = FlexibleSUSY`SUSYScale               /. rules;
          FlexibleSUSY`SUSYScaleFirstGuess     = FlexibleSUSY`SUSYScaleFirstGuess     /. rules;
          FlexibleSUSY`SUSYScaleInput          = FlexibleSUSY`SUSYScaleInput          /. rules;
         ];

EvaluateUserInput[] :=
    Block[{},
          FlexibleSUSY`HighScaleInput          = Map[Evaluate, FlexibleSUSY`HighScaleInput, {0,Infinity}];
          FlexibleSUSY`LowScaleInput           = Map[Evaluate, FlexibleSUSY`LowScaleInput , {0,Infinity}];
          FlexibleSUSY`SUSYScaleInput          = Map[Evaluate, FlexibleSUSY`SUSYScaleInput, {0,Infinity}];
         ];

GUTNormalization[coupling_] :=
    Parameters`GetGUTNormalization[coupling];

GeneralReplacementRules[] :=
    { "@VectorZ@"     -> ToValidCSymbolString[SARAH`VectorZ],
      "@VectorZ_" ~~ num_ ~~ "@" /; IntegerQ[ToExpression[num]] :> ToValidCSymbolString[SARAH`VectorZ] <> If[TreeMasses`GetDimension[SARAH`VectorZ] > 1, "(" <> num <> ")", ""],
      "@VectorZ()_" ~~ num_ ~~ "@" /; IntegerQ[ToExpression[num]] :> ToValidCSymbolString[SARAH`VectorZ] <> "()" <> If[TreeMasses`GetDimension[SARAH`VectorZ] > 1, "(" <> num <> ")", ""],
      "@VectorP@"     -> ToValidCSymbolString[SARAH`VectorP],
      "@VectorW@"     -> ToValidCSymbolString[SARAH`VectorW],
      "@VectorW_" ~~ num_ ~~ "@" /; IntegerQ[ToExpression[num]] :> ToValidCSymbolString[SARAH`VectorW] <> If[TreeMasses`GetDimension[SARAH`VectorW] > 1, "(" <> num <> ")", ""],
      "@VectorW()_" ~~ num_ ~~ "@" /; IntegerQ[ToExpression[num]] :> ToValidCSymbolString[SARAH`VectorW] <> "()" <> If[TreeMasses`GetDimension[SARAH`VectorW] > 1, "(" <> num <> ")", ""],
      "@VectorG@"     -> ToValidCSymbolString[SARAH`VectorG],
      "@TopQuark@"    -> ToValidCSymbolString[SARAH`TopQuark],
      "@BottomQuark@" -> ToValidCSymbolString[SARAH`BottomQuark],
      "@Electron@"    -> ToValidCSymbolString[SARAH`Electron],
      "@Neutrino@"    -> ToValidCSymbolString[SARAH`Neutrino],
      "@HiggsBoson@"  -> ToValidCSymbolString[SARAH`HiggsBoson],
      "@HiggsBoson_" ~~ num_ ~~ "@" /; IntegerQ[ToExpression[num]] :> ToValidCSymbolString[SARAH`HiggsBoson] <> If[TreeMasses`GetDimension[SARAH`HiggsBoson] > 1, "(" <> num <> ")", ""],
      "@HiggsBoson()_" ~~ num_ ~~ "@" /; IntegerQ[ToExpression[num]] :> ToValidCSymbolString[SARAH`HiggsBoson] <> "()"<> If[TreeMasses`GetDimension[SARAH`HiggsBoson] > 1, "(" <> num <> ")", ""],
      "@PseudoScalarBoson@" -> ToValidCSymbolString[SARAH`PseudoScalarBoson],
      "@ChargedHiggs@"   -> ToValidCSymbolString[SARAH`ChargedHiggs],
      "@TopSquark@"      -> ToValidCSymbolString[SARAH`TopSquark],
      "@TopSquark_" ~~ num_ ~~ "@" /; IntegerQ[ToExpression[num]] :> ToValidCSymbolString[SARAH`TopSquark] <> If[TreeMasses`GetDimension[SARAH`TopSquark] > 1, "(" <> num <> ")", ""],
      "@BottomSquark@"   -> ToValidCSymbolString[SARAH`BottomSquark],
      "@BottomSquark_" ~~ num_ ~~ "@" /; IntegerQ[ToExpression[num]] :> ToValidCSymbolString[SARAH`BottomSquark] <> If[TreeMasses`GetDimension[SARAH`BottomSquark] > 1, "(" <> num <> ")", ""],
      "@Sneutrino@"      -> ToValidCSymbolString[SARAH`Sneutrino],
      "@Selectron@"      -> ToValidCSymbolString[SARAH`Selectron],
      "@Gluino@"         -> ToValidCSymbolString[SARAH`Gluino],
      "@UpYukawa@"       -> ToValidCSymbolString[SARAH`UpYukawa],
      "@DownYukawa@"     -> ToValidCSymbolString[SARAH`DownYukawa],
      "@ElectronYukawa@" -> ToValidCSymbolString[SARAH`ElectronYukawa],
      "@LeftUpMixingMatrix@"   -> ToValidCSymbolString[SARAH`UpMatrixL],
      "@LeftDownMixingMatrix@" -> ToValidCSymbolString[SARAH`DownMatrixL],
      "@RightUpMixingMatrix@"  -> ToValidCSymbolString[SARAH`UpMatrixR],
      "@RightDownMixingMatrix@"-> ToValidCSymbolString[SARAH`DownMatrixR],
      "@hyperchargeCoupling@" -> ToValidCSymbolString[SARAH`hyperchargeCoupling],
      "@leftCoupling@"        -> ToValidCSymbolString[SARAH`leftCoupling],
      "@strongCoupling@"      -> ToValidCSymbolString[SARAH`strongCoupling],
      "@hyperchargeCouplingGutNormalization@"  -> RValueToCFormString[Parameters`GetGUTNormalization[SARAH`hyperchargeCoupling]],
      "@leftCouplingGutNormalization@"  -> RValueToCFormString[Parameters`GetGUTNormalization[SARAH`leftCoupling]],
      "@hyperchargeCouplingInverseGutNormalization@" -> RValueToCFormString[1/Parameters`GetGUTNormalization[SARAH`hyperchargeCoupling]],
      "@leftCouplingInverseGutNormalization@" -> RValueToCFormString[1/Parameters`GetGUTNormalization[SARAH`leftCoupling]],
      "@perturbativityThreshold@" -> ToString[N[FlexibleSUSY`FSPerturbativityThreshold]],
      "@ModelName@"           -> FlexibleSUSY`FSModelName,
      "@numberOfModelParameters@" -> ToString[numberOfModelParameters],
      "@InputParameter_" ~~ num_ ~~ "@" /; IntegerQ[ToExpression[num]] :> CConversion`ToValidCSymbolString[Parameters`GetInputParameters[][[ToExpression[num]]]],
      "@DateAndTime@"         -> DateString[],
      "@SARAHVersion@"        -> SA`Version,
      "@FlexibleSUSYVersion@" -> FS`Version,
      "@HiggsGen@" -> ToValidCSymbolString[GetDimension[SARAH`HiggsBoson]],
      "@FlexibleSUSYGitCommit@" -> FS`GitCommit
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
          calcTraces           = {Traces`CreateSARAHTraceCalculation[sarahTraces, "TRACE_STRUCT"],
                                  Sequence @@ Traces`CreateTraceCalculation[betaFun, "TRACE_STRUCT"] };
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
                   "@calc1LTraces@"         -> IndentText @ IndentText[WrapLines[calcTraces[[1]] <> "\n" <> calcTraces[[2]]]],
                   "@calc2LTraces@"         -> IndentText @ IndentText[WrapLines[calcTraces[[3]]]],
                   "@calc3LTraces@"         -> IndentText @ IndentText[WrapLines[calcTraces[[4]]]],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          singleBetaFunctionsDefsFiles = BetaFunction`CreateSingleBetaFunctionDefs[betaFun, templateFile, sarahTraces];
          Print["Creating makefile module for the two-scale method ..."];
          WriteMakefileModule[singleBetaFunctionsDefsFiles,
                              makefileModuleTemplates];
         ];

WriteInputParameterClass[inputParameters_List, files_List] :=
   Module[{defineInputParameters, defaultInputParametersInit, printInputParameters, get, set},
          defineInputParameters = Constraint`DefineInputParameters[inputParameters];
          defaultInputParametersInit = Constraint`InitializeInputParameters[inputParameters];
          printInputParameters = WriteOut`PrintInputParameters[inputParameters,"ostr"];
          get = Parameters`CreateInputParameterArrayGetter[inputParameters];
          set = Parameters`CreateInputParameterArraySetter[inputParameters];
          WriteOut`ReplaceInFiles[files,
                         { "@defineInputParameters@" -> IndentText[defineInputParameters],
                           "@defaultInputParametersInit@" -> WrapLines[defaultInputParametersInit],
                           "@printInputParameters@"       -> IndentText[printInputParameters],
                           "@get@"                        -> IndentText[get],
                           "@set@"                        -> IndentText[set],
                           Sequence @@ GeneralReplacementRules[]
                         } ];
          ];

WriteConstraintClass[condition_, settings_List, scaleFirstGuess_,
                     {minimumScale_, maximumScale_}, files_List] :=
   Module[{applyConstraint = "", calculateScale, scaleGuess,
           restrictScale,
           setDRbarYukawaCouplings,
           calculateDRbarMasses,
           calculateDeltaAlphaEm, calculateDeltaAlphaS,
           calculateGaugeCouplings,
           calculateThetaW,
           recalculateMWPole,
           checkPerturbativityForDimensionlessParameters = "",
           saveEwsbOutputParameters, restoreEwsbOutputParameters},
          Constraint`SetBetaFunctions[GetBetaFunctions[]];
          applyConstraint = Constraint`ApplyConstraints[settings];
          calculateScale  = Constraint`CalculateScale[condition, "scale"];
          scaleGuess      = Constraint`CalculateScale[scaleFirstGuess, "initial_scale_guess"];
          restrictScale   = Constraint`RestrictScale[{minimumScale, maximumScale}];
          calculateDeltaAlphaEm   = ThresholdCorrections`CalculateDeltaAlphaEm[FlexibleSUSY`FSRenormalizationScheme];
          calculateDeltaAlphaS    = ThresholdCorrections`CalculateDeltaAlphaS[FlexibleSUSY`FSRenormalizationScheme];
          calculateThetaW         = ThresholdCorrections`CalculateThetaW[FSWeakMixingAngleOptions,SARAH`SupersymmetricModel];
          calculateGaugeCouplings = ThresholdCorrections`CalculateGaugeCouplings[];
          recalculateMWPole       = ThresholdCorrections`RecalculateMWPole[FSWeakMixingAngleOptions];
          setDRbarYukawaCouplings = {
              ThresholdCorrections`SetDRbarYukawaCouplingTop[settings],
              ThresholdCorrections`SetDRbarYukawaCouplingBottom[settings],
              ThresholdCorrections`SetDRbarYukawaCouplingElectron[settings]
          };
          calculateDRbarMasses = {
              LoopMasses`CallCalculateDRbarMass["Up Quark"         , "Up-Quarks"  , 1, "topDRbar"     , "qedqcd.displayMass(softsusy::mUp)"      ],
              LoopMasses`CallCalculateDRbarMass["Charmed Quark"    , "Up-Quarks"  , 2, "topDRbar"     , "qedqcd.displayMass(softsusy::mCharm)"   ],
              LoopMasses`CallCalculateDRbarMass["Top Quark"        , "Up-Quarks"  , 3, "topDRbar"     , "qedqcd.displayPoleMt()"                 ],
              LoopMasses`CallCalculateDRbarMass["Down Quark"       , "Down-Quarks", 1, "bottomDRbar"  , "qedqcd.displayMass(softsusy::mDown)"    ],
              LoopMasses`CallCalculateDRbarMass["Strange Quark"    , "Down-Quarks", 2, "bottomDRbar"  , "qedqcd.displayMass(softsusy::mStrange)" ],
              LoopMasses`CallCalculateDRbarMass["Bottom Quark"     , "Down-Quarks", 3, "bottomDRbar"  , "qedqcd.displayMass(softsusy::mBottom)"  ],
              LoopMasses`CallCalculateDRbarMass["Electron"         , "Leptons"    , 1, "electronDRbar", "qedqcd.displayMass(softsusy::mElectron)"],
              LoopMasses`CallCalculateDRbarMass["Muon"             , "Leptons"    , 2, "electronDRbar", "qedqcd.displayMass(softsusy::mMuon)"    ],
              LoopMasses`CallCalculateDRbarMass["Tau"              , "Leptons"    , 3, "electronDRbar", "qedqcd.displayMass(softsusy::mTau)"     ],
              LoopMasses`CallCalculateDRbarMass["Electron Neutrino", "Neutrinos"  , 1, "neutrinoDRbar", "qedqcd.displayNeutrinoPoleMass(1)"      ],
              LoopMasses`CallCalculateDRbarMass["Muon Neutrino"    , "Neutrinos"  , 2, "neutrinoDRbar", "qedqcd.displayNeutrinoPoleMass(2)"      ],
              LoopMasses`CallCalculateDRbarMass["Tau Neutrino"     , "Neutrinos"  , 3, "neutrinoDRbar", "qedqcd.displayNeutrinoPoleMass(3)"      ]
          };
          saveEwsbOutputParameters    = Parameters`SaveParameterLocally[FlexibleSUSY`EWSBOutputParameters, "old_", "MODELPARAMETER"];
          restoreEwsbOutputParameters = Parameters`RestoreParameter[FlexibleSUSY`EWSBOutputParameters, "old_", "model"];
          If[FSCheckPerturbativityOfDimensionlessParameters,
             checkPerturbativityForDimensionlessParameters =
                 Constraint`CheckPerturbativityForParameters[
                     Parameters`GetModelParametersWithMassDimension[0],
                     FlexibleSUSY`FSPerturbativityThreshold
                 ];
            ];
          WriteOut`ReplaceInFiles[files,
                 { "@applyConstraint@"      -> IndentText[WrapLines[applyConstraint]],
                   "@calculateScale@"       -> IndentText[WrapLines[calculateScale]],
                   "@scaleGuess@"           -> IndentText[WrapLines[scaleGuess]],
                   "@restrictScale@"        -> IndentText[WrapLines[restrictScale]],
                   "@calculateGaugeCouplings@" -> IndentText[WrapLines[calculateGaugeCouplings]],
                   "@calculateDeltaAlphaEm@" -> IndentText[WrapLines[calculateDeltaAlphaEm]],
                   "@calculateDeltaAlphaS@"  -> IndentText[WrapLines[calculateDeltaAlphaS]],
                   "@calculateThetaW@"       -> IndentText[WrapLines[calculateThetaW]],
                   "@recalculateMWPole@"     -> IndentText[WrapLines[recalculateMWPole]],
                   "@calculateDRbarMassUp@"      -> IndentText[IndentText[calculateDRbarMasses[[1]]]],
                   "@calculateDRbarMassCharm@"   -> IndentText[IndentText[calculateDRbarMasses[[2]]]],
                   "@calculateDRbarMassTop@"     -> IndentText[IndentText[calculateDRbarMasses[[3]]]],
                   "@calculateDRbarMassDown@"    -> IndentText[IndentText[calculateDRbarMasses[[4]]]],
                   "@calculateDRbarMassStrange@" -> IndentText[IndentText[calculateDRbarMasses[[5]]]],
                   "@calculateDRbarMassBottom@"  -> IndentText[IndentText[calculateDRbarMasses[[6]]]],
                   "@calculateDRbarMassElectron@"-> IndentText[IndentText[calculateDRbarMasses[[7]]]],
                   "@calculateDRbarMassMuon@"    -> IndentText[IndentText[calculateDRbarMasses[[8]]]],
                   "@calculateDRbarMassTau@"     -> IndentText[IndentText[calculateDRbarMasses[[9]]]],
                   "@calculateDRbarMassElectronNeutrino@"-> IndentText[IndentText[calculateDRbarMasses[[10]]]],
                   "@calculateDRbarMassMuonNeutrino@"    -> IndentText[IndentText[calculateDRbarMasses[[11]]]],
                   "@calculateDRbarMassTauNeutrino@"     -> IndentText[IndentText[calculateDRbarMasses[[12]]]],
                   "@setDRbarUpQuarkYukawaCouplings@"   -> IndentText[WrapLines[setDRbarYukawaCouplings[[1]]]],
                   "@setDRbarDownQuarkYukawaCouplings@" -> IndentText[WrapLines[setDRbarYukawaCouplings[[2]]]],
                   "@setDRbarElectronYukawaCouplings@"  -> IndentText[WrapLines[setDRbarYukawaCouplings[[3]]]],
                   "@saveEwsbOutputParameters@"    -> IndentText[saveEwsbOutputParameters],
                   "@restoreEwsbOutputParameters@" -> IndentText[restoreEwsbOutputParameters],
                   "@checkPerturbativityForDimensionlessParameters@" -> IndentText[checkPerturbativityForDimensionlessParameters],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteInitialGuesserClass[lowScaleGuess_List, susyScaleGuess_List, highScaleGuess_List, files_List] :=
   Module[{initialGuessAtLowScale, initialGuessAtHighScale, setDRbarYukawaCouplings,
           allSettings},
          initialGuessAtLowScale  = Constraint`ApplyConstraints[lowScaleGuess];
          initialGuessAtSUSYScale = Constraint`ApplyConstraints[susyScaleGuess];
          initialGuessAtHighScale = Constraint`ApplyConstraints[highScaleGuess];
          allSettings             = Join[lowScaleGuess, highScaleGuess];
          setDRbarYukawaCouplings = {
              ThresholdCorrections`SetDRbarYukawaCouplingTop[allSettings],
              ThresholdCorrections`SetDRbarYukawaCouplingBottom[allSettings],
              ThresholdCorrections`SetDRbarYukawaCouplingElectron[allSettings]
          };
          WriteOut`ReplaceInFiles[files,
                 { "@initialGuessAtLowScale@"  -> IndentText[WrapLines[initialGuessAtLowScale]],
                   "@initialGuessAtSUSYScale@" -> IndentText[WrapLines[initialGuessAtSUSYScale]],
                   "@initialGuessAtHighScale@" -> IndentText[WrapLines[initialGuessAtHighScale]],
                   "@setDRbarUpQuarkYukawaCouplings@"   -> IndentText[WrapLines[setDRbarYukawaCouplings[[1]]]],
                   "@setDRbarDownQuarkYukawaCouplings@" -> IndentText[WrapLines[setDRbarYukawaCouplings[[2]]]],
                   "@setDRbarElectronYukawaCouplings@"  -> IndentText[WrapLines[setDRbarYukawaCouplings[[3]]]],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteConvergenceTesterClass[parameters_, files_List] :=
   Module[{compareFunction},
          compareFunction = ConvergenceTester`CreateCompareFunction[parameters];
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

(* Returns a list of three-component lists where the information is
   stored which Higgs corresponds to which EWSB eq. and whether the
   corresponding tadpole is real or imaginary (only in models with CP
   violation).

   Example: MRSSM
   In[] := CreateHiggsToEWSBEqAssociation[]
   Out[] = {{hh, 1, Re}, {hh, 2, Re}, {hh, 4, Re}, {hh, 3, Re}}

   This result means:

   EWSB eq. 1 corresponds to hh[1], the 1L tadpole[1] is real
   EWSB eq. 2 corresponds to hh[2], the 1L tadpole[2] is real
   EWSB eq. 3 corresponds to hh[4], the 1L tadpole[3] is real
   EWSB eq. 4 corresponds to hh[3], the 1L tadpole[4] is real
 *)
CreateHiggsToEWSBEqAssociation[] :=
    Module[{association = {}, v, phi, sigma, higgs, numberOfVEVs, numberOfHiggses, vevs,
            vev, dimVev},
           vevs = Cases[SARAH`DEFINITION[FlexibleSUSY`FSEigenstates][SARAH`VEVs],
                        {_,{v_,_},{s_,_},{p_,_},___} :> {v,s,p}];
           If[Length[vevs] == 1,
              Return[{{SARAH`HiggsBoson, 1, Re}}];
             ];
           FindPositions[es_] :=
               Module[{gaugeES, higgsGaugeES},
                      gaugeES = ExpandGaugeIndices[es];
                      (* list of gauge eigenstate fields, ordered according to Higgs mixing *)
                      higgsGaugeES = Cases[SARAH`DEFINITION[FlexibleSUSY`FSEigenstates][SARAH`MatterSector],
                                           {gauge_List, {SARAH`HiggsBoson, _}} :> gauge][[1]];
                      higgsGaugeES = ExpandGaugeIndices[higgsGaugeES];
                      (* find positions of gaugeES in higgsGaugeES *)
                      {SARAH`HiggsBoson,#}& /@ (Flatten[Position[higgsGaugeES, #]& /@ gaugeES])
                     ];
           Join[Append[#,Re]& /@ FindPositions[Transpose[vevs][[3]]],
                Append[#,Re]& /@ FindPositions[Transpose[vevs][[2]]]]
          ];

WriteModelSLHAClass[massMatrices_List, files_List] :=
    Module[{k,
            slhaYukawaDef = "",
            slhaYukawaGetter = "",
            convertYukawaCouplingsToSLHA = "",
            slhaTrilinearCouplingsDef = "",
            slhaTrilinearCouplingsGetter = "",
            convertTrilinearCouplingsToSLHA = "",
            slhaSoftSquaredMassesDef = "",
            slhaSoftSquaredMassesGetter = "",
            convertSoftSquaredMassesToSLHA = "",
            slhaFerimonMixingMatricesDef = "",
            slhaFerimonMixingMatricesGetters = "",
            slhaPoleMassGetters = "",
            slhaPoleMixingMatrixGetters = "",
            calculateCKMMatrix = "",
            calculatePMNSMatrix = ""
           },
           slhaYukawaDef        = WriteOut`CreateSLHAYukawaDefinition[];
           slhaYukawaGetter     = WriteOut`CreateSLHAYukawaGetters[];
           convertYukawaCouplingsToSLHA = WriteOut`ConvertYukawaCouplingsToSLHA[];
           slhaTrilinearCouplingsDef    = WriteOut`CreateSLHATrilinearCouplingDefinition[];
           slhaTrilinearCouplingsGetter = WriteOut`CreateSLHATrilinearCouplingGetters[];
           convertTrilinearCouplingsToSLHA = WriteOut`ConvertTrilinearCouplingsToSLHA[];
           slhaSoftSquaredMassesDef    = WriteOut`CreateSLHASoftSquaredMassesDefinition[];
           slhaSoftSquaredMassesGetter = WriteOut`CreateSLHASoftSquaredMassesGetters[];
           convertSoftSquaredMassesToSLHA = WriteOut`ConvertSoftSquaredMassesToSLHA[];
           slhaFerimonMixingMatricesDef = WriteOut`CreateSLHAFermionMixingMatricesDef[];
           slhaFerimonMixingMatricesGetters = WriteOut`CreateSLHAFermionMixingMatricesGetters[];
           calculateCKMMatrix = WriteOut`CalculateCKMMatrix[];
           calculatePMNSMatrix = WriteOut`CalculatePMNSMatrix[];
           For[k = 1, k <= Length[massMatrices], k++,
               slhaPoleMassGetters         = slhaPoleMassGetters <> TreeMasses`CreateSLHAPoleMassGetter[massMatrices[[k]]];
               slhaPoleMixingMatrixGetters = slhaPoleMixingMatrixGetters <> TreeMasses`CreateSLHAPoleMixingMatrixGetter[massMatrices[[k]]];
              ];
           WriteOut`ReplaceInFiles[files,
                          { "@slhaYukawaDef@"                  -> IndentText[slhaYukawaDef],
                            "@slhaYukawaGetter@"               -> IndentText[slhaYukawaGetter],
                            "@convertYukawaCouplingsToSLHA@"   -> IndentText[convertYukawaCouplingsToSLHA],
                            "@slhaFerimonMixingMatricesDef@"   -> IndentText[slhaFerimonMixingMatricesDef],
                            "@slhaFerimonMixingMatricesGetters@" -> IndentText[slhaFerimonMixingMatricesGetters],
                            "@slhaTrilinearCouplingsDef@"      -> IndentText[slhaTrilinearCouplingsDef],
                            "@slhaTrilinearCouplingsGetter@"   -> IndentText[slhaTrilinearCouplingsGetter],
                            "@convertTrilinearCouplingsToSLHA@"-> IndentText[convertTrilinearCouplingsToSLHA],
                            "@slhaSoftSquaredMassesDef@"       -> IndentText[slhaSoftSquaredMassesDef],
                            "@slhaSoftSquaredMassesGetter@"    -> IndentText[slhaSoftSquaredMassesGetter],
                            "@convertSoftSquaredMassesToSLHA@" -> IndentText[convertSoftSquaredMassesToSLHA],
                            "@slhaPoleMassGetters@"            -> IndentText[slhaPoleMassGetters],
                            "@slhaPoleMixingMatrixGetters@"    -> IndentText[slhaPoleMixingMatrixGetters],
                            "@calculateCKMMatrix@"             -> IndentText[calculateCKMMatrix],
                            "@calculatePMNSMatrix@"             -> IndentText[calculatePMNSMatrix],
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

(* Returns a list of three-component lists where the information is
   stored which VEV corresponds to which Tadpole eq.

   Example: MRSSM
   It[] := CreateVEVToTadpoleAssociation[]
   Out[] = {{hh, 1, vd}, {hh, 2, vu}, {hh, 4, vS}, {hh, 3, vT}}
 *)
CreateVEVToTadpoleAssociation[] :=
    Module[{association, vev},
           vevs = Cases[SARAH`DEFINITION[FlexibleSUSY`FSEigenstates][SARAH`VEVs],
                        {_,{v_,_},{s_,_},{p_,_},___} :> {v,s,p}];
           association = CreateHiggsToEWSBEqAssociation[];
           {#[[1]], #[[2]], vevs[[#[[2]],1]]}& /@ association
          ];

WriteMatchingClass[files_List] :=
    Module[ {scheme, userMatching},
		scheme = If[SARAH`SupersymmetricModel, FlexibleSUSY`DRbar, FlexibleSUSY`MSbar];
        userMatching="";
        If[SMTower,
           If[Length[SUSYScaleUserMatching] > 0,
              For[i = 1, i <= Length[SUSYScaleUserMatching], i++,
                 userMatching = userMatching <> "model.set_" <> CConversion`ToValidCSymbolString[SUSYScaleUserMatching[[i,1]]]
                                <> "(" <> CConversion`RValueToCFormString[SUSYScaleUserMatching[[i, 2]]]
                                <> ");\n";
              ];
           ];
        ];
        WriteOut`ReplaceInFiles[files, {
          "@gauge1Linit@"       -> IndentText[WrapLines[ Parameters`CreateLocalConstRefs[
                                                                    ThresholdCorrections`CalculateColorCoupling[scheme]
                                                                  + ThresholdCorrections`CalculateElectromagneticCoupling[scheme]
                                                                  ]]],
          "@alphaS1Lmatching@"  ->  IndentText[WrapLines["const double delta_alpha_s = alpha_s/(2.*Pi)*(" <>
                                                         CConversion`RValueToCFormString[ThresholdCorrections`CalculateColorCoupling[scheme]] <> ");\n"]],
          "@alphaEM1Lmatching@" ->  IndentText[WrapLines["const double delta_alpha_em = alpha_em/(2.*Pi)*(" <>
                                                         CConversion`RValueToCFormString[ThresholdCorrections`CalculateElectromagneticCoupling[scheme]] <> ");\n"]],
          "@setYukawas@" -> IndentText[WrapLines[ThresholdCorrections`SetDRbarYukawaCouplings[]]],
          "@applyUserMatching@" -> IndentText[WrapLines[userMatching]],
          Sequence @@ GeneralReplacementRules[]}
                               ];
             ];


WriteModelClass[massMatrices_List, ewsbEquations_List,
                parametersFixedByEWSB_List, ewsbSolution_List, freePhases_List,
                nPointFunctions_List, vertexRules_List, phases_List,
                files_List, diagonalizationPrecision_List] :=
    Module[{ewsbEquationsTreeLevel, independentEwsbEquationsTreeLevel,
            independentEwsbEquations,
            massGetters = "", k,
            mixingMatrixGetters = "",
            slhaPoleMassGetters = "", slhaPoleMixingMatrixGetters = "",
            higgsMassGetters = "",
            tadpoleEqPrototypes = "", tadpoleEqFunctions = "",
            numberOfEWSBEquations = Length[ewsbEquations],
            numberOfIndependentEWSBEquations,
            calculateTreeLevelTadpoles = "",
            ewsbInitialGuess = "", physicalMassesDef = "", mixingMatricesDef = "",
            physicalMassesInit = "", physicalMassesInitNoLeadingComma = "", mixingMatricesInit = "",
            massCalculationPrototypes = "", massCalculationFunctions = "",
            calculateAllMasses = "",
            calculateOneLoopTadpoles = "", calculateTwoLoopTadpoles = "",
            calculateOneLoopTadpolesNoStruct = "", calculateTwoLoopTadpolesNoStruct = "",
            selfEnergyPrototypes = "", selfEnergyFunctions = "",
            twoLoopTadpolePrototypes = "", twoLoopTadpoleFunctions = "",
            twoLoopSelfEnergyPrototypes = "", twoLoopSelfEnergyFunctions = "",
            threeLoopSelfEnergyPrototypes = "", threeLoopSelfEnergyFunctions = "",
            thirdGenerationHelperPrototypes = "", thirdGenerationHelperFunctions = "",
            phasesDefinition = "", phasesGetterSetters = "",
            phasesInit = "",
            loopMassesPrototypes = "", loopMassesFunctions = "",
            runningDRbarMassesPrototypes = "", runningDRbarMassesFunctions = "",
            callAllLoopMassFunctions = "",
            callAllLoopMassFunctionsInThreads = "",
            printMasses = "", printMixingMatrices = "",
            getPhysical = "", setPhysical = "",
            getMasses = "", setMasses = "",
            masses, mixingMatrices, oneLoopTadpoles,
            dependenceNumPrototypes, dependenceNumFunctions,
            clearOutputParameters = "", solveEwsbTreeLevel = "",
            clearPhases = "",
            saveEwsbOutputParameters, restoreEwsbOutputParameters,
            softScalarMasses, treeLevelEWSBOutputParameters,
            saveSoftHiggsMasses, restoreSoftHiggsMasses,
            solveTreeLevelEWSBviaSoftHiggsMasses,
            solveEWSBTemporarily,
            copyDRbarMassesToPoleMasses = "",
            reorderDRbarMasses = "", reorderPoleMasses = "",
            checkPoleMassesForTachyons = "",
            higgsToEWSBEqAssociation,
            twoLoopHiggsHeaders = "", threeLoopHiggsHeaders = "",
            lspGetters = "", lspFunctions = "",
            EWSBSolvers = "",
            setEWSBSolution = "",
            fillArrayWithEWSBParameters = "",
            solveEwsbWithTadpoles = "",
            getEWSBParametersFromGSLVector = "",
            setEWSBParametersFromLocalCopies = "",
            ewsbParametersInitializationList = "",
            setEWSBParametersFromGSLVector = "",
            convertMixingsToSLHAConvention = "",
            convertMixingsToHKConvention = "",
            enablePoleMassThreads = True
           },
           convertMixingsToSLHAConvention = WriteOut`ConvertMixingsToSLHAConvention[massMatrices];
           convertMixingsToHKConvention   = WriteOut`ConvertMixingsToHKConvention[massMatrices];
           independentEwsbEquations = Parameters`FilterOutLinearDependentEqs[ewsbEquations, parametersFixedByEWSB];
           numberOfIndependentEWSBEquations = Length[independentEwsbEquations];
           ewsbEquationsTreeLevel = ewsbEquations /. FlexibleSUSY`tadpole[_] -> 0;
           independentEwsbEquationsTreeLevel = independentEwsbEquations /. FlexibleSUSY`tadpole[_] -> 0;
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
           higgsMassGetters =
               Utils`StringZipWithSeparator[
                   TreeMasses`CreateHiggsMassGetters[SARAH`HiggsBoson,""],
                   TreeMasses`CreateHiggsMassGetters[SARAH`ChargedHiggs,""],
                   TreeMasses`CreateHiggsMassGetters[SARAH`PseudoScalar,""],
                   "\n"
               ];
           clearPhases = Phases`ClearPhases[phases];
           calculateAllMasses = TreeMasses`CallMassCalculationFunctions[massMatrices];
           tadpoleEqPrototypes = EWSB`CreateEWSBEqPrototype[SARAH`HiggsBoson];
           tadpoleEqFunctions  = EWSB`CreateEWSBEqFunction[SARAH`HiggsBoson, ewsbEquationsTreeLevel];
           If[ewsbEquations =!= Table[0, {numberOfEWSBEquations}] &&
              Length[parametersFixedByEWSB] != numberOfIndependentEWSBEquations,
              Print["Error: There are ", numberOfIndependentEWSBEquations, " independent EWSB ",
                    "equations, but you want to fix ", Length[parametersFixedByEWSB],
                    " parameters: ", parametersFixedByEWSB];
             ];
           higgsToEWSBEqAssociation     = CreateHiggsToEWSBEqAssociation[];
           oneLoopTadpoles              = Cases[nPointFunctions, SelfEnergies`Tadpole[___]];
           calculateOneLoopTadpoles     = SelfEnergies`FillArrayWithOneLoopTadpoles[higgsToEWSBEqAssociation, "tadpole", "-"];
           calculateOneLoopTadpolesNoStruct = SelfEnergies`FillArrayWithOneLoopTadpoles[higgsToEWSBEqAssociation, "tadpole", "+"];
           If[SARAH`UseHiggs2LoopMSSM === True ||
              FlexibleSUSY`UseHiggs2LoopNMSSM === True,
              calculateTwoLoopTadpoles  = SelfEnergies`FillArrayWithTwoLoopTadpoles[SARAH`HiggsBoson, "tadpole", "-"];
              calculateTwoLoopTadpolesNoStruct = SelfEnergies`FillArrayWithTwoLoopTadpoles[SARAH`HiggsBoson, "tadpole", "+"];
              {thirdGenerationHelperPrototypes, thirdGenerationHelperFunctions} = TreeMasses`CreateThirdGenerationHelpers[];
             ];
           If[FlexibleSUSY`UseHiggs2LoopSM === True,
              {twoLoopSelfEnergyPrototypes, twoLoopSelfEnergyFunctions} = SelfEnergies`CreateTwoLoopSelfEnergiesSM[{SARAH`HiggsBoson}];
              twoLoopHiggsHeaders = "#include \"sm_twoloophiggs.hpp\"\n";
             ];
           If[FlexibleSUSY`UseHiggs3LoopSplit === True,
              {threeLoopSelfEnergyPrototypes, threeLoopSelfEnergyFunctions} = SelfEnergies`CreateThreeLoopSelfEnergiesSplit[{SARAH`HiggsBoson}];
              threeLoopHiggsHeaders = "#include \"split_threeloophiggs.hpp\"\n";
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
           setEWSBParametersFromGSLVector = EWSB`SetEWSBParametersFromGSLVector[parametersFixedByEWSB, freePhases, "x"];
           calculateTreeLevelTadpoles   = EWSB`FillArrayWithEWSBEqs[SARAH`HiggsBoson, "tadpole"];
           ewsbInitialGuess             = EWSB`FillInitialGuessArray[parametersFixedByEWSB];
           solveEwsbTreeLevel           = EWSB`CreateTreeLevelEwsbSolver[ewsbSolution /. FlexibleSUSY`tadpole[_] -> 0];
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
           runningDRbarMassesFunctions  = LoopMasses`CreateRunningDRbarMassFunctions[FlexibleSUSY`FSRenormalizationScheme];
           enablePoleMassThreads = False;
           callAllLoopMassFunctions     = LoopMasses`CallAllPoleMassFunctions[FlexibleSUSY`FSEigenstates, enablePoleMassThreads];
           enablePoleMassThreads = True;
           callAllLoopMassFunctionsInThreads = LoopMasses`CallAllPoleMassFunctions[FlexibleSUSY`FSEigenstates, enablePoleMassThreads];
           masses                       = FlexibleSUSY`M[TreeMasses`GetMassEigenstate[#]]& /@ massMatrices;
           {lspGetters, lspFunctions}   = LoopMasses`CreateLSPFunctions[FlexibleSUSY`PotentialLSPParticles];
           printMasses                  = WriteOut`PrintParameters[masses, "ostr"];
           getPhysical                  = TreeMasses`CreatePhysicalArrayGetter[massMatrices];
           setPhysical                  = TreeMasses`CreatePhysicalArraySetter[massMatrices, "pars"];
           getMasses                    = TreeMasses`CreateMassArrayGetter[massMatrices];
           setMasses                    = TreeMasses`CreateMassArraySetter[massMatrices, "pars"];
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
           If[Head[softHiggsMasses] === List && Length[softHiggsMasses] > 0,
              saveSoftHiggsMasses       = Parameters`SaveParameterLocally[softHiggsMasses, "old_", ""];
              restoreSoftHiggsMasses    = Parameters`RestoreParameter[softHiggsMasses, "old_", ""];
              solveTreeLevelEWSBviaSoftHiggsMasses = EWSB`SolveTreeLevelEwsbVia[independentEwsbEquationsTreeLevel, softHiggsMasses];
              solveEWSBTemporarily = "solve_ewsb_tree_level_via_soft_higgs_masses();";
              ,
              saveSoftHiggsMasses       = Parameters`SaveParameterLocally[FlexibleSUSY`EWSBOutputParameters, "old_", ""];
              restoreSoftHiggsMasses    = Parameters`RestoreParameter[FlexibleSUSY`EWSBOutputParameters, "old_", ""];
              solveTreeLevelEWSBviaSoftHiggsMasses = "";
              solveEWSBTemporarily = "solve_ewsb_tree_level();";
             ];
           EWSBSolvers                  = EWSB`CreateEWSBRootFinders[FlexibleSUSY`FSEWSBSolvers];
           setEWSBSolution              = EWSB`SetEWSBSolution[parametersFixedByEWSB, freePhases, "solver->get_solution"];
           fillArrayWithEWSBParameters  = EWSB`FillArrayWithParameters["ewsb_parameters", parametersFixedByEWSB];
           solveEwsbWithTadpoles        = EWSB`CreateEwsbSolverWithTadpoles[ewsbSolution, treeLevelEWSBOutputParameters];
           getEWSBParametersFromGSLVector = EWSB`GetEWSBParametersFromGSLVector[parametersFixedByEWSB, freePhases, "x"];
           setEWSBParametersFromLocalCopies = EWSB`SetEWSBParametersFromLocalCopies[parametersFixedByEWSB, "model"];
           ewsbParametersInitializationList = EWSB`CreateEWSBParametersInitializationList[parametersFixedByEWSB];
           reorderDRbarMasses           = TreeMasses`ReorderGoldstoneBosons[""];
           reorderPoleMasses            = TreeMasses`ReorderGoldstoneBosons["PHYSICAL"];
           checkPoleMassesForTachyons   = TreeMasses`CheckPoleMassesForTachyons["PHYSICAL"];
           WriteOut`ReplaceInFiles[files,
                          { "@lspGetters@"           -> IndentText[lspGetters],
                            "@lspFunctions@"         -> lspFunctions,
                            "@massGetters@"          -> IndentText[massGetters],
                            "@mixingMatrixGetters@"  -> IndentText[mixingMatrixGetters],
                            "@slhaPoleMassGetters@"  -> IndentText[slhaPoleMassGetters],
                            "@slhaPoleMixingMatrixGetters@" -> IndentText[slhaPoleMixingMatrixGetters],
                            "@higgsMassGetterPrototypes@"   -> IndentText[higgsMassGetters[[1]]],
                            "@higgsMassGetters@"     -> higgsMassGetters[[2]],
                            "@tadpoleEqPrototypes@"  -> IndentText[tadpoleEqPrototypes],
                            "@tadpoleEqFunctions@"   -> tadpoleEqFunctions,
                            "@numberOfEWSBEquations@"-> ToString[TreeMasses`GetDimension[SARAH`HiggsBoson]],
                            "@calculateTreeLevelTadpoles@" -> IndentText[calculateTreeLevelTadpoles],
                            "@calculateOneLoopTadpoles@"   -> IndentText[calculateOneLoopTadpoles],
                            "@calculateTwoLoopTadpoles@"   -> IndentText[calculateTwoLoopTadpoles],
                            "@calculateOneLoopTadpolesNoStruct@" -> IndentText[calculateOneLoopTadpolesNoStruct],
                            "@calculateTwoLoopTadpolesNoStruct@" -> IndentText[calculateTwoLoopTadpolesNoStruct],
                            "@clearOutputParameters@"  -> IndentText[clearOutputParameters],
                            "@clearPhases@"            -> IndentText[clearPhases],
                            "@copyDRbarMassesToPoleMasses@" -> IndentText[copyDRbarMassesToPoleMasses],
                            "@reorderDRbarMasses@"     -> IndentText[reorderDRbarMasses],
                            "@reorderPoleMasses@"      -> IndentText[reorderPoleMasses],
                            "@checkPoleMassesForTachyons@" -> IndentText[checkPoleMassesForTachyons],
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
                            "@threeLoopSelfEnergyPrototypes@" -> IndentText[threeLoopSelfEnergyPrototypes],
                            "@threeLoopSelfEnergyFunctions@"  -> threeLoopSelfEnergyFunctions,
                            "@threeLoopHiggsHeaders@"         -> threeLoopHiggsHeaders,
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
                            "@getPhysical@"                  -> IndentText[getPhysical],
                            "@setPhysical@"                  -> IndentText[setPhysical],
                            "@getMasses@"                    -> IndentText[getMasses],
                            "@setMasses@"                    -> IndentText[setMasses],
                            "@printMixingMatrices@"          -> IndentText[printMixingMatrices],
                            "@dependenceNumPrototypes@"      -> IndentText[dependenceNumPrototypes],
                            "@dependenceNumFunctions@"       -> WrapLines[dependenceNumFunctions],
                            "@solveEwsbTreeLevel@"           -> IndentText[WrapLines[solveEwsbTreeLevel]],
                            "@saveEwsbOutputParameters@"     -> IndentText[saveEwsbOutputParameters],
                            "@restoreEwsbOutputParameters@"  -> IndentText[restoreEwsbOutputParameters],
                            "@saveSoftHiggsMasses@"          -> IndentText[saveSoftHiggsMasses],
                            "@restoreSoftHiggsMasses@"       -> IndentText[restoreSoftHiggsMasses],
                            "@solveTreeLevelEWSBviaSoftHiggsMasses@" -> IndentText[WrapLines[solveTreeLevelEWSBviaSoftHiggsMasses]],
                            "@solveEWSBTemporarily@"         -> IndentText[solveEWSBTemporarily],
                            "@EWSBSolvers@"                  -> IndentText[IndentText[EWSBSolvers]],
                            "@fillArrayWithEWSBParameters@"  -> IndentText[IndentText[fillArrayWithEWSBParameters]],
                            "@solveEwsbWithTadpoles@"        -> IndentText[WrapLines[solveEwsbWithTadpoles]],
                            "@getEWSBParametersFromGSLVector@" -> IndentText[getEWSBParametersFromGSLVector],
                            "@setEWSBParametersFromLocalCopies@" -> IndentText[setEWSBParametersFromLocalCopies],
                            "@setEWSBParametersFromGSLVector@"   -> IndentText[setEWSBParametersFromGSLVector],
                            "@ewsbParametersInitializationList@" -> ewsbParametersInitializationList,
                            "@setEWSBSolution@"              -> IndentText[setEWSBSolution],
                            "@convertMixingsToSLHAConvention@" -> IndentText[convertMixingsToSLHAConvention],
                            "@convertMixingsToHKConvention@"   -> IndentText[convertMixingsToHKConvention],
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

(* Write the observables files *)
WriteObservables[extraSLHAOutputBlocks_, files_List] :=
    Module[{calculateObservables},
           calculateObservables = Observables`CalculateObservables[extraSLHAOutputBlocks, "observables"];
           WriteOut`ReplaceInFiles[files,
                                   {   "@calculateObservables@" -> IndentText[calculateObservables],
                                       Sequence @@ GeneralReplacementRules[]
                                   } ];
           ];

WriteUserExample[inputParameters_List, files_List] :=
    Module[{parseCmdLineOptions, printCommandLineOptions, spectrumGen},
           parseCmdLineOptions = WriteOut`ParseCmdLineOptions[inputParameters];
           printCommandLineOptions = WriteOut`PrintCmdLineOptions[inputParameters];
           spectrumGen = If[SMTower, CConversion`ToValidCSymbolString[FlexibleSUSY`FSModelName] <> "_SM",
                              CConversion`ToValidCSymbolString[FlexibleSUSY`FSModelName],
                              CConversion`ToValidCSymbolString[FlexibleSUSY`FSModelName]
                           ];
           WriteOut`ReplaceInFiles[files,
                          { "@parseCmdLineOptions@" -> IndentText[IndentText[parseCmdLineOptions]],
                            "@printCommandLineOptions@" -> IndentText[IndentText[printCommandLineOptions]],
                            "@SpectrumGenerator@" -> spectrumGen,
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
           concatenatedFileList = "\t" <> Utils`StringJoinWithSeparator[rgeFile, " \\\n\t"];
           WriteOut`ReplaceInFiles[files,
                          { "@generatedBetaFunctionModules@" -> concatenatedFileList,
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteUtilitiesClass[massMatrices_List, betaFun_List, minpar_List, extpar_List,
                    inputParameters_List,
                    lesHouchesInputParameters_List, extraSLHAOutputBlocks_List,
                    files_List] :=
    Module[{k, particles, susyParticles, smParticles,
            fillSpectrumVectorWithSusyParticles = "",
            fillSpectrumVectorWithSMParticles = "",
            particleLaTeXNames = "",
            particleNames = "", particleEnum = "", particleMultiplicity = "",
            particleMixingEnum = "", particleMixingNames = "",
            parameterNames = "", parameterEnum = "", numberOfParameters = 0,
            inputParameterEnum = "", inputParameterNames = "",
            isLowEnergyModel = "false",
            isSupersymmetricModel = "false",
            fillInputParametersFromMINPAR = "", fillInputParametersFromEXTPAR = "",
            writeSLHAMassBlock = "", writeSLHAMixingMatricesBlocks = "",
            writeSLHAModelParametersBlocks = "", writeSLHAMinparBlock = "",
            writeSLHAExtparBlock = "", readLesHouchesInputParameters,
            writeExtraSLHAOutputBlock = "",
            readLesHouchesOutputParameters, readLesHouchesPhysicalParameters,
            gaugeCouplingNormalizationDecls = "",
            gaugeCouplingNormalizationDefs = "",
            numberOfDRbarBlocks, drBarBlockNames},
           particles = GetMassEigenstate /@ massMatrices;
           susyParticles = Select[particles, (!SARAH`SMQ[#])&];
           smParticles   = Complement[particles, susyParticles];
           particleEnum       = TreeMasses`CreateParticleEnum[particles];
           particleMixingEnum = TreeMasses`CreateParticleMixingEnum[massMatrices];
           particleMultiplicity = TreeMasses`CreateParticleMultiplicity[particles];
           particleNames      = TreeMasses`CreateParticleNames[particles];
           particleLaTeXNames = TreeMasses`CreateParticleLaTeXNames[particles];
           particleMixingNames= TreeMasses`CreateParticleMixingNames[massMatrices];
           inputParameterEnum  = Parameters`CreateInputParameterEnum[inputParameters];
           inputParameterNames = Parameters`CreateInputParameterNames[inputParameters];
           fillSpectrumVectorWithSusyParticles = TreeMasses`FillSpectrumVector[susyParticles];
           fillSpectrumVectorWithSMParticles   = TreeMasses`FillSpectrumVector[smParticles];
           numberOfParameters = BetaFunction`CountNumberOfParameters[betaFun];
           parameterEnum      = BetaFunction`CreateParameterEnum[betaFun];
           parameterNames     = BetaFunction`CreateParameterNames[betaFun];
           isLowEnergyModel = If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY === True, "true", "false"];
           isSupersymmetricModel = If[SARAH`SupersymmetricModel === True, "true", "false"];
           fillInputParametersFromMINPAR = Parameters`FillInputParametersFromTuples[minpar, "MINPAR"];
           fillInputParametersFromEXTPAR = Parameters`FillInputParametersFromTuples[extpar, "EXTPAR"];
           readLesHouchesInputParameters = WriteOut`ReadLesHouchesInputParameters[lesHouchesInputParameters];
           readLesHouchesOutputParameters = WriteOut`ReadLesHouchesOutputParameters[];
           readLesHouchesPhysicalParameters = WriteOut`ReadLesHouchesPhysicalParameters["LOCALPHYSICAL", "DEFINE_PHYSICAL_PARAMETER"];
           writeSLHAMassBlock = WriteOut`WriteSLHAMassBlock[massMatrices];
           writeSLHAMixingMatricesBlocks  = WriteOut`WriteSLHAMixingMatricesBlocks[];
           writeSLHAModelParametersBlocks = WriteOut`WriteSLHAModelParametersBlocks[];
           writeSLHAMinparBlock = WriteOut`WriteSLHAMinparBlock[minpar];
           writeSLHAExtparBlock = WriteOut`WriteSLHAExtparBlock[extpar];
           writeExtraSLHAOutputBlock = WriteOut`WriteExtraSLHAOutputBlock[extraSLHAOutputBlocks];
           numberOfDRbarBlocks  = WriteOut`GetNumberOfDRbarBlocks[];
           drBarBlockNames      = WriteOut`GetDRbarBlockNames[];
           gaugeCouplingNormalizationDecls = WriteOut`GetGaugeCouplingNormalizationsDecls[SARAH`Gauge];
           gaugeCouplingNormalizationDefs  = WriteOut`GetGaugeCouplingNormalizationsDefs[SARAH`Gauge];
           WriteOut`ReplaceInFiles[files,
                          { "@fillSpectrumVectorWithSusyParticles@" -> IndentText[fillSpectrumVectorWithSusyParticles],
                            "@fillSpectrumVectorWithSMParticles@"   -> IndentText[IndentText[fillSpectrumVectorWithSMParticles]],
                            "@particleEnum@"       -> IndentText[WrapLines[particleEnum]],
                            "@particleMixingEnum@" -> IndentText[WrapLines[particleMixingEnum]],
                            "@particleMultiplicity@" -> IndentText[WrapLines[particleMultiplicity]],
                            "@particleNames@"      -> IndentText[WrapLines[particleNames]],
                            "@particleLaTeXNames@" -> IndentText[WrapLines[particleLaTeXNames]],
                            "@parameterEnum@"     -> IndentText[WrapLines[parameterEnum]],
                            "@parameterNames@"     -> IndentText[WrapLines[parameterNames]],
                            "@particleMixingNames@"-> IndentText[WrapLines[particleMixingNames]],
                            "@inputParameterEnum@" -> IndentText[WrapLines[inputParameterEnum]],
                            "@inputParameterNames@"-> IndentText[WrapLines[inputParameterNames]],
                            "@isLowEnergyModel@"   -> isLowEnergyModel,
                            "@isSupersymmetricModel@" -> isSupersymmetricModel,
                            "@fillInputParametersFromMINPAR@" -> IndentText[fillInputParametersFromMINPAR],
                            "@fillInputParametersFromEXTPAR@" -> IndentText[fillInputParametersFromEXTPAR],
                            "@readLesHouchesInputParameters@" -> IndentText[readLesHouchesInputParameters],
                            "@readLesHouchesOutputParameters@" -> IndentText[readLesHouchesOutputParameters],
                            "@readLesHouchesPhysicalParameters@" -> IndentText[readLesHouchesPhysicalParameters],
                            "@writeSLHAMassBlock@" -> IndentText[writeSLHAMassBlock],
                            "@writeSLHAMixingMatricesBlocks@"  -> IndentText[writeSLHAMixingMatricesBlocks],
                            "@writeSLHAModelParametersBlocks@" -> IndentText[writeSLHAModelParametersBlocks],
                            "@writeSLHAMinparBlock@"           -> IndentText[writeSLHAMinparBlock],
                            "@writeSLHAExtparBlock@"           -> IndentText[writeSLHAExtparBlock],
                            "@writeExtraSLHAOutputBlock@"      -> IndentText[writeExtraSLHAOutputBlock],
                            "@gaugeCouplingNormalizationDecls@"-> IndentText[gaugeCouplingNormalizationDecls],
                            "@gaugeCouplingNormalizationDefs@" -> IndentText[gaugeCouplingNormalizationDefs],
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
           fileNames = { "BetaYijk.m", "BetaGauge.m", "BetaMuij.m",
                         "BetaTijk.m", "BetaBij.m", "BetaVEV.m" };
           If[SARAH`AddDiracGauginos === True,
              AppendTo[fileNames, "BetaDGi.m"];
             ];
           If[SARAH`SupersymmetricModel === False,
              AppendTo[fileNames, "BetaLijkl.m"];
             ];
           If[SARAH`SupersymmetricModel === True,
              fileNames = Join[fileNames,
                               { "BetaWijkl.m", "BetaQijkl.m", "BetaLSi.m",
                                 "BetaLi.m", "Betam2ij.m", "BetaMi.m" }];
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
    Module[{particles, particle, i, precisionList = {}, higgs},
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
           higgs = Cases[precisionList, {SARAH`HiggsBoson | SARAH`PseudoScalar | SARAH`ChargedHiggs, LowPrecision}];
           Message[ReadPoleMassPrecisions::ImpreciseHiggs, #[[1]], #[[2]]]& /@ higgs;
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

FindUnfixedParameters[parameters_List, fixed_List] :=
    Module[{fixedParameters, autoInput},
           autoInput = { SARAH`hyperchargeCoupling, SARAH`leftCoupling,
                                            SARAH`strongCoupling };
           If[FlexibleSUSY`SMTower,
           autoInput = Join[autoInput, {SARAH`UpYukawa, SARAH`DownYukawa, SARAH`ElectronYukawa},
                            #[[1]] & /@ SUSYScaleUserMatching
                           ];
           ];
           fixedParameters = DeleteDuplicates[Flatten[Join[fixed,
                                          autoInput]]];
           Complement[parameters, fixedParameters]
          ];

GuessInputParameterType[FlexibleSUSY`Sign[par_]] :=
    CConversion`ScalarType[CConversion`integerScalarCType];
GuessInputParameterType[FlexibleSUSY`Phase[par_]] :=
    CConversion`ScalarType[CConversion`complexScalarCType];
GuessInputParameterType[par_] :=
    CConversion`ScalarType[CConversion`realScalarCType];

(* returns beta functions of VEV phases *)
GetVEVPhases[eigenstates_:FlexibleSUSY`FSEigenstates] :=
    Flatten @ Cases[DEFINITION[eigenstates][SARAH`VEVs], {_,_,_,_, p_} :> p];

AddSM3LoopRGE[beta_List, couplings_List] :=
    Module[{rules, MakeRule},
           MakeRule[coupling_] := {
               RuleDelayed[{coupling         , b1_, b2_}, {coupling       , b1, b2, Last[ThreeLoopSM`BetaSM[coupling]]}],
               RuleDelayed[{coupling[i1_,i2_], b1_, b2_}, {coupling[i1,i2], b1, b2, Last[ThreeLoopSM`BetaSM[coupling]] CConversion`PROJECTOR}]
           };
           rules = Flatten[MakeRule /@ couplings];
           beta /. rules
          ];

AddSM3LoopRGEs[] := Module[{
    gauge = { SARAH`hyperchargeCoupling,
              SARAH`leftCoupling,
              SARAH`strongCoupling },
    yuks  = { SARAH`UpYukawa,
              SARAH`DownYukawa,
              SARAH`ElectronYukawa },
    quart = { Parameters`GetParameterFromDescription["SM Higgs Selfcouplings"] },
    bilin = { Parameters`GetParameterFromDescription["SM Mu Parameter"] }
    },
    SARAH`BetaGauge = AddSM3LoopRGE[SARAH`BetaGauge, gauge];
    SARAH`BetaYijk  = AddSM3LoopRGE[SARAH`BetaYijk , yuks];
    SARAH`BetaLijkl = AddSM3LoopRGE[SARAH`BetaLijkl, quart];
    SARAH`BetaBij   = AddSM3LoopRGE[SARAH`BetaBij  , bilin];
    ];

AddMSSM3LoopRGE[beta_List, couplings_List] :=
    Module[{rules, MakeRule},
           MakeRule[coupling_] := {
               RuleDelayed[{coupling         , b1_, b2_}, {coupling       , b1, b2, Last[ThreeLoopMSSM`BetaMSSM[coupling]]}],
               RuleDelayed[{coupling[i1_,i2_], b1_, b2_}, {coupling[i1,i2], b1, b2, Last[ThreeLoopMSSM`BetaMSSM[coupling]]}]
           };
           rules = Flatten[MakeRule /@ couplings];
           beta /. rules
          ];

AddMSSM3LoopRGEs[] := Module[{
    gauge = { SARAH`hyperchargeCoupling,
              SARAH`leftCoupling,
              SARAH`strongCoupling },
    yuks  = { SARAH`UpYukawa,
              SARAH`DownYukawa,
              SARAH`ElectronYukawa },
    gaugi = { Parameters`GetParameterFromDescription["Bino Mass parameter"],
              Parameters`GetParameterFromDescription["Wino Mass parameter"],
              Parameters`GetParameterFromDescription["Gluino Mass parameter"] },
    trili = { SARAH`TrilinearUp, SARAH`TrilinearDown, SARAH`TrilinearLepton },
    mass2 = { SARAH`SoftSquark, SARAH`SoftUp, SARAH`SoftDown,
              SARAH`SoftLeftLepton, SARAH`SoftRightLepton,
              Parameters`GetParameterFromDescription["Softbreaking Down-Higgs Mass"],
              Parameters`GetParameterFromDescription["Softbreaking Up-Higgs Mass"] },
    mu    = { Parameters`GetParameterFromDescription["Mu-parameter"] },
    bmu   = { Parameters`GetParameterFromDescription["Bmu-parameter"] }
    },
    SARAH`BetaGauge = AddMSSM3LoopRGE[SARAH`BetaGauge, gauge];
    SARAH`BetaYijk  = AddMSSM3LoopRGE[SARAH`BetaYijk , yuks];
    SARAH`BetaMi    = AddMSSM3LoopRGE[SARAH`BetaMi   , gaugi];
    SARAH`BetaMuij  = AddMSSM3LoopRGE[SARAH`BetaMuij , mu   ];
    SARAH`BetaBij   = AddMSSM3LoopRGE[SARAH`BetaBij  , bmu  ];
    SARAH`BetaTijk  = AddMSSM3LoopRGE[SARAH`BetaTijk , trili];
    SARAH`Betam2ij  = AddMSSM3LoopRGE[SARAH`Betam2ij , mass2];
    ];

SelectRenormalizationScheme::UnknownRenormalizationScheme = "Unknown\
 renormalization scheme `1`.";

SelectRenormalizationScheme[renormalizationScheme_] :=
    Switch[renormalizationScheme,
           FlexibleSUSY`DRbar, 0,
           FlexibleSUSY`MSbar, 1,
           _, Message[SelectRenormalizationScheme::UnknownRenormalizationScheme, renormalizationScheme];
              Quit[1];
          ];

Options[MakeFlexibleSUSY] :=
    {
        InputFile -> "FlexibleSUSY.m",
        OutputDirectory -> "",
        DebugOutput -> False
    };

MakeFlexibleSUSY[OptionsPattern[]] :=
    Module[{nPointFunctions, runInputFile, initialGuesserInputFile,
            susyBetaFunctions, susyBreakingBetaFunctions,
            numberOfSusyParameters, anomDim,
            inputParameters (* list of 2-component lists of the form {name, type} *),
            haveEWSB = True,
            ewsbEquations, independentEwsbEquations,
            massMatrices, phases,
            diagonalizationPrecision,
            allInputParameterIndexReplacementRules = {},
            allParticles, allParameters,
            freePhases = {}, ewsbSolution = {},
            fixedParameters,
            treeLevelEwsbSolutionOutputFile, treeLevelEwsbEqsOutputFile,
            lesHouchesInputParameters, lesHouchesInputParameterReplacementRules,
            extraSLHAOutputBlocks,
	    vertexRules, vertexRuleFileName,
	    Lat$massMatrices},
           (* check if SARAH`Start[] was called *)
           If[!ValueQ[Model`Name],
              Print["Error: Model`Name is not defined.  Did you call SARAH`Start[\"Model\"]?"];
              Quit[1];
             ];
           FSDebugOutput = OptionValue[DebugOutput];
           FSOutputDir = OptionValue[OutputDirectory];
           If[!DirectoryQ[FSOutputDir],
              Print["Error: OutputDirectory ", FSOutputDir, " does not exist."];
              Print["   Please run ./createmodel first."];
              Quit[1]];
           CheckSARAHVersion[];
           (* load model file *)
           LoadModelFile[OptionValue[InputFile]];
           Print["FlexibleSUSY model file loaded"];
           Print["  Model: ", FlexibleSUSY`FSModelName];
           Print["  Model file: ", OptionValue[InputFile]];
           Print["  Model output directory: ", FSOutputDir];

           PrintHeadline["Reading SARAH output files"];
           (* get RGEs *)
           FSPrepareRGEs[];
           FSCheckLoopCorrections[FSEigenstates];
           nPointFunctions = EnforceCpColorStructures @ StripInvalidFieldIndices @
	      Join[PrepareSelfEnergies[FSEigenstates], PrepareTadpoles[FSEigenstates]];
           PrepareUnrotatedParticles[FSEigenstates];

           FlexibleSUSY`FSRenormalizationScheme = If[SARAH`SupersymmetricModel,
                                                     FlexibleSUSY`DRbar, FlexibleSUSY`MSbar];

           (* adapt SARAH`Conj to our needs *)
           (* Clear[Conj]; *)
           SARAH`Conj[(B_)[b__]] = .;
           SARAH`Conj /: SARAH`Conj[SARAH`Conj[x_]] := x;
           RXi[_] = 1;
           SARAH`Xi = 1;
           SARAH`Xip = 1;
           SARAH`rMS = SelectRenormalizationScheme[FlexibleSUSY`FSRenormalizationScheme];

           inputParameters = DeleteDuplicates[{#, GuessInputParameterType[#]}& /@ ((#[[2]])& /@ Utils`ForceJoin[SARAH`MINPAR, SARAH`EXTPAR])];
           Parameters`SetInputParameters[(#[[1]])& /@ inputParameters];

           If[FlexibleSUSY`UseSM3LoopRGEs,
              Print["Adding SM 3-loop beta-functions from ",
                    "[arxiv:1303.4364v2, arXiv:1307.3536v4,",
                    " arXiv:1504.05200 (SUSYHD v1.0.1)]"];
              AddSM3LoopRGEs[];
             ];

           If[FlexibleSUSY`UseMSSM3LoopRGEs,
              Print["Adding MSSM 3-loop beta-functions from ",
                    "[arxiv:hep-ph/0308231, arxiv:hep-ph/0408128]"];
              AddMSSM3LoopRGEs[];
             ];

           If[SARAH`SupersymmetricModel,
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
              ,
              (* pick beta functions of dimensionless parameters *)
              susyBetaFunctions = { SARAH`BetaGauge,
                                    SARAH`BetaLijkl, (* quartic scalar interactions *)
                                    SARAH`BetaYijk };

              (* pick beta functions of dimensionfull parameters *)
              susyBreakingBetaFunctions = { SARAH`BetaTijk, (* cubic scalar interactions *)
                                            SARAH`BetaMuij, (* bilinear fermion term *)
                                            SARAH`BetaBij , (* bilinear scalar term *)
                                            SARAH`BetaLi  , (* linear scalar term *)
                                            SARAH`BetaVEV };
             ];

           (* store all model parameters *)
           allParameters = ((#[[1]])& /@ Join[Join @@ susyBetaFunctions, Join @@ susyBreakingBetaFunctions]) /.
                               a_[Susyno`LieGroups`i1] :> a /.
                               a_[Susyno`LieGroups`i1,SARAH`i2] :> a;
           allIndexReplacementRules = Parameters`CreateIndexReplacementRules[allParameters];
           Parameters`SetModelParameters[allParameters];

           (* collect all phases from SARAH *)
           phases = DeleteDuplicates @ Join[
               ConvertSarahPhases[SARAH`ParticlePhases],
               Exp[I #]& /@ GetVEVPhases[FlexibleSUSY`FSEigenstates]];
           Parameters`SetPhases[phases];

           susyBetaFunctions = BetaFunction`ConvertSarahRGEs[susyBetaFunctions];
           susyBetaFunctions = Select[susyBetaFunctions, (BetaFunction`GetAllBetaFunctions[#]!={})&];

           numberOfSusyParameters = BetaFunction`CountNumberOfParameters[susyBetaFunctions];
           anomDim = AnomalousDimension`ConvertSarahAnomDim[SARAH`Gij];

           susyBreakingBetaFunctions = ConvertSarahRGEs[susyBreakingBetaFunctions];
           susyBreakingBetaFunctions = Select[susyBreakingBetaFunctions, (BetaFunction`GetAllBetaFunctions[#]!={})&];

           If[Head[SARAH`RealParameters] === List,
              Parameters`AddRealParameter[SARAH`RealParameters];
             ];

           allBetaFunctions = Join[susyBetaFunctions, susyBreakingBetaFunctions];

           FlexibleSUSY`FSLesHouchesList = SA`LHList;

           (* search for unfixed parameters *)
           Constraint`CheckConstraint[FlexibleSUSY`LowScaleInput, "LowScaleInput"];
           Constraint`CheckConstraint[FlexibleSUSY`SUSYScaleInput, "SUSYScaleInput"];
           Constraint`CheckConstraint[FlexibleSUSY`HighScaleInput, "HighScaleInput"];
           Constraint`CheckConstraint[FlexibleSUSY`InitialGuessAtLowScale, "InitialGuessAtLowScale"];
           Constraint`CheckConstraint[FlexibleSUSY`InitialGuessAtSUSYScale, "InitialGuessAtSUSYScale"];
           Constraint`CheckConstraint[FlexibleSUSY`InitialGuessAtHighScale, "InitialGuessAtHighScale"];
           Constraint`SanityCheck[Join[ If[SMTower, FlexibleSUSY`InitialGuessAtSUSYScale,
                                                                 FlexibleSUSY`InitialGuessAtLowScale,
                                                                 FlexibleSUSY`InitialGuessAtLowScale],
                                       FlexibleSUSY`InitialGuessAtHighScale],
                                  "initial guess"
                                 ];

           (* add EWSB constraint to SUSY-scale constraint if not set *)
           If[FlexibleSUSY`SMTower,
            If[FreeQ[Join[FlexibleSUSY`SUSYScaleInput, FlexibleSUSY`HighScaleInput],
                     FlexibleSUSY`FSSolveEWSBFor[___]],
               AppendTo[FlexibleSUSY`SUSYScaleInput,
                        FlexibleSUSY`FSSolveEWSBFor[FlexibleSUSY`EWSBOutputParameters]
                        ];
               ];
            fixedParameters = Join[ Constraint`FindFixedParametersFromConstraint[FlexibleSUSY`SUSYScaleInput],
                                    Constraint`FindFixedParametersFromConstraint[FlexibleSUSY`HighScaleInput]
                                    ];,
            If[FreeQ[Join[FlexibleSUSY`LowScaleInput, FlexibleSUSY`SUSYScaleInput, FlexibleSUSY`HighScaleInput],
                     FlexibleSUSY`FSSolveEWSBFor[___]],
               AppendTo[FlexibleSUSY`SUSYScaleInput,
                        FlexibleSUSY`FSSolveEWSBFor[FlexibleSUSY`EWSBOutputParameters]
                        ];
               ];
            fixedParameters = Join[Constraint`FindFixedParametersFromConstraint[FlexibleSUSY`LowScaleInput],
                                    Constraint`FindFixedParametersFromConstraint[FlexibleSUSY`SUSYScaleInput],
                                    Constraint`FindFixedParametersFromConstraint[FlexibleSUSY`HighScaleInput]
                                    ];
             ];
           FlexibleSUSY`FSUnfixedParameters = FindUnfixedParameters[allParameters, fixedParameters];
           If[FlexibleSUSY`FSUnfixedParameters =!= {} &&
              FlexibleSUSY`AutomaticInputAtMSUSY =!= True,
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
           If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY === True &&
              FlexibleSUSY`AutomaticInputAtMSUSY,
              FlexibleSUSY`SUSYScaleInput = Join[FlexibleSUSY`SUSYScaleInput,
                                                 {#[[1]],#[[2]]}& /@ FlexibleSUSY`FSUnfixedParameters];
              inputParameters = DeleteDuplicates @ Join[inputParameters,
                                                        {#[[2]], #[[3]]}& /@ FlexibleSUSY`FSUnfixedParameters];
              Parameters`AddInputParameters[(#[[1]])& /@ inputParameters];
             ];

           lesHouchesInputParameters = DeleteDuplicates[
               Flatten[
                   Cases[
                       Join[FlexibleSUSY`LowScaleInput,
                            FlexibleSUSY`SUSYScaleInput,
                            FlexibleSUSY`HighScaleInput,
                            FlexibleSUSY`InitialGuessAtLowScale,
                            FlexibleSUSY`InitialGuessAtHighScale,
                            {FlexibleSUSY`LowScaleFirstGuess,
                             FlexibleSUSY`SUSYScaleFirstGuess,
                             FlexibleSUSY`HighScaleFirstGuess}
                           ],
                       SARAH`LHInput[p_] :> Parameters`StripIndices[p],
                       Infinity
                        ]
                      ]
           ];

           lesHouchesInputParameters = Select[{BetaFunction`GetName[#],
                                               Symbol[ToValidCSymbolString[BetaFunction`GetName[#]] <> "Input"],
                                               Parameters`GetRealTypeFromDimension @ SARAH`getDimParameters @ Parameters`StripIndices @ BetaFunction`GetName[#]}& /@
                                                  Join[susyBetaFunctions, susyBreakingBetaFunctions] /.
                                              a_[Susyno`LieGroups`i1] :> a /.
                                              a_[Susyno`LieGroups`i1,SARAH`i2] :> a,
                                              MemberQ[lesHouchesInputParameters,#[[1]]]&];

           (* determine type of extra input parameters *)
           FlexibleSUSY`FSExtraInputParameters = {#[[1]], #[[2]], Parameters`GetRealTypeFromDimension[#[[3]]]}& /@ FlexibleSUSY`FSExtraInputParameters;

           inputParameters = DeleteDuplicates @ Join[inputParameters,
                                                     {#[[1]], #[[3]]}& /@ FlexibleSUSY`FSExtraInputParameters,
                                                     {#[[2]], #[[3]]}& /@ lesHouchesInputParameters];
           Parameters`AddInputParameters[(#[[1]])& /@ inputParameters];

           FlexibleSUSY`FSLesHouchesList = Join[FlexibleSUSY`FSLesHouchesList, {#[[1]], #[[2]]}& /@ FlexibleSUSY`FSExtraInputParameters];

           allInputParameterIndexReplacementRules = Parameters`CreateIndexReplacementRules[
               (* {parameter, type} *)
               {#[[1]], #[[3]]}& /@ FlexibleSUSY`FSExtraInputParameters
            ];

           (* replace all indices in the user-defined model file variables *)
           EvaluateUserInput[];
           ReplaceIndicesInUserInput[allIndexReplacementRules];
           ReplaceIndicesInUserInput[allInputParameterIndexReplacementRules];

           (* replace LHInput[p] by pInput in the constraints *)

           lesHouchesInputParameterReplacementRules = Flatten[{
               Rule[SARAH`LHInput[#[[1]]], #[[2]]],
               Rule[SARAH`LHInput[#[[1]][p__]], #[[2]][p]]
           }& /@ lesHouchesInputParameters];

           FlexibleSUSY`LowScaleInput = FlexibleSUSY`LowScaleInput /.
               lesHouchesInputParameterReplacementRules;
           FlexibleSUSY`SUSYScaleInput = FlexibleSUSY`SUSYScaleInput /.
               lesHouchesInputParameterReplacementRules;
           FlexibleSUSY`HighScaleInput = FlexibleSUSY`HighScaleInput /.
               lesHouchesInputParameterReplacementRules;

           FlexibleSUSY`InitialGuessAtLowScale = FlexibleSUSY`InitialGuessAtLowScale /.
               lesHouchesInputParameterReplacementRules;
           FlexibleSUSY`InitialGuessAtSUSYScale = FlexibleSUSY`InitialGuessAtSUSYScale /.
               lesHouchesInputParameterReplacementRules;
           FlexibleSUSY`InitialGuessAtHighScale = FlexibleSUSY`InitialGuessAtHighScale /.
               lesHouchesInputParameterReplacementRules;

           FlexibleSUSY`LowScaleFirstGuess = FlexibleSUSY`LowScaleFirstGuess /.
               lesHouchesInputParameterReplacementRules;
           FlexibleSUSY`SUSYScaleFirstGuess = FlexibleSUSY`SUSYScaleFirstGuess /.
               lesHouchesInputParameterReplacementRules;
           FlexibleSUSY`HighScaleFirstGuess = FlexibleSUSY`HighScaleFirstGuess /.
               lesHouchesInputParameterReplacementRules;

           If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY === True,
              lesHouchesInputParameters = Join[FlexibleSUSY`FSUnfixedParameters,
                                               lesHouchesInputParameters];
             ];

           numberOfSusyBreakingParameters = BetaFunction`CountNumberOfParameters[susyBreakingBetaFunctions];
           numberOfModelParameters = numberOfSusyParameters + numberOfSusyBreakingParameters;

           PrintHeadline["Creating model parameter classes"];
           Print["Creating class for susy parameters ..."];
           WriteRGEClass[susyBetaFunctions, anomDim,
                         {{FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_susy_parameters.hpp.in"}],
                           FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_susy_parameters.hpp"}]},
                          {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_susy_parameters.cpp.in"}],
                           FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_susy_parameters.cpp"}]}},
                         "two_scale_susy_beta_.cpp.in",
                         {{FileNameJoin[{$flexiblesusyTemplateDir, "two_scale.mk.in"}],
                           FileNameJoin[{FSOutputDir, "two_scale_susy.mk"}]}}
                        ];

           Print["Creating class for soft parameters ..."];
           WriteRGEClass[susyBreakingBetaFunctions, {},
                         {{FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_soft_parameters.hpp.in"}],
                           FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_soft_parameters.hpp"}]},
                          {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_soft_parameters.cpp.in"}],
                           FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_soft_parameters.cpp"}]}},
                         "two_scale_soft_beta_.cpp.in",
                         {{FileNameJoin[{$flexiblesusyTemplateDir, "two_scale.mk.in"}],
                           FileNameJoin[{FSOutputDir, "two_scale_soft.mk"}]}},
                         If[Head[SARAH`TraceAbbr] === List, SARAH`TraceAbbr, {}],
                         numberOfSusyParameters];

           (********************* EWSB *********************)
           ewsbEquations = SARAH`TadpoleEquations[FSEigenstates] /.
                           Parameters`ApplyGUTNormalization[] /.
                           allIndexReplacementRules /.
                           SARAH`sum[idx_, start_, stop_, expr_] :> Sum[expr, {idx,start,stop}];
           If[Head[ewsbEquations] =!= List,
              Print["Error: Could not find EWSB equations for eigenstates ",
                    FSEigenstates];
              Quit[1];
             ];

           (* filter out trivial EWSB eqs. *)
           ewsbEquations = Select[ewsbEquations, (#=!=0)&];

           haveEWSB = ewsbEquations =!= {};

           If[haveEWSB,
              ewsbEquations = Parameters`ExpandExpressions[ewsbEquations];
              FlexibleSUSY`EWSBOutputParameters = Parameters`DecreaseIndexLiterals[FlexibleSUSY`EWSBOutputParameters];

              (* adding tadpoles to the EWSB eqs. *)
              ewsbEquations = MapIndexed[#1 - tadpole[First[#2]]&, ewsbEquations];
              treeLevelEwsbSolutionOutputFile = FileNameJoin[{FSOutputDir,
                                                              FlexibleSUSY`FSModelName <> "_EWSB_solution.m"}];
              treeLevelEwsbEqsOutputFile      = FileNameJoin[{FSOutputDir,
                                                              FlexibleSUSY`FSModelName <> "_EWSB_equations.m"}];
              Print["Writing EWSB equations to ", treeLevelEwsbEqsOutputFile];
              Put[ewsbEquations, treeLevelEwsbEqsOutputFile];
              Print["Searching for independent EWSB equations ..."];
              independentEwsbEquations = Parameters`FilterOutLinearDependentEqs[ewsbEquations, FlexibleSUSY`EWSBOutputParameters];

              If[FlexibleSUSY`TreeLevelEWSBSolution === {},
                 (* trying to find an analytic solution for the EWSB eqs. *)
                 Print["Solving ", Length[independentEwsbEquations],
                       " independent EWSB equations for ",
                       FlexibleSUSY`EWSBOutputParameters," ..."];
                 {ewsbSolution, freePhases} = EWSB`FindSolutionAndFreePhases[independentEwsbEquations,
                                                                             FlexibleSUSY`EWSBOutputParameters,
                                                                             treeLevelEwsbSolutionOutputFile];
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
                    Print["   The EWSB solution was written to the file:"];
                    Print["      ", treeLevelEwsbSolutionOutputFile];
                   ];
                 ,
                 If[Length[FlexibleSUSY`TreeLevelEWSBSolution] != Length[independentEwsbEquations],
                    Print["Error: not enough EWSB solutions given!"];
                    Print["   You provided solutions for ", Length[FlexibleSUSY`TreeLevelEWSBSolution],
                          " parameters."];
                    Print["   However, there are ", Length[independentEwsbEquations],
                          " independent EWSB eqs."];
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
              ,
              Print["Note: There are no EWSB equations."];
             ];
           If[freePhases =!= {},
              Print["Note: adding free phases: ", freePhases];
              inputParameters = DeleteDuplicates @ Join[inputParameters,
                                                        {#, GuessInputParameterType[#]}& /@ freePhases];
              Parameters`AddInputParameters[(#[[1]])& /@ inputParameters];
             ];

           (* Fixed-point iteration can only be used if an analytic EWSB solution exists *)
           If[ewsbSolution === {} && MemberQ[FlexibleSUSY`FSEWSBSolvers, FlexibleSUSY`FPIRelative],
              Print["Warning: FPIRelative was selected, but no analytic"];
              Print["   solution to the EWSB eqs. is provided."];
              Print["   FPIRelative will be removed from the list of EWSB solvers."];
              FlexibleSUSY`FSEWSBSolvers = Cases[FlexibleSUSY`FSEWSBSolvers, Except[FlexibleSUSY`FPIRelative]];
             ];
           If[ewsbSolution === {} && MemberQ[FlexibleSUSY`FSEWSBSolvers, FlexibleSUSY`FPIAbsolute],
              Print["Warning: FPIAbsolute was selected, but no analytic"];
              Print["   solution to the EWSB eqs. is provided."];
              Print["   FPIAbsolute will be removed from the list of EWSB solvers."];
              FlexibleSUSY`FSEWSBSolvers = Cases[FlexibleSUSY`FSEWSBSolvers, Except[FlexibleSUSY`FPIAbsolute]];
             ];

           Print["Input parameters: ", InputForm[Parameters`GetInputParameters[]]];

           Print["Creating class for input parameters ..."];
           WriteInputParameterClass[inputParameters,
                                    {{FileNameJoin[{$flexiblesusyTemplateDir, "input_parameters.hpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_input_parameters.hpp"}]},
                                     {FileNameJoin[{$flexiblesusyTemplateDir, "input_parameters.cpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_input_parameters.cpp"}]}
                                    }
                                   ];

           lesHouchesInputParameters = Join[lesHouchesInputParameters, FlexibleSUSY`FSExtraInputParameters];

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

           (* check weak mixing angle parameters *)
           FlexibleSUSY`FSWeakMixingAngleOptions =
               Utils`FSSetOption[FlexibleSUSY`FSWeakMixingAngleOptions,
                                 FlexibleSUSY`FSWeakMixingAngleInput ->
                                 CheckWeakMixingAngleInputRequirements[Utils`FSGetOption[
                                     FlexibleSUSY`FSWeakMixingAngleOptions,
                                     FlexibleSUSY`FSWeakMixingAngleInput]
                                 ]
               ];

           PrintHeadline["Creating utilities"];
           Print["Creating class for convergence tester ..."];
           WriteConvergenceTesterClass[FlexibleSUSY`FSConvergenceCheck,
               {{FileNameJoin[{$flexiblesusyTemplateDir, "convergence_tester.hpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_convergence_tester.hpp"}]},
                {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_convergence_tester.hpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_convergence_tester.hpp"}]},
                {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_convergence_tester.cpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_convergence_tester.cpp"}]}
               }
                                      ];

           Print["Creating utilities class ..."];
           WriteUtilitiesClass[massMatrices, Join[susyBetaFunctions, susyBreakingBetaFunctions],
                               MINPAR, EXTPAR, inputParameters,
                               lesHouchesInputParameters, extraSLHAOutputBlocks,
               {{FileNameJoin[{$flexiblesusyTemplateDir, "info.hpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_info.hpp"}]},
                {FileNameJoin[{$flexiblesusyTemplateDir, "info.cpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_info.cpp"}]},
                {FileNameJoin[{$flexiblesusyTemplateDir, "utilities.hpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_utilities.hpp"}]},
                {FileNameJoin[{$flexiblesusyTemplateDir, "utilities.cpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_utilities.cpp"}]},
                {FileNameJoin[{$flexiblesusyTemplateDir, "slha_io.hpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_slha_io.hpp"}]},
                {FileNameJoin[{$flexiblesusyTemplateDir, "slha_io.cpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_slha_io.cpp"}]}
               }
                              ];

           Print["Creating plot scripts ..."];
           WritePlotScripts[{{FileNameJoin[{$flexiblesusyTemplateDir, "plot_spectrum.gnuplot.in"}],
                              FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_plot_spectrum.gnuplot"}]},
                             {FileNameJoin[{$flexiblesusyTemplateDir, "plot_rgflow.gnuplot.in"}],
                              FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_plot_rgflow.gnuplot"}]}}
                           ];

           PrintHeadline["Creating constraints"];
           Print["Creating class for high-scale constraint ..."];
           WriteConstraintClass[FlexibleSUSY`HighScale,
                                FlexibleSUSY`HighScaleInput,
                                FlexibleSUSY`HighScaleFirstGuess,
                                {FlexibleSUSY`HighScaleMinimum, FlexibleSUSY`HighScaleMaximum},
                                {{FileNameJoin[{$flexiblesusyTemplateDir, "high_scale_constraint.hpp.in"}],
                                  FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_high_scale_constraint.hpp"}]},
                                 {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_high_scale_constraint.hpp.in"}],
                                  FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_high_scale_constraint.hpp"}]},
                                 {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_high_scale_constraint.cpp.in"}],
                                  FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_high_scale_constraint.cpp"}]}
                                }
                               ];

           Print["Creating class for susy-scale constraint ..."];
           WriteConstraintClass[FlexibleSUSY`SUSYScale,
                                FlexibleSUSY`SUSYScaleInput,
                                FlexibleSUSY`SUSYScaleFirstGuess,
                                {FlexibleSUSY`SUSYScaleMinimum, FlexibleSUSY`SUSYScaleMaximum},
                                {{FileNameJoin[{$flexiblesusyTemplateDir, "susy_scale_constraint.hpp.in"}],
                                  FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_susy_scale_constraint.hpp"}]},
                                 {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_susy_scale_constraint.hpp.in"}],
                                  FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_susy_scale_constraint.hpp"}]},
                                 {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_susy_scale_constraint.cpp.in"}],
                                  FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_susy_scale_constraint.cpp"}]}
                                }
                               ];

           Print["Creating class for low-scale constraint ..."];
           WriteConstraintClass[FlexibleSUSY`LowScale,
                                FlexibleSUSY`LowScaleInput,
                                FlexibleSUSY`LowScaleFirstGuess,
                                {FlexibleSUSY`LowScaleMinimum, FlexibleSUSY`LowScaleMaximum},
                                {{FileNameJoin[{$flexiblesusyTemplateDir, "low_scale_constraint.hpp.in"}],
                                  FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_low_scale_constraint.hpp"}]},
                                 {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_low_scale_constraint.hpp.in"}],
                                  FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_low_scale_constraint.hpp"}]},
                                 {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_low_scale_constraint.cpp.in"}],
                                  FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_low_scale_constraint.cpp"}]}
                                }
                               ];

           Print["Creating class for initial guesser ..."];
           initialGuesserInputFile = "initial_guesser";
           If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY,
              initialGuesserInputFile = "initial_guesser_low_scale_model";
             ];
           WriteInitialGuesserClass[FlexibleSUSY`InitialGuessAtLowScale,
                                    FlexibleSUSY`InitialGuessAtSUSYScale,
                                    FlexibleSUSY`InitialGuessAtHighScale,
                                    {{FileNameJoin[{$flexiblesusyTemplateDir, "initial_guesser.hpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_initial_guesser.hpp"}]},
                                     {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_" <> initialGuesserInputFile <> ".hpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_initial_guesser.hpp"}]},
                                     {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_" <> initialGuesserInputFile <> ".cpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_initial_guesser.cpp"}]},
                                     {FileNameJoin[{$flexiblesusyTemplateDir, "SM_two_scale_" <> initialGuesserInputFile <> ".hpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_SM_two_scale_initial_guesser.hpp"}]},
                                     {FileNameJoin[{$flexiblesusyTemplateDir, "SM_two_scale_" <> initialGuesserInputFile <> ".cpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_SM_two_scale_initial_guesser.cpp"}]}
                                    }
                                   ];

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

           PrintHeadline["Creating SLHA model"];
           Print["Creating class for SLHA model ..."];
           WriteModelSLHAClass[massMatrices,
                               {{FileNameJoin[{$flexiblesusyTemplateDir, "model_slha.hpp.in"}],
                                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_model_slha.hpp"}]},
                                {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_model_slha.hpp.in"}],
                                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_model_slha.hpp"}]},
                                {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_model_slha.cpp.in"}],
                                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_model_slha.cpp"}]}
                               }];

           PrintHeadline["Creating model"];
           Print["Creating class for model ..."];
           WriteModelClass[massMatrices, ewsbEquations,
                           FlexibleSUSY`EWSBOutputParameters, ewsbSolution, freePhases,
                           nPointFunctions, vertexRules, Parameters`GetPhases[],
                           {{FileNameJoin[{$flexiblesusyTemplateDir, "model.hpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_model.hpp"}]},
                            {FileNameJoin[{$flexiblesusyTemplateDir, "mass_eigenstates.hpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_mass_eigenstates.hpp"}]},
                            {FileNameJoin[{$flexiblesusyTemplateDir, "mass_eigenstates.cpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_mass_eigenstates.cpp"}]},
                            {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_model.hpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_model.hpp"}]},
                            {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_model.cpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_model.cpp"}]},
                            {FileNameJoin[{$flexiblesusyTemplateDir, "physical.hpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_physical.hpp"}]},
                            {FileNameJoin[{$flexiblesusyTemplateDir, "physical.cpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_physical.cpp"}]}
                           },
                           diagonalizationPrecision];

            Print["Creating matching class ..."];
            WriteMatchingClass[ {{FileNameJoin[{$flexiblesusyTemplateDir, "SM_matching.hpp.in"}],
                                  FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_SM_matching.hpp"}]},
                                 {FileNameJoin[{$flexiblesusyTemplateDir, "SM_matching.cpp.in"}],
                                  FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_SM_matching.cpp"}]},
                                 {FileNameJoin[{$flexiblesusyTemplateDir, "SM_two_scale_matching.hpp.in"}],
                                  FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_SM_two_scale_matching.hpp"}]}}
                              ];

           Print["Creating observables"];
           WriteObservables[extraSLHAOutputBlocks,
                            {{FileNameJoin[{$flexiblesusyTemplateDir, "observables.hpp.in"}],
                              FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_observables.hpp"}]},
                             {FileNameJoin[{$flexiblesusyTemplateDir, "observables.cpp.in"}],
                              FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_observables.cpp"}]}}];

           Print["Creating user example spectrum generator program ..."];
           spectrumGeneratorInputFile = "high_scale_spectrum_generator.hpp.in";
           If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY,
              spectrumGeneratorInputFile = "low_scale_spectrum_generator.hpp.in";];
           WriteUserExample[inputParameters,
                            {{FileNameJoin[{$flexiblesusyTemplateDir, spectrumGeneratorInputFile}],
                              FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_spectrum_generator.hpp"}]},
                             {FileNameJoin[{$flexiblesusyTemplateDir, "SM_" <> spectrumGeneratorInputFile}],
                              FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_SM_spectrum_generator.hpp"}]},
                             {FileNameJoin[{$flexiblesusyTemplateDir, "spectrum_generator_interface.hpp.in"}],
                              FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_spectrum_generator_interface.hpp"}]},
                             {FileNameJoin[{$flexiblesusyTemplateDir, "run.cpp.in"}],
                              FileNameJoin[{FSOutputDir, "run_" <> FlexibleSUSY`FSModelName <> ".cpp"}]},
                             {FileNameJoin[{$flexiblesusyTemplateDir, "run_cmd_line.cpp.in"}],
                              FileNameJoin[{FSOutputDir, "run_cmd_line_" <> FlexibleSUSY`FSModelName <> ".cpp"}]},
                             {FileNameJoin[{$flexiblesusyTemplateDir, "scan.cpp.in"}],
                              FileNameJoin[{FSOutputDir, "scan_" <> FlexibleSUSY`FSModelName <> ".cpp"}]}
                            }];

           PrintHeadline["FlexibleSUSY has finished"];
          ];

End[];

EndPackage[];

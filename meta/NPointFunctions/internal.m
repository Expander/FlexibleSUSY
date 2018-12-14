(* ::Package:: *)

BeginPackage["NPointFunctions`",{"FeynArts`","FormCalc`","Utils`"}];
Utils`AssertWithMessage[
	MemberQ[$Packages, "FeynArts`"],
  "NPointFunctions`: Unable to load FeynArts` package"];
Utils`AssertWithMessage[
	MemberQ[$Packages, "FormCalc`"],
  "NPointFunctions`: Unable to load FormCalc` package"];

SetFAFCPaths::usage="";

LoopLevel::usage="";
Regularize::usage="";
DimensionalReduction::usage="";
DimensionalRegularization::usage="";
ZeroExternalMomenta::usage="";
NPointFunctionFAFC::usage="";

GenericS::usage="";
GenericF::usage="";
GenericV::usage="";
GenericU::usage="";
GenericT::usage="";
GenericSum::usage="";
GenericIndex::usage="";
LorentzIndex::usage="";
LorenzDot::usage="";
lt::usage="";

subexpressionToFSRules::usage="";

Begin["Private`"];

feynArtsDir = "";
formCalcDir = "";
feynArtsModel = "";
particleNamesFile = "";
substitutionsFile = "";
particleNamespaceFile = "";

SetFAFCPaths[feynArtsDirS_String, formCalcDirS_String, feynArtsModelS_String,
           particleNamesFileS_String, substitutionsFileS_String,
           particleNamespaceFileS_String] :=
  (feynArtsDir = feynArtsDirS;
   formCalcDir = formCalcDirS;
   feynArtsModel = feynArtsModelS;
   particleNamesFile = particleNamesFileS;
   substitutionsFile = substitutionsFileS;
   particleNamespaceFile = particleNamespaceFileS;

   SetFSConventionRules[];)

NPointFunctionFAFC[inFields_List,outFields_List,
    OptionsPattern[{LoopLevel -> 1,
                    Regularize -> DimensionalReduction,
                    ZeroExternalMomenta -> False}]]:=
  Module[{loopLevel = OptionValue[LoopLevel],
          regularizationScheme = OptionValue[Regularize],
          zeroExternalMomenta = OptionValue[ZeroExternalMomenta],
          topologies,diagrams,amplitudes,genericInsertions,
          symmetryFactors,fsFields, fsInFields,fsOutFields,
          externalMomentumRules, nPointFunction},
    If[loopLevel =!= 1,
      Return["NPointFunctions`NPointFunction[]: Only LoopLevel 1 is supported"]];
    If[regularizationScheme =!= DimensionalReduction &&
       regularizationScheme =!= DimensionalRegularization,
       Return["NPointFunctions`NPointFunction[]: Unknown regularization scheme: " <>
              ToString[regularizationScheme]]];
    If[BooleanQ[zeroExternalMomenta] === False,
       Return["NPointFunctions`NPointFunction[]: ZeroExternalMomenta must \
be either True or False"]];

    If[DirectoryQ[formCalcDir] === False,
       CreateDirectory[formCalcDir]];
    SetDirectory[formCalcDir];

    topologies = FeynArts`CreateTopologies[loopLevel,
      Length[inFields] -> Length[outFields],
      ExcludeTopologies -> Internal];
      
    diagrams = FeynArts`InsertFields[topologies,
      inFields -> outFields,
      InsertionLevel -> Classes,
      Model -> feynArtsModel];
    amplitudes = FeynArts`CreateFeynAmp[diagrams];

    genericInsertions = Flatten[
      GenericInsertionsForDiagram /@ (List @@ diagrams), 1];

    fsInFields = (List @@ Head[amplitudes][[1,2,1,All,1]]) //.
      fieldNameToFSRules;
    fsOutFields = (List @@ Head[amplitudes][[1,2,2,All,1]]) //.
      fieldNameToFSRules;

    fsFields = Join[fsInFields,fsOutFields];

    externalMomentumRules = {
      If[zeroExternalMomenta,
         SARAH`Mom[i_Integer, lorentzIndex_] :> 0,
         SARAH`Mom[i_Integer, lorentzIndex_] :> SARAH`Mom[fsFields[[i]], lorentzIndex]]
    };
    
    nPointFunction = {{fsInFields, fsOutFields},
      CalculateAmplitudes[amplitudes, genericInsertions,
        regularizationScheme, zeroExternalMomenta] /. externalMomentumRules};
    
    ResetDirectory[];
    nPointFunction
  ]

StripParticleIndices[Times[-1,field_]] := Times[-1, StripParticleIndices[field]]
StripParticleIndices[genericType_[classIndex_, ___]] := genericType[classIndex]

FindGenericInsertions[insertions_List]:=
  Module[{toGenericIndexConventionRules, genericFields, genericInsertions},
    toGenericIndexConventionRules =
      Cases[insertions[[1]], Rule[FeynArts`Field[index_Integer],type_Symbol] :>
        Rule[FeynArts`Field[index], type[FeynArts`Index[Generic,index]]]];

    genericFields = toGenericIndexConventionRules[[All,1]];
    genericInsertions = Cases[#, (Rule[genericField_,classesField_] /;
        MemberQ[genericFields, genericField]) :>
      Rule[genericField, StripParticleIndices[classesField]]] &
      /@ insertions[[2]];

    (List @@ genericInsertions) /. toGenericIndexConventionRules
  ]

GenericInsertionsForDiagram[diagram_Rule]:=
  List @@ (FindGenericInsertions /@ (List @@@ diagram[[2]]))

CombinatorialFactorsForAmplitudeInsertions[amplitude_FeynAmp]:=
  Module[{combinatorialPosition},
    combinatorialPosition = Position[amplitude[[-1,1]], FeynArts`RelativeCF][[1,1]];
    (List @@ amplitude[[-1, 2, All, combinatorialPosition]]) /.
      {FeynArts`SumOver[__] -> 1,
       FeynArts`IndexDelta[__] -> 1} (* FIXME: Can we really remove the IndexDelta? *)
  ]

AtomHead[x_] := If[AtomQ[x], x, AtomHead[Head[x]]]

CalculateAmplitudes[classesAmplitudes_, genericInsertions_List,
    regularizationScheme_, zeroExternalMomenta_] :=
  Module[{genericAmplitudes,calculatedAmplitudes,
          abbreviations,subexpressions,combinatorialFactors,
          prefixRules, dimensionParameter, pairs, zeroedRules, pair,
          onShellFlag = zeroExternalMomenta},
    combinatorialFactors = CombinatorialFactorsForAmplitudeInsertions /@
      (List @@ classesAmplitudes);
    genericAmplitudes = FeynArts`PickLevel[Generic][classesAmplitudes];

    dimensionParameter = Switch[regularizationScheme,
      DimensionalReduction, 4,
      DimensionalRegularization, D];

    genericAmplitudes = If[zeroExternalMomenta,
      Evaluate @ OffShell[genericAmplitudes, Sequence @@ Table[k -> 0, {k,
        Total[Length /@ Head[genericAmplitudes][[1,2]]]}]],
      genericAmplitudes];

    calculatedAmplitudes =
      (FormCalc`CalcFeynAmp[Head[genericAmplitudes][#],
                            Dimension -> dimensionParameter,
                            OnShell -> onShellFlag,
                            Invariants -> False] & /@
        genericAmplitudes) //. FormCalc`GenericList[];

    calculatedAmplitudes = SumOverAllFieldIndices /@ (List @@ calculatedAmplitudes);

    pairs = If[zeroExternalMomenta,
      Cases[Select[FormCalc`Abbr[],
              StringMatchQ[ToString @ #[[1]], "Pair"~~__] &],
            Rule[_, Pattern[pair, FormCalc`Pair[FormCalc`k[_], FormCalc`k[_]]]] :> pair],
      {}];
    zeroedRules = (Rule[#, 0] & /@ pairs);

    abbreviations = Complement[FormCalc`Abbr[], pairs] //. FormCalc`GenericList[];
    subexpressions = FormCalc`Subexpr[] //. FormCalc`GenericList[];

    {abbreviations, zeroedRules} = RecursivelyZeroRules[abbreviations, zeroedRules];
    {subexpressions, zeroedRules} = RecursivelyZeroRules[subexpressions, zeroedRules];

    calculatedAmplitudes = calculatedAmplitudes /. zeroedRules;

    FCAmplitudesToFSConvention[
        {calculatedAmplitudes, genericInsertions, combinatorialFactors},
      abbreviations, subexpressions]
  ]

RecursivelyZeroRules[nonzeroRules_List, zeroRules_List] :=
  Module[{nextNonzero, nextZero},
    nextNonzero = Rule @@@ Transpose[
      {nonzeroRules[[All,1]], nonzeroRules[[All,2]] /. zeroRules}];

    If[nextNonzero === nonzeroRules,
       Return[{nonzeroRules, zeroRules}]];

    nextZero = Join[zeroRules, Cases[nextNonzero, HoldPattern[Rule[_, 0]]]];
    nextNonzero = Complement[nextNonzero, nextZero];

    RecursivelyZeroRules[nextNonzero, nextZero]
  ]

SumOverAllFieldIndices[genericAmplitude_] :=
  Module[{genericIndices, fieldType, genericRules},
    Utils`AssertWithMessage[Length[genericAmplitude] === 1,
      "NPointFunctions`SumOverAllFieldIndices[]: \
Length of generic amplitude != 1 not implemented (incompatible FormCalc change?)"];

    genericIndices = DeleteDuplicates[
      Cases[List @@ genericAmplitude,
            Pattern[fieldType,Alternatives[S,F,V,U,T]][
              FeynArts`Index[Generic,number_Integer]] :> {fieldType,number},
            Infinity]];
    GenericSum[genericAmplitude[[1]], genericIndices]
  ]

FCAmplitudesToFSConvention[amplitudes_, abbreviations_, subexpressions_] :=
  Module[{fsAmplitudes, fsAbbreviations, fsSubexpressions},
    fsAmplitudes = amplitudes //. amplitudeToFSRules;
    fsAbbreviations = abbreviations //. subexpressionToFSRules;
    fsSubexpressions = subexpressions //. subexpressionToFSRules;

    {fsAmplitudes, Join[fsAbbreviations,fsSubexpressions]}
  ]

subexpressionToFSRules = {};
SetFSConventionRules[] :=
  Module[{substitutionsHandle, substitutionsContents,
          substitutions, invertibleSubstitutions,
          invertedSubstitutions, indexRules, field,
          fieldsHandle, fieldContents, fieldNames, massRules,
          fieldNamespaces, couplingRules, generalFCRules, fieldType},
    substitutionsHandle = OpenRead[substitutionsFile,BinaryFormat -> True];
    substitutionsContents = ReadString[substitutionsHandle];
    Close[substitutionsHandle];

    substitutions = ReleaseHold[List @@@ Map[Hold,
      ToExpression[substitutionsContents,StandardForm,Hold],
      {2}]];

    invertibleSubstitutions = Cases[substitutions,
      Hold[ass_[f_[fsExpr_],fcExpr_]] /;
        (MemberQ[{UpSet,Set},ass] &&
           !(f === Conjugate &&
             ToString[AtomHead[fsExpr]] === ToString[AtomHead[fcExpr]] <> "C"))];

    invertedSubstitutions =
      InvertFCConventionSubstitution /@ invertibleSubstitutions;

    fieldsHandle = OpenRead[particleNamesFile,BinaryFormat -> True];
    fieldContents = ReadString[fieldsHandle];
    Close[fieldsHandle];

    fieldNames = Cases[StringSplit[fieldContents,"\n"],
                       line_String /; (Length[StringSplit[line,": "]] == 2) :> 
                         StringSplit[line,": "]];
    fieldNamespaces = Get[particleNamespaceFile];
    
    fieldNames = 
      {Cases[fieldNamespaces, {#[[1]], _}][[1,2]] <>
       #[[1]], #[[2]]} & /@ fieldNames;
       
    massRules = Join[
      (Symbol["Mass" <> #][indices__] :> 
       SARAH`Mass[ToExpression[#][{indices}]]) & /@ fieldNames[[All,1]],
      (Symbol["Mass" <> #] ->  
       SARAH`Mass[ToExpression[#]]) & /@ fieldNames[[All,1]],
      {FeynArts`Mass[field_,faSpec_] :> SARAH`Mass[field]}
    ];

    couplingRules = {
        FeynArts`G[_][0][fields__][1] :> SARAH`Cp[fields][1],
        FeynArts`G[_][0][fields__][1] :> SARAH`Cp[fields][1],
        FeynArts`G[_][0][fields__][
            FeynArts`NonCommutative[Global`ChiralityProjector[-1]]] :>
          SARAH`Cp[fields][SARAH`PL],
        FeynArts`G[_][0][fields__][
            FeynArts`NonCommutative[Global`ChiralityProjector[1]]] :>
          SARAH`Cp[fields][SARAH`PR],
        FeynArts`G[_][0][fields__][
            Global`MetricTensor[KI1[i1_Integer],KI1[i2_Integer]]] :>
          SARAH`Cp[fields][SARAH`g[LorentzIndex[{fields}[[i1]]],
                                   LorentzIndex[{fields}[[i2]]]]],
        FeynArts`G[_][0][fields__][
            FeynArts`Mom[i1_Integer] - FeynArts`Mom[i2_Integer]] :>
          SARAH`Cp[fields][SARAH`Mom[{fields}[[i1]]] - SARAH`Mom[{fields}[[i2]]]]

    };

    generalFCRules = {
        FormCalc`Finite -> 1,
        FormCalc`Den[a_,b_] :> 1/(a - b),
        FormCalc`Pair[a_,b_] :> SARAH`sum[lt, 1, 4,
          SARAH`g[lt, lt] * Append[a, lt] * Append[b, lt]],
        Pattern[fieldType,Alternatives[S,F,V,U,T]][
          FeynArts`Index[Generic,number_Integer]] :>  
            fieldType[GenericIndex[number]],
        FormCalc`k[i_Integer, index___] :> SARAH`Mom[i, index]
    };

    indexRules = {
      FeynArts`Index[name_,index_Integer] :> 
        Symbol[If[StringTake[ToString[name],1] === "I",
                  ToString[name],
                  "I" <> ToString[name]] <> ToString[index]]
    };

    fieldNameToFSRules = Join[
      Rule @@@ Transpose[
        {Append[#,{indices___}] & /@ (ToExpression /@ fieldNames[[All,2]]),
         Through[(ToExpression /@ fieldNames[[All,1]])[{indices}]]}],
      Rule @@@ Transpose[
        {ToExpression /@ fieldNames[[All,2]], ToExpression /@ fieldNames[[All,1]]}],
      {S -> GenericS, F -> GenericF, V -> GenericV, U -> GenericU, T -> GenericT},
      {
        Times[-1, field_GenericS | field_GenericV] :>
           Susyno`LieGroups`conj[field],
        Times[-1, field_GenericF | field_GenericU] :>
           SARAH`bar[field]
      },
      (Times[-1, Pattern[field, # | Blank[#]]] :> 
           CXXDiagrams`LorentzConjugate[field]) & /@
         (ToExpression /@ fieldNames[[All,1]]),
      indexRules
    ];

    subexpressionToFSRules = Join[
      invertedSubstitutions,
      massRules,
      fieldNameToFSRules,
      couplingRules,
      generalFCRules
    ];
    amplitudeToFSRules = Join[
      subexpressionToFSRules,
      sumOverRules,
      {FeynArts`IndexSum -> Sum}
    ];
  ]

InvertFCConventionSubstitution[substitution_Hold] :=
  Module[{function,fsSymbol,fcSymbol,pattern,arg},
    function = Head[substitution[[1,1]]];
    fsSymbol = substitution[[1,1,1]];
    fcSymbol = substitution[[1,2]];

    If[Length[fsSymbol] =!= 0,
       pattern = fsSymbol[[1]]; fsSymbol = Head[fsSymbol],
       pattern = Null];
    If[Length[fcSymbol] =!= 0,
       arg = fcSymbol[[1]]; fcSymbol = Head[fcSymbol],
       arg = Null];

    If[(pattern === Null && arg =!= Null) ||
       (pattern =!= Null && arg === Null),
       Return["No heuristic to invert " <> ToString[substitution]]];

    If[pattern =!= Null,
       If[Head[pattern] =!= Pattern,
          Return["No heuristic to invert " <> ToString[substitution]]];
       If[pattern[[1]] =!= arg,
          Return["No heuristic to invert " <> ToString[substitution]]]];

    If[function === Conjugate,
       function = Susyno`LieGroups`conj];

    If[pattern =!= Null,
       fcSymbol[Pattern[Evaluate[arg],___]] :> Evaluate[function[fsSymbol[arg]]],
       fcSymbol -> function[fsSymbol]]
  ]

sumOverRules = {
  FeynArts`SumOver[_,_,FeynArts`External] :> Sequence[],
  Times[expr_,FeynArts`SumOver[index_,max_Integer]]
    :> SARAH`sum[index,1,max,expr],
  Times[expr_,FeynArts`SumOver[index_,{min_Integer,max_Integer}]]
    :> SARAH`sum[index,min,max,expr],
  SARAH`sum[index_,min_Integer,max_Integer,FeynArts`SumOver[_,max2_Integer]]
    :> SARAH`sum[index,1,max,max2],
  SARAH`sum[index_,min_Integer,max_Integer,
      FeynArts`SumOver[_,{min2_Integer,max2_Integer}]]
    :> SARAH`sum[index,1,max,max2-min2]
};

End[];
EndPackage[];




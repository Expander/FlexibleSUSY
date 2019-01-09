(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

BeginPackage["NPointFunctions`",{"FlexibleSUSY`","SARAH`","CXXDiagrams`","Vertices`","Parameters`","Utils`"}];

LoopLevel::usage="";
Regularize::usage="";
DimensionalReduction::usage="";
DimensionalRegularization::usage="";
ZeroExternalMomenta::usage="";
LoopFunctions::usage="Option that controls whether to use FlexibleSUSY or LoopTools for loop functions."

NPointFunction::usage="";
VerticesForNPointFunction::usage="";
CreateCXXFunctions::usage="";
CreateCXXHeaders::usage="";

A0i::usage="";
B0i::usage="";
C0i::usage="";
D0i::usage="";

aa0::usage="";

bb0::usage="";
bb1::usage="";

cc0::usage="";
cc1::usage="";
cc2::usage="";

dd0::usage="";
dd1::usage="";
dd2::usage="";
dd3::usage="";

GenericS::usage="";
GenericF::usage="";
GenericV::usage="";
GenericU::usage="";
GenericT::usage="";
GenericSum::usage="";
GenericIndex::usage="";
LorentzIndex::usage="";

Begin["`Private`"];

SARAHModelName[] := If[SARAH`submodeldir =!= False,
                       SARAH`modelDir <> "-" <> SARAH`submodeldir,
                       SARAH`modelDir]

ClassesModelFileName[] := SARAH`ModelName <> ToString[FlexibleSUSY`FSEigenstates]
SubstitutionsFileName[] :=
  "Substitutions-" <> SARAH`ModelName <>
  ToString[FlexibleSUSY`FSEigenstates] <> ".m"

NPointFunction[inFields_List,outFields_List,
    OptionsPattern[{LoopLevel -> 1,
                    Regularize -> Switch[FlexibleSUSY`FSRenormalizationScheme,
                      FlexibleSUSY`DRbar, DimensionalReduction,
                      FlexibleSUSY`MSbar, DimensionalRegularization],
                    ZeroExternalMomenta -> False}]]:=
  Module[{loopLevel, zeroExternalMomenta, regularizationScheme, nPointMeta,
          sarahOutputDir = SARAH`$sarahCurrentOutputMainDir,
          fsMetaDir = $flexiblesusyMetaDir,
          outputDir,
          currentPath, currentDirectory,
          feynArtsDir,formCalcDir,nPointFunctionsDir,
          cachedNPointFunction,
          feynArtsModel,substitutionsFile,particleNamesFile,
          inFANames,outFANames,
          subKernel,calculationCommand,particleNamespaceFile,
          fileHandle,nPointFunction},
    loopLevel = OptionValue[LoopLevel];
    Utils`AssertWithMessage[loopLevel === 1,
			"NPointFunctions`NPointFunction[]: Only loop level 1 is supported"];
		
    regularizationScheme = OptionValue[Regularize];
    Utils`AssertWithMessage[
			regularizationScheme === DimensionalReduction ||
      regularizationScheme === DimensionalRegularization,
			"NPointFunctions`NPointFunction[]: Unknown regularization scheme: " <>
			ToString[regularizationScheme]];
		
    zeroExternalMomenta = OptionValue[ZeroExternalMomenta];
    Utils`AssertWithMessage[zeroExternalMomenta === True || zeroExternalMomenta === False,
			"NPointFunctions`NPointFunction[]: Option ZeroExternalMomenta must \
be either True or False"];

    Utils`AssertWithMessage[And @@
			TreeMasses`IsScalar /@ Join[inFields, outFields],
			"NPointFunctions`NPointFunction[]: Only external scalars are \
supported (for now)."];

    inFields = Vertices`StripFieldIndices[inFields];
    outFields = Vertices`StripFieldIndices[outFields];
    nPointMeta = {loopLevel, regularizationScheme, zeroExternalMomenta};

    outputDir = FileNameJoin[{sarahOutputDir, ToString[FlexibleSUSY`FSEigenstates]}];

    nPointFunctionsDir = FileNameJoin[{outputDir, "NPointFunctions"}];
    If[DirectoryQ[nPointFunctionsDir] == False,
       CreateDirectory[nPointFunctionsDir]];

    nPointFunctionFile = FileNameJoin[{nPointFunctionsDir, "temp"}];
    nPointFunction = CachedNPointFunction[
      inFields, outFields, nPointMeta];
    If[nPointFunction =!= Null,
       Return[nPointFunction]];

    feynArtsDir = FileNameJoin[{outputDir, "FeynArts"}];
    formCalcDir = FileNameJoin[{outputDir, "FormCalc"}];
    
    feynArtsModel = FileNameJoin[{feynArtsDir, ClassesModelFileName[]}];
    particleNamesFile = FileNameJoin[{feynArtsDir, "ParticleNamesFeynArts.dat"}];
    particleNamespaceFile = FileNameJoin[{feynArtsDir, "ParticleNamespaces.m"}];
    substitutionsFile = FileNameJoin[{feynArtsDir, SubstitutionsFileName[]}];
    
    subKernels = LaunchKernels[2];

    If[FileExistsQ[feynArtsModel <> ".mod"] === False,
       GenerateFAModelFileOnKernel[subKernels[[1]]];
       WriteParticleNamespaceFile[particleNamespaceFile]];

    inFANames = FeynArtsNamesForFields[inFields, particleNamesFile];
    outFANames = FeynArtsNamesForFields[outFields, particleNamesFile];
    
    currentPath = $Path;
    currentDirectory = Directory[];

    (* Unfortunately, there seems to be no way to restrict
      this to a specific kernel *)
    DistributeDefinitions[currentPath, currentDirectory,
      fsMetaDir, feynArtsDir, formCalcDir, feynArtsModel,
      particleNamesFile, substitutionsFile, particleNamespaceFile,
      inFANames, outFANames,
      loopLevel, regularizationScheme, zeroExternalMomenta];

    nPointFunction = ParallelEvaluate[
      $Path = currentPath;
      SetDirectory[currentDirectory];
      
      Get[FileNameJoin[{fsMetaDir, "NPointFunctions", "internal.m"}]];
      
      NPointFunctions`SetFAFCPaths[
        feynArtsDir, formCalcDir, feynArtsModel,
        particleNamesFile, substitutionsFile,
        particleNamespaceFile];
      
      NPointFunctions`NPointFunctionFAFC[
        ToExpression[inFANames], ToExpression[outFANames],
        LoopLevel -> loopLevel,
        Regularize -> regularizationScheme,
        ZeroExternalMomenta -> zeroExternalMomenta],
      subKernels[[2]]];
    
    CloseKernels[subKernels];
    
    Utils`AssertWithMessage[nPointFunction =!= $Failed,
			"NPointFunctions`NPointFunction[]: Calculation failed"];
    
    CacheNPointFunction[nPointFunction, nPointMeta];
    nPointFunction
  ]

GenerateFAModelFileOnKernel[kernel_] :=
  Module[{currentPath, currentDirectory,
          fsMetaDir = $flexiblesusyMetaDir,
          sarahInputDirectories, sarahOutputDirectory,
          sarahModelname, eigenstates},
    currentPath = $Path;
    currentDirectory = Directory[];
    sarahInputDirectories = SARAH`SARAH[InputDirectories];
    sarahOutputDirectory = SARAH`SARAH[OutputDirectory];
    sarahModelName = SARAHModelName[];
    eigenstates = FlexibleSUSY`FSEigenstates;
    
    (* Unfortunately, there seems to be no way to restrict
      this to a specific kernel *)
    DistributeDefinitions[currentPath, currentDirectory,
      fsMetaDir, sarahInputDirectories, sarahOutputDirectory,
      sarahModelName, eigenstates];
    
    ParallelEvaluate[
      $Path = currentPath;
      SetDirectory[currentDirectory];
      
      Get[FileNameJoin[{fsMetaDir, "NPointFunctions", "createFAModelFile.m"}]];
      
      NPointFunctions`CreateFAModelFile[sarahInputDirectories,
				sarahOutputDirectory, sarahModelName, eigenstates],
      kernel];
  ]

WriteParticleNamespaceFile[fileName_String] :=
  Module[{fileHandle = OpenWrite[fileName]},
    Write[fileHandle, {ToString[#], Context[#]} &
      /@ TreeMasses`GetParticles[]];
    Close[fileHandle];
  ]
  
VerticesForNPointFunction[nPointFunction_] := 
  Module[{genericVertices, genericSumPositions,
          genericInsertions, vertices},
    genericVertices = DeleteDuplicates[Cases[nPointFunction,
      SARAH`Cp[fields___] :> {fields}, Infinity, Heads -> True]];

    genericSumPositions = Position[nPointFunction[[2,1,1]], GenericSum[__]];
    genericInsertions = DeleteDuplicates[
      Flatten @ Extract[nPointFunction[[2,1,2]], genericSumPositions]];

    vertices = UniquelyInstantiateGenericFields[genericVertices,
      GatherBy[genericInsertions, First]];
    Map[StripFieldIndices, vertices, {2}]
  ]

UniquelyInstantiateGenericFields[exprs_List, {}] := exprs
UniquelyInstantiateGenericFields[exprs_List, {fieldInsertions_List, next___}] :=
  Module[{nextExprs},
    nextExprs = Flatten[exprs /. (List /@ fieldInsertions), 1];
    UniquelyInstantiateGenericFields[DeleteDuplicates[nextExprs], {next}]
  ]

CreateCXXFunctions[nPointFunctions_List, names_List,
    OptionsPattern[{LoopFunctions -> "FlexibleSUSY"}]] :=
  Module[{loopFunctionRules, hasExternalMomenta, prototypes,
          definitionHeads, definitionBodies,
          auxilliaryClasses, definitions},
    loopFunctionRules = Switch[OptionValue[LoopFunctions],
      "LoopTools", {},
      "FlexibleSUSY",
         Print["Warning: Using FlexibleSUSY loop functions will only remap A0, B0, C0 and D0."];
         Print["Warning: FlexibleSUSY loop functions C0 and D0 require zero external momenta."];
         {
           LoopTools`A0i[LoopTools`aa0, args__] :> "softsusy::a0"[Sequence @@ ("std::sqrt" /@ List[args]),
                                             "context.scale()"],
           LoopTools`A0[arg_] :> "softsusy::a0"[Sequence @@ ("std::sqrt" /@ List[arg]),
                                             "context.scale()"],
           LoopTools`B0i[LoopTools`bb0, args__] :> "softsusy::b0"[Sequence @@ ("std::sqrt" /@ List[args]),
                                             "context.scale()"],
           LoopTools`C0i[LoopTools`cc0, 0, 0, 0, args__] :> "softsusy::c0"[Sequence @@ ("std::sqrt" /@ List[args])],
           LoopTools`D0i[LoopTools`dd0, 0, 0, 0, 0, 0, 0, args__] :> "softsusy::d0"[Sequence @@ ("std::sqrt" /@ List[args])]
         },
       _, Return["Option LoopFunctions must be either LoopTools or FlexibleSUSY"]];

    hasExternalMomenta =
      FreeQ[#, SARAH`Mom[_Integer, ___]] & /@ nPointFunctions;

    prototypes = StringJoin[Riffle[
      "std::complex<double> " <> #[[2]] <>
        CXXArgStringForNPointFunctionPrototype[#[[1]]] <> ";" & /@
      Transpose[{nPointFunctions, names}], "\n"]];
      
    definitionHeads = "std::complex<double> " <> #[[2]] <>
        CXXArgStringForNPointFunctionDefinition[#[[1]]] & /@
      Transpose[{nPointFunctions, names}];
      
    definitionBodies = CXXBodyForNPointFunction /@ nPointFunctions;
    auxilliaryClasses = CXXClassForNPointFunction /@ (nPointFunctions /. loopFunctionRules);

    definitions = StringJoin[Riffle[auxilliaryClasses,"\n\n"]] <> "\n\n" <>
      StringJoin[Riffle[#[[1]] <> "\n{\n" <> #[[2]] <> "\n}" & /@
         Transpose[{definitionHeads, definitionBodies}], "\n\n"]];

    {prototypes, definitions}
  ]

CXXArgStringForNPointFunctionPrototype[nPointFunction_] :=
  Module[{numberOfIndices, numberOfMomenta},
    numberOfIndices = Length[ExternalIndicesForNPointFunction[nPointFunction]];
    numberOfMomenta = If[FreeQ[nPointFunction, SARAH`Mom[_Integer, ___]],
      0, Length[nPointFunction[[1,1]]] + Length[nPointFunction[[1,2]]]];

    "( " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates model, " <>
       "const std::array<int, " <>
      ToString[numberOfIndices] <> "> &indices, const std::array<Eigen::Vector4d, " <>
      ToString[numberOfMomenta] <> "> &momenta = { " <>
        StringJoin[Riffle[Table["Eigen::Vector4d::Zero()", {k,numberOfMomenta}],
                          ", "]] <> " } )"
  ]

CXXArgStringForNPointFunctionDefinition[nPointFunction_] :=
  Module[{numberOfIndices, numberOfMomenta},
    numberOfIndices = Length[ExternalIndicesForNPointFunction[nPointFunction]];
    numberOfMomenta = If[FreeQ[nPointFunction, SARAH`Mom[_Integer, ___]],
      0, Length[nPointFunction[[1,1]]] + Length[nPointFunction[[1,2]]]];

    "( " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates model, const std::array<int, " <>
      ToString[numberOfIndices] <> "> &indices, const std::array<Eigen::Vector4d, " <>
      ToString[numberOfMomenta] <> "> &momenta )"
  ]

CXXBodyForNPointFunction[nPointFunction_] :=
  Module[{className},
    CXXClassNameForNPointFunction[nPointFunction] <>
      " helper{ model, indices, momenta };\n" <>
    "return helper.calculate();"
  ]

CXXClassForNPointFunction[nPointFunction_] :=
  Module[{className = CXXClassNameForNPointFunction[nPointFunction],
          externalIndices, cxxCorrelationContext,
          numberOfIndices, numberOfMomenta, genericSumPositions,
          genericIndices, genericFields, genericSumNames,
          genericSumCode, preCXXRules, cxxExpr,
          subexpressions, cxxSubexpressions},
    externalIndices = ExternalIndicesForNPointFunction[nPointFunction];
    numberOfIndices = Length[externalIndices];
    numberOfMomenta = If[FreeQ[nPointFunction, SARAH`Mom[_Integer, ___]],
      0, Length[nPointFunction[[1,1]]] + Length[nPointFunction[[1,2]]]];

    genericSumPositions = Position[nPointFunction[[2,1,1]], GenericSum[__]];
    genericIndices = DeleteDuplicates[Flatten[
      Extract[nPointFunction[[2,1,1]], genericSumPositions][[All,2]], 1]];
    genericFields = #[[1]][GenericIndex[#[[2]]]] & /@ genericIndices;
    genericSumNames = Table["genericSum" <> ToString[k],
                            {k,Length[genericSumPositions]}];

    subexpressions = nPointFunction[[2,2]];
    preCXXRules = ToCXXPreparationRules[
      externalIndices, genericFields, subexpressions];

    cxxSubexpressions = 
      CXXCodeForSubexpressions[subexpressions, preCXXRules];

    genericSumCode = StringJoin[Riffle[
      CXXCodeForGenericSum[Sequence @@ #, subexpressions, preCXXRules] & /@
        Transpose[{
          Extract[nPointFunction[[2,1,1]], genericSumPositions],
          Extract[nPointFunction[[2,1,2]], genericSumPositions],
          Extract[nPointFunction[[2,1,3]], genericSumPositions],
          genericSumNames}],
      "\n\n"]];

    cxxExpr = Plus @@ ReplacePart[nPointFunction[[2,1,1]],
      Rule @@@ Transpose[{genericSumPositions, # <> "()" & /@ genericSumNames}]];
    cxxExpr = Parameters`ExpressionToString[cxxExpr];
    cxxExpr = StringReplace[cxxExpr, "\"" -> ""];

    cxxCorrelationContext = "correlation_function_context<" <>
       ToString[numberOfIndices] <> ", " <> ToString[numberOfMomenta] <>
    ">";
    
    "class " <> className <> "\n" <>
    ": public " <> cxxCorrelationContext <> "\n{\n" <>
    "using generic_sum_base = " <> cxxCorrelationContext <> ";\n" <>
    "template<class GenericFieldMap> struct subexpression_base\n" <>
    ": generic_sum_base,\n" <>
    "  index_map_interface<GenericFieldMap>\n" <>
    "{\n" <> 
    "subexpression_base( const subexpression_base & ) = default;\n" <>
    "subexpression_base( const generic_sum_base &gsb, \n" <>
    "const typename field_index_map<GenericFieldMap>::type &fim )\n" <>
    ": generic_sum_base( gsb ), index_map_interface<GenericFieldMap>( fim )\n" <>
    "{}\n" <>
    "};\n\n" <>

    StringJoin[Riffle[
        "struct " <> CXXKeyOfGenericField[#] <> " {};" & /@ genericFields,
      "\n"]] <> "\n\n" <>

    cxxSubexpressions <> "\n\n" <> 
    genericSumCode <> "\n\n" <>
    "public:\n" <>
    className <> CXXArgStringForNPointFunctionDefinition[nPointFunction] <> "\n" <>
    ": " <> cxxCorrelationContext <> "{ model, indices, momenta }\n{}\n\n" <>
    
    "std::complex<double> calculate( void )\n{\nreturn " <>
    cxxExpr <>
    ";\n}" <>
    "\n};"
  ]

ToCXXPreparationRules[externalIndices_List, 
    genericFields_List, subexpressions_List] :=
  Module[{externalIndexRules, genericRules,
          subexprRules, massRules, couplingRules},
    externalIndexRules = Rule[#[[1]],
        "this->external_indices(" <> ToString[#[[2]]] <> ")"] & /@
      Transpose[{externalIndices, Table[k, {k, Length[externalIndices]}] - 1}];

    genericRules = Join[
      Rule[Susyno`LieGroups`conj[#], CXXNameOfGenericField[Susyno`LieGroups`conj[#]][
        CXXIndicesForField[#]]] & /@ genericFields,
      Rule[SARAH`bar[#], CXXNameOfGenericField[SARAH`bar[#]][
        CXXIndicesForField[#]]] & /@ genericFields,
      Rule[#, CXXNameOfGenericField[#][
        CXXIndicesForField[#]]] & /@ genericFields
    ];

    couplingRules = {
      SARAH`Cp[fields___][1] :>
      I * "context.vertex<" <> StringJoin[Riffle[
        If[IsGenericField[#], Head[# /. genericRules],
           CXXNameOfField[Vertices`StripFieldIndices[#],
             prefixNamespace -> "fields"]] & /@ {fields},
        ", "]] <>
      ">( lorentz_scalar{}, concatenate( " <>
        StringJoin[Riffle[CXXIndicesForField /@ {fields}, ", "]] <> 
      " ) )",
      SARAH`Cp[fields___][SARAH`PL] :>
      I * "context.vertex<" <> StringJoin[Riffle[
        If[IsGenericField[#], Head[# /. genericRules],
           CXXNameOfField[Vertices`StripFieldIndices[#],
             prefixNamespace -> "fields"]] & /@ {fields},
        ", "]] <>
      ">( lorentz_left{}, concatenate( " <>
        StringJoin[Riffle[CXXIndicesForField /@ {fields}, ", "]] <> 
      " ) )",
      SARAH`Cp[fields___][SARAH`PR] :>
      I * "context.vertex<" <> StringJoin[Riffle[
        If[IsGenericField[#], Head[# /. genericRules],
           CXXNameOfField[Vertices`StripFieldIndices[#],
             prefixNamespace -> "fields"]] & /@ {fields},
        ", "]] <>
      ">( lorentz_right{}, concatenate( " <>
        StringJoin[Riffle[CXXIndicesForField /@ {fields}, ", "]] <> 
      " ) )",
      SARAH`Cp[fields___][SARAH`Mom[f1_] - SARAH`Mom[f2_]] :>
      I * "context.vertex<" <> StringJoin[Riffle[
        If[IsGenericField[#], Head[# /. genericRules],
           CXXNameOfField[Vertices`StripFieldIndices[#],
             prefixNamespace -> "fields"]] & /@ {fields},
        ", "]] <>
      ">( lorentz_momentum_diff{" <>
            ToString[Position[{fields},f1,{1},Heads->False][[1,1]]-1] <> ", " <>
            ToString[Position[{fields},f2,{1},Heads->False][[1,1]]-1] <> "}, " <>
        "concatenate( " <>
        StringJoin[Riffle[CXXIndicesForField /@ {fields}, ", "]] <> 
      " ) )",
      SARAH`Cp[fields___][SARAH`g[_,_]] :>
      I * "context.vertex<" <> StringJoin[Riffle[
        If[IsGenericField[#], Head[# /. genericRules],
           CXXNameOfField[Vertices`StripFieldIndices[#],
             prefixNamespace -> "fields"]] & /@ {fields},
        ", "]] <>
      ">( lorentz_inverse_metric{}, concatenate( " <>
        StringJoin[Riffle[CXXIndicesForField /@ {fields}, ", "]] <> 
      " ) )"
    };
    massRules = SARAH`Mass[field_String[indices_String]] :>
      "context.mass<" <> field <> ">( " <> indices <> " )";
    subexprRules = Rule[#[[1]], ToString[#[[1]]] <> "_()"] & /@
      subexpressions;

    {externalIndexRules, couplingRules,
     genericRules, massRules, subexprRules}
  ]

CXXCodeForSubexpressions[subexpressions_List, preCXXRules_List] :=
  Module[{names = ToString /@ subexpressions[[All,1]],
          exprs = subexpressions[[All,2]], subexpr,
          relevantSubexpressions, relevantGenericFields,
          needsContexts, cxxExprs},
    relevantGenericFields = DeleteDuplicates /@
      (Cases[#,_[GenericIndex[_]], Infinity, Heads -> True] & /@ exprs);

    relevantSubexpressions = DeleteDuplicates /@
      (Cases[#,
        Pattern[subexpr, Alternatives @@ subexpressions[[All,1]]] :> subexpr,
          Infinity,
            Heads -> True] & /@ exprs);

    needsContexts = !(FreeQ[#, SARAH`Cp] && FreeQ[#, SARAH`Mass]) & /@ exprs;

    cxxExprs = Parameters`ExpressionToString[
      Fold[ReplaceAll, #, preCXXRules]] & /@ exprs;
    cxxExprs = StringReplace[#, "\"" -> ""] & /@ cxxExprs;

    StringJoin[Riffle[
    Module[{name = #[[1]], genericFields = #[[2]],
       subs = #[[3]], needsContext = #[[4]], cxxExpr = #[[5]]},
    "template<class GenericFieldMap> struct " <> name <>
    "\n: subexpression_base<GenericFieldMap>\n{\n" <>
    "template<class ...Args> " <> name <> "( Args &&...args )\n" <>
    ": subexpression_base<GenericFieldMap>( std::forward<Args>( args )... )\n{}\n\n" <>
    "std::complex<double> operator()( void ) const\n{\n" <>

    If[Length[subs] =!= 0,
      StringJoin[ToString[#] <> "<GenericFieldMap> " <>
        ToString[#] <> "_{ *this };\n" & /@
        subs] <> "\n",
      ""] <>

    If[Length[genericFields] =!= 0,
      "using boost::mpl::at;\n" <>
      "using boost::fusion::at_key;\n\n" <>

      StringJoin[Module[{genericField = #},
        "using " <> CXXNameOfGenericField[genericField] <>
        " = typename at<GenericFieldMap, " <>
          CXXKeyOfGenericField[genericField] <>
        ">::type;\n"
        ] & /@ genericFields] <> "\n" <>

      StringJoin[Module[{genericField = #},
        "const auto &" <> CXXIndicesForField[genericField] <>
        " = at_key<" <> CXXKeyOfGenericField[genericField] <>
          ">( this->index_map() );\n"
        ] & /@ genericFields] <> "\n",
      ""] <>

    If[needsContext,
       "const context_with_vertices &context =  *this;\n\n",
       ""] <>
    
    "return " <> cxxExpr <> ";\n}\n};"] & /@
      Transpose[{names, relevantGenericFields,
                 relevantSubexpressions, needsContexts, cxxExprs}],
    "\n\n"]]
  ]

CXXCodeForGenericSum[sum_GenericSum, genericInsertions_List,
    combinatorialFactors_List, functionName_String,
    subexpressions_List, preCXXRules_List] :=
  Module[{expr = sum[[1]], indices = sum[[2]],
          sortedGenericInsertions, genericFields, relevantSubexpressions,
          subexpr, needsContext, cxxExpr},

    relevantSubexpressions = DeleteDuplicates[Cases[expr,
      Pattern[subexpr, Alternatives @@ subexpressions[[All,1]]] :> subexpr,
      Infinity, Heads -> True]];

    needsContext = !(FreeQ[expr, SARAH`Cp] && FreeQ[expr, SARAH`Mass]);

    cxxExpr = Parameters`ExpressionToString[
      Fold[ReplaceAll, expr, preCXXRules]];
    cxxExpr = StringReplace[cxxExpr, "\"" -> ""];

    genericFields = Sort[indices[[#,1]][GenericIndex[indices[[#,2]]]] & /@
      Table[k, {k,Length[indices]}]];
    sortedGenericInsertions = SortBy[#, First] & /@ genericInsertions;

    "template<class GenericFieldMap> struct " <> functionName <> "_impl\n" <>
    ": generic_sum_base\n{\n" <>
    functionName <> "_impl( const generic_sum_base &base )\n" <>
    ": generic_sum_base( base ) {} \n\n" <>
    "std::complex<double> operator()( void )\n{\n" <>
    "using boost::mpl::at;\n" <>
    "using boost::fusion::at_key;\n\n" <>

    StringJoin[
      "using " <> CXXNameOfGenericField[#] <>
      " = typename at<GenericFieldMap, " <>
        CXXKeyOfGenericField[#] <>
      ">::type;\n" &
      /@ genericFields] <> "\n" <>

    "typename field_index_map<GenericFieldMap>::type index_map;\n\n" <>

    If[Length[relevantSubexpressions] =!= 0,
       StringJoin[ToString[#] <> "<GenericFieldMap> " <>
         ToString[#] <> "_{ *this, index_map };\n" & /@
         relevantSubexpressions] <> "\n",
       ""] <>

    If[needsContext,
       "const context_with_vertices &context = *this;\n",
       ""] <>
    "std::complex<double> value = 0.0;\n\n" <>

    StringJoin[Riffle[
      "for( const auto &" <> CXXIndicesForField[#] <> " : " <>
        "index_range<" <> CXXNameOfGenericField[#] <> ">() ) {\n" <>
      "at_key<" <> CXXKeyOfGenericField[#] <> ">( index_map ) = " <>
        CXXIndicesForField[#] <> ";" & /@
      genericFields, "\n"]] <> "\n\n" <>
    "value += " <> cxxExpr <> ";\n" <>
    StringJoin[Table["}",{k,Length[indices]}]] <> "\n\n" <> 

    "return value;" <>
    "\n}\n};\n\n" <>

    "std::complex<double>" <> functionName <> "( void )\n{\n" <>
    "using GenericKeys = boost::mpl::vector<\n" <>
    StringJoin[Riffle[CXXKeyOfGenericField /@ genericFields, ",\n"]] <>
    "\n>;\n\n" <>

    "using GenericInsertions" <> " = boost::mpl::vector<\n" <>
    StringJoin[Riffle[
      "boost::mpl::vector<" <> StringJoin[Riffle[
          CXXDiagrams`CXXNameOfField[#, prefixNamespace -> "fields"] & /@ #, ", "]] <>
        ">" & /@ sortedGenericInsertions[[All,All,2]],
      ",\n"]] <> "\n>;\n\n" <>

    "using combinatorial_factor" <> " = boost::mpl::vector<\n" <>
    StringJoin[Riffle[
      "boost::mpl::int_<" <> ToString[#] <> ">" & /@ combinatorialFactors,
      ", "]] <> "\n>;\n\n" <>

    "return accumulate_generic<GenericKeys, GenericInsertions,\n" <> 
      "combinatorial_factor, " <> functionName <> 
    "_impl>( *this );\n}"
  ]

ExternalIndicesForNPointFunction[nPointFunction_] :=
  Flatten[Cases[Join[nPointFunction[[1,1]], nPointFunction[[1,2]]],
            _[indices_List] :> indices]]

IndicesForGenericField[_GenericS] :=
  Table[k, {k,Length[
    Select[TreeMasses`GetParticles[],TreeMasses`IsScalar]]}]
IndicesForGenericField[_GenericF] :=
  Table[k, {k,Length[
    Select[TreeMasses`GetParticles[],TreeMasses`IsFermion]]}]
IndicesForGenericField[_GenericV] :=
  Table[k, {k,Length[
    Select[TreeMasses`GetParticles[],TreeMasses`IsVector]]}]
IndicesForGenericField[_GenericU] :=
  Table[k, {k,Length[
    Select[TreeMasses`GetParticles[],TreeMasses`IsGhost]]}]

CXXIndicesForField[SARAH`bar[field_]] := CXXIndicesForField[field]
CXXIndicesForField[Susyno`LieGroups`conj[field_]] := CXXIndicesForField[field]
CXXIndicesForField[head_[GenericIndex[index_Integer]]] := 
  "indices" <> StringTake[SymbolName[head],-1] <> ToString[index]
CXXIndicesForField[field_] :=
  If[Length[field] === 0, "std::array<int, 0>()",
     "std::array<int, " <> ToString[Length[field[[1]]]] <> ">{ " <>
       StringJoin[Riffle[ToString /@ field[[1]], ", "]] <> " }"]

CXXIndexTypeOfGenericField[field_] :=
  "typename field_indices<" <> CXXNameOfGenericField[field] <> ">::type"

IsGenericField[field_] :=
  Module[{head = Head[CXXDiagrams`RemoveLorentzConjugation[field]]},
    Switch[head,
           GenericS, True, 
           GenericF, True, 
           GenericV, True, 
           GenericU, True, 
           _, False]
  ] 

CXXKeyOfGenericField[head_[GenericIndex[index_Integer]]] :=
  ToString[head] <> ToString[index] <> "Key"

CXXNameOfGenericField[head_[GenericIndex[index_Integer]]] :=
  ToString[head] <> ToString[index]
CXXNameOfGenericField[SARAH`bar[field_]] :=
  "typename bar<" <> CXXNameOfGenericField[field] <> ">::type";
CXXNameOfGenericField[Susyno`LieGroups`conj[field_]] :=
  "typename conj<" <> CXXNameOfGenericField[field] <> ">::type";

CXXGenericFieldVector[GenericS[_]] := "fields::scalars"
CXXGenericFieldVector[GenericF[_]] := "fields::fermions"
CXXGenericFieldVector[GenericV[_]] := "fields::vectors"
CXXGenericFieldVector[GenericU[_]] := "fields::ghosts"

InstancesOfGenericField[GenericS[_]] :=
  Select[TreeMasses`GetParticles[], TreeMasses`IsScalar]
InstancesOfGenericField[GenericF[_]] :=
  Select[TreeMasses`GetParticles[], TreeMasses`IsFermion]
InstancesOfGenericField[GenericV[_]] :=
  Select[TreeMasses`GetParticles[], TreeMasses`IsVector]
InstancesOfGenericField[GenericU[_]] :=
  Select[TreeMasses`GetParticles[], TreeMasses`IsGhost]

CXXClassNameForNPointFunction[nPointFunction_] :=
  Module[{fields},
    fields = Vertices`StripFieldIndices[
      Join[nPointFunction[[1,1]], nPointFunction[[1,2]]]];
    "nPoint" <> StringJoin[ToString /@ Flatten[fields //. a_[b_] :> {a,b}]]
  ]

CacheNPointFunction[nPointFunction_, nPointMeta_List] := 
  Module[{sarahOutputDir = SARAH`$sarahCurrentOutputMainDir,
          outputDir, nPointFunctionsDir, cacheName,
          nPointFunctionsFile,fileHandle,nPointFunctions,
          position},
    cacheName = CacheNameForMeta[nPointMeta];

    outputDir = FileNameJoin[{sarahOutputDir, ToString[FlexibleSUSY`FSEigenstates]}];
    nPointFunctionsDir = FileNameJoin[{outputDir, "NPointFunctions"}];
    nPointFunctionsFile = FileNameJoin[{nPointFunctionsDir, cacheName}];

    If[!FileExistsQ[nPointFunctionsFile],
       fileHandle = OpenWrite[nPointFunctionsFile];
       Write[fileHandle,{}];
       Close[fileHandle]];

    nPointFunctions = Get[nPointFunctionsFile];

    position = Position[nPointFunctions[[All,1]], nPointFunction[[1]]];
    If[Length[position] =!= 0,
       nPointFunctions[[position[[1,1]]]] = nPointFunction,

       AppendTo[nPointFunctions, nPointFunction]
    ];

    fileHandle = OpenWrite[nPointFunctionsFile];
    Write[fileHandle,nPointFunctions];
    Close[fileHandle];
  ]

CachedNPointFunction[inFields_List, outFields_List, nPointMeta_List] :=
  Module[{sarahOutputDir = SARAH`$sarahCurrentOutputMainDir,
          outputDir, nPointFunctionsDir, cacheName,
          nPointFunctionsFile, nPointFunctions,
          unindexedExternalFields, position},
    cacheName = CacheNameForMeta[nPointMeta];

    outputDir = FileNameJoin[{sarahOutputDir, ToString[FlexibleSUSY`FSEigenstates]}];
    nPointFunctionsDir = FileNameJoin[{outputDir, "NPointFunctions"}];
    nPointFunctionsFile = FileNameJoin[{nPointFunctionsDir, cacheName}];

    If[!FileExistsQ[nPointFunctionsFile],
       Return[Null]];

    nPointFunctions = Get[nPointFunctionsFile];

    unindexedExternalFields = {Vertices`StripFieldIndices[#[[1]]],
                               Vertices`StripFieldIndices[#[[2]]]} & /@
      nPointFunctions[[All,1]];

    position = Position[unindexedExternalFields,
      {inFields, outFields}, {1}, Heads -> False];
    If[Length[position] === 0,
       Null,
       nPointFunctions[[position[[1,1]]]]
    ]
  ]

CacheNameForMeta[nPointMeta_List] :=
  "cache_" <> StringJoin[Riffle[ToString /@ nPointMeta, "_"]] <> ".m"

FeynArtsNamesForFields[fields_List,particleNamesFile_String] :=
  Module[{lines, unindexedBaseFields, fieldLines,
          faFieldNames, faNameRules},
    lines = Utils`ReadLinesInFile[particleNamesFile];
    
    unindexedBaseFields = DeleteDuplicates[CXXDiagrams`AtomHead[
      CXXDiagrams`RemoveLorentzConjugation[#]] & /@ fields];

    fieldLines = Module[{fieldName = ToString[#]},
      Select[lines, StringMatchQ[#,___~~fieldName~~":"~~___] &][[1]]] & /@
        unindexedBaseFields;

    faFieldNames = ("FeynArts`" <> StringSplit[#][[2]]) & /@ fieldLines;
    
    faNameRules = Join[
      Sequence @@ {#[[1]][indices_List] :>
         StringDrop[#[[2]], -1] <> ", " <> ToString[indices] <> "]",
       #[[1]] -> #[[2]]} & /@ Transpose[{unindexedBaseFields, faFieldNames}],
      {SARAH`bar[field_String] :> "-" <> field,
       Susyno`LieGroups`conj[field_String] :> "-" <> field}
    ];
    
    fields //. faNameRules
  ]

GenericFieldType[field_] :=
  Module[{head = CXXDiagrams`AtomHead[CXXDiagrams`RemoveLorentzConjugation[field]]},
    Switch[head,
           GenericS, GenericS,
           GenericF, GenericF,
           GenericV, GenericV,
           GenericU, GenericU,
           GenericT, GenericT,
           _, If[TreeMasses`IsScalar[head],GenericS,
                 If[TreeMasses`IsFermion[head],GenericF,
                    If[TreeMasses`IsVector[head],GenericV,
                      If[TreeMasses`IsGhost[head],GenericU,
                         "unrecognized field type " <> ToString[head]]]]]]

  ]

CreateCXXHeaders[OptionsPattern[{LoopFunctions -> "FlexibleSUSY"}]] :=
  "#include \"cxx_qft/" <> FlexibleSUSY`FSModelName <> "_npointfunctions.hpp\"\n" <>
  "#include \"concatenate.hpp\"\n\n" <>
  "#include <boost/fusion/include/at_key.hpp>\n\n" <>
  Switch[OptionValue[LoopFunctions],
         "LoopTools", "#include <clooptools.h>",
         "FlexibleSUSY", "#include \"numerics.h\""]
     

End[];
EndPackage[];




BeginPackage["EDM`", {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "Parameters`"}];

(* This module generates c++ code that calculates electric dipole moments of fields *)

Initialize::usage="Initialize the EDM module.";

SetEDMFields::usage="Set the fields for which the EDMs shall be calculated.";
CreateFields::usage="Returns the c++ code that contains all field fields";
CreateDiagrams::usage="Returns the c++ code that contains all relevant diagram classes";
CreateInterfaceFunctions::usage="Returns the c++ code containing the interface functions {prototypeCode, definitionCode}."
CreateVertexFunctionData::usage="Returns the c++ code that contains all relevant vertex function data";
CreateDefinitions::usage="Returns the c++ that contains all function definitions"

NPointFunctions::usage="Returns a list of all n point functions that are needed. Actually it is a list of fake functions to extract vertex functions...";

(******** TODO: IMPORTANT NOTES:
 If you add new kinds of vertices (e.g for new diagram types):
 - Add the new types to vertexTypes
 - Expand CouplingsForFields[] and VertexTypeForFields[] accordingly
 - Write the c++ class for the new vertex type
 
 When adding support for new diagram types, do the following:
 - Add the new types to diagramTypes
 - Write new overloads for CreateDiagramEvaluatorClass[], ContributingDiagramsOfType[] and VerticesForDiagram[]
 - Write the necessary c++ code: loop functions, DiagramEvaluator<> specialisations
 **********)

(************* Begin public interface *******************)

Initialize[] := (subIndexPattern = (Alternatives @@ SARAH`subIndizes[[All, 1]] -> ___);)

edmFields = Null;
SetEDMFields[fields_List] := (edmFields = fields;)

CreateFields[] :=
Module[{fields, code},
       fields = TreeMasses`GetParticles[];
       
       code = ("struct Field {};\n\n" <>
               
               StringJoin @ Riffle[("struct " <> CXXNameOfField[#] <>
                                    ": public Field {\n" <>
                                    TextFormatting`IndentText["static const unsigned numberOfGenerations = " <>
                                                              ToString @ TreeMasses`GetDimension[#] <> ";\n"] <>
                                    "};\n" &) /@ fields, "\n"] <> "\n\n" <>
               "// Special field families\n" <>
               "using Photon = " <> CXXNameOfField @ SARAH`Photon <> ";\n" <>
               "using Electron = " <> CXXNameOfField @ SARAH`Electron <> ";\n\n" <>
               
               "// Anti fields\n" <>
               "template<class P> struct anti : public Field\n" <>
               "{\n" <>
               IndentText @
               ("static const unsigned numberOfGenerations = P::numberOfGenerations;\n" <>
                "using type = anti<P>;\n") <>
               "};\n" <>
               "template<class P> struct anti<anti<P>> { using type = P; };\n\n" <>
               
               "// Fields that are their own anti fields\n" <>
               StringJoin @ Riffle[("template<> struct " <>
                                    "anti<" <> CXXNameOfField[#] <> ">" <>
                                    " { using type = " <> CXXNameOfField[#] <> "; };"
                                    &) /@ Select[fields, (# == SARAH`AntiField[#] &)],
                                   "\n"] <> "\n\n" <>
               
               "template<class Field> struct field_indices;\n\n" <>
               
               StringJoin @ Riffle[("template<> struct field_indices<" <>
                                    CXXNameOfField[#] <> ">\n" <>
                                    "{\n" <>
                                    IndentText @
                                    ("using type = std::array<unsigned, " <>
                                     ToString @ Length @ CleanFieldInfo[#][[5]] <>
                                     ">;\n"
                                    ) <>
                                    "};\n" &) /@ fields, "\n"]
               );
       
       Return[code];
       ];

CreateDiagrams[] :=
Module[{code},
       code = StringJoin @ Riffle[(Module[{diagramType = #},
                                  "template<unsigned> class " <> SymbolName[diagramType] <> ";\n" <>
                                   StringJoin @ Riffle[("template<> class " <> SymbolName[diagramType] <>
                                                        "<" <> ToString @ # <> "> {};"
                                                        &) /@ diagramSubTypes[diagramType], "\n"]
                                   ] &) /@ diagramTypes, "\n\n"];
       
       code = (code <> "\n\n" <>
               StringJoin @ Riffle[(Module[{diagramType = #},
                                    StringJoin @ Riffle[
                                    ("template<class EDMField, class PhotonEmitter, class ExchangeField>\n" <>
                                    "struct DiagramEvaluator<" <> SymbolName[diagramType] <>
                                    "<" <> ToString @ # <>
                                    ">, EDMField, PhotonEmitter, ExchangeField>\n" <>
                                    "{ static double value(const typename field_indices<EDMField>::type &indices, EvaluationContext& context); };"
                                    &) /@ diagramSubTypes[diagramType], "\n\n"]] &) /@ diagramTypes, "\n\n"]
               );
       
       Return[code];
       ];

CreateInterfaceFunctions[] :=
Module[{prototypes, definitions, evaluators},
       evaluators = ConcreteDiagramEvaluators[];
       
       prototypes = ("namespace " <> FlexibleSUSY`FSModelName <> "_edm {\n" <>
                     StringJoin @ Riffle[("double calculate_edm_" <> CXXNameOfField[#] <>
                                          "(" <>
                                          If[TreeMasses`GetDimension[#] =!= 1,
                                             " unsigned generationIndex, ",
                                             " "] <>
                                          "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model );"
                                          &) /@ edmFields, "\n"] <>
                     "\n}");
       
       definitions = StringJoin @ Riffle[
                         Module[{field = #[[1]], fieldEvaluators = #[[2]],
                                 numberOfIndices},
                                numberOfIndices = Length @ CleanFieldInfo[field][[5]];
                                
                                "double " <> FlexibleSUSY`FSModelName <> "_edm::calculate_edm_" <> CXXNameOfField[field] <>
                                "(" <>
                                If[TreeMasses`GetDimension[field] =!= 1,
                                   " unsigned generationIndex, ",
                                   " "] <>
                                "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model )\n" <>
                                "{\n" <>
                                IndentText @
                                ("CMSSM_mass_eigenstates model_ = model;\n" <>
                                 "EvaluationContext context{ model_ };\n" <>
                                 "std::array<unsigned, " <>
                                 ToString @ numberOfIndices <>
                                 "> indices = {" <>
                                 If[TreeMasses`GetDimension[field] =!= 1,
                                    " generationIndex" <>
                                    If[numberOfIndices =!= 1,
                                       StringJoin @ Table[", 1", {numberOfIndices-1}],
                                       ""
                                       ] <> " ",
                                    If[numberOfIndices =!= 0,
                                       StringJoin @ Riffle[Table[" 1", {numberOfIndices}], ","] <> " ",
                                       ""
                                       ]
                                    ] <>
                                 "};\n\n" <>
                                 
                                 "double val = 0.0;\n\n" <>
                                 StringJoin @ Riffle[("val += " <> ToString @ # <>
                                                      "::value(indices, context);"
                                                      &) /@ fieldEvaluators, "\n"] <>
                                 "\n\n" <>
                                 "return val;"
                                 ) <>
                                "\n}"] & /@ evaluators, "\n\n"];
       
       Return[{prototypes, definitions}];
       ];

CreateVertexFunctionData[vertexRules_List] := CreateVertices[vertexRules][[1]];

CreateDefinitions[vertexRules_List] :=
(CreateVertices[vertexRules][[2]] <> "\n\n" <>
 CreateEvaluationContextSpecializations[]);

NPointFunctions[] :=
Module[{contributingDiagrams, vertices},
       contributingDiagrams = ContributingDiagrams[];
       
       vertices = DeleteDuplicates @ Flatten[VerticesForDiagram /@
                                             Flatten @ contributingDiagrams[[All, 2]], 1];
       
       Flatten[(Null[Null, #] &) /@ ((CouplingsForFields[#] &) /@ vertices)]
       ];

(**************** End public interface *****************)
Begin["`Private`"];

(* The supported vertex types.
 They have the same names as their c++ counterparts. *)
vertexTypes = {
    SingleComponentedVertex,
    LeftAndRightComponentedVertex
};

(* The supported diagram types.
 They have the same names as their c++ counterparts. *)
diagramTypes = {
    OneLoopDiagram
};

(* The supported diagram types.
 Indexed by the diagram type, gives a set of (c++-compatible) unsigned integer indices. *)
diagramSubTypes[OneLoopDiagram] = { 0, 1 }; (* 0: fermion emits photon, exchange field is a scalar
                                             1: scalar emits photon, exchange field is a fermion *)

(**************** CXX conversion routines ***************)

(* Return a string corresponding to the c++ class name of the field.
 Note that "bar" and "conj" get turned into anti<...>::type! *)
CXXNameOfField[p_] := SymbolName[p];
CXXNameOfField[SARAH`bar[p_]] := "anti<" <> SymbolName[p] <> ">::type";
CXXNameOfField[Susyno`LieGroups`conj[p_]] := "anti<" <> SymbolName[p] <> ">::type";

(**************** Other Functions ***************)

GetElectronIndex[] := If[TreeMasses`GetDimension[SARAH`Electron] =!= 1, 1, Null];

IsLorentzIndex[index_] := StringMatchQ[ToString @ index, "lt" ~~ __];

StripLorentzIndices[p_Symbol] := p;
StripLorentzIndices[SARAH`bar[p_]] := SARAH`bar[StripLorentzIndices[p]];
StripLorentzIndices[Susyno`LieGroups`conj[p_]] := Susyno`LieGroups`conj[StripLorentzIndices[p]];
StripLorentzIndices[p_] := Module[{remainingIndices},
                                  remainingIndices = Select[p[[1]], (!IsLorentzIndex[#] &)];
                                  If[Length[remainingIndices] === 0, Head[p],
                                     Head[p][remainingIndices]]
                                  ];

CleanFieldInfo[field_] := Module[{fieldInfo = FirstCase[SARAH`Particles[FlexibleSUSY`FSEigenstates],
                                                        {SARAH`getParticleName @ field, ___}]},
                                 fieldInfo = DeleteCases[fieldInfo, {SARAH`generation, 1}, {2}];
                                 DeleteCases[fieldInfo, {SARAH`lorentz, _}, {2}]
                                 ];

CreateEvaluationContextSpecializations[] :=
Module[{fields, contributingDiagrams, photonEmitters,
        electronDimension = TreeMasses`GetDimension[SARAH`Electron]},
       fields = Select[TreeMasses`GetParticles[], (! TreeMasses`IsGhost[#] &)];
       
       contributingDiagrams = ContributingDiagrams[];
       photonEmitters = DeleteDuplicates @ Flatten[contributingDiagrams[[All, 2]], 1][[All, 3]];

       StringJoin @ Riffle[(Module[{fieldInfo = CleanFieldInfo[#], numberOfIndices},
                                   numberOfIndices = Length @ fieldInfo[[5]];
                                   
                                   "template<> double EvaluationContext::mass<" <> ToString[#] <>
                                   ">( const std::array<unsigned, " <> ToString @ numberOfIndices <>
                                   "> &indices ) const\n" <>
                                   "{ return model.get_M" <> CXXNameOfField[#] <>
                                   If[TreeMasses`GetDimension[#] === 1, "()", "( indices[0] )"] <> "; }"
                                   ] &) /@ fields, "\n\n"] <> "\n\n" <>
       
       StringJoin @ Riffle[(Module[{fieldInfo = CleanFieldInfo[#],
                                    photonVertexType = VertexTypeForFields[{SARAH`Photon, #, SARAH`AntiField @ #}],
                                    numberOfIndices},
                                   numberOfIndices = Length @ fieldInfo[[5]];
                                   
                                   "template<>\n" <>
                                   "double EvaluationContext::charge<" <> CXXNameOfField[#] <>
                                   ">( const std::array<unsigned, " <> ToString @ numberOfIndices <>
                                   "> &indices ) const\n" <>
                                   "{\n" <>
                                   IndentText @
                                   ("using PhotonVertex = VertexFunction<Photon, " <>
                                    CXXNameOfField[#] <> ", " <> CXXNameOfField[SARAH`AntiField @ #] <>
                                    ">;\n\n" <>
                                    "double fieldCharge = PhotonVertex::vertex( concatenate( indices, indices ), *this )" <>
                                    If[photonVertexType === SingleComponentedVertex,
                                       ".value().real();\n",
                                       ".left().real();\n"] <>
                                    "return fieldCharge;\n"
                                    ) <>
                                   "}"] &) /@ photonEmitters, "\n\n"]
       ];

(* Find all diagrams of the type type_, testing all corresponding combinations of fields *)
(* IMPORTANT: Return value should have the format
 {{edmField1, {Diagram[DIAGRAMTYPENAME[_Integer], Fields___], Diagram[...], ...}},
  {edmField2, {...}},
  ...} *)
ContributingDiagramsOfType[OneLoopDiagram] :=
    Module[{edmField = #, diagrams = SARAH`InsFields[
                           {{C[#, SARAH`AntiField[SARAH`FieldToInsert[1]],
                               SARAH`AntiField[SARAH`FieldToInsert[2]]],
                             C[SARAH`FieldToInsert[1], SARAH`Photon,
                               SARAH`AntiField[SARAH`FieldToInsert[1]]],
                             C[SARAH`FieldToInsert[1], SARAH`FieldToInsert[2],
                               SARAH`AntiField[#]]},
                            {SARAH`FieldToInsert[1], SARAH`FieldToInsert[2]}}],
            subtypedDiagrams, uniqueDiagrams},
           
           If[TreeMasses`IsFermion[edmField] =!= True,
              Return[{edmField,{}}]];
           
           subtypedDiagrams = (Module[{photonEmitter = #[[2,1]],
                                       exchangeField = #[[2,2]],
                                       subType},
                                      subType = If[TreeMasses`IsFermion[photonEmitter] &&
                                                   TreeMasses`IsScalar[exchangeField],
                                                   0,
                                                   If[TreeMasses`IsScalar[photonEmitter] &&
                                                      TreeMasses`IsFermion[exchangeField],
                                                      1]];
                                      If[subType === Null,
                                         Null,
                                         Diagram[OneLoopDiagram[subType], edmField, photonEmitter, exchangeField]]
                                      ]
                               &) /@ diagrams;
           uniqueDiagrams = DeleteDuplicates @ Cases[subtypedDiagrams, Except[Null]];
           
           {edmField, uniqueDiagrams}] & /@ edmFields;

(* Returns the necessary c++ code corresponding to the vertices that need to be calculated.
 The returned value is a list {prototypes, definitions}. *)
CreateVertices[vertexRules_List] :=
    Module[{contributingDiagrams, vertices,
            vertexClassesPrototypes, vertexClassesDefinitions},
           contributingDiagrams = ContributingDiagrams[];

           vertices = DeleteDuplicates @ Flatten[VerticesForDiagram /@
                                                 Flatten @ contributingDiagrams[[All, 2]], 1];
           
           {vertexClassesPrototypes, vertexClassesDefinitions} = Transpose @
               ((CreateVertexFunction[#, vertexRules] &) /@ vertices);

           (StringJoin @ Riffle[#, "\n\n"] &) /@ {vertexClassesPrototypes, vertexClassesDefinitions}
          ];

(* Returns the vertices that are present in the specified diagram.
 This function should be overloaded for future diagram types. *)
VerticesForDiagram[Diagram[loopDiagram_OneLoopDiagram, edmField_, photonEmitter_, exchangeField_]] :=
    Module[{edmVertex, photonVertex},
           edmVertex = {edmField, SARAH`AntiField[photonEmitter], SARAH`AntiField[exchangeField]};
           photonVertex = {SARAH`Photon, photonEmitter, SARAH`AntiField[photonEmitter]};
           Return[{edmVertex, photonVertex}];
           ];

(* Returns the vertex type for a vertex with a given list of fields *)
VertexTypeForFields[fields_List] :=
    Module[{fermions, scalarCount, vectorCount, fermionCount},
           fermions = Select[fields, TreeMasses`IsFermion];
           
           scalarCount = Length @ Select[fields, TreeMasses`IsScalar];
           vectorCount = Length @ Select[fields, TreeMasses`IsVector];
           fermionCount = Length @ fermions;
           
           If[fermionCount === 2 && scalarCount === 1 && vectorCount === 0,
              Return[LeftAndRightComponentedVertex]];
           If[fermionCount === 2 && scalarCount === 0 && vectorCount === 1,
              If[fermions[[1]] === SARAH`AntiField[fermions[[2]]],
                 Return[LeftAndRightComponentedVertex]]];
           If[fermionCount === 0 && scalarCount === 2 && vectorCount === 1,
              Return[SingleComponentedVertex]];
           
           Return[Null];
           ];

(* Returns the different SARAH`Cp coupling parts for a vertex with a given list of fields *)
CouplingsForFields[fields_List] :=
    Module[{vertexType, couplings},
           vertexType = VertexTypeForFields[fields];
           couplings = {SARAH`Cp @@ fields};

           Switch[vertexType,
                  SingleComponentedVertex, couplings,
                  LeftAndRightComponentedVertex, {couplings[[1]][SARAH`PL], couplings[[1]][SARAH`PR]}]
           ];

(* Creates the actual c++ code for a vertex with given fields.
 This involves creating the VertexFunctionData<> code as well as
 the VertexFunction<> code. You should never need to change this code! *)
CreateVertexFunction[fields_List, vertexRules_List] :=
    Module[{prototype, definition,
        parsedVertex, dataClassName, functionClassName, fieldIndexStartF,
        fieldIndexStart, indexBounds},
           parsedVertex = ParseVertex[fields, vertexRules];
           
           dataClassName = "VertexFunctionData<" <> StringJoin @ Riffle[CXXNameOfField /@ fields, ", "] <> ">";
           functionClassName = "VertexFunction<" <> StringJoin @ Riffle[CXXNameOfField /@ fields, ", "] <> ">";
           
           fieldIndexStartF[1] = 0;
           fieldIndexStartF[pIndex_] := fieldIndexStartF[pIndex-1] + NumberOfIndices[parsedVertex, pIndex-1];
           fieldIndexStartF[Length[fields]+1] = NumberOfIndices[parsedVertex];
           
           fieldIndexStart = Table[fieldIndexStartF[i], {i, 1, Length[fields] + 1}];
           
           indexBounds = IndexBounds[parsedVertex];
           
           prototype = ("template<> struct " <> dataClassName <> "\n" <>
                        "{\n" <>
                        IndentText @
                        ("static constexpr IndexBounds<" <> ToString @ NumberOfIndices[parsedVertex] <>
                         "> index_bounds" <>
                         If[NumberOfIndices[parsedVertex] =!= 0,
                            " = { " <>
                            "{ " <> StringJoin @ Riffle[ToString /@ indexBounds[[1]], ", "] <> " }, " <>
                            "{ " <> StringJoin @ Riffle[ToString /@ indexBounds[[2]], ", "] <> " } };\n"
                            ,
                            ";/n"
                            ] <>
                         "static constexpr std::array<unsigned, " <> ToString @ Length[fieldIndexStart] <>
                         "> fieldIndexStart = { { " <> StringJoin @ Riffle[ToString /@ fieldIndexStart, ", "] <>
                         " } };\n" <>
                         "using vertex_type = " <> VertexClassName[parsedVertex] <> ";\n"
                         ) <>
                        "};");
           
           definition = ("template<> " <> functionClassName <> "::vertex_type\n" <>
                         functionClassName <> "::vertex(const indices_type &indices, const EvaluationContext &context)\n" <>
                         "{\n" <>
                         IndentText @ VertexFunctionBody[parsedVertex] <> "\n" <>
                         "}");
           
           Return[{prototype, definition}];
          ];

(* Creates local declarations of field indices, whose values are taken
   from the elements of `arrayName'.
 *)
DeclareIndices[indexedFields_List, arrayName_String] :=
    Module[{p, total = 0, fieldIndexList, decl = "", idx},
           DeclareIndex[idx_, num_Integer, an_String] := (
               "const unsigned " <> CConversion`ToValidCSymbolString[idx] <>
               " = " <> an <> "[" <> ToString[num] <> "];\n");
           For[p = 1, p <= Length[indexedFields], p++,
               fieldIndexList = Vertices`FieldIndexList[indexedFields[[p]]];
               decl = decl <> StringJoin[DeclareIndex[#, total++, arrayName]& /@ fieldIndexList];
              ];
           Assert[total == Total[Length[Vertices`FieldIndexList[#]]& /@ indexedFields]];
           decl
          ];

GetComplexScalarCType[] :=
    CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];

(* ParsedVertex structure:
 ParsedVertex[
              {numP1Indices, numP2Indices, ...},
              {{minIndex1, minIndex2, ...}, {maxIndex1+1, maxIndex2+1, ...}},
              VertexClassName,
              VertexFunctionBody
              ]

 Getters are available! Given below ParseVertex[]
 *)

(* The heart of the algorithm! From the field content, determine all
 necessary information. *)
ParseVertex[fields_List, vertexRules_List] :=
    Module[{indexedFields, numberOfIndices, declareIndices,
        parsedVertex, vertexClassName, vertexFunctionBody,
        fieldInfo, trIndexBounds, indexBounds,
        expr, exprL, exprR},
           indexedFields = MapIndexed[(Module[{field = #1,
                                               index = #2[[1]]},
                                              SARAH`getFull[#1] /. SARAH`subGC[index] /. SARAH`subIndFinal[index,index]
                                              ] &), fields];
           indexedFields = StripLorentzIndices /@ indexedFields;
           
           
           numberOfIndices = ((Length @ Vertices`FieldIndexList[#] &) /@ indexedFields);
           declareIndices = DeclareIndices[indexedFields, "indices"];
           
           vertexClassName = SymbolName[VertexTypeForFields[fields]];
           vertexFunctionBody = Switch[vertexClassName,
                                       "SingleComponentedVertex",
                                       expr = (SARAH`Cp @@ fields) /. vertexRules;
                                       expr = TreeMasses`ReplaceDependenciesReverse[expr];
                                       declareIndices <>
                                       Parameters`CreateLocalConstRefs[expr] <> "\n" <>
                                       "const " <> GetComplexScalarCType[] <> " result = " <>
                                       Parameters`ExpressionToString[expr] <> ";\n\n" <>
                                       "return vertex_type(result);",
                                       
                                       "LeftAndRightComponentedVertex",
                                       exprL = SARAH`Cp[Sequence @@ fields][SARAH`PL] /. vertexRules;
                                       exprR = SARAH`Cp[Sequence @@ fields][SARAH`PR] /. vertexRules;
                                       exprL = TreeMasses`ReplaceDependenciesReverse[exprL];
                                       exprR = TreeMasses`ReplaceDependenciesReverse[exprR];
                                       declareIndices <>
                                       Parameters`CreateLocalConstRefs[exprL + exprR] <> "\n" <>
                                       "const " <> GetComplexScalarCType[] <> " left = " <>
                                       Parameters`ExpressionToString[exprL] <> ";\n\n" <>
                                       "const " <> GetComplexScalarCType[] <> " right = " <>
                                       Parameters`ExpressionToString[exprR] <> ";\n\n" <>
                                       "return vertex_type(left, right);"];

           fieldInfo = CleanFieldInfo /@ fields;

           trIndexBounds = Cases[Flatten[(With[{fieldIndex = #},
                                             (If[#[[1]] === SARAH`generation,
                                                 {fieldInfo[[fieldIndex, 2]]-1, fieldInfo[[fieldIndex, 3]]},
                                                 {1, #[[2]]}]
                                              &) /@ fieldInfo[[fieldIndex, 5]]]
                                        &) /@ Table[i, {i, Length[fields]}],
                                       1],
                                 Except[{}]];
           
           If[trIndexBounds === {},
              indexBounds = {{},{}},
              indexBounds = Transpose @ trIndexBounds];

           parsedVertex = ParsedVertex[numberOfIndices,
                                       indexBounds,
                                       vertexClassName,
                                       vertexFunctionBody];

           Return[parsedVertex];
           ];

(** Getters to the ParsedVertex structure **)
NumberOfIndices[parsedVertex_ParsedVertex] := Total @ parsedVertex[[1]];
NumberOfIndices[parsedVertex_ParsedVertex, pIndex_Integer] := parsedVertex[[1, pIndex]];

IndexBounds[parsedVertex_ParsedVertex] := parsedVertex[[2]];

VertexClassName[parsedVertex_ParsedVertex] := parsedVertex[[3]];
VertexFunctionBody[parsedVertex_ParsedVertex] := parsedVertex[[4]];
(** End getters **)

(* Find all contributing diagrams *)
cachedContributingDiagrams = Null;
ContributingDiagrams[] :=
       Module[{diagrams},
           If[cachedContributingDiagrams =!= Null, Return[cachedContributingDiagrams]];
           LoadVerticesIfNecessary[];

           diagrams = Flatten[(ContributingDiagramsOfType[#] &)
                                 /@ diagramTypes
                                 , 1];
              
           cachedContributingDiagrams = ({#, Union @
                   (Sequence @@ Cases[diagrams,
                         {#, diags_List} -> diags])} &) /@ edmFields;
           Return[cachedContributingDiagrams];
          ];

LoadVerticesIfNecessary[] :=
    Module[{},
           If[SARAH`VertexList3 =!= List || Length[SARAH`VertexList3] === 0,
              SA`CurrentStates = FlexibleSUSY`FSEigenstates;
              SARAH`InitVertexCalculation[FlexibleSUSY`FSEigenstates, False];
              SARAH`ReadVertexList[FlexibleSUSY`FSEigenstates, False, False, True];
              SARAH`MakeCouplingLists;
              ];
           ];

ConcreteDiagramEvaluators[] :=
     ({#[[1]],
         (("DiagramEvaluator<" <> SymbolName @ Head @ #[[1]] <> "<" <>
           ToString @ #[[1,1]] <> ">, " <>
           StringJoin @ (Riffle[CXXNameOfField /@ List @@ #[[2;;]], ", "]) <>
           ">" &)
          /@ #[[2]]) } &) /@ ContributingDiagrams[];

End[];

EndPackage[];

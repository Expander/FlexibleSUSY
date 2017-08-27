(* ::Package:: *)

BeginPackage["CXXDiagrams`", {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "Parameters`"}];

(* This module generates c++ code intended to be used similarly to SARAH's fields and Vertex[] function *)

Initialize[] := LoadVerticesIfNecessary[]

LoadVerticesIfNecessary[] :=
   If[SARAH`VertexList3 =!= List || Length[SARAH`VertexList3] === 0,
        SA`CurrentStates = FlexibleSUSY`FSEigenstates; 
        SARAH`InitVertexCalculation[FlexibleSUSY`FSEigenstates, False];
        SARAH`ReadVertexList[FlexibleSUSY`FSEigenstates, False, False, True];
        SARAH`MakeCouplingLists;
   ]

(* The supported vertex types.
 They have the same names as their c++ counterparts. *)
vertexTypes = {
    SingleComponentedVertex,
    LeftAndRightComponentedVertex
};

(* Return a string corresponding to the c++ class name of the field.
 Note that "bar" and "conj" get turned into bar<...>::type and
 conj<...>::type respectively! *)
CXXNameOfField[p_] := SymbolName[p];
CXXNameOfField[SARAH`bar[p_]] := "bar<" <> SymbolName[p] <> ">::type";
CXXNameOfField[Susyno`LieGroups`conj[p_]] := "conj<" <> SymbolName[p] <> ">::type";

(* TODO: Better name than LorentzConjugate *)
LorentzConjugateOperation[field_] := If[FermionQ[field] || GhostQ[field],
                                        "bar",
                                        "conj"];
LorentzConjugate[field_] := SARAH`AntiField[field]

NumberOfFieldIndices[field_] := Length @ FieldInfo[field][[5]]

CreateFields[] :=
  Module[{fields},
       fields = TreeMasses`GetParticles[];
       
       StringJoin @ Riffle[
         ("struct " <> CXXNameOfField[#] <> ": public Field {\n" <>
            TextFormatting`IndentText["static constexpr unsigned numberOfGenerations = " <>
                                         ToString @ TreeMasses`GetDimension[#] <> ";\n" <>
                                      "static constexpr unsigned numberOfFieldIndices = " <>
                                         ToString @ NumberOfFieldIndices[#] <> ";\n"] <>
            "};\n" &) /@ fields, "\n"] <> "\n\n" <>
              
       "// Special field families\n" <>
       "using Photon = " <> CXXNameOfField @ GetPhoton[] <> ";\n\n" <>
       
       "// Fields that contract with themselves to form Lorentz scalars.\n" <>
       StringJoin @ Riffle[
         ("template<> struct " <> LorentzConjugateOperation[#] <> "<" <> CXXNameOfField[#] <> ">" <>
            " { using type = " <> CXXNameOfField[#] <> "; };"
            &) /@ Select[fields, (# == LorentzContractionPartner[#] &)],
          "\n"] <> "\n\n" <>
       
       "// Declare a type that can hold the field indices for any given field\n" <>
       "template<class Field> struct field_indices\n" <>
       "{\n" <>
          TextFormatting`IndentText[
            "using type = std::array<unsigned, Field::numberOfFieldIndices>;\n"] <>
       "};\n"
  ]

(* adjacencyMatrix must be undirected (i.e symmetric) *)
FeynmanDiagramsOfType[adjacencyMatrix_List,externalFields_List] :=
  Module[{externalVertices = externalFields[[All,1]],
          internalVertices,externalRules,
          internalFieldCouplings,
          unspecifiedEdgesLess,unspecifiedEdgesEqual,
          insertFieldRulesLess,insertFieldRulesGreater,insertFieldRulesEqual,
          fieldsToInsert,
          unresolvedFieldCouplings,resolvedFields,resolvedFieldCouplings},
   internalVertices = Complement[Table[k,{k,Length[adjacencyMatrix]}],externalVertices];
   externalRules = Flatten @ ({{_,#,_} :> SARAH`AntiField[# /. externalFields],
                               {#,_,_} :> SARAH`AntiField[# /. externalFields]} & /@ externalVertices);

   internalFieldCouplings = (Flatten[(Flatten @ Position[adjacencyMatrix[[#]],Except[0],{1},Heads -> False]
                                /. {i_Integer :> Table[{#,i,k},{k,adjacencyMatrix[[#,i]]}]}),1] &
                             /@ internalVertices) /. externalRules;

   unspecifiedEdgesLess = Cases[internalFieldCouplings,{i_,j_,_} /; i < j,{2}];
   unspecifiedEdgesEqual = Cases[internalFieldCouplings,{i_,i_,_},{2}];

   insertFieldRulesLess = MapIndexed[#1 -> SARAH`FieldToInsert[#2[[1]]] &,unspecifiedEdgesLess];
   insertFieldRulesGreater = (insertFieldRulesLess /. {Rule[{i_,j_,k_},field_] :> Rule[{j,i,k},SARAH`AntiField[field]]});
   insertFieldRulesEqual = MapIndexed[#1 -> {SARAH`FieldToInsert[#2[[1]]+Length[insertFieldRulesLess]],
                                            SARAH`AntiField[SARAH`FieldToInsert[#2[[1]]+Length[insertFieldRulesLess]]]} &,
                                      unspecifiedEdgesEqual];
   fieldsToInsert = Table[SARAH`FieldToInsert[k],
             {k,Length[insertFieldRulesLess] + Length[insertFieldRulesEqual]}];
   
   unresolvedFieldCouplings = internalFieldCouplings
     /. insertFieldRulesLess /. insertFieldRulesGreater /. insertFieldRulesEqual;
   resolvedFields = SARAH`InsFields[{C @@@ unresolvedFieldCouplings,
                                     fieldsToInsert}][[All,2]];
   resolvedFieldCouplings = unresolvedFieldCouplings /.
     ((Rule @@@ Transpose[{fieldsToInsert,#}]) & /@ resolvedFields);
   
   Table[k,{k,Length[adjacencyMatrix]}] /. externalFields /. 
     ((Rule @@@ Transpose[{internalVertices,#}]) & /@ resolvedFieldCouplings)
  ]

VerticesForDiagram[diagram_] := Select[diagram,Length[#] > 1 &]

CreateVertexData[vertices_List,vertexRules_List] := 
  Module[{dataClassName,indexBounds,parsedVertex,fieldIndexStartF,fieldIndexStart},
    parsedVertex = ParseVertex[fields, vertexRules];
    indexBounds = IndexBounds[parsedVertex];
    
    fieldIndexStartF[1] = 0;
    fieldIndexStartF[pIndex_] := fieldIndexStartF[pIndex-1] + NumberOfIndices[parsedVertex, pIndex-1];
    fieldIndexStartF[Length[fields]+1] = NumberOfIndices[parsedVertex];
           
    fieldIndexStart = Table[fieldIndexStartF[i], {i, 1, Length[fields] + 1}];
    
    dataClassName = "VertexFunctionData<" <> StringJoin @ Riffle[CXXNameOfField /@ fields, ", "] <> ">";
    
    "template<> struct " <> dataClassName <> "\n" <>
    "{\n" <>
    TextFormatting`IndentText[
      "static constexpr IndexBounds<" <> ToString @ NumberOfIndices[parsedVertex] <>
        "> index_bounds" <>
        If[NumberOfIndices[parsedVertex] =!= 0,
           " = { " <>
             "{ " <> StringJoin @ Riffle[ToString /@ indexBounds[[1]], ", "] <> " }, " <>
             "{ " <> StringJoin @ Riffle[ToString /@ indexBounds[[2]], ", "] <> " } }",
           "{}"
        ] <> ";\n" <>
      "static constexpr unsigned fieldIndexStart[" <> ToString @ Length[fieldIndexStart] <>
         "] = { " <> StringJoin @ Riffle[ToString /@ fieldIndexStart, ", "] <> " };\n" <>
      "using vertex_type = " <> VertexClassName[parsedVertex] <> ";"] <> "\n" <>
    "};"
  ]

(* Returns the necessary c++ code corresponding to the vertices that need to be calculated.
 The returned value is a list {prototypes, definitions}. *)
CreateVertexFunctions[vertices_List,vertexRules_List] :=
  StringJoin @\[NonBreakingSpace]Riffle[CreateVertexFunction[#, vertexRules] & /@ vertices,
                      "\n\n"]

(* Creates the actual c++ code for a vertex with given fields.
 You should never need to change this code! *)
CreateVertexFunction[fields_List, vertexRules_List] :=
  Module[{parsedVertex, functionClassName},
         parsedVertex = ParseVertex[fields, vertexRules];
         functionClassName = "VertexFunction<" <> StringJoin @ Riffle[CXXNameOfField /@ fields, ", "] <> ">";
         
         "template<> " <> functionClassName <> "::vertex_type\n" <>
         functionClassName <> "::vertex(const indices_type &indices, const EvaluationContext &context)\n" <>
         "{\n" <>
         TextFormatting`IndentText @ VertexFunctionBody[parsedVertex] <> "\n" <>
         "}"
  ]

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
           indexedFields = IndexFields[fields];
           
           numberOfIndices = ((Length @ Vertices`FieldIndexList[#] &) /@ indexedFields);
           declareIndices = DeclareIndices[indexedFields, "indices"];
           
           vertexClassName = SymbolName[VertexTypeForFields[fields]];
           vertexFunctionBody = Switch[vertexClassName,
                                       "SingleComponentedVertex",
                                       expr = Vertices`SortCp[SARAH`Cp @@ indexedFields] /. vertexRules;
                                       expr = TreeMasses`ReplaceDependenciesReverse[expr];
                                       declareIndices <>
                                       Parameters`CreateLocalConstRefs[expr] <> "\n" <>
                                       "const " <> GetComplexScalarCType[] <> " result = " <>
                                       Parameters`ExpressionToString[expr] <> ";\n\n" <>
                                       "return vertex_type(result);",
                                       
                                       "LeftAndRightComponentedVertex",
                                       exprL = Vertices`SortCp @ SARAH`Cp[Sequence @@ indexedFields][SARAH`PL] /. vertexRules;
                                       exprR = Vertices`SortCp @ SARAH`Cp[Sequence @@ indexedFields][SARAH`PR] /. vertexRules;
                                       exprL = TreeMasses`ReplaceDependenciesReverse[exprL];
                                       exprR = TreeMasses`ReplaceDependenciesReverse[exprR];
                                       declareIndices <>
                                       Parameters`CreateLocalConstRefs[exprL + exprR] <> "\n" <>
                                       "const " <> GetComplexScalarCType[] <> " left = " <>
                                       Parameters`ExpressionToString[exprL] <> ";\n\n" <>
                                       "const " <> GetComplexScalarCType[] <> " right = " <>
                                       Parameters`ExpressionToString[exprR] <> ";\n\n" <>
                                       "return vertex_type(left, right);"];

           fieldInfo = FieldInfo /@ fields;

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

           parsedVertex
           ];

IndexFields[fields_List] :=
    MapIndexed[
	Module[{field = #1,
		index = #2[[1]]},
	       StripLorentzIndices[
		   SARAH`getFull[field] /. SARAH`subGC[index] /.
		   SARAH`subIndFinal[index,index]]
               ] &, fields];

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

(** Getters to the ParsedVertex structure **)
NumberOfIndices[parsedVertex_ParsedVertex] := Total @ parsedVertex[[1]];
NumberOfIndices[parsedVertex_ParsedVertex, pIndex_Integer] := parsedVertex[[1, pIndex]];

IndexBounds[parsedVertex_ParsedVertex] := parsedVertex[[2]];

VertexClassName[parsedVertex_ParsedVertex] := parsedVertex[[3]];
VertexFunctionBody[parsedVertex_ParsedVertex] := parsedVertex[[4]];
(** End getters **)

(* Returns the vertex type for a vertex with a given list of fields *)
VertexTypeForFields[fields_List] :=
  Module[{fermions, scalarCount, vectorCount, fermionCount, vertexType = "UnknownVertexType"},
    fermions = Select[fields, TreeMasses`IsFermion];
      
    scalarCount = Length @ Select[fields, TreeMasses`IsScalar];
    vectorCount = Length @ Select[fields, TreeMasses`IsVector];
    fermionCount = Length @ fermions;
       
    If[fermionCount === 2 && scalarCount === 1 && vectorCount === 0,
       vertexType = LeftAndRightComponentedVertex];
    If[fermionCount === 2 && scalarCount === 0 && vectorCount === 1,
       If[fermions[[1]] === SARAH`AntiField[fermions[[2]]],
          vertexType = LeftAndRightComponentedVertex]];
    If[fermionCount === 0 && scalarCount === 2 && vectorCount === 1,
       vertexType = SingleComponentedVertex];
    vertexType
  ]

(* Returns the different SARAH`Cp coupling parts for a vertex with a given list of fields *)
CouplingsForFields[fields_List] :=
    Module[{vertexType, couplings},
      vertexType = VertexTypeForFields[fields];
      couplings = {SARAH`Cp @@ fields};
      
      Switch[vertexType,
             SingleComponentedVertex, couplings,
             LeftAndRightComponentedVertex, {couplings[[1]][SARAH`PL], couplings[[1]][SARAH`PR]},
             "UnknownVertexType",{}]
   ]

IsLorentzIndex[index_] := StringMatchQ[ToString @ index, "lt" ~~ __];

StripLorentzIndices[p_Symbol] := p;
StripLorentzIndices[SARAH`bar[p_]] := SARAH`bar[StripLorentzIndices[p]];
StripLorentzIndices[Susyno`LieGroups`conj[p_]] := Susyno`LieGroups`conj[StripLorentzIndices[p]];
StripLorentzIndices[p_] := Module[{remainingIndices},
                                  remainingIndices = Select[p[[1]], (!IsLorentzIndex[#] &)];
                                  If[Length[remainingIndices] === 0, Head[p],
                                     Head[p][remainingIndices]]
                                  ];

FieldInfo[field_,OptionsPattern[{includeLorentzIndices -> False}]] := 
    Module[{fieldInfo = Cases[SARAH`Particles[FlexibleSUSY`FSEigenstates],
                                {SARAH`getParticleName @ field, ___}][[1]]},
            fieldInfo = DeleteCases[fieldInfo, {SARAH`generation, 1}, {2}];
            If[!OptionValue[includeLorentzIndices],
               DeleteCases[fieldInfo, {SARAH`lorentz, _}, {2}],
               fieldInfo]
           ];

NPointFunctions[vertices_List] :=
  Flatten[(Null[Null, #] &) /@ ((CouplingsForFields[#] &) /@ vertices)]

EndPackage[];

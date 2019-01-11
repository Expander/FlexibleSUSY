(* ::Package:: *)

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

BeginPackage["CXXDiagrams`", {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "Parameters`","CConversion`"}];

(* This module generates c++ code intended to be used similarly to SARAH's fields and Vertex[] function *)

(* Vertex types *)
ScalarVertex::usage="";
ChiralVertex::usage="";
MomentumDifferenceVertex::usage="";
InverseMetricVertex::usage="";

VertexTypes::usage="";
CXXNameOfField::usage="";
LorentzConjugateOperation::usage="";
LorentzConjugate::usage="";
RemoveLorentzConjugation::usage="";
AtomHead::usage="";
CreateFields::usage="";
FeynmanDiagramsOfType::usage="";
VerticesForDiagram::usage="";
ContractionsBetweenVerticesForDiagramFromGraph::usage="";
CreateVertexData::usage="";
CreateVertices::usage="";
VertexRulesForVertices::usage="";
CreateMassFunctions::usage="";
CreateUnitCharge::usage="";
NumberOfFieldIndices::usage="";
FieldInfo::usage="";
includeLorentzIndices::usage="";

Begin["`Private`"];

VertexTypes[] := {
    ScalarVertex,
    ChiralVertex,
    MomentumDifferenceVertex,
    InverseMetricVertex
};

(* Return a string corresponding to the c++ class name of the field.
 Note that "bar" and "conj" get turned into bar<...>::type and
 conj<...>::type respectively! *)
CXXNameOfField[p_, OptionsPattern[{prefixNamespace -> False}]] :=
  If[StringQ[OptionValue[prefixNamespace]],
     OptionValue[prefixNamespace] <> "::",
     ""] <> SymbolName[p];
CXXNameOfField[SARAH`bar[p_], OptionsPattern[{prefixNamespace -> False}]] :=
  "typename bar<" <> CXXNameOfField[p, prefixNamespace -> OptionValue[prefixNamespace]] <>
  ">::type";
CXXNameOfField[Susyno`LieGroups`conj[p_],
               OptionsPattern[{prefixNamespace -> False}]] :=
  "typename conj<" <> CXXNameOfField[p, prefixNamespace -> OptionValue[prefixNamespace]] <>
  ">::type";

CXXBoolValue[True] = "true"
CXXBoolValue[False] = "false"

(* TODO: Better name than LorentzConjugate *)
(* FIXME: We decide this on our own which is bad, but 
	SARAH`AntiField[] only works after one has called some
	SARAH routines like in LoadVerticesIfNecessary[].
	But CXXDiagrams should be stateless, hence we avoid relying
	SARAH here. *)
LorentzConjugateOperation[field_] :=
	If[TreeMasses`IsFermion[field] || TreeMasses`IsGhost[field],
		"bar", "conj"]
LorentzConjugate[field_] :=
	If[TreeMasses`IsFermion[field] || TreeMasses`IsGhost[field],
		SARAH`bar[field], Susyno`LieGroups`conj[field]]

RemoveLorentzConjugation[p_] := p
RemoveLorentzConjugation[SARAH`bar[p_]] := p
RemoveLorentzConjugation[Susyno`LieGroups`conj[p_]] := p

AtomHead[x_] := If[AtomQ[x], x, AtomHead[Head[x]]]

CreateFields[] :=
  Module[{fields, scalars, fermions, vectors, ghosts},
       fields = TreeMasses`GetParticles[];
       scalars = Select[fields, TreeMasses`IsScalar];
       fermions = Select[fields, TreeMasses`IsFermion];
       vectors = Select[fields, TreeMasses`IsVector];
       ghosts = Select[fields, TreeMasses`IsGhost];
       
       StringJoin @ Riffle[
         ("struct " <> CXXNameOfField[#] <> " {\n" <>
            TextFormatting`IndentText[
              "using index_bounds = boost::mpl::pair<\n" <>
              "  boost::mpl::vector_c<int" <>
                   StringJoin[", " <> ToString[#] & /@ 
                     (IndexBoundsForField[#][[1]] - 1)] <> ">,\n" <>
              "  boost::mpl::vector_c<int" <>
                   StringJoin[", " <> ToString[#] & /@
                     IndexBoundsForField[#][[2]]] <> ">\n" <>
              ">;\n" <>
              "static constexpr int numberOfGenerations = " <>
                   ToString @ TreeMasses`GetDimension[#] <> ";\n" <>
                     "using sm_flags = boost::mpl::vector_c<bool, " <>
                        If[TreeMasses`GetDimension[#] === 1,
                           CXXBoolValue @ TreeMasses`IsSMParticle[#],
                           StringJoin @ Riffle[CXXBoolValue /@
                             (TreeMasses`IsSMParticle[#] & /@ Table[#[{k}],{k,TreeMasses`GetDimension[#]}]),
                                               ", "]] <>
                        ">;\n" <>
                "static constexpr int numberOfFieldIndices = " <>
                   ToString @ NumberOfFieldIndices[#] <> ";\n" <>
                "static constexpr double electric_charge = " <>
                   CConversion`RValueToCFormString[TreeMasses`GetElectricCharge[#]] <> ";\n" <>
                "using lorentz_conjugate = " <>
                   CXXNameOfField[LorentzConjugate[#]] <> ";\n"] <>
              "};" &) /@ fields, "\n\n"] <> "\n\n" <>
       
       "// Named fields\n" <>
       "using Photon = " <> CXXNameOfField[SARAH`Photon] <> ";\n" <>
       "using Electron = " <> CXXNameOfField[AtomHead @ TreeMasses`GetSMElectronLepton[]] <> ";\n\n" <>
       
       "// Fields that are their own Lorentz conjugates.\n" <>
       StringJoin @ Riffle[
         ("template<> struct " <> LorentzConjugateOperation[#] <> "<" <> CXXNameOfField[#] <> ">" <>
            " { using type = " <> CXXNameOfField[#] <> "; };"
            &) /@ Select[fields, (# == LorentzConjugate[#] &)],
          "\n"] <> "\n\n" <>

       "using scalars = boost::mpl::vector<" <>
         StringJoin[Riffle[CXXNameOfField /@ scalars, ", "]] <> ">;\n" <>
       "using fermions = boost::mpl::vector<" <>
         StringJoin[Riffle[CXXNameOfField /@ fermions, ", "]] <> ">;\n" <>
       "using vectors = boost::mpl::vector<" <>
         StringJoin[Riffle[CXXNameOfField /@ vectors, ", "]] <> ">;\n" <>
       "using ghosts = boost::mpl::vector<" <>
         StringJoin[Riffle[CXXNameOfField /@ ghosts, ", "]] <> ">;"
  ]

(* adjacencyMatrix must be undirected (i.e symmetric) *)
FeynmanDiagramsOfType[adjacencyMatrix_List,externalFields_List] :=
  Module[{externalVertices = externalFields[[All,1]],
          internalVertices,externalRules,
          internalFieldCouplings,
          unspecifiedEdgesLess,unspecifiedEdgesEqual,
          insertFieldRulesLess,insertFieldRulesGreater,insertFieldRulesEqual,
          fieldsToInsert,
          unresolvedFieldCouplings,resolvedFields,resolvedFieldCouplings,
          diagrams},
   LoadVerticesIfNecessary[];
   
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
   
   diagrams = Table[k,{k,Length[adjacencyMatrix]}] /. externalFields /. 
     ((Rule @@@ Transpose[{internalVertices,#}]) & /@ resolvedFieldCouplings);
   
   DeleteDuplicates[diagrams,
     AllTrue[Cases[Transpose[{#1,#2}],{{___},{___}}], (* External lines *)
             (Sort[#[[1]]] === Sort[#[[2]]]&)] &]
  ]

VerticesForDiagram[diagram_] := Select[diagram,Length[#] > 1 &]

ContractionsBetweenVerticesForDiagramFromGraph[v1_Integer, v2_Integer,
		diagram_List, graph_List] :=
	Module[{fields1 = diagram[[v1]], fields2 = diagram[[v2]],
			preceedingNumberOfFields1 = Total[graph[[v1, ;;v2]]] - graph[[v1,v2]],
			preceedingNumberOfFields2 = Total[graph[[v2, ;;v1]]] - graph[[v2,v1]],
			contractedFieldIndices1, contractedFieldIndices2},
		contractedFieldIndices1 = Table[k, {k, preceedingNumberOfFields1 + 1,
			preceedingNumberOfFields1 + graph[[v1,v2]]}];
		contractedFieldIndices2 = Table[k, {k, preceedingNumberOfFields2 + 1,
			preceedingNumberOfFields2 + graph[[v2,v1]]}];
		
		Transpose[{contractedFieldIndices1, contractedFieldIndices2}]
	]

CreateVertexData[fields_List] := 
  Module[{dataClassName},
    dataClassName = "VertexData<" <> StringJoin[Riffle[
      CXXNameOfField[#, prefixNamespace -> "fields"] & /@ fields,
    ", "]] <> ">";
    
    "template<> struct " <> dataClassName <> "\n" <>
    "{\n" <>
    TextFormatting`IndentText[
      "using vertex_type = " <> SymbolName[VertexTypeForFields[fields]] <>
         ";"] <> "\n" <>
    "};"
  ]

(* Returns the necessary c++ code corresponding to the vertices that need to be calculated.
 The returned value is a list {prototypes, definitions}. *)
CreateVertices[vertices_List, OptionsPattern[{StripColorStructure -> False}]] :=
  StringJoin @ Riffle[
		CreateVertex[#, StripColorStructure -> OptionValue[StripColorStructure]] & /@
			DeleteDuplicates[vertices], "\n\n"]

(* Creates the actual c++ code for a vertex with given fields.
 You should never need to change this code! *)
CreateVertex[fields_List, OptionsPattern[{StripColorStructure -> False}]] :=
  Module[{functionClassName},
    LoadVerticesIfNecessary[];
    functionClassName = "Vertex<" <> StringJoin @ Riffle[
    CXXNameOfField[#, prefixNamespace -> "fields"] & /@ fields, ", "] <> ">";

    "template<> inline\n" <> 
    functionClassName <> "::vertex_type\n" <>
    functionClassName <> "::evaluate(const indices_type& indices, const context_base& context)\n" <>
    "{\n" <>
    TextFormatting`IndentText @ VertexFunctionBodyForFields[fields,
			StripColorStructure -> OptionValue[StripColorStructure]] <> "\n" <>
    "}"
  ]

VertexFunctionBodyForFields[fields_List, OptionsPattern[{StripColorStructure -> False}]] := 
	Switch[Length[fields],
		3, VertexFunctionBodyForFieldsImpl[fields, SARAH`VertexList3,
			StripColorStructure -> OptionValue[StripColorStructure]],
		4, VertexFunctionBodyForFieldsImpl[fields, SARAH`VertexList4,
			StripColorStructure -> OptionValue[StripColorStructure]],
		_, "non-(3,4)-point vertex"]

VertexFunctionBodyForFieldsImpl[fields_List, vertexList_List,
		OptionsPattern[{StripColorStructure -> False}]] :=
  Module[{sortedFields, sortedIndexedFields, indexedFields,
          fieldsOrdering, sortedFieldsOrdering, inverseFOrdering,
          fOrderingWRTSortedF, vertex, vExpression, vertexIsZero = False,
          vertexType = VertexTypeForFields[fields], expr, exprL, exprR,
          vertexRules, incomingScalar, outgoingScalar},
    sortedFields = Vertices`SortFieldsInCp[fields];
    
    vertex = Select[vertexList, StripFieldIndices[#[[1]]] === sortedFields &, 1];
    If[vertex === {}, vertexIsZero = True, vertex = vertex[[1]]];

    If[vertexIsZero,
       Return[Switch[vertexType,
         ScalarVertex,
         "return vertex_type(0);",
         
         ChiralVertex,
         "return vertex_type(0, 0);",
         
         MomentumDifferenceVertex,
         "return vertex_type(0, " <> StringJoin[Riffle[
            ToString /@ Flatten[Position[fields,
               field_ /; TreeMasses`IsScalar[field] || TreeMasses`IsGhost[field],
               {1}, Heads -> False] - 1],
            ", "]] <> ");",
         
         InverseMetricVertex,
         "return vertex_type(0);"]]];

    sortedIndexedFields = vertex[[1]];
    
    (* Mathematica 7 does not know about permutations... :'-( *)
    fieldsOrdering = Ordering[fields];
    sortedFieldsOrdering = Ordering[sortedFields];

    inverseFOrdering = Ordering[fieldsOrdering];
    fOrderingWRTSortedF = sortedFieldsOrdering[[inverseFOrdering]];

    indexedFields = sortedIndexedFields[[fOrderingWRTSortedF]];

    Switch[vertexType,
      ScalarVertex,
      vertexRules = {(SARAH`Cp @@ sortedIndexedFields) ->
        Vertices`FindVertexWithLorentzStructure[Rest[vertex], 1][[1]]};
      
      expr = CanonicalizeCoupling[SARAH`Cp @@ fields,
        sortedFields, sortedIndexedFields] /. vertexRules;
      
      expr = Vertices`SarahToFSVertexConventions[sortedFields, expr,
				StripColorStructure -> OptionValue[StripColorStructure]];
      expr = TreeMasses`ReplaceDependenciesReverse[expr];
      DeclareIndices[StripLorentzIndices /@ indexedFields, "indices"] <>
      Parameters`CreateLocalConstRefs[expr] <> "\n" <>
      "const " <> GetComplexScalarCType[] <> " result = " <>
      Parameters`ExpressionToString[expr] <> ";\n\n" <>
      "return vertex_type(result);",
      
      ChiralVertex,
      vertexRules = {
        (SARAH`Cp @@ sortedIndexedFields)[SARAH`PL] ->
          Vertices`FindVertexWithLorentzStructure[Rest[vertex], SARAH`PL][[1]],
        (SARAH`Cp @@ sortedIndexedFields)[SARAH`PR] ->
          Vertices`FindVertexWithLorentzStructure[Rest[vertex], SARAH`PR][[1]]};

      exprL = CanonicalizeCoupling[(SARAH`Cp @@ fields)[SARAH`PL],
        sortedFields, sortedIndexedFields] /. vertexRules;
      exprR = CanonicalizeCoupling[(SARAH`Cp @@ fields)[SARAH`PR],
        sortedFields, sortedIndexedFields] /. vertexRules;

      exprL = Vertices`SarahToFSVertexConventions[sortedFields, exprL,
				StripColorStructure -> OptionValue[StripColorStructure]];
      exprR = Vertices`SarahToFSVertexConventions[sortedFields, exprR,
				StripColorStructure -> OptionValue[StripColorStructure]];
      exprL = TreeMasses`ReplaceDependenciesReverse[exprL];
      exprR = TreeMasses`ReplaceDependenciesReverse[exprR];
      DeclareIndices[StripLorentzIndices /@ indexedFields, "indices"] <>
      Parameters`CreateLocalConstRefs[{exprL, exprR}] <> "\n" <>
      "const " <> GetComplexScalarCType[] <> " left = " <>
      Parameters`ExpressionToString[exprL] <> ";\n\n" <>
      "const " <> GetComplexScalarCType[] <> " right = " <>
      Parameters`ExpressionToString[exprR] <> ";\n\n" <>
      "return vertex_type(left, right);",

      MomentumDifferenceVertex,
      {incomingScalar, outgoingScalar} = Replace[vertex[[2,2]],
        SARAH`Mom[is_,_] - SARAH`Mom[os_,_] :> {is, os}];
      vertexRules = {(SARAH`Cp @@ sortedIndexedFields) -> vertex[[2,1]]};
      
      expr = CanonicalizeCoupling[SARAH`Cp @@ fields,
        sortedFields, sortedIndexedFields] /. vertexRules;
      
      expr = Vertices`SarahToFSVertexConventions[sortedFields, expr,
				StripColorStructure -> OptionValue[StripColorStructure]];
      expr = TreeMasses`ReplaceDependenciesReverse[expr];
      "int minuend_index = " <> 
        ToString[Position[indexedFields, incomingScalar][[1,1]] - 1] <> ";\n" <>
      "int subtrahend_index = " <>
        ToString[Position[indexedFields, outgoingScalar][[1,1]] - 1] <> ";\n\n" <>
      DeclareIndices[StripLorentzIndices /@ indexedFields, "indices"] <>
      Parameters`CreateLocalConstRefs[expr] <> "\n" <>
      "const " <> GetComplexScalarCType[] <> " result = " <>
      Parameters`ExpressionToString[expr] <> ";\n\n" <>
      "return vertex_type(result, minuend_index, subtrahend_index);",

      InverseMetricVertex,
      vertexRules = {(SARAH`Cp @@ sortedIndexedFields) -> vertex[[2,1]]};
      
      expr = CanonicalizeCoupling[SARAH`Cp @@ fields,
        sortedFields, sortedIndexedFields] /. vertexRules;
      
      expr = Vertices`SarahToFSVertexConventions[sortedFields, expr,
				StripColorStructure -> OptionValue[StripColorStructure]];
      expr = TreeMasses`ReplaceDependenciesReverse[expr];
      DeclareIndices[StripLorentzIndices /@ indexedFields, "indices"] <>
      Parameters`CreateLocalConstRefs[expr] <> "\n" <>
      "const " <> GetComplexScalarCType[] <> " result = " <>
      Parameters`ExpressionToString[expr] <> ";\n\n" <>
      "return vertex_type(result);"
    ]
  ]

DeclareIndices[indexedFields_List, arrayName_String] :=
    Module[{p, total = 0, fieldIndexList, decl = ""},
           DeclareIndex[idx_, num_Integer, an_String] := (
               "const int " <> CConversion`ToValidCSymbolString[idx] <>
               " = " <> an <> "[" <> ToString[num] <> "];\n");
           For[p = 1, p <= Length[indexedFields], p++,
               fieldIndexList = Vertices`FieldIndexList[indexedFields[[p]]];
               decl = decl <> StringJoin[DeclareIndex[#, total++, arrayName]& /@ fieldIndexList];
              ];
           Assert[total == Total[Length[Vertices`FieldIndexList[#]]& /@ indexedFields]];
           decl
          ]

GetComplexScalarCType[] :=
    CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]]

(* Get a mathematical expression of the requested vertex
   in terms of its canonically ordered vertex. *)
CanonicalizeCoupling[
    coupling_, sortedFields_List, sortedIndexedFields_List] :=
  Module[{expr},
    (* Note: Cannot use indexed fields, because their internal
    SARAH ordering is different... *)
    expr = Vertices`SortCp[coupling];
    
    (* Index the fields *)
    expr = expr /. {SARAH`Cp @@ sortedFields -> SARAH`Cp @@ sortedIndexedFields}
  ]
  
CreateMassFunctions[] :=
  Module[{massiveFields,
          ghostMappings = SelfEnergies`ReplaceGhosts[FlexibleSUSY`FSEigenstates]},
    massiveFields = TreeMasses`GetParticles[];

    StringJoin @ Riffle[
      Module[{fieldInfo = FieldInfo[#], numberOfIndices},
             numberOfIndices = Length @ fieldInfo[[5]];
                                   
             "template<> inline\n" <>
             "double context_base::mass_impl<" <>
               CXXNameOfField[#, prefixNamespace -> "fields"] <>
             ">(const std::array<int, " <> ToString @ numberOfIndices <>
             ">& indices) const\n" <>
             "{ return model.get_M" <> CXXNameOfField[# /. ghostMappings] <>
             If[TreeMasses`GetDimension[#] === 1, "()", "(indices[0])"] <> "; }"
            ] & /@ massiveFields, "\n\n"]
        ]

CreateUnitCharge[] :=
  Module[{electron,photon,vertex,vertexBody,
          numberOfElectronIndices,numberOfPhotonIndices},
         electron = AtomHead @ TreeMasses`GetSMElectronLepton[];
         photon = SARAH`Photon;
         vertex = {SARAH`bar[electron], electron, photon};
         vertexBody = VertexFunctionBodyForFields[vertex];
         numberOfElectronIndices = NumberOfFieldIndices[electron];
         numberOfPhotonIndices = NumberOfFieldIndices[photon];

         "static ChiralVertex unit_charge(const context_base& context)\n" <>
         "{\n" <>
         TextFormatting`IndentText["using vertex_type = ChiralVertex;"] <> "\n\n" <>
         TextFormatting`IndentText @ 
           ("std::array<int, " <> ToString @ numberOfElectronIndices <> "> electron_indices = {" <>
              If[TreeMasses`GetDimension[electron] =!= 1,
                 " " <> ToString @ (FieldInfo[electron][[2]]-1) <> (* Electron has the lowest index *)
                 If[numberOfElectronIndices =!= 1,
                    StringJoin @ Table[", 0", {numberOfElectronIndices-1}],
                    ""] <> " ",
                 If[numberOfElectronIndices =!= 0,
                    StringJoin @ Riffle[Table[" 0", {numberOfElectronIndices}], ","] <> " ",
                    ""]
                ] <>
            "};\n") <>
         TextFormatting`IndentText @ 
           ("std::array<int, " <> ToString @ numberOfPhotonIndices <> "> photon_indices = {" <>
               If[TreeMasses`GetDimension[photon] =!= 1,
                 " " <> ToString @ (FieldInfo[photon][[2]]-1) <>
                 If[numberOfPhotonIndices =!= 1,
                    StringJoin @ Table[", 0", {numberOfPhotonIndices-1}],
                    ""] <> " ",
                 If[numberOfPhotonIndices =!= 0,
                    StringJoin @ Riffle[Table[" 0", {numberOfPhotonIndices}], ","] <> " ",
                    ""]
                ] <>
            "};\n") <>
         TextFormatting`IndentText @ 
           ("std::array<int, " <> ToString[
              Total[NumberOfFieldIndices /@ {photon,electron,electron}]] <>
            "> indices = " <>
              "concatenate(photon_indices, electron_indices, electron_indices);\n\n") <>
           
         TextFormatting`IndentText @ vertexBody <> "\n" <>
         "}"
  ]

NumberOfFieldIndices[field_] := Length @ FieldInfo[field][[5]]

IndexBoundsForField[field_] :=
  Module[{fieldInfo = FieldInfo[field]},
    If[NumberOfFieldIndices[field] === 0,
       Return[{{},{}}]];
    If[Length @ Cases[fieldInfo[[5]],{SARAH`generation,_}] === 0,
       Transpose[{1,#[[2]]} & /@ fieldInfo[[5]]],
       Transpose[Prepend[
         {1,#[[2]]} & /@ DeleteCases[fieldInfo[[5]],{SARAH`generation,_}],
         {fieldInfo[[2]],fieldInfo[[3]]}]]]]

LoadVerticesIfNecessary[] :=
   If[Head[SARAH`VertexList3] === Symbol || Length[SARAH`VertexList3] === 0,
        SA`CurrentStates = FlexibleSUSY`FSEigenstates; 
        SARAH`InitVertexCalculation[FlexibleSUSY`FSEigenstates, False];
        SARAH`partDefinition = ParticleDefinitions[FlexibleSUSY`FSEigenstates];
        SARAH`Particles[SARAH`Current] = SARAH`Particles[FlexibleSUSY`FSEigenstates];
        SARAH`ReadVertexList[FlexibleSUSY`FSEigenstates, False, False, True];
        SARAH`MakeCouplingLists;
   ]

(* Returns the vertex type for a vertex with a given list of fields *)
VertexTypeForFields[fields_List] :=
  Module[{scalarCount, vectorCount, fermionCount, ghostCount},
    scalarCount = Length @ Select[fields, TreeMasses`IsScalar];
    vectorCount = Length @ Select[fields, TreeMasses`IsVector];
    fermionCount = Length @ Select[fields, TreeMasses`IsFermion];
    ghostCount = Length @ Select[fields, TreeMasses`IsGhost];
    
    Switch[{fermionCount, scalarCount, vectorCount, ghostCount},
      {0, 3, 0, 0}, ScalarVertex,
      {0, 1, 0, 2}, ScalarVertex,
      {0, 4, 0, 0}, ScalarVertex,
      {2, 1, 0, 0}, ChiralVertex,
      {2, 0, 1, 0}, ChiralVertex,
      {0, 2, 1, 0}, MomentumDifferenceVertex,
      {0, 1, 2, 0}, InverseMetricVertex,
      {0, 2, 2, 0}, InverseMetricVertex,
      _, "(UnknownVertexType: " <> ToString[fields] <> ")"]
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
          ]

End[];
EndPackage[];

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

BeginPackage["CXXDiagrams`", {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "Parameters`","CConversion`", "SelfEnergies`"}];

(* This module generates c++ code intended to be used similarly to SARAH's fields and Vertex[] function *)

(* Vertex types *)
ScalarVertex::usage="";
ChiralVertex::usage="";
MomentumVertex::usage="";
TripleVectorVertex::usage="";
QuadrupleVectorVertex::usage="";
MomentumDifferenceVertex::usage="";
InverseMetricVertex::usage="";

VertexTypes::usage="";
VertexTypeForFields::usage="";
CXXNameOfField::usage="";
CXXNameOfVertex::usage="";
LorentzConjugateOperation::usage="";
LorentzConjugate::usage="";
RemoveLorentzConjugation::usage="";
AtomHead::usage="";
CreateFields::usage="";
FeynmanDiagramsOfType::usage="Obtain all instantiations of Feynman \
diagrams of a given topology with given external fields.";
VerticesForDiagram::usage="";
ContractionsBetweenVerticesForDiagramFromGraph::usage="";
CreateVertices::usage="Creates c++ code that makes functions available that \
numerically evaluate any of the given vertices.";
MaximumVerticesLimit::usage"";
VertexRulesForVertices::usage="";
CreateMassFunctions::usage="";
CreateUnitCharge::usage="";
StripLorentzIndices::usage="";
NumberOfFieldIndices::usage="";
FieldInfo::usage="";
includeLorentzIndices::usage="";
LoadVerticesIfNecessary::usage="";
IsLorentzIndex::usage = "";
IsGenerationIndex::usage = "";

Begin["`Private`"];

VertexTypes[] := {
    ScalarVertex,
    ChiralVertex,
    MomentumVertex,
    TripleVectorVertex,
    QuadrupleVectorVertex,
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

CXXNameOfVertex[fields_List] := "Vertex<" <> StringJoin[Riffle[
		CXXNameOfField[#, prefixNamespace -> "fields"] & /@ fields,
	", "]] <> ">"

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

(** \brief Get the lorentz index of a given field
 * \param[in] field The indexed field
 * \returns The lorentz index of \a field
 * For any given indexed field or antifield that has a lorentz index,
 * return it.
 *)
LorentzIndexOfField[field_]:=
	Module[{lorentzIndices = Select[Vertices`FieldIndexList[field], IsLorentzIndex]},
		Utils`AssertWithMessage[Length[lorentzIndices] === 1,
			"CXXDiagrams`LorentzIndexOfField[]: Argument " <>
			ToString[field] <> " does not have exactly one color index."];
    lorentzIndices[[1]]
  ]

(** \brief Obtain all instantiations of Feynman diagrams of a given
 * topology with given external fields.
 * \param adjacencyMatrix the symmetric adjacency matrix determining the
 * topology of the diagrams
 * \param externalFields a list of rules of the form
 * `{index1 -> externalField1, ... }`
 * where outgoing fields should be specified as conjugate fields and
 * incoming fields should be specified as unconjugate fields.
 * \returns a list of Feynman diagrams of the given topology with
 * given external fields.
 * All returned diagrams are of the form:
 * `{v1, v2, ..., vn}`
 * where `n` is the dimension of the given adjacency matrix and
 * `vk` for \f$k \in \{ 1, ..., n \} \f$ corresponds to the vertex
 * with index `k` as specified by the adjacency matrix.
 * A vertex in the returned diagram is either external or internal,
 * where external vertices are by definition taken to be those
 * specified by `externalFields`. Thus an external vertex `vk` is simply
 * represented as the external field `externalFieldk` as specified by
 * `externalFields`.
 * The internal vertices are given as lists of fields corresponding to 
 * terms in the interaction Lagrangian. They are obtained in
 * correspondence with the standard perturbative QFT approach when
 * calculating
 * \f[ \langle \left( \text{outgoing} \right)^* \vert
 * 	\text{incoming} \rangle \f]
 * where the asterisk denotes Lorentz conjugation. Of course only
 * diagrams of the given topology are taken into account.
 * The order of fields in the internal vertices is given by the 
 * adjacency matrix as follows:
 * Let `k` be the index of an internal vertex and `l` be a
 * non-zero entry in `adjacencyMatrix[[k, l]]`. Then, if `k =!= l`
 * `vk` will contain exactly `adjacencyMatrix[[k, l]]` entries of fields
 * that get contracted with `adjacencyMatrix[[k, l]]` fields in the
 * vertex `vl`. The order of those `adjacencyMatrix[[k, l]]` fields is
 * unspecified. If `k === l` `vk` will contain exactly 
 * `2 * adjacencyMatrix[[k, l]]` entries of fields where
 * the first one gets contracted with the second one, the third one with
 * the fourth one and so on. Otherwise the order of these
 * `2 * adjacencyMatrix[[k, l]]` fields is unspecified.
 * \note Here, 'contracted' means that contracted fields will be
 * Lorentz conjugates to each other - unless one of the fields is
 * external, then they will be the same.
 * Finally, the vertex `vk` is given by concatenating the thusly
 * constructed lists of fields for all `l` where 
 * `adjacencyMatrix[[k, l]]` is non-zero in the order given by
 * the index `l`.
 * Furthermore, two diagrams are considered equal if they only differ by
 * a reordering of fields in internal vertices and thus any appearing
 * duplicate diagrams are removed.
 **)
FeynmanDiagramsOfType[adjacencyMatrix_List,externalFields_List] :=
	Module[{externalVertices = externalFields[[All,1]],
			internalVertices,externalRules, internalFieldCouplings,
			unspecifiedEdgesLess,unspecifiedEdgesEqual,
			insertFieldRulesLess,insertFieldRulesGreater,insertFieldRulesEqual,
			fieldsToInsert,
			unresolvedFieldCouplings,resolvedFields,resolvedFieldCouplings,
			diagrams},
	LoadVerticesIfNecessary[];
	
	internalVertices = Complement[Table[k,{k,Length[adjacencyMatrix]}],externalVertices];
	externalRules = Flatten @ ({{_,#,_} :> (# /. externalFields)} & /@
		externalVertices);

	internalFieldCouplings = (Flatten[(Flatten @ Position[adjacencyMatrix[[#]],Except[0],{1},Heads -> False]
		/. {i_Integer :> Table[{#,i,k},{k,adjacencyMatrix[[#,i]]}]}),1] &
		/@ internalVertices) /. externalRules;

	unspecifiedEdgesLess = Cases[internalFieldCouplings,{i_,j_,_} /; i < j,{2}];
	unspecifiedEdgesEqual = Cases[internalFieldCouplings,{i_,i_,_},{2}];

	insertFieldRulesLess = MapIndexed[#1 -> SARAH`FieldToInsert[#2[[1]]] &,unspecifiedEdgesLess];
	insertFieldRulesGreater = (insertFieldRulesLess /. {Rule[{i_,j_,k_},field_] :> Rule[{j,i,k},SARAH`AntiField[field]]});
	insertFieldRulesEqual = MapIndexed[#1 -> Sequence @@ {SARAH`FieldToInsert[#2[[1]]+Length[insertFieldRulesLess]],
		SARAH`AntiField[SARAH`FieldToInsert[#2[[1]]+Length[insertFieldRulesLess]]]} &,
			unspecifiedEdgesEqual];
	fieldsToInsert = Table[SARAH`FieldToInsert[k],
		{k,Length[insertFieldRulesLess] + Length[insertFieldRulesEqual]}];
	
	unresolvedFieldCouplings = internalFieldCouplings
		/. insertFieldRulesLess /. insertFieldRulesGreater /. insertFieldRulesEqual;

	resolvedFields = If[fieldsToInsert === {}, {{}},
		SARAH`InsFields[{C @@@ unresolvedFieldCouplings,
			fieldsToInsert}][[All,2]]];

	If[resolvedFields === {}, Return[{}]];

	resolvedFieldCouplings = unresolvedFieldCouplings /.
		((Rule @@@ Transpose[{fieldsToInsert,#}]) & /@ resolvedFields);

	diagrams = Table[k,{k,Length[adjacencyMatrix]}] /. externalFields /. 
		((Rule @@@ Transpose[{internalVertices,#}]) & /@ resolvedFieldCouplings);
  
	(* Prevent overcounting of diagrams by removing diagrams that only
	 * differ by permutations of internal fields within the internal
	 * vertices. This is automatically performed by SARAH if 
	 * SA`CheckSameVertices === True, but it is better to not depend on
	 * some internal SARAH state. *)
	DeleteDuplicates[diagrams,
		(And @@ ((Sort[#[[1]]] === Sort[#[[2]]] &) /@
			Cases[Transpose[{#1, #2}],{{___},{___}}] (* Only check internal vertices *)
		) &)]
  ]

VerticesForDiagram[diagram_] := Select[diagram,Length[#] > 1 &]

ContractionsBetweenVerticesForDiagramFromGraph[v1_Integer, v2_Integer,
		diagram_List, graph_List] :=
	Module[{fields1 = diagram[[v1]], fields2 = diagram[[v2]],
			preceedingNumberOfFields1 = Total[graph[[v1, ;;v2]]] - graph[[v1,v2]],
			preceedingNumberOfFields2 = Total[graph[[v2, ;;v1]]] - graph[[v2,v1]],
			contractedFieldIndices1, contractedFieldIndices2,
			stepSize = If[v1 === v2, 2, 1]},
		If[v1 === v2,
			preceedingNumberOfFields2 = preceedingNumberOfFields2 + 1];
		If[v1 < v2,
			preceedingNumberOfFields1 =
				preceedingNumberOfFields1 + graph[[v1,v1]]];
		If[v2 < v1,
			preceedingNumberOfFields2 =
				preceedingNumberOfFields2 + graph[[v2,v2]]];
		
		contractedFieldIndices1 = Table[k, {k, preceedingNumberOfFields1 + 1,
			preceedingNumberOfFields1 + stepSize * graph[[v1,v2]], stepSize}];
		contractedFieldIndices2 = Table[k, {k, preceedingNumberOfFields2 + 1,
			preceedingNumberOfFields2 + stepSize * graph[[v2,v1]], stepSize}];
		
		Transpose[{contractedFieldIndices1, contractedFieldIndices2}]
	]

(** \brief Creates c++ code that makes functions available that
 * numerically evaluate any of the given vertices.
 * \param vertices a list of vertices
 * \param StripColorStructure A boolean option to specify whether
 * to strip away parts of vertices possessing colour structures other
 * than a single ``SARAH`Delta[]``. The default value is `False`.
 * \param MaximumVerticesLimit An integer option that specify an upper
 * limit of vertices that shall go into a single block of code.
 * \returns a list `{{prototypes1, definitions1}, ...}` containing the
 * corresponding c++ code where no sublist contains more than
 * `MaximumVerticesLimit` number of vertices.
 **)
CreateVertices[vertices_List,
		OptionsPattern[{StripColorStructure -> False,
			MaximumVerticesLimit -> 500}]] :=
	Module[{cxxVertices, vertexPartition},
		cxxVertices = CreateVertex[#, StripColorStructure -> OptionValue[StripColorStructure]] & /@
			DeleteDuplicates[vertices];
		
		(* Mathematica 7 does not support the `UpTo[n]` notation *)
		vertexPartition = Partition[cxxVertices, OptionValue[MaximumVerticesLimit]];
		If[vertexPartition === {},
			vertexPartition = {cxxVertices}];
		
		Map[StringJoin[Riffle[#, "\n\n"]] &, Transpose /@ vertexPartition, {2}]
	]

CreateVertex[fields_List, OptionsPattern[{StripColorStructure -> False}]] :=
  Module[{fieldSequence},
		LoadVerticesIfNecessary[];
		
		fieldSequence = StringJoin @ Riffle[
			CXXNameOfField[#, prefixNamespace -> "fields"] & /@ fields, ", "];

		{
		"template<> struct VertexImpl<" <> fieldSequence <> ">" <> "\n" <>
		"{\n" <> TextFormatting`IndentText[
			"static " <> SymbolName[VertexTypeForFields[fields]] <>
				" evaluate(const std::array<int, " <>
				ToString[Total[NumberOfFieldIndices /@ fields]] <>
			">& indices, const context_base& context);"] <> "\n" <>
		"};"
		,
		SymbolName[VertexTypeForFields[fields]] <>
			" VertexImpl<" <> fieldSequence <> ">::evaluate(\n" <>
				TextFormatting`IndentText["const std::array<int, " <>
					ToString[Total[NumberOfFieldIndices /@ fields]] <> ">& indices, " <>
				"const context_base& context)"] <> "\n" <>
		"{\n" <> TextFormatting`IndentText[
				VertexFunctionBodyForFields[fields,
					StripColorStructure -> OptionValue[StripColorStructure]] <> "\n"] <>
		"};"
		}
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
  Module[{sortedFields, sortedIndexedFields, indexedFields, ghosts,
          fieldsOrdering, sortedFieldsOrdering, inverseFOrdering,
          fOrderingWRTSortedF, vertex, vExpression, vertexIsZero = False,
          vertexType = VertexTypeForFields[fields], expr, exprL, exprR,
          vertexRules, incomingScalar, outgoingScalar, incomingGhost,
					lIndex1, lIndex2, lIndex3, lIndex4, expr1, expr2, expr3},
    sortedFields = Vertices`SortFieldsInCp[fields];
    
    vertex = Select[vertexList, StripFieldIndices[#[[1]]] === sortedFields &, 1];
    If[vertex === {}, vertexIsZero = True, vertex = vertex[[1]]];

    If[vertexIsZero,
       Return[Switch[vertexType,
         ScalarVertex,
         "return {0};",
         
         ChiralVertex,
         "return {0, 0};",
         
         MomentumVertex,
         ghosts = Select[fields, TreeMasses`IsGhost];
         ghosts = Select[ghosts, RemoveLorentzConjugation[#] =!= 
           LorentzConjugate[#] &];
         
         Utils`AssertWithMessage[Length[ghosts] === 1,
				   "CXXDiagrams`VertexFunctionBodyForFieldsImpl[]: \
Invalid ghosts in vertex: " <> ToString[fields]];
				   
         "return {0, " <>
           ToString[Position[fields, ghosts[[1]], {1}][[1,1]] - 1] <>
         "};",
         
         TripleVectorVertex,
         "return {0, TripleVectorVertex::even_permutation{}};",
         
         QuadrupleVectorVertex,
         "return {0, 0, 0};",
         
         MomentumDifferenceVertex,
         "return {0, " <> StringJoin[Riffle[
            ToString /@ Flatten[Position[fields,
               field_ /; TreeMasses`IsScalar[field] || TreeMasses`IsGhost[field],
               {1}, Heads -> False] - 1],
            ", "]] <> "};",
         
         InverseMetricVertex,
         "return {0};"]]];

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
      "return {result};",
      
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
      "return {left, right};",
      
      MomentumVertex,
      incomingGhost = Replace[vertex[[2,2]], SARAH`Mom[ig_,_] :> ig];
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
      "return {result, " <> 
				ToString[Position[indexedFields, incomingGhost,{1}][[1,1]] - 1] <>
			"};",
         
			TripleVectorVertex,
			{lIndex1, lIndex2, lIndex3} = LorentzIndexOfField /@ 
				sortedIndexedFields;
			
      vertexRules = {(SARAH`Cp @@ sortedIndexedFields) -> vertex[[2,1]]};
      
      expr = CanonicalizeCoupling[SARAH`Cp @@ fields,
        sortedFields, sortedIndexedFields] /. vertexRules;
      
      expr = Switch[Thread[Expand /@ {vertex[[2,2]], - vertex[[2,2]]} == Expand[
				SARAH`g[lIndex1, lIndex2] * (SARAH`Mom[sortedIndexedFields[[1]], lIndex3] -
					SARAH`Mom[sortedIndexedFields[[2]], lIndex3]) +
				SARAH`g[lIndex2, lIndex3] * (SARAH`Mom[sortedIndexedFields[[2]], lIndex1] -
					SARAH`Mom[sortedIndexedFields[[3]], lIndex1]) + 
				SARAH`g[lIndex1, lIndex3] * (SARAH`Mom[sortedIndexedFields[[3]], lIndex2] -
					SARAH`Mom[sortedIndexedFields[[1]], lIndex2])]],
				
				{True, _}, expr,
				{_, True}, - expr,
				_, Print["CXXDiagrams`VertexFunctionBodyForFieldsImpl[]: Unknown Lorentz \
structure in vertex: " <> ToString[vertex]]; Quit[1]];
      
      expr = Vertices`SarahToFSVertexConventions[sortedFields, expr,
				StripColorStructure -> OptionValue[StripColorStructure]];
      expr = TreeMasses`ReplaceDependenciesReverse[expr];
      DeclareIndices[StripLorentzIndices /@ indexedFields, "indices"] <>
      Parameters`CreateLocalConstRefs[expr] <> "\n" <>
      "const " <> GetComplexScalarCType[] <> " result = " <>
      Parameters`ExpressionToString[expr] <> ";\n\n" <>
      "return {result, TripleVectorVertex::even_permutation{}};",
         
			QuadrupleVectorVertex,
			{lIndex1, lIndex2, lIndex3, lIndex4} = LorentzIndexOfField /@ 
				sortedIndexedFields;
      
      vertexRules = {
				(SARAH`Cp @@ sortedIndexedFields)[
					SARAH`g[lIndex1, lIndex2] * SARAH`g[lIndex3, lIndex4]] ->
        Vertices`FindVertexWithLorentzStructure[Rest[vertex],
					SARAH`g[lIndex1, lIndex2] * SARAH`g[lIndex3, lIndex4]][[1]],
				(SARAH`Cp @@ sortedIndexedFields)[
					SARAH`g[lIndex1, lIndex3] * SARAH`g[lIndex2, lIndex4]] ->
        Vertices`FindVertexWithLorentzStructure[Rest[vertex],
					SARAH`g[lIndex1, lIndex3] * SARAH`g[lIndex2, lIndex4]][[1]],
				(SARAH`Cp @@ sortedIndexedFields)[
					SARAH`g[lIndex1, lIndex4] * SARAH`g[lIndex2, lIndex3]] ->
        Vertices`FindVertexWithLorentzStructure[Rest[vertex],
					SARAH`g[lIndex1, lIndex4] * SARAH`g[lIndex2, lIndex3]][[1]]
			};
      
      expr1 = CanonicalizeCoupling[(SARAH`Cp @@ fields)[
					SARAH`g[lIndex1, lIndex2] * SARAH`g[lIndex3, lIndex4]],
        sortedFields, sortedIndexedFields] /. vertexRules;
      expr2 = CanonicalizeCoupling[(SARAH`Cp @@ fields)[
					SARAH`g[lIndex1, lIndex3] * SARAH`g[lIndex2, lIndex4]],
        sortedFields, sortedIndexedFields] /. vertexRules;
      expr3 = CanonicalizeCoupling[(SARAH`Cp @@ fields)[
					SARAH`g[lIndex1, lIndex4] * SARAH`g[lIndex2, lIndex3]],
        sortedFields, sortedIndexedFields] /. vertexRules;

      expr1 = Vertices`SarahToFSVertexConventions[sortedFields, expr1,
				StripColorStructure -> OptionValue[StripColorStructure]];
      expr2 = Vertices`SarahToFSVertexConventions[sortedFields, expr2,
				StripColorStructure -> OptionValue[StripColorStructure]];
      expr3 = Vertices`SarahToFSVertexConventions[sortedFields, expr3,
				StripColorStructure -> OptionValue[StripColorStructure]];
			
      expr1 = TreeMasses`ReplaceDependenciesReverse[expr1];
      expr2 = TreeMasses`ReplaceDependenciesReverse[expr2];
      expr3 = TreeMasses`ReplaceDependenciesReverse[expr3];
      
      DeclareIndices[StripLorentzIndices /@ indexedFields, "indices"] <>
      Parameters`CreateLocalConstRefs[{expr1, expr2, expr3}] <> "\n" <>
      "const " <> GetComplexScalarCType[] <> " part1 = " <>
      Parameters`ExpressionToString[expr1] <> ";\n\n" <>
      "const " <> GetComplexScalarCType[] <> " part2 = " <>
      Parameters`ExpressionToString[expr2] <> ";\n\n" <>
      "const " <> GetComplexScalarCType[] <> " part3 = " <>
      Parameters`ExpressionToString[expr3] <> ";\n\n" <>
      "return {part1, part2, part3};",

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
      "return {result, minuend_index, subtrahend_index};",

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
      "return {result};"
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

         "ChiralVertex unit_charge(const context_base& context)\n" <>
         "{\n" <>
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
      {0, 0, 1, 2}, MomentumVertex,
      {0, 0, 3, 0}, TripleVectorVertex,
      {0, 0, 4, 0}, QuadrupleVectorVertex,
      {0, 4, 0, 0}, ScalarVertex,
      {2, 1, 0, 0}, ChiralVertex,
      {2, 0, 1, 0}, ChiralVertex,
      {0, 2, 1, 0}, MomentumDifferenceVertex,
      {0, 1, 2, 0}, InverseMetricVertex,
      {0, 2, 2, 0}, InverseMetricVertex,
      _, "(UnknownVertexType: " <> ToString[fields] <> ")"]
  ]

IsLorentzIndex[index_] := StringMatchQ[ToString @ index, "lt" ~~ __];
IsGenerationIndex[index_] := StringMatchQ[ToString @ index, "gt" ~~ __];

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

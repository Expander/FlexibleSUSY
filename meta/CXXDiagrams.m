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

(*** This module generates c++ code capable of representing fields, vertices and diagrams. ***)
BeginPackage["CXXDiagrams`", {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "Parameters`", "CConversion`", "ColorMath`"}];


(*** Debug only ***)
SortedVertex::usage="";
GaugeStructureOfVertexLorentzPart::usage="";
FullGaugeStructureFromParts::usage="";
GaugeStructureOfVertex::usage="";
IsLorentzIndex::usage="";
IsColourIndex::usage="";
IndexFields::usage="";

(*** Public interfaces that model fields ***)
CreateFields::usage="";
CXXNameOfField::usage="";
LorentzConjugateOperation::usage="";
LorentzConjugate::usage="";
RemoveLorentzConjugation::usage="";
NumberOfFieldIndices::usage="";
CreateMassFunctions::usage="";
includeLorentzIndices::usage="";

(*** Public interfaces that model vertices ***)
(** We need to encode the vertex structure for every unbroken gauge group:
- Identity component of Poincare group
- SU(3)) **)
(* A basis for the currently implemented Lorentz structures of vertices *)
ScalarVertex::usage="A Lorentz scalar (invariant) vertex";
ChiralVertex::usage="A vertex decomposing into left/right projectors";
MomentumDifferenceVertex::usage="A vertex with a momentum difference (p^\\mu - q^\\nv)";
InverseMetricVertex::usage="A vertex proportional to g^{\\mu \\nu}";

(* A basis for the currently implemented colour structures of vertices *)
UncolouredVertex::usage="A SU(3) invariant vertex";
KroneckerDeltaColourVertex::usage="A vertex proportional to \\delta( ct1, ct2 )";
GellMannVertex::usage="A vertex proportional to \\Lambda^t_{a b}";
AdjointlyColouredVertex::usage="A vertex proportional to f^{a b c}";

CreateVertices::usage="";
VertexRulesForVertices::usage="";

(*** Public interfaces that model diagrams ***)
FeynmanDiagramsOfType::usage="";
VerticesForDiagram::usage="";
ColourFactorForDiagramFromGraph::usage="";
ContractionsBetweenVerticesForDiagramFromGraph::usage="";

(*** Public utility functions ***)
AtomHead::usage="";
CreateUnitCharge::usage="";

Begin["`Private`"];

(* Internally, we need to represent left chiral and right chiral
 vertices separately. *)
LeftChiralVertex::usage="A left projector part of a vertex";
RightChiralVertex::usage="A right projector part of a vertex";

(* Return a string corresponding to the c++ class name of the field.
 Note that "bar" and "conj" get turned into "typename bar<...>::type" and
 "typename conj<...>::type" respectively! *)
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

(* Convert color structures to the ColorMath convention *)
ConvertColourStructureToColorMathConvention[fields_List,
	GellMannVertex[cIndex1_, cIndex2_, cIndex3_]] :=
(* FIXME: Factor of two? *)
	2 ColorMath`CMt[{cIndex1}, cIndex2, cIndex3] 

ConvertColourStructureToColorMathConvention[fields_List,
	AdjointlyColouredVertex[cIndex1_, cIndex2_, cIndex3_]] :=
	ColorMath`CMf[cIndex1, cIndex2, cIndex3] 

ConvertColourStructureToColorMathConvention[indexedFields_List,
	KroneckerDeltaColourVertex[cIndex1_, cIndex2_]] := 
	Module[{colouredField1, colouredField2, colourRep1, colourRep2},
		(* If the result has a color Delta, we need to find out if it's adj. or fundamental
			because they are represented by different symbols in ColorMath.
			Also, in ColorMath the order of indices matters. As stated in  ColorMath tutorial notebook:

			"Note that lower and upper indices are different. One position is taken to define incoming
			quarks/outgoing anti-quarks, and the other position is used for outgoing quarks and incoming
			anti-quarks. What is taken to be what is a question of definition."
		*)
    
    colouredField1 = Select[indexedFields, !FreeQ[#, cIndex1] &][[1]];
    colouredField2 = Select[indexedFields, !FreeQ[#, cIndex2] &][[1]];
    
    colourRep1 = SARAH`getColorRep[colouredField1];
    colourRep2 = SARAH`getColorRep[colouredField2];
		
    Utils`AssertWithMessage[colourRep1 == colourRep2,
			"CXXDiagrams`ConvertColourStructureToColorMathConvention[]: " <>
			"Two colour indices in Kronecker delta that come from fields " <>
			"of incompatible representations."];
		
		(* FIXME: Are these orderings correct? *)
		Switch[{colourRep1},
			{SARAH`T}, If[RemoveLorentzConjugation[colouredField1] === colouredField1,
				ColorMath`CMdelta[cIndex2, cIndex1],
				ColorMath`CMdelta[cIndex1, cIndex2]
			],
			{O}, ColorMath`CMDelta[cIndex2, cIndex1],
			_,
			Print["CXXDiagrams`ConvertColourStructureToColorMathConvention[]: " <>
				"Two colour indices in Kronecker delta that come from fields " <>
				"or representations we cannot handle: " <>
				ToString[{colourRep1, colourRep2}]];
			Quit[1];
		]
	]

(** \brief Calculate the colour factor of a given diagram with a given
 * topology.
 * \param[in] diagram The diagram specifiying all occurring fields
 * \param[in] graph The topology of the diagram as adjacency matrix
 * \returns The color factor for \a diagram.
 *)
ColourFactorForDiagramFromGraph[diagram_, graph_] :=
	Module[{indexedFields = IndexFields[Flatten[diagram]],
			vertexFields = Cases[diagram, {__}],
			indexedVertexFields,
			indexedDiagram, sortedVertices, vertices,
			verticesOrderings, sortedVerticesOrderings,
			fOrderingWRTSortedVs, defaultIndexedVertices,
			colourStructures, colorMathExpression},
		indexedDiagram = diagram /. Rule @@@ Transpose[{
			#[Unique[]] & /@ Flatten[diagram], indexedFields}];
		indexedVertexFields = Cases[indexedDiagram, {__}];
		
		(** TODO: Connect the color indices **)
		
		sortedVertices = SortedVertex /@ vertexFields;
		
    (* Mathematica 7 does not know about permutations... :'-( *)
    verticesOrderings = Ordering /@ vertexFields;
    sortedVerticesOrderings = Ordering /@ sortedVertices[[All,1]];

    inverseVOrderings = Ordering /@ verticesOrderings;

    fOrderingWRTSortedVs = Part @@@ Transpose[
			{sortedVerticesOrderings, inverseVOrderings}];

    defaultIndexedVertices = Part @@@ Transpose[
			{sortedVerticesOrderings, fOrderingWRTSortedVs}];
		
		vertices = (#[[1]] /. (Rule @@@
			Map[Vertices`FieldIndexList,
				Transpose[{#[[2]], #[[3]]}], {2}]
			) &) /@
			Transpose[{sortedVertices,
				defaultIndexedVertices, indexedVertexFields}];
		
		colourStructures = Cases[
			Transpose[{
				vertices[[All,1]],
				(GaugeStructureOfVertex /@ vertices)[[All,2]]
			}],
			Except[UncolouredVertex]];
		
		Print[colourStructures];
		
		colorMathExpression = Times @@
			(ConvertColourStructureToColorMathConvention @@@ colourStructures);
		
		If[colorMathExpression === 1, 1,
			ColorMath`CSimplify[colorMathExpression]]
	]

IndexPrefixForType[SARAH`generation] := "gt"
IndexPrefixForType[SARAH`lorentz] := "lt"
IndexPrefixForType[SARAH`color] := "ct"
IndexPrefixForType[indexType_] := 
	(Print["Unknown index type: " <> ToString[indexType] <> " encountered."]; Quit[1])

IndexFields[fields_List] :=
	Module[{indexSpecifications, indexTypes, indexNames},
		indexSpecifications =
			GetParticleIndices[SARAH`getParticleName[#]] & /@ fields;
		indexSpecifications =
			Cases[#, Except[{SARAH`generation, 1}]] & /@ indexSpecifications;
		
		indexTypes = Map[First, indexSpecifications, {2}];
		
		indexNames = MapIndexed["SARAH`" <> #1 <> ToString[#2[[1]]] &,
			Map[IndexPrefixForType, indexTypes, {2}], {2}];
		
		(#[[1]][Symbol /@ #[[2]]] & /@ Transpose[{fields, indexNames}]) /.
			field_[{}] :> field
	]

SortedVertex[fields_List] :=
	Module[{sortedFields, indexedSortedFields, vertexList, vertex, types,
			similarVertexList, similarVertex, fieldReplacementRules},
    LoadVerticesIfNecessary[];
		sortedFields = Vertices`SortFieldsInCp[fields];
    
    vertexList = Switch[Length[fields],
			3, SARAH`VertexList3,
			4, SARAH`VertexList4,
			_, Print["CXXDiagrams`SortedVertex[]: Only three- and four-point vertices are supported"];
				Quit[1]];
    
    vertex = Select[vertexList, Vertices`StripFieldIndices[#[[1]]] === sortedFields &, 1];
    vertex = If[vertex =!= {}, vertex[[1]],
			types = SARAH`getType[#, False, FlexibleSUSY`FSEigenstates] & /@ sortedFields;
			similarVertexList = SA`VertexList[Symbol[StringJoin[SymbolName /@ types]]];
			
			Utils`AssertWithMessage[similarVertexList =!= {},
				"CXXDiagrams`SortedVertex[]: Cannot determine Lorentz structure of " <> ToString[fields] <> " vertex."];
			similarVertex = similarVertexList[[1]];
			
			indexedSortedFields = IndexFields[sortedFields];
			fieldReplacementRules = Rule @@@ Transpose[{similarVertex[[1]], indexedSortedFields}];
			
			{indexedSortedFields, Sequence @@ Transpose[{
				Table[0, {Length[similarVertex] - 1}],
					Rest[similarVertex][[All,2]] /. fieldReplacementRules}]}
		];
		
		(* Be consistent with the conventions in Vertices.m *)
		{vertex[[1]], Sequence @@ Transpose[
			{-I * Rest[vertex][[All,1]], Rest[vertex][[All,2]]}]}
	]

LabelLorentzPart[{scalar_, 1}] := {scalar, ScalarVertex}

LabelLorentzPart[{scalar_,
		SARAH`PL | SARAH`LorentzProduct[
			SARAH`gamma[lIndex_] /; IsLorentzIndex[lIndex], SARAH`PL]}] :=
	{scalar, LeftChiralVertex}
	
LabelLorentzPart[{scalar_,
		SARAH`PR | SARAH`LorentzProduct[
			SARAH`gamma[lIndex_] /; IsLorentzIndex[lIndex], SARAH`PR]}] :=
	{scalar, RightChiralVertex}
	
LabelLorentzPart[{scalar_, SARAH`Mom[in_, lIndex1_] - SARAH`Mom[out_, lIndex2_]}] :=
	{scalar, MomentumDifferenceVertex[in, out]} /;
	IsLorentzIndex[lIndex1] && IsLorentzIndex[lIndex2]

LabelLorentzPart[{scalar_, SARAH`g[lIndex1_, lIndex2_]}] :=
	{scalar, InverseMetricVertex} /;
	IsLorentzIndex[lIndex1] && IsLorentzIndex[lIndex2]

LabelLorentzPart[vertexPart_] := 
	(Print["Unknown Lorentz structure in vertex ", vertexPart]; Quit[1])


GaugeStructureOfVertexLorentzPart[{scalar_, lorentzStructure_}] :=
	{scalar, UncolouredVertex, lorentzStructure} /;
	FreeQ[scalar, atom_ /; IsColourIndex[atom], -1]

GaugeStructureOfVertexLorentzPart[
	{scalar_ * SARAH`Delta[cIndex1_, cIndex2_], lorentzStructure_}] :=
	{scalar, KroneckerDeltaColourVertex[cIndex1, cIndex2], lorentzStructure} /;
	FreeQ[scalar, atom_ /; IsColourIndex[atom], -1]
	
GaugeStructureOfVertexLorentzPart[
	{scalar_ * SARAH`Lam[cIndex1_, cIndex2_, cIndex3_], lorentzStructure_}] :=
	{scalar, GellMannVertex[cIndex1, cIndex2, cIndex3], lorentzStructure} /;
	FreeQ[scalar, atom_ /; IsColourIndex[atom], -1]
	
GaugeStructureOfVertexLorentzPart[
	{scalar_ * SARAH`fSU3[cIndex1_, cIndex2_, cIndex3_], lorentzStructure_}] :=
	{scalar, AdjointlyColouredVertex[cIndex1, cIndex2, cIndex3], lorentzStructure} /;
	FreeQ[scalar, atom_ /; IsColourIndex[atom], -1]

GaugeStructureOfVertexLorentzPart[vertexPart_] := 
	(Print["Unknown colour structure in vertex ", vertexPart]; Quit[1])

(* Assuming all Gauge structure combinations are different *)
FullGaugeStructureFromParts[gaugeParts_List] :=
	Module[{leftChiralParts, rightChiralParts, nonChiralParts, chiralParts,
			gaugeStructures},
		leftChiralParts = Cases[gaugeParts, {_, _, LeftChiralVertex}];
		rightChiralParts = Cases[gaugeParts, {_, _, RightChiralVertex}];
		nonChiralParts = Complement[gaugeParts, leftChiralParts, rightChiralParts];
		
		chiralParts = Flatten[Cases[rightChiralParts,
			{rightScalar_, #[[2]], RightChiralVertex} :>
			{{#[[1]], rightScalar}, #[[2]], ChiralVertex}] &
			/@ leftChiralParts, 1];
			
		leftChiralParts = Complement[leftChiralParts, chiralParts,
			SameTest -> (#1[[2]] == #2[[2]] &)];
		rightChiralParts = Complement[rightChiralParts, chiralParts,
			SameTest -> (#1[[2]] == #2[[2]] &)];
		
		gaugeStructures = Sequence @@@
			{nonChiralParts, chiralParts, leftChiralParts, rightChiralParts};
		Utils`AssertWithMessage[Length[gaugeStructures] === 1,
			"Vertex with composite gauge structure encountered."];
		
		gaugeStructures[[1]]
	]

GaugeStructureOfVertex[vertex_] :=
	Module[{lorentzParts, gaugeStructures},
		(* This is arguably a bug in SARAH causing e.g {hh, hh, VZ} to be mapped to {0, 0} *)
		lorentzParts = If[MatchQ[sortedVertex, {_, {_, 0}}],
			Utils`AssertWithMessage[sortedVertex[[2,1]] === 0,
				"CXXDiagrams`SortedVertex[]: ill-formed SARAH vertex: " <> ToString[sortedVertex]];
			{{0, MomentumDifferenceVertex[sortedVertex[[1,1]], sortedVertex[[1,2]]]}},
			LabelLorentzPart /@ vertex[[2;;]]];
		
		gaugeStructures = GaugeStructureOfVertexLorentzPart /@ lorentzParts;
		FullGaugeStructureFromParts[gaugeStructures]
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

CreateVertices[vertices_List] :=
  StringJoin @ Riffle[CreateVertex /@ DeleteDuplicates[vertices], "\n\n"]

CreateVertex[fields_List] :=
  Module[{fieldSequence},
		fieldSequence = StringJoin @ Riffle[
			CXXNameOfField[#, prefixNamespace -> "fields"] & /@ fields, ", "];

    "template<> struct VertexImpl<" <> fieldSequence <> ">" <> "\n" <>
    "{\n" <> TextFormatting`IndentText[
			"using vertex_type = " <> SymbolName[VertexTypeForFields[fields]] <> ";\n\n" <>
			"static vertex_type evaluate(const typename Vertex<" <> fieldSequence <> 
				">::indices_type& indices, \n" <>
				TextFormatting`IndentText["const context_base& context)\n"] <>
			"{\n" <> TextFormatting`IndentText[
				VertexFunctionBodyForFields[fields] <> "\n"] <>
			"}" <> "\n"] <>
		"};"
  ]

VertexTypeForFields[fields_List] :=
	AtomHead[GaugeStructureOfVertex[SortedVertex[fields]][[3]]]

VertexFunctionBodyForFields[fields_List] := 
	Module[{sortedFields, sortedIndexedFields, indexedFields, indexFields,
          fieldsOrdering, sortedFieldsOrdering, inverseFOrdering,
          fOrderingWRTSortedF, expr, exprL, exprR,
          vertexRules, incomingScalar, outgoingScalar},
    sortedVertex = SortedVertex[fields];
    sortedIndexedFields = sortedVertex[[1]];
    sortedFields = Vertices`StripFieldIndices /@ sortedIndexedFields;
    
		indexFields = {SARAH`Cp @@ sortedFields -> SARAH`Cp @@ sortedIndexedFields};
    
    (* Mathematica 7 does not know about permutations... :'-( *)
    fieldsOrdering = Ordering[fields];
    sortedFieldsOrdering = Ordering[sortedFields];

    inverseFOrdering = Ordering[fieldsOrdering];
    fOrderingWRTSortedF = sortedFieldsOrdering[[inverseFOrdering]];

    indexedFields = sortedIndexedFields[[fOrderingWRTSortedF]];
    
    gaugeStructure = GaugeStructureOfVertex[sortedVertex];
    Switch[gaugeStructure[[3]],
      ScalarVertex | InverseMetricVertex,
      vertexRules = {(SARAH`Cp @@ sortedIndexedFields) -> gaugeStructure[[1]]};
      
      expr = Vertices`SortCp[SARAH`Cp @@ fields] /. indexFields /. vertexRules;
      
      DeclareIndices[StripUnbrokenGaugeIndices /@ indexedFields, "indices"] <>
      Parameters`CreateLocalConstRefs[expr] <> "\n" <>
      "const " <> GetComplexScalarCType[] <> " result = " <>
      Parameters`ExpressionToString[expr] <> ";\n\n" <>
      "return vertex_type(result);",
      
      ChiralVertex,
      vertexRules = {
        (SARAH`Cp @@ sortedIndexedFields)[SARAH`PL] -> gaugeStructure[[1,2]],
        (SARAH`Cp @@ sortedIndexedFields)[SARAH`PR] -> gaugeStructure[[1,2]]
      };

      exprL = Vertices`SortCp[(SARAH`Cp @@ fields)[SARAH`PL]] /.
				indexFields /. vertexRules;
      exprR = Vertices`SortCp[(SARAH`Cp @@ fields)[SARAH`PR]] /.
				indexFields /. vertexRules;

      DeclareIndices[StripUnbrokenGaugeIndices /@ indexedFields, "indices"] <>
      Parameters`CreateLocalConstRefs[{exprL, exprR}] <> "\n" <>
      "const " <> GetComplexScalarCType[] <> " left = " <>
				Parameters`ExpressionToString[exprL] <> ";\n\n" <>
      "const " <> GetComplexScalarCType[] <> " right = " <>
				Parameters`ExpressionToString[exprR] <> ";\n\n" <>
      "return vertex_type(left, right);",

      MomentumDifferenceVertex[__],
      {incomingScalar, outgoingScalar} = List @@ gaugeStructure[[3]];
      vertexRules = {(SARAH`Cp @@ sortedIndexedFields) -> gaugeStructure[[1]]};
      
      expr = Vertices`SortCp[SARAH`Cp @@ fields] /. indexFields /. vertexRules;
      
      "int minuend_index = " <> 
        ToString[Position[indexedFields, incomingScalar][[1,1]] - 1] <> ";\n" <>
      "int subtrahend_index = " <>
        ToString[Position[indexedFields, outgoingScalar][[1,1]] - 1] <> ";\n\n" <>
      DeclareIndices[StripUnbrokenGaugeIndices /@ indexedFields, "indices"] <>
      Parameters`CreateLocalConstRefs[expr] <> "\n" <>
      "const " <> GetComplexScalarCType[] <> " result = " <>
      Parameters`ExpressionToString[expr] <> ";\n\n" <>
      "return vertex_type(result, minuend_index, subtrahend_index);",
      
      _,
      Print["Unrecognized gauge structure: " <> ToString[gaugeStructure[[3]]]];
      Quit[1]
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
  
CreateMassFunctions[] :=
  Module[{massiveFields,
          ghostMappings = SelfEnergies`ReplaceGhosts[FlexibleSUSY`FSEigenstates]},
    massiveFields = TreeMasses`GetParticles[];

    StringJoin @ Riffle[
      Module[{fieldInfo = TreeMasses`FieldInfo[#], numberOfIndices},
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
                 " " <> ToString @ (TreeMasses`FieldInfo[electron][[2]]-1) <> (* Electron has the lowest index *)
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
                 " " <> ToString @ (TreeMasses`FieldInfo[photon][[2]]-1) <>
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

NumberOfFieldIndices[field_] := Length @ TreeMasses`FieldInfo[field][[5]]

IndexBoundsForField[field_] :=
  Module[{fieldInfo = TreeMasses`FieldInfo[field]},
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

IsColourIndex[index_] := StringMatchQ[ToString @ index, "ct" ~~ __]
IsLorentzIndex[index_] := StringMatchQ[ToString @ index, "lt" ~~ __]

StripUnbrokenGaugeIndices[p_] := StripColourIndices[StripLorentzIndices[p]]

StripLorentzIndices[p_Symbol] := p
StripLorentzIndices[SARAH`bar[p_]] := SARAH`bar[StripLorentzIndices[p]]
StripLorentzIndices[Susyno`LieGroups`conj[p_]] := Susyno`LieGroups`conj[StripLorentzIndices[p]]
StripLorentzIndices[p_] := 
	Module[{remainingIndices},
		remainingIndices = Select[p[[1]], (!IsLorentzIndex[#] &)];
		If[Length[remainingIndices] === 0, Head[p],
			Head[p][remainingIndices]]
	]

StripColourIndices[p_Symbol] := p
StripColourIndices[SARAH`bar[p_]] := SARAH`bar[StripColourIndices[p]]
StripColourIndices[Susyno`LieGroups`conj[p_]] := Susyno`LieGroups`conj[StripColourIndices[p]]
StripColourIndices[p_] := 
	Module[{remainingIndices},
		remainingIndices = Select[p[[1]], (!IsColourIndex[#] &)];
		If[Length[remainingIndices] === 0, Head[p],
			Head[p][remainingIndices]]
	]

End[];
EndPackage[];

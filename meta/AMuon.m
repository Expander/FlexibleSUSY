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

BeginPackage["AMuon`", {"SARAH`", "CXXDiagrams`", "TextFormatting`", "TreeMasses`", "LoopMasses`"}];

AMuonContributingGraphs::usage="";
AMuonGetMuon::usage="";
AMuonContributingDiagramsForGraph::usage="";
AMuonCreateMuonPhysicalMass::usage="";
AMuonCreateCalculation::usage="";
AMuonGetMSUSY::usage="";

Begin["`Private`"];

(* The graphs that contribute to the EDM are precisely those with three
   external lines given by the field in question, its Lorentz conjugate
   and a photon.
   They are given as a List of undirected adjacency matrices where
    1 is the field itself
    2 is its Lorentz conjugate
    3 is the photon
   and all other indices unspecified. *)
vertexCorrectionGraph = {{0,0,0,1,0,0},
                         {0,0,0,0,1,0},
                         {0,0,0,0,0,1},
                         {1,0,0,0,1,1},
                         {0,1,0,1,0,1},
                         {0,0,1,1,1,0}};
contributingGraphs = {vertexCorrectionGraph};

AMuonContributingGraphs[] := contributingGraphs

AMuonGetMuon[] := If[TreeMasses`GetDimension[TreeMasses`GetSMMuonLeptonMultiplet[]] =!= 1,
                TreeMasses`GetSMMuonLeptonMultiplet[],
                Cases[SARAH`ParticleDefinitions[FlexibleSUSY`FSEigenstates],
                      {p_, {Description -> "Muon", ___}} -> p, 1][[1]]
               ]
GetCXXMuonIndex[] := If[TreeMasses`GetDimension[TreeMasses`GetSMMuonLeptonMultiplet[]] =!= 1,
                        1,
                        Null]

AMuonContributingDiagramsForGraph[graph_] :=
  Module[{diagrams},
    diagrams = CXXDiagrams`FeynmanDiagramsOfType[graph,
         {1 -> AMuonGetMuon[], 2 -> CXXDiagrams`LorentzConjugate[AMuonGetMuon[]], 3 -> TreeMasses`GetPhoton[]}];
         
    Select[diagrams,IsDiagramSupported[graph,#] &]
 ]
 
IsDiagramSupported[vertexCorrectionGraph,diagram_] :=
  Module[{photonEmitter,exchangeParticle},
    photonEmitter = diagram[[4,3]]; (* Edge between vertices 4 and 6 (3rd edge of vertex 4) *)
    exchangeParticle = diagram[[4,2]]; (* Edge between vertices 4 and 5 (2nd edge of vertex 4) *)
    
    If[diagram[[6]] =!= {TreeMasses`GetPhoton[],CXXDiagrams`LorentzConjugate[photonEmitter],photonEmitter},
       Return[False]];
    If[TreeMasses`IsFermion[photonEmitter] && TreeMasses`IsScalar[exchangeParticle],
       Return[True]];
    If[TreeMasses`IsFermion[exchangeParticle] && TreeMasses`IsScalar[photonEmitter],
       Return[True]];
    
    Return[False];
  ]

AMuonCreateMuonPhysicalMass[] := "return context.model.get_physical().M" <>
                             CXXDiagrams`CXXNameOfField[AMuonGetMuon[]] <>
                             If[GetCXXMuonIndex[] =!= Null,
                                "( " <> ToString @ GetCXXMuonIndex[] <> " )",
                                ""] <>
                             ";"

AMuonCreateCalculation[gTaggedDiagrams_List] :=
  Module[{muon = AMuonGetMuon[], cxxMuonIndex = GetCXXMuonIndex[],
          calculation,numberOfIndices},
    numberOfIndices = CXXDiagrams`NumberOfFieldIndices[muon];
    
    "std::array<int, " <> ToString @ numberOfIndices <>
    "> indices = {" <>
      If[cxxMuonIndex =!= Null,
         " " <> ToString @ cxxMuonIndex <>
         If[numberOfIndices =!= 1,
            StringJoin @ Table[", 0", {numberOfIndices-1}],
            ""] <>
         " ",
         If[numberOfIndices =!= 0,
            StringJoin @ Riffle[Table[" 0", {numberOfIndices}], ","] <> " ",
            ""]
        ] <> "};\n\n" <>
                                 
    StringJoin @ Riffle[Module[{graph = #[[1]], diagrams = #[[2]]},
			StringJoin @ Riffle[Module[{diagram = #, indexedDiagram},
				indexedDiagram = CXXDiagrams`IndexDiagramFromGraph[diagram, graph];
				
				"val += " <> 
				ToString @ ProjectColourFactor[
					CXXDiagrams`ColourFactorForIndexedDiagramFromGraph[indexedDiagram, graph]] <>
				" * " <> 
				CXXEvaluatorForDiagramFromGraph[diagram, graph] <>
				"::value(indices, context);"
			] & /@ diagrams, "\n"]
		] & /@ gTaggedDiagrams, "\n"]
  ];

ProjectColourFactor[colourFactor_] := (
	Utils`AssertWithMessage[!TreeMasses`ColorChargedQ[AMuonGetMuon[]],
		"AMuon::ProjectColourFactor[]: The muon has a colour charge (unsupported)"];
	colourFactor)

CXXEvaluatorForDiagramFromGraph[diagram_, vertexCorrectionGraph] := 
  Module[{photonEmitter, exchangeParticle},
    photonEmitter = diagram[[4,3]]; (* Edge between vertices 4 and 6 (3rd edge of vertex 4) *)
    exchangeParticle = diagram[[4,2]]; (* Edge between vertices 4 and 5 (2nd edge of vertex 4) *)
    
    If[diagram[[6]] =!= {TreeMasses`GetPhoton[],CXXDiagrams`LorentzConjugate[photonEmitter],photonEmitter},
       Return["(unknown diagram)"]];
    If[TreeMasses`IsFermion[photonEmitter] && TreeMasses`IsScalar[exchangeParticle],
       Return[CXXEvaluatorFS[photonEmitter,exchangeParticle]]];
    If[TreeMasses`IsFermion[exchangeParticle] && TreeMasses`IsScalar[photonEmitter],
       Return[CXXEvaluatorSF[photonEmitter,exchangeParticle]]];
    
    Return["(unknown diagram)"];
  ]

CXXEvaluatorFS[photonEmitter_,exchangeParticle_] :=
  "AMuonVertexCorrectionFS<" <>
  CXXDiagrams`CXXNameOfField[photonEmitter] <> ", " <>
  CXXDiagrams`CXXNameOfField[exchangeParticle] <> ">"

CXXEvaluatorSF[photonEmitter_,exchangeParticle_] :=
  "AMuonVertexCorrectionSF<" <> 
  CXXDiagrams`CXXNameOfField[photonEmitter] <> ", " <>
  CXXDiagrams`CXXNameOfField[exchangeParticle] <> ">"
  
GetMinMass[particle_] :=
    Module[{dim = TreeMasses`GetDimension[particle],
            mStr = CConversion`ToValidCSymbolString[FlexibleSUSY`M[particle]],
            tail},
           If[dim == 1,
              "model.get_" <> mStr <> "()",
              tail = ToString[GetDimension[particle] - GetDimensionStartSkippingGoldstones[particle] + 1];
              "model.get_" <> mStr <> "().tail<" <> tail <> ">().minCoeff()"
             ]
          ];

AMuonGetMSUSY[] :=
    Module[{susyParticles},
           susyParticles = Select[TreeMasses`GetSusyParticles[], IsElectricallyCharged];
           If[susyParticles === {},
              "return 0.;",
              "return Min(" <>
                 StringJoin[Riffle[GetMinMass /@ susyParticles, ", "]] <>
              ");"
             ]
          ];

End[];
EndPackage[];


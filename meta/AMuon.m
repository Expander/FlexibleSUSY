(* ::Package:: *)

BeginPackage["AMuon`", {"SARAH`", "CXXDiagrams`", "TextFormatting`", "TreeMasses`", "LoopMasses`"}];

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

ContributingGraphs[] := contributingGraphs

GetPhoton[] := SARAH`Photon
GetMuon[] := If[TreeMasses`GetDimension[TreeMasses`GetSMMuonLeptonMultiplet[]] =!= 1,
                TreeMasses`GetSMMuonLeptonMultiplet[],
                Cases[SARAH`ParticleDefinitions[FlexibleSUSY`FSEigenstates],
                      {p_, {Description -> "Muon", ___}} -> p, 1][[1]]
               ]
GetMuonIndex[] := If[TreeMasses`GetDimension[TreeMasses`GetSMMuonLeptonMultiplet[]] =!= 1,
                     2,
                     Null]

ContributingDiagramsForGraph[graph_] :=
  Module[{diagrams},
    diagrams = CXXDiagrams`FeynmanDiagramsOfType[graph,
         {1 -> GetMuon[], 2 -> SARAH`AntiField[GetMuon[]], 3 -> GetPhoton[]}];
         
    Select[diagrams,IsDiagramSupported[graph,#] &]
 ]
 
IsDiagramSupported[vertexCorrectionGraph,diagram_] :=
  Module[{photonEmitter,exchangeParticle},
    photonEmitter = diagram[[4,3]]; (* Edge between vertices 4 and 6 (3rd edge of vertex 4) *)
    exchangeParticle = diagram[[4,2]]; (* Edge between vertices 4 and 5 (2nd edge of vertex 4) *)
    
    If[diagram[[6]] =!= {GetPhoton[],CXXDiagrams`LorentzConjugate[photonEmitter],photonEmitter},
       Return[False]];
    If[TreeMasses`IsFermion[photonEmitter] && TreeMasses`IsScalar[exchangeParticle],
       Return[True]];
    If[TreeMasses`IsFermion[exchangeParticle] && TreeMasses`IsScalar[photonEmitter],
       Return[True]];
    
    Return[False];
  ]
 
CreateCalculation[gTaggedDiagrams_List] :=
  Module[{muon = GetMuon[], muonIndex = GetMuonIndex[],
          calculation,numberOfIndices},
    numberOfIndices = CXXDiagrams`NumberOfFieldIndices[muon];
    
    "EvaluationContext context{ model };\n" <>
    "std::array<unsigned, " <> ToString @ numberOfIndices <>
    "> indices = {" <>
      If[muonIndex =!= Null,
         " " <> ToString @ muonIndex <>
         If[numberOfIndices =!= 1,
            StringJoin @ Table[", 1", {numberOfIndices-1}],
            ""] <>
         " ",
         If[numberOfIndices =!= 0,
            StringJoin @ Riffle[Table[" 1", {numberOfIndices}], ","] <> " ",
            ""]
        ] <> "};\n\n" <>
                                 
    "double val = 0.0;\n\n" <>
                   
    StringJoin @ Riffle[("val += " <> ToString @ # <> "::value(indices, context);") & /@ 
      Flatten[CXXEvaluatorsForDiagramsFromGraph[#[[2]],#[[1]]] & /@ gTaggedDiagrams],
                                       "\n"] <> "\n\n" <>
                   
    "return val;"
  ];

CXXEvaluatorsForDiagramsFromGraph[diagrams_,graph_] :=
  CXXEvaluatorForDiagramFromGraph[#,graph] & /@ diagrams
CXXEvaluatorForDiagramFromGraph[diagram_,vertexCorrectionGraph] := 
  Module[{photonEmitter,exchangeParticle},
    photonEmitter = diagram[[4,3]]; (* Edge between vertices 4 and 6 (3rd edge of vertex 4) *)
    exchangeParticle = diagram[[4,2]]; (* Edge between vertices 4 and 5 (2nd edge of vertex 4) *)
    
    If[diagram[[6]] =!= {GetPhoton[],CXXDiagrams`LorentzConjugate[photonEmitter],photonEmitter},
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

GetMSUSY[] :=
    Module[{susyParticles},
           susyParticles = Select[TreeMasses`GetSusyParticles[], IsElectricallyCharged];
           If[susyParticles === {},
              "return 0.;",
              "return Min(" <>
                 StringJoin[Riffle[GetMinMass /@ susyParticles, ", "]] <>
              ");"
             ]
          ];

GetQED2L[] :=
    "const double MSUSY = Abs(get_MSUSY(context.model));\n" <>
    "const double m_muon = muonPhysicalMass(context);\n" <>
    "const double alpha_em = Sqr(muonCharge(context))/(4*Pi);\n" <>
    "const double qed_2L = alpha_em/(4*Pi) * 16 * FiniteLog(m_muon/MSUSY);\n\n" <>
    "return qed_2L;";

EndPackage[];


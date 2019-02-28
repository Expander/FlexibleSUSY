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

BeginPackage["FFVFormFactors`",
  {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "CXXDiagrams`", "ColorMathInterface`"}
];

FFVFormFactorsCreateInterfaceFunction::usage = "";
FFVGraphs::usage = "";
FFVContributingDiagramsForGraph::usage = "";
ffff::usage = "";

Begin["Private`"];

vertexCorrectionGraph = {{0,0,0,1,0,0},
                         {0,0,0,0,1,0},
                         {0,0,0,0,0,1},
                         {1,0,0,0,1,1},
                         {0,1,0,1,0,1},
                         {0,0,1,1,1,0}};
contributingGraphs = {vertexCorrectionGraph};
FFVGraphs[] := contributingGraphs;

FFVContributingDiagramsForGraph[graph_, Fj_ -> {Fi_, V_}] :=
   Module[{diagrams},
      diagrams =
         CXXDiagrams`FeynmanDiagramsOfType[
            graph,
            {1 -> Fj, 2 -> CXXDiagrams`LorentzConjugate[Fi], 3 -> CXXDiagrams`LorentzConjugate[V]}
         ];
      Print["From FeynmanDiagramsOfType", diagrams];
         
      Select[diagrams, IsDiagramSupported[graph, #]&]
   ];

IsDiagramSupported[vertexCorrectionGraph, diagram_] :=
  Module[{photonEmitter,exchangeParticle},
    photonEmitter = diagram[[4,3]]; (* Edge between vertices 4 and 6 (3rd edge of vertex 4) *)
    exchangeParticle = diagram[[4,2]]; (* Edge between vertices 4 and 5 (2nd edge of vertex 4) *)
    
    (*If[diagram[[6]] =!= {TreeMasses`GetPhoton[],CXXDiagrams`LorentzConjugate[photonEmitter],photonEmitter},*)
       (*Return[False]];*)
    If[TreeMasses`IsFermion[photonEmitter] && TreeMasses`IsScalar[exchangeParticle],
       Return[True]];
    If[TreeMasses`IsFermion[exchangeParticle] && TreeMasses`IsScalar[photonEmitter],
       Return[True]];

    Return[False];
  ]

MyReIm[z_] := If[$VersionNumber >= 10.1,
   ReIm[z],
   {Re[z], Im[z]}
];

FFVFormFactorsCreateInterfaceFunction::field = "Field `1` is not Gluon or Photon.";

(* given external particles Fj, Fi & V and a list of lists like
{
   particles in the loop
   {Fv, Hpm},
   {
      color factors
      {1, 1},
      vertices
      {{bar[Fv], Fe, conj[Hpm]}, {bar[Fe], Fv, Hpm},  {conj[Hpm], Hpm, VP}, {bar[Fv], Fv, VP}}
   }
}

*)
FFVFormFactorsCreateInterfaceFunction[Fj_ -> {Fi_, V_}, topology_, diagrams_] :=
   Module[{prototype, definition,
           numberOfIndices1 = CXXDiagrams`NumberOfFieldIndices[Fj],
           numberOfIndices2 = CXXDiagrams`NumberOfFieldIndices[Fi],
           numberOfIndices3 = CXXDiagrams`NumberOfFieldIndices[V]},

      Print["Topology ", topology];
      Print["Diagrams ", diagrams];
      prototype =
         "std::valarray<std::complex<double>> calculate_" <> CXXNameOfField[Fj] <>
            "_" <> CXXNameOfField[Fi] <> "_" <> CXXNameOfField[V] <> "_form_factors (\n" <>
            IndentText[
               If[TreeMasses`GetDimension[Fj] =!= 1,
                  "int generationIndex1, ",
                  " "
               ] <>
               If[TreeMasses`GetDimension[Fi] =!= 1,
                  "int generationIndex2,\n",
                  " "
               ] <>
               "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model)"
            ];
                 
      definition =
         prototype <> " {\n\n" <>
            IndentText[
               FlexibleSUSY`FSModelName <> "_mass_eigenstates model_ = model;\n" <>
               "context_base context {model_};\n" <>
               "std::array<int, " <> ToString @ numberOfIndices1 <> "> indices1 = {" <>
                     (* TODO: Specify indices correctly *)
                       If[TreeMasses`GetDimension[Fj] =!= 1,
                          "generationIndex1" <>
                          If[numberOfIndices1 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices1-1}],
                             ""] <> " ",
                          If[numberOfIndices1 =!= 0,
                             StringJoin @ Riffle[Table[" 0", {numberOfIndices1}], ", "],
                             ""]
                         ] <> "};\n" <>
                   "std::array<int, " <> ToString @ numberOfIndices2 <>
                     "> indices2 = {" <>
                       If[TreeMasses`GetDimension[Fi] =!= 1,
                          "generationIndex2" <>
                          If[numberOfIndices2 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices2-1}],
                             ""] <> " ",
                          If[numberOfIndices2 =!= 0,
                             StringJoin @ Riffle[Table[" 0", {numberOfIndices2}], ", "] <> " ",
                             ""]
                         ] <> "};\n\n" <>

               "std::valarray<std::complex<double>> val {0.0, 0.0, 0.0, 0.0};\n\n" <>

                  (*
               Switch[V,
                  SARAH`Photon,
                  StringJoin[(
                     (* emit photon from the scalar *)
                     If[IsElectricallyCharged[#[[1,2]]],
                        CreateCall[#[[2,1,1]], "FFVEmitterS", Fj, Fi, V, #[[1,1]], #[[1,2]]],
                        ""
                     ] <>
                     (* emit photon from the fermion *)
                     If[IsElectricallyCharged[#[[1,1]]],
                        CreateCall[#[[2,1,2]], "FFVEmitterF", Fj, Fi, V, #[[1,1]], #[[1,2]]],
                        ""
                     ]
                  )& /@ gTaggedDiagrams
                  ],
                  SARAH`Gluon,
                  StringJoin[(
                     (* emit gluon from the scalar *)
                     If[ColorChargedQ[#[[1,2]]],
                        CreateCall[#[[2,1,1]], "FFVEmitterS", Fj, Fi, V, #[[1,1]], #[[1,2]]],
                        ""
                     ] <>
                     (* emit gluon from the fermion *)
                     If[ColorChargedQ[#[[1,1]]],
                        CreateCall[#[[2,1,2]], "FFVEmitterF", Fj, Fi, V, #[[1,1]], #[[1,2]]],
                        ""
                     ]
                  )& /@ gTaggedDiagrams
                  ],
                  (* we assume that there are no unbroken gauge groups other than U(1)_em and SU(3)_C *)
                  _, Message[FFVFormFactorsCreateInterfaceFunction::field, V]; Abort[];
               ] <>
                  *)

               "return val;"
            ] <> "\n}\n\n";

      {prototype <> ";", definition}
   ];

CreateCall[color_, type_, Fj_, Fi_, V_, F_, S_] :=
   "val += std::complex<double> " <> (ToString @ N[MyReIm@ColorN[color], 16]) <> " * " <> type <> "<" <>
                     StringJoin @ Riffle[CXXDiagrams`CXXNameOfField /@ {Fj, Fi, V, F, S}, ","]  <>
                     ">::value(indices1, indices2, context);\n";

(* create a list of insertion that will get passed to the FFVCreateInterface *)
ffff[inFermion_ -> {outFermion_, spectator_}] :=
   Module[{scalars, fermions, internalParticles = {}, temp},

      scalars = getParticlesOfType[TreeMasses`IsScalar];
      fermions = getParticlesOfType[TreeMasses`IsFermion];

      Map[
         (temp = singleDiagram[inFermion, outFermion, spectator, #[[1]], #[[2]]];
         If[temp =!= {},
            AppendTo[internalParticles, {#, temp}]
         ])&,
         Tuples[{fermions, scalars}]
      ];

      internalParticles
   ];

(* if a diagram exists, return a color factor and a list of particles in vertices,
   otherwise return an empty list *)
singleDiagram[inFermion_, outFermion_, spectator_, F_?TreeMasses`IsFermion, S_?TreeMasses`IsScalar] :=
   Module[{FBarFjSBar, FiBarFS, SBarSVBar, FBarFVBar, v1, v2, v3, v4,colorIndexAssociation, p, sortColorFacRep,
   colorFacwithSEmit, colorFacWithFEmit},

      On[Assert];
      Assert[IsMassless[spectator]];

      (* if the electric charge of an incoming particle doesn't equal to the sum of charges of outgoing ones,
         return an {} *)
      If[TreeMasses`GetElectricCharge[inFermion] =!= Plus @@ (TreeMasses`GetElectricCharge /@ {S,F}),
         Return[{}]
      ];
      (* if an incoming fermion is color charged, at least on of the particles in the loop has to be colored *)
      If[TreeMasses`ColorChargedQ[inFermion] && !(TreeMasses`ColorChargedQ[F] || TreeMasses`ColorChargedQ[S]),
         Return[{}]
      ];

      colorIndexAssociation =
         If[TreeMasses`ColorChargedQ[#],
            If[TreeMasses`GetDimension[#] === 1,
               #[{Unique["ct"]}],
               #[{Unique["gt"], Unique["ct"]}]
            ], #
         ]& /@ {outFermion, S, F, F, S, inFermion, spectator};

      p = colorIndexAssociation;
      p = {p[[6]], p[[1]], p[[7]], p[[4]], p[[3]], p[[5]], p[[2]]};

      (* ColorMath package distinquishes between first and second index in the delta function *)
      v1 = {CXXDiagrams`LorentzConjugate[F], inFermion, CXXDiagrams`LorentzConjugate[S]};
      FBarFjSBar = SARAHToColorMathSymbols @ SARAH`Vertex[
         {CXXDiagrams`LorentzConjugate[p[[6]]], CXXDiagrams`LorentzConjugate[p[[4]]], p[[1]]}
      ];

      v2 = {CXXDiagrams`LorentzConjugate[outFermion], F, S};
      FiBarFS = SARAHToColorMathSymbols @ SARAH`Vertex[
         {p[[7]], CXXDiagrams`LorentzConjugate[p[[2]]], p[[5]]}
      ];

      v3 = {CXXDiagrams`LorentzConjugate[S], S, CXXDiagrams`LorentzConjugate[spectator]};
      SBarSVBar = SARAHToColorMathSymbols @ SARAH`Vertex[
         {CXXDiagrams`LorentzConjugate[p[[6]]], p[[7]], CXXDiagrams`LorentzConjugate[p[[3]]]}
      ];

      v4 = {CXXDiagrams`LorentzConjugate[F], F, CXXDiagrams`LorentzConjugate[spectator]};
      FBarFVBar = SARAHToColorMathSymbols @ SARAH`Vertex[
         {CXXDiagrams`LorentzConjugate[p[[5]]], p[[4]], CXXDiagrams`LorentzConjugate[p[[3]]]}
      ];
      sortColorFacRep = SortColorDeltas @@ p;

      If[vertexNonZero[FBarFjSBar] && vertexNonZero[FiBarFS],
         If[vertexNonZeroS[SBarSVBar] && !vertexNonZero[FBarFVBar],
            Return[{
               {StripSU3Generators[p[[1]], p[[2]], p[[3]],
                  ColorMath`CSimplify[
                     CalculateColorFactor[{FBarFjSBar, FiBarFS, SBarSVBar}//.sortColorFacRep]
                     ConnectColorLines[p[[5]], p[[4]]]//.sortColorFacRep
                  ]
                  ],
                  0},
               (*{v1, v2, v3}*){v1, v2, v3, v4}}
            ]
         ];
         If[vertexNonZero[FBarFVBar] && !vertexNonZeroS[SBarSVBar],
            Return[{
                  {0,
                     StripSU3Generators[p[[1]], p[[2]], p[[3]],
                        ColorMath`CSimplify[
                           CalculateColorFactor[{FBarFjSBar, FiBarFS, FBarFVBar}//.sortColorFacRep]
                           ConnectColorLines[p[[7]], p[[6]]]//.sortColorFacRep
                        ]
                     ]
                  },
                  (*{v1,v2, v4}*){v1, v2, v3, v4}
               }
            ]
         ];
         If[vertexNonZero[FBarFVBar] && vertexNonZeroS[SBarSVBar],
            Return[
               {
                  StripSU3Generators[p[[1]], p[[2]], p[[3]], #]& /@ {
                     ColorMath`CSimplify[
                        CalculateColorFactor[{FBarFjSBar, FiBarFS, SBarSVBar}//.sortColorFacRep]
                        ConnectColorLines[p[[5]], p[[4]]]//.sortColorFacRep
                     ],
                     ColorMath`CSimplify[
                        CalculateColorFactor[{FBarFjSBar, FiBarFS, FBarFVBar}//.sortColorFacRep]
                        ConnectColorLines[p[[7]], p[[6]]]//.sortColorFacRep
                     ]
                  },
                  {v1, v2, v3, v4}
               }
            ]
         ],
         Return[{}]
      ];

      Return[{}];
   ];
(* TODO: add other topologies? *)

ColorFactorForDiagramFromGraph[graph_, diagram_] := 1;

End[];
EndPackage[];

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
   {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "CXXDiagrams`", "ColorMathInterface`", "Utils`"}
];

FFVFormFactorsCreateInterfaceFunction::usage = "";
FFVGraphs::usage = "";
FFVContributingDiagramsForGraph::usage = "";

Begin["Private`"];

vertexCorrectionGraph = {
   {0, 1, 0, 0, 0, 0},
   {1, 0, 1, 0, 1, 0},
   {0, 1, 0, 0, 1, 1},
   {0, 0, 0, 0, 1, 0},
   {0, 1, 1, 1, 0, 0},
   {0, 0, 1, 0, 0, 0}
};
contributingGraphs = {vertexCorrectionGraph};
FFVGraphs[] := contributingGraphs;

FFVContributingDiagramsForGraph[graph_, Fj_ -> {Fi_, V_}] :=
   Module[{diagrams},
      diagrams =
         CXXDiagrams`FeynmanDiagramsOfType[
            graph,
            {4 -> Fj, 1 -> CXXDiagrams`LorentzConjugate[Fi], 6 -> CXXDiagrams`LorentzConjugate[V]}
         ];

      (* for now we only support selected classes of diagrams *)
      Select[diagrams, IsDiagramSupported[graph, #]&]
   ];

IsDiagramSupported[vertexCorrectionGraph, diagram_] :=
   Module[{photonEmitter, exchangeParticle},

      photonEmitter = diagram[[3,2]]; (* Edge between vertices ? and ? (3rd edge of vertex 4) *)
      exchangeParticle = diagram[[2,3]]; (* Edge between vertices ? and ? (2nd edge of vertex 4) *)

      Print["emit: ", photonEmitter, ", exchange: ", exchangeParticle];
    
      (*If[diagram[[6]] =!= {TreeMasses`GetPhoton[],CXXDiagrams`LorentzConjugate[photonEmitter],photonEmitter},*)
       (*Return[False]];*)
      If[TreeMasses`IsFermion[photonEmitter] && TreeMasses`IsScalar[exchangeParticle],
         Return[True]];
      If[TreeMasses`IsFermion[exchangeParticle] && TreeMasses`IsScalar[photonEmitter], Return[True]];

      Return[False];
   ];

FFVFormFactorsCreateInterfaceFunction::field = "Field `1` is not Gluon or Photon.";

FFVFormFactorsCreateInterfaceFunction[Fj_ -> {Fi_, V_}, topologies_, diagrams_] :=
   Module[{prototype, definition,
           numberOfIndices1 = CXXDiagrams`NumberOfFieldIndices[Fj],
           numberOfIndices2 = CXXDiagrams`NumberOfFieldIndices[Fi],
           numberOfIndices3 = CXXDiagrams`NumberOfFieldIndices[V], temp = {}},

      Utils`AssertWithMessage[Length[topologies] === Length[diagrams],
         "Length of diagrams should be the same as length of topologies"];

      Print[CXXDiagrams`IndexDiagramFromGraph[diagrams[[1,1]], topologies[[1]]]];
      Print[topologies];
      Print[diagrams];
      For[i = 1, i <= Length[topologies], i++,
         For[j = 1, j <= Length[diagrams[[i]]], j++,
         Print[
            ColourFactorForIndexedDiagramFromGraph[
               CXXDiagrams`IndexDiagramFromGraph[diagrams[[i,j]], topologies[[i]]], topologies[[i]]
            ]
         ];
               If[TreeMasses`IsScalar[diagrams[[i,j,3,2]]],
                 temp = temp <> "val += std::complex<double> " <> (ToString @ N[FSReIm@
               ColourFactorForIndexedDiagramFromGraph[
               CXXDiagrams`IndexDiagramFromGraph[diagrams[[i,j]], topologies[[i]]], topologies[[i]]
                  ]
                     ]) <> " * FFVEmitterS<" <>
                        StringJoin@Riffle[CXXDiagrams`CXXNameOfField /@ {Fi, Fj, V, diagrams[[i,j,2,3]], diagrams[[i,j,3,2]]}, ", "] <>
">::value(indices1, indices2, context);\n",
                  ""
               ] <>
                  If[TreeMasses`IsFermion[diagrams[[i,j,3,2]]],
                     temp = temp <> "val += std::complex<double> " <> (ToString @ N[FSReIm@
                        ColourFactorForIndexedDiagramFromGraph[
                           CXXDiagrams`IndexDiagramFromGraph[diagrams[[i,j]], topologies[[i]]], topologies[[i]]
                        ]
                     ]) <> " * FFVEmitterF<" <>
                        StringJoin@Riffle[CXXDiagrams`CXXNameOfField /@ {Fi, Fj, V, diagrams[[i,j,3,2]], diagrams[[i,j,2,3]]},
                           ", "]
                     <> ">::value(indices1, indices2, context);\n",
                     ""
                  ]

            ];
         ];
      Print[temp];
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
               "context_base context {model};\n" <>
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

               temp <>
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
   "val += std::complex<double> " <> (ToString @ N[FSReIm@ColorN[color], 16]) <> " * " <> type <> "<" <>
                     StringJoin @ Riffle[CXXDiagrams`CXXNameOfField /@ {Fj, Fi, V, F, S}, ","]  <>
                     ">::value(indices1, indices2, context);\n";

End[];
EndPackage[];

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
  {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "CXXDiagrams`"}
];

FFVFormFactorsCreateInterfaceFunctionForField::usage="";
FFVFormFactorsContributingDiagramsForFieldAndGraph::usage="";
FFVFormFactorsContributingGraphs::usage="";
FFVFormFactorsCreateInterfaceFunctionForLeptonPair::usage="";

Begin["Private`"];

FFVFormFactorsCreateInterfaceFunctionForLeptonPair[inFermion_, outFermion_, spectator_, gTaggedDiagrams_List] :=
   Module[
      {prototype, definition, numberOfIndices1 = CXXDiagrams`NumberOfFieldIndices[inFermion],
         numberOfIndices2 = CXXDiagrams`NumberOfFieldIndices[outFermion],
         numberOfIndices3 = CXXDiagrams`NumberOfFieldIndices[spectator]
      },

      prototype =
         "std::valarray<std::complex<double>> calculate_" <> CXXNameOfField[inFermion] <>
            "_" <> CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[spectator] <> "_form_factors" <>
            " (\n" <>
            If[TreeMasses`GetDimension[inFermion] =!= 1,
               "   int generationIndex1, ",
               " "
            ] <>
            If[TreeMasses`GetDimension[outFermion] =!= 1,
               " int generationIndex2, ",
               " "
            ] <>
            "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model )";
                 
      definition =
          prototype <> "{\n" <>
            IndentText[
               FlexibleSUSY`FSModelName <> "_mass_eigenstates model_ = model;\n" <>
               "EvaluationContext context{ model_ };\n" <>
               "std::array<int, " <> ToString @ numberOfIndices1 <> "> indices1 = {" <>
                     (* TODO: Specify indices correctly *)
                       If[TreeMasses`GetDimension[inFermion] =!= 1,
                          " generationIndex1" <>
                          If[numberOfIndices1 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices1-1}],
                             ""] <> " ",
                          If[numberOfIndices1 =!= 0,
                             StringJoin @ Riffle[Table[" 0", {numberOfIndices1}], ","] <> " ",
                             ""]
                         ] <> "};\n" <>
                   "std::array<int, " <> ToString @ numberOfIndices2 <>
                     "> indices2 = {" <>
                       If[TreeMasses`GetDimension[outFermion] =!= 1,
                          " generationIndex2" <>
                          If[numberOfIndices2 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices2-1}],
                             ""] <> " ",
                          If[numberOfIndices2 =!= 0,
                             StringJoin @ Riffle[Table[" 0", {numberOfIndices2}], ","] <> " ",
                             ""]
                         ] <> "};\n\n" <>

               "std::valarray<std::complex<double>> val {0.0, 0.0, 0.0, 0.0};\n\n" <>

                  StringJoin[
                    If[spectator === SARAH`Photon,
                       (* emit photon from the scalar *)
                      If[IsElectricallyCharged[#[[1,2]]],
                    ("val += std::complex<double> " <> (ToString @ N[ReIm@ColorN[#[[2,1,1]]], 16]) <> " * FFVEmitterS<" <>
                     StringRiffle[CXXDiagrams`CXXNameOfField /@ {inFermion, outFermion, spectator, #[[1,1]], #[[1,2]]}, ","]  <>
                     ">::value(indices1, indices2, context);\n"), ""] <>
                      (* emit photon from the fermion *)
                          If[IsElectricallyCharged[#[[1,1]]],
                            ("val += std::complex<double> " <> (ToString @ N[ReIm@ColorN[#[[2,1,2]]], 16]) <> " * FFVEmitterF<" <>
                                StringRiffle[CXXDiagrams`CXXNameOfField /@ {inFermion, outFermion, spectator, #[[1,1]], #[[1,2]]}, ","]  <>
                                ">::value(indices1, indices2, context);\n"), ""],
                       ""
                      ]& /@ gTaggedDiagrams
                  ] <> "\n" <>

                   StringJoin[
                      If[spectator === SARAH`Gluon,
                         If[ColorChargedQ[#[[1,2]]],
                            ("val += std::complex<double> " <> (ToString @ N[ReIm@ColorN[#[[2,1,1]]], 16]) <> " * FFVEmitterS<" <>
                                StringRiffle[CXXDiagrams`CXXNameOfField /@ {inFermion, outFermion, spectator, #[[1,1]], #[[1,2]]}, ","]  <>
                                ">::value(indices1, indices2, context);\n"), ""] <>
                             If[ColorChargedQ[#[[1,1]]],
                                ("val += std::complex<double> " <> (ToString @ N[ReIm@ColorN[#[[2,1,2]]], 16]) <> " * FFVEmitterF<" <>
                                    StringRiffle[CXXDiagrams`CXXNameOfField /@ {inFermion, outFermion, spectator, #[[1,1]], #[[1,2]]}, ","]  <>
                                    ">::value(indices1, indices2, context);\n"), ""],
                         ""
                      ]& /@ gTaggedDiagrams
                   ] <> "\n" <>

               "return val;"
            ] <> "\n}\n\n";

    {prototype <> ";", definition}
  ];

(* TODO: add other topologies? *)

End[];
EndPackage[];

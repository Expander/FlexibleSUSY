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

BeginPackage["FToFConversionInNucleus`", {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "CXXDiagrams`"}];

FToFConversionInNucleusCreateInterfaceFunctionForField::usage="";
FToFConversionInNucleusContributingDiagramsForFieldAndGraph::usage="";
FToFConversionInNucleusContributingGraphs::usage="";

(* TODO: uncomment this in the end *)
(*Begin["Private`"];*)

FToFConversionInNucleus`FToFConversionInNucleusCreateInterface[{inFermion_, outFermion_, spectator_}] :=
  Module[
    {prototype, definition, 
      nucleus = Au,
     numberOfIndices1 = CXXDiagrams`NumberOfFieldIndices[inFermion],
     numberOfIndices2 = CXXDiagrams`NumberOfFieldIndices[outFermion],
     numberOfIndices3 = CXXDiagrams`NumberOfFieldIndices[spectator]},
   
      prototype =
          "double calculate_" <> CXXNameOfField[inFermion] <> "_to_" <>
            CXXNameOfField[outFermion] <> "_in_nucleus (" <>
            If[TreeMasses`GetDimension[inFermion] =!= 1,
               "int generationIndex1, ",
               " "
            ] <>
            If[TreeMasses`GetDimension[outFermion] =!= 1,
               " int generationIndex2, ",
               " "
            ] <>
            "const " <> FlexibleSUSY`FSModelName <> "_f_to_f_conversion::Nucleus nucleus, " <>
            "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model)";
                 
      definition = 
        prototype <> " {\n" <>
        IndentText[
FlexibleSUSY`FSModelName <> "_mass_eigenstates model_ = model;\n" <>
                  "EvaluationContext context{ model_ };\n" <>
                  "std::array<int, " <> ToString @ numberOfIndices1 <>
                     "> indices1 = {" <>
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
                    "const auto form_factors = calculate_" <> CXXNameOfField[inFermion] <> "_"
                   <> CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[VP] <> "_form_factors "<>
                   "(" <> If[TreeMasses`GetDimension[inFermion] =!= 1,
                            "generationIndex1, ",
                            " "] <>
                         If[TreeMasses`GetDimension[outFermion] =!= 1,
                            " generationIndex2, ",
                            " "] <>
                  "model);\n" <>
                  "const auto nuclear_form_factors = get_overlap_integrals(flexiblesusy::" <> ToString[FlexibleSUSY`FSModelName] <> "_f_to_f_conversion::Nucleus::" <> ToString[nucleus] <> ");\n" <>
                  "constexpr auto GF {1.1667e-5};\n" <>
                     "auto A2L {form_factors[2]};" <>
   "auto A2R {form_factors[3]};" <>
   "// translate from the convention of Hisano, Moroi & Tobe to Kitano, Koike & Okada" <>
   "A2L = A2L/(4.*GF/sqrt(2.));" <>
   "A2R = A2R/(4.*GF/sqrt(2.));" <>
   "const auto left {A2L*nuclear_form_factors.D};" <>
   "const auto right {A2R*nuclear_form_factors.D};" <>
        "// eq. 14 of Kitano, Koike and Okada" <>
        "return 2.*pow(GF,2)*(std::norm(left) + std::norm(right));\n"
        ] <>
        "}\n";
    
    {prototype <> ";", definition}
  ];
EndPackage[];

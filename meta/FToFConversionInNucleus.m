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

FToFConversionInNucleus`FToFConversionInNucleusCreateInterfaceFunctionForLeptonPair[{inFermion_, outFermion_, spectator_}] :=
  Module[
    {prototype, definition, 
     numberOfIndices1 = CXXDiagrams`NumberOfFieldIndices[inFermion],
     numberOfIndices2 = CXXDiagrams`NumberOfFieldIndices[outFermion],
     numberOfIndices3 = CXXDiagrams`NumberOfFieldIndices[spectator]},
   
    Print["Inside of FToFConvertionInNucles"];
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
            "const " <> FlexibleSUSY`FSModelName <> "_f_to_f_conversion::Nucleus, " <>
            "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model)";
                 
      definition = 
        prototype <> "{\n" <>
        "return 0.;\n" <>
        "}\n";
    
    {prototype <> ";", definition}
  ];
EndPackage[];

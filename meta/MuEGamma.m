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

BeginPackage["MuEGamma`", 
   {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "CXXDiagrams`"}
];

MuEGammaCreateInterfaceFunctionForLeptonPair::usage = "";

Begin["`Private`"];

MuEGammaCreateInterfaceFunctionForLeptonPair[{inFermion_, outFermion_, spectator_}] :=
    Module[{prototype, definition,
            numberOfIndices1 = CXXDiagrams`NumberOfFieldIndices[inFermion],
            numberOfIndices2 = CXXDiagrams`NumberOfFieldIndices[outFermion],
            numberOfIndices3 = CXXDiagrams`NumberOfFieldIndices[spectator]},

        prototype =
            "double calculate_" <> CXXNameOfField[inFermion] <> "_to_" <>
            CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[spectator] <> "(" <>
            If[TreeMasses`GetDimension[inFermion] =!= 1,
               " int generationIndex1, ",
               " "
            ] <>
            If[TreeMasses`GetDimension[outFermion] =!= 1,
               " int generationIndex2, ",
               " "
            ] <>
            "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Physical_input& physical_input);";

                 definition =
            (* calculate observable using formfactors *)
            "double calculate_" <> CXXNameOfField[inFermion] <> "_to_" <> CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[spectator] <> " (" <>
               IndentText[
                  If[TreeMasses`GetDimension[inFermion] =!= 1, "int generationIndex1, ", " "] <>
                  If[TreeMasses`GetDimension[outFermion] =!= 1, "int generationIndex2, ", " "] <> "\n" <>
                  "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Physical_input& physical_input) {\n\n"
               ] <>
            (* choose which observable to compute from form factors *)
            Switch[{inFermion, spectator}, {SARAH`Fe, SARAH`VP},
               (* write routine for mu to e gamma *)
               IndentText[
                  FlexibleSUSY`FSModelName <> "_mass_eigenstates model_ = model;\n" <>
                  "context_base context{ model_ };\n" <>
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
                   <> CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[spectator] <> "_form_factors "<>
                   "(" <> If[TreeMasses`GetDimension[inFermion] =!= 1,
                            "generationIndex1, ",
                            ""] <>
                         If[TreeMasses`GetDimension[outFermion] =!= 1,
                            "generationIndex2, ",
                            ""] <>
                  "model);\n" <>
                  (* Dominik suggest that the phase space prefactor should use pole masses  so we get them from the input file *)
                  "double leptonInMassOS;\n" <>
                  "switch (generationIndex1) {\n" <> 
                  IndentText[
                     "case 0: leptonInMassOS = qedqcd.displayMass(softsusy::mElectron); break;\n" <> 
                     "case 1: leptonInMassOS = qedqcd.displayMass(softsusy::mMuon);     break;\n" <> 
                     "case 2: leptonInMassOS = qedqcd.displayMass(softsusy::mTau);      break;\n" <> 
                     "default: exit;\n"
                  ] <>
                  "}\n" <>
                  "\n" <>
                  "// eq. 51 of arXiv:hep-ph/9510309 (note that we include 'e' in the definition of form_factor)\n" <>
                  "const double partial_width = pow(leptonInMassOS,5)/(16.0*Pi) * (std::norm(form_factors[2]) + std::norm(form_factors[3]));\n" <>

                  "const double total_width = lepton_total_decay_width(indices1, indices2, model, qedqcd);\n" <>
                  "return partial_width/total_width;\n"
               ], {SARAH`Fd, SARAH`VP},
               (* write routine for b -> s gamma *)
               IndentText[
                  FlexibleSUSY`FSModelName <> "_mass_eigenstates model_ = model;\n" <>
                  "context_base context{ model_ };\n" <>
                  "std::array<int, " <> ToString @ numberOfIndices1 <>
                     "> indices1 = {" <>
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
                   <> CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[spectator] <> "_form_factors "<>
                   "(" <> If[TreeMasses`GetDimension[inFermion] =!= 1,
                            " generationIndex1, ",
                            " "] <>
                         If[TreeMasses`GetDimension[outFermion] =!= 1,
                            " generationIndex2, ",
                            " "] <>
                  "model );\n" <>
                  "std::valarray<std::complex<double>> c7NP {0.0, 0.0};\n" <>
                  "c7NP[0] = -1/(2*unit_charge(context)) * std::complex<double>(form_factors[3]);\n" <>
                  "c7NP[1] = -1/(2*unit_charge(context)) * std::complex<double>(form_factors[2]);\n" <>
                  "return std::real(c7NP[0]);\n"
               ], {SARAH`Fd, SARAH`VG},
               (* write routine for b -> s gluon *)
               IndentText[
                  FlexibleSUSY`FSModelName <> "_mass_eigenstates model_ = model;\n" <>
                  "context_base context{ model_ };\n" <>
                  "std::array<int, " <> ToString @ numberOfIndices1 <>
                     "> indices1 = {" <>
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
                   <> CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[spectator] <> "_form_factors "<>
                   "(" <> If[TreeMasses`GetDimension[inFermion] =!= 1,
                            " generationIndex1, ",
                            " "] <>
                         If[TreeMasses`GetDimension[outFermion] =!= 1,
                            " generationIndex2, ",
                            " "] <>
                  "model );\n" <>
                  "std::valarray<std::complex<double>> c8NP {0.0, 0.0};\n" <>
                  (* TODO: Get correct g3 *)
                  "const double g3 = 1.04407069; \n" <>
                  "c8NP[0] = -1/(2*g3) * std::complex<double>(form_factors[3]);\n" <>
                  "c8NP[1] = -1/(2*g3) * std::complex<double>(form_factors[2]);\n" <>
                  "return std::real(c8NP[0]);\n"
               ],
            _, ""

            ] <> "}";

        {prototype, definition}
    ];

End[];
EndPackage[];

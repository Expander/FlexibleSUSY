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

BeginPackage["FToFConversionInNucleus`",
    {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "CXXDiagrams`"}
];

FToFConversionInNucleusCreateInterface::usage = "";

Begin["Private`"];

FToFConversionInNucleusCreateInterface[{inFermion_, outFermion_, nucleus_}] :=
    Module[{prototype, definition,
            numberOfIndices1 = CXXDiagrams`NumberOfFieldIndices[inFermion],
            numberOfIndices2 = CXXDiagrams`NumberOfFieldIndices[outFermion]},

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
                "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)";

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

                "const auto photon_exchange = calculate_" <> CXXNameOfField[inFermion] <> "_" <>
                    CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[SARAH`VP] <> "_form_factors (" <>
                    If[TreeMasses`GetDimension[inFermion] =!= 1, "generationIndex1, ", " "] <>
                    If[TreeMasses`GetDimension[outFermion] =!= 1, " generationIndex2, ", " "] <>
                    "model);\n" <>

                "const auto nuclear_form_factors = get_overlap_integrals(flexiblesusy::" <>
                    ToString[FlexibleSUSY`FSModelName] <> "_f_to_f_conversion::Nucleus::" <> SymbolName[nucleus] <>
                    ", qedqcd" <> ");\n" <>

                "// get Fermi constant from Les Houches input file\n" <>
                "const auto GF {qedqcd.displayFermiConstant()};\n" <>
                (* TODO: replace with the value from qedqcd *)
                "const auto e = sqrt(1./137. * 4. * Pi);\n" <>

                "\n// translate from the convention of Hisano, Moroi & Tobe to Kitano, Koike & Okada\n" <>
                (* TODO: check the statement below *)
                "// Hisano defines form factors A2 through a matrix element in eq. 14\n" <>
                "// Kitano uses a lagrangian with F_munu. There is a factor of 2 from translation\n" <>
                "// because Fmunu = qeps - eps q\n" <>
                "const auto A2L = 0.5 * photon_exchange[2]/(4.*GF/sqrt(2.));\n" <>
                "const auto A2R = 0.5 * photon_exchange[3]/(4.*GF/sqrt(2.));\n" <>

                "\n// construct 4-fermion operators from A1 form factors\n" <>
                "const auto gpLV = photon_exchange[0]/(GF/sqrt(2.0)) * e * (2.*2./3.-1./3.);\n" <>
                "const auto gpRV = photon_exchange[1]/(GF/sqrt(2.0)) * e * (2.*2./3.-1./3.);\n" <>
                "const auto gnLV = photon_exchange[0]/(GF/sqrt(2.0)) * e * (2.*(-1./3.)+2./3.);\n" <>
                "const auto gnRV = photon_exchange[1]/(GF/sqrt(2.0)) * e * (2.*(-1./3.)+2./3.);\n" <>

                "\nconst auto left {A2L*nuclear_form_factors.D + gpLV*nuclear_form_factors.Vp + gnLV*nuclear_form_factors.Vn};\n" <>
                "const auto right {A2R*nuclear_form_factors.D + gpRV*nuclear_form_factors.Vp + gnRV*nuclear_form_factors.Vn};\n" <>

                "\n// eq. 14 of Kitano, Koike and Okada\n" <>
                "return 2.*pow(GF,2)*(std::norm(left) + std::norm(right));\n"
            ] <>
            "}\n";
    
        {prototype <> ";", definition}
    ];

End[];
EndPackage[];

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

BeginPackage["Himalaya`", {"SARAH`", "TextFormatting`", "CConversion`",
                           "Parameters`", "TreeMasses`", "Utils`"}];

FillHimalayaInput::usage = "Fills parameters into Himalaya input
 struct.";

{ RenormalizationScheme, Lambda3LEFT, Lambda3LH3m, Lambda3LUncertainty,
  \[Mu], g3, vu, vd, MSQ2, MSD2, MSU2, At, Ab, mG, mW, mZ, Global`mt, mb, mA };

Begin["`Private`"];

AppendAtEnd[str_, char_] :=
    StringReplace[str, EndOfLine|EndOfString -> char];

FillHimalayaInput[inputPars_, struct_:"pars"] :=
    Module[{result = "",
            inPars = Parameters`DecreaseIndexLiterals[inputPars],
            pars = First /@ Reverse /@ inputPars
           },
           result = Parameters`CreateLocalConstRefs[pars] <> "

pars.scale = MODELPARAMETER(scale);
pars.mu = Re(" <> CConversion`RValueToCFormString[\[Mu] /. inPars] <> ");
pars.g3 = Re(" <> CConversion`RValueToCFormString[g3 /. inPars] <> ");
pars.vd = Re(" <> CConversion`RValueToCFormString[vd /. inPars] <> ");
pars.vu = Re(" <> CConversion`RValueToCFormString[vu /. inPars] <> ");
pars.mq2 = Re(" <> CConversion`RValueToCFormString[MSQ2 /. inPars] <> ");
pars.md2 = Re(" <> CConversion`RValueToCFormString[MSD2 /. inPars] <> ");
pars.mu2 = Re(" <> CConversion`RValueToCFormString[MSU2 /. inPars] <> ");
pars.At = Re(" <> CConversion`RValueToCFormString[At /. inPars] <> ");
pars.Ab = Re(" <> CConversion`RValueToCFormString[Ab /. inPars] <> ");
pars.MG = " <> CConversion`RValueToCFormString[mG /. inPars] <> ";
pars.MW = " <> CConversion`RValueToCFormString[mW /. inPars] <> ";
pars.MZ = " <> CConversion`RValueToCFormString[mZ /. inPars] <> ";
pars.Mt = " <> CConversion`RValueToCFormString[Global`mt /. inPars] <> ";
pars.Mb = " <> CConversion`RValueToCFormString[mb /. inPars] <> ";
pars.MA = " <> CConversion`RValueToCFormString[mA /. inPars] <> ";

const double msbar_scheme = " <> ToString[If[(RenormalizationScheme /. inPars) === DRbar, 0, 1]] <> ";
const double lambda_3L_eft = " <> CConversion`RValueToCFormString[Lambda3LEFT /. inPars] <> ";
const double lambda_3L_h3m = " <> CConversion`RValueToCFormString[Lambda3LH3m /. inPars] <> ";
const double lambda_3L_uncertainty = " <> CConversion`RValueToCFormString[Lambda3LUncertainty /. inPars] <> ";
";
           AppendAtEnd[result, " \\"]
          ];

End[];

EndPackage[];

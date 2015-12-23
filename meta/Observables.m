BeginPackage["Observables`", {"FlexibleSUSY`", "Utils`", "TextFormatting`"}];

(* observables *)
Begin["FlexibleSUSYObservable`"];
FSObservables = { aMuon, aMuonGM2Calc, aMuonGM2CalcUncertainty };
End[];

CalculateObservables::usage="";

Begin["`Private`"];

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuon, structName_String] :=
    structName <> ".AMU = " <> FlexibleSUSY`FSModelName <> "_GMuonMinus2::calculate_amuon(MODEL);";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc, structName_String] :=
    structName <> ".AMUGM2CALC = gm2calc_calculate_amu(gm2calc_data);";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty, structName_String] :=
    structName <> ".AMUGM2CALCUNCERTAINTY = gm2calc_calculate_amu_uncertainty(gm2calc_data);";

FillGM2CalcInterfaceData[struct_String] :=
    Module[{filling, mwStr,
            w, pseudoscalar, smuon, muonsneutrino, chargino, neutralino,
            mu, m1, m2, m3, mq2, mu2, md2, ml2, me2, tu, td, te, yu, yd, ye},
           w             = Parameters`GetParticleFromDescription["W-Boson"];
           pseudoscalar  = Parameters`GetParticleFromDescription["Pseudo-Scalar Higgs"];
           smuon         = Parameters`GetParticleFromDescription["Smuon"];
           muonsneutrino = Parameters`GetParticleFromDescription["Muon Sneutrino"];
           chargino      = Parameters`GetParticleFromDescription["Charginos"];
           neutralino    = Parameters`GetParticleFromDescription["Neutralinos"];
           mu            = Parameters`GetParameterFromDescription["Mu-parameter"];
           m1            = Parameters`GetParameterFromDescription["Bino Mass parameter"];
           m2            = Parameters`GetParameterFromDescription["Wino Mass parameter"];
           m3            = Parameters`GetParameterFromDescription["Gluino Mass parameter"];
           mq2           = Parameters`GetParameterFromDescription["Softbreaking left Squark Mass"];
           mu2           = Parameters`GetParameterFromDescription["Softbreaking right Up-Squark Mass"];
           md2           = Parameters`GetParameterFromDescription["Softbreaking right Down-Squark Mass"];
           ml2           = Parameters`GetParameterFromDescription["Softbreaking left Slepton Mass"];
           me2           = Parameters`GetParameterFromDescription["Softbreaking right Slepton Mass"];
           tu            = Parameters`GetParameterFromDescription["Trilinear-Up-Coupling"];
           td            = Parameters`GetParameterFromDescription["Trilinear-Down-Coupling"];
           te            = Parameters`GetParameterFromDescription["Trilinear-Lepton-Coupling"];
           yu            = Parameters`GetParameterFromDescription["Up-Yukawa-Coupling"];
           yd            = Parameters`GetParameterFromDescription["Down-Yukawa-Coupling"];
           ye            = Parameters`GetParameterFromDescription["Lepton-Yukawa-Coupling"];
           mwStr         = "MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[w]];
           filling = \
           struct <> ".alpha_s_MZ = ALPHA_S_MZ;\n" <>
           struct <> ".MZ    = MZPole;\n" <>
           "if (!is_zero(" <> mwStr <> "))\n" <>
              TextFormatting`IndentText[struct <> ".MW = " <> mwStr <> ";"] <> "\n" <>
           "else if (!is_zero(MWPole))\n" <>
              TextFormatting`IndentText[struct <> ".MW = MWPole;"] <> "\n" <>
           struct <> ".mb_mb = MBMB;\n" <>
           struct <> ".MT    = MTPole;\n" <>
           struct <> ".MTau  = MTauPole;\n" <>
           struct <> ".MM    = MMPole;\n" <>
           struct <> ".MA0   = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[pseudoscalar][1]] <> ";\n" <>
           struct <> ".MSvm  = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[muonsneutrino]] <> ";\n" <>
           struct <> ".MSm   = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[smuon]] <> ";\n" <>
           struct <> ".MCha  = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[chargino]] <> ";\n" <>
           struct <> ".MChi  = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[neutralino]] <> ";\n" <>
           struct <> ".scale = MODEL.get_scale();\n" <>
           struct <> ".TB    = MODEL.get_" <> CConversion`RValueToCFormString[SARAH`VEVSM2] <> "() / " <>
                              "MODEL.get_" <> CConversion`RValueToCFormString[SARAH`VEVSM1] <> "();\n" <>
           struct <> ".Mu    = MODEL.get_" <> CConversion`RValueToCFormString[mu] <> "();\n" <>
           struct <> ".M1    = MODEL.get_" <> CConversion`RValueToCFormString[m1] <> "();\n" <>
           struct <> ".M2    = MODEL.get_" <> CConversion`RValueToCFormString[m2] <> "();\n" <>
           struct <> ".M3    = MODEL.get_" <> CConversion`RValueToCFormString[m3] <> "();\n" <>
           struct <> ".mq2   = MODEL.get_" <> CConversion`RValueToCFormString[mq2] <> "();\n" <>
           struct <> ".mu2   = MODEL.get_" <> CConversion`RValueToCFormString[mu2] <> "();\n" <>
           struct <> ".md2   = MODEL.get_" <> CConversion`RValueToCFormString[md2] <> "();\n" <>
           struct <> ".ml2   = MODEL.get_" <> CConversion`RValueToCFormString[ml2] <> "();\n" <>
           struct <> ".me2   = MODEL.get_" <> CConversion`RValueToCFormString[me2] <> "();\n" <>
           struct <> ".Au    = div_save(MODEL.get_" <> CConversion`RValueToCFormString[tu] <>
                               "(), MODEL.get_" <> CConversion`RValueToCFormString[yu] <> "());\n" <>
           struct <> ".Ad    = div_save(MODEL.get_" <> CConversion`RValueToCFormString[td] <>
                               "(), MODEL.get_" <> CConversion`RValueToCFormString[yd] <> "());\n" <>
           struct <> ".Ae    = div_save(MODEL.get_" <> CConversion`RValueToCFormString[te] <>
                               "(), MODEL.get_" <> CConversion`RValueToCFormString[ye] <> "());\n";
           "GM2Calc_data " <> struct <> ";\n" <> filling
          ];

FillInterfaceData[{}] := "";

FillInterfaceData[obs_List] :=
    Module[{},
           If[MemberQ[obs,FlexibleSUSYObservable`aMuonGM2Calc] ||
              MemberQ[obs,FlexibleSUSYObservable`aMuonGM2CalcUncertainty],
              FillGM2CalcInterfaceData["gm2calc_data"]
             ]
          ];

CalculateObservables[something_, structName_String] :=
    Module[{observables},
           observables = Cases[something, a_?(MemberQ[FlexibleSUSYObservable`FSObservables,#]&) :> a, {0, Infinity}];
           FillInterfaceData[observables] <> "\n" <>
           Utils`StringJoinWithSeparator[CalculateObservable[#,structName]& /@ observables, "\n"]
          ];

End[];

EndPackage[];

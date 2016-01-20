BeginPackage["Observables`", {"FlexibleSUSY`", "SARAH`", "TreeMasses`", "Utils`", "CConversion`", "TextFormatting`"}];

(* observables *)
Begin["FlexibleSUSYObservable`"];
FSObservables = { aMuonGM2Calc, aMuonGM2CalcUncertainty,
                  CpHiggsPhotonPhoton, CpHiggsGluonGluon,
                  CpPseudoScalarPhotonPhoton, CpPseudoScalarGluonGluon };
End[];

GetRequestedObservables::usage="";
GetNumberOfObservables::usage="";
CreateObservablesDefinitions::usage="";
CreateObservablesInitialization::usage="";
CreateSetAndDisplayObservablesFunctions::usage="";
CreateClearObservablesFunction::usage="";
CalculateObservables::usage="";

Begin["`Private`"];

(* @todo proper error handling, e.g. when no Higgs or pseudo-scalar defined *)
(* @todo properly count complex parameters *)

GetRequestedObservables[blocks_] :=
    DeleteDuplicates[Cases[blocks, a_?(MemberQ[FlexibleSUSYObservable`FSObservables,#]&) :> a, {0, Infinity}]];

GetObservableName[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc] := "a_muon_gm2_calc";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty] := "a_muon_gm2_calc_uncertainty";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`CpHiggsPhotonPhoton] := "eff_cp_higgs_photon_photon";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`CpHiggsGluonGluon] := "eff_cp_higgs_gluon_gluon";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] := "eff_cp_pseudoscalar_photon_photon";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarGluonGluon] := "eff_cp_pseudoscalar_gluon_gluon";

GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc] := "a_muon = (g-2)/2 of the muon (calculated with GM2Calc)";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty] := "uncertainty of (g-2)/2 of the muon (calculated with GM2Calc)";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`CpHiggsPhotonPhoton] := "effective H-Photon-Photon coupling";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`CpHiggsGluonGluon] := "effective H-Gluon-Gluon coupling";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] := "effective A-Photon-Photon coupling";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarGluonGluon] := "effective A-Gluon-Gluon coupling";

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty] := CConversion`ScalarType[CConversion`realScalarCType];

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`CpHiggsPhotonPhoton] :=
    Module[{dim, type},
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson];
           If[dim == 1,
              type = CConversion`ScalarType[CConversion`complexScalarCType],
              type = CConversion`ArrayType[CConversion`complexScalarCType, dim]
             ];
           type
          ];

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`CpHiggsGluonGluon] :=
    Module[{dim, type},
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson];
           If[dim == 1,
              type = CConversion`ScalarType[CConversion`complexScalarCType],
              type = CConversion`ArrayType[CConversion`complexScalarCType, dim]
             ];
           type
          ];

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] :=
    Module[{dim, type},
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar];
           If[dim == 1,
              type = CConversion`ScalarType[CConversion`complexScalarCType],
              type = CConversion`ArrayType[CConversion`complexScalarCType, dim]
             ];
           type
          ];

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarGluonGluon] :=
    Module[{dim, type},
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar];
           If[dim == 1,
              type = CConversion`ScalarType[CConversion`complexScalarCType],
              type = CConversion`ArrayType[CConversion`complexScalarCType, dim]
             ];
           type
          ];

GetNumberOfObservables[observables_List] :=
    Module[{i, number = 0},
           For[i = 1, i <= Length[observables], i++,
               Which[observables[[i]] === FlexibleSUSYObservable`aMuonGM2Calc ||
                     observables[[i]] === FlexibleSUSYObservable`aMuonGM2CalcUncertainty,
                     number = number + 1,
                     observables[[i]] === FlexibleSUSYObservable`CpHiggsPhotonPhoton ||
                     observables[[i]] === FlexibleSUSYObservable`CpHiggsGluonGluon,
                     number = number + TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson],
                     observables[[i]] === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton ||
                     observables[[i]] === FlexibleSUSYObservable`CpPseudoScalarGluonGluon,
                     number = number + TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar],
                     True,
                     Print["Warning: ignoring invalid observable ", observables[[i]]];
                    ];
              ];
           number
          ];

CreateObservablesDefinitions[observables_List] :=
    Module[{i, type, name, description, dim, definitions = ""},
           For[i = 1, i <= Length[observables], i++,
               If[MemberQ[FlexibleSUSYObservable`FSObservables, observables[[i]]],
                  name = GetObservableName[observables[[i]]];
                  description = GetObservableDescription[observables[[i]]];
                  type = CConversion`CreateCType[GetObservableType[observables[[i]]]];
                  definitions = definitions <> type <> " " <> name <> "; ///< " <> description <> "\n";,
                  Print["Warning: ignoring invalid observable ", observables[[i]]];
              ];
           definitions
          ];

CreateObservablesInitialization[observables_List] :=
    Module[{i, name, value, init = ": "},
           For[i = 1, i <= Length[observables], i}},
               name = GetObservableName[observables[[i]]];
               Which[observables[[i]] === FlexibleSUSYObservable`aMuonGM2Calc ||
                     observables[[i]] === FlexibleSUSYObservable`aMuonGM2CalcUncertainty,
                     value = "(0)",
                     observables[[i]] === FlexibleSUSYObservable`CpHiggsPhotonPhoton ||
                     observables[[i]] === FlexibleSUSYObservable`CpHiggsGluonGluon,
                     dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson];
                     If[dim == 1,
                        value = "(0,0)";,
                        value = "(" <> CConversion`CreateCType[CConversion`ArrayType[CConversion`complexScalarCType, dim]]
                                <> "::Zero())";
                       ];,
                     observables[[i]] === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton ||
                     observables[[i]] === FlexibleSUSYObservable`CpPseudoScalarGluonGluon,
                     dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar];
                     If[dim == 1,
                        value = "(0,0)";,
                        value = "(" <> CConversion`CreateCType[CConversion`ArrayType[CConversion`complexScalarCType, dim]]
                                <> "::Zero())";
                       ];,
                     True,
                     Print["Warning: ignoring invalid observable ", observables[[i]]];
                     ];
               If[init == ": ", init = init <> name <> value <> "\n", init = init <> ", " <> name <> value <> "\n"];
              ];
           init
          ];

CreateSetAndDisplayObservablesFunctions[observables_List] :=
    Module[{i, name, dim, display = "", displayNames = "", set = ""},
           For[i = 1, i <= Length[observables], i++,
               If[MemberQ[FlexibleSUSYObservable`FSObservables, observables[[i]]],
                  name = GetObservableName[observables[[i]]];
                  Which[observables[[i]] === FlexibleSUSYObservable`aMuonGM2Calc ||
                        observables[[i]] === FlexibleSUSYObservable`aMuonGM2CalcUncertainty,
                        display = display <> "vec(" <> ToString[i-1] <> ") = " <> name <> ";\n";
                        displayNames = displayNames <> "names[" <> ToString[i-1] <> "] = \""
                                       <> name <> "\";\n";
                        set = set <> name <> " = vec(" <> ToString[i-1] <> ");\n";,
                        observables[[i]] === FlexibleSUSYObservable`CpHiggsPhotonPhoton ||
                        observables[[i]] === FlexibleSUSYObservable`CpHiggsGluonGluon,
                        dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson];
                        If[dim == 1,
                           (* continue here *)
                       ];
                  ,
                  Print["Warning: ignoring invalid observable ", observables[[i]]];
                 ];
              ];
           {display, displayNames, set}
          ];

CreateClearObservablesFunction[observables_List] :=
    Module[{i, name, value, dim, result = ""},
           For[i = 1, i <= Length[observables], i++,
               If[MemberQ[FlexibleSUSYObservable`FSObservables, observables[[i]]],
                  name = GetObservableName[observables[[i]]];
                  Switch[observables[[i]],
                         FlexibleSUSYObservable`aMuonGM2Calc, value = "0.",
                         FlexibleSUSYObservable`aMuonGM2CalcUncertainty, value = "0.",
                         FlexibleSUSYObservable`CpHiggsPhotonPhoton,
                         dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson];
                         If[dim == 1,
                            value = "(0.,0.)",
                            value = CConversion`CreateCType[CConversion`ArrayType[CConversion`complexScalarCType, dim]]
                                    <> "::Zero()"
                           ];,
                         FlexibleSUSYObservable`CpHiggsGluonGluon,
                         dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson];
                         If[dim == 1,
                            value = "(0.,0.)",
                            value = CConversion`CreateCType[CConversion`ArrayType[CConversion`complexScalarCType, dim]]
                                    <> "::Zero()"
                           ];,
                         FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton,
                         dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar];
                         If[dim == 1,
                            value = "(0.,0.)",
                            value = CConversion`CreateCType[CConversion`ArrayType[CConversion`complexScalarCType, dim]]
                                    <> "::Zero()"
                           ];,
                         FlexibleSUSYObservable`CpPseudoScalarGluonGluon,
                         dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar];
                         If[dim == 1,
                            value = "(0.,0.)",
                            value = CConversion`CreateCType[CConversion`ArrayType[CConversion`complexScalarCType, dim]]
                                    <> "::Zero()"
                           ];
                        ];
                  result = result <> name <> " = " <> value <> ";\n" ,
                  Print["Warning: ignoring invalid observable ", observables[[i]]];
                 ];
              ];
          ];

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc, structName_String] :=
    structName <> ".AMUGM2CALC = gm2calc_calculate_amu(gm2calc_data);";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty, structName_String] :=
    structName <> ".AMUGM2CALCUNCERTAINTY = gm2calc_calculate_amu_uncertainty(gm2calc_data);";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`CpHiggsPhotonPhoton, structName_String] :=
    structName <> ".EFFCPHIGGSPHOTONPHOTON = 0.;";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`CpHiggsGluonGluon, structName_String] :=
    structName <> ".EFFCPHIGGSGLUONGLUON = 0.;";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton, structName_String] :=
    structName <> ".EFFCPPSEUDOSCALARPHOTONPHOTON = 0.;";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarGluonGluon, structName_String] :=
    structName <> ".EFFCPPSEUDOSCALARGLUONGLUON = 0.;";

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

FillEffectiveCouplingsInterfaceData[struct_String] :=
    FlexibleSUSY`FSModelName <> "_effective_couplings " <> struct <> "(MODEL);\n";

FillInterfaceData[{}] := "";

FillInterfaceData[obs_List] :=
    Module[{filled = ""},
           If[MemberQ[obs,FlexibleSUSYObservable`aMuonGM2Calc] ||
              MemberQ[obs,FlexibleSUSYObservable`aMuonGM2CalcUncertainty],
              filled = filled <> FillGM2CalcInterfaceData["gm2calc_data"];
             ];
           If[MemberQ[obs,FlexibleSUSYObservable`CpHiggsPhotonPhoton]         ||
              MemberQ[obs,FlexibleSUSYObservable`CpHiggsGluonGluon]           ||
              MemberQ[obs, FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] ||
              MemberQ[obs, FlexibleSUSYObservable`CpPseudoScalarGluonGluon],
              filled = filled <> FillEffectiveCouplingsInterfaceData["effective_couplings"];
             ];
           filled
          ];

CalculateObservables[something_, structName_String] :=
    Module[{observables},
           observables = Cases[something, a_?(MemberQ[FlexibleSUSYObservable`FSObservables,#]&) :> a, {0, Infinity}];
           FillInterfaceData[observables] <> "\n" <>
           Utils`StringJoinWithSeparator[CalculateObservable[#,structName]& /@ observables, "\n"]
          ];

End[];

EndPackage[];

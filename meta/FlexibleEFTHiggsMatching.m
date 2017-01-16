BeginPackage["FlexibleEFTHiggsMatching`", {"CConversion`", "TreeMasses`", "LoopMasses`", "Constraint`", "ThresholdCorrections`", "Parameters`", "Utils`"}];

CalculateMHiggsPoleOneMomentumIteration::usage = "";
CalculateRunningUpQuarkMasses::usage = "";
CalculateRunningDownQuarkMasses::usage = "";
CalculateRunningDownLeptonMasses::usage = "";
FillSMFermionPoleMasses::usage = "";
SetBSMParameters::usage = "";

Begin["`Private`"];

CalculateMHiggsPoleOneMomentumIteration[particle_] :=
    If[GetDimension[particle] == 1,
       "Mh2_pole = mh2_tree - self_energy - tadpole(0);"
       ,
"const auto M_loop = (mh2_tree - self_energy - " <> CreateCType[TreeMasses`GetMassMatrixType[particle]] <> "(tadpole.asDiagonal())).eval();

" <> CreateCType[CConversion`ArrayType[CConversion`realScalarCType, GetDimension[particle]]] <> " M_pole;
fs_diagonalize_hermitian(M_loop, M_pole);

Mh2_pole = M_pole(idx);"
      ];

CalculateRunningUpQuarkMasses[] :=
    Module[{result = "", i, iStr, mq, mqFun},
           For[i = 0, i < 3, i++,
               iStr = ToString[i];
               mq = mqFun = CConversion`RValueToCFormString[TreeMasses`GetUpQuark[i + 1, True]];
               If[Length[TreeMasses`GetSMUpQuarks[]] == 3, mqFun = mq <> "()"];
               result = result <>
                        "upQuarksDRbar(" <> iStr <> "," <> iStr <> ") = " <>
                        "sm.get_physical().MFu(" <> iStr <> ") - " <>
                        "model.get_physical().M" <> mq <> " + model.get_M" <> mqFun <> ";\n";
              ];
           result
          ];

CalculateRunningDownQuarkMasses[] :=
    Module[{result = "", i, iStr, mq, mqFun},
           For[i = 0, i < 3, i++,
               iStr = ToString[i];
               mq = mqFun = CConversion`RValueToCFormString[TreeMasses`GetDownQuark[i + 1, True]];
               If[Length[TreeMasses`GetSMDownQuarks[]] == 3, mqFun = mq <> "()"];
               result = result <>
                        "downQuarksDRbar(" <> iStr <> "," <> iStr <> ") = " <>
                        "sm.get_physical().MFd(" <> iStr <> ") - " <>
                        "model.get_physical().M" <> mq <> " + model.get_M" <> mqFun <> ";\n";
              ];
           result
          ];

CalculateRunningDownLeptonMasses[] :=
    Module[{result = "", i, iStr, mq, mqFun},
           For[i = 0, i < 3, i++,
               iStr = ToString[i];
               mq = mqFun = CConversion`RValueToCFormString[TreeMasses`GetDownLepton[i + 1, True]];
               If[Length[TreeMasses`GetSMChargedLeptons[]] == 3, mqFun = mq <> "()"];
               result = result <>
                        "downLeptonsDRbar(" <> iStr <> "," <> iStr <> ") = " <>
                        "sm.get_physical().MFe(" <> iStr <> ") - " <>
                        "model.get_physical().M" <> mq <> " + model.get_M" <> mqFun <> ";\n";
              ];
           result
          ];

FillSMFermionPoleMasses[] :=
    Module[{result = "", i, mq},
           For[i = 0, i < 3, i++,
               mq = CConversion`RValueToCFormString[TreeMasses`GetUpQuark[i + 1, True]];
               result = result <>
                        "this->model.get_physical().M" <> mq <> " = " <>
                        "eft.get_physical().MFu(" <> ToString[i] <> ");\n";
              ];
           For[i = 0, i < 3, i++,
               mq = CConversion`RValueToCFormString[TreeMasses`GetDownQuark[i + 1, True]];
               result = result <>
                        "this->model.get_physical().M" <> mq <> " = " <>
                        "eft.get_physical().MFd(" <> ToString[i] <> ");\n";
              ];
           For[i = 0, i < 3, i++,
               mq = CConversion`RValueToCFormString[TreeMasses`GetDownLepton[i + 1, True]];
               result = result <>
                        "this->model.get_physical().M" <> mq <> " = " <>
                        "eft.get_physical().MFe(" <> ToString[i] <> ");\n";
              ];
           result
          ];

(* find loop order at which parameter needs to be determined given the
   expression which determines the Higgs mass

   par = parameter

   exprs = list of expressions which determine Mh.
   exprs[[1]] contains the tree-level expression.
   exprs[[2]] contains the 1-loop expression.
   exprs[[3]] contains the 2-loop expression.

   towerLoopOrder = requested loop order of the tower
 *)
FindMHiggsLoopOrderFor[par_, exprs_List, towerLoopOrder_] :=
    Module[{i},
           For[i = 0, i < Length[exprs], i++,
               If[!FreeQ[exprs[[i+1]], par],
                  Return[Max[{0, towerLoopOrder - i}]];
                 ];
              ];
           0
          ];

SetBSMParameterAtLoopOrder[par_, lo_, struct_String] :=
    Module[{parName = CConversion`ToValidCSymbolString[par]},
           struct <> "set_" <> parName <> "(model_" <> ToString[lo] <> "l.get_" <> parName <> "());\n"
          ];

SetBSMParameters[susyScaleMatching_List, higgsMassMatrix_, struct_String:""] :=
    Module[{pars, loopOrder},
           (* BSM parameters fixed by matching SM -> BSM *)
           pars = Intersection[Join[{SARAH`hyperchargeCoupling, SARAH`leftCoupling, SARAH`strongCoupling,
                                     SARAH`UpYukawa, SARAH`DownYukawa, SARAH`ElectronYukawa},
                                    First /@ susyScaleMatching],
                               Parameters`GetModelParameters[]
                              ];
           (* find matching loop orders for parameters for FlexibleEFTHiggs-1L *)
           loopOrder = FindMHiggsLoopOrderFor[#, {higgsMassMatrix, 0}, 1]& /@ pars;
           StringJoin[SetBSMParameterAtLoopOrder[#[[1]], #[[2]], struct]& /@ Utils`Zip[pars, loopOrder]]
          ];

End[];

EndPackage[];

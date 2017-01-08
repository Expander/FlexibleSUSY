BeginPackage["FlexibleEFTHiggsMatching`", {"CConversion`", "TreeMasses`", "LoopMasses`", "Constraint`", "ThresholdCorrections`", "Parameters`"}];

CallSMPoleMassFunctions::usage = "";
CalculateRunningUpQuarkMasses::usage = "";
CalculateRunningDownQuarkMasses::usage = "";
CalculateRunningDownLeptonMasses::usage = "";
FillSMFermionPoleMasses::usage = "";

Begin["`Private`"];

CallSMPoleMassFunctions[states_, enablePoleMassThreads_] :=
    Module[{particles, result, Selector},
           Selector[p_] := SARAH`SMQ[p] && !IsMassless[p] && (IsVector[p] || IsFermion[p]);
           particles = Select[LoopMasses`GetLoopCorrectedParticles[states], Selector];
           If[enablePoleMassThreads =!= True,
              result = StringJoin[LoopMasses`CallPoleMassFunction[#,"model."]& /@ particles];
              ,
              result = StringJoin[LoopMasses`CallThreadedPoleMassFunction[#,"model_ptr"]& /@ particles] <> "\n";
             ];
           result
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

End[];

EndPackage[];

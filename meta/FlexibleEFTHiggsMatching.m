BeginPackage["FlexibleEFTHiggsMatching`", {"CConversion`", "TreeMasses`", "LoopMasses`", "Constraint`", "ThresholdCorrections`", "Parameters`"}];

CallSMPoleMassFunctions::usage = "";
CalculateRunningFermionMasses::usage = "";
FillSMFermionPoleMasses::usage = "";

Begin["`Private`"];

DefineFutureAndCallThreadedPoleMassFunction[particle_Symbol, ptr_:"this"] :=
    Module[{massStr},
           massStr = ToValidCSymbolString[FlexibleSUSY`M[particle]];
          "auto fut_" <> massStr <> " = run_async([obj_ptr] () { obj_ptr->" <>
           CreateLoopMassFunctionName[particle] <> "(); });\n"
          ];

JoinFutures[particle_Symbol] :=
    "fut_" <> ToValidCSymbolString[FlexibleSUSY`M[particle]] <> ".get();\n";

CallSMPoleMassFunctions[states_, enablePoleMassThreads_] :=
    Module[{particles, result, Selector},
           Selector[p_] := SARAH`SMQ[p] && !IsMassless[p] && (IsVector[p] || IsFermion[p]);
           particles = Select[LoopMasses`GetLoopCorrectedParticles[states], Selector];
           If[enablePoleMassThreads =!= True,
              result = StringJoin[LoopMasses`CallPoleMassFunction[#,"model."]& /@ particles];
              ,
              result = StringJoin[DefineFutureAndCallThreadedPoleMassFunction[#,"&model"]& /@ particles] <> "\n" <>
                       StringJoin[JoinFutures /@ particles];
             ];
           result
          ];

CalculateRunningFermionMasses[] :=
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
           For[i = 0, i < 3, i++,
               iStr = ToString[i];
               mq = mqFun = CConversion`RValueToCFormString[TreeMasses`GetDownQuark[i + 1, True]];
               If[Length[TreeMasses`GetSMDownQuarks[]] == 3, mqFun = mq <> "()"];
               result = result <>
                        "downQuarksDRbar(" <> iStr <> "," <> iStr <> ") = " <>
                        "sm.get_physical().MFd(" <> iStr <> ") - " <>
                        "model.get_physical().M" <> mq <> " + model.get_M" <> mqFun <> ";\n";
              ];
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

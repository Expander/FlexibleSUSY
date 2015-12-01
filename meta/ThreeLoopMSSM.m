BeginPackage["ThreeLoopMSSM`", {"SARAH`", "Parameters`"}];

BetaMSSM::usage = "beta functions of the MSSM from SUSYHD v1.0";

Begin["`Private`"];

subDir = FileNameJoin[{FlexibleSUSY`$flexiblesusyMetaDir, "ThreeLoopMSSM"}];

IsDefinedAndEqual[descr_String, c_] :=
    Parameters`GetParameterFromDescription[descr] =!= Null &&
    Parameters`GetParameterFromDescription[descr] === c

G1GUTNormalization[] :=
    Parameters`GetGUTNormalization[SARAH`hyperchargeCoupling] / Sqrt[3/5];

BetaMSSM[gc_] :=
    Switch[gc,
           SARAH`hyperchargeCoupling, Get[FileNameJoin[{subDir, "beta_g1.m"}]] / G1GUTNormalization[],
           SARAH`leftCoupling       , Get[FileNameJoin[{subDir, "beta_g2.m"}]],
           SARAH`strongCoupling     , Get[FileNameJoin[{subDir, "beta_g3.m"}]],
           (* SARAH`UpYukawa           , Get[FileNameJoin[{subDir, "beta_gt.m"  }]], *)
           (* SARAH`DownYukawa         , Get[FileNameJoin[{subDir, "beta_gb.m"  }]], *)
           (* SARAH`ElectronYukawa     , Get[FileNameJoin[{subDir, "beta_gtau.m"}]], *)
           (* \[Lambda]                , Get[FileNameJoin[{subDir, "beta_lambda.m"}]], *)
           (* m2                       , Get[FileNameJoin[{subDir, "beta_m2.m"}]], *)
           _, Which[IsDefinedAndEqual["MSSM Mu Parameter", gc],
                    Get[FileNameJoin[{subDir, "beta_m2.m"}]],
                    IsDefinedAndEqual["MSSM Higgs Selfcouplings", gc],
                    Get[FileNameJoin[{subDir, "beta_lambda.m"}]],
                    True, Print["Error: unknown coupling: ", gc]; {0,0,0}
                    ]
          ] /. ThreeLoopMSSM`ToSARAHNamingConvention[];

(* Note:
   g1, g2, g3, Ye are global variables in SARAH
 *)
ToSARAHNamingConvention[] := {
    g1 -> SARAH`hyperchargeCoupling G1GUTNormalization[],
    g2 -> SARAH`leftCoupling,
    g3 -> SARAH`strongCoupling,
    Global`Yt -> SARAH`UpYukawa,
    Global`Yb -> SARAH`DownYukawa,
    Ye -> SARAH`ElectronYukawa
    (* , *)
    (* \[Lambda] -> Parameters`GetParameterFromDescription["MSSM Higgs Selfcouplings"], *)
    (* m2        -> Parameters`GetParameterFromDescription["MSSM Mu Parameter"] *)
};

End[];

EndPackage[];

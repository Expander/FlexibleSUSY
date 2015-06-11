BeginPackage["ThreeLoopSM`", {"SARAH`", "Parameters`"}];

BetaSM::usage = "beta functions of the SM from SUSYHD v1.0";

Begin["`Private`"];

subDir = FileNameJoin[{Global`$flexiblesusyMetaDir, "ThreeLoopSM"}];

BetaSM[gc_] :=
    Switch[gc,
           SARAH`hyperchargeCoupling, Get[FileNameJoin[{subDir, "beta_g1.m"}]],
           SARAH`leftCoupling       , Get[FileNameJoin[{subDir, "beta_g2.m"}]],
           SARAH`strongCoupling     , Get[FileNameJoin[{subDir, "beta_g3.m"}]],
           SARAH`UpYukawa           , Get[FileNameJoin[{subDir, "beta_gt.m"  }]],
           SARAH`DownYukawa         , Get[FileNameJoin[{subDir, "beta_gb.m"  }]],
           SARAH`ElectronYukawa     , Get[FileNameJoin[{subDir, "beta_gtau.m"}]],
           \[Lambda]                , Get[FileNameJoin[{subDir, "beta_lambda.m"}]],
           _, Print["Error: unknown coupling: ", gc]; {0,0,0}
          ] /. ThreeLoopSM`ToSARAHNamingConvention[];

ToSARAHNamingConvention[] := {
    g1 -> SARAH`hyperchargeCoupling,
    g2 -> SARAH`leftCoupling,
    g3 -> SARAH`strongCoupling,
    Global`gt -> SARAH`UpYukawa[3,3],
    gb -> SARAH`DownYukawa[3,3],
    Global`g\[Tau] -> SARAH`ElectronYukawa[3,3],
    \[Lambda] -> Parameters`GetParameterFromDescription["SM Higgs Selfcouplings"]
};

End[];

EndPackage[];

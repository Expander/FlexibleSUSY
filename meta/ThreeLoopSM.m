BeginPackage["ThreeLoopSM`", {"SARAH`"}];

BetaSM::usage = "gauge coupling beta function";

Begin["`Private`"];

subDir = FileNameJoin[{Global`$flexiblesusyMetaDir, "ThreeLoopSM"}];

BetaSM[gc_] :=
    Switch[gc,
           SARAH`hyperchargeCoupling, Get[FileNameJoin[{subDir, "beta_g1.m"}]],
           SARAH`leftCoupling       , Get[FileNameJoin[{subDir, "beta_g2.m"}]],
           SARAH`strongCoupling     , Get[FileNameJoin[{subDir, "beta_g3.m"}]],
           SARAH`UpYukawa           , Get[FileNameJoin[{subDir, "beta_gt.m"}]],
           SARAH`DownYukawa         , Get[FileNameJoin[{subDir, "beta_gb.m"}]],
           SARAH`ElectronYukawa     , Get[FileNameJoin[{subDir, "beta_gtau.m"}]],
           \[Lambda]                , Get[FileNameJoin[{subDir, "beta_lambda.m"}]],
           _, Print["Error: unknown coupling: ", gc]; Null
          ];

End[];

EndPackage[];

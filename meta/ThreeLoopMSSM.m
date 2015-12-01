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
           SARAH`UpYukawa           , Get[FileNameJoin[{subDir, "beta_Yu.m"}]],
           SARAH`DownYukawa         , Get[FileNameJoin[{subDir, "beta_Yd.m"}]],
           SARAH`ElectronYukawa     , Get[FileNameJoin[{subDir, "beta_Ye.m"}]],
           _, Which[IsDefinedAndEqual["Bino Mass parameter", gc]  , Get[FileNameJoin[{subDir, "beta_M1.m"}]],
                    IsDefinedAndEqual["Wino Mass parameter", gc]  , Get[FileNameJoin[{subDir, "beta_M2.m"}]],
                    IsDefinedAndEqual["Gluino Mass parameter", gc], Get[FileNameJoin[{subDir, "beta_M3.m"}]],
                    True, Print["Error: unknown coupling: ", gc]; {0,0,0}
                    ]
          ] /. ThreeLoopMSSM`ToSARAHNamingConvention[] /. Zeta[s_] :> N[Zeta[s]];

(* Note:
   g1, g2, g3, Ye, M1, M2, M3 are in SARAH` context
 *)
ToSARAHNamingConvention[] := {
    g1 -> SARAH`hyperchargeCoupling G1GUTNormalization[],
    g2 -> SARAH`leftCoupling,
    g3 -> SARAH`strongCoupling,
    Global`Yt -> SARAH`UpYukawa,
    Global`Yb -> SARAH`DownYukawa,
    Ye -> SARAH`ElectronYukawa,
    M1 -> Parameters`GetParameterFromDescription["Bino Mass parameter"],
    M2 -> Parameters`GetParameterFromDescription["Wino Mass parameter"],
    M3 -> Parameters`GetParameterFromDescription["Gluino Mass parameter"],
    Global`ht -> SARAH`TrilinearUp,
    Global`hb -> SARAH`TrilinearDown,
    Global`he -> SARAH`TrilinearLepton
};

End[];

EndPackage[];

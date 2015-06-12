(*
   This script extracts the SM beta-functions from 1303.4364v2.

   Author: Alexander Voigt

   Run it like this:

   math -run "<< extract_SM_beta_functions_m2_from_1303.4364.m; Quit[]"
*)

Get["beta_m2_1303.4364.m"];

rules = {
    al1 -> g1^2,
    al2 -> g2^2,
    al3 -> g3^2,
    at -> gt^2,
    ab -> gb^2,
    atau -> g\[Tau]^2,
    lam -> \[Lambda]/4,
    NR -> 3,
    cR -> 4/3,
    NG -> 3
};

beta = {
    Simplify[2 m2 Coefficient[bms, h, 1] /. rules],
    Simplify[2 m2 Coefficient[bms, h, 2] /. rules],
    Simplify[2 m2 Coefficient[bms, h, 3] /. rules]
};

outputDir = ".";
filename = FileNameJoin[{outputDir, "beta_m2.m"}];

Put[beta, filename];

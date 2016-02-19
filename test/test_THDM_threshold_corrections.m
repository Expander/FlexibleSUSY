Needs["TestSuite`", "TestSuite.m"];
Needs["ThreeLoopSM`", "ThreeLoopSM.m"];

FlexibleSUSY`$flexiblesusyMetaDir = FileNameJoin[{Directory[], "meta"}];

gRules = {
    gbar -> (g2^2 + gY^2)/4,
    gbarm -> (-g2^2 + gY^2)/4,
    g1 -> gY Sqrt[5/3],
    GUTNormalization[g1] -> Sqrt[3/5]
};

quarks  = {AU -> At, AD -> Ab  , hD -> hb  , hU -> ht, Nc -> 3};
leptons = {AU -> 0 , AD -> Atau, hD -> htau, hU -> 0 , Nc -> 1};

(* Eq. (6.16) from arxiv:hep-ph/9307201 *)
Aprime = -Nc/(96 Pi^2 MSUSY^2) (
    hU^2 {{Mu^2, -Mu AU}, {-Mu AU, AU^2}}
    + hD^2 {{AD^2, -Mu AD}, {-Mu AD, Mu^2}}
);

(* vertex corrections from arxiv:hep-ph/9307201 *)
lam1Ver = (Nc/(16 Pi^2 MSUSY^2) (
    AD^2 hD^2 (2 hD^2 - gbar) + Mu^2 hU^2 gbar));

lam2Ver = (Nc/(16 Pi^2 MSUSY^2) (
    AU^2 hU^2 (2 hU^2 - gbar) + Mu^2 hD^2 gbar));

lam3Ver = (Nc/(32 Pi^2 MSUSY^2) (
    Mu^2 (hU^2 - hD^2)^2 + hU^2 hD^2 (AU + AD)^2
    + gbarm ((AD^2 - Mu^2) hD^2 + (AU^2 - Mu^2) hU^2)
));

lam4Ver = (Nc/(32 Pi^2 MSUSY^2) (
    Mu^2 (hU^2 + hD^2)^2 - hU^2 hD^2 (AU + AD)^2
    + g2^2/2 ((AD^2 - Mu^2) hD^2 + (AU^2 - Mu^2) hU^2)
));

lam5Ver = 0;

lam6Ver = (Nc Mu /(32 Pi^2 MSUSY^2) (
    AD hD^2 (gbar - 2 hD^2) - AU hU^2 gbar));

lam7Ver = (Nc Mu /(32 Pi^2 MSUSY^2) (
    AU hU^2 (gbar - 2 hU^2) - AD hD^2 gbar));

(* box corrections from arxiv:hep-ph/9307201 *)
lam1Box = -(Nc/(96 Pi^2 MSUSY^4) (AD^4 hD^4 + Mu^4 hU^4));
lam2Box = -(Nc/(96 Pi^2 MSUSY^4) (AU^4 hU^4 + Mu^4 hD^4));
lam3Box = -(Nc/(96 Pi^2 MSUSY^4) (
      Mu^2 AU^2 hU^4 + Mu^2 AD^2 hD^4 + hU^2 hD^2 (Mu^2 - AU AD)^2));
lam4Box = -(Nc/(96 Pi^2 MSUSY^4) (
      Mu^2 AU^2 hU^4 + Mu^2 AD^2 hD^4 - hU^2 hD^2 (Mu^2 - AU AD)^2));
lam5Box = -(Nc Mu^2/(96 Pi^2 MSUSY^4) (AD^2 hD^4 + AU^2 hU^4));
lam6Box = (Nc Mu /(96 Pi^2 MSUSY^4) (Mu^2 AU hU^4 + AD^3 hD^4));
lam7Box = (Nc Mu /(96 Pi^2 MSUSY^4) (Mu^2 AD hD^4 + AU^3 hU^4));

(* field renormalization from arxiv:hep-ph/9307201 *)
lam1Field = 2 gbar Aprime[[1, 1]];
lam2Field = 2 gbar Aprime[[2, 2]];
lam3Field = -gbarm (Aprime[[1, 1]] + Aprime[[2, 2]]);
lam4Field = -g2^2/2 (Aprime[[1, 1]] + Aprime[[2, 2]]);
lam5Field = 0;
(* the following two terms are 0 in the PhysRev.D.48.4280 *)
lam6Field = gbar/2 Aprime[[1, 2]];
lam7Field = gbar/2 Aprime[[1, 2]];

(* list with vertex + box corrections *)
lamTh = {
   lam1Ver + lam1Box,
   lam2Ver + lam2Box,
   lam3Ver + lam3Box,
   lam4Ver + lam4Box,
   lam5Ver + lam5Box,
   lam6Ver + lam6Box,
   lam7Ver + lam7Box
};

lamTh = (lamTh //. quarks) + (lamTh //. leptons);

(* list with field renormalizations *)
lamPhi = {lam1Field, lam2Field, lam3Field, lam4Field,
          lam5Field, lam6Field, lam7Field};

lamPhi = (lamPhi //. quarks) + (lamPhi //. leptons);

Print["testing THDM threshold corrections[] ..."];

Get["model_files/THDMIIMSSMBC/FlexibleSUSY.m.in"];

renameRules = {
    AtauInput -> Atau,
    AtInput -> At,
    AbInput -> Ab,
    Yu[3,3] -> ht,
    Yd[3,3] -> hb,
    Ye[3,3] -> htau,
    MuInput -> Mu,
    \[Mu] -> Mu
};

deltaLambdaTh = {
    deltaLambda1th1L, deltaLambda2th1L, deltaLambda3th1L, 
    deltaLambda4th1L, deltaLambda5th1L, deltaLambda6th1L, 
    deltaLambda7th1L
} //. renameRules;

deltaLambdaPhi = {
    deltaLambda1Phi1L, deltaLambda2Phi1L, deltaLambda3Phi1L, 
    deltaLambda4Phi1L, deltaLambda5Phi1L, deltaLambda6Phi1L, 
    deltaLambda7Phi1L
} //. renameRules;

lamThDiff  = Simplify[lamTh  - deltaLambdaTh  //. gRules];
lamPhiDiff = Simplify[lamPhi - deltaLambdaPhi //. gRules];

TestEquality[lamThDiff , Table[0, {i,1,7}]];
TestEquality[lamPhiDiff, Table[0, {i,1,7}]];

Get["model_files/HTHDMIIMSSMBC/FlexibleSUSY.m.in"];

deltaLambdaTh = {
    deltaLambda1th1L, deltaLambda2th1L, deltaLambda3th1L, 
    deltaLambda4th1L, deltaLambda5th1L, deltaLambda6th1L, 
    deltaLambda7th1L
} //. renameRules;

deltaLambdaPhi = {
    deltaLambda1Phi1L, deltaLambda2Phi1L, deltaLambda3Phi1L, 
    deltaLambda4Phi1L, deltaLambda5Phi1L, deltaLambda6Phi1L, 
    deltaLambda7Phi1L
} //. renameRules;

lamThDiff  = Simplify[lamTh  - deltaLambdaTh  //. gRules];
lamPhiDiff = Simplify[lamPhi - deltaLambdaPhi //. gRules];

TestEquality[lamThDiff , Table[0, {i,1,7}]];
TestEquality[lamPhiDiff, Table[0, {i,1,7}]];

Get["model_files/HGTHDMIIMSSMBC/FlexibleSUSY.m.in"];

deltaLambdaTh = {
    deltaLambda1th1L, deltaLambda2th1L, deltaLambda3th1L, 
    deltaLambda4th1L, deltaLambda5th1L, deltaLambda6th1L, 
    deltaLambda7th1L
} //. renameRules;

deltaLambdaPhi = {
    deltaLambda1Phi1L, deltaLambda2Phi1L, deltaLambda3Phi1L, 
    deltaLambda4Phi1L, deltaLambda5Phi1L, deltaLambda6Phi1L, 
    deltaLambda7Phi1L
} //. renameRules;

lamThDiff  = Simplify[lamTh  - deltaLambdaTh  //. gRules];
lamPhiDiff = Simplify[lamPhi - deltaLambdaPhi //. gRules];

TestEquality[lamThDiff , Table[0, {i,1,7}]];
TestEquality[lamPhiDiff, Table[0, {i,1,7}]];

PrintTestSummary[];

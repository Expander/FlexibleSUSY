Needs["TestSuite`", "TestSuite.m"];
Needs["EWSB`", "EWSB.m"];

Print["testing MSSM-like EWSB for Mu and BMu ..."];

FlexibleSUSY`FSSolveEWSBTimeConstraint = 120;

mssmEwsbEqs = {
    mu^2 + x^2 + x y + z + 5,
    Bmu  - x^2 + x y + z + 5
};

mssmEwsbOutputParameters = { mu, Bmu };

Parameters`AddRealParameter[mssmEwsbOutputParameters];

TestEquality[EWSB`FindFreePhasesInEWSB[mssmEwsbEqs, mssmEwsbOutputParameters],
             {FlexibleSUSY`Sign[mu]}];

mssmFullSolution = EWSB`Private`FindSolution[mssmEwsbEqs, mssmEwsbOutputParameters];

TestEquality[mssmFullSolution,
             { {B[mu] -> -5 + x^2 - x*y - z, mu -> -Sqrt[-5 - x^2 - x*y - z]}, 
               {B[mu] -> -5 + x^2 - x*y - z, mu -> Sqrt[-5 - x^2 - x*y - z]}
             }];

mssmSigns = Cases[EWSB`Private`FindFreePhase /@ mssmEwsbOutputParameters,
              FlexibleSUSY`Sign[_]];

TestEquality[mssmSigns, {FlexibleSUSY`Sign[mu]}];

TestEquality[EWSB`Private`CanReduceSolution[mssmFullSolution, mssmSigns], True];

mssmReducedSolution = EWSB`Private`ReduceSolution[mssmFullSolution, mssmSigns];

TestEquality[mssmReducedSolution,
             { B[mu] -> -5 + x^2 - x*y - z,
               mu -> LOCALINPUT[FlexibleSUSY`Signmu] Sqrt[-5 - x^2 - x*y - z]
             }];

Print["testing NMSSM-like EWSB for Kappa, vS and mS2 ..."];

muEff = lambda s;

BmuEff = Alambda + kappa s;

nmssmEwsbEqs = {
    vu (mHu2 + muEff^2 + lambda^2 vd^2 + g^2 (vu^2 - vd^2)) - vd muEff BmuEff,
    vd (mHd2 + muEff^2 + lambda^2 vu^2 + g^2 (vd^2 - vu^2)) - vu muEff BmuEff,
    s (mS2 + X + kappa s Akappa + kappa^2 s^2) + Y
};

nmssmEwsbOutputParameters = { s, kappa, mS2 };

Parameters`AddRealParameter[nmssmEwsbOutputParameters];

TestEquality[EWSB`FindFreePhasesInEWSB[nmssmEwsbEqs, nmssmEwsbOutputParameters],
             {}];

nmssmFullSolution = EWSB`Private`FindSolution[nmssmEwsbEqs, nmssmEwsbOutputParameters];

TestEquality[Length[nmssmFullSolution], 2];

PrintTestSummary[];

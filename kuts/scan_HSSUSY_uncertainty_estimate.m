Get["models/HSSUSY/HSSUSY_librarylink.m"];
Get["kuts/HSSUSY_uncertainty_estimate.m"];

settings = {
    precisionGoal -> 1.*^-5,
    maxIterations -> 100,
    poleMassLoopOrder -> 4,
    ewsbLoopOrder -> 4,
    betaFunctionLoopOrder -> 5,
    thresholdCorrectionsLoopOrder -> 4,
    thresholdCorrections -> 124111421
};

smpars = {
    alphaEmMZ -> 0.00775526,(* SMINPUTS[1] *)
    GF -> 1.16639*^-5,      (* SMINPUTS[2] *)
    alphaSMZ -> 0.118,      (* SMINPUTS[3] *)
    MZ -> 91.1876,          (* SMINPUTS[4] *)
    mbmb -> 4.2,            (* SMINPUTS[5] *)
    Mt -> 173.34,           (* SMINPUTS[6] *)
    Mtau -> 1.77703,        (* SMINPUTS[7] *)
    Mv3 -> 0,               (* SMINPUTS[8] *)
    MW -> 80.385,           (* SMINPUTS[9] *)
    Me -> 0.000510998902,   (* SMINPUTS[11] *)
    Mv1 -> 0,               (* SMINPUTS[12] *)
    Mm -> 0.1056583715,     (* SMINPUTS[13] *)
    Mv2 -> 0,               (* SMINPUTS[14] *)
    md2GeV -> 0.00475,      (* SMINPUTS[21] *)
    mu2GeV -> 0.0024,       (* SMINPUTS[22] *)
    ms2GeV -> 0.104,        (* SMINPUTS[23] *)
    mcmc -> 1.27,           (* SMINPUTS[24] *)
    CKMTheta12 -> 0,
    CKMTheta13 -> 0,
    CKMTheta23 -> 0,
    CKMDelta -> 0,
    PMNSTheta12 -> 0,
    PMNSTheta13 -> 0,
    PMNSTheta23 -> 0,
    PMNSDelta -> 0,
    PMNSAlpha1 -> 0,
    PMNSAlpha2 -> 0,
    alphaEm0 -> 1/137.035999074,
    Mh -> 125.09
};

HSSUSYCalcMh[MS_, TB_, Xtt_] :=
    CalcHSSUSYDMh[
        fsSettings -> settings,
        fsSMParameters -> smpars,
        fsModelParameters -> {
            TanBeta -> TB,
            MEWSB -> 173.34,
            MSUSY -> MS,
            M1Input -> MS,
            M2Input -> MS,
            M3Input -> MS,
            MuInput -> MS,
            mAInput -> MS,
            AtInput -> (Xtt + 1/TB) MS,
            AbInput -> 0,
            AtauInput -> 0,
            msq2 -> MS^2 IdentityMatrix[3],
            msu2 -> MS^2 IdentityMatrix[3],
            msd2 -> MS^2 IdentityMatrix[3],
            msl2 -> MS^2 IdentityMatrix[3],
            mse2 -> MS^2 IdentityMatrix[3],
            LambdaLoopOrder -> 2, (* may use 3 here *)
            TwoLoopAtAs -> 1,
            TwoLoopAbAs -> 1,
            TwoLoopAtAb -> 1,
            TwoLoopAtauAtau -> 1,
            TwoLoopAtAt -> 1,
            ThreeLoopAtAsAs -> 1
        }
   ];

LinearRange[start_, stop_, steps_] :=
    Range[start, stop, (stop - start)/steps];

Xtt = -Sqrt[6];
TB  = 20;

data = ParallelMap[
    { N[#], Sequence @@ HSSUSYCalcMh[#, TB, Xtt] }&,
    LogRange[300, 10^5, 4]
];

Export["HSSUSY_uncertainty.dat", data];
Quit[];

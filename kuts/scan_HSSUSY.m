Get["models/HSSUSY/HSSUSY_librarylink.m"];

(* get digit of [num] at position [pos] *)
GetDigit[num_, pos_, base_:10] :=
    IntegerPart[Mod[num / base^pos, base]];

(* set digit of [num] at position [pos] to [val] *)
SetDigit[num_, pos_, val_, base_:10] :=
    num + (val - GetDigit[num,pos,base]) base^pos;

(* generate logarithmically spaced range [start, stop] *)
LogRange[start_, stop_, steps_] :=
    Exp /@ Range[Log[start], Log[stop], (Log[stop] - Log[start])/steps];

(* calculate Higgs mass *)
CalcHSSUSYMh[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    CalcHSSUSYMh[a, Sequence @@ s, r];

CalcHSSUSYMh[ytLoops_?NumericQ, Qpole_?NumericQ, Qm_?NumericQ, eft_?NumericQ, ytMSSM_?NumericQ, args__] :=
    Module[{handle, spec, tc},
           tc = thresholdCorrections /. { args };
           tc = If[IntegerQ[tc], tc,
                   thresholdCorrections /. Options[FSHSSUSYOpenHandle]];
           If[ytLoops >= 0, tc = SetDigit[tc, 6, ytLoops]];
           handle = FSHSSUSYOpenHandle[args];
           FSHSSUSYSet[handle,
               fsSettings -> {
                   calculateStandardModelMasses -> 1,
                   thresholdCorrectionsLoopOrder -> 4,
                   poleMassScale -> Qpole,
                   thresholdCorrections -> tc
               },
               fsModelParameters -> {
                   DeltaEFT -> eft,
                   DeltaYt -> ytMSSM,
                   Qmatch -> Qm
               }
           ];
           spec = FSHSSUSYCalculateSpectrum[handle];
           FSHSSUSYCloseHandle[handle];
           If[spec === $Failed, $Failed,
              Pole[M[hh]] /. (HSSUSY /. spec)]
          ];

(* calculate Higgs mass and uncertainty estimate *)
CalcHSSUSYDMh[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    CalcHSSUSYDMh[a, Sequence @@ s, r];

CalcHSSUSYDMh[args__] :=
    Module[{Mh0, Mh, MhYt3L, MhEFT, MhYtMSSM, varyQpole, varyQmatch,
            DMhSM, DMhEFT, DMhSUSY,
            MS = MSUSY /. { args }, Mlow = MEWSB /. { args }},
           Mh0        = CalcHSSUSYMh[-1, 0, 0, 0, 0, args];
           If[Mh0 === $Failed, Return[{$Failed, $Failed}]];
           Mh         = CalcHSSUSYMh[3, 0, 0, 0, 0, args];
           If[Mh === $Failed, Return[{Mh0, $Failed}]];
           MhYt3L     = CalcHSSUSYMh[4, 0, 0, 0, 0, args];
           If[MhYt3L === $Failed, Return[{Mh0, $Failed}]];
           MhEFT      = CalcHSSUSYMh[3, 0, 0, 1, 0, args];
           If[MhEFT === $Failed, Return[{Mh0, $Failed}]];
           MhYtMSSM   = CalcHSSUSYMh[3, 0, 0, 0, 1, args];
           If[MhYtMSSM === $Failed, Return[{Mh0, $Failed}]];
           varyQpole  = CalcHSSUSYMh[3, #, 0, 0, 0, args]& /@
                        LogRange[Mlow/2, 2 Mlow, 10];
           varyQmatch = CalcHSSUSYMh[3, 0, #, 0, 0, args]& /@
                        LogRange[MS/2, 2 MS, 10];
           varyQpole  = Select[varyQpole , NumericQ];
           varyQmatch = Select[varyQmatch, NumericQ];
           If[varyQmatch === {} || varyQpole === {}, Return[{Mh0, $Failed}]];
           (* combine uncertainty estimates *)
           DMhSM   = Max[Abs[Max[varyQpole] - Mh],
                         Abs[Min[varyQpole] - Mh]] +
                     Abs[Mh - MhYt3L];
           DMhEFT  = Abs[Mh - MhEFT];
           DMhSUSY = Max[Abs[Max[varyQmatch] - Mh],
                         Abs[Min[varyQmatch] - Mh]] +
                     Abs[Mh - MhYtMSSM];
           (* { Mh0, DMhSM + DMhEFT + DMhSUSY } *)
           {
               (* 1 Mh   *) Mh0,
               (* 2 SM   *) Max[Abs[Max[varyQpole] - Mh], Abs[Min[varyQpole] - Mh]],
               (* 3 SM   *) Abs[Mh - MhYt3L],
               (* 4 SUSY *) Max[Abs[Max[varyQmatch] - Mh], Abs[Min[varyQmatch] - Mh]],
               (* 5 SUSY *) Abs[Mh - MhYtMSSM],
               (* 6 EFT  *) DMhEFT
           }
          ];

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

Xtt = -Sqrt[6];
TB  = 20;

data = ParallelMap[
    { N[#], Sequence @@ HSSUSYCalcMh[#, TB, Xtt] }&,
    LogRange[200, 10^5, 100]
];

Export["HSSUSY_uncertainty.dat", data];

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

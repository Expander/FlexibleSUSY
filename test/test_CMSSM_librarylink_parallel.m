Needs["TestSuite`", "TestSuite.m"];
Get["models/CMSSM/CMSSM_librarylink.m"];

Off[FSCMSSM::info];
Off[FSCMSSMCalculateSpectrum::error];
Off[FSCMSSMCalculateSpectrum::warning];
Off[FSCMSSMCalculateObservables::error];
Off[FSCMSSMCalculateObservables::warning];

CalcMh[TB_] :=
    Module[{spec, handle},
           handle = FSCMSSMOpenHandle[
               fsSettings -> {calculateStandardModelMasses -> 1,
                              betaFunctionLoopOrder -> 2},
               fsModelParameters -> {m0 -> 125, m12 -> 500, TanBeta -> TB,
                                     SignMu -> 1, Azero -> 0}];
           spec = FSCMSSMCalculateSpectrum[handle];
           FSCMSSMCloseHandle[handle];
           If[spec === $Failed, 0,
              (Pole[M[hh]] /. spec)[[1]]]
          ];

kernels = LaunchKernels[];
Print["Using ", Length[kernels], " kernels."];

DistributeDefinitions[CalcMh];

range = Range[1, 50, 1];
resultMap = AbsoluteTiming[Map[CalcMh, range]];
resultParallelMap = AbsoluteTiming[ParallelMap[CalcMh, range]];

Print["time for sequential runs: ", resultMap[[1]]];
Print["time for parallel runs  : ", resultParallelMap[[1]]];

If[Length[kernels] > 1,
   TestLowerThan[resultParallelMap[[1]], resultMap[[1]]];
  ];

TestEquality[resultMap[[2]], resultParallelMap[[2]]];

PrintTestSummary[];

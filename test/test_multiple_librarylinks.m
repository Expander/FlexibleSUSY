Needs["TestSuite`", "TestSuite.m"];

Get["models/SM/SM_librarylink.m"];
Get["models/SMtower/SMtower_librarylink.m"];

hSMtower = FSSMtowerOpenHandle[betaFunctionLoopOrder -> 3, Mt -> 173.34, LambdaIN -> 0.2, Qin -> 2000];

TestEquality[LambdaIN              /. FSSMtowerGetInputParameters[hSMtower]  , 0.2];
TestEquality[Mt                    /. FSSMtowerGetSMInputParameters[hSMtower], 173.34];
TestEquality[betaFunctionLoopOrder /. FSSMtowerGetSettings[hSMtower]         , 3];

TestEquality[FSSMGetInputParameters[hSMtower]  , $Failed];
TestEquality[FSSMGetSMInputParameters[hSMtower], $Failed];
TestEquality[FSSMGetSettings[hSMtower]         , $Failed];
TestEquality[FSSMCalculateSpectrum[hSMtower]   , $Failed];
TestEquality[FSSMCalculateObservables[hSMtower], $Failed];

hSM = FSSMOpenHandle[betaFunctionLoopOrder -> 2, Mt -> 173, LambdaIN -> 0.1, Qin -> 1000];

TestEquality[LambdaIN              /. FSSMtowerGetInputParameters[hSMtower]  , 0.2];
TestEquality[Mt                    /. FSSMtowerGetSMInputParameters[hSMtower], 173.34];
TestEquality[betaFunctionLoopOrder /. FSSMtowerGetSettings[hSMtower]         , 3];

TestEquality[LambdaIN              /. FSSMGetInputParameters[hSMtower]  , 0.1];
TestEquality[Mt                    /. FSSMGetSMInputParameters[hSMtower], 173.];
TestEquality[betaFunctionLoopOrder /. FSSMGetSettings[hSMtower]         , 2];

FSSMCloseHandle[hSMtower];
FSSMCloseHandle[hSM];

PrintTestSummary[];

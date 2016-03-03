Needs["TestSuite`", "TestSuite.m"];
Needs["LoopFunctions`", "LoopFunctions.m"];

Print["testing divergences ..."];

loopFunctions = {A0, B0, B1, B22, B22tilde, C0, D0, D27, F, G, H};

For[i = 1, i <= Length[loopFunctions], i++,
    TestEquality[
        Coefficient[loopFunctions[[i]][m,Q] /. LoopFunctions`FullMomentum[], LoopFunctions`Delta],
        Coefficient[loopFunctions[[i]][m,Q] /. LoopFunctions`ZeroMomentum[], LoopFunctions`Delta]
    ];
   ];

Print["testing B0[] function ..."];

TestEquality[Simplify[
    (B0[1, 2, 3, 4] /. LoopFunctions`FullMomentum[]) -
    LoopFunctions`Private`B0integral[1, 2, 3, 4]],
    0
];

PrintTestSummary[];

Needs["TestSuite`", "TestSuite.m"];
Needs["LoopFunctions`", "LoopFunctions.m"];

Print["testing divergences ..."];

assumptions = m > 0 && p > 0 && Q > 0;

TestEquality[
    Print["   testing A0 ..."];
    Coefficient[A0[m,Q] /. FullMomentum[], Delta],
    Coefficient[A0[m,Q] /. ZeroMomentum[], Delta]
];

loopFunctions = {B0, B1, B22, B22tilde, F, G, H};

For[i = 1, i <= Length[loopFunctions], i++,
    Print["   testing ", loopFunctions[[i]], "[p,m,m,Q] ..."];
    expr1 = Limit[loopFunctions[[i]][p,m,m,Q] /. FullMomentum[], p -> 0,
                  Assumptions :> assumptions];
    expr2 = loopFunctions[[i]][0,m,m,Q] /. ZeroMomentum[];
    TestEquality[FullSimplify[expr1 - expr2, assumptions], 0];
   ];

Print["testing B0[] function ..."];

TestEquality[Simplify[
    (B0[1, 2, 3, 4] /. FullMomentum[]) -
    LoopFunctions`Private`B0integral[1, 2, 3, 4]],
    0
];

PrintTestSummary[];

Needs["TestSuite`", "TestSuite.m"];
Needs["LoopFunctions`", "LoopFunctions.m"];

Print["testing p -> 0 limit ..."];

assumptions = m > 0 && m1 > 0 && m2 > 0 && p > 0 && Q > 0;

TestEquality[
    Print["   testing A0 ..."];
    Coefficient[A0[m,Q] /. FullMomentum[], Delta],
    Coefficient[A0[m,Q] /. ZeroMomentum[], Delta]
];

loopFunctions = {B0, B1, B00, B22, B22tilde, F, G, H};

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

Print["testing divergences ..."];

loopFunctions = {
    A0[m,Q],
    B0[p,m1,m2,Q],
    B1[p,m1,m2,Q],
    B00[p,m1,m2,Q],
    B11[p,m1,m2,Q],
    B22[p,m1,m2,Q],
    B22tilde[p,m1,m2,Q],
    F[p,m1,m2,Q],
    G[p,m1,m2,Q],
    H[p,m1,m2,Q]
};

For[i = 1, i <= Length[loopFunctions], i++,
    Print["   testing ", loopFunctions[[i]], " ..."];
    expr1 = 1/2 Q D[loopFunctions[[i]] /. FullMomentum[], Q];
    expr2 = loopFunctions[[i]] /. Divergence[] /. Delta -> 1;
    TestEquality[FullSimplify[expr1 - expr2, assumptions], 0];
   ];

loopFunctions = {
    C0[p1,p2,m1,m2,m3,Q],
    D0[p1,p2,p3,m1,m2,m3,m4,Q],
    D27[p1,p2,p3,m1,m2,m3,m4,Q]
};

For[i = 1, i <= Length[loopFunctions], i++,
    Print["   testing ", loopFunctions[[i]], " ..."];
    expr2 = loopFunctions[[i]] /. Divergence[] /. Delta -> 1;
    TestEquality[expr2, 0];
   ];

PrintTestSummary[];

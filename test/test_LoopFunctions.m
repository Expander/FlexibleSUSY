Needs["TestSuite`", "TestSuite.m"];
Needs["LoopFunctions`", "LoopFunctions.m"];

assumptions = m > 0 && m1 > 0 && m2 > 0 && p > 0 && Q > 0;

RunTests[] := (
    Print["running all tests with $BPMZSign = ", $BPMZSign];
    Print["testing p -> 0 limit ..."];

    TestEquality[
        Print["   testing A0 ..."];
        Coefficient[A0[m,Q] /. LFFull[], Delta],
        Coefficient[A0[m,Q] /. LFZeroMomentum[], Delta]
    ];

    loopFunctions = {B0, B1, B00, B22, B22tilde, F, G, H};

    For[i = 1, i <= Length[loopFunctions], i++,
        Print["   testing ", loopFunctions[[i]], "[p,m,m,Q] ..."];
        expr1 = Limit[loopFunctions[[i]][p,m,m,Q] /. LFFull[], p -> 0,
                      Assumptions :> assumptions];
        expr2 = loopFunctions[[i]][0,m,m,Q] /. LFZeroMomentum[];
        TestEquality[FullSimplify[expr1 - expr2, assumptions], 0];
       ];

    Print["testing m -> 0 limits ..."];

    loopFunctions = {
        {B0[p,m1,m2,Q], m1 -> 0},
        {B0[p,m1,m2,Q], m2 -> 0},
        {B0[0,m1,m2,Q], m1 -> 0},
        {B0[0,m1,m2,Q], m2 -> 0},
        {B1[p,m1,m2,Q], m1 -> 0},
        {B1[p,m1,m2,Q], m2 -> 0},
        {B1[0,m1,m2,Q], m1 -> 0},
        {B1[0,m1,m2,Q], m2 -> 0},
        {B00[p,m1,m2,Q], m1 -> 0},
        {B00[p,m1,m2,Q], m2 -> 0},
        {B00[0,m1,m2,Q], m1 -> 0},
        {B00[0,m1,m2,Q], m2 -> 0},
        {B11[p,m1,m2,Q], m1 -> 0},
        {B11[p,m1,m2,Q], m2 -> 0},
        {B22[p,m1,m2,Q], m1 -> 0},
        {B22[p,m1,m2,Q], m2 -> 0},
        {B22[0,m1,m2,Q], m1 -> 0},
        {B22[0,m1,m2,Q], m2 -> 0},
        {B22tilde[p,m1,m2,Q], m1 -> 0},
        {B22tilde[p,m1,m2,Q], m2 -> 0},
        {B22tilde[0,m1,m2,Q], m1 -> 0},
        {B22tilde[0,m1,m2,Q], m2 -> 0},
        {F[p,m1,m2,Q], m1 -> 0},
        {F[p,m1,m2,Q], m2 -> 0},
        {F[0,m1,m2,Q], m1 -> 0},
        {F[0,m1,m2,Q], m2 -> 0},
        {G[p,m1,m2,Q], m1 -> 0},
        {G[p,m1,m2,Q], m2 -> 0},
        {G[0,m1,m2,Q], m1 -> 0},
        {G[0,m1,m2,Q], m2 -> 0},
        {H[p,m1,m2,Q], m1 -> 0},
        {H[p,m1,m2,Q], m2 -> 0},
        {H[0,m1,m2,Q], m1 -> 0},
        {H[0,m1,m2,Q], m2 -> 0}
    };

    For[i = 1, i <= Length[loopFunctions], i++,
        lf = First[loopFunctions[[i]]];
        limit = Drop[loopFunctions[[i]], 1];
        Print["   testing ", lf, " in the limit ", limit," ..."];
        expr1 = Limit[lf /. LFFull[], Sequence @@ limit, Assumptions :> assumptions];
        expr2 = (lf /. limit) /. LFFull[];
        TestEquality[FullSimplify[expr1 - expr2, assumptions], 0];
       ];

    Print["testing B0[] function ..."];

    TestEquality[Simplify[
        (B0[1, 2, 3, 4] /. LFFull[]) -
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
        expr1 = 1/2 Q D[loopFunctions[[i]] /. LFFull[], Q];
        expr2 = Coefficient[loopFunctions[[i]] /. LFDivergence[], Delta];
        TestEquality[FullSimplify[expr1 - expr2, assumptions], 0];
       ];

    loopFunctions = {
        C0[p1,p2,m1,m2,m3,Q],
        D0[p1,p2,p3,m1,m2,m3,m4,Q],
        D27[p1,p2,p3,m1,m2,m3,m4,Q]
    };

    For[i = 1, i <= Length[loopFunctions], i++,
        Print["   testing ", loopFunctions[[i]], " ..."];
        expr2 = loopFunctions[[i]] /. LFDivergence[];
        TestEquality[expr2, 0];
       ];

    Print["testing logarithms ..."];

    loopFunctions = {
        A0[m, mu],
        B0[p, m1, m2, mu],
        B1[p, m1, m2, mu],
        B00[p, m1, m2, mu],
        B11[p, m1, m2, mu],
        B22[p, m1, m2, mu],
        B22tilde[p, m1, m2, mu],
        C0[p1, p2, m1, m2, m3, mu],
        D0[p1, p2, p3, m1, m2, m3, m4, mu],
        D27[p1, p2, p3, m1, m2, m3, m4, mu],
        F[p, m1, m2, mu],
        G[p, m1, m2, mu],
        H[p, m1, m2, mu]
    };

    divs = loopFunctions /. LFDivergence[] /. Delta -> Log[mu^2];
    logs = loopFunctions /. LFScaleDependence[];

    divs = (mu D[#, mu]) & /@ divs;
    logs = (mu D[#, mu]) & /@ logs;

    TestEquality[Simplify[divs - logs], Table[0,{i,1,Length[loopFunctions]}]];
);

$BPMZSign = 1;
RunTests[];

$BPMZSign = -1;
RunTests[];

PrintTestSummary[];

Needs["TestSuite`", "TestSuite.m"];
Needs["LoopFunctions`", "LoopFunctions.m"];

assumptions = m > 0 && m1 > 0 && m2 > 0 && p > 0 && Q > 0;

RunTests[] := (
    Print["running all tests with $BPMZSign = ", $BPMZSign];
    Print["testing p -> 0 limit ..."];

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
        expr2 = loopFunctions[[i]] /. Divergence[];
        TestEquality[expr2, 0];
       ];

    Print["testing logaritms ..."];

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

    divs = loopFunctions /. Divergence[] /. Delta -> Log[mu^2];
    logs = loopFunctions /. Logarithms[];

    divs = (mu D[#, mu]) & /@ divs;
    logs = (mu D[#, mu]) & /@ logs;

    TestEquality[Simplify[divs - logs], Table[0,{i,1,Length[loopFunctions]}]];
);

$BPMZSign = 1;
RunTests[];

$BPMZSign = -1;
RunTests[];

PrintTestSummary[];

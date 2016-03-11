Needs["TestSuite`", "TestSuite.m"];
Needs["RGIntegrator`", "RGIntegrator.m"];

On[Assert];

betas = {
    {},
    {{g}},
    {{g}, {l}},
    {{g, h a g}},
    {{g, h a g, h^2 b g^2}},
    {{g, h a g, h^2 b g^2 l^2},
     {l, h (c l + d g), h^2 (e l^2 + f g^2 + k l g)}}
};

Print["Testing running up and back down ..."]

For[i = 1, i <= Length[betas], i++,
    For[p = 1, p <= Length[betas[[i]]], p++,
        par = First[betas[[i, p]]];
        Print["  Running ", par, " using ", InputForm[betas[[i]]]];
        runningPar = FullSimplify[
            Normal@Series[(par[Q0] /. RGIntegrate[betas[[i]], Q0, Q1] /. 
                           RGIntegrate[betas[[i]], Q1, Q0]),
                          {h, 0, Length[betas[[i]] - 1]}], Assumptions :> Q0 > 0 && Q1 > 0];
        TestEquality[runningPar, par[Q0]];
       ];
   ];

Print["Testing scale dependence of solution ..."]

For[i = 1, i <= Length[betas], i++,
    pars = #[[1]]& /@ betas[[i]];
    parsQ0 = #[[1]][Q0]& /@ betas[[i]];
    sols = parsQ0 /. RGIntegrate[betas[[i]], Q0, Q];
    deriv = (Q0 D[#, Q0]& /@ sols) /. Q -> Q0;
    beta = Total /@ (Drop[#,1]& /@ (betas[[i]] /. (Rule[#,#[Q0]]& /@ pars)));

    Print["  Testing ", InputForm[beta]];

    TestEquality[Simplify[deriv - beta], Table[0,{i,1,Length[betas[[i]]]}]];
   ];

PrintTestSummary[];

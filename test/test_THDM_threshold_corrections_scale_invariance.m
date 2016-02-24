Needs["TestSuite`", "TestSuite.m"];

FlexibleSUSY`$flexiblesusyMetaDir = FileNameJoin[{Directory[], "meta"}];

mssmRGEs = FileNameJoin[{Directory[], "Output", "MSSM", "RGEs", "BetaGauge.m"}];
thdmRGEs = FileNameJoin[{Directory[], "Output", "THDM-II", "RGEs", "BetaLijkl.m"}];
thdmThresholds = FileNameJoin[{Directory[], "model_files", "THDMIIMSSMBC", "full_1L_thresholds.m"}];

gRules = {g1 -> Sqrt[5/3] gY};

approx = {
    Yu[i_, j_] :> 0 /; i < 3 || j < 3 || i != j,
    Yd[i_, j_] :> 0 /; i < 3 || j < 3 || i != j,
    Ye[i_, j_] :> 0 /; i < 3 || j < 3 || i != j,
    Tu[i_, j_] :> 0 /; i < 3 || j < 3 || i != j,
    Td[i_, j_] :> 0 /; i < 3 || j < 3 || i != j,
    Te[i_, j_] :> 0 /; i < 3 || j < 3 || i != j,
    Nc -> 3,
    Yu[3, 3] -> ht,
    Yd[3, 3] -> hb,
    Ye[3, 3] -> htau,
    Tu[3, 3] -> At ht,
    Td[3, 3] -> Ab hb,
    Te[3, 3] -> Atau htau,
    mse[_] :> MSUSY,
    msu[_] :> MSUSY,
    msd[_] :> MSUSY,
    msq[_] :> MSUSY,
    msl[_] :> MSUSY,
    M1 -> Mu,
    M2 -> Mu,
    Abs[p_] :> p,
    Conjugate[p_] :> p,
    Re[p_] :> p,
    THRESHOLD -> 1,
    flagSferm -> 1,
    flagIno -> 1,
    trace[Yu, Adj[Yu]] -> Yu[3, 3] Conjugate[Yu[3, 3]], 
    trace[Yd, Adj[Yd]] -> Yd[3, 3] Conjugate[Yd[3, 3]],
    trace[Ye, Adj[Ye]] -> Ye[3, 3] Conjugate[Ye[3, 3]],
    trace[Yu, Adj[Yu], Yu, Adj[Yu]] -> Yu[3, 3]^2 Conjugate[Yu[3, 3]]^2,
    trace[Yd, Adj[Yd], Yd, Adj[Yd]] -> Yd[3, 3]^2 Conjugate[Yd[3, 3]]^2,
    trace[Ye, Adj[Ye], Ye, Adj[Ye]] -> Ye[3, 3]^2 Conjugate[Ye[3, 3]]^2,
    trace[Yd, Adj[Yu], Yu, Adj[Yd]] -> Abs[Yd[3, 3]]^2 Abs[Yu[3, 3]]^2
};

betag1MSSM = Cases[Get[mssmRGEs], {g1, b_, __} :> b][[1]];
betagYMSSM = (Sqrt[3/5] betag1MSSM /. gRules);
betag2MSSM = Cases[Get[mssmRGEs], {g2, b_, __} :> b][[1]] /. gRules;

(* tree-level in SARAH convention *)
lambdaTree = {
    1/8 (gY^2 + g2^2),
    1/8 (gY^2 + g2^2),
    1/4 (-gY^2 + g2^2),
    -1/2 g2^2,
    0, 0, 0
};

lambdaTreeRules = {
    Lambda1 -> lambdaTree[[1]],
    Lambda2 -> lambdaTree[[2]],
    Lambda3 -> lambdaTree[[3]],
    Lambda4 -> lambdaTree[[4]],
    Lambda5 -> lambdaTree[[5]],
    Lambda6 -> lambdaTree[[6]],
    Lambda7 -> lambdaTree[[7]]
};

(* beta functions of lambda_i in the MSSM *)

betaLambdaMSSM = (Dt[#] & /@ lambdaTree) /. {
    Dt[gY] -> betagYMSSM,
    Dt[g2] -> betag2MSSM
} // Simplify;

(* beta functions of lambda_i in the THDM *)
betaLambdaTHDM = {
    Cases[Get[thdmRGEs], {Lambda1, b_, __} :> b][[1]],
    Cases[Get[thdmRGEs], {Lambda2, b_, __} :> b][[1]],
    Cases[Get[thdmRGEs], {Lambda3, b_, __} :> b][[1]],
    Cases[Get[thdmRGEs], {Lambda4, b_, __} :> b][[1]],
    Cases[Get[thdmRGEs], {Lambda5, b_, __} :> b][[1]],
    Cases[Get[thdmRGEs], {Lambda6, b_, __} :> b][[1]],
    Cases[Get[thdmRGEs], {Lambda7, b_, __} :> b][[1]]
} /. gRules;

(* difference of the beta functions in the two models *)

betaDiff =
  Expand[(betaLambdaTHDM - betaLambdaMSSM) /. lambdaTreeRules //. approx];

(* keep only Yukawa terms *)
betaDiff = Simplify[
    If[AtomQ[#], #,
       Plus @@ Select[MonomialList[#], 
                      Function[x, !FreeQ[x, ht] || !FreeQ[x, hb] || !FreeQ[x, htau]]]]]& /@ betaDiff;

(* load threshold corrections *)
Get[thdmThresholds];

(* convert to SARAH convention *)
lamSARAH = lamWagnerLee;
lamSARAH[[1]] = lamWagnerLee[[1]]/2;
lamSARAH[[2]] = lamWagnerLee[[2]]/2;

thresh = lamSARAH /. coefficients /. Summation -> Sum //. approx //. loopFunctions;

(* mu-dependence of threshold corrections *)
threshMuDep = Expand[Q D[thresh, Q] 16 Pi^2];

(* keep only Yukawa terms *)
threshMuDep = Simplify[
    If[AtomQ[#], #, 
       Plus @@ Select[MonomialList[#], 
                      Function[x, !FreeQ[x, ht] || !FreeQ[x, hb] || !FreeQ[x, htau]]]]]& /@ threshMuDep;

TestEquality[Simplify[betaDiff - threshMuDep], Table[0, {i,1,7}]];

PrintTestSummary[];

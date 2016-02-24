Needs["TestSuite`", "TestSuite.m"];

FlexibleSUSY`$flexiblesusyMetaDir = FileNameJoin[{Directory[], "meta"}];

mssmRGEs = FileNameJoin[{Directory[], "Output", "MSSM", "RGEs", "BetaGauge.m"}];
thdmRGEs = FileNameJoin[{Directory[], "Output", "THDM-II", "RGEs", "BetaLijkl.m"}];
thdmThresholds = FileNameJoin[{Directory[], "meta", "THDM", "Thresholds_1L_full.m"}];

gRules = {g1 -> Sqrt[5/3] gY};

(* Note: Yukawa couplings are transposed, compared to SARAH *)
YuMat = Table[Yu[j, i], {i, 1, 3}, {j, 1, 3}];
YdMat = Table[Yd[j, i], {i, 1, 3}, {j, 1, 3}];
YeMat = Table[Ye[j, i], {i, 1, 3}, {j, 1, 3}];

approx = {
    Nc -> 3,
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
    trace[Yu, Adj[Yu]] -> Tr[YuMat.ConjugateTranspose[YuMat]],
    trace[Yd, Adj[Yd]] -> Tr[YdMat.ConjugateTranspose[YdMat]],
    trace[Ye, Adj[Ye]] -> Tr[YeMat.ConjugateTranspose[YeMat]],
    trace[Yu, Adj[Yu], Yu, Adj[Yu]] -> Tr[YuMat.ConjugateTranspose[YuMat].YuMat.ConjugateTranspose[YuMat]],
    trace[Yd, Adj[Yd], Yd, Adj[Yd]] -> Tr[YdMat.ConjugateTranspose[YdMat].YdMat.ConjugateTranspose[YdMat]],
    trace[Ye, Adj[Ye], Ye, Adj[Ye]] -> Tr[YeMat.ConjugateTranspose[YeMat].YeMat.ConjugateTranspose[YeMat]],
    trace[Yd, Adj[Yu], Yu, Adj[Yd]] -> Tr[YdMat.ConjugateTranspose[YuMat].YuMat.ConjugateTranspose[YdMat]]
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

(* load threshold corrections *)
Get[thdmThresholds];

(* convert to SARAH convention *)
lamSARAH = lamWagnerLee;
lamSARAH[[1]] = lamWagnerLee[[1]]/2;
lamSARAH[[2]] = lamWagnerLee[[2]]/2;

thresh = lamSARAH /. flags /. coefficients /. Summation -> Sum //. approx //. loopFunctions;

(* mu-dependence of threshold corrections *)
threshMuDep = Expand[Q D[thresh, Q] 16 Pi^2];

TestEquality[Simplify[betaDiff - threshMuDep], Table[0, {i,1,7}]];

PrintTestSummary[];

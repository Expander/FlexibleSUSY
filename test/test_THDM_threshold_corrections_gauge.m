(* This test compares the implementation of arxiv:0901.2065 with the
   implementation of arxiv:hep-ph/9307201 and arxiv:1508.00576 *)

Needs["TestSuite`", "TestSuite.m"];

FlexibleSUSY`$flexiblesusyMetaDir = FileNameJoin[{Directory[], "meta"}];

Get[FileNameJoin[{Directory[], "model_files", "THDMIIMSSMBC", "full_1L_thresholds.m"}]];

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
    THRESHOLD -> 1
};

result = lamWagnerLee /. coefficients /. Summation -> Sum;

TestEquality[FreeQ[result, Undef] && FreeQ[result, Null], True];

result = result /. flags //. approx //. loopFunctions /. Q -> MSUSY // Expand;

(* collect only terms involving Yukawa couplings *)
result = Simplify[Plus @@ Select[List @@ Expand[#],
                                 Function[x,!FreeQ[x,ht] || !FreeQ[x,hb] || !FreeQ[x,htau]]]]& /@ result

Print["testing THDM threshold corrections[] ..."];

Get["model_files/THDMIIMSSMBC/FlexibleSUSY.m.in"];

renameRules = {
    AtauInput -> Atau,
    AtInput -> At,
    AbInput -> Ab,
    Yu[3,3] -> ht,
    Yd[3,3] -> hb,
    Ye[3,3] -> htau,
    MuInput -> Mu,
    \[Mu] -> Mu
};

gRules = {
    gbar -> (g2^2 + gY^2)/4,
    gbarm -> (-g2^2 + gY^2)/4,
    g1 -> gY Sqrt[5/3],
    GUTNormalization[g1] -> Sqrt[3/5]
};

deltaLambdaTh = {
    deltaLambda1th1L, deltaLambda2th1L, deltaLambda3th1L, 
    deltaLambda4th1L, deltaLambda5th1L, deltaLambda6th1L, 
    deltaLambda7th1L
} //. renameRules //. gRules //. approx;

deltaLambdaPhi = {
    deltaLambda1Phi1L, deltaLambda2Phi1L, deltaLambda3Phi1L, 
    deltaLambda4Phi1L, deltaLambda5Phi1L, deltaLambda6Phi1L, 
    deltaLambda7Phi1L
} //. renameRules //. gRules //. approx;

deltaLambdaFull = deltaLambdaTh + deltaLambdaPhi;

lamThDiff  = Simplify[result - deltaLambdaFull];

TestEquality[lamThDiff, Table[0, {i,1,7}]];

PrintTestSummary[];

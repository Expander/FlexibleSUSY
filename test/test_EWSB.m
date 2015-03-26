Needs["TestSuite`", "TestSuite.m"];
Needs["EWSB`", "EWSB.m"];

Print["testing MSSM-like EWSB for Mu and BMu ..."];

FlexibleSUSY`FSSolveEWSBTimeConstraint = 120;

mssmEwsbEqs = {
    mu^2 + x^2 + x y + z + 5,
    Bmu  - x^2 + x y + z + 5
};

mssmEwsbOutputParameters = { mu, Bmu };

Parameters`AddRealParameter[mssmEwsbOutputParameters];

{mssmSolution, mssmFreePhases} = EWSB`FindSolutionAndFreePhases[mssmEwsbEqs, mssmEwsbOutputParameters];

TestEquality[mssmFreePhases, {FlexibleSUSY`Sign[mu]}];

mssmFullSolution = EWSB`Private`FindSolution[mssmEwsbEqs, mssmEwsbOutputParameters];

TestEquality[Sort /@ mssmFullSolution,
             Sort /@ { {{mu -> -Sqrt[-5 - x^2 - x*y - z]},
                        {mu -> Sqrt[-5 - x^2 - x*y - z]}},
                       {{B[mu] -> -5 + x^2 - x*y - z}}
                     }];

TestEquality[Sort[mssmSolution],
             Sort[{ B[mu] -> -5 + x^2 - x*y - z,
                    mu -> FlexibleSUSY`Sign[mu] Sqrt[-5 - x^2 - x*y - z]
                  }]];

Print["testing NMSSM-like EWSB for Kappa, vS and mS2 ..."];

muEff = lambda s;

BmuEff = Alambda + kappa s;

nmssmEwsbEqs = {
    vu (mHu2 + muEff^2 + lambda^2 vd^2 + g^2 (vu^2 - vd^2)) - vd muEff BmuEff,
    vd (mHd2 + muEff^2 + lambda^2 vu^2 + g^2 (vd^2 - vu^2)) - vu muEff BmuEff,
    s (mS2 + X + kappa s Akappa + kappa^2 s^2) + Y
};

nmssmEwsbOutputParameters = { s, kappa, mS2 };

Parameters`AddRealParameter[nmssmEwsbOutputParameters];

{nmssmSolution, nmssmFreePhases} = EWSB`FindSolutionAndFreePhases[nmssmEwsbEqs, nmssmEwsbOutputParameters];

TestEquality[nmssmFreePhases, {FlexibleSUSY`Sign[s]}];

nmssmFullSolution = EWSB`Private`FindSolution[nmssmEwsbEqs, nmssmEwsbOutputParameters];

TestEquality[Sort /@ nmssmFullSolution,
             Sort /@ {{{s -> -(Sqrt[-(mHd2*vd^2) - g^2*vd^4 + mHu2*vu^2 + g^2*vu^4]/
                               Sqrt[lambda^2*vd^2 - lambda^2*vu^2])},
                       {s -> Sqrt[-(mHd2*vd^2) - g^2*vd^4 + mHu2*vu^2 + g^2*vu^4]/
                        Sqrt[lambda^2*vd^2 - lambda^2*vu^2]}},
                      {{kappa ->
                        (-(Alambda*lambda*s*vd) + mHu2*vu + lambda^2*s^2*vu - g^2*vd^2*vu +
                         lambda^2*vd^2*vu + g^2*vu^3)/(lambda*s^2*vd)}},
                      {{mS2 -> (-(Akappa*kappa*s^2) - kappa^2*s^3 - s*X - Y)/s}}
                     }];

TestEquality[Sort[nmssmSolution],
             Sort[{ kappa -> (-(Alambda*lambda*s*vd) + mHu2*vu + lambda^2*s^2*vu -
                              g^2*vd^2*vu + lambda^2*vd^2*vu + g^2*vu^3)/(lambda*s^2*vd),
                    mS2 -> (-(Akappa*kappa*s^2) - kappa^2*s^3 - s*X - Y)/s,
                    s -> (Sqrt[-(mHd2*vd^2) - g^2*vd^4 + mHu2*vu^2 + g^2*vu^4]*
                          FlexibleSUSY`Sign[s])/Sqrt[lambda^2*vd^2 - lambda^2*vu^2]
                  }]];

Print["testing UMSSM-like EWSB for mHu2, mHd2 and mS2 ..."];

umssmEwsbEqs = {
    vu mHu2 + g^2 (vu^2 - vd^2),
    vd mHd2 + g^2 (vd^2 - vu^2),
    s mS2 + X
};

umssmEwsbOutputParameters = { mHu2, mHd2, mS2 };

Parameters`AddRealParameter[umssmEwsbOutputParameters];

{umssmSolution, umssmFreePhases} = EWSB`FindSolutionAndFreePhases[umssmEwsbEqs, umssmEwsbOutputParameters];

TestEquality[umssmFreePhases, {}];

umssmFullSolution = EWSB`Private`FindSolution[umssmEwsbEqs, umssmEwsbOutputParameters];

TestEquality[Sort /@ umssmFullSolution,
             Sort /@ { {{mHu2 -> (g^2*(vd^2 - vu^2))/vu}},
                       {{mHd2 -> -(g^2*(vd^2 - vu^2))/vd}},
                       {{mS2 -> -(X/s)}}
                     }];

TestEquality[Sort /@ umssmSolution,
             Sort /@ { mHu2 -> (g^2*(vd^2 - vu^2))/vu,
                       mHd2 -> -(g^2*(vd^2 - vu^2))/vd,
                       mS2 -> -(X/s)
                     }];

Print["testing EWSB for vu, vd, s ..."];

(* This test ensures that the algorith works even in cases where only
   one independent sub-equation can be found. *)

oneIndependentSubeq = {
    vu mHu2 + X vd + s,
    vd mHd2 + Y vu + s,
    s mS2 + Z           (* this Eq. is independent of all the others*)
};

oneIndependentSubeqEwsbOutputParameters = { vu, vd, s };

Parameters`AddRealParameter[oneIndependentSubeqEwsbOutputParameters];

{oneIndependentSubeqSolution, oneIndependentSubeqFreePhases} =
    EWSB`FindSolutionAndFreePhases[oneIndependentSubeq, oneIndependentSubeqEwsbOutputParameters];

TestEquality[oneIndependentSubeqFreePhases, {}];

TestEquality[Sort /@ oneIndependentSubeqSolution,
             Sort /@ { s -> -(Z/mS2),
                       vu -> -((s*(mHd2 - X))/(mHd2*mHu2 - X*Y)),
                       vd -> (-s - mHu2*vu)/X
                     }];

Print["testing NSM EWSB for mH2, mS2 ..."];

nsmEwsbOutputParameters = {mH2, mS2};

Parameters`AddRealParameter[nsmEwsbOutputParameters];

nsmEwsbEqs = {
    mH2*vH - vH^3*l1 - vH*vS^2*l3 - vH*vS*l4 - tadpole[1],
    2*mS2*vS - 4*vS^3*l2 - vH^2*vS*l3 - (vH^2*l4)/2 - 3*vS^2*l5 - tadpole[2]
};

nsmFullSolution = EWSB`Private`FindSolution[nsmEwsbEqs, nsmEwsbOutputParameters];

TestEquality[Sort /@ nsmFullSolution,
             Sort /@ { {{mH2 -> (vH^3*l1 + vH*vS^2*l3 + vH*vS*l4 + tadpole[1])/vH}},
                       {{mS2 -> (8*vS^3*l2 + 2*vH^2*vS*l3 + vH^2*l4 + 6*vS^2*l5 + 2*tadpole[2])/(4*vS)}}
                     }];

Print["testing cMSSM-like EWSB for |Mu| and BMu ..."];

cmssmEwsbEqs = {
    Susyno`LieGroups`conj[Mu] Mu + x^2 + x y + z + 5,
    BMu - x^2 + x y + z + 5
};

cmssmEwsbOutputParameters = { Mu, BMu };

TestEquality[Parameters`IsRealParameter[Mu], False];
TestEquality[Parameters`IsRealParameter[BMu], False];

Print["\t calling FindSolution[] ..."];

cmssmFullSolution = EWSB`Private`FindSolution[cmssmEwsbEqs, cmssmEwsbOutputParameters];

TestEquality[Sort /@ cmssmFullSolution,
             Sort /@ { {{Mu -> -Sqrt[-5 - x^2 - x*y - z]},
                        {Mu -> Sqrt[-5 - x^2 - x*y - z]}},
                       {{B[Mu] -> -5 + x^2 - x*y - z}}
                     }];

Print["\t calling FindSolutionAndFreePhases[] ..."];

{cmssmSolution, cmssmFreePhases} = EWSB`FindSolutionAndFreePhases[cmssmEwsbEqs, cmssmEwsbOutputParameters];

TestEquality[cmssmFreePhases, {FlexibleSUSY`Phase[Mu]}];

TestEquality[Sort[cmssmSolution],
             Sort[{ B[Mu] -> -5 + x^2 - x*y - z,
                    Mu -> FlexibleSUSY`Phase[Mu] Sqrt[-5 - x^2 - x*y - z]
                  }]];

PrintTestSummary[];

Needs["TestSuite`", "TestSuite.m"];
Needs["Constraint`", "Constraint.m"];

Constraint`SetInputParameters[{g0, m0}];
Constraint`SetModelParameters[{g1, g2, Yu, Yd}];
Constraint`SetOutputParameters[{VZ, Se}];
Constraint`SetBetaFunctions[
    { BetaFunction`BetaFunction[g1, CConversion`ScalarType["double"], 2 g1],
      BetaFunction`BetaFunction[g2, CConversion`ScalarType["double"], 3 g2],
      BetaFunction`BetaFunction[Yu[i1,i2], CConversion`MatrixType["DoubleMatrix", 3, 3], 4 Yu],
      BetaFunction`BetaFunction[Yd[i1,i2], CConversion`MatrixType["DoubleMatrix", 3, 3], 5 Yd]
    }];

Print["testing GetListOfIndexedParameters[] ..."];

TestEquality[Private`GetListOfIndexedParameters[],
             {g1, g2,
              Yu[1,1], Yu[1,2], Yu[1,3],
              Yu[2,1], Yu[2,2], Yu[2,3],
              Yu[3,1], Yu[3,2], Yu[3,3],
              Yd[1,1], Yd[1,2], Yd[1,3],
              Yd[2,1], Yd[2,2], Yd[2,3],
              Yd[3,1], Yd[3,2], Yd[3,3]
             }];

Print["testing CalculateScaleFromExpr[] ..."];

TestEquality[Private`CalculateScaleFromExprSymb[1 == 1],
             currentScale
            ];

TestEquality[Private`CalculateScaleFromExprSymb[1 == 2],
             Null
            ];

TestEquality[Private`CalculateScaleFromExprSymb[g1 == g2],
             currentScale*E^((-g1 + g2)/(BETA[g1] - BETA[g2]))
            ];

TestEquality[Private`CalculateScaleFromExprSymb[g1 == 2 g2],
             Simplify[currentScale*E^((-g1 + 2*g2)/(BETA[g1] - 2*BETA[g2]))]
            ];

TestEquality[Private`CalculateScaleFromExprSymb[g1 == g2^2],
             Simplify[currentScale*E^((-g1 + g2^2)/(BETA[g1] - 2*g2*BETA[g2]))]
            ];

TestEquality[Private`CalculateScaleFromExprSymb[g1 == Yu[3,3]],
             Simplify[currentScale*E^((-g1 + Yu[3,3])/(BETA[g1] - BETA[Yu[3,3]]))]
            ];

TestEquality[Private`CalculateScaleFromExprSymb[g1 == g2 Yu[3,3]],
             Simplify[currentScale*E^((-g1 + g2 Yu[3,3])/(BETA[g1] - g2 BETA[Yu[3,3]] - Yu[3,3] BETA[g2]))]
            ];

TestEquality[Private`CalculateScaleFromExprSymb[Yd[3,3] == Yu[3,3]],
             Simplify[currentScale*E^((-Yd[3,3] + Yu[3,3])/(BETA[Yd[3,3]] - BETA[Yu[3,3]]))]
            ];

TestEquality[Private`CalculateScaleFromExprSymb[Sqrt[g1] == g2],
             FullSimplify[currentScale/E^((-Sqrt[g1] + g2)/(-BETA[g1]/(2 Sqrt[g1]) + BETA[g2]))]
            ];

TestEquality[Private`CalculateScaleFromExprSymb[g1 == 2],
             FullSimplify[currentScale*E^((-g1 + 2)/(BETA[g1]))]
            ];

TestEquality[Private`CalculateScaleFromExprSymb[g1 == MZ],
             FullSimplify[currentScale*E^((-g1 + MZ)/(BETA[g1]))]
            ];

(* Print["testing for Mathematica 9 bug where the GUT condition cannot be solved ..."]; *)
(* TestEquality[Solve[(BETA[g1]-BETA[g2]) Log[x/a] == g1-g2, x, Reals], *)
(*              {{x -> a*E^((g1 - g2)/(BETA[g1] - BETA[g2]))}} *)
(*             ]; *)

PrintTestSummary[];

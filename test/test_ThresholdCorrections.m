Needs["TestSuite`", "TestSuite.m"];
Needs["ThresholdCorrections`", "ThresholdCorrections.m"];

Print["testing InvertRelation[] ..."];

TestEquality[Private`InvertRelation[A,B,A],
             {A,B}
            ];

TestEquality[Private`InvertRelation[A,B + C,A],
             {A,B + C}
            ];

TestEquality[Private`InvertRelation[Transpose[A],B,A],
             {A,SARAH`Tp[B]}
            ];

TestEquality[Private`InvertRelation[ConjugateTranspose[A],B,A],
             {A,SARAH`Adj[B]}
            ];

TestEquality[Private`InvertRelation[FlexibleSUSY`Diag[A],B,A],
             {A,FlexibleSUSY`Diag[B]}
            ];


PrintTestSummary[];

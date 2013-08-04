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

TestEquality[Private`InvertRelation[SARAH`Tp[A],B,A],
             {A,SARAH`Tp[B]}
            ];

TestEquality[Private`InvertRelation[ConjugateTranspose[A],B,A],
             {A,SARAH`Adj[B]}
            ];

TestEquality[Private`InvertRelation[SARAH`Adj[A],B,A],
             {A,SARAH`Adj[B]}
            ];

TestEquality[Private`InvertRelation[FlexibleSUSY`Diag[A],B,A],
             {A,FlexibleSUSY`Diag[B]}
            ];

TestEquality[Private`InvertRelation[SARAH`MatMul[A,V],C,A],
             {A,SARAH`MatMul[C,SARAH`Adj[V]]}
            ];

TestEquality[Private`InvertRelation[SARAH`MatMul[V,A],C,A],
             {A,SARAH`MatMul[SARAH`Adj[V],C]}
            ];

TestEquality[Private`InvertRelation[SARAH`MatMul[U,A,V],C,A],
             {A,SARAH`MatMul[SARAH`Adj[U],C,SARAH`Adj[V]]}
            ];

TestEquality[Private`InvertRelation[SARAH`MatMul[SARAH`Adj[U],A,V],C,A],
             {A,SARAH`MatMul[U,C,SARAH`Adj[V]]}
            ];

TestEquality[Private`InvertRelation[SARAH`MatMul[U,A,SARAH`Adj[V]],C,A],
             {A,SARAH`MatMul[SARAH`Adj[U],C,V]}
            ];

TestEquality[Private`InvertRelation[SARAH`MatMul[SARAH`Adj[U],A,SARAH`Adj[V]],C,A],
             {A,SARAH`MatMul[U,C,V]}
            ];


Print["testing ToMatrixExpression[] ..."];

TestEquality[Private`ToMatrixExpression[{}],
             Null
            ];

TestEquality[Private`ToMatrixExpression[{{Y[1,1]}}],
             Y
            ];

TestEquality[Private`ToMatrixExpression[{{Y[1,1],Y[1,2]},{Y[2,1],Y[2,2]}}],
             Y
            ];

TestEquality[Private`ToMatrixExpression[{{Y[1,1],Y[2,1]},{Y[1,2],Y[2,2]}}],
             SARAH`Tp[Y]
            ];

TestEquality[Private`ToMatrixExpression[{{Y[1,1],0},{0,Y[2,2]}}],
             FlexibleSUSY`Diag[Y]
            ];


PrintTestSummary[];

Needs["TestSuite`", "TestSuite.m"];
Needs["BetaFunction`", "BetaFunction.m"];

Print["testing ProtectTensorProducts[] ..."];

TestEquality[BetaFunction`Private`ProtectTensorProducts[a,i,j],
             a];

TestEquality[BetaFunction`Private`ProtectTensorProducts[f[i],i,j],
             f[i]];

TestEquality[BetaFunction`Private`ProtectTensorProducts[f[i] f[1],i,j],
             f[i] f[1]];

TestEquality[BetaFunction`Private`ProtectTensorProducts[f[i] f[j],i,j],
             CConversion`TensorProd[f,f][i,j]];

TestEquality[BetaFunction`Private`ProtectTensorProducts[f[i] f[j],j,i],
             CConversion`TensorProd[f,f][j,i]];

TestEquality[BetaFunction`Private`ProtectTensorProducts[f[i] F[a+b][j],i,j],
             CConversion`TensorProd[f,F[a+b]][i,j]];

TestEquality[BetaFunction`Private`ProtectTensorProducts[f[y][i] F[a+b][j],i,j],
             CConversion`TensorProd[f[y],F[a+b]][i,j]];

TestEquality[BetaFunction`Private`ProtectTensorProducts[(f+g)[i] F[a+b][j],i,j],
             CConversion`TensorProd[f+g,F[a+b]][i,j]];

PrintTestSummary[];

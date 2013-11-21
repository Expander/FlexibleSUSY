Needs["TestSuite`", "TestSuite.m"];
Needs["BetaFunction`", "BetaFunction.m"];

Print["testing ProtectTensorProducts[] ..."];

TestEquality[CConversion`ProtectTensorProducts[a,i,j],
             a];

TestEquality[CConversion`ProtectTensorProducts[f[i],i,j],
             f[i]];

TestEquality[CConversion`ProtectTensorProducts[f[i] f[1],i,j],
             f[i] f[1]];

TestEquality[CConversion`ProtectTensorProducts[f[i] f[j],i,j],
             CConversion`TensorProd[f,f][i,j]];

TestEquality[CConversion`ProtectTensorProducts[f[i] f[j],j,i],
             CConversion`TensorProd[f,f][j,i]];

TestEquality[CConversion`ProtectTensorProducts[f[i] F[a+b][j],i,j],
             CConversion`TensorProd[f,F[a+b]][i,j]];

TestEquality[CConversion`ProtectTensorProducts[f[y][i] F[a+b][j],i,j],
             CConversion`TensorProd[f[y],F[a+b]][i,j]];

TestEquality[CConversion`ProtectTensorProducts[(f+g)[i] F[a+b][j],i,j],
             CConversion`TensorProd[f+g,F[a+b]][i,j]];

TestEquality[CConversion`ProtectTensorProducts[x[i](f[j]+g[j]),i,j],
             CConversion`TensorProd[x,f+g][i,j]];

TestEquality[CConversion`ProtectTensorProducts[(x[i]+y[i])f[j],i,j],
             CConversion`TensorProd[x+y,f][i,j]];

TestEquality[CConversion`ProtectTensorProducts[conj[x[i]] f[j],i,j],
             CConversion`TensorProd[conj[x],f][i,j]];

TestEquality[CConversion`ProtectTensorProducts[x[i] conj[f[j]],i,j],
             CConversion`TensorProd[x,conj[f]][i,j]];

PrintTestSummary[];

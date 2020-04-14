(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

Needs["TestSuite`", "TestSuite.m"];
Needs["MatMulSimplify`", "MatMulSimplify.m"];

(* MatMul symbol *)
MM = SARAH`MatMul;

MMS[expr_] := MatMulSimplify`MatMulSimplify[expr, MM];
MMR[expr_] := MatMulSimplify`MatMulRefine[expr, MM];

Print["testing number of replacement rules ..."];

{simpEx, simpRep} = MMS[1];
TestEquality[Length[simpRep], 0];

{simpEx, simpRep} = MMS[a];
TestEquality[Length[simpRep], 0];

{simpEx, simpRep} = MMS[MM[a]];
TestEquality[Length[simpRep], 0];

{simpEx, simpRep} = MMS[MM[a,b]];
TestEquality[Length[simpRep], 1];

{simpEx, simpRep} = MMS[MM[a,b,c]];
TestEquality[Length[simpRep], 1];

{simpEx, simpRep} = MMS[MM[a,b,c,d]];
TestEquality[Length[simpRep], 1];

{simpEx, simpRep} = MMS[MM[a,b] + MM[c,b]];
TestEquality[Length[simpRep], 2];

{simpEx, simpRep} = MMS[MM[a,b] + MM[a,b,c]];
TestEquality[Length[simpRep], 2];

{simpEx, simpRep} = MMS[MM[a,b] + MM[a,b,c] + MM[a,b,c,d]];
TestEquality[Length[simpRep], 3];

Print["testing equality of full and simplified expressions ..."];

CheckExpr[expr_] :=
    Module[{simpEx, simpRep},
           {simpEx, simpRep} = MMS[expr];
           TestEquality[MatMulRefine[Expand[expr - (simpEx //. simpRep)], MM], 0];
    ];

CheckExpr[1];
CheckExpr[a];
CheckExpr[MM[a]];
CheckExpr[MM[a,b]];
CheckExpr[MM[a,b,c]];
CheckExpr[MM[a,b,c,d]];
CheckExpr[MM[a,b] + MM[c,d]];
CheckExpr[1 + MM[a,b] + MM[a,b,c] + MM[a,b,c,d]];

Print["testing equality of full and simplified MSSM 3-loop beta functions ..."];

betaDir = FileNameJoin[{Directory[], "meta", "MSSM"}];

For[l = 1, l <= 3, l++,
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_BMu.m" }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_g1.m"  }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_g2.m"  }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_g3.m"  }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_M1.m"  }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_M2.m"  }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_M3.m"  }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_md2.m" }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_me2.m" }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_mHd2.m"}]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_mHu2.m"}]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_ml2.m" }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_mq2.m" }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_mu2.m" }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_Mu.m"  }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_TYd.m" }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_TYe.m" }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_TYu.m" }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_Yd.m"  }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_Ye.m"  }]][[l]]];
    CheckExpr[Get[FileNameJoin[{betaDir, "beta_Yu.m"  }]][[l]]];
];

Print["testing MatMulRefine ..."];

TestEquality[MMR[MM[a]], a];
TestEquality[MMR[MM[2]], 2];
TestEquality[MMR[MM[a] + MM[a,2]], 3 a];
TestEquality[MMR[MM[a] + MM[a,2,c]],a + 2 MM[a,c]];

PrintTestSummary[];

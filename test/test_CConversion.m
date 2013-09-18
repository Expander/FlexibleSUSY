Needs["TestSuite`", "TestSuite.m"];
Needs["CConversion`", "CConversion.m"];

Print["testing ConvertGreekLetters[] ..."];

TestEquality[CConversion`Private`ConvertGreekLetters[\[Alpha]], Alpha];
TestEquality[CConversion`Private`ConvertGreekLetters[\[Beta]], Betax]; (* Beta is already defined by Mathematica *)
TestEquality[CConversion`Private`ConvertGreekLetters[\[Mu]]   , Mu];
TestEquality[CConversion`Private`ConvertGreekLetters[\[Zeta]], Zetax]; (* Zeta is already defined by Mathematica *)
SARAH`Delta;
TestEquality[CConversion`Private`ConvertGreekLetters[\[Delta]], Deltax];

Print["testing ToValidCSymbol[] ..."];

TestEquality[ToValidCSymbol[1], 1];
TestEquality[ToValidCSymbol[1.1], 1.1];
TestEquality[ToValidCSymbol[a], a];
TestEquality[ToValidCSymbol[a[1]], a1];
TestEquality[ToValidCSymbol[a[b]], ab];
TestEquality[ToValidCSymbol[a[b,c]], abc];
TestEquality[ToValidCSymbol[a[b,c][d]], abcd];

(* SARAH sometimes appends indices [i1,i2] to express that a symbol is
   of type matrix. These indices should be stripped before
   conversion. *)
TestEquality[ToValidCSymbol[a[Susyno`LieGroups`i1,SARAH`i2]], a];
TestEquality[ToValidCSymbol[a[b,c][Susyno`LieGroups`i1,SARAH`i2]], abc];

(* Test that symbols, which consist of greek letters are correctly
   converted. *)
TestEquality[ToValidCSymbol[\[Lambda]3], Lambda3];
TestEquality[ToValidCSymbol[T[\[Lambda]3]], TLambda3];
TestEquality[ToValidCSymbol[\[Lambda]], Lambdax];
(* no unique name necessary *)
(* TestEquality[ToValidCSymbol[T[\[Lambda]]], TLambda]; *)

Print["testing ToValidCSymbolString[] ..."];

TestEquality[ToValidCSymbolString[1], "1"];
TestEquality[ToValidCSymbolString[1.1], "1.1"];
TestEquality[ToValidCSymbolString[a], "a"];
TestEquality[ToValidCSymbolString[a[b]], "ab"];
TestEquality[ToValidCSymbolString[a[b,c]], "abc"];
TestEquality[ToValidCSymbolString[a[b,c][d]], "abcd"];

(* SARAH sometimes appends indices [i1,i2] to express that a type is
   of type matrix. These indices should be stripped before
   conversion. *)
TestEquality[ToValidCSymbolString[a[Susyno`LieGroups`i1,SARAH`i2]], "a"];
TestEquality[ToValidCSymbolString[a[b][Susyno`LieGroups`i1,SARAH`i2]], "ab"];
TestEquality[ToValidCSymbolString[a[b,c][Susyno`LieGroups`i1,SARAH`i2]], "abc"];

Print["testing RValueToCFormString[] ..."];

TestEquality[RValueToCFormString[1], "1"];
TestEquality[RValueToCFormString[a], "a"];
TestEquality[RValueToCFormString[a+b], "a + b"];
TestEquality[RValueToCFormString[a^1], "a"];
TestEquality[RValueToCFormString[a^2], "Sqr(a)"];
TestEquality[RValueToCFormString[f[x]], "f(x)"];
TestEquality[RValueToCFormString[a MatMul[A]], "a*A"];
TestEquality[RValueToCFormString[a MatMul[A,B]], "a*(A*B)"];
TestEquality[RValueToCFormString[a MatMul[A,B,A]], "a*(A*B*A)"];
TestEquality[RValueToCFormString[a trace[A]], "a*(A).trace()"];
TestEquality[RValueToCFormString[a trace[A,B]], "a*(A*B).trace()"];
TestEquality[RValueToCFormString[a trace[A,B,A]], "a*(A*B*A).trace()"];

TestEquality[RValueToCFormString[MatMul[Adj[A]]], "A.adjoint()"];
TestEquality[RValueToCFormString[MatMul[A Adj[A]]], "A*A.adjoint()"];
TestEquality[RValueToCFormString[trace[Adj[A]]], "(A.adjoint()).trace()"];
TestEquality[RValueToCFormString[trace[A Adj[A]]], "(A*A.adjoint()).trace()"];

(* SARAH sometimes appends indices [i1,i2] to express that a type is
   of type matrix. These indices should be stripped before
   conversion. *)
TestEquality[RValueToCFormString[a[Susyno`LieGroups`i1,SARAH`i2]], "a"];
TestEquality[RValueToCFormString[a[b][Susyno`LieGroups`i1,SARAH`i2]], "a(b)"];
TestEquality[RValueToCFormString[a[b,c][Susyno`LieGroups`i1,SARAH`i2]], "a(b,c)"];

(* test greek symbol conversion *)
TestEquality[RValueToCFormString[\[Mu]], "Mu"];
TestEquality[RValueToCFormString[SARAH`B[\[Mu]]], "BMu"];
TestEquality[RValueToCFormString[A[\[Mu]]], "A(Mu)"];
TestEquality[RValueToCFormString[\[Mu][1]], "Mu(1)"];
TestEquality[RValueToCFormString[\[Mu][1,2]], "Mu(1,2)"];
TestEquality[RValueToCFormString[\[Mu][1,2,3]], "Mu(1,2,3)"];

Print["testing GetHead[] ..."];

TestEquality[GetHead[1], 1];
TestEquality[GetHead[a], a];
TestEquality[GetHead[f[x]], f];
TestEquality[GetHead[f[x][y]], f];
TestEquality[GetHead[f[x[y]]], f];

Print["testing ExpandSums[] ..."];

TestCPPCode[{"", "int result = 0;" <> ExpandSums[1, "result"]}, "1", "int", "1"];
TestCPPCode[{"", "int a = 9;" <>
             "int result = 0;" <>
             ExpandSums[a, "result", "int"]}, "result", "int", "9"];
TestCPPCode[{"", "int a = 9, b = 2;" <>
             "int result = 0;" <>
             ExpandSums[a+b, "result", "int"]}, "result", "int", "11"];
TestCPPCode[{"", "int a = 9, b = 2;" <>
             "int result = 0;" <>
             ExpandSums[a b, "result", "int"]}, "result", "int", "18"];
TestCPPCode[{"", "int a = 2;" <>
             "int result = 0;" <>
             ExpandSums[sum[i,1,3,a], "result", "int"]}, "result", "int", "6"];
TestCPPCode[{"", "int a = 2, b = 3;" <>
             "int result = 0;" <>
             ExpandSums[b sum[i,1,3,a], "result", "int"]}, "result", "int", "18"];
TestCPPCode[{"", "int a = 2, b = 3, c = 4;" <>
             "int result = 0;" <>
             ExpandSums[(a+b) sum[i,1,2,c], "result", "int"]}, "result", "int", "40"];
TestCPPCode[{"", "int a = 2, b = 3;" <>
             "int result = 0;" <>
             ExpandSums[sum[i,1,3,sum[k,1,2,b]], "result", "int"]}, "result", "int", "18"];
TestCPPCode[{"", "int a = 2, b = 3;" <>
             "int result = 0;" <>
             ExpandSums[sum[i,1,3,a sum[k,1,2,b]], "result", "int"]}, "result", "int", "36"];

(* use function as sum argument to avoid automatic simplification *)
TestCPPCode[{"int sqr(int n) { return n*n; }",
             "int result = 0;" <>
             ExpandSums[sum[k,1,4,sqr[k]], "result", "int", " = 0"]}, "result", "int", "30"];
TestCPPCode[{"int sqr(int n) { return n*n; }",
             "int a = 2; int result = 0;" <>
             ExpandSums[a + sum[k,1,4,sqr[k]], "result", "int", " = 0"]}, "result", "int", "32"];
TestCPPCode[{"int sqr(int n) { return n*n; }",
             "int a = 2; int result = 0;" <>
             ExpandSums[a sum[k,1,4,sqr[k]], "result", "int", " = 0"]}, "result", "int", "60"];
TestCPPCode[{"int sqr(int n) { return n*n; }",
             "int a = 2; int result = 0;" <>
             ExpandSums[a sum[k,1,4,sqr[k]], "result", "int", " = 0"]}, "result", "int", "60"];
TestCPPCode[{"int sqr(int n) { return n*n; }",
             "int a = 2, b = 3; int result = 0;" <>
             ExpandSums[(a + b) sum[k,1,4,sqr[k]], "result", "int", " = 0"]}, "result", "int", "150"];
TestCPPCode[{"int sqr(int n) { return n*n; }",
             "int result = 0;" <>
             ExpandSums[sum[i,1,2,sum[k,1,2,sqr[i k]]], "result", "int", " = 0"]}, "result", "int", "25"];
TestCPPCode[{"int sqr(int n) { return n*n; }",
             "int result = 0;" <>
             ExpandSums[sum[i,1,2,(i+1) sum[k,1,4,sqr[k]]], "result", "int", " = 0"]}, "result", "int", "150"];
TestCPPCode[{"int sqr(int n) { return n*n; }",
             "int result = 0;" <>
             ExpandSums[sum[i,1,2,(i+1) sum[k,1,4,sqr[i k]]], "result", "int", " = 0"]}, "result", "int", "420"];

Clear[SARAH`sum];

(* Test ThetaStep expansion *)
TestCPPCode[{"",
             "int result = 0;\n" <>
             ExpandSums[ThetaStep[1,3], "result", "int", " = 0"]}, "result", "int", "1"];
TestCPPCode[{"",
             "int result = 0;\n" <>
             ExpandSums[3 ThetaStep[1,3], "result", "int", " = 0"]}, "result", "int", "3"];
TestCPPCode[{"",
             "int result = 0, i1 = 1;\n" <>
             ExpandSums[2 ThetaStep[i1,3], "result", "int", " = 0"]}, "result", "int", "2"];
TestCPPCode[{"",
             "int result = 0, i1 = 3;\n" <>
             ExpandSums[2 ThetaStep[i1,3], "result", "int", " = 0"]}, "result", "int", "2"];
TestCPPCode[{"",
             "int result = 0, i1 = 4;\n" <>
             ExpandSums[2 ThetaStep[i1,3], "result", "int", " = 0"]}, "result", "int", "0"];

TestCPPCode[{"#include<cassert>\nint Aborts() { assert(false && \"this should not happen\"); return 1; }",
             "int result = 0, i1 = 4;\n" <>
             ExpandSums[Aborts[] ThetaStep[i1,3], "result", "int", " = 0"]}, "result", "int", "0"];

TestCPPCode[{"#include<cassert>\nint Aborts() { assert(false && \"this should not happen\"); return 1; }",
             "int result = 0, i1 = 4;\n" <>
             ExpandSums[Aborts[] (ThetaStep[i1,3] + ThetaStep[i1,2]), "result", "int", " = 0"]}, "result", "int", "0"];

TestCPPCode[{"",
             "int result = 0, i1 = 4;\n" <>
             ExpandSums[2 (ThetaStep[i1,3] + ThetaStep[i1,4]), "result", "int", " = 0"]}, "result", "int", "2"];

TestCPPCode[{"",
             "int result = 0, i1 = 4;\n" <>
             ExpandSums[2 ThetaStep[i1,3] ThetaStep[i1,4], "result", "int", " = 0"]}, "result", "int", "0"];

TestCPPCode[{"",
             "int result = 0, i1 = 3;\n" <>
             ExpandSums[2 ThetaStep[i1,3] ThetaStep[i1,4], "result", "int", " = 0"]}, "result", "int", "2"];

TestCPPCode[{"",
             "int result = 0, i1 = 4;\n" <>
             ExpandSums[2 ThetaStep[i1,5] (ThetaStep[i1,3] + ThetaStep[i1,4]), "result", "int", " = 0"]}, "result", "int", "2"];

TestCPPCode[{"",
             "int result = 0, i1 = 4;\n" <>
             ExpandSums[Power[ThetaStep[i1,5],2], "result", "int", " = 0"]}, "result", "int", "1"];

TestCPPCode[{"",
             "int result = 0, i1 = 4;\n" <>
             ExpandSums[Power[ThetaStep[i1,2],2], "result", "int", " = 0"]}, "result", "int", "0"];

PrintTestSummary[];

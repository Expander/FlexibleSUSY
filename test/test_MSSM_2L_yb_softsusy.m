Needs["TestSuite`", "TestSuite.m"];

Get[FileNameJoin[{"meta", "MSSM", "twoloopbubble.m"}]];

b2l  = Get[FileNameJoin[{"meta", "MSSM", "bquark_2loop_sqcd_decoupling.m"}]];
b2ls = Get[FileNameJoin[{"meta", "MSSM", "dmbas2.m"}]];

colorCA = 3; colorCF = 4/3; GS = g3;
zt2 = Zeta[2]; zt3 = Zeta[3];
MGl = mgl; MT = mt; MB = mb;
mmgl = mgl^2; mmt = mt^2;
mmst1 = mst1^2; mmst2 = mst2^2;
mmsb1 = msb1^2; mmsb2 = msb2^2;
msd1 = msusy; msd2 = msusy;
mmsusy = msusy^2; mmu = scale^2;
Xt = At - MUE CB/SB;
Xb = Ab - MUE SB/CB;
s2t = 2 mt Xt/(mmst1 - mmst2);
s2b = 2 mb Xb/(mmsb1 - mmsb2);
snb = Sin[ArcSin[s2b]/2];

Delta[mass1_, mass2_, mass3_, x_] := (mass1^2 + mass2^2 + mass3^2 - 2*(mass1*mass2 + mass1*mass3 + mass2*mass3))^x;
fpart[0, m1_, m2_, scale_]        := Fin20[m1^2, m2^2, scale^2];
fpart[m1_, m2_, m3_, scale_]      := Fin3[m1^2, m2^2, m3^2, scale^2];
delta3[m1_, m2_, m3_]             := Delta[m1^2, m2^2, m3^2, -1]; 
fin[0, mass1_, mass2_]            := Fin20[mass1, mass2, mmu];
fin[mass1_, mass2_, mass3_]       := Fin3[mass1, mass2, mass3, mmu];

(* check for specific parameter point *)
point = {
   g3 -> 1, mst1 -> 100, mst2 -> 200, msb1 -> 300, msb2 -> 400, 
   msusy -> 500, scale -> 90,
   mt -> 173, mb -> 3, SB -> TB/Sqrt[1 + TB^2], 
   CB -> 1/Sqrt[1 + TB^2], Ab -> 100, At -> 110,
   MUE -> 400, mgl -> 600, TB -> 20
};

TestEquality[PossibleZeroQ[(b2l - b2ls) //. point // N], True];

PrintTestSummary[];

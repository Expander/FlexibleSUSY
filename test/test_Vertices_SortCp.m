Needs["TestSuite`", "TestSuite.m"];
Needs["Vertices`", "Vertices.m"];

workingDirectory = Directory[];
SARAH`SARAH[OutputDirectory] = CreateDirectory[];
Print["Current working directory: ", workingDirectory];
Print["SARAH output directory: ", SARAH`SARAH[OutputDirectory]];

Start["MSSM"];

FlexibleSUSY`FSEigenstates = SARAH`EWSB;

Print["testing SortCp[] ..."];

TestEquality[SortCp @ Cp[conj[USd[{gO2}]], VG, Sd[{gI2}]],
		      Cp[Sd[{gI2}], conj[USd[{gO2}]], VG]];

TestEquality[SortCp @ Cp[conj[VWm], Hpm[{gI1}], hh[{gI2}]],
		      Cp[hh[{gI2}], Hpm[{gI1}], conj[VWm]]];

TestEquality[SortCp @ Cp[conj[VWm], Hpm[{gI1}], Ah[{gI2}]],
		      Cp[Ah[{gI2}], Hpm[{gI1}], conj[VWm]]];

TestEquality[SortCp @ Cp[conj[VWm], hh[{gI2}], Hpm[{gI1}]],
		      Cp[hh[{gI2}], Hpm[{gI1}], conj[VWm]]];

TestEquality[SortCp @ Cp[Ah[{gI2}], conj[VWm], Hpm[{gI1}]],
		      Cp[Ah[{gI2}], Hpm[{gI1}], conj[VWm]]];

TestEquality[SortCp @ Cp[VWm, conj[VWm], conj[VWm], VWm][1],
		      Cp[conj[VWm], conj[VWm], VWm, VWm][2]];

TestEquality[SortCp @ Cp[VWm, conj[VWm], conj[VWm], VWm][2],
		      Cp[conj[VWm], conj[VWm], VWm, VWm][3]];

TestEquality[SortCp @ Cp[VWm, conj[VWm], conj[VWm], VWm][3],
		      Cp[conj[VWm], conj[VWm], VWm, VWm][1]];

TestEquality[SortCp @ Cp[UChi[{gO1}], bar[Cha[{gI1}]], VWm][PL],
		     -Cp[bar[Cha[{gI1}]], UChi[{gO1}], VWm][PR]];

TestEquality[SortCp @ Cp[UChi[{gO1}], bar[Cha[{gI1}]], VWm][PR],
		     -Cp[bar[Cha[{gI1}]], UChi[{gO1}], VWm][PL]];

TestEquality[SortCp @ Cp[UChi[{gO1}], bar[Cha[{gI1}]], Hpm][PL],
		      Cp[bar[Cha[{gI1}]], UChi[{gO1}], Hpm][PL]];

TestEquality[SortCp @ Cp[UChi[{gO1}], bar[Cha[{gI1}]], Hpm][PR],
		      Cp[bar[Cha[{gI1}]], UChi[{gO1}], Hpm][PR]];

Print["testing Vertices`Private`RestoreBarOnMajorana[] ..."];

RBOM := Vertices`Private`RestoreBarOnMajorana;

Block[{SARAH`bar},

TestEquality[RBOM[{UChi[{gO1}], Cha[{gI2}], conj[Hpm[{gI1}]]}, PR],
	     {bar @ UChi[{gO1}], Cha[{gI2}], conj[Hpm[{gI1}]]}];

TestEquality[RBOM[{bar[Cha[{gI1}]], UChi[{gO1}], Hpm}, PR],
	     {bar[Cha[{gI1}]], bar @ UChi[{gO1}], Hpm}];

TestEquality[RBOM[{UChi[{gO2}], Cha[{gI2}], conj[VWm]},
		  LorentzProduct[gamma[lt3], PR]],
	     {bar @ UChi[{gO2}], Cha[{gI2}], conj[VWm]}];

TestEquality[RBOM[{bar[Cha[{gI1}]], UChi[{gO1}], VWm},
		  LorentzProduct[gamma[lt3], PR]],
	     {bar[Cha[{gI1}]], (* w/o bar *) UChi[{gO1}], VWm}];

];

DeleteDirectory[SARAH`SARAH[OutputDirectory], DeleteContents -> True];

PrintTestSummary[];

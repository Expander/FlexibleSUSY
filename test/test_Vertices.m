Needs["TestSuite`", "TestSuite.m"];
Needs["Vertices`", "Vertices.m"];

workingDirectory = Directory[];
$sarahOutputDir = CreateDirectory[];
Print["Current working directory: ", workingDirectory];
Print["SARAH output directory: ", $sarahOutputDir];

Start["MRSSM"];

On[Assert];
Off[EnforceCpColorStructures::cpext];

FlexibleSUSY`FSEigenstates = SARAH`EWSB;
Parameters`SetModelParameters[{
    Yd, Ye, LamTD, LamTU, LamSD, LamSU, Yu, \[Mu], MuD, MuU, g1, g2, g3,
    vd, vu, vT, vS, B[\[Mu]], B[MuD], B[MuU],
    mq2, ml2, mHd2, mHu2, md2, mu2, me2, mS2, mT2, moc2, mRd2, mRu2,
    MDBS, MDWBT, MDGoc}];

FH := Vertices`Private`FieldHead;
SFI := Vertices`Private`StripFieldIndices;

CheckExternalField[extField_, cp:Cp[a_, b_, __]] /;
    FH[a] === FH[b] === extField &&
    SFI[a] === Susyno`LieGroups`conj@SFI[b] := (
    Print[extField, " correctly appears at the first two positions in ", cp];
    True
);

CheckExternalField[extField_, cp_Cp] := (
    Print[extField, " does not appear at the first two positions in ", cp];
    False
);

CheckExternalField[extField_, lst_List] :=
    And @@ (CheckExternalField[extField, #]& /@ Cases[lst, _Cp, Infinity]);

(* external fields are at wrong positions in following Cp[]'s *)

selfEnergySOc = {
    FSSelfEnergy[SOc,
     -C A0[Mass2[SOc]] Cp[conj[SOc], conj[SOc], SOc, SOc] -
      C sum[gI1, 1, 6,
        A0[Mass2[Sd[{gI1}]]] Cp[conj[Sd[{gI1}]], SOc, conj[SOc], Sd[{gI1}]]] -
      C sum[gI1, 1, 6,
        A0[Mass2[Su[{gI1}]]] Cp[conj[Su[{gI1}]], SOc, Su[{gI1}], conj[SOc]]]
    ]};

selfEnergySd = {
    FSSelfEnergy[Sd[gO1, gO2],
     -C A0[Mass2[SOc]]
	Cp[USd[{gO1}], conj[SOc], conj[USd[{gO2}]], SOc] -
      C sum[gI1, 1, 6,
        A0[Mass2[Sd[{gI1}]]]
	Cp[USd[{gO1}], conj[Sd[{gI1}]], Sd[{gI1}], conj[USd[{gO2}]]]] -
      C sum[gI1, 1, 6,
        A0[Mass2[Su[{gI1}]]]
	Cp[Su[{gI1}], conj[USd[{gO2}]], conj[Su[{gI1}]], USd[{gO1}]]]
    ]};

selfEnergySu = {
    FSSelfEnergy[Su[gO1, gO2],
     -C A0[Mass2[SOc]]
	Cp[conj[SOc], USu[{gO1}], conj[USu[{gO2}]], SOc] -
      C sum[gI1, 1, 6,
        A0[Mass2[Sd[{gI1}]]]
	Cp[Sd[{gI1}], conj[USu[{gO2}]], USu[{gO1}], conj[Sd[{gI1}]]]] -
      C sum[gI1, 1, 6,
        A0[Mass2[Su[{gI1}]]]
	Cp[Su[{gI1}], conj[Su[{gI1}]], conj[USu[{gO2}]], USu[{gO1}]]]
    ]};

massMatrices = {
    FSMassMatrix[dummy, SOc, Null],
    FSMassMatrix[dummy, Sd, ZD],
    FSMassMatrix[dummy, Su, ZU]
};

Print["testing EnforceCpColorStructures[] ..."];

selfEnergySOc = EnforceCpColorStructures[selfEnergySOc];
selfEnergySd  = EnforceCpColorStructures[selfEnergySd ];
selfEnergySu  = EnforceCpColorStructures[selfEnergySu ];

TestEquality[CheckExternalField[SOc, selfEnergySOc], True];
TestEquality[CheckExternalField[USd, selfEnergySd ], True];
TestEquality[CheckExternalField[USu, selfEnergySu ], True];

Print["testing color summation ..."];

vertexRulesSOc = VertexRules[selfEnergySOc, massMatrices];
vertexRulesSd  = VertexRules[selfEnergySd , massMatrices];
vertexRulesSu  = VertexRules[selfEnergySu , massMatrices];

TestEquality[Cp[SOc, conj[SOc], conj[SOc], SOc] /. vertexRulesSOc,
    -128*g3^2];

TestEquality[Cp[SOc, conj[SOc], conj[Sd[{gI1}]], Sd[{gI2}]] /. vertexRulesSOc,
    (-24*I)*g3^2*(sum[j1, 1, 3, conj[ZD[gI2, j1]]*ZD[gI1, j1]] -
		  sum[j1, 1, 3, conj[ZD[gI2, 3 + j1]]*ZD[gI1, 3 + j1]])];

TestEquality[Cp[SOc, conj[SOc], conj[Su[{gI1}]], Su[{gI2}]] /. vertexRulesSOc,
    (-24*I)*g3^2*(sum[j1, 1, 3, conj[ZU[gI2, j1]]*ZU[gI1, j1]] -
		  sum[j1, 1, 3, conj[ZU[gI2, 3 + j1]]*ZU[gI1, 3 + j1]])];

TestEquality[Cp[USd[{gO1}], conj[USd[{gO2}]], conj[SOc], SOc] /. vertexRulesSd,
    (-I)*(64*g3^2*sum[j1, 1, 3, Delta[gO1, 3 + j1]*
          Delta[gO2, 3 + j1]] - 64*g3^2*Delta[gO1, gO2]*ThetaStep[gO1, 3])];

TestEquality[Cp[USu[{gO1}], conj[USu[{gO2}]], conj[SOc], SOc] /. vertexRulesSu,
    (-I)*(64*g3^2*sum[j1, 1, 3, Delta[gO1, 3 + j1]*
	  Delta[gO2, 3 + j1]] - 64*g3^2*Delta[gO1, gO2]*ThetaStep[gO1, 3])];

DeleteDirectory[$sarahOutputDir, DeleteContents -> True];

PrintTestSummary[];

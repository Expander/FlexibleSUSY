Needs["TestSuite`", "TestSuite.m"];
Needs["SelfEnergies`", "SelfEnergies.m"];

Print["testing CreateMixingMatrixReplacementRulesFor[] ..."];

Field;
FindMixingMatrixSymbolFor[Field] := ZField;
FindMixingMatrixSymbolFor[SVDField] := {U, V};
GetDimension[Field] := 3;

TestEquality[Private`CreateMixingMatrixReplacementRulesFor[Cp[Field[{i}]]], {}];

TestEquality[Private`CreateMixingMatrixReplacementRulesFor[Cp[UField[{i}]]],
             {ZField[i,Private`idx_] :> SARAH`Delta[i,Private`idx],
              ZField[Private`idx_,i] :> SARAH`Delta[Private`idx,i]}];

TestEquality[Private`CreateMixingMatrixReplacementRulesFor[Cp[USVDField[{i}]]],
             {U[i,Private`idx_] :> SARAH`Delta[i,Private`idx],
              U[Private`idx_,i] :> SARAH`Delta[Private`idx,i],
              V[i,Private`idx_] :> SARAH`Delta[i,Private`idx],
              V[Private`idx_,i] :> SARAH`Delta[Private`idx,i]}];

Print["testing ReplaceMixingMatrixByIdentityIn[] ..."];

(* check that mixing matrix is replaced by Delta *)
couplingWithMixingMatrix = {Cp[UField[{idx}]], g ZField[idx, i] };

TestEquality[Private`ReplaceMixingMatrixByIdentityIn[couplingWithMixingMatrix[[2]],
                                                     couplingWithMixingMatrix[[1]]],
             g SARAH`Delta[idx,i]
            ];

(* check that if no mixing matrix exists, a sum over the mixing matrix
   is introduced *)
couplingWithDelta = {Cp[Susyno`LieGroups`conj[UField[{idx1}]], Field[{idx2}]],
                     g SARAH`Delta[idx1, idx2] };

TestEquality[Private`ReplaceMixingMatrixByIdentityIn[couplingWithDelta[[2]],
                                                     couplingWithDelta[[1]]],
             SARAH`sum[idx, 1, 3, g ZField[idx,idx1] SARAH`Delta[idx,idx2]]
            ];

(* check case where the conj is on the other field *)
couplingWithDelta = {Cp[UField[{idx1}], Susyno`LieGroups`conj[Field[{idx2}]]],
                     g SARAH`Delta[idx1, idx2] };

TestEquality[Private`ReplaceMixingMatrixByIdentityIn[couplingWithDelta[[2]],
                                                     couplingWithDelta[[1]]],
             SARAH`sum[idx, 1, 3, g Susyno`LieGroups`conj[ZField[idx,idx1]] SARAH`Delta[idx,idx2]]
            ];

PrintTestSummary[];

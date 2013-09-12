Needs["TestSuite`", "TestSuite.m"];
Needs["SelfEnergies`", "SelfEnergies.m"];

Print["testing CreateMixingMatrixReplacementRulesFor[] ..."];

Field;
VField;
FindMixingMatrixSymbolFor[Field] := ZField;
FindMixingMatrixSymbolFor[SVDField] := {U, V};
GetDimension[Field] := 3;
IsVector[VField] = True;

TestEquality[SelfEnergies`Private`CreateMixingMatrixReplacementRulesFor[Cp[Field[{i}]]], {}];

TestEquality[SelfEnergies`Private`CreateMixingMatrixReplacementRulesFor[Cp[UField[{i}]]],
             {ZField[i,SelfEnergies`Private`idx_] :> SARAH`Delta[i,SelfEnergies`Private`idx],
              ZField[SelfEnergies`Private`idx_,i] :> SARAH`Delta[SelfEnergies`Private`idx,i]}];

TestEquality[SelfEnergies`Private`CreateMixingMatrixReplacementRulesFor[Cp[USVDField[{i}]]],
             {U[i,SelfEnergies`Private`idx_] :> SARAH`Delta[i,SelfEnergies`Private`idx],
              U[SelfEnergies`Private`idx_,i] :> SARAH`Delta[SelfEnergies`Private`idx,i],
              V[i,SelfEnergies`Private`idx_] :> SARAH`Delta[i,SelfEnergies`Private`idx],
              V[SelfEnergies`Private`idx_,i] :> SARAH`Delta[SelfEnergies`Private`idx,i]}];

Print["testing ReplaceMixingMatrixByIdentityIn[] ..."];

(* check that mixing matrix is replaced by Delta *)
couplingWithMixingMatrix = {Cp[UField[{idx}]], g ZField[idx, i] };

TestEquality[SelfEnergies`Private`ReplaceMixingMatrixByIdentityIn[couplingWithMixingMatrix[[2]],
                                                     couplingWithMixingMatrix[[1]]],
             g SARAH`Delta[idx,i]
            ];

(* check that if no mixing matrix exists, a sum over the mixing matrix
   is introduced *)
couplingWithDelta = {Cp[Susyno`LieGroups`conj[UField[{idx1}]], VField, Field[{idx2}]],
                     g SARAH`Delta[idx1, idx2] };

TestEquality[SelfEnergies`Private`ReplaceMixingMatrixByIdentityIn[couplingWithDelta[[2]],
                                                     couplingWithDelta[[1]]],
             SARAH`sum[idx, 1, 3, g ZField[idx,idx1] SARAH`Delta[idx,idx2]]
            ];

(* check case where the conj is on the other field *)
couplingWithDelta = {Cp[UField[{idx1}], VField, Susyno`LieGroups`conj[Field[{idx2}]]],
                     g SARAH`Delta[idx1, idx2] };

TestEquality[SelfEnergies`Private`ReplaceMixingMatrixByIdentityIn[couplingWithDelta[[2]],
                                                     couplingWithDelta[[1]]],
             SARAH`sum[idx, 1, 3, g Susyno`LieGroups`conj[ZField[idx,idx1]] SARAH`Delta[idx,idx2]]
            ];

PrintTestSummary[];

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

BeginPackage["Vertices`", {
    "SARAH`",
    "Parameters`",
    "TreeMasses`",
    "LatticeUtils`"}]

VertexRules::usage;
ToCpPattern::usage="ToCpPattern[cp] converts field indices inside cp to patterns, e.g. ToCpPattern[Cp[bar[UFd[{gO1}]], Sd[{gI1}], Glu[{1}]][PL]] === Cp[bar[UFd[{gO1_}]], Sd[{gI1_}], Glu[{1}]][PL].";
ToCp::usage="ToCp[cpPattern] converts field index patterns inside cpPattern to symbols, e.g. ToCp@Cp[bar[UFd[{gO1_}]], Sd[{gI1_}], Glu[{1}]][PL] === Cp[bar[UFd[{gO1}]], Sd[{gI1}], Glu[{1}]][PL].";
FieldIndexList::usage;

GetLorentzStructure::usage;
GetParticleList::usage;
IsUnrotated::usage;
ToRotatedField::usage;
ReplaceUnrotatedFields::usage;

Begin["`Private`"]

VertexRules[nPointFunctions_, massMatrices_] := Block[{
	UnitaryMatrixQ,
	nCpPatterns,
	cpPatterns = DeleteRedundantCpPatterns[
	    ToCpPattern /@ RenumberCpIndices /@ Cases[
		nPointFunctions, _SARAH`Cp|_SARAH`Cp[_], {0, Infinity}]]
    },
    (UnitaryMatrixQ[#] = True)& /@
        Flatten@DeleteCases[massMatrices[[All,3]], Null];
    UnitaryMatrixQ[_] := False;
    nCpPatterns = Length[cpPatterns];
    MapIndexed[
	DoneLn[#1 -> VertexExp[#1, nPointFunctions, massMatrices],
	       "[",First[#2],"/",nCpPatterns,"] calculating ", #1, "... "]&,
	cpPatterns]
];

DeleteRedundantCpPatterns[cpPatterns_] :=
    First @ Sort[#, MatchQ[#2, #1]&]& /@
    Gather[cpPatterns, MatchQ[#1, #2] || MatchQ[#2, #1]&];

VertexExp[cpPattern_, nPointFunctions_, massMatrices_] := Module[{
	cp = ToCp[cpPattern],
	rotatedCp, fieldsInRotatedCp,
	sarahVertex,
	fields, vertices,
	lorentzTag, lorentz, vertex,
	strippedIndices,
	contraction
    },
    rotatedCp = ReplaceUnrotatedFields[cp];
    fieldsInRotatedCp = GetParticleList[rotatedCp];
    sarahVertex = SARAH`Vertex[fieldsInRotatedCp];
    Assert[MatchQ[sarahVertex, {_, __}]];
    fields = First[sarahVertex];
    vertices = Rest[sarahVertex];
    lorentzTag = GetLorentzStructure[rotatedCp];
    {vertex, lorentz} = FindVertexWithLorentzStructure[vertices, lorentzTag];
    strippedIndices = Complement[Flatten[FieldIndexList /@ fields],
				 Flatten[FieldIndexList /@ fieldsInRotatedCp]];
    vertex = StripGroupStructure[
	ResolveColorFactor[
	    vertex, fields, cpPattern, nPointFunctions[[All,2]]],
	strippedIndices];
    contraction = Block[{
	    SARAH`sum
	    (* corrupts a polynomial (monomial + monomial + ...) summand *)
	},
	ExpandSarahSum @ SimplifyContraction @
	InTermsOfRotatedVertex[
	    vertex, lorentz,
	    GetParticleList[cp], massMatrices]];
    (* Q: is the factor -I right? *)
    -I TreeMasses`ReplaceDependencies[contraction] /.
	Parameters`ApplyGUTNormalization[]
];

StripGroupStructure[expr_, indices_List] := Module[{
	indexPattern = Alternatives@@indices
    },
    expr /. {
	PatternWithDistinctIndices[SARAH`Delta, 2, indexPattern] -> 1,
	SARAH`Lam[__] -> 2,
	SARAH`Sig[__] -> 2,
	SARAH`fSU2[__] -> 1,
	SARAH`fSU3[__] -> 1
    }
];

PatternWithDistinctIndices[head_, nIndices_, indexPattern_] :=
ReleaseHold[Block[{
	patternNames = Table[Unique[], {nIndices}],
	UnsameQ
    },
    With[{condition = UnsameQ @@ patternNames},
	 head @@ (Pattern[#, indexPattern]& /@ patternNames) /;
	 Hold[condition]]]];

FindVertexWithLorentzStructure[vertices_List, lorentz_Integer] :=
    vertices[[lorentz]];

FindVertexWithLorentzStructure[vertices_List, lorentz_] :=
    SingleCase[vertices, {_, str_ /; !FreeQ[str, lorentz]}];

simplifyContractionDispatch = Dispatch[{
    HoldPattern[
	SARAH`sum[i_, 1, d_, x_ z_?UnitaryMatrixQ[i_, j_]]
	/; HasContractionQ[x, Susyno`LieGroups`conj[z[i, _]]] &&
	SARAH`getDim[z] === d] :>
    (x /. Susyno`LieGroups`conj[z[i, k_]] :> SARAH`Delta[j, k]),

    HoldPattern[
	SARAH`sum[i_,1,d_, x_ Susyno`LieGroups`conj[z_?UnitaryMatrixQ[i_,j_]]]
	/; HasContractionQ[x, z[i, _]] &&
	SARAH`getDim[z] === d] :>
    (x /. z[i, k_] :> SARAH`Delta[j, k]),

    HoldPattern[
	SARAH`sum[i_, 1, d_, x_ z_?UnitaryMatrixQ[j_, i_]]
	/; HasContractionQ[x, Susyno`LieGroups`conj[z[_, i]]] &&
	SARAH`getDim[z] === d] :>
    (x /. Susyno`LieGroups`conj[z[k_, i]] :> SARAH`Delta[j, k]),

    HoldPattern[
	SARAH`sum[i_,1,d_, x_ Susyno`LieGroups`conj[z_?UnitaryMatrixQ[j_,i_]]]
	/; HasContractionQ[x, z[_, i]] &&
	SARAH`getDim[z] === d] :>
    (x /. z[k_, i] :> SARAH`Delta[j, k]),

    s:HoldPattern[
    (SARAH`sum[i_, 1, d_, x_ z_?UnitaryMatrixQ[i_, j_]]
     /; HasMixedContractionQ[x, Susyno`LieGroups`conj[z[i, _]]]) |
    (SARAH`sum[i_, 1, d_, x_ Susyno`LieGroups`conj[z_?UnitaryMatrixQ[i_, j_]]]
     /; HasMixedContractionQ[x, z[i, _]]) |
    (SARAH`sum[i_, 1, d_, x_ z_?UnitaryMatrixQ[j_, i_]]
     /; HasMixedContractionQ[x, Susyno`LieGroups`conj[z[_, i]]]) |
    (SARAH`sum[i_, 1, d_, x_ Susyno`LieGroups`conj[z_?UnitaryMatrixQ[j_, i_]]]
     /; HasMixedContractionQ[x, z[_, i]]) /;
    SARAH`getDim[z] === d] :> ExpandSarahSum[s]
}];

SimplifyContraction[expr_] := expr //. simplifyContractionDispatch;

HasContractionQ[x_, form_] :=
    !MatchQ[FindContraction[x, form], Unindexed|Indexed|Mixed];

HasMixedContractionQ[x_, form_] := FindContraction[x, form] === Mixed;

FindContraction[x_, form:_?UnitaryMatrixQ[___, i_Symbol,___]] :=
    FindContraction[x, form, i];

FindContraction[x_, form:HoldPattern@Susyno`LieGroups`conj
		[_?UnitaryMatrixQ[___, i_Symbol, ___]]] :=
    FindContraction[x, form, i];

FindContraction[x_Times, form_, index_] := Module[{
	structure = FindContraction[#, form, index]& /@ List@@x,
	match
    },
    Which[MatchQ[structure, {Unindexed..}], Unindexed,
	  (match =
	   Replace[structure,
		   {{Unindexed... , p:Except[Unindexed], Unindexed...} -> p,
		    _ -> False}]) =!= False, match,
	  MatchQ[structure, {___, Mixed, ___}], Mixed,
	  MatchQ[structure, {___, Except[Unindexed], ___}], Indexed,
	  True, Print["Vertices`Private`FindContraction[",
		      x, ", ", form, ", ", index, "] failed."]; Abort[]]
];

FindContraction[x_Plus, form_, index_] := Module[{
	structure = FindContraction[#, form, index]& /@ List@@x,
	match
    },
    Which[MatchQ[structure, {Unindexed..}], Unindexed,
	  MatchQ[structure, {(f:Except[Unindexed])..}], First[structure],
	  MatchQ[structure, {___, Except[Unindexed], ___}], Mixed,
	  True, Print["Vertices`Private`FindContraction[",
		      x, ", ", form, ", ", index, "] failed."]; Abort[]]
];

FindContraction[HoldPattern@SARAH`sum[_, _, _, x_], form_, index_] :=
    FindContraction[x, form, index];

FindContraction[x_, form_, index_] /; FreeQ[x, index] := Unindexed;

FindContraction[x_, form_, index_] /; MatchQ[x, form] := x;

FindContraction[x_, form_, index_] := Indexed;

ExpandSarahSum[expr_] := expr //.
    SARAH`sum[a_, b_, c_, x_] /; Head[Expand[x]] === Plus :>
    (SARAH`sum[a, b, c, #]& /@ Expand[x]);

InTermsOfRotatedVertex[vertex_, lorentz_, uFields_List, massMatrices_] :=
Block[{
	fixedFields,
	SARAH`bar
    },
    (* reconstruct hidden bar applied on a Majorana spinor *)
    fixedFields = MapIndexed[
	If[MajoranaQ@FieldHead@ToRotatedField[#1] && First[#2] <= 2,
	   SARAH`bar[#1], #1]&,
	uFields];
    Fold[If[IsUnrotated[#2],
	    RewriteUnrotatedField[
		#1, lorentz, #2,
		SingleCase[massMatrices,
			   _[_,FieldHead@ToRotatedField[#2],z_] :> z]],
	    #1]&,
	 vertex, fixedFields]
];

RewriteUnrotatedField[
    expr_, _,
    uField:_Symbol[{__}], z_Symbol] :=
    ContractMixingMatrix[expr, uField, z, Identity];

RewriteUnrotatedField[
    expr_, _,
    uField:HoldPattern@Susyno`LieGroups`conj[_Symbol[{__}]], z_Symbol] :=
    ContractMixingMatrix[expr, uField, z, Susyno`LieGroups`conj];

RewriteUnrotatedField[
    expr_, LorentzProduct[gamma[_], PL],
    uField:_Symbol[{__}], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, u, Identity];

RewriteUnrotatedField[
    expr_, LorentzProduct[gamma[_], PL],
    uField:HoldPattern@SARAH`bar[_Symbol[{__}]], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, u, Susyno`LieGroups`conj];

RewriteUnrotatedField[
    expr_, LorentzProduct[gamma[_], PR],
    uField:_Symbol[{__}], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, v, Susyno`LieGroups`conj];

RewriteUnrotatedField[
    expr_, LorentzProduct[gamma[_], PR],
    uField:HoldPattern@SARAH`bar[_Symbol[{__}]], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, v, Identity];

RewriteUnrotatedField[
    expr_, PL,
    uField:_Symbol[{__}], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, u, Identity];

RewriteUnrotatedField[
    expr_, PL,
    uField:HoldPattern@SARAH`bar[_Symbol[{__}]], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, v, Identity];

RewriteUnrotatedField[
    expr_, PR,
    uField:_Symbol[{__}], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, v, Susyno`LieGroups`conj];

RewriteUnrotatedField[
    expr_, PR,
    uField:HoldPattern@SARAH`bar[_Symbol[{__}]], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, u, Susyno`LieGroups`conj];

ContractMixingMatrix[expr_, uField_, z_, op_] := Module[{
	uIndices = FieldIndexList[uField],
	field = ToRotatedField[uField],
	uIndex, i = Unique["gl"]
    },
    (* CHECK: does the first index always denote the flavour? *)
    uIndex = First[uIndices];
    SARAH`sum[i, SARAH`getGenStart[field], SARAH`getGen[field],
	      (expr /. uIndex -> i) op @ z[i, uIndex]]
];

FieldHead[Susyno`LieGroups`conj[field_]] := FieldHead[field];

FieldHead[SARAH`bar[field_]] := FieldHead[field];

FieldHead[field_Symbol[{__}]] := field;

FieldHead[field_Symbol] := field;

ToCpPattern[cp : _SARAH`Cp|_SARAH`Cp[_]] := cp /.
    ((# -> (# /. Thread[(# -> If[Head[#] === Symbol, Pattern[#, _], #])& /@
			FieldIndexList[#]]))& /@
     GetParticleList[cp]);

ToCp[cpPattern : _SARAH`Cp|_SARAH`Cp[_]] := cpPattern /. p_Pattern :> First[p];

CpType[cp : _SARAH`Cp|_SARAH`Cp[_]] := RotatedCpType @
    ReplaceUnrotatedFields[cp];

RotatedCpType[SARAH`Cp[fields__]] := SARAH`getVertexType[{fields}];

RotatedCpType[SARAH`Cp[fields__][_]] := SARAH`getVertexType[{fields}];

RenumberCpIndices[SARAH`Cp[fields__]] :=
    SARAH`Cp @@ RenumberFieldIndices[{fields}]

RenumberCpIndices[SARAH`Cp[fields__][l_]] :=
    (SARAH`Cp @@ RenumberFieldIndices[{fields}])[l]

RenumberFieldIndices[fields_List] := Block[{
	UsedSarahIndexQ
    },
    UsedSarahIndexQ[_] := False;
    RenumberFieldIndices /@ fields
];

RenumberFieldIndices[field_] :=
    field /. ((# -> RenumberSarahIndex[#])& /@ FieldIndexList[field]);

RenumberSarahIndex[index_Symbol?UsedSarahIndexQ] :=
    RenumberSarahIndex @ Symbol @
	StringReplace[ToString[index],
		      RegularExpression["[[:digit:]]+$"] :>
		      ToString[ToExpression["$0"]+1]];

RenumberSarahIndex[index_Symbol] := (UsedSarahIndexQ[index] = True; index);

RenumberSarahIndex[index_] := index;

(* see SelfEnergies`Private`CreateCouplingFunctions[] *)
ResolveColorFactor[vertex_, fields_, cpPattern_, exprs_] :=
    If[UnresolvedColorFactorFreeQ[cpPattern, exprs],
       vertex,
       Module[{loopArgs,
	       internalColorIndices = InternalColorIndices[fields],
	       externalColorIndices = ExternalColorIndices[fields]},
	      (* Q: does one need to sum also over external color indices
		 as in SelfEnergies`Private`CreateCouplingFunctions[]?
		 A: it is a way to strip the color structure of this class
		 of vertices *)
	   loopArgs = Join[{#, 3}& /@ internalColorIndices,
			   {#, 1}& /@ externalColorIndices];
	   Sum @@ Prepend[loopArgs, vertex]
       ]];

InternalColorIndices[fields_List] :=
    Union@Cases[DeleteCases[FieldIndexList /@ fields,
			    {___,_?SarahExternalGenerationIndexQ,___}],
		_?SarahColorIndexQ, Infinity];

ExternalColorIndices[fields_List] :=
    Union@Cases[Cases[FieldIndexList /@ fields,
		      {___,_?SarahExternalGenerationIndexQ,___}],
		_?SarahColorIndexQ, Infinity];

FieldIndexList[field_] := Flatten@Cases[field, _?VectorQ, {0, Infinity}];

UnresolvedColorFactorFreeQ[cpPattern_, exprs_] := Module[{
	fstPos = First@Position[exprs, cpPattern],
	cpInstance,
	exprInstance
    },
    cpInstance = Extract[exprs, fstPos];
    exprInstance = Extract[exprs, Take[fstPos, 1]];
    FreeQ[Coefficient[exprInstance //. SARAH`sum[__, ex_] :> ex, cpInstance],
	  C]
];

(* CHECK: are the following right semantics of SARAH indices? *)
SarahExternalGenerationIndexQ[index_Symbol] :=
    StringMatchQ[ToString[index], RegularExpression["gO[[:digit:]]+"]];

SarahInternalGenerationIndexQ[index_Symbol] :=
    StringMatchQ[ToString[index], RegularExpression["gI[[:digit:]]+"]];

SarahColorIndexQ[index_Symbol] :=
    StringMatchQ[ToString[index], RegularExpression["ct[[:digit:]]+"]];

GetLorentzStructure[SARAH`Cp[__]] := 1;

GetLorentzStructure[SARAH`Cp[__][a_]] := a;

GetParticleList[SARAH`Cp[a__]] := {a};

GetParticleList[SARAH`Cp[a__][_]] := {a};

IsUnrotated[SARAH`bar[field_]] := IsUnrotated[field];

IsUnrotated[Susyno`LieGroups`conj[field_]] := IsUnrotated[field];

IsUnrotated[field_[__]] := IsUnrotated[field];

IsUnrotated[field_Symbol] := StringTake[ToString[field],1] === "U";

ToRotatedField[field_Symbol] :=
    Symbol[StringReplace[ToString[field], StartOfString ~~ "U" ~~ rest_ :> rest]];

ToRotatedField[SARAH`bar[field_]] := SARAH`bar[ToRotatedField[field]];

ToRotatedField[Susyno`LieGroups`conj[field_]] := Susyno`LieGroups`conj[ToRotatedField[field]];

ToRotatedField[field_List] := ToRotatedField /@ field;

ToRotatedField[field_[indices__]] := ToRotatedField[field][indices];

ReplaceUnrotatedFields[SARAH`Cp[p__]] :=
    SARAH`Cp @@ ToRotatedField[{p}];

ReplaceUnrotatedFields[SARAH`Cp[p__][lorentz_]] :=
    ReplaceUnrotatedFields[SARAH`Cp[p]][lorentz];

End[] (* `Private` *)

EndPackage[]

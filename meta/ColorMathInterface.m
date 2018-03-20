(* ::Package:: *)

(* Outline
1.  I assume that the input for the calculation is a list from Jobst 
	having the following format:
	{list of external fields, list of lists of fields in vertices}
	e.g.
	{Fe,bar[Fe],VP,{bar[Fe],Ah,Fe},{Fe,Ah,bar[Fe]},{VP,bar[Fe],Fe}}
	Sizes of those lists are arbitrary.
2.  For each external color charged field I generate a random name for 
	the color index. This is done through
	GenerateUniqueColorAssociationsForExternalParticles function
	which takes list from Jobst and returns a association list,
	e.g. \[LeftAssociation]1\[Rule]c22\[RightAssociation] meaning that external particle at position on in
	Jobst list will have color index c22.	
*)

BeginPackage["ColorMathInterface`", {"SARAH`", "ColorMath`"}]

CalculateColorFactor::usage =
	"dasdas"
   
(* Begin["`Private`"]; *)
ccParticles={Fd,bar[Fd],VG};
ColorChargedQ[particle_] :=
 Module[{p = particle /. bar[x__] -> x},
  If[Head[p] === Symbol, MemberQ[ccParticles, p],
   MemberQ[ccParticles, Head[p]]
   ]
  ]
  
(* give a field, e.g. Fd[{a,b}] or bar[Fd[{a,b}] will return {a,b} *)
GetFieldIndices[field_] :=
  field //. bar[x__] -> x /. _[x_List] :> x 
  
GetFieldColorIndex[field_/;ColorChargedQ[field]]:=
  Module[{res},
    res=GetFieldIndices[field];
    res = Select[res,ColorIndexQ];
    Print[res];
    Assert[Length[res]==1];
    res[[1]]
  ]
CalculateColorFactor[vertex_List, rule_] :=
   Module[{return},
      return = 
         vertex // 
         DropColorles //
         TakeOnlyColor // 
         SARAHToColorMathSymbols;
         return = Times @@ return;
         Print[return];
         (*return = (#/.(x___ SARAH`Delta[col1_,col2_]\[RuleDelayed] x/.col2\[Rule]col1)&)/@return;*)
         return=return //. (x___ SARAH`Delta[col1_, col2_]y___:> (x y /. col2 -> col1));
      CSimplify[return]
   ];
    
ColorIndexQ::notes="Checks if a field index is a color index. Color indices start with 'c'"
ColorIndexQ[x_Symbol] :=
   (Characters@SymbolName[x])[[1]] == "c"
   
GenerationIndexQ[x_Symbol] :=
   (Characters@SymbolName[x])[[1]] == "g"

LorentzIndexQ[x_Symbol] :=
   (Characters@SymbolName[x])[[1]] == "l"
DropColorles::notes = "Drop colorles vertices from the list of Vertex objets  "
DropColorles[vertices_List] :=  
   DeleteCases[vertices, el_ /; 
      FreeQ[el, 
         SARAH`Lam[__] |
         SARAH`fSU3[__] | 
         SARAH`Delta[c1_/;ColorIndexQ[c1], c2_/;ColorIndexQ[c2]]
      ]
   ];

TakeOnlyColor[v__] :=
(Take[
   #, {2}][[1, 1]] & /@ v) /. ___ Lam[colIdx__] :> Lam[colIdx]/. ___ fSU3[colIdx__] :> fSU3[colIdx] /. ___ SARAH`Delta[c1_/;ColorIndexQ[c1], c2_/;ColorIndexQ[c2]] :> SARAH`Delta[c1,c2];

SARAHToColorMathSymbols[s_] := 
   s /.SARAH`Lam[colIdx1_, colIdx2_, colIdx3_] :> 2 ColorMath`t[{colIdx1}, colIdx2, colIdx3]/. SARAH`fSU3[colSeq__] :> ColorMath`f[colSeq];

(* input
    {Fe,bar[Fe],VP,{bar[Fe],Ah,Fe},{Fe,Ah,bar[Fe]},{VP,bar[Fe],Fe}}
    *)
GenerateUniqueColorAssociationsForExternalParticles::notes=
  "Generates unique color indices for external particles"    
GenerateUniqueColorAssociationsForExternalParticles[v_List]:=
  Module[{inOutParticles,inOutColoredParticles,a},
    inOutParticles=TakeWhile[v,(Head[#]=!=List)&];
    (* generate a unique color index for every external particle *)
    a = Association[{}];
    inOutParticlesWithColorIndices = 
    MapIndexed[
      If[ColorChargedQ[#1],AssociateTo[a,#2[[1]]->Unique["c"]] ]&,
inOutParticles
];
a
    ]

Print[GenerateColorIndices[
  {Fe,bar[Fe],VP,{bar[Fe],Ah,Fe},{Fe,Ah,bar[Fe]},{VP,bar[Fe],Fe}}
  ]]
(* End[] *)

EndPackage[];










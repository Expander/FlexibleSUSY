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
3.  Final step. Take fields in vertices and pass them to SARAH`Vertex function.
    Each vertex will get automatically generated indices starting from c1. 
    Using info the info from Jobst adjacency matrix and my association list 
    create new dummny indices and connect them between vertices.		
*)

BeginPackage["ColorMathInterface`", {"SARAH`", "CXXDiagrams`", "ColorMath`"}];

RegenerateIndices[l_List, graph_]:=
   Module[{keys, extFields},
      keys = GenerateUniqueColorAssociationsForExternalParticles[l];
      extFields = TakeWhile[l,(Head[#]=!=List)&];
      Print[extFields, " ", CXXDiagrams`ColorChargedQ /@ extFields];
      vertices = Drop[l,Length@extFields];
      ll = SARAH`Vertex[#]& /@ vertices;
      (* loop over external particles *)
      For[extIdx=1, extIdx <= Length[extFields], extIdx++,
         (* skip if uncollored *)
         If[!CXXDiagrams`ColorChargedQ[extFields[[extIdx]]], Continue[]];
(* loop over vertices *)
For[vertIdx=1,vertIdx<=Length[Complement[l,extFields]],vertIdx++,
(* check graph if enternal field is connected to the vertex at all *)
If[graph[[extIdx,vertIdx+Length[extFields]]]==0,Continue[]];
(* loop over particles in the vertex *)
For[vertFieldIdx=1,vertFieldIdx<=Length[ll[[vertIdx,1]]],vertFieldIdx++,
pInV=l[[vertIdx+Length[extFields],vertFieldIdx]];
If[!CXXDiagrams`ColorChargedQ[pInV],Continue[]];
If[AntiField[extFields[[extIdx]]]=!= pInV,Continue[]];
ll=MapAt[(
(#//.GetFieldColorIndex[#[[1,vertFieldIdx]]]->keys[extIdx]))&,ll,vertIdx]];
]
]
(* loop over vertices pairs *)
For[vertIdx1 = 1,vertIdx1<=Length[vertices],vertIdx1++,
For[vertIdx2=vertIdx1+1, vertIdx2<=Length[vertices], vertIdx2++,
(* if two vertices are not conneted at all *)
If[graph[[vertIdx1+Length[extFields],vertIdx2+Length[extFields]]]==0,Continue[]];
(* loop over fields in the vertex *)
For[v1i=1,v1i<=Length[vertices[[vertIdx1,1]]],v1i++,
If[!CXXDiagrams`ColorChargedQ[ll[[vertIdx1,1,v1i]]],Continue[]];
For[v2i=1,v2i<=Length[vertices[[vertIdx2,1]]],v2i++,
If[!CXXDiagrams`ColorChargedQ[ll[[vertIdx2,1,v2i]]],Continue[]];
If[ll[[v1,1,v1i]]=!=AntiField[ll[[v2,1,v2i]]],Continue[]];
ll=MapAt[(#//.GetFieldColorIndex[#[[1,v2i]]]:>GetFieldColorIndex[ll[[vertIdx1,1,v1i]]])&,ll,vertIdx2];
]
]
]
];
ll
];

CalculateColorFactor::usage =
	"dasdas"
   
(* Begin["`Private`"]; *)

(* give a field, e.g. Fd[{a,b}] or bar[Fd[{a,b}] will return {a,b} *)
GetFieldIndices[field_] :=
  field //. bar[x__] -> x /. _[x_List] :> x 
  
GetFieldColorIndex[field_/;CXXDiagrams`ColorChargedQ[field]]:=
  Module[{res},
    res=GetFieldIndices[field];
    res = Select[res,ColorIndexQ];
    Assert[Length[res]==1];
    res[[1]]
  ]
CalculateColorFactor[vertex_List,graph_] :=
   Module[{return},
      return = 
         RegenerateIndices[vertex,graph] // DropColorles;
      If[ return === {}, Return[1]];
      return = 
         return //  TakeOnlyColor // 
         SARAHToColorMathSymbols;
      Print[return];
      return = Times @@ return;
      return = return //. (x___ SARAH`Delta[col1_, col2_] y___ :> (x y /. col2 -> col1));
      return = return //. x___ SARAH`Delta[col1_, col2_] y___ :> x y ColorMath`delta[col1, col2];
      (* CSimplify[1] doesn't evaluate *)
      If[ return === 1, 1, Return[AllSimpleRules[return]]];
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
   Module[{vert},
   vert = DeleteCases[vertices, el_ /; 
      FreeQ[el, 
         SARAH`Lam[__] |
         SARAH`fSU3[__] | 
         SARAH`Delta[c1_/;ColorIndexQ[c1], c2_/;ColorIndexQ[c2]]
      ]
   ];
   vert
   ]
   

TakeOnlyColor[v__] :=
    Module[{result},
      (* the generic structure of the Vertex "object" is 
         {{ParticleList},{{Coefficient 1, Lorentz 1},{Coefficient 2, Lorentz 2},...} *)
      (* drop ParticleList *) 
      (*Print["start --------------------------------------------------------------------------------------------------------------"];*)
      result = Drop[#, 1]& /@ v;
      (*Print["1: ", result];*)
      result = (Transpose @ Drop[Transpose[#], -1])& /@ result;
      (*Print["2: ", result];*)
      result = result //. 
         ___ SARAH`Lam[colIdx__] :> SARAH`Lam[colIdx] //. 
         ___ SARAH`fSU3[colIdx__] :> SARAH`fSU3[colIdx] //. 
         ___ SARAH`Delta[c1_/;ColorIndexQ[c1], c2_/;ColorIndexQ[c2]] :> SARAH`Delta[c1,c2];
      (*Print["3: ", result];*)
      result = DeleteCases[#, {0}]& /@ result;
      Assert[CountDistinct[#] === 1]& /@ result;
      result = DeleteDuplicates[#]& /@ result;
      (*Print["4: ", result];*)
      result = Flatten[result, 2];
      (*Print["5: ", result];*)
      (*Print["end --------------------------------------------------------------------------------------------------------------"];*)
      result
    ];

SARAHToColorMathSymbols[s__] := s //.
   SARAH`Lam[colIdx1_, colIdx2_, colIdx3_] :> 2 ColorMath`t[{colIdx1}, colIdx2, colIdx3] //. 
   SARAH`fSU3[colSeq__] :> ColorMath`f[colSeq];

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
      If[CXXDiagrams`ColorChargedQ[#1],AssociateTo[a,#2[[1]]->Unique["c"]] ]&,
inOutParticles
];
a
    ]

(* End[] *)

EndPackage[];










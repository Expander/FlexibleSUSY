(* ::Package:: *)

(*
Outline
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

Limitations:
1.  We asumme that there are only 3-particle vertices.
*)

BeginPackage["ColorMathInterface`",
   {"SARAH`", "TreeMasses`", (*IsColorIndex is there*)"CXXDiagrams`", "ColorMath`"}
];

CalculateColorFactor::usage = "";
SARAHToColorMathSymbols::usage = "";
GetFieldColorIndex::usage = "";

Begin["`Private`"];

(*
RegenerateIndices[l_List, graph_]:=
    Module[{keys, extFields, particlesInVertices,
            vertices, vertex, vertex1, vertex2,
            extField, inField,
            fieldsInVertex, fieldsInVertex1, fieldsInVertex2},

        extFields = TakeWhile[l, (Head[#]=!=List)&];
        Print["Are external fields charged: ", extFields, " ", TreeMasses`ColorChargedQ /@ extFields];
        keys = GenerateUniqueColorAssociationsForExternalParticles[l];
        Print["Generated color indices for external particles: ", keys];
        particlesInVertices = Drop[l, Length@extFields];
        (* change to vertices = SARAH`Vertex /@ particlesInVertices; *)
        vertices = SARAH`Vertex[#]& /@ particlesInVertices;

        (* STEP1: set colors of external particles in vertices to values in keys *)
        (* loop over external particles *)
        For[extIdx=1, extIdx <= Length[extFields], extIdx++,
            extField = extFields[[extIdx]];
            (* skip if un-colored *)
            If[!TreeMasses`ColorChargedQ[extField], Continue[]];
            (* loop over vertices *)
            For[vertIdx=1, vertIdx<=Length[vertices], vertIdx++,
                vertex = vertices[[vertIdx]];
                fieldsInVertex = vertex[[1]];
                (* check if enternal field is connected to the vertex *)
                If[graph[[extIdx, vertIdx+Length[extFields]]] === 0, Continue[]];
                (* loop over particles in the vertex *)
                For[vertFieldIdx=1, vertFieldIdx<=Length[fieldsInVertex], vertFieldIdx++,
                    inField = fieldsInVertex[[vertFieldIdx]];
                    If[!TreeMasses`ColorChargedQ[inField], Continue[]];
                    If[AntiField[extField] =!= (inField /. f_[_List] -> f), Continue[]];
                    If[MemberQ[Values[keys], GetFieldColorIndex[inField]], Continue[]];
                    vertices = MapAt[
                        (# //. GetFieldColorIndex[inField] :> keys[extIdx])&,
                        vertices,
                        vertIdx
                    ];
                ];
            ];
        ];

        symbol = {};
        (* loop over vertices pairs *)
        vertex1 = vertices[[1]];
        fieldsInVertex1 = vertex1[[1]];
        field1ColorIndexOld = GetFieldColorIndex[fieldsInVertex1[[2]]];
        field1ColorIndexNew = Unique["c"];
        vertices = MapAt[
        (# //. field1ColorIndexOld -> field1ColorIndexNew)&,
        vertices, 1
        ];
        vertices = MapAt[
          (# //. GetFieldColorIndex[vertices[[2,1,2]]] -> field1ColorIndexNew)&,
          vertices, 2
        ];
        field1ColorIndexOld = GetFieldColorIndex[fieldsInVertex1[[3]]];
        field1ColorIndexNew = Unique["c"];
        vertices = MapAt[
          (# //. field1ColorIndexOld -> field1ColorIndexNew)&,
          vertices, 1
        ];
        vertices = MapAt[
          (# //. GetFieldColorIndex[vertices[[3,1,2]]] -> field1ColorIndexNew)&,
          vertices, 3
        ];
        vertices = MapAt[
          (# //. GetFieldColorIndex[vertices[[2,1,3]]] -> GetFieldColorIndex[vertices[[3,1,3]]])&,
          vertices, 3
        ];

        (*For[vertIdx1 = 1, vertIdx1<=Length[vertices], vertIdx1++,*)
            (*vertex1 = vertices[[vertIdx1]];*)
            (*fieldsInVertex1 = vertex1[[1]];*)
            (** loop over fields in vertex1 *)
            (*For[v1i=1, v1i<=Length[fieldsInVertex1], v1i++,*)
                (*field1 = fieldsInVertex1[[v1i]];*)
                (*If[!TreeMasses`ColorChargedQ[field1], Continue[]];*)
                (*field1ColorIndexOld = GetFieldColorIndex[field1];*)
                (** cycle if external particle *)
                (*If[MemberQ[Values[keys], field1ColorIndexOld], Continue[]];*)
                (*field1ColorIndexNew = Unique["c"];*)
                (*Print["generate new index ", field1ColorIndexNew, " for ", field1];*)
                (*vertices = MapAt[*)
                    (*(# //. field1ColorIndexOld -> field1ColorIndexNew)&,*)
                    (*vertices, vertIdx1*)
                (*];*)
                (*AppendTo[symbol, field1ColorIndexNew];*)
                (*For[vertIdx2=vertIdx1+1, vertIdx2<=Length[vertices], vertIdx2++,*)
                    (** cycle if two vertices are not connected *)
                    (*If[graph[[vertIdx1+Length[extFields], vertIdx2+Length[extFields]]] === 0, Continue[]];*)
                    (*vertex2 = vertices[[vertIdx2]];*)
                    (*fieldsInVertex2 = vertex2[[1]];*)
                    (*For[v2i=1, v2i<=Length[fieldsInVertex2], v2i++,*)
                        (*field2 = fieldsInVertex2[[v2i]];*)
                        (*Print["Second field ", field2];*)
                        (*If[!TreeMasses`ColorChargedQ[field2], Continue[]];*)
                        (*field2ColorIndex = GetFieldColorIndex[field2];*)
                        (*If[MemberQ[Values[keys], field2ColorIndex], Continue[]];*)
                        (*If[MemberQ[symbol, field2ColorIndex], Print["Ania......."];Continue[]];*)
                        (** we want to make a propagator < field bar[field]> *)
                        (*If[(field1 /. f_[_List] -> f) =!= (AntiField[field2] /. f_[_List] -> f), Continue[]];*)
                        (*Print["Overriding ", field2ColorIndex,  " ", field1ColorIndexNew];*)
                        (*If[vertIdx1 === 4 && v1i === 2 && vertIdx2 =!= 5, Continue[]];*)
                        (*vertices = MapAt[*)
                            (*(# //. field2ColorIndex -> field1ColorIndexNew)&,*)
                            (*vertices, vertIdx2*)
                        (*];*)
                    (*];*)
                (*];*)
            (*];*)
        (*];*)
        vertices
   ];
*)

(*  given a field will return it's indices,
    e.g. Fd[{a,b}] or bar[Fd[{a,b}] or conj[Sd[{a,b}]] will return {a,b} *)
GetFieldIndices[field_] :=
    field /. SARAH`bar | Susyno`LieGroups`conj -> Identity /. _[x_List] :> x;

(*  given a field with indices will discard indices,
    e.g. bar[Fd[{a,b}] -> bar[Fd] *)
DropFieldIndices[field_] :=
    field /. f_[_List] -> f;

NumberOfExternalParticles[l_List] :=
    Count[l, Except[_List], 1];

NumberOfVertices[l_List] :=
    Length[l] - NumberOfExternalParticles[l];

GetFieldColorIndex[field_/;TreeMasses`ColorChargedQ[field]]:=
  Module[{res},
    res = GetFieldIndices[field];
    Print[IsColorIndex /@ res];
    res = Select[res, IsColorIndex];
    Assert[Length[res]==1];
    res[[1]]
  ];
  
CalculateColorFactor[vertex_List] :=
   Module[{return},
      return =
         vertex // DropColorles;
      If[ return === {}, Return[1]];
      return =
         TakeOnlyColor@return;
      return = Times @@ return;
      (* CSimplify[1] doesn't evaluate *)
      If[return === 1,
         1,
         CSimplify[return]
      ]
   ];

GenerationIndexQ[x_Symbol] :=
   (Characters@SymbolName[x])[[1]] == "g";

LorentzIndexQ[x_Symbol] :=
   (Characters@SymbolName[x])[[1]] == "l";

ColorStructureFreeQ[el_] :=
   FreeQ[el,
      Subscript[Superscript[Superscript[ColorMath`t,List[__]],_],_] |
         Superscript[ColorMath`f,List[__]] |
         Subscript[Superscript[ColorMath`\[Delta], _], _] |
         Superscript[ColorMath`\[CapitalDelta],List[_, _]]
   ];

DropColorles::notes = "Drop colorles vertices from the list of Vertex objets  ";
DropColorles[vertices_List] :=  
   DeleteCases[vertices,
      el_ /; ColorStructureFreeQ[el]
   ];


TakeOnlyColor[vvvv__] :=
    Module[{result},
      (* the generic structure of the Vertex "object" is 
         {{ParticleList},{{Coefficient 1, Lorentz 1},{Coefficient 2, Lorentz 2},...} *)
      (* drop ParticleList *) 
      result = Drop[#, 1]& /@ vvvv;
      result = (Transpose @ Drop[Transpose[#], -1])& /@ result;
      result = result //.
         __?ColorStructureFreeQ el_ /; !ColorStructureFreeQ[el] :> el;
      result = DeleteCases[#, {0}]& /@ result;
      Assert[CountDistinct[#] === 1]& /@ result;
      result = DeleteDuplicates[#]& /@ result;
      result = Flatten[result, 2];
      result
    ];

SARAHToColorMathSymbols[vertex_List] :=
    Module[{result},

   result = vertex //.
      SARAH`Lam[colIdx1_, colIdx2_, colIdx3_] :> 2 ColorMath`t[{colIdx1}, colIdx2, colIdx3] //.
      SARAH`fSU3[colSeq__] :> ColorMath`f[colSeq];

   (* if the result has Delta, we need to find out if it's adj. or fundamental *)
    result = result /. SARAH`Delta[c1_?IsColorIndex, c2_?IsColorIndex] :>
        If[getColorRep[Select[vertex[[1]], ! FreeQ[#, c1] &][[1]]] === T,
           ColorMath`delta[c1, c2]
        ];
   result = result /. SARAH`Delta[c1_?IsColorIndex, c2_?IsColorIndex] :>
       If[getColorRep[Select[vertex[[1]], ! FreeQ[#, c1] &][[1]]] === O,
          ColorMath`Delta[c1, c2]
       ];

       result
   ];

(* input
    {Fe,bar[Fe],VP,{bar[Fe],Ah,Fe},{Fe,Ah,bar[Fe]},{VP,bar[Fe],Fe}}
    *)
(*
GenerateUniqueColorAssociationsForExternalParticles::notes=
  "Generates unique color indices for external particles"    
GenerateUniqueColorAssociationsForExternalParticles[vvvv_List]:=
  Module[{inOutParticles,inOutColoredParticles,a},
    inOutParticles=TakeWhile[vvvv,(Head[#]=!=List)&];
    (* generate a unique color index for every external particle *)
    a = Association[{}];
    inOutParticlesWithColorIndices = 
    MapIndexed[
      If[TreeMasses`ColorChargedQ[#1],AssociateTo[a,#2[[1]]->Unique["c"]] ]&,
inOutParticles
];
a
    ]
*)

End[];

EndPackage[];

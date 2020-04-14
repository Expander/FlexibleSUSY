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

BeginPackage["MatMulSimplify`"];

MatMulSimplify::usage =
"MatMulSimplify[expr] returns the tuple { simpExpr, simpRules }, where
simpExpr is a simplified expression where all (nested) appearances of
MatMul[___] are replaced by unique symbols.  simpRules contains the
replacement rules for these symbols to restore the expression

Example:

    MatMulSimplify[MatMul[a,b] + MatMul[a,b,c] + MatMul[a,b,c,d], MatMul]
";

MatMulRefine::usage =
"MatMulRefine[expr] refines all MatMul[___] expressions in expr by
pulling out numerical terms.

Example:

    MatMulRefine[MatMul[a] + MatMul[a,2], MatMul]
";

Begin["`Private`"];

(* returns expressions matching pattern *)
MatMulFindPattern[expr_, patt_] :=
    DeleteDuplicates @ Cases[expr, patt, {0, Infinity}]

(* replace rules in expressions, including nestings of the form head[___, term, ___] *)
MatMulReplaceSubExpr[expr_, rules_, head_] :=
    Module[{a, b},
           expr //. rules //. MatMulRefine[Rule[head[a___, PatternSequence @@ First[#], b___], head[a, Last[#], b]]& /@ rules, head]
    ]

(* creates pattern of the form head[_,...,_] with n blanks *)
MatMulMakePattern[n_Integer, head_] :=
    head[PatternSequence @@ Array[Blank[]&, n]]

(* returns tuple of simplified expr and replacement rules *)
MatMulReplaceNested[expr_, head_, {nMin_, nMax_}, ruleSymb_] :=
    Module[{simpRep = {}, simpExpr = expr, rules, patt, m},
           While[!FreeQ[simpExpr, head],
                 For[m = nMin, !FreeQ[simpExpr, head] && m <= nMax, m++,
                     patt = MatMulMakePattern[m, head];
                     rules = Rule[#, Unique[ruleSymb]]& /@ MatMulFindPattern[simpExpr, patt];
                     simpExpr = MatMulReplaceSubExpr[simpExpr, rules, head];
                     simpRep = Join[simpRep, rules];
                     If[rules =!= {},
                        Break[];
                     ];
                 ];
           ];
           {MatMulRefine[simpExpr, head], Reverse /@ simpRep}
    ]

MatMulRefine[expr_, head_] :=
    expr //. {
	head[a_] :> a,
	head[a__, n_?NumericQ, b___] :> head[n, a, b],
	head[n_?NumericQ, a___] :> n head[a],
	head[a___, head[m__], b___] :> head[a, m, b]
    }

MatMulSimplify[expr_, head_, ruleSymb_:"mat"] :=
    MatMulReplaceNested[MatMulRefine[expr, head], head, {2, Infinity}, ruleSymb]

End[];

EndPackage[];

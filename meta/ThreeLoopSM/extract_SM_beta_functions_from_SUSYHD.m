(*
   This script extracts the SM beta-functions from SUSYHD.

   Author: Alexander Voigt

   Run it like this:

   math -run "Get[extract_SM_beta_functions_from_SUSYHD.m]; Quit"
*)

Needs["SUSYHD`"];

Begin["Private`"]

outputDir = ".";

NoSubscript[sym_] := 
  Simplify[sym /. 
    Subscript[a_, b_] :> Symbol[ToString[a] <> ToString[b]]];

WriteOut[sym_, name_, factor_, repl_:{}] :=
  Module[{l1, l2, l3, \[Kappa] = 1/(4 \[Pi])^2, filename},
    l1 = SM\[Beta][1, sym];
    l2 = SM\[Beta][2, sym] - l1;
    l3 = SM\[Beta][3, sym] - l2;
    filename = FileNameJoin[{outputDir, "beta_" <> name <> ".m"}];
    Print["Writing beta-function for ", NoSubscript[sym], " to ", filename];
    Put[{\[Kappa]^-1 l1, \[Kappa]^-2 l2, \[Kappa]^-3 l3}*factor /. 
      repl // NoSubscript, filename];
 ];

(*
   Convert to SARAH convention.
   SARAH : L = - \[Lambda]/2 H^4
           d g / d Log[Q] = \[Beta]_SARARH

   SUSYHD: L = - \[Lambda]   H^4
           d g^2 / d Log[Q^2] = \[Beta]_SUSYHD
*)

WriteOut[Subscript[g, 1]     , "g1"    , 1/Subscript[g, 1]     , \[Lambda] -> \[Lambda]/2];
WriteOut[Subscript[g, 2]     , "g2"    , 1/Subscript[g, 2]     , \[Lambda] -> \[Lambda]/2];
WriteOut[Subscript[g, 3]     , "g3"    , 1/Subscript[g, 3]     , \[Lambda] -> \[Lambda]/2];
WriteOut[\[Lambda]           , "lambda", 4/(2 \[Lambda])       , \[Lambda] -> \[Lambda]/2];
WriteOut[Subscript[g, t]     , "gt"    , 1/Subscript[g, t]     , \[Lambda] -> \[Lambda]/2];
WriteOut[Subscript[g, b]     , "gb"    , 1/Subscript[g, b]     , \[Lambda] -> \[Lambda]/2];
WriteOut[Subscript[g, \[Tau]], "gtau"  , 1/Subscript[g, \[Tau]], \[Lambda] -> \[Lambda]/2];

End[]

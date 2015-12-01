(*
   This script extracts the 3-loop MSSM beta-functions from
   http://www.liv.ac.uk/~dij/betas/allgennb.log .

   Author: Alexander Voigt

   Run it like this:

   math -run "<< extract_MSSM_beta_functions_from_hep-ph-0308231.m; Quit[]"
*)

betaStr = Import["allgennb.log"];

(* fix syntax error *)
betaStr = StringReplace[betaStr, "beta_htau^(3)) =" -> "beta_htau^(3) ="];
    
(* wrap expressions around brackets *)
betaStr = StringReplace[betaStr, {"=" -> "= (", ";" -> ");"}];

(* convert symbols into valid Wolfram symbols *)
betaStr = StringReplace[betaStr,
   {"_" -> "", "^(1)" -> "1", "^(2)" -> "2", "^(3)" -> "3"}];

matrixPattern = (
    (WordBoundary ~~ "Ye" ~~ WordBoundary) | "Yec" |
    (WordBoundary ~~ "Yt" ~~ WordBoundary) | "Ytc" |
    (WordBoundary ~~ "Yb" ~~ WordBoundary) | "Ybc" |
    (WordBoundary ~~ "he" ~~ WordBoundary) | "hec" |
    (WordBoundary ~~ "ht" ~~ WordBoundary) | "htc" |
    (WordBoundary ~~ "hb" ~~ WordBoundary) | "hbc" |
    (WordBoundary ~~ "mq" ~~ WordBoundary) |
    (WordBoundary ~~ "mt" ~~ WordBoundary) |
    (WordBoundary ~~ "mb" ~~ WordBoundary) |
    (WordBoundary ~~ "ml" ~~ WordBoundary) |
    (WordBoundary ~~ "me" ~~ WordBoundary));

(* convert matrix products to non-commutative products *)
betaStr = StringReplace[betaStr,
   {"tr(" ~~ p:(Except[")"]..) ~~ ")" :>
    "trace[" <> StringReplace[p, "*" -> ","] <> "]",
    p:(matrixPattern ~~ ("*" ~~ matrixPattern)..) :>
    "MatMul[" <> StringReplace[p, "*" -> ","] <> "]"}];

(* convert adjoint symbols to SARAH convention *)
betaStr = StringReplace[betaStr,
   {"Yec" -> "Adj[Ye]",
    "Ytc" -> "Adj[Yt]",
    "Ybc" -> "Adj[Yb]",
    "hec" -> "Adj[he]",
    "htc" -> "Adj[ht]",
    "hbc" -> "Adj[hb]" }];

ToExpression[betaStr];

repl = {
    a1 -> g1^2,
    a2 -> g2^2,
    a3 -> g3^2,
    k -> 6 Zeta[3],
    n5 -> 0,
    n10 -> 0
};

{g1 bg11  , g1 bg12  , g1 bg13  } /. repl >> "beta_g1.m";
{g2 bg21  , g2 bg22  , g2 bg23  } /. repl >> "beta_g2.m";
{g3 bg31  , g3 bg3   , g3 bg33  } /. repl >> "beta_g3.m";
{betat1   , betat2   , betat3   } /. repl >> "beta_Yu.m";
{betab1   , betab2   , betab3   } /. repl >> "beta_Yd.m";
{betae1   , betae2   , betae3   } /. repl >> "beta_Ye.m";
{betamq1  , betamq2  , betamq3  } /. repl >> "beta_mq2.m";
{betamt1  , betamt2  , betamt3  } /. repl >> "beta_mu2.m";
{betamb1  , betamb2  , betamb3  } /. repl >> "beta_md2.m";
{betamL1  , betamL2  , betamL3  } /. repl >> "beta_ml2.m";
{betamtau1, betamtau2, betamtau3} /. repl >> "beta_me2.m";
{betamH21 , betamH22 , betamH23 } /. repl >> "beta_mHu2.m";
{betamH11 , betamH12 , betamH13 } /. repl >> "beta_mHd2.m";
{betaht1  , betaht2  , betaht3  } /. repl >> "beta_TYu.m";
{betahb1  , betahb2  , betahb3  } /. repl >> "beta_TYd.m";
{betahtau1, betahtau2, betahtau3} /. repl >> "beta_TYe.m";
{betaM11  , betaM12  , betaM13  } /. repl >> "beta_M1.m";
{betaM21  , betaM22  , betaM23  } /. repl >> "beta_M2.m";
{betaM31  , betaM32  , betaM33  } /. repl >> "beta_M3.m";

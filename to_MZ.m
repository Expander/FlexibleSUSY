(* ****** NMSSM ****** *)

tempRules = {
    T[\[Lambda]] -> Alam lam,
    T[\[Kappa]]  -> Akap kap,
    \[Lambda]    -> lam,
    \[Kappa]     -> kap
};

(* EWSB solution, obtained with FlexibleSUSY *)
solution = {vS -> (Sign[vS]*Sqrt[(-40*mHd2*vd^2 - 3*g1^2*vd^4 - 5*g2^2*vd^4 + 
      40*mHu2*vu^2 + 3*g1^2*vu^4 + 5*g2^2*vu^4 + 40*vd*tadpole[1] - 
      40*vu*tadpole[2])/(vd^2*\[Lambda]*conj[\[Lambda]] - 
      vu^2*\[Lambda]*conj[\[Lambda]])])/(2*Sqrt[5]), 
 \[Kappa] -> (40*mHd2*vd + 3*g1^2*vd^3 + 5*g2^2*vd^3 - 3*g1^2*vd*vu^2 - 
    5*g2^2*vd*vu^2 + 20*vd*vS^2*\[Lambda]*conj[\[Lambda]] + 
    20*vd*vu^2*\[Lambda]*conj[\[Lambda]] - 10*Sqrt[2]*vS*vu*
     conj[T[\[Lambda]]] - 10*Sqrt[2]*vS*vu*T[\[Lambda]] - 40*tadpole[1])/
   (10*vS^2*vu*(\[Lambda] + conj[\[Lambda]])), 
 ms2 -> (-4*vS^3*\[Kappa]^2 + 2*vd*vS*vu*\[Kappa]*\[Lambda] + 
    2*vd*vS*vu*\[Kappa]*conj[\[Lambda]] - 2*vd^2*vS*\[Lambda]*
     conj[\[Lambda]] - 2*vS*vu^2*\[Lambda]*conj[\[Lambda]] - 
    Sqrt[2]*vS^2*conj[T[\[Kappa]]] + Sqrt[2]*vd*vu*conj[T[\[Lambda]]] - 
    Sqrt[2]*vS^2*T[\[Kappa]] + Sqrt[2]*vd*vu*T[\[Lambda]] + 4*tadpole[3])/
   (4*vS)};

(* simplify solution *)
solution = solution /. conj[x_] :> x /. tempRules;

(* combine to 1 equation *)
mS2 = ms2 /. solution[[3]] /. solution[[2]] /. solution[[1]];

(* simplify expression for mS2 *)
mS2 = Refine[mS2,
             vu > 0 && vd > 0 && g1 > 0 && g2 > 0 && tanBeta > 1 && Element[vu, Reals] && Element[vS, Reals] && Element[vd, Reals] && Element[lam, Reals] && Element[tanBeta, Reals] && MZ2 > 0
            ];

simplificationRules = {
    tadpole[_] -> 0
    , g1^2       -> (MZ2/((vu^2 + vd^2)/4) - g2^2 ) 5/3
    , vu         -> tanBeta vd
    , Sign[_]^2  -> 1
    , Sign[_]^-2 -> 1
    (* expressions (1+tanBeta^n) and (-1+tanBeta^n) appear often -> abbreviate *)
    , Plus[n_, Power[tanBeta, p_]] :> TB[n,p]
    };

mS2 = mS2 //. simplificationRules;
mS2 = Simplify[mS2];
mS2 = mS2 //. simplificationRules;
mS2 = Simplify[mS2];
mS2 = mS2 //. simplificationRules;

mS2 >> "NMSSM_EWSB_solution.m";

Print["mS2 == ", InputForm[mS2]];

Quit[0];


(* ******* MSSM ******* *)

solution = {B[\[Mu]] -> (-20*mHd2*vd*vu + 20*mHu2*vd*vu - 3*g1^2*vd^3*vu - 
    5*g2^2*vd^3*vu + 3*g1^2*vd*vu^3 + 5*g2^2*vd*vu^3 + 20*vu*tadpole[1] - 
    20*vd*tadpole[2])/(20*(vd^2 - vu^2)), 
 \[Mu] -> (Sign[\[Mu]]*Sqrt[-(mHd2*vd) - (3*g1^2*vd^3)/40 - (g2^2*vd^3)/8 + 
      (3*g1^2*vd*vu^2)/40 + (g2^2*vd*vu^2)/8 + vu*B[\[Mu]] + tadpole[1]])/
   Sqrt[vd]}

(* combine to 1 equation *)
mu2 = \[Mu]^2 /. solution[[2]] /. solution[[1]];

(* simplify solution *)
simplificationRules = {
    tadpole[_] -> 0,
    conj[x_]   :> x,
    g1^2       -> (MZ2/((vu^2 + vd^2)/4) - g2^2 ) 5/3,
    vu         -> Tan[beta] vd,
    Sign[_]^2  -> 1
    };

mu2 = mu2 //. simplificationRules;
mu2 = FullSimplify[mu2];

mu2 >> "MSSM_EWSB_solution.m";

Print["mu^2 == ", InputForm[mu2]];

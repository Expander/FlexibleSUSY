BeginPackage["RGIntegrator`"];
EndPackage[];

RGIntegrate::usage = "Integrates a given list of renormalization group
 equations from a starting scale to a destination scale.

 The first parameter of RGIntegrate[] is the list of beta functions.
 Each list element is again a list, where the first entry is the name
 of the parameter and the rest are the 1-loop, 2-loop, ... beta
 functions.

 The loopOrder of the integration can be given as an optional
 argument, see Example 3, (default: Automatic).

 Example 1
 =========

 betas = { {g, h a g^2} }

 RGIntegrate[betas, Q1, Q2]

 Result:

 {g[Q1] -> g[Q2] - a h g[Q2]^2 Log[Q2/Q1]}

 Example 2
 =========

 betas = { {g, h a g^2, h^2 b g^2 l^2},
           {l, h (c l^2 + d g^2), h^2 (e l^4 + f g^2 l^2)} }

 RGIntegrate[betas, Q1, Q2]

 Result:

 {g[Q1] -> g[Q2] - a h g[Q2]^2 Log[Q2/Q1] +
    h^2 (a^2 g[Q2]^3 Log[Q1/Q2]^2 - b g[Q2]^2 l[Q2]^2 Log[Q2/Q1]),
  l[Q1] -> l[Q2] + h (-(d g[Q2]^2 Log[Q2/Q1]) - c l[Q2]^2 Log[Q2/Q1]) +
    h^2 ((a d g[Q2]^3 + c d g[Q2]^2 l[Q2] + c^2 l[Q2]^3) Log[Q1/Q2]^2 -
      f g[Q2]^2 l[Q2]^2 Log[Q2/Q1] - e l[Q2]^4 Log[Q2/Q1])}

 Example 3
 =========

 Integrate the beta functions up to the 2nd order:

 betas = { {g, h a g^2} }

 RGIntegrate[betas, Q1, Q2, loopOrder -> 2]

 Result:

 {g[Q1] -> g[Q2] - a h g[Q2]^2 Log[Q2/Q1] + a^2 h^2 g[Q2]^3 Log[Q1/Q2]^2}
";

Begin["RGIntegrator`Private`"];

(* loop factor for expansion *)
{h};

AddScale[par_, Q_] := Rule[par, par[Q]];

collectLogs := {
    Plus[a___, x_ Log[Q1_], b___, y_ Log[Q2_], c___] /; x === -y :> Plus[a, b, c, x Log[Q1/Q2]]
};

MultiplyLoopFactor[{par_, betas___}, h_] :=
    Join[{par}, MapIndexed[#1 h^(First[#2])&, {betas}]];

PerformIntegrals[expr_] := expr /. Integral -> Integrate;

IntegrateSingleRHS[{par_, betas___}, Q1_, Q2_, Qp_, addScales_, sol_] :=
    Module[{loopOrder = Length[{betas}], integrand},
           integrand = Total[{betas} /. addScales /. (sol /. Q1 -> Qp)]/Qp;
           integrand = Normal @ Series[integrand, {h,0,loopOrder}];
           par[Q1] -> par[Q2] - Integral[integrand, {Qp, Q1, Q2},
                                         Assumptions :> Q1 > 0 && Q2 > 0 && Qp > 0 && Q2 > Q1]
          ];

IntegrateRHS[{}, Q1_, Q2_, sol_] := {};

IntegrateRHS[betas_List, Q1_, Q2_, sol_] :=
    Module[{Qp, addScales, ints},
           addScales = AddScale[First[#], Qp]& /@ betas;
           ints = IntegrateSingleRHS[#, Q1, Q2, Qp, addScales, sol]& /@ betas;
           PerformIntegrals /@ ints
          ];

(* integrate beta functions of order h^N,
   given a solution of the (N-1)th order beta functions *)
RGIntegrateRecursively[betas_List, Q1_, Q2_, sol_] :=
    IntegrateRHS[betas, Q1, Q2, sol];

(* solution of empty beta functions *)
RGIntegrateRecursively[{}, Q1_, Q2_] := {};

(* solution of trivial beta functions *)
RGIntegrateRecursively[betas : {{par_}, ___}, Q1_, Q2_] :=
    Rule[First[#][Q1], First[#][Q2]]& /@ betas;

RGIntegrateRecursively[betas_List, Q1_, Q2_] :=
    Module[{sol, betasLower},
           betasLower = Drop[#, -1]& /@ betas;
           sol = RGIntegrateRecursively[betasLower, Q1, Q2];
           RGIntegrateRecursively[betas, Q1, Q2, sol]
      ];

Options[RGIntegrate] = {
    loopOrder -> Automatic
};

RGIntegrate[beta_List, Q1_, Q2_, OptionsPattern[]] :=
    Module[{lbeta, lo = OptionValue[loopOrder]},
           lbeta = MultiplyLoopFactor[#, h]& /@ beta;
           If[IntegerQ[lo],
              (* chop or fill beta functions with zeros *)
              lbeta = PadRight[#, lo + 1]& /@ lbeta;
             ];
           RGIntegrateRecursively[lbeta, Q1, Q2] /. h -> 1 /. collectLogs
          ];

End[];

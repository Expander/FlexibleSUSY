(* Global symbols *)
{A0, B0, B1, B22, B22tilde, C0, D0, D27, F, G, H};

BeginPackage["LoopFunctions`"];

FullMomentum::usage = "Returns loop functions with full momentum
dependence in BPMZ convention (arxiv:hep-ph/9606211)";

ZeroMomentum::usage = "Returns loop functions at zero momentum in BPMZ
 convention (arxiv:hep-ph/9606211)";

Delta::usage = "1/eps - \[Gamma]_E + Log[4 Pi]";

Begin["`Private`"];

ZeroMomentum[] := {
                                        Global`A0[m_,mu_]               :> A0impl[m,mu],
    Global`B0[0,m1_,m2_,mu_]          | Global`B0[m1_,m2_,mu_]          :> B0zero[m1,m2,mu],
    Global`B1[0,m1_,m2_,mu_]          | Global`B1[m1_,m2_,mu_]          :> B1zero[m1,m2,mu],
    Global`B00[0,m1_,m2_,mu_]         | Global`B00[m1_,m2_,mu_]         :> B22zero[m1,m2,mu],
    Global`B11[0,m1_,m2_,mu_]         | Global`B11[m1_,m2_,mu_]         :> B11zero[m1,m2,mu],
    Global`B22[0,m1_,m2_,mu_]         | Global`B22[m1_,m2_,mu_]         :> B22zero[m1,m2,mu],
    Global`B22tilde[0,m1_,m2_,mu_]    | Global`B22tilde[m1_,m2_,mu_]    :> B22tildezero[m1,m2,mu],
    Global`C0[0,m1_,m2_,m2_,mu_]      | Global`C0[m1_,m2_,m2_,mu_]      :> C0zero[m1,m2,m3,mu],
    Global`D0[0,m1_,m2_,m3_,m4_,mu_]  | Global`D0[m1_,m2_,m3_,m4_,mu_]  :> D0zero[m1,m2,m3,m4,mu],
    Global`D27[0,m1_,m2_,m3_,m4_,mu_] | Global`D27[m1_,m2_,m3_,m4_,mu_] :> D27zero[m1,m2,m3,m4,mu],
    Global`F[0,m1_,m2_,mu_]           | Global`F[m1_,m2_,mu_]           :> Fzero[m1,m2,mu],
    Global`G[0,m1_,m2_,mu_]           | Global`G[m1_,m2_,mu_]           :> Gzero[m1,m2,mu],
    Global`H[0,m1_,m2_,mu_]           | Global`H[m1_,m2_,mu_]           :> Hzero[m1,m2,mu]
};

FullMomentum[] := {
    Global`A0[m_,mu_]                  :> A0impl[m,mu],
    Global`B0[p_,m1_,m2_,mu_]          :> B0[p,m1,m2,mu],
    Global`B1[p_,m1_,m2_,mu_]          :> B1[p,m1,m2,mu],
    Global`B00[p_,m1_,m2_,mu_]         :> B22[p,m1,m2,mu],
    Global`B11[p_,m1_,m2_,mu_]         :> B11[p,m1,m2,mu],
    Global`B22[p_,m1_,m2_,mu_]         :> B22[p,m1,m2,mu],
    Global`B22tilde[p_,m1_,m2_,mu_]    :> B22tilde[p,m1,m2,mu],
    Global`C0[p_,m1_,m2_,m2_,mu_]      :> C0[p,m1,m2,m3,mu],
    Global`D0[p_,m1_,m2_,m3_,m4_,mu_]  :> D0[p,m1,m2,m3,m4,mu],
    Global`D27[p_,m1_,m2_,m3_,m4_,mu_] :> D27[p,m1,m2,m3,m4,mu],
    Global`F[p_,m1_,m2_,mu_]           :> F[p,m1,m2,mu],
    Global`G[p_,m1_,m2_,mu_]           :> G[p,m1,m2,mu],
    Global`H[p_,m1_,m2_,mu_]           :> H[p,m1,m2,mu]
};

(* A0, Eq. (B.5) *)
A0impl[m_, mu_] := m^2 (Delta + 1 + Log[mu^2/m^2]);

(* A0, Eq. (B.5) *)
B0zero[m1_, m2_, mu_] := If[PossibleZeroQ[m1 - m2],
   Delta + Log[mu^2/m2^2],
   Delta + 1 + (m1^2 Log[mu^2/m1^2] - m2^2 Log[mu^2/m2^2])/(m1^2 - m2^2)
   ];

(* Eq. (B.6) *)
B0[p_, m1_, m2_, mu_] :=
    If[PossibleZeroQ[p],
       B0zero[m1,m2,mu],
       B0analytic[p,m1,m2,mu]
      ];

B0analytic[p_, m1_, m2_, mu_] :=
    Module[{eps, fB, s, xp, xm},
           s = p^2 - m2^2 + m1^2;
           xp = (s + Sqrt[s^2 - 4 p^2 (m1^2 - I eps)]) / (2 p^2);
           xm = (s - Sqrt[s^2 - 4 p^2 (m1^2 - I eps)]) / (2 p^2);
           fB[x_] := Log[1-x] - x Log[1 - 1/x] - 1;
           Limit[Delta - Log[p^2/mu^2] - fB[xp] - fB[xm], eps -> 0]
          ];

B0integral[p_, m1_, m2_, mu_] :=
    Module[{eps},
           Limit[
               Delta - Integrate[Log[((1-x) m1^2 + x m2^2 - x (1-x) p^2 - I eps)/mu^2], {x,0,1}],
               eps -> 0]
          ];

(* Eq. (B.9) *)
B1[p_, m1_, m2_, mu_] :=
    If[PossibleZeroQ[p],
       B1zero[m1,m2,mu],
       1/(2 p^2) (A0impl[m2,mu] - A0impl[m1,mu] + (p^2 + m1^2 - m2^2) B0[p,m1,m2])
      ];

B1zero[m1_, m2_, mu_] :=
    If[PossibleZeroQ[m1 - m2],
       (Delta + Log[mu^2/m2^2])/2,
       1/2 (Delta + 1 + Log[mu^2/m2^2] + (m1^2/(m1^2 - m2^2))^2 Log[m2^2/m1^2]
            + 1/2 (m1^2 + m2^2)/(m1^2 - m2^2))
      ];

B11[p_, m1_, m2_, mu_] :=
    If[PossibleZeroQ[p],
       B11zero[m1,m2,mu],
       1/(6 p^2) (
           2 A0impl[m2,mu]
           - 2 m1^2 B0[p,m1,m2,mu]
           - 4 (p^2 - m2^2 + m1^2) B1[p,m1,m2,mu]
           - m1^2 - m2^2 - p^2/3)
      ];

B11zero[m1_, m2_, mu_] := 0;

(* Eq. (B.10), identical to B00[p,m1,m2,mu] *)
B22[p_, m1_, m2_, mu_] :=
    If[PossibleZeroQ[p],
       B22zero[m1,m2,mu],
       1/6 (1/2 (A0impl[m1,mu] + A0impl[m2,mu])
            + (m1^2 + m2^2 - p^2/2) B0[p,m1,m2,mu]
            + (m2^2 - m1^2)/(2 p^2) (A0impl[m2,mu] - A0impl[m1,mu]
                                     - (m2^2 - m1^2) B0[p,m1,m2,mu])
            + m1^2 + m2^2 - p^2/3)
      ];

B22zero[m1_, m2_, mu_] :=
    If[PossibleZeroQ[m1 - m2],
       (m2^2*(1 + Delta + Log[mu^2/m2^2]))/2,
       ((5 + 3*Delta)*(m1^4 - m2^4) + m1^2*(3*m1^2 + m2^2)*Log[mu^2/m1^2] - 
        m2^2*(m1^2 + 3*m2^2)*Log[mu^2/m2^2])/(12*(m1^2 - m2^2))
      ];

(* Eq. (C.19) *)
C0zero[m1_, m2_, m3_] := Which[
    PossibleZeroQ[m1 - m2] && PossibleZeroQ[m1 - m3],
    -1/(2*m3^2),
    PossibleZeroQ[m1 - m2],
    (-m2^2 + m3^2 - m3^2*Log[m3^2/m2^2])/(m2^2 - m3^2)^2
    PossibleZeroQ[m1 - m3],
    (m2^2 - m3^2 - m2^2*Log[m2^2/m3^2])/(m2^2 - m3^2)^2,
    PossibleZeroQ[m2 - m3],
    (m1^2 - m3^2 + m1^2*Log[m3^2/m1^2])/(m1^2 - m3^2)^2,
    (* general case *)
    True,
    (+ m2^2 / (m1^2 - m2^2) Log[m2^2/m1^2]
     - m3^2 / (m1^2 - m3^2) Log[m3^2/m1^2]) / (m2^2 - m3^2)
];

(* Eq. (C.21) *)
D0zero[m1_, m2_, m3_, m4_] := Which[
   PossibleZeroQ[m1 - m2] && PossibleZeroQ[m1 - m3] && 
    PossibleZeroQ[m1 - m4],
   1/(6 m1^4),
   PossibleZeroQ[m1 - m2] && PossibleZeroQ[m2 - m3],
   (8*(m2^2 - m4^2)^2*(m2^2 + m4^2) + (-23*m2^4*m4^2 - 10*m2^2*m4^4 + 
         m4^6)*Log[m2^2/m4^2] + (-7*m2^4*m4^2 - 26*m2^2*m4^4 + m4^6)*
       Log[m4^2/m2^2])/(16*m2^2*(m2^2 - m4^2)^4),
   PossibleZeroQ[m1 - m2] && PossibleZeroQ[m2 - m4],
   (8*(m2^2 - m3^2)^2*(m2^2 + m3^2) + (7*m2^4*m3^2 + 26*m2^2*m3^4 - 
         m3^6)*Log[m2^2/m3^2] + (23*m2^4*m3^2 + 10*m2^2*m3^4 - m3^6)*
       Log[m3^2/m2^2])/(16*m2^2*(m2^2 - m3^2)^4),
   PossibleZeroQ[m1 - m2] && PossibleZeroQ[m3 - m4],
   (-8*(m2^2 - m4^2)^2 + (m2^4 - 8*m2^2*m4^2 - 5*m4^4)*
       Log[m2^2/m4^2] - (3*m2^4 + 8*m2^2*m4^2 + m4^4)*
       Log[m4^2/m2^2])/(4*(m2^2 - m4^2)^4),
   PossibleZeroQ[m1 - m2],
   (-((m2 - m3)*(m2 + m3)*(m2 - m4)*(m3 - m4)*(m2 + m4)*(m3 + 
           m4)) - (m2^4 - m3^2*m4^2)*(m3^2*Log[m3^2/m2^2] + 
         m4^2*Log[m2^2/m4^2]) + 
      m3^2*m4^2*(-2*m2^2 + m3^2 + m4^2)*
       Log[m4^2/m3^2])/((m2^2 - m3^2)^2*(m2^2 - m4^2)^2*(m3^2 - m4^2)),
   True,
   (C0[m1, m3, m4] - C0[m2, m3, m4])/(m1^2 - m2^2)];

(* Eq. (C.22) *)
D27zero[m1_, m2_, m3_, m4_] :=
   (m1^2 C0[m1, m3, m4] - m2^2 C0[m2, m3, m4])/(4 (m1^2 - m2^2));

(* Eq. (B.11) *)
F[p_, m1_, m2_, mu_] :=
    A0impl[m1,mu] - 2 A0impl[m2,mu] - (2 p^2 + 2 m1^2 - m2^2) B0[p,m1,m2,mu];

Fzero[m1_, m2_, mu_] := F[0,m1,m2,mu];

(* Eq. (B.12) *)
G[p_, m1_, m2_, mu_] :=
    (p^2 - m1^2 - m2^2) B0[p,m1,m2,mu] - A0impl[m1,mu] - A0impl[m2,mu];

Gzero[m1_, m2_, mu_] := G[0,m1,m2,mu];

(* Eq. (B.13) *)
H[p_, m1_, m2_, mu_] := 4 B22[p,m1,m2,mu] + G[p,m1,m2,mu];

Hzero[m1_, m2_, mu_] := H[0,m1,m2,mu];

(* Eq. (B.14) *)
B22tilde[p_, m1_, m2_, mu_] := B22[p,m1,m2,mu] - A0impl[m1,mu]/4 - A0impl[m2,mu]/4;

B22tildezero[m1_, m2_, mu_] := B22tilde[0,m1,m2,mu];

End[];

EndPackage[];

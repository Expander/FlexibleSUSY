BeginPackage["LoopFunctions`"];
EndPackage[];

(* loop functions *)
{A0, B0, B1, B00, B11, B22, B22tilde, C0, D0, D27, F, G, H};

LFFull::usage = "Returns explicit form of loop functions with full
 momentum dependence.  The convention depends on the value of
 $BPMZSign .";

LFZeroMomentum::usage = "Returns loop functions at zero momentum.  The
 convention depends on the value of $BPMZSign .";

LFDivergence::usage = "Returns divergent part of loop functions.  The
 convention depends on the value of $BPMZSign .";

LFScaleDependence::usage = "Returns the renormalization scale
 dependent part of loop functions.  The convention depends on the
 value of $BPMZSign .";

Delta::usage = "1/eps - \[Gamma]_E + Log[4 Pi]";

$BPMZSign::usage = "If set to 1, loop functions are returned in BPMZ
 convention (arxiv:hep-ph/9606211) (default).  If set to -1, loop
 functions are returned in Denner convention
 (arxiv:hep-ph/0709.1075).";

$BPMZSign = 1;

Begin["LoopFunctions`Private`"];

LFZeroMomentum[] := {
    A0[m_,mu_]                                               :> A0impl[m,mu],
    B0[0,m1_,m2_,mu_]             | B0[m1_,m2_,mu_]          :> B0zero[m1,m2,mu],
    B1[0,m1_,m2_,mu_]             | B1[m1_,m2_,mu_]          :> B1zero[m1,m2,mu],
    B00[0,m1_,m2_,mu_]            | B00[m1_,m2_,mu_]         :> B22zero[m1,m2,mu],
    B11[0,m1_,m2_,mu_]            | B11[m1_,m2_,mu_]         :> B11zero[m1,m2,mu],
    B22[0,m1_,m2_,mu_]            | B22[m1_,m2_,mu_]         :> B22zero[m1,m2,mu],
    B22tilde[0,m1_,m2_,mu_]       | B22tilde[m1_,m2_,mu_]    :> B22tildezero[m1,m2,mu],
    C0[0,0,m1_,m2_,m3_,mu_]       | C0[m1_,m2_,m3_,mu_]      :> C0zero[m1,m2,m3],
    D0[0,0,0,m1_,m2_,m3_,m4_,mu_] | D0[m1_,m2_,m3_,m4_,mu_]  :> D0zero[m1,m2,m3,m4],
    D27[0,0,0,m1_,m2_,m3_,m4_,mu_]| D27[m1_,m2_,m3_,m4_,mu_] :> D27zero[m1,m2,m3,m4],
    F[0,m1_,m2_,mu_]              | F[m1_,m2_,mu_]           :> Fzero[m1,m2,mu],
    G[0,m1_,m2_,mu_]              | G[m1_,m2_,mu_]           :> Gzero[m1,m2,mu],
    H[0,m1_,m2_,mu_]              | H[m1_,m2_,mu_]           :> Hzero[m1,m2,mu]
};

LFFull[] := {
    A0[m_,mu_]                          :> A0impl[m,mu],
    B0[p_,m1_,m2_,mu_]                  :> B0impl[p,m1,m2,mu],
    B1[p_,m1_,m2_,mu_]                  :> B1impl[p,m1,m2,mu],
    B00[p_,m1_,m2_,mu_]                 :> B22impl[p,m1,m2,mu],
    B11[p_,m1_,m2_,mu_]                 :> B11impl[p,m1,m2,mu],
    B22[p_,m1_,m2_,mu_]                 :> B22impl[p,m1,m2,mu],
    B22tilde[p_,m1_,m2_,mu_]            :> B22tildeimpl[p,m1,m2,mu],
    C0[p1_,p2_,m1_,m2_,m3_,mu_]         :> C0impl[p1,p2,m1,m2,m3,mu],
    D0[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_] :> D0impl[p1,p2,p3,m1,m2,m3,m4,mu],
    D27[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_]:> D27impl[p1,p2,p3,m1,m2,m3,m4,mu],
    F[p_,m1_,m2_,mu_]                   :> Fimpl[p,m1,m2,mu],
    G[p_,m1_,m2_,mu_]                   :> Gimpl[p,m1,m2,mu],
    H[p_,m1_,m2_,mu_]                   :> Himpl[p,m1,m2,mu]
};

LFDivergence[] := {
    A0[m_,mu_]                          :> DivA0[m,mu],
    B0[p_,m1_,m2_,mu_]                  :> DivB0[p,m1,m2,mu],
    B1[p_,m1_,m2_,mu_]                  :> DivB1[p,m1,m2,mu],
    B00[p_,m1_,m2_,mu_]                 :> DivB22[p,m1,m2,mu],
    B11[p_,m1_,m2_,mu_]                 :> DivB11[p,m1,m2,mu],
    B22[p_,m1_,m2_,mu_]                 :> DivB22[p,m1,m2,mu],
    B22tilde[p_,m1_,m2_,mu_]            :> DivB22tilde[p,m1,m2,mu],
    C0[p1_,p2_,m1_,m2_,m3_,mu_]         :> 0,
    D0[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_] :> 0,
    D27[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_]:> 0,
    F[p_,m1_,m2_,mu_]                   :> DivF[p,m1,m2,mu],
    G[p_,m1_,m2_,mu_]                   :> DivG[p,m1,m2,mu],
    H[p_,m1_,m2_,mu_]                   :> DivH[p,m1,m2,mu]
};

LFScaleDependence[] := {
    A0[m_,mu_]                          :> LogA0[m,mu],
    B0[p_,m1_,m2_,mu_]                  :> LogB0[p,m1,m2,mu],
    B1[p_,m1_,m2_,mu_]                  :> LogB1[p,m1,m2,mu],
    B00[p_,m1_,m2_,mu_]                 :> LogB22[p,m1,m2,mu],
    B11[p_,m1_,m2_,mu_]                 :> LogB11[p,m1,m2,mu],
    B22[p_,m1_,m2_,mu_]                 :> LogB22[p,m1,m2,mu],
    B22tilde[p_,m1_,m2_,mu_]            :> LogB22tilde[p,m1,m2,mu],
    C0[p1_,p2_,m1_,m2_,m3_,mu_]         :> 0,
    D0[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_] :> 0,
    D27[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_]:> 0,
    F[p_,m1_,m2_,mu_]                   :> LogF[p,m1,m2,mu],
    G[p_,m1_,m2_,mu_]                   :> LogG[p,m1,m2,mu],
    H[p_,m1_,m2_,mu_]                   :> LogH[p,m1,m2,mu]
};

(* A0, Eq. (B.5) *)
A0impl[m_, mu_] :=
    If[PossibleZeroQ[m],
       0,
       m^2 (Delta + 1 + Log[mu^2/m^2])
      ];

DivA0[m_, mu_] := m^2 Delta;

LogA0[m_, mu_] :=
    If[PossibleZeroQ[m],
       0,
       m^2 Log[mu^2/m^2]
      ];

(* A0, Eq. (B.5) *)
B0zero[m1_, m2_, mu_] :=
    Which[PossibleZeroQ[m1 - m2],
          Delta + Log[mu^2/m2^2],
          PossibleZeroQ[m1],
          Delta + 1 + Log[mu^2/m2^2],
          PossibleZeroQ[m2],
          Delta + 1 + Log[mu^2/m1^2],
          True,
          Delta + 1 + (m1^2 Log[mu^2/m1^2] - m2^2 Log[mu^2/m2^2])/(m1^2 - m2^2)
   ];

(* Eq. (B.6) *)
B0impl[p_, m1_, m2_, mu_] :=
    Which[PossibleZeroQ[p],
          B0zero[m1,m2,mu],
          PossibleZeroQ[p - m2] && PossibleZeroQ[m1],
          Delta + Log[mu^2/m2^2] + 2,
          PossibleZeroQ[p - m1] && PossibleZeroQ[m2],
          Delta + Log[mu^2/m1^2] + 2,
          PossibleZeroQ[m1] && PossibleZeroQ[m2],
          Delta - Log[(-p^2 - I eps)/mu^2] + 2,
          True,
          B0analytic[p,m1,m2,mu]
      ];

B0analytic[p_, m1_, m2_, mu_] :=
    Module[{eps, fB, s, xp, xm},
           s = p^2 - m2^2 + m1^2;
           xp = (s + Sqrt[s^2 - 4 p^2 (m1^2 - I eps)]) / (2 p^2);
           xm = (s - Sqrt[s^2 - 4 p^2 (m1^2 - I eps)]) / (2 p^2);
           fB[x_] := Log[1-x] - x Log[1 - 1/x] - 1;
           Limit[Delta - Log[p^2/mu^2] - fB[xp] - fB[xm], eps -> 0,
                 Direction -> -1,
                 Assumptions :> p > 0 && m1 >= 0 && m2 >= 0 && mu > 0]
          ];

B0integral[p_, m1_, m2_, mu_] :=
    Module[{eps},
           Limit[
               Delta - Integrate[Log[((1-x) m1^2 + x m2^2 - x (1-x) p^2 - I eps)/mu^2], {x,0,1}],
               eps -> 0]
          ];

DivB0[_, _, _, _] := Delta;

LogB0[p_, _, _, mu_] := Log[mu^2/p^2];

(* Eq. (B.9) *)
B1impl[p_, m1_, m2_, mu_] :=
    If[PossibleZeroQ[p],
       B1zero[m1,m2,mu],
       $BPMZSign 1/(2 p^2) (A0impl[m2,mu] - A0impl[m1,mu]
                            + (p^2 + m1^2 - m2^2) B0impl[p,m1,m2,mu])
      ];

B1zero[m1_, m2_, mu_] :=
    Which[PossibleZeroQ[m1],
          $BPMZSign (1/4 + Delta/2 + Log[mu^2/m2^2]/2),
          PossibleZeroQ[m2],
          $BPMZSign (3/4 + Delta/2 + Log[mu^2/m1^2]/2),
          PossibleZeroQ[m1 - m2],
          $BPMZSign (Delta + Log[mu^2/m2^2])/2,
          True,
          $BPMZSign 1/2 (Delta + 1 + Log[mu^2/m2^2]
                         + (m1^2/(m1^2 - m2^2))^2 Log[m2^2/m1^2]
                         + 1/2 (m1^2 + m2^2)/(m1^2 - m2^2))
         ];

DivB1[_, _, _, _] := $BPMZSign Delta / 2;

LogB1[p_, m1_, m2_, mu_] :=
    If[PossibleZeroQ[p],
       $BPMZSign Log[mu^2/m2^2]/2,
       $BPMZSign 1/(2 p^2) (LogA0[m2,mu] - LogA0[m1,mu]
                            + (p^2 + m1^2 - m2^2) LogB0[p,m1,m2,mu])
      ];

B11impl[p_, m1_, m2_, mu_] :=
    If[PossibleZeroQ[p],
       B11zero[m1,m2,mu],
       1/(6 p^2) (
           2 A0impl[m2,mu]
           - 2 m1^2 B0impl[p,m1,m2,mu]
           + $BPMZSign 4 (p^2 - m2^2 + m1^2) B1impl[p,m1,m2,mu]
           - m1^2 - m2^2 + p^2/3)
      ];

(* DirectedInfinity[(m1^2 - m2^2)*(-2*m1^2 + 2*m2^2)] *)
B11zero[m1_, m2_, mu_] := Undefined;

DivB11[p_, m1_, m2_, _] := Delta / 3;

LogB11[p_, m1_, m2_, mu_] :=
    If[PossibleZeroQ[p],
       Undefined,
       1/(6 p^2) (
           2 LogA0[m2,mu]
           - 2 m1^2 LogB0[p,m1,m2,mu]
           + $BPMZSign 4 (p^2 - m2^2 + m1^2) LogB1[p,m1,m2,mu])
      ];

(* Eq. (B.10), identical to B00[p,m1,m2,mu] *)
B22impl[p_, m1_, m2_, mu_] :=
    If[PossibleZeroQ[p],
       B22zero[m1,m2,mu],
       1/6 (1/2 (A0impl[m1,mu] + A0impl[m2,mu])
            + (m1^2 + m2^2 - p^2/2) B0impl[p,m1,m2,mu]
            + (m2^2 - m1^2)/(2 p^2) (A0impl[m2,mu] - A0impl[m1,mu]
                                     - (m2^2 - m1^2) B0impl[p,m1,m2,mu])
            + m1^2 + m2^2 - p^2/3)
      ];

B22zero[m1_, m2_, mu_] :=
    Which[PossibleZeroQ[m1] && PossibleZeroQ[m2],
          0,
          PossibleZeroQ[m1],
          (m2^2*(5 + 3*Delta + 3*Log[mu^2/m2^2]))/12,
          PossibleZeroQ[m2],
          (m1^2*(5 + 3*Delta + 3*Log[mu^2/m1^2]))/12,
          PossibleZeroQ[m1 - m2],
          (m2^2*(1 + Delta + Log[mu^2/m2^2]))/2,
          True,
          ((5 + 3*Delta)*(m1^4 - m2^4) + m1^2*(3*m1^2 + m2^2)*Log[mu^2/m1^2] -
           m2^2*(m1^2 + 3*m2^2)*Log[mu^2/m2^2])/(12*(m1^2 - m2^2))
      ];

DivB22[p_, m1_, m2_, _] := Delta (3*m1^2 + 3*m2^2 - p^2)/12;

LogB22[p_, m1_, m2_, mu_] :=
    Which[PossibleZeroQ[p] && PossibleZeroQ[m1 - m2],
          m2^2 Log[mu^2/m2^2]/2,
          PossibleZeroQ[p],
          (m1^2*(3*m1^2 + m2^2)*Log[mu^2/m1^2] -
           m2^2*(m1^2 + 3*m2^2)*Log[mu^2/m2^2])/(12*(m1^2 - m2^2)),
          True,
          1/6 (1/2 (LogA0[m1,mu] + LogA0[m2,mu])
               + (m1^2 + m2^2 - p^2/2) LogB0[p,m1,m2,mu]
               + (m2^2 - m1^2)/(2 p^2) (LogA0[m2,mu] - LogA0[m1,mu]
                                        - (m2^2 - m1^2) LogB0[p,m1,m2,mu]))
         ];

(* Eq. (C.19) *)
C0zero[m1_, m2_, m3_] := Which[
    PossibleZeroQ[m1 - m2] && PossibleZeroQ[m1 - m3],
    -1/(2*m3^2),
    PossibleZeroQ[m1 - m2],
    (-m2^2 + m3^2 - m3^2*Log[m3^2/m2^2])/(m2^2 - m3^2)^2,
    PossibleZeroQ[m1 - m3],
    (m2^2 - m3^2 - m2^2*Log[m2^2/m3^2])/(m2^2 - m3^2)^2,
    PossibleZeroQ[m2 - m3],
    (m1^2 - m3^2 + m1^2*Log[m3^2/m1^2])/(m1^2 - m3^2)^2,
    (* general case *)
    True,
    (+ m2^2 / (m1^2 - m2^2) Log[m2^2/m1^2]
     - m3^2 / (m1^2 - m3^2) Log[m3^2/m1^2]) / (m2^2 - m3^2)
];

C0impl[0, 0, m1_, m2_, m3_, mu_] := C0zero[m1,m2,m3];

C0impl[p1_, p2_, m1_, m2_, m3_, mu_] := NotImplemented;

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
   (C0zero[m1, m3, m4] - C0zero[m2, m3, m4])/(m1^2 - m2^2)];

D0impl[0, 0, 0, m1_, m2_, m3_, m4_, mu_] := D0zero[m1, m2, m3, m4];

D0impl[p1_, p2_, p3_, m1_, m2_, m3_, m4_, mu_] := NotImplemented;

(* Eq. (C.22) *)
D27zero[m1_, m2_, m3_, m4_] :=
   (m1^2 C0zero[m1, m3, m4] - m2^2 C0zero[m2, m3, m4])/(4 (m1^2 - m2^2));

D27impl[0, 0, 0, m1_, m2_, m3_, m4_, mu_] := D27zero[m1, m2, m3, m4];

D27impl[p1_, p2_, p3_, m1_, m2_, m3_, m4_, mu_] := NotImplemented;

(* Eq. (B.11) *)
Fimpl[p_, m1_, m2_, mu_] :=
    A0impl[m1,mu] - 2 A0impl[m2,mu] - (2 p^2 + 2 m1^2 - m2^2) B0impl[p,m1,m2,mu];

Fzero[m1_, m2_, mu_] := Fimpl[0,m1,m2,mu];

DivF[p_, m1_, m2_, mu_] := Delta (-m1^2 - m2^2 - 2*p^2);

LogF[p_, m1_, m2_, mu_] :=
    LogA0[m1,mu] - 2 LogA0[m2,mu] - (2 p^2 + 2 m1^2 - m2^2) LogB0[p,m1,m2,mu];

(* Eq. (B.12) *)
Gimpl[p_, m1_, m2_, mu_] :=
    (p^2 - m1^2 - m2^2) B0impl[p,m1,m2,mu] - A0impl[m1,mu] - A0impl[m2,mu];

Gzero[m1_, m2_, mu_] := Gimpl[0,m1,m2,mu];

DivG[p_, m1_, m2_, mu_] := Delta (-2*m1^2 - 2*m2^2 + p^2);

LogG[p_, m1_, m2_, mu_] :=
    (p^2 - m1^2 - m2^2) LogB0[p,m1,m2,mu] - LogA0[m1,mu] - LogA0[m2,mu];

(* Eq. (B.13) *)
Himpl[p_, m1_, m2_, mu_] := 4 B22impl[p,m1,m2,mu] + Gimpl[p,m1,m2,mu];

Hzero[m1_, m2_, mu_] := Himpl[0,m1,m2,mu];

DivH[p_, m1_, m2_, mu_] := Delta (-m1^2 - m2^2 + (2*p^2)/3);

LogH[p_, m1_, m2_, mu_] := 4 LogB22[p,m1,m2,mu] + LogG[p,m1,m2,mu];

(* Eq. (B.14) *)
B22tildeimpl[p_, m1_, m2_, mu_] := B22impl[p,m1,m2,mu] - A0impl[m1,mu]/4 - A0impl[m2,mu]/4;

B22tildezero[m1_, m2_, mu_] := B22tildeimpl[0,m1,m2,mu];

DivB22tilde[p_, m1_, m2_, mu_] := Delta (-p^2/12);

LogB22tilde[p_, m1_, m2_, mu_] := LogB22[p,m1,m2,mu] - LogA0[m1,mu]/4 - LogA0[m2,mu]/4;

End[];

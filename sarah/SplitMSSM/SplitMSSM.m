Off[General::spell]

Model`Name = "SplitMSSM";
Model`NameLaTeX ="SplitMSSM";
Model`Authors = "F.Staub,A.Voigt";
Model`Date = "2015-02-16";

(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Global symmetries *)

Global[[1]] = {Z[2],RParity};
RpM = {-1,-1,1};
RpP = {1,1,-1};

(* Vector Superfields *)

Gauge[[1]]={B,   U[1], hypercharge, g1,False,RpM};
Gauge[[2]]={WB, SU[2], left,        g2,True, RpM};
Gauge[[3]]={G,  SU[3], color,       g3,False,RpM};

(* Matter fields *)

FermionFields[[1]] = {q, 3, {uL, dL},     1/6, 2,  3, RpM};
FermionFields[[2]] = {l, 3, {vL, eL},    -1/2, 2,  1, RpM};
FermionFields[[3]] = {d, 3, conj[dR],     1/3, 1, -3, RpM};
FermionFields[[4]] = {u, 3, conj[uR],    -2/3, 1, -3, RpM};
FermionFields[[5]] = {e, 3, conj[eR],       1, 1,  1, RpM};

FermionFields[[6]]  = {FG, 1, fG          ,   0  , 1, 8, RpP};
FermionFields[[7]]  = {WB, 1, {{fW0/Sqrt[2],fWp},
                               {fWm,-fW0/Sqrt[2]}}, 0 , 3, 1, RpP};
FermionFields[[8]]  = {B , 1, fB          ,   0  , 1, 1, RpP};
FermionFields[[9]]  = {Hd, 1, {FHd0, FHdm},  -1/2, 2, 1, RpP};
FermionFields[[10]] = {Hu, 1, {FHup, FHu0},   1/2, 2, 1, RpP};

ScalarFields[[1]]  = {H, 1, {Hp, H0},     1/2, 2,  1, RpP};

(*----------------------------------------------*)
(*   ROTATIONS                                  *)
(*----------------------------------------------*)

NameOfStates = {GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][LagrangianInput] = {
    {LagHC  , {AddHC->True }},
    {LagNoHC, {AddHC->False}},
    {LagSplit, {AddHC->True}}
};

LagNoHC = mu2 conj[H].H - 1/2 \[Lambda] conj[H].H.conj[H].H;
LagHC = Yd conj[H].d.q + Ye conj[H].e.l + Yu H.u.q;
LagSplit = - MassG/2 FG.FG - MassWB/2 WB.WB - MassB/2 B.B - \[Mu] Hu.Hd \
    - g2u conj[H].WB.Hu - (g1u/Sqrt[2]) conj[H].B.Hu \
    + g2d H.WB.Hd - (g1d/Sqrt[2]) H.B.Hd;

DEFINITION[GaugeES][DiracSpinors] = {
    Bino -> {fB , conj[fB] },
    FWB0 -> {fW0, 0},
    FWBp -> {fWp, conj[fWm]},
    Glu  -> {fG , conj[fG] },
    FH0  -> {FHd0, conj[FHu0]},
    FHC  -> {FHdm, conj[FHup]},
    Fd1  -> {dL, 0},
    Fd2  -> {0, dR},
    Fu1  -> {uL, 0},
    Fu2  -> {0, uR},
    Fe1  -> {eL, 0},
    Fe2  -> {0, eR},
    Fv   -> {vL,0}
};

(* ----- After EWSB ----- *)

(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] = {
  {{VB,VWB[3]},{VP,VZ},ZZ},
  {{VWB[1],VWB[2]},{VWp,conj[VWp]},ZW}
};

(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs] = {
    {H0, {v,1}, {Ah,\[ImaginaryI]/Sqrt[2]}, {hh,1/Sqrt[2]}}
};

(* ---- Mixings ---- *)

DEFINITION[EWSB][MatterSector] = {
    {{{dL}, {conj[dR]}}, {{DL,Vd}, {DR,Ud}}},
    {{{uL}, {conj[uR]}}, {{UL,Vu}, {UR,Uu}}},
    {{{eL}, {conj[eR]}}, {{EL,Ve}, {ER,Ue}}},
    {{fB, fW0, FHd0, FHu0}, {L0, ZN}},
    {{{fWm, FHdm}, {fWp, FHup}}, {{Lm,UM}, {Lp,UP}}}
};

DEFINITION[EWSB][Phases] = {
    {fG, PhaseGlu}
};

(*------------------------------------------------------*)
(* Dirac-Spinors *)
(*------------------------------------------------------*)

DEFINITION[EWSB][DiracSpinors] = {
    Fd -> {DL, conj[DR]},
    Fe -> {EL, conj[ER]},
    Fu -> {UL, conj[UR]},
    Fv -> {vL, 0},
    Chi -> {L0, conj[L0]},
    Cha -> {Lm, conj[Lp]},
    Glu -> {fG, conj[fG]}
};

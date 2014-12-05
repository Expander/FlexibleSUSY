Off[General::spell]

Model`Name = "NSM";
Model`NameLaTeX ="Next-to Standard Model";
Model`Authors = "F.Staub, A.Voigt";
Model`Date = "2014-11-27";

(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Gauge Groups *)

Gauge[[1]] = {B,   U[1], hypercharge, g1, False};
Gauge[[2]] = {WB, SU[2], left,        g2, True};
Gauge[[3]] = {G,  SU[3], color,       g3, False};

(* Matter Fields *)

FermionFields[[1]] = {q, 3, {uL, dL},     1/6, 2,  3};
FermionFields[[2]] = {l, 3, {vL, eL},    -1/2, 2,  1};
FermionFields[[3]] = {d, 3, conj[dR],     1/3, 1, -3};
FermionFields[[4]] = {u, 3, conj[uR],    -2/3, 1, -3};
FermionFields[[5]] = {e, 3, conj[eR],       1, 1,  1};

ScalarFields[[1]]  = {H, 1, {Hp, H0},     1/2, 2,  1};
ScalarFields[[2]]  = {S, 1, S0      ,       0, 1,  1};

RealScalars = {S};

(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)

NameOfStates = {GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][LagrangianInput] = {
    {LagHC  , {AddHC->True }},
    {LagNoHC, {AddHC->False}}
};

LagZ2 = \
    - mH2 conj[H].H \
    - mS2 S.S \
    + \[Lambda]1 conj[H].H.conj[H].H \
    + \[Lambda]2 S.S.S.S \
    + \[Lambda]3 S.S.conj[H].H;

LagZ2v = \
    + \[Lambda]4 S.conj[H].H \
    + \[Lambda]5 S.S.S;

LagNoHC = LagZ2 + LagZ2v;

LagHC = Yd conj[H].d.q + Ye conj[H].e.l + Yu H.u.q;

(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] = {
    {{VB,VWB[3]}, {VP,VZ}, ZZ},
    {{VWB[1],VWB[2]}, {VWp,conj[VWp]}, ZW}
};

(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs] = {
    {H0, {vH, 1/Sqrt[2]}, {Ah, \[ImaginaryI]/Sqrt[2]}, {phiH, 1/Sqrt[2]}},
    {S0, {vS, 1}, {0, 0}, {phiS, 1}}
};

DEFINITION[EWSB][MatterSector] = {
    {{{dL}, {conj[dR]}}, {{DL,Vd}, {DR,Ud}}},
    {{{uL}, {conj[uR]}}, {{UL,Vu}, {UR,Uu}}},
    {{{eL}, {conj[eR]}}, {{EL,Ve}, {ER,Ue}}},
    {{phiH, phiS}      , {hh, ZH}}
};

(*------------------------------------------------------*)
(* Dirac-Spinors *)
(*------------------------------------------------------*)

DEFINITION[EWSB][DiracSpinors] = {
    Fd -> {DL, conj[DR]},
    Fe -> {EL, conj[ER]},
    Fu -> {UL, conj[UR]},
    Fv -> {vL, 0}
};

DEFINITION[EWSB][GaugeES] = {
    Fd1 -> {FdL, 0},
    Fd2 -> {0, FdR},
    Fu1 -> {Fu1, 0},
    Fu2 -> {0, Fu2},
    Fe1 -> {Fe1, 0},
    Fe2 -> {0, Fe2}
};

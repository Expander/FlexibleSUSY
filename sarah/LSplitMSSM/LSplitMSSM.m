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

Off[General::spell]

Model`Name = "LSplitMSSM";
Model`NameLaTeX ="LSplitMSSM";
Model`Authors = "A.Voigt";
Model`Date = "2018-12-04";

(* Implementation of the split-MSSM with extra light 1st and 2nd
   generation sleptons.  The convention of arXiv:1407.4081 is used. *)

(*-------------------------------------------*)
(*   Particle Content                        *)
(*-------------------------------------------*)

(* Global symmetries *)

Global[[1]] = {Z[2],RParity};
RpP = 1;
RpM = -1;

(* Vector Superfields *)

Gauge[[1]]={B,   U[1], hypercharge, g1,False,RpP};
Gauge[[2]]={WB, SU[2], left,        g2,True, RpP};
Gauge[[3]]={G,  SU[3], color,       g3,False,RpP};

(* Matter fields *)

FermionFields[[1]] = {q, 3, {uL, dL},     1/6, 2,  3, RpP};
FermionFields[[2]] = {l, 3, {vL, eL},    -1/2, 2,  1, RpP};
FermionFields[[3]] = {d, 3, conj[dR],     1/3, 1, -3, RpP};
FermionFields[[4]] = {u, 3, conj[uR],    -2/3, 1, -3, RpP};
FermionFields[[5]] = {e, 3, conj[eR],       1, 1,  1, RpP};

FermionFields[[7]]  = {FW, 1, {{fW0/Sqrt[2],fWp},
                               {fWm,-fW0/Sqrt[2]}}, 0 , 3, 1, RpM};
FermionFields[[8]]  = {FB, 1, fB          ,   0  , 1, 1, RpM};
FermionFields[[9]]  = {Hd, 1, {FHd0, FHdm},  -1/2, 2, 1, RpM};
FermionFields[[10]] = {Hu, 1, {FHup, FHu0},   1/2, 2, 1, RpM};

ScalarFields[[1]]  = {H  , 1, {Hp, H0}  ,   1/2, 2,  1, RpP};
ScalarFields[[2]]  = {SEL, 2, {SvL, SeL},  -1/2, 2,  1, RpM};
ScalarFields[[3]]  = {SER, 2, conj[SeR] ,     1, 1,  1, RpM};

(*----------------------------------------------*)
(*   ROTATIONS                                  *)
(*----------------------------------------------*)

NameOfStates = {GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][LagrangianInput] = {
    {LagHC    , {AddHC->True }},
    {LagNoHC  , {AddHC->False}},
    {LagSplit , {AddHC->True }},
    {LagReSlep, {AddHC->False}},
    {LagSlep  , {AddHC->True }}
};

LagNoHC = mu2 conj[H].H - 1/2 \[Lambda] conj[H].H.conj[H].H;

LagHC = - Yd conj[H].d.q - Ye conj[H].e.l + Yu H.u.q;

(* arXiv:1407.4081, Eq. (4), w/o gluino *)
LagSplit = \
    - MassWB/2 FW.FW - MassB/2 FB.FB - \[Mu] Hu.Hd \
    - g2u conj[H].FW.Hu - (gYu/Sqrt[2]) conj[H].FB.Hu \
    + g2d H.FW.Hd - (gYd/Sqrt[2]) H.FB.Hd;

(* The following terms are in principle also allowed, but omitted in
   this model *)
LagSplit2 = \
    - g2u Hu.conj[FW].conj[H] \
    + g2d FW.conj[H].conj[Hd] \
    - (gYu/Sqrt[2]) Hu.conj[FB].conj[H] \
    - (gYd/Sqrt[2]) FB.conj[H].conj[Hd] ;

(* real part of slepton Lagrangian *)
LagReSlep = (
    - ml2 conj[SEL].SEL - me2 conj[SER].SER
    (* using completeness relation of pauli matrices *)
    + g2llll (2 Delta[lef1,lef4] Delta[lef2,lef3] - Delta[lef1,lef2] Delta[lef3,lef4]) conj[SEL].SEL.SEL.conj[SEL]
    + g2llhh (2 Delta[lef1,lef4] Delta[lef2,lef3] - Delta[lef1,lef2] Delta[lef3,lef4]) conj[H].H.SEL.conj[SEL]
    + g2hhhh (2 Delta[lef1,lef4] Delta[lef2,lef3] - Delta[lef1,lef2] Delta[lef3,lef4]) conj[H].H.H.conj[H]
    + g1llll Delta[lef1,lef1] Delta[lef2,lef2] conj[SEL].SEL.SEL.conj[SEL]
    + g1eeee Delta[lef1,lef1] Delta[lef2,lef2] conj[SER].SER.SER.conj[SER]
    + g1hhhh Delta[lef1,lef1] Delta[lef2,lef2] conj[H].H.H.conj[H]
    + g1llhh Delta[lef1,lef1] Delta[lef2,lef2] conj[SEL].SEL.H.conj[H]
    + g1eehh Delta[lef1,lef1] Delta[lef2,lef2] conj[SER].SER.H.conj[H]
    + g1llee Delta[lef1,lef1] Delta[lef2,lef2] conj[SEL].SEL.SER.conj[SER]
    + gylele Delta[lef1,lef1] Delta[lef2,lef2] conj[SEL].SER.conj[SER].SEL
    + gylehh Delta[lef1,lef1] Delta[lef2,lef2] conj[H].SER.SEL
);

(* part of slepton Lagrangian that needs to be conjugated *)
LagSlep = (
    - Te conj[H].SER.SEL
    - gydsle Hd.SEL.e
    - gydlse Hd.l.SER
    + g2slwl conj[SEL].FW.l
    + g1slbl conj[SEL].FB.l
    + g1sebe conj[SER].FB.e
);

DEFINITION[GaugeES][DiracSpinors] = {
    Bino -> {fB , conj[fB] },
    FW0  -> {fW0, 0},
    FWp  -> {fWp, conj[fWm]},
    FH0  -> {FHd0, conj[FHu0]},
    FHC  -> {FHdm, conj[FHup]},
    Fd1  -> {dL, 0},
    Fd2  -> {0, dR},
    Fu1  -> {uL, 0},
    Fu2  -> {0, uR},
    Fe1  -> {eL, 0},
    Fe2  -> {0, eR},
    Fv   -> {vL, 0}
};

(* ----- After EWSB ----- *)

(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] = {
    {{VB,VWB[3]}, {VP,VZ}, ZZ},
    {{VWB[1],VWB[2]}, {VWp,conj[VWp]}, ZW}
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
    {{{fWm, FHdm}, {fWp, FHup}}, {{Lm,UM}, {Lp,UP}}},
    {{SvL}, {Sv, ZV}},
    {{SeL, SeR}, {Se, ZE}}
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
    Cha -> {Lm, conj[Lp]}
};

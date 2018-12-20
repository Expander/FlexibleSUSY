(* ::Package:: *)

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
(*   Parameters in matrix form               *)
(*-------------------------------------------*)

me2 = DiagonalMatrix[{me21, me22}];
ml2 = DiagonalMatrix[{ml21, ml22}];
g1llll = {{g1llll11, g1llll12}, {g1llll12, g1llll22}};
g1eeee = {{g1eeee11, g1eeee12}, {g1eeee12, g1eeee22}};
g1llhh = DiagonalMatrix[{g1llhh1, g1llhh2}];
g1eehh = DiagonalMatrix[{g1eehh1, g1eehh2}];
g1llee = {{g1llee11, g1llee12}, {g1llee21, g1llee22}};
gylele = {{gylele11, gylele12}, {gylele12, gylele22}};
gyllhh = DiagonalMatrix[{gyllhh1, gyllhh2}];
gyeehh = DiagonalMatrix[{gyeehh1, gyeehh2}];
Te = {{Te11, Te12},{Te21, Te22}};
gyleh = DiagonalMatrix[{gyleh1, gyleh2}];
gydsle = DiagonalMatrix[{gydsle1, gydsle2}];
gydlse = DiagonalMatrix[{gydlse1, gydsle2}];
g2slwl = DiagonalMatrix[{g2slwl1, g2slwl2}];
g1slbl = DiagonalMatrix[{g1slbl1, g1slbl2}];
g1sebe = DiagonalMatrix[{g1sebe1, g1sebe2}];


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
ScalarFields[[2]]  = {SEL1, 1, {SvL1,SeL1}, -1/2, 2,  1, RpM};
ScalarFields[[3]]  = {SEL2, 1, {SvL2,SeL2}, -1/2, 2,  1, RpM};
ScalarFields[[4]]  = {SER1, 1, conj[SeR1],    1, 1,  1, RpM};
ScalarFields[[5]]  = {SER2, 1, conj[SeR2],    1, 1, 1, RpM};

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
    - ml2[[1,1]] conj[SEL1].SEL1 - ml2[[2,2]] conj[SEL2].SEL2 - me2[[1,1]] conj[SER1].SER1 - me2[[2,2]] conj[SER2].SER2
    (* using completeness relation of pauli matrices, assuming diagonal and real Yukawa couplings and mass matrices *)
    + g2llll (2 Delta[lef1,lef4] Delta[lef2,lef3] - Delta[lef1,lef2] Delta[lef3,lef4])conj[SEL1].SEL1.conj[SEL1].SEL1
    + 2 g2llll (2 Delta[lef1,lef4] Delta[lef2,lef3] - Delta[lef1,lef2] Delta[lef3,lef4])conj[SEL1].SEL1.conj[SEL2].SEL2
    + g2llll (2 Delta[lef1,lef4] Delta[lef2,lef3] - Delta[lef1,lef2] Delta[lef3,lef4])conj[SEL2].SEL2.conj[SEL2].SEL2
    + g2llhh (2 Delta[lef1,lef4] Delta[lef2,lef3] - Delta[lef1,lef2] Delta[lef3,lef4]) conj[H].H.SEL1.conj[SEL1] 
    + g2llhh (2 Delta[lef1,lef4] Delta[lef2,lef3] - Delta[lef1,lef2] Delta[lef3,lef4]) conj[H].H.SEL2.conj[SEL2]
    + g2hhhh (2 Delta[lef1,lef4] Delta[lef2,lef3] - Delta[lef1,lef2] Delta[lef3,lef4]) conj[H].H.H.conj[H]
    + g1llll[[1,1]] Delta[lef1,lef2] Delta[lef3,lef4] conj[SEL1].SEL1.conj[SEL1].SEL1
    + 2 g1llll[[1,2]] Delta[lef1,lef2] Delta[lef3,lef4] conj[SEL1].SEL1.conj[SEL2].SEL2 
    + g1llll[[2,2]] Delta[lef1,lef2] Delta[lef3,lef4] conj[SEL2].SEL2.conj[SEL2].SEL2
    + g1eeee[[1,1]] conj[SER1].SER1.SER1.conj[SER1]
    + 2 g1eeee[[1,2]] conj[SER1].SER1.SER2.conj[SER2] 
    + g1eeee[[2,2]] conj[SER2].SER2.SER2.conj[SER2]
    + g1hhhh Delta[lef1,lef2] Delta[lef3,lef4] conj[H].H.H.conj[H]
    + g1llhh[[1,1]] Delta[lef1,lef2] Delta[lef3,lef4] conj[SEL1].SEL1.H.conj[H] 
    + g1llhh[[2,2]] Delta[lef1,lef2] Delta[lef3,lef4] conj[SEL2].SEL2.H.conj[H]
    + g1eehh[[1,1]] Delta[lef3,lef4] conj[SER1].SER1.H.conj[H] 
    + g1eehh[[2,2]] Delta[lef3,lef4] conj[SER2].SER2.H.conj[H]
    + g1llee[[1,1]] Delta[lef1,lef2] conj[SEL1].SEL1.SER1.conj[SER1]
    + g1llee[[1,2]] Delta[lef1,lef2] conj[SEL1].SEL1.SER2.conj[SER2]
    + g1llee[[2,1]] Delta[lef1,lef2] conj[SEL2].SEL2.SER1.conj[SER2] 
    + g1llee[[2,2]] Delta[lef1,lef2] conj[SEL2].SEL2.SER2.conj[SER2]
    + gylele[[1,1]] Delta[lef1,lef4] conj[SEL1].SER1.conj[SER1].SEL1 
    + 2 gylele[[1,2]] Delta[lef1,lef4] conj[SEL1].SER1.conj[SER2].SEL2 
    + gylele[[2,2]] Delta[lef1,lef4] conj[SEL2].SER2.conj[SER2].SEL2
    + gyllhh[[1,1]] Delta[lef1,lef2] Delta[lef3, lef4] conj[H].H.conj[SEL1].SEL1 
    + gyllhh[[2,2]] Delta[lef1,lef2] Delta[lef3, lef4] conj[H].H.conj[SEL2].SEL2
    + gyeehh[[1,1]] conj[H].H.conj[SER1].SER1 
    + gyeehh[[2,2]] conj[H].H.conj[SER2].SER2
);

(* part of slepton Lagrangian that needs to be conjugated *)
LagSlep = (
    - Te[[1,1]] H.conj[SER1].SEL1 - Te[[1,2]] H.conj[SER1].SEL2 - Te[[2,1]] H.conj[SER2].SEL1 - Te[[2,2]] H.conj[SER2].SEL2
    + gyleh[[1,1]] conj[H].conj[SER1].SEL1 + gyleh[[2,2]] conj[H].conj[SER2].SEL2
    - gydsle[[1,1]] Hd.SEL1.e - gydsle[[2,2]] Hd.SEL2.e
    - gydlse[[1,1]] Hd.l.SER1 - gydlse[[2,2]] Hd.l.SER2 
    + g2slwl[[1,1]] conj[SEL1].FW.l + g2slwl[[2,2]] conj[SEL2].FW.l
    + g1slbl[[1,1]] conj[SEL1].FB.l + g1slbl[[2,2]] conj[SEL2].FB.l
    + g1sebe[[1,1]] conj[SER1].FB.e + g1sebe[[2,2]] conj[SER2].FB.e
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
    {H0, {v,1}, {Ah,I/Sqrt[2]}, {hh,1/Sqrt[2]}}
};

(* ---- Mixings ---- *)

DEFINITION[EWSB][MatterSector] = {
    {{{dL}, {conj[dR]}}, {{DL,Vd}, {DR,Ud}}},
    {{{uL}, {conj[uR]}}, {{UL,Vu}, {UR,Uu}}},
    {{{eL}, {conj[eR]}}, {{EL,Ve}, {ER,Ue}}},
    {{fB, fW0, FHd0, FHu0}, {L0, ZN}},
    {{{fWm, FHdm}, {fWp, FHup}}, {{Lm,UM}, {Lp,UP}}},
    {{SvL1, SvL2 }, {Sv, ZV}},
    {{SeL1, SeL2, SeR1, SeR2}, {Se, ZE}}
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

Off[General::spell]
Print["Model file for the E6SSM loaded"];

Model`Name = "E6SSM";
Model`NameLaTeX = "E6SSM";
Model`Authors = "G.Hellwig, P.Diessner";
Model`Date = "2013-21-05";

(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Gauge Superfields *)

Gauge[[1]]={B,        U[1], hypercharge, g1, False};
Gauge[[2]]={WB,      SU[2], left,        g2, True};
Gauge[[3]]={G,       SU[3], color,       g3, False};
Gauge[[4]]={U,        U[1], Ncharge,     gN, False};

(* Chiral Superfields *)

(* Quarks & Leptons *)
Fields[[1]] = {{uL0, dL0},  3, q,  1/6, 2,  3, 1};
Fields[[2]] = {{vL0, eL0},  3, l, -1/2, 2,  1, 2};
Fields[[3]] = {conj[dR0],   3, d,  1/3, 1, -3, 2};
Fields[[4]] = {conj[uR0],   3, u, -2/3, 1, -3, 1};
Fields[[5]] = {conj[eR0],   3, e,    1, 1,  1, 1};
(*Exotic-Quarks*)
Fields[[6]] = {xL,       3, X,    -1/3, 1,  3, -2};
Fields[[7]] = {conj[xR], 3, XBar,  1/3, 1, -3, -3};
(*splitting of the Higgses into 2 Inert-generations and 1 VEV-generation*)
Fields[[8]]  = {{Hd0Inert, HdmInert}, 2, HdInert, -1/2, 2, 1, -3};
Fields[[9]]  = {{HupInert, Hu0Inert}, 2, HuInert,  1/2, 2, 1, -2};
Fields[[10]] = {{Hd0, Hdm},           1, Hd,      -1/2, 2, 1, -3};
Fields[[11]] = {{Hup, Hu0},           1, Hu,       1/2, 2, 1, -2};
(*splitting of the Higgses into 2 Inert-generations and 1 VEV-generation*)
Fields[[12]] = {SRInert, 2, sInert, 0, 1, 1, 5};
Fields[[13]] = {SR,      1, s,      0, 1, 1, 5};
(*Survival-Higgs*)
Fields[[14]] = {{HPrime0,HPrimem},       1, HPrime,   -1/2, 2, 1,  2};
Fields[[15]] = {{HBarPrimep,HBarPrime0}, 1, HBarPrime, 1/2, 2, 1, -2};

NoU1Mixing = True;
AddMixedSofts = False;

(*------------------------------------------------------*)
(* Superpotential *)
(*------------------------------------------------------*)

SuperPotential = { {{1, Yu},{q,Hu,u}}, {{-1,Yd},{q,Hd,d}},
                   {{-1,Ye},{l,Hd,e}},
                   {{1,\[Lambda]3},{Hu,Hd,s}},
                   {{1,\[Kappa]},{X,XBar,s}},
                   {{1,\[Lambda]},{HuInert,HdInert,s}},
                   {{1,muPrime},{HBarPrime,HPrime}}
};

(*-------------------------------------------*)
(* Integrate Out or Delete Particles         *)
(*-------------------------------------------*)

IntegrateOut = {};

(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)

NameOfStates = {GaugeES, EWSB};

(* ----- After EWSB ----- *)

DEFINITION[EWSB][VEVs] = {
    {SHd0,{vd, 1/Sqrt[2]}, {sigmad, \[ImaginaryI]/Sqrt[2]}, {phid, 1/Sqrt[2]}},
    {SHu0,{vu, 1/Sqrt[2]}, {sigmau, \[ImaginaryI]/Sqrt[2]}, {phiu, 1/Sqrt[2]}},
    {SSR, {vS, 1/Sqrt[2]}, {sigmaS, \[ImaginaryI]/Sqrt[2]}, {phiS, 1/Sqrt[2]}}
};

(*Mixings GaugeSector*)
DEFINITION[EWSB][GaugeSector] = {
    {{VB,VWB[3],VU},{VP,VZ,VZprime},ZZ},
    {{VWB[1],VWB[2]},{VWm,conj[VWm]},ZW},
    {{fWB[1],fWB[2],fWB[3]},{fWm,fWp,fW0},ZfW}
};

(* Mixings MatterSector*)
DEFINITION[EWSB][MatterSector] = {
    {{SdL0, SdR0}, {Sd, ZD}},
    {{SuL0, SuR0}, {Su, ZU}},
    {{SeL0, SeR0}, {Se, ZE}},
    {{SvL0}, {Sv, ZV}},
    {{SxL, SxR}, {SX, ZX}},
    {{phid, phiu, phiS}, {hh, ZH}},
    {{sigmad, sigmau,sigmaS}, {Ah, ZA}},
    {{SHdm,conj[SHup]},{Hpm,ZP}},
    {{fB, fW0, FHd0, FHu0,FSR,fU}, {L0, ZN}},
    {{{fWm, FHdm}, {fWp, FHup}}, {{Lm,UM}, {Lp,UP}}},

    {{{FeL0},{conj[FeR0]}},{{FEL,ZEL},{FER,ZER}}},
    {{{FdL0},{conj[FdR0]}},{{FDL,ZDL},{FDR,ZDR}}},
    {{{FuL0},{conj[FuR0]}},{{FUL,ZUL},{FUR,ZUR}}},
    {{{FxL},{conj[FxR]}},{{FXL,ZXL},{FXR,ZXR}}},
    {{SHd0Inert,conj[SHu0Inert]},{SH0Inert,UH0Inert}},
    {{SHdmInert,conj[SHupInert]},{SHpInert,UHpInert}},
    {{FHd0Inert,FHu0Inert},{L0Inert, ZNInert}},
    {{{FHdmInert},{FHupInert}},{{LmInert,UMInert},{LpInert,UPInert}}},
    {{SSRInert},{SS0Inert,ZSSInert}},
    {{SHPrime0,conj[SHBarPrime0]},{SH0Prime,UH0Prime}},
    {{SHPrimem,conj[SHBarPrimep]},{SHpPrime,UHpPrime}},
    {{FHPrime0,FHBarPrime0},{L0Prime,ZNPrime}}
    (* singlet particles should not mix *)
};

DEFINITION[EWSB][Phases] = {
    {fG, PhaseGlu}
};

(*------------------------------------------------------*)
(* Dirac-Spinors *)
(*------------------------------------------------------*)

DEFINITION[GaugeES][DiracSpinors] = {
    Bino -> {fB, conj[fB]},
    Wino -> {fWB, conj[fWB]},
    Glu -> {fG, conj[fG]},
    H01 -> {FHd0, 0},
    H02 -> {0, conj[FHu0]},
    HC1 -> {FHdm, 0},
    HC2 -> {0, conj[FHup]},
    Fd1 -> {FdL0, 0},
    Fd2 -> {0, FdR0},
    Fu1 -> {FuL0, 0},
    Fu2 -> {0, FuR0},
    Fe1 -> {FeL0, 0},
    Fe2 -> {0, FeR0},
    FUprime -> {fU, conj[fU]},
    FS1 -> {FSR, 0},
    FS2 -> {0, conj[FSR]},
    FSInert1 -> {FSRInert, 0},
    FSInert2 -> {0, conj[FSRInert]},
    Fv -> {FvL0, 0},
    FX1 -> {FxL, 0},
    FX2 -> {0, FxR},
    H0I1 -> {FHd0Inert, 0},
    H0I2 -> {0, conj[FHu0Inert]},
    HCI1 -> {FHdmInert, 0},
    HCI2 -> {0, conj[FHupInert]},
    H0P1 -> {FHPrime0, 0},
    H0P2 -> {0, conj[FHBarPrime0]},
    HCP1 -> {FHPrimem, 0},
    HCP2 -> {0, conj[FHBarPrimep]}
};

DEFINITION[EWSB][DiracSpinors] = {
    Chi -> {L0, conj[L0]},
    Cha -> {Lm, conj[Lp]},
    Glu -> {fG, conj[fG]},
    Fd -> {FDL, conj[FDR]},
    Fe -> {FEL, conj[FER]},
    Fu -> {FUL, conj[FUR]},
    Fv -> {FvL0, 0},
    FX -> {FXL, conj[FXR]},
    ChiInert -> {L0Inert, conj[L0Inert]},
    ChaInert -> {LmInert, conj[LpInert]},
    FSInert -> {FSRInert, conj[FSRInert]},
    ChiPrime -> {L0Prime, conj[L0Prime]},
    ChaPrime -> {FHPrimem, conj[FHBarPrimep]}
};

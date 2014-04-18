Off[General::spell]

Model`Name = "EZSSM";
Model`NameLaTeX ="EZSSM";
Model`Authors = "Sophie Underwood";
Model`Date = "2013-11-10";

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
Gauge[[4]]={Bp,  U[1], Ncharge,    g1p,False,RpM};

(* Chiral Superfields *)

SuperFields[[1]] = {q, 3, {uL,  dL},     1/6, 2, 3, 1, RpM};  
SuperFields[[2]] = {l, 3, {vL,  eL},    -1/2, 2, 1, 2, RpM};
SuperFields[[3]] = {Hd,1, {Hd0, Hdm},   -1/2, 2, 1, -3, RpP};
SuperFields[[4]] = {Hu,1, {Hup, Hu0},    1/2, 2, 1, -2, RpP};

SuperFields[[5]] = {d, 3, conj[dR],    1/3, 1, -3, 2, RpM};
SuperFields[[6]] = {u, 3, conj[uR],   -2/3, 1, -3, 1, RpM};
SuperFields[[7]] = {e, 3, conj[eR],      1, 1,  1, 1, RpM};
SuperFields[[8]] = {s, 1, sR,     0, 1,  1, 5, RpP};

SuperFields[[9]]  = {H11I, 1, {H11I0, H11Im},  -1/2, 2, 1, -3, RpP};
SuperFields[[10]] = {H21I, 1, {H21Ip, H21I0},   1/2, 2, 1, -2, RpP};

SuperFields[[11]]  = {H12I, 1, {H12I0, H12Im},  -1/2, 2, 1, -3, RpP};
SuperFields[[12]] = {H22I, 1, {H22Ip, H22I0},   1/2, 2, 1, -2, RpP};

SuperFields[[13]] = {sI, 2, sIR,    0, 1,  1, 5, RpP};
SuperFields[[14]] = {Dx, 3, DxL,  -1/3, 1, 3, -2, RpP};
SuperFields[[15]] = {Dxbar, 3, conj[DxbarR],  1/3, 1, -3, -3, RpP};

SuperFields[[16]] = {Hp, 1, {Hpd0, Hpdm},  -1/2, 2,  1, 2, RpP};
SuperFields[[17]] = {Hpbar, 1, {Hpup, Hpu0}, 1/2, 2,  1, -2, RpP};
NoU1Mixing=True;
AddMixedSofts = False;

(*------------------------------------------------------*)
(*Z2H exact Superpotential *)
(*------------------------------------------------------*)

SuperPotential = Yu u.q.Hu - Yd d.q.Hd - Ye e.l.Hd + \[Lambda] s.Hu.Hd +  \[Lambda]1 s.H21I.H11I +  \[Lambda]2 s.H22I.H12I + \[Kappa] s.Dx.Dxbar + \[Mu]Pr Hpbar.Hp + Xu1 s.Hu.H11I +  Xu2 s.Hu.H12I +  Xd1 s.H21I.Hd + Xd2 s.H22I.Hd ;

(*-------------------------------------------*)
(* Integrate Out or Delete Particles         *)
(*-------------------------------------------*)

IntegrateOut={};
DeleteParticles={};


(*----------------------------------------------*)
(*   ROTATIONS                                  *)
(*----------------------------------------------*)

NameOfStates={GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][DiracSpinors]={
  Bino ->{fB, conj[fB]},
  Wino -> {fWB, conj[fWB]},
  Glu -> {fG, conj[fG]},
  H01 -> {FHd0, 0},
  H02 -> {0, conj[FHu0]},
  HC1 -> {FHdm, 0},
  HC2 -> {0, conj[FHup]},
  Fd1 -> {FdL, 0},
  Fd2 -> {0, FdR},
  Fu1 -> {FuL, 0},
  Fu2 -> {0, FuR},
  Fe1 -> {FeL, 0},
  Fe2 -> {0, FeR},
  Fv -> {FvL,0},
  FS1 -> {FsR,0},
  FS2 -> {0,conj[FsR]},
  FBp -> {fBp,conj[fBp]},
  H01I1 -> {FH11I0, 0},
  H02I1 -> {FH12I0, 0},
  H01I2 -> {0, conj[FH21I0]},
  H02I2 -> {0, conj[FH22I0]}, 
  HC1I1 -> {FH11Im, 0},
  HC2I1 -> {FH12Im, 0},
  HC1I2 -> {0, conj[FH21Ip]},
  HC2I2 -> {0, conj[FH22Ip]},
  FSI1 -> {FsIR, 0},
  FSI2 -> {0, conj[FsIR]},
  FDx1 -> {FDxL, 0},
  FDx2 -> {0, FDxbarR},
  Hp01 -> {FHpd0, 0},
  Hp02 -> {0, conj[FHpu0]},
  HpC1 -> {FHpdm, 0} ,
  HpC2 -> {0, conj[FHpup]} 
};


(* ----- After EWSB ----- *)


(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] =
{ 
   {{VB,VWB[3],VBp},{VP,VZ,VZp},ZZ},
  {{VWB[1],VWB[2]},{VWm,conj[VWm]},ZW},
  {{fWB[1],fWB[2],fWB[3]},{fWm,fWp,fW0},ZfW}
};      
          	

(* ----- VEVs ---- *)

(* DEFINITION[EWSB][VEVs]= 
  { {SHd0, {vd, 1/Sqrt[2]}, {sigmad, \[ImaginaryI]/Sqrt[2]},{phid,1/Sqrt[2]}},
    {SHu0, {vu, 1/Sqrt[2]}, {sigmau, \[ImaginaryI]/Sqrt[2]},{phiu,1/Sqrt[2]}},
    {SsR,  {vs, 1/Sqrt[2]}, {sigmaS, \[ImaginaryI]/Sqrt[2]},{phiS,1/Sqrt[2]}},
    {SH1I0, {vI1, 1/Sqrt[2]}, {sigmaI1, \[ImaginaryI]/Sqrt[2]},{phiI1,1/Sqrt[2]}},
    {SH120,{vI2, 1/Sqrt[2]}, {sigmaI2, \[ImaginaryI]/Sqrt[2]},{phiI2,1/Sqrt[2]}}
   
}; *)

DEFINITION[EWSB][VEVs]= 
  { {SHd0, {vd, 1/Sqrt[2]}, {sigmad, \[ImaginaryI]/Sqrt[2]},{phid,1/Sqrt[2]}},
    {SHu0, {vu, 1/Sqrt[2]}, {sigmau, \[ImaginaryI]/Sqrt[2]},{phiu,1/Sqrt[2]}},
    {SsR,  {vs, 1/Sqrt[2]}, {sigmaS, \[ImaginaryI]/Sqrt[2]},{phiS,1/Sqrt[2]}},
    {SH11I0, {vI11, 1/Sqrt[2]}, {sigmaI11, \[ImaginaryI]/Sqrt[2]},{phiI11,1/Sqrt[2]}},
    {SH12I0, {vI12, 1/Sqrt[2]}, {sigmaI12, \[ImaginaryI]/Sqrt[2]},{phiI12,1/Sqrt[2]}},
    {SH21I0, {vI21, 1/Sqrt[2]}, {sigmaI21, \[ImaginaryI]/Sqrt[2]},{phiI21,1/Sqrt[2]}},
    {SH22I0, {vI22, 1/Sqrt[2]}, {sigmaI22, \[ImaginaryI]/Sqrt[2]},{phiI22,1/Sqrt[2]}}         
};

(*DEFINITION[EWSB][VEVs]= 
  { {SHd0, {vd, 1/Sqrt[2]}, {sigmad, \[ImaginaryI]/Sqrt[2]},{phid,1/Sqrt[2]}},
    {SHu0, {vu, 1/Sqrt[2]}, {sigmau, \[ImaginaryI]/Sqrt[2]},{phiu,1/Sqrt[2]}},
    {SsR,  {vs, 1/Sqrt[2]}, {sigmaS, \[ImaginaryI]/Sqrt[2]},{phiS,1/Sqrt[2]}}
}; *)        


 
(* ---- Mixings ---- *)

DEFINITION[EWSB][MatterSector]= 
{    {{SdL, SdR}, {Sd, ZD}},
     {{SvL}, {Sv, ZV}},
     {{SuL, SuR}, {Su, ZU}},
     {{SeL, SeR}, {Se, ZE}},
     {{SDxL, SDxbarR}, {SDX, ZDX}},
     {{phid, phiu, phiS, phiI11,phiI12, phiI21,phiI22}, {hh, ZH}},
     {{sigmad, sigmau,sigmaS, sigmaI11,sigmaI12,sigmaI21,sigmaI22}, {Ah,ZA}},
     {{SHdm, SH11Im,SH12Im,conj[SHup], conj[SH21Ip], conj[SH22Ip]},{Hpm,ZP}},
     {{fB, fW0, FHd0, FHu0, FsR, fBp, FH11I0, FH21I0,  FH12I0, FH22I0}, {L0, ZN}}, 
     {{{fWm, FHdm, FH11Im, FH12Im}, {fWp, FHup,FH21Ip,FH22Ip}}, {{Lm,UM}, {Lp,UP}}}, 
     {{{FeL},{conj[FeR]}},{{FEL,ZEL},{FER,ZER}}},
     {{{FdL},{conj[FdR]}},{{FDL,ZDL},{FDR,ZDR}}},
     {{{FuL},{conj[FuR]}},{{FUL,ZUL},{FUR,ZUR}}},
     {{{FDxL},{conj[FDxbarR]}},{{FDXL,ZDXL},{FDXR,ZDXR}}}, 
     {{SsIR},{SSI0,ZSSI}},
     {{SHpd0,conj[SHpu0]},{SHp0,UHp0}},
     {{SHpdm,conj[SHpup]},{SHpp,UHpp}},
     {{FHpd0,FHpu0},{L0p,ZNp}}
}; 
       
DEFINITION[EWSB][Phases]= 
   {    {fG, PhaseGlu},
        {ChaP, PhaseChaP}
}; 

	
	
DEFINITION[EWSB][DiracSpinors]={
 Fd -> {FDL, conj[FDR]},
 Fe -> {FEL, conj[FER]},
 Fu -> {FUL, conj[FUR]},
 Fv -> {FvL, 0},
 Chi -> {L0, conj[L0]},
 Cha -> {Lm, conj[Lp]},
 Glu -> {fG, conj[fG]},
 (* ChiI -> {L0I, conj[L0I]},
  ChaI -> {LmI, conj[LpI]},*)
 ChiP -> {L0p, conj[L0p]},
 ChaP -> {FHpdm, conj[FHpup]}, 
 FDX -> {FDXL, conj[FDXR]},
 FSI -> {FsIR, conj[FsIR]}
};	

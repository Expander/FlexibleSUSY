Off[General::spell]

Model`Name = "MSSM";
Model`NameLaTeX ="MSSM";
Model`Authors = "F.Staub";
Model`Date = "2012-09-01";

(* 2013-09-01: changing to new conventions for Superfields, Superpotential and global symmetries *)



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


(* Chiral Superfields *)

SuperFields[[1]] = {q,  3, {uL,  dL},    1/6, 2, 3, RpM};  
SuperFields[[2]] = {l,  3, {vL,  eL},   -1/2, 2, 1, RpM};
SuperFields[[3]] = {Hd, 1, {Hd0, Hdm},  -1/2, 2, 1, RpP};
SuperFields[[4]] = {Hu, 1, {Hup, Hu0},   1/2, 2, 1, RpP};

SuperFields[[5]] = {d, 3, conj[dR],  1/3, 1, -3, RpM};
SuperFields[[6]] = {u, 3, conj[uR], -2/3, 1, -3, RpM};
SuperFields[[7]] = {e, 3, conj[eR],    1, 1,  1, RpM};


(*------------------------------------------------------*)
(* Superpotential *)
(*------------------------------------------------------*)

SuperPotential = Yu u.q.Hu - Yd d.q.Hd - Ye e.l.Hd + \[Mu] Hu.Hd;


(*----------------------------------------------*)
(*   ROTATIONS                                  *)
(*----------------------------------------------*)

NameOfStates={GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][DiracSpinors]={
  Bino ->{fB, conj[fB]},
  Wino -> {fWB, conj[fWB]},
  Glu -> {fG, conj[fG]},
  H0 -> {FHd0, conj[FHu0]},
  HC -> {FHdm, conj[FHup]},
  Fd1 -> {FdL, 0},
  Fd2 -> {0, FdR},
  Fu1 -> {FuL, 0},
  Fu2 -> {0, FuR},
  Fe1 -> {FeL, 0},
  Fe2 -> {0, FeR},
  Fv -> {FvL,0}
};


(* ----- After EWSB ----- *)


(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] =
{ 
  {{VB,VWB[3]},{VP,VZ},ZZ},
  {{VWB[1],VWB[2]},{VWm,conj[VWm]},ZW},
  {{fWB[1],fWB[2],fWB[3]},{fWm,fWp,fW0},ZfW}
};      
          	

(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs]= 
  {{SHd0, {vd, 1/Sqrt[2]}, {sigmad, \[ImaginaryI]/Sqrt[2]},{phid,1/Sqrt[2]}},
   {SHu0, {vu, 1/Sqrt[2]}, {sigmau, \[ImaginaryI]/Sqrt[2]},{phiu,1/Sqrt[2]}}};


 
(* ---- Mixings ---- *)

DEFINITION[EWSB][MatterSector]= 
{    {{SdL, SdR}, {Sd, ZD}},
	 {{SvL}, {Sv, ZV}},
     {{SuL, SuR}, {Su, ZU}},
     {{SeL, SeR}, {Se, ZE}},
     {{phid, phiu}, {hh, ZH}},
     {{sigmad, sigmau}, {Ah, ZA}},
     {{SHdm,conj[SHup]},{Hpm,ZP}},
     {{fB, fW0, FHd0, FHu0}, {L0, ZN}}, 
     {{{fWm, FHdm}, {fWp, FHup}}, {{Lm,UM}, {Lp,UP}}},
     {{{FeL},{conj[FeR]}},{{FEL,ZEL},{FER,ZER}}},
     {{{FdL},{conj[FdR]}},{{FDL,ZDL},{FDR,ZDR}}},
     {{{FuL},{conj[FuR]}},{{FUL,ZUL},{FUR,ZUR}}}       }; 
       
DEFINITION[EWSB][Phases]= 
{    {fG, PhaseGlu}
    }; 

	
	
DEFINITION[EWSB][DiracSpinors]={
 Fd ->{  FDL, conj[FDR]},
 Fe ->{  FEL, conj[FER]},
 Fu ->{  FUL, conj[FUR]},
 Fv ->{  FvL, 0},
 Chi ->{ L0, conj[L0]},
 Cha ->{ Lm, conj[Lp]},
 Glu ->{ fG, conj[fG]}
};	

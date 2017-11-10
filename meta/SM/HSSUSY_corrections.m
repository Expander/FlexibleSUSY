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

(*** High-scale MSSM boundary conditions to the MSSM ***)
(* Taken from arXiv:1407.4081, arxiv:1504.05200, arxiv:1703.08166 *)

(* abbreviations *)
xQU = Sqrt[Abs[msq2[3,3]/msu2[3,3]]];
xQD = Sqrt[Abs[msq2[3,3]/msd2[3,3]]];
xLE = Sqrt[Abs[msl2[3,3]/mse2[3,3]]];

At = AtInput;
xt = At - MuInput/TanBeta;
yt = At + MuInput*TanBeta;
xtt = xt^2/Sqrt[Abs[msq2[3,3] msu2[3,3]]];

Ab = AbInput;
xb = Ab - MuInput*TanBeta;
yb = Ab + MuInput/TanBeta;
xbb = xb^2/Sqrt[Abs[msq2[3,3] msd2[3,3]]];

Atau = AtauInput;
xtau = Atau - MuInput*TanBeta;
xtaut = xtau^2/Sqrt[Abs[msl2[3,3] mse2[3,3]]];
ytau = Atau + MuInput/TanBeta;

(* arXiv:1407.4081, Eq. (5) *)
gYu = Sqrt[3/5] g1 Sin[ArcTan[TanBeta]];
gYd = Sqrt[3/5] g1 Cos[ArcTan[TanBeta]];
g2u = g2 Sin[ArcTan[TanBeta]];
g2d = g2 Cos[ArcTan[TanBeta]];

(* arXiv:1407.4081, Eq. (3) *)
lambdaTree = 1/4 (g2^2 + 3/5 g1^2) Cos[2 ArcTan[TanBeta]]^2;

(* arXiv:1407.4081, Eq. (9) *)
lambda1LReg = 1/(4 Pi)^2 (
    - 9/100 g1^4 - 3/10 g1^2 g2^2
    - (3/4 - Cos[2 ArcTan[TanBeta]]^2/6) * g2^4
    );

(* arXiv:1407.4081, Eq. (10) *)
lambda1LPhi = 1/(4 Pi)^2 (
    3 Yu[3,3]^2 (Yu[3,3]^2 + 1/2 (g2^2-g1^2/5) Cos[2 ArcTan[TanBeta]]) Log[msq2[3,3]/SCALE^2] (1 + DeltaEFT v^2/Abs[msq2[3,3]])
    + 3 Yu[3,3]^2 (Yu[3,3]^2 + 2/5 g1^2 Cos[2 ArcTan[TanBeta]]) Log[msu2[3,3]/SCALE^2] (1 + DeltaEFT v^2/Abs[msu2[3,3]])
    + Cos[2 ArcTan[TanBeta]]^2/300 (
        3 (g1^4 + 25 g2^4) (
            + Log[msq2[1,1]/SCALE^2] (1 + DeltaEFT v^2/Abs[msq2[1,1]])
            + Log[msq2[2,2]/SCALE^2] (1 + DeltaEFT v^2/Abs[msq2[2,2]])
            + Log[msq2[3,3]/SCALE^2] (1 + DeltaEFT v^2/Abs[msq2[3,3]])
        )
        + 24 g1^4 (
            + Log[msu2[1,1]/SCALE^2] (1 + DeltaEFT v^2/Abs[msu2[1,1]])
            + Log[msu2[2,2]/SCALE^2] (1 + DeltaEFT v^2/Abs[msu2[2,2]])
            + Log[msu2[3,3]/SCALE^2] (1 + DeltaEFT v^2/Abs[msu2[3,3]])
        )
        + 6 g1^4 (
            + Log[msd2[1,1]/SCALE^2] (1 + DeltaEFT v^2/Abs[msd2[1,1]])
            + Log[msd2[2,2]/SCALE^2] (1 + DeltaEFT v^2/Abs[msd2[2,2]])
            + Log[msd2[3,3]/SCALE^2] (1 + DeltaEFT v^2/Abs[msd2[3,3]])
        )
        + (9 g1^4 + 25 g2^4) (
            + Log[msl2[1,1]/SCALE^2] (1 + DeltaEFT v^2/Abs[msl2[1,1]])
            + Log[msl2[2,2]/SCALE^2] (1 + DeltaEFT v^2/Abs[msl2[2,2]])
            + Log[msl2[3,3]/SCALE^2] (1 + DeltaEFT v^2/Abs[msl2[3,3]])
        )
        + 18 g1^4 (
            + Log[mse2[1,1]/SCALE^2] (1 + DeltaEFT v^2/Abs[mse2[1,1]])
            + Log[mse2[2,2]/SCALE^2] (1 + DeltaEFT v^2/Abs[mse2[2,2]])
            + Log[mse2[3,3]/SCALE^2] (1 + DeltaEFT v^2/Abs[mse2[3,3]])
        )
    )
    + (
       + 1/4800 (261 g1^4 + 630 g1^2 g2^2 + 1325 g2^4
                 -4 Cos[4 ArcTan[TanBeta]] (9 g1^4 + 90 g1^2 g2^2 + 175 g2^4)
                 -9 Cos[8 ArcTan[TanBeta]] (3 g1^2 + 5 g2^2)^2) Log[mAInput^2/SCALE^2]
       - 3/16 (3/5 g1^2 + g2^2)^2 Sin[4 ArcTan[TanBeta]]^2
    ) (1 + DeltaEFT v^2/mAInput^2)
    + (
       + 6 Yu[3,3]^4 xtt (TCF[1][xQU] - xtt/12 TCF[2][xQU])
       + 3/4 Yu[3,3]^2 xtt Cos[2 ArcTan[TanBeta]] (3/5 g1^2 TCF[3][xQU] + g2^2 TCF[4][xQU])
       - 1/4 Yu[3,3]^2 xtt Cos[2 ArcTan[TanBeta]]^2 (3/5 g1^2 + g2^2) TCF[5][xQU]
    ) (1 + DeltaEFT v^2/Min[Abs[msq2[3,3]],Abs[msu2[3,3]]])
    );

(* arXiv:1407.4081, Eq. (11) *)
lambda1LChi1 = With[
    { r1 = M1Input / MuInput, r2 = M2Input / MuInput },
    1/(4 Pi)^2 (
        1/2 betaLambda Log[MuInput^2/SCALE^2]                       (1 + DeltaEFT v^2/MuInput^2)
        - 7/12 TCf[1][r1] (gYd^4 + gYu^4)                           (1 + DeltaEFT v^2/Min[M1Input^2,MuInput^2])
        - 9/4 TCf[2][r2] (g2d^4 + g2u^4)                            (1 + DeltaEFT v^2/Min[M2Input^2,MuInput^2])
        - 3/2 TCf[3][r1] gYd^2 gYu^2                                (1 + DeltaEFT v^2/Min[M1Input^2,MuInput^2])
        - 7/2 TCf[4][r2] g2d^2 g2u^2                                (1 + DeltaEFT v^2/Min[M2Input^2,MuInput^2])
        - 8/3 TCf[5][r1,r2] gYd gYu g2d g2u                         (1 + DeltaEFT v^2/Min[M1Input^2,M2Input^2,MuInput^2])
        - 7/6 TCf[6][r1,r2] (gYd^2 g2d^2 + gYu^2 g2u^2)             (1 + DeltaEFT v^2/Min[M1Input^2,M2Input^2,MuInput^2])
        - 1/6 TCf[7][r1,r2] (gYd^2 g2u^2 + gYu^2 g2d^2)             (1 + DeltaEFT v^2/Min[M1Input^2,M2Input^2,MuInput^2])
        - 4/3 TCf[8][r1,r2] (gYd g2u + gYu g2d) (gYd g2d + gYu g2u) (1 + DeltaEFT v^2/Min[M1Input^2,M2Input^2,MuInput^2])
        + 2/3 TCf0[r1] gYd gYu (lambdaTree - 2 (gYd^2 + gYu^2))     (1 + DeltaEFT v^2/Min[M1Input^2,MuInput^2])
        + 2 TCf0[r2] g2d g2u (lambdaTree - 2 (g2d^2 + g2u^2))       (1 + DeltaEFT v^2/Min[M2Input^2,MuInput^2])
        + 1/3 TCg0[r1] lambdaTree (gYd^2 + gYu^2)                   (1 + DeltaEFT v^2/Min[M1Input^2,MuInput^2])
        + TCg0[r2] lambdaTree (g2d^2 + g2u^2)                       (1 + DeltaEFT v^2/Min[M2Input^2,MuInput^2])
    )
];

(* arXiv:1407.4081, Eq. (13) *)
lambda1LChi2 = 1/(4 Pi)^2 (
    -1/6 Cos[2 ArcTan[TanBeta]]^2 (
        + 2 g2^4 Log[M2Input^2/SCALE^2]             (1 + DeltaEFT v^2/M2Input^2)
        + (9/25 g1^4 + g2^4) Log[MuInput^2/SCALE^2] (1 + DeltaEFT v^2/MuInput^2)
    )
);

(* complete lambda 1-loop threshold correction O(gb^4); as discussed
   in arxiv:1703.08166.

   We express the correction in terms of the bottom yukawa of the MSSM
   to resum tanb enhanced corrections
*)
lambda1Lbottom = With[{
    Nc = 3,
    mQ3  = Sqrt[msq2[3,3]],
    mQ32 = msq2[3,3],
    mU3 = Sqrt[msu2[3,3]],
    mD3 = Sqrt[msd2[3,3]],
    mD32 = msd2[3,3],
    M3 = M3Input,
    Mu = MuInput,
    Xt = xt,
    gb = Yd[3,3], (* SM Yukawa coupling *)
    gt = Yu[3,3], (* SM Yukawa coupling *)
    Q = SCALE,
    Q2 = SCALE^2,
    MA = mAInput,
    cosb = Cos[ArcTan[TanBeta]],
    cos2beta = Cos[2*ArcTan[TanBeta]],
    sinb = Sin[ArcTan[TanBeta]],
    Xb = xb,
    Xtildet = xtt,
    Xtildeb = xbb,
    k = 1/(4 Pi)^2,
    CF = 4/3
    },
      ybMSSM[mQ3_,mU3_,mD3_,M3_,Mu_,TanBeta_,Xt_,Xb_] := Module[{deltagsb, deltagbyL1, deltagbyL2, deltagbyL3, deltagbyL4},
        deltagsb   = - g3^2*CF*k*(1+Log[M3^2/Q2]+TCF[6][mQ3/M3]+TCF[6][mD3/M3]-Xb/M3*TCF[9][mQ3/M3,mD3/M3]);
        deltagbyL1 = - gb^2/cosb^2*k*(3/4*Log[Mu^2/Q2]+3/8*sinb^2*(2*Log[MA^2/Q2]-1)+TCF[6][mQ3/Mu]+1/2*TCF[6][mD3/Mu]);
        deltagbyL2 = - gt^2/sinb^2*k*(1/4*Log[Mu^2/Q2]+1/8*cosb^2*(2*Log[MA^2/Q2]-1)+sinb^2*(Log[MA^2/Q2]-1));
        deltagbyL3 = - gt^2/sinb^2*k*(1/2*TCF[6][mU3/Mu]+(Xt*TanBeta)/(2*Mu)*TCF[9][mQ3/Mu,mU3/Mu] );
        deltagbyL4 = - 1/2*(-gt^2*Nc*k*Xtildet/6*TCF[5][xQU]-gb^2*Nc*k*Xtildeb/6*TCF[5][xQD]);
        (gb/(1-deltagsb-(deltagbyL1+deltagbyL2+deltagbyL3+deltagbyL4)))
       ];
      k*(6 Xtildeb ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^4 (TCF[1][xQD] - 1/12 Xtildeb TCF[2][xQD])
      + 3/4 Xtildeb ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2 cos2beta (-(3/10) g1^2 TCF[3][xQD] + (-((3 g1^2)/10) - g2^2) TCF[4][xQD])
      - 1/4 ((3 g1^2)/5 + g2^2) Xtildeb ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2 cos2beta^2 TCF[5][xQD]
      + 3*ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2*(ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2 - 1/2*(g2^2+g1^2/5)*cos2beta)*Log[mQ32/Q2]
      + 3 ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2 (ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2 - 1/5 g1^2 cos2beta) Log[mD32/Q2]
      ) (1 + DeltaEFT v^2/Min[Abs[msq2[3,3]], Abs[msu2[3,3]], Abs[msd2[3,3]], M3Input^2, MuInput^2])
];

(* lambda 1-loop threshold correction O(alpha_tau),
   taken from SUSYHD 1.0.2, arxiv:1504.05200 *)
lambda1Ltau = With[{
    ML32 = msl2[3,3],
    ME32 = mse2[3,3],
    Q2 = SCALE^2,
    gtau = Ye[3,3], (* SM Yukawa coupling *)
    XtildeTau = xtaut,
    cos2beta = Cos[2*ArcTan[TanBeta]],
    k = 1/(4*Pi)^2
    },
     k*(2*XtildeTau*gtau^4*(TCF[1][xLE] - (XtildeTau*TCF[2][xLE])/12)
     + (XtildeTau*gtau^2*cos2beta*((-9*g1^2*TCF[3][xLE])/10 + ((3*g1^2)/10 - g2^2)*TCF[4][xLE]))/4
     -  (((3*g1^2)/5 + g2^2)*XtildeTau*gtau^2*cos2beta^2*TCF[5][xLE])/12
     + gtau^2*(gtau^2 - (3*g1^2*cos2beta)/5)*Log[ME32/Q2] + gtau^2*(gtau^2 - 1/2*(g2^2 - 3/5*g1^2)*cos2beta)*Log[ML32/Q2]
     ) (1 + DeltaEFT v^2/Min[Abs[msl2[3,3]], Abs[mse2[3,3]], MuInput^2])
    ];

(* lambda 2-loop threshold correction O(alpha_s alpha_t^2),
   arXiv:1407.4081, Eq. (36)
   Valid in the limit M3 = MQ[3,3] = MU[3,3] = MSUSY *)
lambda2LPhiHSSAlphaTAlphaSDegenerate = With[{ r = xt / SCALE },
    (g3^2 Yu[3,3]^4)/(96 Pi^4) (
        -12 r - 6 r^2 + 14 r^3 + 1/2 r^4 - r^5
    )
];

(* lambda 2-loop threshold correction O(alpha_s alpha_t^2),
   arXiv:1407.4081, provided by Pietro Slavich,
   in the limit M3 = MQ[3,3] = MU[3,3] != MSUSY *)
lambda2LPhiHSSAlphaTAlphaSMQMUM3Degenerate = With[{
    htSM = Yu[3,3],
    CF = 4/3,
    Nc = 3,
    as = g3^2 / (4 Pi),
    g = M3Input^2,
    x1 = msq2[3,3]/M3Input^2,
    x2 = msu2[3,3]/M3Input^2,
    q2 = SCALE^2,
    Xt = xt / M3Input
    },
    (
      -(as*CF*htSM^4*Nc*(Xt*(24 + 12*Xt - 28*Xt^2 - Xt^3 + 2*Xt^4) -
       2*Xt*(24 - 24*Xt - 4*Xt^2 + Xt^3)*Log[g/q2] +
       36*Log[g/q2]^2))/(192*Pi^3)
    )
];

(* lambda 2-loop threshold correction O(alpha_s alpha_t^2),
   arXiv:1407.4081, provided by Pietro Slavich,
    General expression. *)
lambda2LPhiHSSAlphaTAlphaSFull = With[{
    htSM = Yu[3,3],
    CF = 4/3,
    Nc = 3,
    as = g3^2 / (4 Pi),
    g = (M3Input)^2,
    x1 = msq2[3,3]*(1-0.005)/(M3Input*(1+0.001))^2,
    x2 = msu2[3,3]/(M3Input*(1-0.005))^2,
    q2 = SCALE^2,
    Xt = xt / M3Input
    },    
    htSM^4*as*CF*Nc/(4*Pi)^3 (
       -12*Log[g/q2]^2 + ((-6 + 8*x1 - 4*x1^2)*Log[x1]^2)/(-1 + x1)^2 +
        Log[g/q2]*(8 - 4/x1 - 4/x2 - 12*Log[x1] - 12*Log[x2]) +
        (2*(-7 + 7*x2 - 3*x2^2 + x1*(6 - 5*x2 + 2*x2^2))*Log[x2])/
         ((-1 + x1)*(-1 + x2)^2) + ((-6 + 8*x2 - 4*x2^2)*Log[x2]^2)/(-1 + x2)^2 +
        Log[x1]*((2*(-7 + x1*(7 - 5*x2) + 6*x2 + x1^2*(-3 + 2*x2)))/
           ((-1 + x1)^2*(-1 + x2)) - (2*((-2 + x2)*x2 + x1*(-2 + 8*x2 - 4*x2^2) +
             x1^2*(1 - 4*x2 + 2*x2^2))*Log[x2])/((-1 + x1)^2*(-1 + x2)^2)) +
        Xt^2*(-8/(x1*x2) - (4*(4*x1^3 - x2 - 2*x1^2*(4 + x2) + x1*(3 + 4*x2))*
            Log[x1]^2)/((-1 + x1)^2*(x1 - x2)^2) -
          (4*(7 + 3*x1*(-2 + x2) - 4*x2)*Log[x2])/((-1 + x1)*(x1 - x2)*(-1 + x2)) +
          (4*(x2*(-3 + 8*x2 - 4*x2^2) + x1*(1 - 4*x2 + 2*x2^2))*Log[x2]^2)/
           ((x1 - x2)^2*(-1 + x2)^2) + Log[g/q2]*(8/(x1*x2) -
            (24*Log[x1])/(x1 - x2) + (24*Log[x2])/(x1 - x2)) +
          Log[x1]*((4*(7 - 6*x2 + x1*(-4 + 3*x2)))/((-1 + x1)*(x1 - x2)*(-1 + x2)) +
            (4*(((-2 + x1)*x1*(x1 - x2))/(-1 + x1)^2 + ((-2 + x2)*x2*(-x1 + x2))/
                (-1 + x2)^2 + 2*(x1 + x2))*Log[x2])/(x1 - x2)^2)) +
        Xt^5*((-8*x1*(x1 + x2)*Log[x1]^2)/((-1 + x1)*(x1 - x2)^4) +
          (16*x2*Log[x2])/((-1 + x2)*(-x1 + x2)^3) - (8*x2*(x1 + x2)*Log[x2]^2)/
           ((x1 - x2)^4*(-1 + x2)) + Log[x1]*((16*x1)/((-1 + x1)*(x1 - x2)^3) +
            (8*(x1 + x2)*(-x2 + x1*(-1 + 2*x2))*Log[x2])/((-1 + x1)*(x1 - x2)^4*
              (-1 + x2)))) + Xt*((8*(-2 + x1)*Log[x1]^2)/((-1 + x1)*(x1 - x2)) +
          (16*Log[x2])/(x1 - x2) - (8*(-2 + x2)*Log[x2]^2)/((x1 - x2)*(-1 + x2)) +
          Log[g/q2]*((16*Log[x1])/(x1 - x2) - (16*Log[x2])/(x1 - x2)) +
          Log[x1]*(-16/(x1 - x2) - (8*Log[x2])/((-1 + x1)*(-1 + x2))) +
          (32*(-((-1 + x2)*PolyLog[2, (-1 + x1)/x1]) +
             (-1 + x1)*PolyLog[2, (-1 + x2)/x2]))/((-1 + x1)*(x1 - x2)*(-1 + x2))) -
        (2*(4*x1*(-1 + x2)^2*x2*PolyLog[2, (-1 + x1)/x1] +
           (-1 + x1)*((-1 + x2)*(2*(-1 + x2)*x2 + x1*(-2 + 9*x2 - 6*x2^2) +
               x1^2*(2 - 6*x2 + 3*x2^2)) + 4*(-1 + x1)*x1*x2*
              PolyLog[2, (-1 + x2)/x2])))/((-1 + x1)^2*x1*(-1 + x2)^2*x2) +
        Xt^3*((8*(-2 + x1^2 - 3*x1*(-1 + x2) + x2)*Log[x1]^2)/
           ((-1 + x1)*(x1 - x2)^3) - (16*(x1 + 3*x2)*Log[x2])/(x1 - x2)^3 +
          (8*(-2 + x1 + 3*x2 - 3*x1*x2 + x2^2)*Log[x2]^2)/((-1 + x2)*(-x1 + x2)^3) +
          Log[g/q2]*(32/(x1 - x2)^2 - (16*(x1 + x2)*Log[x1])/(x1 - x2)^3 +
            (16*(x1 + x2)*Log[x2])/(x1 - x2)^3) +
          Log[x1]*((16*(3*x1 + x2))/(x1 - x2)^3 + (16*(x1 + x2 - 2*x1*x2)*Log[x2])/
             ((-1 + x1)*(x1 - x2)^2*(-1 + x2))) -
          (16*(4*(x1 - x2) + (-2 + x1 + x2)*PolyLog[2, (-1 + x1)/x1] -
             (-2 + x1 + x2)*PolyLog[2, (-1 + x2)/x2]))/(x1 - x2)^3) +
        Xt^4*((2*(5*x1^4 + 4*x1^3*(-2 + x2) - 2*x2 - x1^2*x2*(10 + x2) +
             2*x1*(1 + 4*x2 + x2^2))*Log[x1]^2)/((-1 + x1)^2*(x1 - x2)^4) +
          (2*(-4 - 7*x2 + 29*x2^2 - 16*x2^3 + x1^2*(6 - 7*x2 + 3*x2^2) +
             x1*(-3 + 15*x2 - 31*x2^2 + 15*x2^3))*Log[x2])/
           ((-1 + x1)*(x1 - x2)^3*(-1 + x2)^2) -
          (2*(-2*x2 + x1^2*(-2 + x2)*x2 + 8*x2^3 - 5*x2^4 -
             2*x1*(-1 + x2)^2*(-1 + 2*x2))*Log[x2]^2)/((x1 - x2)^4*(-1 + x2)^2) +
          Log[g/q2]*((-4*(x1 + x2 + 6*x1*x2))/(x1*(x1 - x2)^2*x2) +
            (4*(2 + 3*x1 + 3*x2)*Log[x1])/(x1 - x2)^3 -
            (4*(2 + 3*x1 + 3*x2)*Log[x2])/(x1 - x2)^3) +
          Log[x1]*((2*(4 + x1^3*(16 - 15*x2) + 3*x2 - 6*x2^2 +
               x1^2*(-29 + 31*x2 - 3*x2^2) + x1*(7 - 15*x2 + 7*x2^2)))/
             ((-1 + x1)^2*(x1 - x2)^3*(-1 + x2)) +
            (2*(-4*x1^2 - 8*x1*x2 - 4*x2^2 - ((-2 + x1)*x1*(x1 - x2)*(x1 + x2))/
                (-1 + x1)^2 - ((-2 + x2)*x2*(-x1 + x2)*(x1 + x2))/(-1 + x2)^2)*
              Log[x2])/(x1 - x2)^4) -
          (4*(-x1^2 + x1^3 - 6*x1^2*x2 + 5*x1^3*x2 + x2^2 + 6*x1*x2^2 -
             5*x1^3*x2^2 - x2^3 - 5*x1*x2^3 + 5*x1^2*x2^3 +
             (-1 + x1)*x1*(-1 + x2)*x2*(-2 + x1 + x2)*PolyLog[2, (-1 + x1)/x1] -
             (-1 + x1)*x1*(-1 + x2)*x2*(-2 + x1 + x2)*PolyLog[2, (-1 + x2)/x2]))/
           ((-1 + x1)*x1*(x1 - x2)^3*(-1 + x2)*x2))
    )
];

(* lambda 2-loop threshold correction O(alpha_s alpha_t^2),
   in the limit MQ[3,3]/M3 = MU[3,3]/M3
 *)
lambda2LPhiHSSAlphaTAlphaSX1X2Degenerate = With[{
    htSM = Yu[3,3],
    CF = 4/3,
    Nc = 3,
    as = g3^2 / (4 Pi),
    g = M3Input^2,
    x2 = msu2[3,3]/M3Input^2,
    q2 = SCALE^2,
    Xt = xt / M3Input
    },
    (as*CF*htSM^4*Nc*(24*x2^2 - 78*x2^3 + 72*x2^4 - 18*x2^5 -
   48*x2^2*Xt + 96*x2^3*Xt - 48*x2^4*Xt - 24*x2*Xt^2 +
   108*x2^2*Xt^2 - 96*x2^3*Xt^2 + 12*x2^4*Xt^2 + 8*x2*Xt^3 -
   72*x2^2*Xt^3 + 64*x2^3*Xt^3 + 4*Xt^4 - 13*x2*Xt^4 +
   10*x2^2*Xt^4 - x2^3*Xt^4 + 4*x2*Xt^5 - 4*x2^2*Xt^5 -
   24*x2^2*Log[g] + 72*x2^3*Log[g] - 72*x2^4*Log[g] +
   24*x2^5*Log[g] + 48*x2^2*Xt*Log[g] - 96*x2^3*Xt*Log[g] +
   48*x2^4*Xt*Log[g] + 24*x2*Xt^2*Log[g] - 120*x2^2*Xt^2*Log[g] +
   168*x2^3*Xt^2*Log[g] - 72*x2^4*Xt^2*Log[g] - 8*x2*Xt^3*Log[g] +
   16*x2^2*Xt^3*Log[g] - 8*x2^3*Xt^3*Log[g] - 4*Xt^4*Log[g] +
   14*x2*Xt^4*Log[g] - 16*x2^2*Xt^4*Log[g] + 6*x2^3*Xt^4*Log[g] -
   36*x2^3*Log[g]^2 + 72*x2^4*Log[g]^2 - 36*x2^5*Log[g]^2 +
   24*x2^2*Log[q2] - 72*x2^3*Log[q2] + 72*x2^4*Log[q2] -
   24*x2^5*Log[q2] - 48*x2^2*Xt*Log[q2] + 96*x2^3*Xt*Log[q2] -
   48*x2^4*Xt*Log[q2] - 24*x2*Xt^2*Log[q2] +
   120*x2^2*Xt^2*Log[q2] - 168*x2^3*Xt^2*Log[q2] +
   72*x2^4*Xt^2*Log[q2] + 8*x2*Xt^3*Log[q2] -
   16*x2^2*Xt^3*Log[q2] + 8*x2^3*Xt^3*Log[q2] + 4*Xt^4*Log[q2] -
   14*x2*Xt^4*Log[q2] + 16*x2^2*Xt^4*Log[q2] -
   6*x2^3*Xt^4*Log[q2] + 72*x2^3*Log[g]*Log[q2] -
   144*x2^4*Log[g]*Log[q2] + 72*x2^5*Log[g]*Log[q2] -
   36*x2^3*Log[q2]^2 + 72*x2^4*Log[q2]^2 - 36*x2^5*Log[q2]^2 +
   96*x2^2*Xt*Log[1 - (-1 + x2)/x2] +
   32*Xt^3*Log[1 - (-1 + x2)/x2] -
   72*x2*Xt^3*Log[1 - (-1 + x2)/x2] +
   48*x2^2*Xt^3*Log[1 - (-1 + x2)/x2] +
   8*Xt^4*Log[1 - (-1 + x2)/x2] -
   18*x2*Xt^4*Log[1 - (-1 + x2)/x2] +
   12*x2^2*Xt^4*Log[1 - (-1 + x2)/x2] + 84*x2^3*Log[x2] -
   72*x2^4*Log[x2] + 24*x2^5*Log[x2] + 96*x2^2*Xt*Log[x2] -
   144*x2^3*Xt*Log[x2] + 48*x2^4*Xt*Log[x2] -
   48*x2^2*Xt^2*Log[x2] + 168*x2^3*Xt^2*Log[x2] -
   72*x2^4*Xt^2*Log[x2] + 32*Xt^3*Log[x2] - 72*x2*Xt^3*Log[x2] -
   8*x2^3*Xt^3*Log[x2] + 8*Xt^4*Log[x2] - 12*x2*Xt^4*Log[x2] -
   4*x2^2*Xt^4*Log[x2] + 6*x2^3*Xt^4*Log[x2] + 4*x2*Xt^5*Log[x2] -
   72*x2^3*Log[g]*Log[x2] + 144*x2^4*Log[g]*Log[x2] -
   72*x2^5*Log[g]*Log[x2] + 72*x2^3*Log[q2]*Log[x2] -
   144*x2^4*Log[q2]*Log[x2] + 72*x2^5*Log[q2]*Log[x2] -
   36*x2^3*Log[x2]^2 + 72*x2^4*Log[x2]^2 - 36*x2^5*Log[x2]^2 -
   48*x2^3*PolyLog[2, (-1 + x2)/x2] +
   96*x2^3*Xt*PolyLog[2, (-1 + x2)/x2]))/(192*Pi^3*(-1 + x2)^2*x2^3)
];

(* lambda 2-loop threshold correction O(alpha_s alpha_t^2),
   in the limit MQ[3,3] = M3 (x1 -> 1), MU[3,3] != M3 (x2 != 1)
 *)
lambda2LPhiHSSAlphaTAlphaSMQM3Degenerate = With[{
    htSM = Yu[3,3],
    CF = 4/3,
    Nc = 3,
    as = g3^2 / (4 Pi),
    g = M3Input^2,
    x1 = msq2[3,3]/M3Input^2,
    x2 = msu2[3,3]/M3Input^2,
    q2 = SCALE^2,
    Xt = xt / M3Input
    },
    (
       -(as*CF*htSM^4*Nc*(4 - 25*x2 + 63*x2^2 - 82*x2^3 + 58*x2^4 -
       21*x2^5 + 3*x2^6 - 32*x2*Xt + 128*x2^2*Xt - 192*x2^3*Xt +
       128*x2^4*Xt - 32*x2^5*Xt - 8*Xt^2 + 32*x2*Xt^2 -
       48*x2^2*Xt^2 + 32*x2^3*Xt^2 - 8*x2^4*Xt^2 - 64*x2*Xt^3 +
       192*x2^2*Xt^3 - 192*x2^3*Xt^3 + 64*x2^4*Xt^3 + 4*Xt^4 +
       14*x2*Xt^4 - 62*x2^2*Xt^4 + 66*x2^3*Xt^4 - 22*x2^4*Xt^4 +
       16*x2*Xt^5 - 32*x2^2*Xt^5 + 16*x2^3*Xt^5 +
       12*(-1 + x2)^5*x2*Log[g/q2]^2 + 11*x2*Log[x2] -
       41*x2^2*Log[x2] + 60*x2^3*Log[x2] - 44*x2^4*Log[x2] +
       17*x2^5*Log[x2] - 3*x2^6*Log[x2] + 24*x2*Xt*Log[x2] -
       96*x2^2*Xt*Log[x2] + 144*x2^3*Xt*Log[x2] -
       96*x2^4*Xt*Log[x2] + 24*x2^5*Xt*Log[x2] -
       22*x2*Xt^2*Log[x2] + 76*x2^2*Xt^2*Log[x2] -
       96*x2^3*Xt^2*Log[x2] + 52*x2^4*Xt^2*Log[x2] -
       10*x2^5*Xt^2*Log[x2] - 32*x2*Xt^3*Log[x2] +
       32*x2^2*Xt^3*Log[x2] + 32*x2^3*Xt^3*Log[x2] -
       32*x2^4*Xt^3*Log[x2] + 19*x2*Xt^4*Log[x2] -
       x2^2*Xt^4*Log[x2] - 47*x2^3*Xt^4*Log[x2] +
       29*x2^4*Xt^4*Log[x2] + 8*x2*Xt^5*Log[x2] +
       16*x2^2*Xt^5*Log[x2] - 24*x2^3*Xt^5*Log[x2] -
       6*x2*Log[x2]^2 + 26*x2^2*Log[x2]^2 - 46*x2^3*Log[x2]^2 +
       42*x2^4*Log[x2]^2 - 20*x2^5*Log[x2]^2 + 4*x2^6*Log[x2]^2 -
       16*x2*Xt*Log[x2]^2 + 56*x2^2*Xt*Log[x2]^2 -
       72*x2^3*Xt*Log[x2]^2 + 40*x2^4*Xt*Log[x2]^2 -
       8*x2^5*Xt*Log[x2]^2 + 4*x2*Xt^2*Log[x2]^2 -
       32*x2^2*Xt^2*Log[x2]^2 + 68*x2^3*Xt^2*Log[x2]^2 -
       56*x2^4*Xt^2*Log[x2]^2 + 16*x2^5*Xt^2*Log[x2]^2 -
       8*x2*Xt^3*Log[x2]^2 + 8*x2^2*Xt^3*Log[x2]^2 +
       8*x2^3*Xt^3*Log[x2]^2 - 8*x2^4*Xt^3*Log[x2]^2 -
       4*x2*Xt^4*Log[x2]^2 + 20*x2^2*Xt^4*Log[x2]^2 -
       2*x2^3*Xt^4*Log[x2]^2 - 10*x2^4*Xt^4*Log[x2]^2 +
       8*x2^2*Xt^5*Log[x2]^2 + 8*x2^3*Xt^5*Log[x2]^2 +
       4*(-1 + x2)^2*Log[g/q2]*(-x2^4 - 2*x2^3*(-2 + Xt^2) -
         (-1 + Xt^2)^2 + x2*(4 - 6*Xt^2 + 8*Xt^3 - 6*Xt^4) +
         x2^2*(-6 + 6*Xt^2 - 8*Xt^3 + 7*Xt^4) +
         x2*(-3 + 3*x2^3 - 4*Xt + 6*Xt^2 + 4*Xt^3 - 5*Xt^4 +
           x2^2*(-9 - 4*Xt + 6*Xt^2) + x2*(9 + 8*Xt - 12*Xt^2 +
             4*Xt^3 - 3*Xt^4))*Log[x2]) + 4*(-1 + x2)^3*x2*
        (2 + 8*Xt + 4*Xt^3 + Xt^4)*PolyLog[2, (-1 + x2)/x2]))/
    (64*Pi^3*(-1 + x2)^5*x2)
    )
];

(* lambda 2-loop threshold correction O(alpha_s alpha_t^2),
   in the limit MU[3,3] = M3 (x2 -> 1), MQ[3,3] != M3 (x1 != 1)
 *)
lambda2LPhiHSSAlphaTAlphaSMUM3Degenerate = With[{
    htSM = Yu[3,3],
    CF = 4/3,
    Nc = 3,
    as = g3^2 / (4 Pi),
    g = M3Input^2,
    x1 = msq2[3,3]/M3Input^2,
    x2 = msu2[3,3]/M3Input^2,
    q2 = SCALE^2,
    Xt = xt / M3Input
    },
    (
       -(as*CF*htSM^4*Nc*(4 - 25*x1 + 63*x1^2 - 82*x1^3 + 58*x1^4 -
       21*x1^5 + 3*x1^6 - 32*x1*Xt + 128*x1^2*Xt - 192*x1^3*Xt +
       128*x1^4*Xt - 32*x1^5*Xt - 8*Xt^2 + 32*x1*Xt^2 -
       48*x1^2*Xt^2 + 32*x1^3*Xt^2 - 8*x1^4*Xt^2 - 64*x1*Xt^3 +
       192*x1^2*Xt^3 - 192*x1^3*Xt^3 + 64*x1^4*Xt^3 + 4*Xt^4 +
       14*x1*Xt^4 - 62*x1^2*Xt^4 + 66*x1^3*Xt^4 - 22*x1^4*Xt^4 +
       16*x1*Xt^5 - 32*x1^2*Xt^5 + 16*x1^3*Xt^5 +
       12*(-1 + x1)^5*x1*Log[g/q2]^2 + 11*x1*Log[x1] -
       41*x1^2*Log[x1] + 60*x1^3*Log[x1] - 44*x1^4*Log[x1] +
       17*x1^5*Log[x1] - 3*x1^6*Log[x1] + 24*x1*Xt*Log[x1] -
       96*x1^2*Xt*Log[x1] + 144*x1^3*Xt*Log[x1] -
       96*x1^4*Xt*Log[x1] + 24*x1^5*Xt*Log[x1] -
       22*x1*Xt^2*Log[x1] + 76*x1^2*Xt^2*Log[x1] -
       96*x1^3*Xt^2*Log[x1] + 52*x1^4*Xt^2*Log[x1] -
       10*x1^5*Xt^2*Log[x1] - 32*x1*Xt^3*Log[x1] +
       32*x1^2*Xt^3*Log[x1] + 32*x1^3*Xt^3*Log[x1] -
       32*x1^4*Xt^3*Log[x1] + 19*x1*Xt^4*Log[x1] -
       x1^2*Xt^4*Log[x1] - 47*x1^3*Xt^4*Log[x1] +
       29*x1^4*Xt^4*Log[x1] + 8*x1*Xt^5*Log[x1] +
       16*x1^2*Xt^5*Log[x1] - 24*x1^3*Xt^5*Log[x1] -
       6*x1*Log[x1]^2 + 26*x1^2*Log[x1]^2 - 46*x1^3*Log[x1]^2 +
       42*x1^4*Log[x1]^2 - 20*x1^5*Log[x1]^2 + 4*x1^6*Log[x1]^2 -
       16*x1*Xt*Log[x1]^2 + 56*x1^2*Xt*Log[x1]^2 -
       72*x1^3*Xt*Log[x1]^2 + 40*x1^4*Xt*Log[x1]^2 -
       8*x1^5*Xt*Log[x1]^2 + 4*x1*Xt^2*Log[x1]^2 -
       32*x1^2*Xt^2*Log[x1]^2 + 68*x1^3*Xt^2*Log[x1]^2 -
       56*x1^4*Xt^2*Log[x1]^2 + 16*x1^5*Xt^2*Log[x1]^2 -
       8*x1*Xt^3*Log[x1]^2 + 8*x1^2*Xt^3*Log[x1]^2 +
       8*x1^3*Xt^3*Log[x1]^2 - 8*x1^4*Xt^3*Log[x1]^2 -
       4*x1*Xt^4*Log[x1]^2 + 20*x1^2*Xt^4*Log[x1]^2 -
       2*x1^3*Xt^4*Log[x1]^2 - 10*x1^4*Xt^4*Log[x1]^2 +
       8*x1^2*Xt^5*Log[x1]^2 + 8*x1^3*Xt^5*Log[x1]^2 +
       4*(-1 + x1)^2*Log[g/q2]*(-x1^4 - 2*x1^3*(-2 + Xt^2) -
         (-1 + Xt^2)^2 + x1*(4 - 6*Xt^2 + 8*Xt^3 - 6*Xt^4) +
         x1^2*(-6 + 6*Xt^2 - 8*Xt^3 + 7*Xt^4) +
         x1*(-3 + 3*x1^3 - 4*Xt + 6*Xt^2 + 4*Xt^3 - 5*Xt^4 +
           x1^2*(-9 - 4*Xt + 6*Xt^2) + x1*(9 + 8*Xt - 12*Xt^2 +
             4*Xt^3 - 3*Xt^4))*Log[x1]) + 4*(-1 + x1)^3*x1*
        (2 + 8*Xt + 4*Xt^3 + Xt^4)*PolyLog[2, (-1 + x1)/x1]))/
    (64*Pi^3*(-1 + x1)^5*x1)
    )
];

(* arXiv:1407.4081, Eq. (12) *)
betaLambda = (
    2 lambdaTree (gYd^2 + gYu^2 + 3 g2d^2 + 3 g2u^2)
    - gYd^4 - gYu^4 - 5 g2d^4 - 5 g2u^4
    - 4 gYd gYu g2d g2u
    - 2 (gYd^2 + g2u^2) (gYu^2 + g2d^2)
);

(* lambda 2-loop threshold correction O(alpha_t^2),
   arxiv:1504.05200, Eq. (21), in the limit of degenerate masses for the stops and the pseudoscalar *)
lambda2LHSSAlphaT2 = With[{
    tan\[Beta] = TanBeta,
    mst = Sqrt[Sqrt[msq2[3,3]*msu2[3,3]]], (* average stop mass *)
    Q = SCALE,
    \[Mu] = MuInput,
    K = -0.1953256,
    Xt = xt,
    gt = Yu[3,3] (* SM Yukawa coupling *)
    },
    (
(3*(1 + tan\[Beta]^2)*(1/2 - 8.34993159891064/(1 + tan\[Beta]^2) +
   (6*\[Mu]^2)/mst^2 +
   (Xt^6*(-1/2 + (1/2 - Log[mst^2/Q^2]/2)/(1 + tan\[Beta]^2) +
      Log[mst^2/Q^2]/2))/mst^6 - 4*Log[mst^2/Q^2] +
   (13*Log[mst^2/Q^2])/(1 + tan\[Beta]^2) - (6*\[Mu]^2*Log[mst^2/Q^2])/
    mst^2 + 3*Log[mst^2/Q^2]^2 - (3*Log[mst^2/Q^2]^2)/(1 + tan\[Beta]^2) +
   (Xt^3*(Xt/mst + (2*\[Mu]*Csc[2*ArcTan[tan\[Beta]]])/mst)*
     (0.8747904000000002/(1 + tan\[Beta]^2) - (2*Log[mst^2/Q^2])/
       (1 + tan\[Beta]^2)))/mst^3 +
   (Xt/mst + (2*\[Mu]*Csc[2*ArcTan[tan\[Beta]]])/mst)^2*
    (0.1252095999999998/(1 + tan\[Beta]^2) + (3*Log[mst^2/Q^2])/
      (1 + tan\[Beta]^2)) +
   (Xt*(Xt/mst + (2*\[Mu]*Csc[2*ArcTan[tan\[Beta]]])/mst)*
     (0.5008383999999992/(1 + tan\[Beta]^2) + (12*Log[mst^2/Q^2])/
       (1 + tan\[Beta]^2)))/mst +
   4*Piecewise[{{-9/4, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
     Re[((-1 + (2*\[Mu]^2)/mst^2 + (2*\[Mu]^4)/mst^4)*
        (-Pi^2/6 - (\[Mu]^2*Log[\[Mu]^2/mst^2])/mst^2 + Log[\[Mu]^2/mst^2]*
          Log[Abs[1 - \[Mu]^2/mst^2]] + PolyLog[2, \[Mu]^2/mst^2]))/
       Abs[(1 - \[Mu]^2/mst^2)^2]]] -
   8*Piecewise[{{-1, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
     (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2))] -
   (2*\[Mu]^2*Piecewise[{{-1, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
      (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2))])/mst^2 +
   (3*\[Mu]^2*Piecewise[{{1/2, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
      (1 + (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2)))/
       (1 - \[Mu]^2/mst^2)])/mst^2 +
   (Xt^2*(-7 + 19.6878144/(1 + tan\[Beta]^2) - (6*\[Mu]^2)/mst^2 +
      27*Log[mst^2/Q^2] - (24*Log[mst^2/Q^2])/(1 + tan\[Beta]^2) +
      (6*\[Mu]^2*Log[mst^2/Q^2])/mst^2 +
      (Xt/mst + (2*\[Mu]*Csc[2*ArcTan[tan\[Beta]]])/mst)^2*
       (-0.021147733333332752/(1 + tan\[Beta]^2) - (3*Log[mst^2/Q^2])/
         (1 + tan\[Beta]^2)) +
      4*Piecewise[{{-1, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
        (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2))] -
      (6*\[Mu]^2*Piecewise[{{-1, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
         (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2))])/mst^2 -
      4*Piecewise[{{1/2, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
        (1 + (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2)))/
         (1 - \[Mu]^2/mst^2)] -
      (6*\[Mu]^2*Piecewise[{{1/2, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
         (1 + (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2)))/
          (1 - \[Mu]^2/mst^2)])/mst^2))/mst^2 +
   (Xt^4*(11/2 - 25/(4*(1 + tan\[Beta]^2)) + \[Mu]^2/mst^2 -
      (13*Log[mst^2/Q^2])/2 + (6*Log[mst^2/Q^2])/(1 + tan\[Beta]^2) -
      (\[Mu]^2*Log[mst^2/Q^2])/mst^2 +
      (Xt/mst + (2*\[Mu]*Csc[2*ArcTan[tan\[Beta]]])/mst)^2*
       (-0.020728533333333354/(1 + tan\[Beta]^2) + Log[mst^2/Q^2]/
         (2*(1 + tan\[Beta]^2))) -
      Piecewise[{{-1, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
        (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2))]/2 +
      (\[Mu]^2*Piecewise[{{-1, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
         (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2))])/mst^2 +
      Piecewise[{{1/2, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
        (1 + (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2)))/
         (1 - \[Mu]^2/mst^2)]/2 +
      (\[Mu]^2*Piecewise[{{1/2, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
         (1 + (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2)))/
          (1 - \[Mu]^2/mst^2)])/(2*mst^2)))/mst^4)*gt^6)/
    (256*Pi^4*tan\[Beta]^2)
    ) /. Piecewise[{{val_, cond_}}, default_] :> If[cond, val, default]
    ];


(* lambda 2-loop threshold correction O(alpha_t^2), computation for generic masses *)
lambda2LHSSAlphaT2Generic = With[{
    k = 1/(4*Pi)^2,
    sbe = Sqrt[TanBeta^2/(1+TanBeta^2)],
    cbe = Sqrt[1/(1+TanBeta^2)],
    Nc = 3, (* number of colors *)
    Q  = msq2[3,3]*(1+0.02),
    U  = msu2[3,3]*(1-0.02),
    mu2 = MuInput^2*(1+0.01),
    q2 = SCALE^2, (* renormalization/matching scale *)
    gt = Yu[3,3], (* SM Yukawa coupling *)
    A0 = mAInput^2*(1-0.01),
    Xt = xt,
    Yt = yt
    },
      gt^6*k^2*(Nc*(5 + 3*Log[Q]^2 + Log[Q]*(-6 - 6*Log[q2]) + 5*Log[q2]^2 +
    Log[q2]*(10 - 4*Log[U]) - 4*Log[U] + 2*Log[U]^2) +
  Nc*Xt^4*((2*(2*Q^2 - 27*Q*U + U^2))/(Q*(Q - U)^2*U) +
    ((-7*Q - 9*U)*Log[Q]^2)/(Q - U)^3 + ((-6*Q^2 - 44*Q*U + 2*U^2)*Log[U])/
     (Q*(Q - U)^3) + ((3*Q + 5*U)*Log[U]^2)/(Q - U)^3 +
    Log[Q]*((-4*(Q^2 - 10*Q*U - 3*U^2))/((Q - U)^3*U) +
      (2*(5*Q + 7*U)*Log[q2])/(Q - U)^3 + (4*(Q + U)*Log[U])/(Q - U)^3) +
    Log[q2]*((2*(2*Q^2 - 15*Q*U + U^2))/(Q*(Q - U)^2*U) -
      (2*(5*Q + 7*U)*Log[U])/(Q - U)^3) - (6*PolyLog[2, 1 - Q/U])/
     (Q - U)^2) + Nc*Xt^2*((-2*Q^2 - 3*Q*U + U^2)/(Q*(Q - U)*U) +
    ((7*Q - 9*U)*Log[Q]^2)/(Q - U)^2 + ((16*Q^2 - 21*Q*U + U^2)*Log[U])/
     (Q*(Q - U)^2) - (9*Log[U]^2)/(Q - U) +
    Log[Q]*((2*Q^2 - 17*Q*U + 19*U^2)/((Q - U)^2*U) -
      (2*(8*Q - 9*U)*Log[q2])/(Q - U)^2 + (2*Q*Log[U])/(Q - U)^2) +
    Log[q2]*((-2*Q^2 - Q*U + U^2)/(Q*(Q - U)*U) + (2*(8*Q - 9*U)*Log[U])/
       (Q - U)^2) - (6*PolyLog[2, 1 - Q/U])/(Q - U)) +
  Nc*Xt^6*((-2*Q^2 - 11*Q*U + U^2)/(Q*(Q - U)^3*U) +
    ((-9*Q - 5*U)*Log[Q]^2)/(Q - U)^4 + ((-10*Q^2 - 3*Q*U + U^2)*Log[U])/
     (Q*(Q - U)^4) + ((-3*Q - 5*U)*Log[U]^2)/(Q - U)^4 +
    Log[q2]*((-2*Q^2 - 5*Q*U + U^2)/(Q*(Q - U)^3*U) -
      (6*Q*Log[U])/(Q - U)^4) +
    Log[Q]*((2*Q^2 + 13*Q*U - 3*U^2)/((Q - U)^4*U) +
      (6*Q*Log[q2])/(Q - U)^4 + (2*(6*Q + 5*U)*Log[U])/(Q - U)^4) -
    (2*PolyLog[2, 1 - Q/U])/(Q - U)^3 + (4*PolyLog[2, (Q - U)/Q])/
     (Q - U)^3) +
  (cbe^2*(Nc*Yt^2*(-Q^(-1) - 2/U +
        ((((Q - U)^4 - A0^3*(Q + U) - A0*(Q - U)^2*(3*Q + 5*U) +
             A0^2*(3*Q^2 + 4*Q*U + 5*U^2))*delta[A0, U, Q] +
           delta[A0, Q, U]*(-(Q*(-A0^2 + (Q - U)^2)*U) +
             (-Q^2 + 3*Q*U - U^2 + A0*(Q + U))*delta[A0, U, Q]))*Log[A0])/
         (2*Q*U^2*delta[A0, Q, U]*delta[A0, U, Q]) +
        (((A0^3 - (Q - U)^3 - 3*A0^2*(Q + U) + 3*A0*(Q^2 - U^2))*
            delta[A0, U, Q] + delta[A0, Q, U]*(-(U*(A0^2 - Q^2 - 2*A0*U +
                U^2)) + (-A0 + Q + 2*U)*delta[A0, U, Q]))*Log[Q])/
         (2*U^2*delta[A0, Q, U]*delta[A0, U, Q]) + (-Q^(-1) - 2/U)*Log[q2] +
        (((A0^3 + (Q - U)^3 - A0^2*(Q + 5*U) - A0*(Q^2 + 4*Q*U - 5*U^2))*
            delta[A0, U, Q] - delta[A0, Q, U]*(2*Q*(A0 + Q - U)*U +
             (A0 + Q - 3*U)*delta[A0, U, Q]))*Log[U])/(2*Q*U*delta[A0, Q, U]*
          delta[A0, U, Q]) + ((A0^4 + (Q - U)^4 - 4*A0^3*(Q + U) -
           4*A0*(Q - U)^2*(Q + U) + 2*A0^2*(3*Q^2 + 2*Q*U + U^2) -
           3*(A0^2 + (Q - U)^2 - 2*A0*(Q + U))*delta[A0, Q, U] +
           2*delta[A0, Q, U]^2)*phi[A0, Q, U])/(2*U^3*delta[A0, Q, U]) +
        ((-(A0 + Q - U)^2 + delta[A0, U, Q])*phi[A0, U, Q])/
         (2*Q*delta[A0, U, Q])) + Nc*(3/2 + Pi^2 - A0/Q - (2*A0)/U +
        3*Log[A0]^2 + 3*Log[Q]^2 + Log[Q]*(-9/2 - 3*Log[q2]) + 2*Log[q2]^2 +
        Log[A0]*(7 + A0*(Q^(-1) + 2/U) - 3*Log[Q] - 3*Log[U]) +
        Log[q2]*(-((A0*(2*Q + U))/(Q*U)) - Log[U]) - (5*Log[U])/2 +
        2*Log[U]^2 + ((A0^2 - 7*A0*Q - delta[A0, Q, Q])*phi[A0, Q, Q])/Q^2 +
        ((A0^2 - 6*A0*U - delta[A0, U, U])*phi[A0, U, U])/U^2) +
      Nc*Xt*Yt*((6*Log[Q]^2)/(Q - U) + (12*Log[U])/(Q - U) +
        (12*Log[q2]*Log[U])/(Q - U) - (4*Log[U]^2)/(Q - U) +
        Log[A0]*((2*Log[Q])/(Q - U) - (2*Log[U])/(Q - U)) +
        Log[Q]*(-12/(Q - U) - (12*Log[q2])/(Q - U) - (2*Log[U])/(Q - U)) -
        (2*(-(A0*(A0 - 7*Q)) + delta[A0, Q, Q])*phi[A0, Q, Q])/
         (Q^2*(Q - U)) + (2*(A0 + Q - U)*phi[A0, U, Q])/(Q*(Q - U)) +
        (2*(-(A0*(A0 - 6*U)) + delta[A0, U, U])*phi[A0, U, U])/
         ((Q - U)*U^2)) + Nc*Xt^3*Yt*(-48/(Q - U)^2 +
        (2*(4*A0 - 5*Q - 3*U)*Log[Q]^2)/(Q - U)^3 - (12*(Q + 3*U)*Log[U])/
         (Q - U)^3 + (4*(-A0 + Q + U)*Log[U]^2)/(Q - U)^3 +
        Log[A0]*((2*(-6*A0 + Q - U)*Log[Q])/(Q - U)^3 +
          (2*(6*A0 - Q + U)*Log[U])/(Q - U)^3) +
        Log[q2]*(-24/(Q - U)^2 - (12*(Q + U)*Log[U])/(Q - U)^3) +
        Log[Q]*((12*(3*Q + U))/(Q - U)^3 + (12*(Q + U)*Log[q2])/(Q - U)^3 +
          (2*(-2*A0 + 3*Q + U)*Log[U])/(Q - U)^3) +
        (2*(A0*(A0 - 7*Q)*(Q - U) + (-5*Q + U)*delta[A0, Q, Q])*
          phi[A0, Q, Q])/(Q^2*(Q - U)^3) -
        (2*((Q - U)*(A0 + Q - U) - 2*delta[A0, U, Q])*phi[A0, U, Q])/
         (Q*(Q - U)^3) + (2*(A0*(A0 - 6*U)*(Q - U) - (Q - 3*U)*
            delta[A0, U, U])*phi[A0, U, U])/((Q - U)^3*U^2)) +
      Xt^2*(Nc*Yt^2*((4*Q - 2*U)/(Q^2*U - Q*U^2) + (3*Log[Q]^2)/(Q - U)^2 +
          ((-((Q - U)*((Q - U)^2*(2*Q - U) + A0^2*(2*Q + U) - 4*A0*Q*
                 (Q + 2*U))*delta[A0, U, Q]) + delta[A0, Q, U]*
              (2*Q*(Q - U)*(A0 + Q - U)*U + (2*Q^2 + Q*U - U^2)*delta[A0, U,
                 Q]))*Log[U])/(Q*(Q - U)^2*U*delta[A0, Q, U]*
            delta[A0, U, Q]) + (2*Log[U]^2)/(Q - U)^2 +
          Log[Q]*(((Q - U)*(-A0^3 + Q*(Q - U)^2 + A0^2*(3*Q + 4*U) -
                A0*(3*Q^2 + 2*Q*U + 7*U^2))*delta[A0, U, Q] +
              delta[A0, Q, U]*(-((Q - U)*U*(-A0^2 + Q^2 + 2*A0*U - U^2)) +
                (-Q^2 + A0*(Q - U) - 2*Q*U + U^2)*delta[A0, U, Q]))/
             ((Q - U)^2*U^2*delta[A0, Q, U]*delta[A0, U, Q]) -
            (2*Log[q2])/(Q - U)^2 - (5*Log[U])/(Q - U)^2) +
          Log[A0]*(((A0^3*Q - (Q - U)^4 + A0*Q*(3*Q^2 - 2*Q*U - U^2) +
                A0^2*(-3*Q^2 - 2*Q*U + U^2))*delta[A0, U, Q] +
              delta[A0, Q, U]*(Q*(-A0^2 + (Q - U)^2)*U + (-(A0*Q) + Q^2 -
                  3*Q*U + U^2)*delta[A0, U, Q]))/(Q*(Q - U)*U^2*
              delta[A0, Q, U]*delta[A0, U, Q]) + Log[Q]/(Q - U)^2 -
            Log[U]/(Q - U)^2) + Log[q2]*((4*Q - 2*U)/(Q^2*U - Q*U^2) +
            (2*Log[U])/(Q - U)^2) + ((A0*(A0 - 7*Q) - delta[A0, Q, Q])*
            phi[A0, Q, Q])/(Q^2*(Q - U)^2) +
          ((-((Q - U)*(A0^4 - 4*A0*Q^2*(Q - U) + (Q - U)^4 - 4*A0^3*(Q + U) +
                A0^2*(6*Q^2 + 4*Q*U + 6*U^2))) + (A0^2*(3*Q - 5*U) +
               (3*Q - 5*U)*(Q - U)^2 + A0*(-6*Q^2 + 4*Q*U + 14*U^2))*
              delta[A0, Q, U] - 2*(Q - 2*U)*delta[A0, Q, U]^2)*phi[A0, Q, U])/
           ((Q - U)^2*U^3*delta[A0, Q, U]) +
          (((Q - U)*(A0 + Q - U)^2 + A0*delta[A0, U, Q])*phi[A0, U, Q])/
           (Q*(Q - U)^2*delta[A0, U, Q]) + ((A0*(A0 - 6*U) - delta[A0, U, U])*
            phi[A0, U, U])/((Q - U)^2*U^2)) +
        Nc*((4*A0*Q - 2*A0*U - 4*Q*U)/(Q^2*U - Q*U^2) +
          ((A0 + Q - 3*U)*Log[Q]^2)/(Q - U)^2 + ((2*A0 + Q - 5*U)*Log[U])/
           (Q - U)^2 - (2*Log[U]^2)/(Q - U) + Log[q2]*
           ((4*A0*Q - 2*A0*U - 2*Q*U)/(Q^2*U - Q*U^2) + (2*(A0 - Q)*Log[U])/
             (Q - U)^2) + Log[Q]*((-2*A0 + Q + 3*U)/(Q - U)^2 -
            (2*(A0 - Q)*Log[q2])/(Q - U)^2 + ((-A0 + Q + U)*Log[U])/
             (Q - U)^2) + Log[A0]*((2*A0*(-2*Q + U))/(Q*(Q - U)*U) +
            ((A0 - 5*Q + 5*U)*Log[Q])/(Q - U)^2 - ((A0 - 5*Q + 5*U)*Log[U])/
             (Q - U)^2) + ((A0*(A0 - 7*Q)*(Q - U) + (-2*Q + U)*
              delta[A0, Q, Q])*phi[A0, Q, Q])/(Q^2*(Q - U)^2) +
          (delta[A0, U, Q]*phi[A0, U, Q])/(Q*(Q - U)^2) +
          ((-(A0*(A0 - 6*U)) + delta[A0, U, U])*phi[A0, U, U])/
           ((Q - U)*U^2))) +
      Xt^4*(Nc*((3*Q*(Q - U)*U + A0*(-2*Q^2 - 5*Q*U + U^2))/(Q*(Q - U)^3*U) +
          Log[Q]*((3*(4*A0*Q - Q^2 + U^2))/(2*(Q - U)^4) +
            (3*(2*A0*Q - Q^2 + U^2)*Log[q2])/(Q - U)^4) +
          (3*(-4*A0*Q + Q^2 - U^2)*Log[U])/(2*(Q - U)^4) +
          Log[q2]*((6*Q*(Q - U)*U + A0*(-2*Q^2 - 5*Q*U + U^2))/
             (Q*(Q - U)^3*U) + (3*(-2*A0*Q + Q^2 - U^2)*Log[U])/(Q - U)^4) +
          Log[A0]*((6*Q*U*(-Q + U) + A0*(2*Q^2 + 5*Q*U - U^2))/
             (Q*(Q - U)^3*U) + (3*(-2*A0*Q + Q^2 - U^2)*Log[Q])/(Q - U)^4 +
            (3*(2*A0*Q - Q^2 + U^2)*Log[U])/(Q - U)^4)) +
        Nc*Yt^2*((-2*Q^2 - 11*Q*U + U^2)/(Q*(Q - U)^3*U) +
          ((7*A0 - 11*Q - 3*U)*Log[Q]^2)/(Q - U)^4 -
          (((Q - U)^2*(A0^3 - 3*Q^3 + 11*Q^2*U - 5*Q*U^2 - 3*U^3 - A0^2*
                (5*Q + 3*U) + A0*(7*Q^2 + 4*Q*U + 5*U^2))*delta[A0, U, Q] +
             delta[A0, Q, U]*(2*Q*(Q - U)^2*(A0 + Q - U)*U + (3*Q^3 -
                 A0*(Q - U)^2 + 7*Q^2*U + 13*Q*U^2 + U^3)*delta[A0, U, Q]))*
            Log[U])/(2*Q*(Q - U)^4*U*delta[A0, Q, U]*delta[A0, U, Q]) -
          (2*(-2*A0 + Q + 3*U)*Log[U]^2)/(Q - U)^4 +
          Log[q2]*((-2*Q^2 - 5*Q*U + U^2)/(Q*(Q - U)^3*U) -
            (6*Q*Log[U])/(Q - U)^4) + Log[A0]*
           (((Q - U)*(-A0^3 + Q^3 - 3*Q^2*U - Q*U^2 + 3*U^3 + 3*A0^2*
                 (Q + U) - A0*(3*Q^2 + 5*U^2))*delta[A0, U, Q] +
              delta[A0, Q, U]*(-(Q*(-A0^2 + (Q - U)^2)*U) +
                (-Q^2 + A0*(Q - U) + 3*Q*U + 3*U^2)*delta[A0, U, Q]))/
             (2*Q*(Q - U)^2*U^2*delta[A0, Q, U]*delta[A0, U, Q]) -
            (3*(A0 - Q + U)*Log[Q])/(Q - U)^4 + (3*(A0 - Q + U)*Log[U])/
             (Q - U)^4) + Log[Q]*(((Q - U)^2*(A0^3 - Q^3 + Q^2*U + 9*Q*U^2 -
                9*U^3 - A0^2*(3*Q + 5*U) + A0*(3*Q^2 + 4*Q*U + 9*U^2))*delta[
                A0, U, Q] + delta[A0, Q, U]*((Q - U)^2*U*(-A0^2 + Q^2 +
                  2*A0*U - U^2) + (Q^3 - A0*(Q - U)^2 + 2*Q^2*U + 17*Q*U^2 +
                  4*U^3)*delta[A0, U, Q]))/(2*(Q - U)^4*U^2*delta[A0, Q, U]*
              delta[A0, U, Q]) + (6*Q*Log[q2])/(Q - U)^4 +
            ((-11*A0 + 13*Q + 9*U)*Log[U])/(Q - U)^4) +
          ((A0*(A0 - 7*Q)*(Q - U) + (-8*Q + U)*delta[A0, Q, Q])*
            phi[A0, Q, Q])/(Q^2*(Q - U)^4) +
          (((Q - U)^2*(A0^4 - 4*A0^3*(Q + U) - 4*A0*(Q - U)^2*(Q + U) +
               (Q - U)^2*(Q^2 - 2*Q*U - 3*U^2) + A0^2*(6*Q^2 + 4*Q*U +
                 6*U^2)) - (Q - U)*(3*Q^3 + 3*A0^2*(Q - 3*U) - 15*Q^2*U + 29*
                Q*U^2 - 17*U^3 - 6*A0*(Q^2 - 2*Q*U - 3*U^2))*
              delta[A0, Q, U] + 2*(Q^2 - 5*Q*U + 12*U^2)*delta[A0, Q, U]^2)*
            phi[A0, Q, U])/(2*(Q - U)^4*U^3*delta[A0, Q, U]) -
          (((Q - U)^2*(A0 + Q - U)^2 + (4*A0 + 3*Q - 3*U)*(Q - U)*
              delta[A0, U, Q] - 6*delta[A0, U, Q]^2)*phi[A0, U, Q])/
           (2*Q*(Q - U)^4*delta[A0, U, Q]) +
          ((A0*(A0 - 6*U)*(-Q + U) + (Q - 5*U)*delta[A0, U, U])*
            phi[A0, U, U])/((Q - U)^4*U^2)))) +
    Nc*((4*mu2^3*(2*Q + U) - Q*U*(4*Q^2 + 3*Q*U + 2*U^2) -
        3*mu2^2*(4*Q^2 + 7*Q*U + 2*U^2) + mu2*(4*Q^3 + 17*Q^2*U + 13*Q*U^2 +
          2*U^3))/(2*(mu2 - Q)*Q*(mu2 - U)*U) +
      ((4*mu2^4 - 4*mu2^3*U - 2*Q^2*U^2 + 2*mu2^2*U*(-4*Q + U) +
         4*mu2*Q*U*(Q + U))*Log[mu2]^2)/((mu2 - Q)^2*(mu2 - U)^2) +
      (2*(2*mu2 - Q)*Q*Log[Q]^2)/(mu2 - Q)^2 - 2*Log[q2]^2 +
      ((-(Q*U^2*(Q + 2*U)) + mu2^3*(13*Q + 2*U) +
         mu2*U*(4*Q^2 + 9*Q*U + 2*U^2) - mu2^2*(9*Q^2 + 14*Q*U + 4*U^2))*
        Log[U])/(2*(mu2 - Q)*Q*(mu2 - U)^2) + ((2*mu2 - U)*U*Log[U]^2)/
       (mu2 - U)^2 + Log[q2]*(-(((-2*mu2 + Q + U)*(2*Q + U))/(Q*U)) +
        Log[U]) + Log[mu2]*((mu2*(-2*Q^2*U^2*(Q + U) - 2*mu2^4*(2*Q + U) +
           mu2*Q*U*(3*Q^2 + 8*Q*U - 2*U^2) + mu2^3*(8*Q^2 + Q*U + 4*U^2) -
           2*mu2^2*(2*Q^3 + 4*Q^2*U - Q*U^2 + U^3)))/((mu2 - Q)^2*Q*
          (mu2 - U)^2*U) + ((-5*mu2^4 - 2*mu2^3*(Q - 4*U) + 2*Q^2*U^2 -
           4*mu2*Q*U*(Q + U) + mu2^2*(Q^2 + 8*Q*U - 4*U^2))*Log[Q])/
         ((mu2 - Q)^2*(mu2 - U)^2) + ((-3*mu2^4 + 2*mu2^3*Q -
           mu2^2*Q*(Q - 8*U) + 2*Q^2*U^2 - 4*mu2*Q*U*(Q + U))*Log[U])/
         ((mu2 - Q)^2*(mu2 - U)^2)) +
      Log[Q]*((-(Q^2*U*(4*Q + 5*U)) + mu2^3*(4*Q + 15*U) +
          mu2*Q*(4*Q^2 + 15*Q*U + 6*U^2) - mu2^2*(8*Q^2 + 14*Q*U + 13*U^2))/
         (2*(mu2 - Q)^2*(mu2 - U)*U) + 3*Log[q2] +
        ((2*mu2^4 - 2*mu2^3*U - Q^2*U^2 + mu2^2*U*(-4*Q + U) +
           2*mu2*Q*U*(Q + U))*Log[U])/((mu2 - Q)^2*(mu2 - U)^2)) +
      (4*mu2^2*PolyLog[2, (mu2 - Q)/mu2])/(mu2 - Q)^2 +
      (2*(mu2^2 + 2*mu2*Q - Q^2)*PolyLog[2, 1 - Q/mu2])/(mu2 - Q)^2 +
      (2*(mu2^2 + 2*mu2*U - U^2)*PolyLog[2, 1 - U/mu2])/(mu2 - U)^2) +
    Xt^2*(Nc^2*((4*Q*Log[Q]^2)/(Q - U)^2 + (4*Log[U])/(Q - U) +
        (4*Log[q2]*Log[U])/(Q - U) + (4*U*Log[U]^2)/(Q - U)^2 +
        Log[Q]*(-4/(Q - U) - (4*Log[q2])/(Q - U) - (4*(Q + U)*Log[U])/
           (Q - U)^2)) + Nc*((-8*mu2*Q + 4*Q^2 + 4*mu2*U + 6*Q*U - 2*U^2)/
         (Q*(Q - U)*U) + (2*Q*(3*mu2^2 - 2*mu2*(Q + 2*U) + Q*(Q + 2*U))*
          Log[Q]^2)/((mu2 - Q)^2*(Q - U)^2) +
        ((-4*mu2^3*Q + Q*U*(7*Q^2 + 3*Q*U - 2*U^2) +
           mu2^2*(Q^2 + 17*Q*U - 2*U^2) - mu2*(Q^3 + 16*Q^2*U + 5*Q*U^2 -
             2*U^3))*Log[U])/((mu2 - Q)*Q*(mu2 - U)*(Q - U)^2) +
        (2*U*(2*mu2^2 - 2*mu2*(Q + U) + U*(Q + U))*Log[U]^2)/
         ((mu2 - U)^2*(Q - U)^2) + Log[q2]*
         ((-8*mu2*Q + 4*Q^2 + 4*mu2*U + 4*Q*U - 2*U^2)/(Q*(Q - U)*U) +
          (2*(-2*mu2 + 2*Q + U)*Log[U])/(Q - U)^2) +
        Log[Q]*((4*mu2^3*U + Q*U*(-4*Q^2 - 7*Q*U + 3*U^2) -
            mu2^2*(4*Q^2 + Q*U + 11*U^2) + mu2*(4*Q^3 + 9*Q^2*U + 2*Q*U^2 +
              5*U^3))/((mu2 - Q)*(mu2 - U)*(Q - U)^2*U) -
          (2*(-2*mu2 + 2*Q + U)*Log[q2])/(Q - U)^2 +
          ((-6*Q)/(Q - U)^2 + (4*Q*(-2*mu2 + Q))/((mu2 - Q)^2*(Q - U)) -
            (4*U)/(Q - U)^2 + (2*U*(-2*mu2 + U))/((mu2 - U)^2*(-Q + U)))*
           Log[U]) + Log[mu2]*((4*mu2*(mu2^2*(2*Q - U) + Q^2*U +
             mu2*(-2*Q^2 - Q*U + U^2)))/((mu2 - Q)*Q*(mu2 - U)*(Q - U)*U) +
          (2*(-2*mu2^5 + 2*Q^3*U^2 - 2*mu2*Q^2*U*(2*Q + 3*U) +
             mu2^4*(3*Q + 7*U) - 2*mu2^3*(2*Q^2 + 5*Q*U + 3*U^2) +
             mu2^2*(Q^3 + 13*Q^2*U + 4*Q*U^2 + 2*U^3))*Log[Q])/
           ((mu2 - Q)^2*(mu2 - U)^2*(Q - U)^2) +
          (2*(2*mu2^5 - 2*Q^3*U^2 + 2*mu2*Q^2*U*(2*Q + 3*U) -
             mu2^4*(3*Q + 7*U) + 2*mu2^3*(2*Q^2 + 5*Q*U + 3*U^2) -
             mu2^2*(Q^3 + 13*Q^2*U + 4*Q*U^2 + 2*U^3))*Log[U])/
           ((mu2 - Q)^2*(mu2 - U)^2*(Q - U)^2)) +
        (4*(mu2 - Q)*PolyLog[2, (mu2 - Q)/mu2])/(Q - U)^2 -
        (4*(mu2 - Q)*PolyLog[2, 1 - U/mu2])/(Q - U)^2)) +
    Xt^4*(Nc^2*(-8/(Q - U)^2 - (4*Q*(Q + U)*Log[Q]^2)/(Q - U)^4 -
        (4*(Q + 3*U)*Log[U])/(Q - U)^3 - (4*U*(Q + U)*Log[U]^2)/(Q - U)^4 +
        Log[q2]*(-8/(Q - U)^2 - (4*(Q + U)*Log[U])/(Q - U)^3) +
        Log[Q]*((4*(3*Q + U))/(Q - U)^3 + (4*(Q + U)*Log[q2])/(Q - U)^3 +
          (4*(Q + U)^2*Log[U])/(Q - U)^4)) +
      Nc*((2*mu2^3*(2*Q^2 + 8*Q*U - U^2) - mu2^2*(6*Q^3 + 38*Q^2*U +
            13*Q*U^2 - 3*U^3) + Q*U*(-2*Q^3 - 24*Q^2*U + 7*Q*U^2 + U^3) +
          mu2*(2*Q^4 + 28*Q^3*U + 31*Q^2*U^2 - 6*Q*U^3 - U^4))/
         ((mu2 - Q)*Q*(mu2 - U)*(Q - U)^3*U) -
        (2*Q*(mu2^2*(4*Q + 3*U) - 2*mu2*(3*Q^2 + 3*Q*U + U^2) +
           Q*(3*Q^2 + 3*Q*U + U^2))*Log[Q]^2)/((mu2 - Q)^2*(Q - U)^4) +
        ((12*mu2^4*Q*(2*Q + U) + Q*U^2*(25*Q^3 + 34*Q^2*U - 21*Q*U^2 -
             2*U^3) - mu2^3*(39*Q^3 + 102*Q^2*U + 5*Q*U^2 - 2*U^3) +
           mu2^2*(19*Q^4 + 130*Q^3*U + 93*Q^2*U^2 - 22*Q*U^3 - 4*U^4) +
           mu2*U*(-48*Q^4 - 113*Q^3*U - 6*Q^2*U^2 + 21*Q*U^3 + 2*U^4))*
          Log[U])/(2*(mu2 - Q)*Q*(mu2 - U)^2*(Q - U)^4) -
        (U*(Q + U)*(4*mu2^2 - 2*mu2*(Q + 3*U) + U*(Q + 3*U))*Log[U]^2)/
         ((mu2 - U)^2*(Q - U)^4) + Log[q2]*
         ((-2*Q^3 - 13*Q^2*U + 2*Q*U^2 + U^3 + 2*mu2*(2*Q^2 + 5*Q*U - U^2))/
           (Q*(Q - U)^3*U) + (3*(4*mu2*Q - 3*Q^2 - 2*Q*U + U^2)*Log[U])/
           (Q - U)^4) + Log[Q]*((-36*mu2^4*Q*U - Q^2*U*(4*Q^3 + 47*Q^2*U +
              4*Q*U^2 - 19*U^3) + mu2^3*(4*Q^3 + 121*Q^2*U + 24*Q*U^2 -
              5*U^3) + mu2*Q*(4*Q^4 + 53*Q^3*U + 130*Q^2*U^2 - 9*Q*U^3 -
              34*U^4) + mu2^2*(-8*Q^4 - 126*Q^3*U - 131*Q^2*U^2 + 42*Q*U^3 +
              7*U^4))/(2*(mu2 - Q)^2*(mu2 - U)*(Q - U)^4*U) -
          (3*(4*mu2*Q - 3*Q^2 - 2*Q*U + U^2)*Log[q2])/(Q - U)^4 +
          ((10*Q)/(Q - U)^3 - 2/(Q - U)^2 + (22*Q*U)/(Q - U)^4 +
            (6*U)/(-Q + U)^3 + (2*(2*mu2 - Q)*Q*(Q + U))/((mu2 - Q)^2*
              (Q - U)^3) + ((2*mu2 - U)*U*(Q + U))/((mu2 - U)^2*(-Q + U)^3))*
           Log[U]) + Log[mu2]*((-2*mu2*(-3*mu2*Q^3*U*(Q + 3*U) +
             mu2^4*(2*Q^2 + 2*Q*U - U^2) + Q^2*U^2*(2*Q^2 + 2*Q*U - U^2) -
             mu2^3*(4*Q^3 + 5*Q^2*U + 5*Q*U^2 - 2*U^3) +
             mu2^2*(2*Q^4 + 8*Q^3*U + 7*Q^2*U^2 + 2*Q*U^3 - U^4)))/
           ((mu2 - Q)^2*Q*(mu2 - U)^2*(Q - U)^3*U) +
          ((12*mu2^5*Q + 2*Q^2*U^2*(-2*Q^2 - 2*Q*U + U^2) -
             mu2^4*(25*Q^2 + 28*Q*U + U^2) + 6*mu2^3*Q*(3*Q^2 + 10*Q*U + 3*
                U^2) - mu2^2*Q*(3*Q^3 + 44*Q^2*U + 41*Q*U^2 - 4*U^3) +
             4*mu2*Q*U*(2*Q^3 + 7*Q^2*U + Q*U^2 - U^3))*Log[Q])/
           ((mu2 - Q)^2*(mu2 - U)^2*(Q - U)^4) +
          ((-12*mu2^5*Q + 2*Q^2*U^2*(2*Q^2 + 2*Q*U - U^2) +
             mu2^4*(25*Q^2 + 28*Q*U + U^2) - 6*mu2^3*Q*(3*Q^2 + 10*Q*U + 3*
                U^2) + mu2^2*Q*(3*Q^3 + 44*Q^2*U + 41*Q*U^2 - 4*U^3) -
             4*mu2*Q*U*(2*Q^3 + 7*Q^2*U + Q*U^2 - U^3))*Log[U])/
           ((mu2 - Q)^2*(mu2 - U)^2*(Q - U)^4)) +
        (2*(3*mu2^2 - 6*mu2*Q + 2*Q^2 + 2*Q*U - U^2)*PolyLog[2,
           (mu2 - Q)/mu2])/(Q - U)^4 + (2*(-3*mu2^2 + 6*mu2*Q - 2*Q^2 -
           2*Q*U + U^2)*PolyLog[2, 1 - U/mu2])/(Q - U)^4)))/sbe^2) /. { delta[x_,y_,z_] :> TDelta[x,y,z], phi[x_,y_,z_] :> TPhi[x,y,z] }
];


    
(* Top and bottom Yukawa lambda 2-loop threshold correction, all masses degenerate and equal to the matching scale *)
lambda2LHSSAlphaTAlphaBAllDegenerate = With[{
    k = 1/(4*Pi)^2,
    sbe = Sqrt[TanBeta^2/(1+TanBeta^2)],
    cbe = Sqrt[1/(1+TanBeta^2)],
    Nc = 3, (* number of colors *)
    CF = 4/3,
    gt = Yu[3,3], (* SM Yukawa coupling *)
    gb = Yd[3,3],
    Xt = xt,
    Yt = yt,
    Xb = xb,
    Yb = yb,
    M  = Sqrt[msq2[3,3]*msu2[3,3]],
    K = -0.1953256,
    sgn = MuInput/Abs[MuInput],
    (* For the bottom Yukawa coupling *)
    Mu = MuInput,
    M3 = M3Input,
    mQ3 = Sqrt[msq2[3,3]],
    mU3 = Sqrt[msu2[3,3]],
    mD3 = Sqrt[msd2[3,3]],
    Xtildet = xtt,
    Xtildeb = xbb,
    Q2 = SCALE^2,
    cosb = Cos[ArcTan[TanBeta]],
    sinb = Sin[ArcTan[TanBeta]],
    MA = mAInput
    },
    ybMSSM[mQ3_,mU3_,mD3_,M3_,Mu_,TanBeta_,Xt_,Xb_] := Module[{deltagsb, deltagbyL1, deltagbyL2, deltagbyL3, deltagbyL4},
          deltagsb = - g3^2*CF*k*(1+Log[M3^2/Q2]+TCF[6][mQ3/M3]+TCF[6][mD3/M3]-Xb/M3*TCF[9][mQ3/M3,mD3/M3]);
          deltagbyL1 = - gb^2/cosb^2*k*(3/4*Log[Mu^2/Q2]+3/8*sinb^2*(2*Log[MA^2/Q2]-1)+TCF[6][mQ3/Mu]+1/2*TCF[6][mD3/Mu]);
          deltagbyL2 = - gt^2/sinb^2*k*(1/4*Log[Mu^2/Q2]+1/8*cosb^2*(2*Log[MA^2/Q2]-1)+sinb^2*(Log[MA^2/Q2]-1));
          deltagbyL3 = - gt^2/sinb^2*k*(1/2*TCF[6][mU3/Mu]+(Xt*TanBeta)/(2*Mu)*TCF[9][mQ3/Mu,mU3/Mu] );
          deltagbyL4 = - 1/2*(-gt^2*Nc*k*Xtildet/6*TCF[5][xQU]-gb^2*Nc*k*Xtildeb/6*TCF[5][xQD]);
          (gb/(1-deltagsb-(deltagbyL1+deltagbyL2+deltagbyL3+deltagbyL4)))
    ];
    k^2*( ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^6*(5*Nc - (14*Nc*Xb^2)/M + ((11*Nc)/(2*M^2) - (2*Nc^2)/(3*M^2))*Xb^4 +
    (-Nc/(2*M^3) + Nc^2/(18*M^3))*Xb^6 +
    (4*Nc - (2*Nc*Xb^2)/M + (Nc*Xb^4)/M^2)/cbe^2 +
    (sbe^2*(Nc*(-3/2 + 60*K + Pi^2) + (-12/M - (64*K)/M)*Nc*Xb*Yb +
       (4/M^2 + (16*K)/M^2)*Nc*Xb^3*Yb + (-3/M - (16*K)/M)*Nc*Yb^2 +
       Xb^4*(-Nc/(2*M^2) + (-19/(12*M^3) - (8*K)/M^3)*Nc*Yb^2) +
       Xb^2*((-2/M - (24*K)/M)*Nc + (14/(3*M^2) + (24*K)/M^2)*Nc*Yb^2)))/
     cbe^2) + ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^4*gt^2*(Nc*(-1 + 48*K + (4*Pi^2)/3) +
    (-10/M - (32*K)/M)*Nc*Xb^2 + (Nc*Xb^4)/(6*M^2) +
    ((4*Nc)/M + (-(Nc/M^2) - (2*Nc^2)/(3*M^2))*Xb^2 +
      (Nc/(6*M^3) + Nc^2/(18*M^3))*Xb^4)*Xt^2 +
    ((-2*Nc*sgn*Xb^3)/(3*M^(3/2)) + (4*Nc*sgn*Xt)/Sqrt[M])/(cbe*sbe) +
    (-8*Nc - (Nc*Xb^2)/M + (-(Nc/M) + (Nc*Xb^2)/(2*M^2))*Xt^2)/cbe^2 +
    (-4/M - (32*K)/M)*Nc*Xb*Yb + (8/(9*M^2) + (16*K)/(3*M^2))*Nc*Xb^3*Yb +
    (sbe^2*(24*K*Nc + (16*K*Nc*Xb*Yb)/M + (8*K*Nc*Yb^2)/M +
       Xt*((16*K*Nc*Xb)/M + (16*K*Nc*Yb)/M + (2*Nc*Xb^2*Yb)/(3*M^2)) +
       Xb^2*((-M^(-1) - (8*K)/M)*Nc + (Nc*Yb^2)/(3*M^2)) +
       Xt^2*((-M^(-1) - (8*K)/M)*Nc + (2*Nc*Xb*Yb)/(3*M^2) +
         (Nc*Yb^2)/(3*M^2) + (-11/(18*M^3) - (8*K)/(3*M^3))*Nc*Xb^2*Yb^2)))/
     cbe^2 + ((-4/M - (32*K)/M)*Nc*Xb + (8/(9*M^2) + (16*K)/(3*M^2))*Nc*Xb^3)*
     Yt + Xt*((16*K*Nc*Xb)/M + (16*K*Nc*Yb)/M + (2*Nc*Xb^2*Yb)/(3*M^2) +
      ((16*K*Nc)/M + (2*Nc*Xb^2)/(3*M^2) + (4*Nc*Xb*Yb)/(3*M^2) +
        (2/(27*M^3) + (16*K)/(9*M^3))*Nc*Xb^3*Yb)*Yt) +
    (4*Nc + (Nc*Xb^4)/(6*M^2) + cbe^2*(Nc*(12*K + (-3 + 2*Pi^2)/6) -
        (8*K*Nc*Xb^2)/M - (Nc*Xb^4)/(6*M^2) + ((-4/M - (32*K)/M)*Nc*Xb +
          (8/(9*M^2) + (16*K)/(3*M^2))*Nc*Xb^3)*Yt +
        ((-M^(-1) - (8*K)/M)*Nc + (4/(3*M^2) + (8*K)/M^2)*Nc*Xb^2 +
          (-35/(108*M^3) - (16*K)/(9*M^3))*Nc*Xb^4)*Yt^2))/sbe^2) +
     + ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2*gt^4*(Nc*(-1 + 48*K + (4*Pi^2)/3) + (4*Nc*Xb^2)/M +
    (-Nc/(2*M^2) + (Nc*Xb^2)/(6*M^3))*Xt^4 + (4*Nc + (Nc*Xt^4)/(6*M^2))/
     cbe^2 + ((4*Nc*sgn*Xb)/Sqrt[M] - (4*Nc*sgn*Xb*Xt^2)/M^(3/2) -
      (2*Nc*sgn*Xt^3)/(3*M^(3/2)) + (Nc*sgn*Xb*Xt^4)/(3*M^(5/2)))/(cbe*sbe) +
    (16*K*Nc*Xb*Yb)/M + (sbe^2*(Nc*(12*K + (-3 + 2*Pi^2)/6) +
       (-4/M - (32*K)/M)*Nc*Xt*Yb + (8/(9*M^2) + (16*K)/(3*M^2))*Nc*Xt^3*Yb +
       (-M^(-1) - (8*K)/M)*Nc*Yb^2 + Xt^4*(-Nc/(4*M^2) +
         (-35/(108*M^3) - (16*K)/(9*M^3))*Nc*Yb^2) +
       Xt^2*((M^(-1) - (8*K)/M)*Nc + (4/(3*M^2) + (8*K)/M^2)*Nc*Yb^2)))/
     cbe^2 + (16*K*Nc*Xb*Yt)/M + Xt^2*((-2/M - (32*K)/M)*Nc - (Nc*Xb^2)/M^2 +
      (2*Nc*Xb*Yb)/(3*M^2) + (2*Nc*Xb*Yt)/(3*M^2)) +
    Xt^3*((8/(9*M^2) + (16*K)/(3*M^2))*Nc*Yb +
      ((8/(9*M^2) + (16*K)/(3*M^2))*Nc + (2/(27*M^3) + (16*K)/(9*M^3))*Nc*Xb*
         Yb)*Yt) + Xt*((16*K*Nc*Xb)/M + (-4/M - (32*K)/M)*Nc*Yb +
      ((-4/M - (32*K)/M)*Nc + (4*Nc*Xb*Yb)/(3*M^2))*Yt) +
    (-8*Nc - (Nc*Xb^2)/M + (-(Nc/M) + (Nc*Xb^2)/(2*M^2))*Xt^2 +
      cbe^2*(24*K*Nc + (-M^(-1) - (8*K)/M)*Nc*Xb^2 + (16*K*Nc*Xb*Yt)/M +
        ((8*K*Nc)/M + (Nc*Xb^2)/(3*M^2))*Yt^2 +
        Xt*((16*K*Nc*Xb)/M + ((16*K*Nc)/M + (2*Nc*Xb^2)/(3*M^2))*Yt) +
        Xt^2*((-M^(-1) - (8*K)/M)*Nc + (2*Nc*Xb*Yt)/(3*M^2) +
          (Nc/(3*M^2) + (-11/(18*M^3) - (8*K)/(3*M^3))*Nc*Xb^2)*Yt^2)))/
     sbe^2))
];

(* Top and bottom Yukawa lambda 2-loop threshold correction, computation with general masses *)
lambda2LHSSAlphaTAlphaBGeneric = With[{
    k = 1/(4*Pi)^2,
    sbe = Sqrt[TanBeta^2/(1+TanBeta^2)],
    cbe = Sqrt[1/(1+TanBeta^2)],
    Nc = 3, (* number of colors *)
    Q  = msq2[3,3]*(1+0.02),
    U  = msu2[3,3]*(1-0.02),
    mD  = msd2[3,3]*(1-0.04),
    mu2 = MuInput^2*(1+0.03),
    mQ3 = Sqrt[msq2[3,3]],
    mU3 = Sqrt[msu2[3,3]],
    mD3 = Sqrt[msd2[3,3]],
    M3 = M3Input,
    Mu = MuInput,
    sgn = MuInput/Abs[MuInput],      
    q2 = SCALE^2, (* renormalization/matching scale *)
    Q2 = SCALE^2,
    Xtildet = xtt,
    Xtildeb = xbb,
    cosb = Cos[ArcTan[TanBeta]],
    sinb = Sin[ArcTan[TanBeta]],
    MA = mAInput,
    CF = 4/3,
    gt = Yu[3,3], (* SM Yukawa coupling *)
    gb = Yd[3,3],
    A0 = mAInput^2*(1-0.02),
    Xt = xt,
    Yt = yt,
    Xb = xb,
    Yb = yb
    },
    ybMSSM[mQ3_,mU3_,mD3_,M3_,Mu_,TanBeta_,Xt_,Xb_] := Module[{deltagsb, deltagbyL1, deltagbyL2, deltagbyL3, deltagbyL4},
          deltagsb = - g3^2*CF*k*(1+Log[M3^2/Q2]+TCF[6][mQ3/M3]+TCF[6][mD3/M3]-Xb/M3*TCF[9][mQ3/M3,mD3/M3]);
          deltagbyL1 = - gb^2/cosb^2*k*(3/4*Log[Mu^2/Q2]+3/8*sinb^2*(2*Log[MA^2/Q2]-1)+TCF[6][mQ3/Mu]+1/2*TCF[6][mD3/Mu]);
          deltagbyL2 = - gt^2/sinb^2*k*(1/4*Log[Mu^2/Q2]+1/8*cosb^2*(2*Log[MA^2/Q2]-1)+sinb^2*(Log[MA^2/Q2]-1));
          deltagbyL3 = - gt^2/sinb^2*k*(1/2*TCF[6][mU3/Mu]+(Xt*TanBeta)/(2*Mu)*TCF[9][mQ3/Mu,mU3/Mu] );
          deltagbyL4 = - 1/2*(-gt^2*Nc*k*Xtildet/6*TCF[5][xQU]-gb^2*Nc*k*Xtildeb/6*TCF[5][xQD]);
          (gb/(1-deltagsb-(deltagbyL1+deltagbyL2+deltagbyL3+deltagbyL4)))
      ];
k^2*(ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^6*(Nc*(5 + 2*Log[mD]^2 + 3*Log[Q]^2 + Log[Q]*(-6 - 6*Log[q2]) + 
      Log[mD]*(-4 - 4*Log[q2]) + 10*Log[q2] + 5*Log[q2]^2) + 
    (sbe^2*(Nc*Yb^2*(-2/mD - Q^(-1) + ((-(A0^3*(mD + 2*Q)) + 
            (mD - Q)^2*(mD^2 - 4*mD*Q + 2*Q^2) + A0^2*(5*mD^2 + 6*mD*Q + 
              6*Q^2) + A0*(-5*mD^3 + 8*mD^2*Q + 3*mD*Q^2 - 6*Q^3) + 
            ((A0 - mD)*mD + 2*(A0 + 2*mD)*Q - 2*Q^2)*delta[A0, Q, mD])*
           Log[A0])/(2*mD^2*Q*delta[A0, Q, mD]) + 
         ((A0^3 + A0^2*(-5*mD + Q) + 5*A0*(mD^2 - 2*mD*Q - Q^2) - 
            (mD - Q)*(mD^2 - 6*mD*Q + 3*Q^2) - (A0 - 3*mD + 3*Q)*
             delta[A0, Q, mD])*Log[mD])/(2*mD*Q*delta[A0, Q, mD]) - 
         ((-2*A0^3 + A0^2*(7*mD + 6*Q) + (mD - Q)*(mD^2 + 3*mD*Q - 2*Q^2) - 
            2*A0*(mD^2 + mD*Q + 3*Q^2) + (2*A0 - 3*mD - 2*Q)*
             delta[A0, Q, mD])*Log[Q])/(2*mD^2*delta[A0, Q, mD]) + 
         (-2/mD - Q^(-1))*Log[q2] + ((2*A0^4 - 8*A0^3*(mD + Q) - 
            2*A0*(mD - Q)*(3*mD^2 - 4*Q^2) + (mD - Q)^2*(mD^2 - 4*mD*Q + 
              2*Q^2) + A0^2*(7*mD^2 + 8*mD*Q + 12*Q^2) + delta[A0, Q, mD]*
             (-6*A0^2 - 5*mD^2 + 12*mD*Q - 6*Q^2 + 12*A0*(mD + Q) + 
              4*delta[A0, Q, mD]))*phi[A0, Q, mD])/
          (2*mD^3*delta[A0, Q, mD])) + Nc*(3/2 - (2*A0)/mD + Pi^2 - A0/Q + 
         3*Log[A0]^2 + 2*Log[mD]^2 + 3*Log[Q]^2 + Log[Q]*(-6 - 6*Log[q2]) + 
         Log[A0]*(7 + A0*(2/mD + Q^(-1)) - 6*Log[q2]) + 
         Log[mD]*(-4 - 4*Log[q2]) + (3 - A0*(2/mD + Q^(-1)))*Log[q2] + 
         8*Log[q2]^2 + ((A0^2 - 6*A0*mD - delta[A0, mD, mD])*phi[A0, mD, mD])/
          mD^2 + ((2*A0^2 - 11*A0*Q - 2*delta[A0, Q, Q])*phi[A0, Q, Q])/
          Q^2) + Nc*Xb*Yb*((4*Log[mD]^2)/(mD - Q) + (6*Log[Q]^2)/(-mD + Q) + 
         Log[A0]*((2*Log[mD])/(mD - Q) + (2*Log[Q])/(-mD + Q)) + 
         Log[Q]*(12/(mD - Q) + (12*Log[q2])/(mD - Q)) + 
         Log[mD]*(12/(-mD + Q) + (2*Log[Q])/(mD - Q) + (12*Log[q2])/
            (-mD + Q)) - (2*(-(A0*(A0 - 6*mD)) + delta[A0, mD, mD])*
           phi[A0, mD, mD])/(mD^2*(mD - Q)) + 
         (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2 - 
            delta[A0, Q, mD])*phi[A0, Q, mD])/(mD^2*(mD - Q)) + 
         ((-4*A0^2 + 22*A0*Q + 4*delta[A0, Q, Q])*phi[A0, Q, Q])/
          ((mD - Q)*Q^2)) + Nc*Xb^3*Yb*(-48/(mD - Q)^2 + 
         (4*(A0 - mD - Q)*Log[mD]^2)/(mD - Q)^3 + 
         (2*(-4*A0 + 3*mD + 5*Q)*Log[Q]^2)/(mD - Q)^3 + 
         Log[A0]*((-2*(6*A0 + mD - Q)*Log[mD])/(mD - Q)^3 + 
           (2*(6*A0 + mD - Q)*Log[Q])/(mD - Q)^3) - (24*Log[q2])/(mD - Q)^2 + 
         Log[Q]*((-12*(mD + 3*Q))/(mD - Q)^3 - (12*(mD + Q)*Log[q2])/
            (mD - Q)^3) + Log[mD]*((12*(3*mD + Q))/(mD - Q)^3 - 
           (2*(-2*A0 + mD + 3*Q)*Log[Q])/(mD - Q)^3 + (12*(mD + Q)*Log[q2])/
            (mD - Q)^3) + (2*(A0*(A0 - 6*mD)*(mD - Q) + (-3*mD + Q)*
             delta[A0, mD, mD])*phi[A0, mD, mD])/(mD^2*(mD - Q)^3) + 
         (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2) + 
            (-3*mD + Q)*delta[A0, Q, mD])*phi[A0, Q, mD])/(mD^2*(mD - Q)^3) + 
         ((2*A0*(2*A0 - 11*Q)*(mD - Q) - 4*(mD - 3*Q)*delta[A0, Q, Q])*
           phi[A0, Q, Q])/((mD - Q)^3*Q^2)) + 
       Xb^4*(Nc*((A0*(-mD^2 + 5*mD*Q + 2*Q^2))/(mD*(mD - Q)^3*Q) + 
           Log[A0]*((A0*(mD^2 - 5*mD*Q - 2*Q^2))/(mD*(mD - Q)^3*Q) + 
             (6*A0*Q*Log[mD])/(mD - Q)^4 - (6*A0*Q*Log[Q])/(mD - Q)^4) + 
           (A0*(-mD^2 + 5*mD*Q + 2*Q^2)*Log[q2])/(mD*(mD - Q)^3*Q) + 
           Log[mD]*((-6*A0*Q)/(mD - Q)^4 - (6*A0*Q*Log[q2])/(mD - Q)^4) + 
           Log[Q]*((6*A0*Q)/(mD - Q)^4 + (6*A0*Q*Log[q2])/(mD - Q)^4)) + 
         Nc*Yb^2*((-mD^2 + 11*mD*Q + 2*Q^2)/(mD*(mD - Q)^3*Q) - 
           (2*(-2*A0 + 3*mD + Q)*Log[mD]^2)/(mD - Q)^4 + 
           ((7*A0 - 3*mD - 11*Q)*Log[Q]^2)/(mD - Q)^4 + 
           Log[A0]*((A0^3*(mD - 2*Q) - (mD - Q)^2*(3*mD^2 + 4*mD*Q - 2*Q^2) + 
               A0*(mD - Q)*(5*mD^2 + mD*Q + 6*Q^2) + A0^2*(-3*mD^2 + 2*mD*Q + 
                 6*Q^2) + (-(A0*mD) + 3*mD^2 + 2*A0*Q + 4*mD*Q - 2*Q^2)*
                delta[A0, Q, mD])/(2*mD^2*(mD - Q)^2*Q*delta[A0, Q, mD]) + 
             (3*(A0 + mD - Q)*Log[mD])/(mD - Q)^4 - (3*(A0 + mD - Q)*Log[Q])/
              (mD - Q)^4) + ((-mD^2 + 5*mD*Q + 2*Q^2)*Log[q2])/
            (mD*(mD - Q)^3*Q) + Log[mD]*(((mD - Q)^2*(-A0^3 + 
                 A0^2*(3*mD + 7*Q) + (mD - Q)*(3*mD^2 + 12*mD*Q - 5*Q^2) - 
                 A0*(5*mD^2 + 10*mD*Q + 11*Q^2)) - (mD^3 - A0*(mD - Q)^2 + 
                 15*mD^2*Q + 3*mD*Q^2 + 5*Q^3)*delta[A0, Q, mD])/
              (2*mD*(mD - Q)^4*Q*delta[A0, Q, mD]) + ((-11*A0 + 9*mD + 13*Q)*
               Log[Q])/(mD - Q)^4 - (6*Q*Log[q2])/(mD - Q)^4) + 
           Log[Q]*(((mD - Q)^2*(2*A0^3 - 3*A0^2*(3*mD + 2*Q) - (mD - Q)*
                  (11*mD^2 + mD*Q - 2*Q^2) + 2*A0*(7*mD^2 + 3*mD*Q + 
                   3*Q^2)) + (5*mD^3 - 2*A0*(mD - Q)^2 + 16*mD^2*Q + mD*Q^2 + 
                 2*Q^3)*delta[A0, Q, mD])/(2*mD^2*(mD - Q)^4*delta[A0, Q, 
                mD]) + (6*Q*Log[q2])/(mD - Q)^4) + 
           ((A0*(A0 - 6*mD)*(mD - Q) + (-5*mD + Q)*delta[A0, mD, mD])*
             phi[A0, mD, mD])/(mD^2*(mD - Q)^4) + 
           (((mD - Q)^2*(2*A0^4 - 8*A0^3*(mD + Q) - 2*A0*(mD - Q)*
                 (3*mD^2 - 4*Q^2) - (mD - Q)^2*(3*mD^2 + 4*mD*Q - 2*Q^2) + 
                A0^2*(11*mD^2 + 8*mD*Q + 12*Q^2)) - (mD - Q)*(2*A0^2*
                 (8*mD - 3*Q) - 4*A0*(9*mD^2 + 5*mD*Q - 3*Q^2) + 
                (mD - Q)*(27*mD^2 - 22*mD*Q + 6*Q^2))*delta[A0, Q, mD] + 
              2*(18*mD^2 - 9*mD*Q + 2*Q^2)*delta[A0, Q, mD]^2)*
             phi[A0, Q, mD])/(2*mD^3*(mD - Q)^4*delta[A0, Q, mD]) + 
           ((A0*(2*A0 - 11*Q)*(-mD + Q) + (2*mD - 9*Q)*delta[A0, Q, Q])*
             phi[A0, Q, Q])/((mD - Q)^4*Q^2))) + 
       Xb^2*(Nc*Yb^2*((2*mD - 4*Q)/(mD^2*Q - mD*Q^2) + (2*Log[mD]^2)/
            (mD - Q)^2 + (3*Log[Q]^2)/(mD - Q)^2 + 
           Log[A0]*(-(((A0 + mD - Q)*(2*A0^2*Q + A0*(mD^2 - 6*mD*Q - 4*Q^2) - 
                  (mD - Q)*(mD^2 - 4*mD*Q + 2*Q^2)) + (mD^2 - 4*mD*Q + 
                  2*Q*(-A0 + Q))*delta[A0, Q, mD])/(mD^2*(mD - Q)*Q*
                delta[A0, Q, mD])) - Log[mD]/(mD - Q)^2 + 
             Log[Q]/(mD - Q)^2) + ((2*mD - 4*Q)*Log[q2])/(mD^2*Q - mD*Q^2) + 
           Log[Q]*((-2*(mD - Q)*(-A0^3 + mD^3 - 2*mD*Q^2 + Q^3 + 
                 A0^2*(4*mD + 3*Q) - A0*(6*mD^2 + 2*mD*Q + 3*Q^2)) - 2*
                (-mD^2 + A0*(mD - Q) + mD*Q + Q^2)*delta[A0, Q, mD])/
              (mD^2*(mD - Q)^2*delta[A0, Q, mD]) - (2*Log[q2])/(mD - Q)^2) + 
           Log[mD]*(((mD - Q)*(A0^2*(mD + 4*Q) - 2*A0*Q*(7*mD + 4*Q) - 
                 (mD - Q)*(mD^2 - 7*mD*Q + 4*Q^2)) - (mD^2 + mD*Q - 4*Q^2)*
                delta[A0, Q, mD])/(mD*(mD - Q)^2*Q*delta[A0, Q, mD]) - 
             (5*Log[Q])/(mD - Q)^2 + (2*Log[q2])/(mD - Q)^2) + 
           ((A0*(A0 - 6*mD) - delta[A0, mD, mD])*phi[A0, mD, mD])/
            (mD^2*(mD - Q)^2) + (((mD - Q)*(2*A0^4 - 8*A0^3*(mD + Q) - 
                2*A0*(mD - 2*Q)*(mD - Q)*(mD + 2*Q) + (mD - Q)^2*
                 (mD^2 - 4*mD*Q + 2*Q^2) + A0^2*(11*mD^2 + 8*mD*Q + 
                  12*Q^2)) + delta[A0, Q, mD]*(-3*(3*mD - 2*Q)*(mD - Q)^2 + 
                A0^2*(-9*mD + 6*Q) + A0*(23*mD^2 + 6*mD*Q - 12*Q^2) + 
                (7*mD - 4*Q)*delta[A0, Q, mD]))*phi[A0, Q, mD])/
            (mD^3*(mD - Q)^2*delta[A0, Q, mD]) + 
           ((A0*(2*A0 - 11*Q) - 2*delta[A0, Q, Q])*phi[A0, Q, Q])/
            ((mD - Q)^2*Q^2)) + Nc*((2*A0*mD - 4*A0*Q + 4*mD*Q)/
            (mD^2*Q - mD*Q^2) + (2*Log[mD]^2)/(mD - Q) + 
           ((A0 - 3*mD + Q)*Log[Q]^2)/(mD - Q)^2 + 
           Log[A0]*((-2*A0*(mD - 2*Q))/(mD*(mD - Q)*Q) - 
             ((A0 - mD + Q)*Log[mD])/(mD - Q)^2 + ((A0 - mD + Q)*Log[Q])/
              (mD - Q)^2) + ((2*A0*mD - 4*A0*Q + 2*mD*Q)*Log[q2])/
            (mD^2*Q - mD*Q^2) + Log[Q]*((-2*(A0 - 3*mD + Q))/(mD - Q)^2 - 
             (2*(A0 - 3*mD + 2*Q)*Log[q2])/(mD - Q)^2) + 
           Log[mD]*((2*(A0 - 4*mD + 2*Q))/(mD - Q)^2 + 
             ((-A0 + mD + Q)*Log[Q])/(mD - Q)^2 + (2*(A0 - 3*mD + 2*Q)*Log[
                q2])/(mD - Q)^2) + ((A0*(A0 - 6*mD) - delta[A0, mD, mD])*
             phi[A0, mD, mD])/(mD^2*(mD - Q)) + (delta[A0, Q, mD]*
             phi[A0, Q, mD])/(mD*(mD - Q)^2) + 
           ((A0*(2*A0 - 11*Q)*(-mD + Q) + (2*mD - 3*Q)*delta[A0, Q, Q])*
             phi[A0, Q, Q])/((mD - Q)^2*Q^2)))))/cbe^2 + 
    Xb^4*(Nc^2*((4*mD*Q*Log[mD]^2)/(mD - Q)^4 + (2*(mD + Q)*Log[Q])/
         (mD - Q)^3 + (4*mD*Q*Log[Q]^2)/(mD - Q)^4 + 
        Log[mD]*((-2*(mD + Q))/(mD - Q)^3 - (8*mD*Q*Log[Q])/(mD - Q)^4)) + 
      Nc*((2*(mD^2 - 27*mD*Q + 2*Q^2))/(mD*(mD - Q)^2*Q) - 
        (2*(mD + 3*Q)*Log[mD]^2)/(mD - Q)^3 + (4*(3*mD + Q)*Log[Q]^2)/
         (mD - Q)^3 + (2*(mD^2 - 15*mD*Q + 2*Q^2)*Log[q2])/
         (mD*(mD - Q)^2*Q) + Log[Q]*((4*(-3*mD^2 - 10*mD*Q + Q^2))/
           (mD*(mD - Q)^3) - (2*(7*mD + 5*Q)*Log[q2])/(mD - Q)^3) + 
        Log[mD]*((-2*mD^2 + 44*mD*Q + 6*Q^2)/((mD - Q)^3*Q) + 
          (2*(-5*mD + Q)*Log[Q])/(mD - Q)^3 + (2*(7*mD + 5*Q)*Log[q2])/
           (mD - Q)^3) + (6*PolyLog[2, 1 - mD/Q])/(mD - Q)^2)) + 
    Xb^6*(Nc^2*((-2*(mD + Q))/(mD - Q)^4 - (2*mD*Q*(mD + Q)*Log[mD]^2)/
         (mD - Q)^6 - ((mD^2 + 6*mD*Q + Q^2)*Log[Q])/(mD - Q)^5 - 
        (2*mD*Q*(mD + Q)*Log[Q]^2)/(mD - Q)^6 + 
        Log[mD]*((mD^2 + 6*mD*Q + Q^2)/(mD - Q)^5 + (4*mD*Q*(mD + Q)*Log[Q])/
           (mD - Q)^6)) + Nc*((-mD^2 + 11*mD*Q + 2*Q^2)/(mD*(mD - Q)^3*Q) - 
        (2*(3*mD + Q)*Log[mD]^2)/(mD - Q)^4 - (2*(3*mD + 4*Q)*Log[Q]^2)/
         (mD - Q)^4 + ((-mD^2 + 5*mD*Q + 2*Q^2)*Log[q2])/(mD*(mD - Q)^3*Q) + 
        Log[mD]*(((mD - 5*Q)*(mD + 2*Q))/((mD - Q)^4*Q) + 
          (2*(6*mD + 5*Q)*Log[Q])/(mD - Q)^4 - (6*Q*Log[q2])/(mD - Q)^4) + 
        Log[Q]*((-3*mD^2 + 13*mD*Q + 2*Q^2)/(mD*(mD - Q)^4) + 
          (6*Q*Log[q2])/(mD - Q)^4) + (6*PolyLog[2, 1 - mD/Q])/
         (-mD + Q)^3)) + Xb^2*(Nc^2*((2*mD*Q*Log[mD]^2)/(mD - Q)^3 - 
        (2*mD*Q*Log[Q]^2)/(mD - Q)^3 + (2*(mD + Q)*Log[q2])/(mD - Q)^2 + 
        Log[mD]*(-((mD + Q)/(mD - Q)^2) - (4*mD*Q*Log[q2])/(mD - Q)^3) + 
        Log[Q]*(-((mD + Q)/(mD - Q)^2) + (4*mD*Q*Log[q2])/(mD - Q)^3)) + 
      Nc*(-2/mD + 4/(mD - Q) - Q^(-1) + (6*Log[mD]^2)/(mD - Q) - 
        (2*(6*mD - 5*Q)*Log[Q]^2)/(mD - Q)^2 + (-2/mD + 2/(mD - Q) - Q^(-1))*
         Log[q2] + Log[mD]*((mD^2 - 21*mD*Q + 16*Q^2)/((mD - Q)^2*Q) + 
          ((6*mD - 4*Q)*Log[Q])/(mD - Q)^2 - (2*(9*mD - 8*Q)*Log[q2])/
           (mD - Q)^2) + Log[Q]*((19*mD^2 - 17*mD*Q + 2*Q^2)/
           (mD*(mD - Q)^2) + (2*(9*mD - 8*Q)*Log[q2])/(mD - Q)^2) + 
        (6*PolyLog[2, 1 - mD/Q])/(-mD + Q))) + 
    (Nc*(-5/2 + mD*((mD - mu2)^(-1) - Q^(-1)) + (4*mu2 - 2*Q)/mD + 
        2*mu2*(Q^(-1) + (-mu2 + Q)^(-1)) + ((mD^2 - 2*mD*mu2 - mu2^2)*
          Log[mD]^2)/(mD - mu2)^2 + (1 - (4*mu2^2)/(mu2 - Q)^2)*Log[Q]^2 + 
        (2 + (4*mu2)/mD + (2*mD)/(-mD + mu2) + (4*mu2)/(mu2 - Q) - 
          (mD - 2*mu2)/Q - (2*Q)/mD)*Log[q2] + 4*Log[q2]^2 + 
        Log[mD]*((-mD^2 + 2*mD*mu2 + 2*mu2^2)/(mD - mu2)^2 + mD/Q + 
          (2*mu2^2*Log[mu2])/(mD - mu2)^2 + 2*Log[Q] - 
          (2*(2*mD^2 - 4*mD*mu2 + mu2^2)*Log[q2])/(mD - mu2)^2) + 
        Log[mu2]*(mu2*((-2*mD^2 + 3*mD*mu2 - 4*mu2^2)/(mD*(mD - mu2)^2) - 
            (6*mu2)/(mu2 - Q)^2 - 2/Q) + (4*mu2^2*Log[Q])/(mu2 - Q)^2 + 
          2*mu2^2*(-(mD - mu2)^(-2) - 2/(mu2 - Q)^2)*Log[q2]) + 
        Log[Q]*((2*Q)/mD + (3*mu2^2 + 2*mu2*Q + Q^2)/(mu2 - Q)^2 + 
          (-4 + (4*mu2^2)/(mu2 - Q)^2)*Log[q2]) + 
        (2*(mD^2 - 2*mD*mu2 - mu2^2)*PolyLog[2, 1 - mu2/mD])/(mD - mu2)^2 + 
        (2 - (8*mu2^2)/(mu2 - Q)^2)*PolyLog[2, 1 - mu2/Q]) + 
      Xb^2*(Nc^2*((4*mD*Log[mD]^2)/(mD - Q)^2 + (4*Q*Log[Q]^2)/(mD - Q)^2 + 
          Log[Q]*(4/(mD - Q) + (4*Log[q2])/(mD - Q)) + 
          Log[mD]*(4/(-mD + Q) - (4*(mD + Q)*Log[Q])/(mD - Q)^2 + 
            (4*Log[q2])/(-mD + Q))) + 
        Nc*((2*mD^2 - 4*mD*mu2 - 6*mD*Q + 8*mu2*Q - 4*Q^2)/
           (mD^2*Q - mD*Q^2) + (2*(2*mD + mu2 - Q)*Log[mD]^2)/(mD - Q)^2 + 
          4*mu2*(1/((mD - mu2)*(mu2 - Q)) + (Q^(-1) + (-mD + Q)^(-1))/mD)*
           Log[mu2] + ((-2*mu2 + 8*Q)*Log[Q]^2)/(mD - Q)^2 + 
          (2*(mD^2 + 2*(2*mu2 - Q)*Q - 2*mD*(mu2 + Q))*Log[q2])/
           (mD*(mD - Q)*Q) + Log[Q]*((-2*mD*mu2*(mD + 2*mu2) + 
              2*mD*(3*mD + 5*mu2)*Q + 2*(-5*mD + 2*mu2)*Q^2 - 4*Q^3)/
             (mD*(mD - Q)^2*(-mu2 + Q)) + ((4*(mD + mu2) - 10*Q)*Log[q2])/
             (mD - Q)^2) + Log[mD]*((-4*(mD^2 + mD*(mu2 - 3*Q) - mu2*
                (mu2 - 2*Q)))/((mD - mu2)*(mD - Q)^2) - 2/Q - 
            (2*(2*mD + 3*Q)*Log[Q])/(mD - Q)^2 + ((-4*(mD + mu2) + 10*Q)*
              Log[q2])/(mD - Q)^2) + (4*(mu2 - Q)*PolyLog[2, 1 - mu2/mD])/
           (mD - Q)^2 + (4*(-mu2 + Q)*PolyLog[2, 1 - mu2/Q])/(mD - Q)^2)) + 
      Xb^4*(Nc^2*(-8/(mD - Q)^2 - (4*mD*(mD + Q)*Log[mD]^2)/(mD - Q)^4 - 
          (4*Q*(mD + Q)*Log[Q]^2)/(mD - Q)^4 - (8*Log[q2])/(mD - Q)^2 + 
          Log[Q]*((-4*(mD + 3*Q))/(mD - Q)^3 - (4*(mD + Q)*Log[q2])/
             (mD - Q)^3) + Log[mD]*((4*(3*mD + Q))/(mD - Q)^3 + 
            (4*(mD + Q)^2*Log[Q])/(mD - Q)^4 + (4*(mD + Q)*Log[q2])/
             (mD - Q)^3)) + Nc*((-mD^3 + 2*mD^2*(mu2 - 5*Q) + 
            2*Q^2*(-2*mu2 + Q) + mD*Q*(-16*mu2 + 27*Q))/(mD*(mD - Q)^3*Q) + 
          ((-5*mD^2 + 3*mu2^2 - 2*(mD + 3*mu2)*Q + 2*Q^2)*Log[mD]^2)/
           (mD - Q)^4 + ((mD^2 - 3*mu2^2 - 8*mD*Q + 6*mu2*Q - 10*Q^2)*
            Log[Q]^2)/(mD - Q)^4 + Log[mu2]*((-2*mu2*(mD^2 - 2*mD*Q - 2*Q^2))/
             (mD*(mD - Q)^3*Q) + (6*mu2^2*Log[Q])/(mD - Q)^4) + 
          ((-mD^3 + 2*mD^2*(mu2 - 4*Q) + 2*Q^2*(-2*mu2 + Q) + 
             mD*Q*(-10*mu2 + 19*Q))*Log[q2])/(mD*(mD - Q)^3*Q) + 
          Log[Q]*((-7*mD^3 - 6*mD*(mD + 3*mu2)*Q + 29*mD*Q^2 + 2*Q^3)/
             (mD*(mD - Q)^4) - (6*(mD^2 - mD*Q + 2*(mu2 - Q)*Q)*Log[q2])/
             (mD - Q)^4) + Log[mD]*((mD^3 + 14*mD^2*Q + 3*mD*(2*mu2 - 7*Q)*
               Q + 12*(mu2 - Q)*Q^2)/((mD - Q)^4*Q) - (6*mu2^2*Log[mu2])/
             (mD - Q)^4 + (2*(2*mD^2 + 5*mD*Q + 4*Q^2)*Log[Q])/(mD - Q)^4 + 
            (6*(mD^2 - mD*Q + 2*(mu2 - Q)*Q)*Log[q2])/(mD - Q)^4) + 
          ((-2*mD^2 + 6*mu2^2 + 4*mD*Q + 4*Q*(-3*mu2 + Q))*
            PolyLog[2, 1 - mu2/mD])/(mD - Q)^4 + 
          (2*(mD^2 - 3*mu2^2 - 2*mD*Q + 6*mu2*Q - 2*Q^2)*
            PolyLog[2, 1 - mu2/Q])/(mD - Q)^4)))/cbe^2) + 
  ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2*gt^4*(Nc*Xb*Yb*((2*Log[mD]*Log[Q])/(mD - Q) + (2*Log[Q]^2)/(-mD + Q) + 
      Log[A0]*((2*Log[mD])/(-mD + Q) + (2*Log[Q])/(mD - Q)) - 
      (2*(A0 + mD - Q)*phi[A0, Q, mD])/(mD*(mD - Q)) + 
      (2*A0*phi[A0, Q, Q])/((mD - Q)*Q)) + 
    Nc*(-1 + (4*Pi^2)/3 + 4*Log[A0]^2 + Log[Q]^2 - 2*Log[q2] - Log[q2]^2 + 
      Log[Q]*(2 + 2*Log[q2]) + Log[A0]*(-2*Log[mD] - 4*Log[Q] - 2*Log[U]) + 
      2*Log[mD]*Log[U] - (2*A0*phi[A0, Q, Q])/Q - 
      (2*(A0 + mD - U)*phi[A0, U, mD])/mD) + 
    (sbe^2*(Nc*Yb^2*(-Q^(-1) + ((A0^2 - (mD - Q)^2 + delta[A0, Q, mD])*
           Log[A0])/(2*Q*delta[A0, Q, mD]) + 
         ((-A0^2 + mD^2 + 2*A0*Q - Q^2 + delta[A0, Q, mD])*Log[mD])/
          (2*Q*delta[A0, Q, mD]) - ((A0 + mD - Q)*Log[Q])/delta[A0, Q, mD] - 
         Log[q2]/Q + ((-(A0 + mD - Q)^2 + delta[A0, Q, mD])*phi[A0, Q, mD])/
          (2*mD*delta[A0, Q, mD])) + Nc*(1/2 + Pi^2/3 - A0/Q + Log[A0]^2 + 
         Log[A0]*((A0 + Q)/Q - Log[mD] - Log[Q]) + Log[Q]*(1/2 + Log[q2]) + 
         Log[q2]*(-(A0/Q) - Log[U]) - (3*Log[U])/2 + Log[mD]*Log[U] - 
         ((A0 + mD - U)*phi[A0, U, mD])/mD) + 
       Nc*Xt*Yb*(Log[Q]*(4/(-Q + U) + (4*Log[q2])/(-Q + U)) + 
         (4*Log[U])/(Q - U) + (4*Log[q2]*Log[U])/(Q - U) + 
         Log[A0]*((2*Log[Q])/(Q - U) + (2*Log[U])/(-Q + U)) + 
         Log[mD]*((2*Log[Q])/(Q - U) + (2*Log[U])/(-Q + U)) + 
         (2*(A0 + mD - Q)*phi[A0, Q, mD])/(mD*(-Q + U)) + 
         (2*(A0 + mD - U)*phi[A0, U, mD])/(mD*(Q - U))) + 
       Nc*Xt^3*Yb*(-16/(Q - U)^2 + Log[Q]*((4*(3*Q + U))/(Q - U)^3 + 
           (4*(Q + U)*Log[q2])/(Q - U)^3) - (4*(Q + 3*U)*Log[U])/(Q - U)^3 + 
         Log[q2]*(-8/(Q - U)^2 - (4*(Q + U)*Log[U])/(Q - U)^3) + 
         Log[A0]*((-2*(2*A0 - 2*mD + Q + U)*Log[Q])/(Q - U)^3 + 
           (2*(2*A0 - 2*mD + Q + U)*Log[U])/(Q - U)^3) + 
         Log[mD]*((-2*(-2*A0 + 2*mD + Q + U)*Log[Q])/(Q - U)^3 + 
           (2*(-2*A0 + 2*mD + Q + U)*Log[U])/(Q - U)^3) - 
         (2*((A0 + mD - Q)*(Q - U) + 2*delta[A0, Q, mD])*phi[A0, Q, mD])/
          (mD*(Q - U)^3) + ((2*(A0 + mD - U)*(-Q + U) + 4*delta[A0, U, mD])*
           phi[A0, U, mD])/(mD*(Q - U)^3)) + 
       Xt^4*(Nc*((Q*(Q - U) + A0*(5*Q + U))/(Q*(Q - U)^3) + 
           Log[Q]*((-(Q*(4*A0 + Q)) - 8*A0*U + U^2)/(2*(Q - U)^4) + 
             ((-(Q*(2*A0 + Q)) - 4*A0*U + U^2)*Log[q2])/(Q - U)^4) + 
           ((4*A0*Q + Q^2 + 8*A0*U - U^2)*Log[U])/(2*(Q - U)^4) + 
           Log[q2]*((2*Q*(Q - U) + A0*(5*Q + U))/(Q*(Q - U)^3) + 
             ((2*A0*Q + Q^2 + 4*A0*U - U^2)*Log[U])/(Q - U)^4) + 
           Log[A0]*((2*Q*(-Q + U) - A0*(5*Q + U))/(Q*(Q - U)^3) + 
             ((2*A0*Q + Q^2 + 4*A0*U - U^2)*Log[Q])/(Q - U)^4 + 
             ((-(Q*(2*A0 + Q)) - 4*A0*U + U^2)*Log[U])/(Q - U)^4)) + 
         Nc*Yb^2*((11*Q + U)/(Q*(Q - U)^3) + Log[Q]*
            ((-4*(2*Q + U))/(Q - U)^4 - (A0 + mD - Q)/((Q - U)^2*delta[A0, Q, 
                mD]) - (2*(Q + 2*U)*Log[q2])/(Q - U)^4) + 
           (2*(Q + 5*U)*Log[U])/(Q - U)^4 + Log[q2]*
            ((5*Q + U)/(Q*(Q - U)^3) + (2*(Q + 2*U)*Log[U])/(Q - U)^4) + 
           Log[A0]*((A0^2 - (mD - Q)^2 + delta[A0, Q, mD])/(2*Q*(Q - U)^
                2*delta[A0, Q, mD]) + ((3*A0 - 3*mD + Q + 2*U)*Log[Q])/
              (Q - U)^4 - ((3*A0 - 3*mD + Q + 2*U)*Log[U])/(Q - U)^4) + 
           Log[mD]*((-A0^2 + mD^2 + 2*A0*Q - Q^2 + delta[A0, Q, mD])/
              (2*Q*(Q - U)^2*delta[A0, Q, mD]) + ((-3*A0 + 3*mD + Q + 2*U)*
               Log[Q])/(Q - U)^4 - ((-3*A0 + 3*mD + Q + 2*U)*Log[U])/
              (Q - U)^4) + ((-((A0 + mD - Q)^2*(Q - U)^2) + (4*A0 + 4*mD - 
                3*Q - U)*(Q - U)*delta[A0, Q, mD] + 6*delta[A0, Q, mD]^2)*
             phi[A0, Q, mD])/(2*mD*(Q - U)^4*delta[A0, Q, mD]) + 
           (((A0 + mD - U)*(Q - U) - 3*delta[A0, U, mD])*phi[A0, U, mD])/
            (mD*(Q - U)^4))) + Xt^2*(Nc*Yb^2*(-2/(Q^2 - Q*U) + 
           Log[Q]*((2*((A0 + mD - Q)*(-Q + U) + delta[A0, Q, mD]))/
              ((Q - U)^2*delta[A0, Q, mD]) + (2*Log[q2])/(Q - U)^2) - 
           (2*Log[U])/(Q - U)^2 + Log[q2]*(-2/(Q^2 - Q*U) - 
             (2*Log[U])/(Q - U)^2) + Log[A0]*((A0^2 - (mD - Q)^2 + delta[A0, 
                Q, mD])/(Q*(Q - U)*delta[A0, Q, mD]) - Log[Q]/(Q - U)^2 + 
             Log[U]/(Q - U)^2) + Log[mD]*((-A0^2 + mD^2 + 2*A0*Q - Q^2 + 
               delta[A0, Q, mD])/(Q*(Q - U)*delta[A0, Q, mD]) - 
             Log[Q]/(Q - U)^2 + Log[U]/(Q - U)^2) + 
           ((-((A0 + mD - Q)^2*(Q - U)) + (A0 + mD - U)*delta[A0, Q, mD])*
             phi[A0, Q, mD])/(mD*(Q - U)^2*delta[A0, Q, mD]) - 
           ((A0 + mD - U)*phi[A0, U, mD])/(mD*(Q - U)^2)) + 
         Nc*((-2*A0 + 4*Q)/(Q*(Q - U)) + Log[Q]*
            (-((-2*A0 + 3*Q + U)/(Q - U)^2) + (2*(A0 - U)*Log[q2])/
              (Q - U)^2) + ((-2*A0 + Q + 3*U)*Log[U])/(Q - U)^2 + 
           Log[mD]*(((-A0 + mD + Q)*Log[Q])/(Q - U)^2 - 
             ((-A0 + mD + Q)*Log[U])/(Q - U)^2) + 
           Log[A0]*((2*A0)/(Q^2 - Q*U) - ((A0 + mD + Q - 2*U)*Log[Q])/
              (Q - U)^2 + ((A0 + mD + Q - 2*U)*Log[U])/(Q - U)^2) + 
           Log[q2]*((2*(-A0 + Q))/(Q*(Q - U)) - (2*(A0 - U)*Log[U])/
              (Q - U)^2) + (delta[A0, Q, mD]*phi[A0, Q, mD])/(mD*(Q - U)^2) + 
           (((A0 + mD - U)*(Q - U) - delta[A0, U, mD])*phi[A0, U, mD])/
            (mD*(Q - U)^2)))))/cbe^2 + 
    Nc*Xb*Yt*(Log[A0]*((2*Log[mD])/(-mD + Q) + (2*Log[Q])/(mD - Q)) + 
      (2*Log[mD]*Log[U])/(mD - Q) + (2*Log[Q]*Log[U])/(-mD + Q) - 
      (2*(A0 + mD - U)*phi[A0, U, mD])/(mD*(mD - Q)) + 
      (2*(A0 + Q - U)*phi[A0, U, Q])/((mD - Q)*Q)) + 
    Xt*(Nc*Yb*(Log[Q]*(4/(-Q + U) + (4*Log[q2])/(-Q + U)) + 
        (4*Log[U])/(Q - U) + (4*Log[q2]*Log[U])/(Q - U) + 
        Log[A0]*((2*Log[Q])/(Q - U) + (2*Log[U])/(-Q + U)) + 
        Log[mD]*((2*Log[Q])/(Q - U) + (2*Log[U])/(-Q + U)) + 
        (2*(A0 + mD - Q)*phi[A0, Q, mD])/(mD*(-Q + U)) + 
        (2*(A0 + mD - U)*phi[A0, U, mD])/(mD*(Q - U))) + 
      Nc*Xb*((2*(A0 - 2*Q)*Log[Q]^2)/((mD - Q)*(Q - U)) + 
        (2*(-A0 + Q + U)*Log[Q]*Log[U])/((mD - Q)*(Q - U)) + 
        Log[A0]*((2*Log[mD])/(-mD + Q) + 2*((mD - Q)^(-1) + (-Q + U)^(-1))*
           Log[Q] + (2*Log[U])/(Q - U)) + 
        Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/((mD - Q)*(Q - U)) - 
          (2*(-A0 + mD + U)*Log[U])/((mD - Q)*(Q - U))) + 
        (2*delta[A0, Q, mD]*phi[A0, Q, mD])/((mD^2 - mD*Q)*(Q - U)) - 
        (2*delta[A0, Q, Q]*phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)) - 
        (2*delta[A0, U, mD]*phi[A0, U, mD])/((mD^2 - mD*Q)*(Q - U)) + 
        (2*delta[A0, U, Q]*phi[A0, U, Q])/((mD - Q)*Q*(Q - U))) + 
      Yt*(Nc*((2*Log[Q]^2)/(Q - U) + (4*Log[U])/(Q - U) + 
          (4*Log[q2]*Log[U])/(Q - U) + Log[A0]*((2*Log[Q])/(Q - U) + 
            (2*Log[U])/(-Q + U)) + Log[Q]*(4/(-Q + U) + (4*Log[q2])/
             (-Q + U) + (2*Log[U])/(-Q + U)) - (2*A0*phi[A0, Q, Q])/
           (Q^2 - Q*U) + (2*(A0 + Q - U)*phi[A0, U, Q])/(Q*(Q - U))) + 
        Nc*Xb*Yb*((2*Log[Q]^2)/((mD - Q)*(-Q + U)) + (2*Log[Q]*Log[U])/
           ((mD - Q)*(Q - U)) + Log[mD]*((2*Log[Q])/((mD - Q)*(Q - U)) + 
            (2*Log[U])/((mD - Q)*(-Q + U))) - 
          (2*(A0 + mD - Q)*phi[A0, Q, mD])/((mD^2 - mD*Q)*(Q - U)) + 
          (2*A0*phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)) + 
          (2*(A0 + mD - U)*phi[A0, U, mD])/(mD*(mD - Q)*(Q - U)) + 
          (2*(A0 + Q - U)*phi[A0, U, Q])/(Q*(-mD + Q)*(Q - U))))) + 
    Xt^3*(Nc*Yb*(-16/(Q - U)^2 + Log[Q]*((4*(3*Q + U))/(Q - U)^3 + 
          (4*(Q + U)*Log[q2])/(Q - U)^3) - (4*(Q + 3*U)*Log[U])/(Q - U)^3 + 
        Log[q2]*(-8/(Q - U)^2 - (4*(Q + U)*Log[U])/(Q - U)^3) + 
        Log[A0]*((-2*(2*A0 - 2*mD + Q + U)*Log[Q])/(Q - U)^3 + 
          (2*(2*A0 - 2*mD + Q + U)*Log[U])/(Q - U)^3) + 
        Log[mD]*((-2*(-2*A0 + 2*mD + Q + U)*Log[Q])/(Q - U)^3 + 
          (2*(-2*A0 + 2*mD + Q + U)*Log[U])/(Q - U)^3) - 
        (2*((A0 + mD - Q)*(Q - U) + 2*delta[A0, Q, mD])*phi[A0, Q, mD])/
         (mD*(Q - U)^3) + ((2*(A0 + mD - U)*(-Q + U) + 4*delta[A0, U, mD])*
          phi[A0, U, mD])/(mD*(Q - U)^3)) + 
      Yt*(Nc*Xb*Yb*((2*(-2*A0 + 3*Q + U)*Log[Q]^2)/((mD - Q)*(Q - U)^3) - 
          (2*(-2*A0 + 3*Q + U)*Log[Q]*Log[U])/((mD - Q)*(Q - U)^3) + 
          Log[A0]*((4*Log[Q])/(Q - U)^3 + (4*Log[U])/(-Q + U)^3) + 
          Log[mD]*((-2*(-2*A0 + 2*mD + Q + U)*Log[Q])/((mD - Q)*(Q - U)^3) + 
            (2*(-2*A0 + 2*mD + Q + U)*Log[U])/((mD - Q)*(Q - U)^3)) - 
          (2*((A0 + mD - Q)*(Q - U) + 2*delta[A0, Q, mD])*phi[A0, Q, mD])/
           (mD*(mD - Q)*(Q - U)^3) + (2*(A0*(Q - U) + 2*delta[A0, Q, Q])*
            phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)^3) + 
          ((2*(A0 + mD - U)*(-Q + U) + 4*delta[A0, U, mD])*phi[A0, U, mD])/
           (mD*(mD - Q)*(Q - U)^3) + (2*((Q - U)*(A0 + Q - U) - 
             2*delta[A0, U, Q])*phi[A0, U, Q])/((mD - Q)*Q*(Q - U)^3)) + 
        Nc*(-16/(Q - U)^2 - (2*(-2*A0 + 3*Q + U)*Log[Q]^2)/(Q - U)^3 - 
          (4*(Q + 3*U)*Log[U])/(Q - U)^3 + Log[A0]*
           ((2*(-2*A0 + Q - U)*Log[Q])/(Q - U)^3 + (2*(2*A0 - Q + U)*Log[U])/
             (Q - U)^3) + Log[q2]*(-8/(Q - U)^2 - (4*(Q + U)*Log[U])/
             (Q - U)^3) + Log[Q]*((4*(3*Q + U))/(Q - U)^3 + 
            (4*(Q + U)*Log[q2])/(Q - U)^3 + (2*(-2*A0 + 3*Q + U)*Log[U])/
             (Q - U)^3) - (2*(A0*(Q - U) + 2*delta[A0, Q, Q])*phi[A0, Q, Q])/
           (Q*(Q - U)^3) + ((-2*(Q - U)*(A0 + Q - U) + 4*delta[A0, U, Q])*
            phi[A0, U, Q])/(Q*(Q - U)^3)))) + 
    Nc*Xb^2*((mD - 5*Q)/(Q*(-mD + Q)) - (2*mD*Log[Q]^2)/(mD - Q)^2 + 
      ((mD - 3*Q)*Log[q2])/(Q*(-mD + Q)) + 
      Log[mD]*((mD*(mD - 5*Q))/((mD - Q)^2*Q) + (2*mD*Log[Q])/(mD - Q)^2 - 
        (2*mD*Log[q2])/(mD - Q)^2) + Log[Q]*((mD + 3*Q)/(mD - Q)^2 + 
        (2*mD*Log[q2])/(mD - Q)^2) + (6*PolyLog[2, 1 - mD/Q])/(-mD + Q)) + 
    (Nc*Xb*Xt^2*((8*Sqrt[mu2]*Q*sgn*Log[Q]^2)/((mu2 - Q)*(-mD + Q)*(Q - U)) + 
        (8*Sqrt[mu2]*Q*sgn*Log[Q]*Log[U])/((mu2 - Q)*(-mD + Q)*(-Q + U)) + 
        Log[mD]*((8*mD*Sqrt[mu2]*sgn*Log[Q])/((mD - mu2)*(mD - Q)*(-Q + U)) + 
          (8*mD*Sqrt[mu2]*sgn*Log[U])/((mD - mu2)*(mD - Q)*(Q - U))) + 
        Log[mu2]*((8*mu2^(3/2)*sgn*Log[Q])/((mD - mu2)*(mu2 - Q)*(Q - U)) + 
          (8*mu2^(3/2)*sgn*Log[U])/((-mD + mu2)*(mu2 - Q)*(Q - U)))) + 
      Nc*Xb*Xt^4*((-4*Sqrt[mu2]*Q*sgn*(Q + U)*Log[Q]^2)/((mu2 - Q)*(-mD + Q)*
          (Q - U)^3) + Log[mD]*((-8*mD*Sqrt[mu2]*sgn)/((mD - mu2)*(mD - Q)*
            (Q - U)^2) + (4*mD*Sqrt[mu2]*sgn*(Q + U)*Log[Q])/
           ((mD - mu2)*(mD - Q)*(Q - U)^3) - (4*mD*Sqrt[mu2]*sgn*(Q + U)*
            Log[U])/((mD - mu2)*(mD - Q)*(Q - U)^3)) + 
        Log[mu2]*((8*mu2^(3/2)*sgn)/((mD - mu2)*(mu2 - Q)*(Q - U)^2) - 
          (4*mu2^(3/2)*sgn*(Q + U)*Log[Q])/((mD - mu2)*(mu2 - Q)*(Q - U)^3) + 
          (4*mu2^(3/2)*sgn*(Q + U)*Log[U])/((mD - mu2)*(mu2 - Q)*
            (Q - U)^3)) + Log[Q]*((8*Sqrt[mu2]*Q*sgn)/((mu2 - Q)*(-mD + Q)*
            (Q - U)^2) + (4*Sqrt[mu2]*Q*sgn*(Q + U)*Log[U])/
           ((mu2 - Q)*(-mD + Q)*(Q - U)^3))) + 
      Nc*Xb*((4*Sqrt[mu2]*(mD + mu2)*sgn*Log[mD]^2)/((mD - mu2)*(mD - Q)) + 
        (4*mu2^(3/2)*sgn*Log[Q]^2)/((mD - Q)*(mu2 - Q)) + 
        (4*Sqrt[mu2]*Q*sgn*Log[Q]*Log[U])/((mu2 - Q)*(-mD + Q)) + 
        Log[mu2]*((-4*mu2^(3/2)*(mD - 2*mu2 + Q)*sgn*Log[Q])/
           ((mD - mu2)*(mD - Q)*(mu2 - Q)) + (4*mu2^(3/2)*sgn*Log[U])/
           ((mD - mu2)*(mu2 - Q))) + Log[mD]*((8*mu2^(3/2)*sgn*Log[mu2])/
           ((-mD + mu2)*(mD - Q)) + (4*mD*Sqrt[mu2]*sgn*Log[Q])/
           ((mD - mu2)*(-mD + Q)) + (4*mD*Sqrt[mu2]*sgn*Log[U])/
           ((mD - mu2)*(-mD + Q))) + (8*Sqrt[mu2]*(mD + mu2)*sgn*
          PolyLog[2, 1 - mu2/mD])/((mD - mu2)*(mD - Q)) - 
        (8*Sqrt[mu2]*(mu2 + Q)*sgn*PolyLog[2, 1 - mu2/Q])/
         ((mu2 - Q)*(-mD + Q))) + Nc*Xt*((4*Sqrt[mu2]*sgn*Log[Q]^2)/
         (-Q + U) + Log[Q]*((8*Sqrt[mu2]*sgn)/(Q - U) + 
          (8*Sqrt[mu2]*sgn*Log[q2])/(Q - U)) + (8*Sqrt[mu2]*sgn*Log[U])/
         (-Q + U) + (8*Sqrt[mu2]*sgn*Log[q2]*Log[U])/(-Q + U) + 
        (4*Sqrt[mu2]*sgn*Log[U]^2)/(Q - U) + 
        (8*Sqrt[mu2]*sgn*PolyLog[2, 1 - mu2/Q])/(-Q + U) + 
        (8*Sqrt[mu2]*sgn*PolyLog[2, 1 - mu2/U])/(Q - U)) + 
      Nc*Xt^3*((32*Sqrt[mu2]*sgn)/(Q - U)^2 + 
        (4*Sqrt[mu2]*sgn*(-2*mu2 + Q + U)*Log[Q]^2)/(Q - U)^3 + 
        Log[Q]*((-8*Sqrt[mu2]*sgn*(3*Q + U))/(Q - U)^3 - 
          (8*Sqrt[mu2]*sgn*(Q + U)*Log[q2])/(Q - U)^3) + 
        (8*Sqrt[mu2]*sgn*(Q + 3*U)*Log[U])/(Q - U)^3 - 
        (4*Sqrt[mu2]*sgn*(-2*mu2 + Q + U)*Log[U]^2)/(Q - U)^3 + 
        Log[mu2]*((16*mu2^(3/2)*sgn*Log[Q])/(Q - U)^3 + 
          (16*mu2^(3/2)*sgn*Log[U])/(-Q + U)^3) + 
        Log[q2]*((16*Sqrt[mu2]*sgn)/(Q - U)^2 + (8*Sqrt[mu2]*sgn*(Q + U)*
            Log[U])/(Q - U)^3) + (8*Sqrt[mu2]*sgn*(-2*mu2 + Q + U)*
          PolyLog[2, 1 - mu2/Q])/(Q - U)^3 - 
        (8*Sqrt[mu2]*sgn*(-2*mu2 + Q + U)*PolyLog[2, 1 - mu2/U])/(Q - U)^3))/
     (cbe*sbe) + (Nc*(-3/2 + mD/(mD - mu2) - (mD - 2*mu2)/Q + 
        ((mD^2 - 2*mD*mu2 - mu2^2)*Log[mD]^2)/(mD - mu2)^2 + Log[Q]^2 + 
        Log[Q]*(1/2 + mD/(-mD + mu2) - Log[q2]) + (3/2 + mD/(-mD + mu2))*
         Log[U] + Log[q2]*(-((mD - 2*mu2 + Q)/Q) + Log[U]) + 
        Log[mD]*(mD*((mD + 2*mu2)/(mD - mu2)^2 + Q^(-1)) + 
          (2*mu2^2*Log[mu2])/(mD - mu2)^2 - (mD*(mD - 2*mu2)*Log[Q])/
           (mD - mu2)^2 - (mD*(mD - 2*mu2)*Log[U])/(mD - mu2)^2) + 
        Log[mu2]*(-((mu2*(2*mD + mu2))/(mD - mu2)^2) - (2*mu2)/Q - 
          (mu2^2*Log[Q])/(mD - mu2)^2 - (mu2^2*Log[U])/(mD - mu2)^2) + 
        (2*(mD^2 - 2*mD*mu2 - mu2^2)*PolyLog[2, 1 - mu2/mD])/(mD - mu2)^2 + 
        2*PolyLog[2, 1 - mu2/Q]) + Nc*Xt^2*((2*(mD - 2*mu2 + 3*Q))/
         (Q*(-Q + U)) + (2*(mu2 - U)*Log[Q]^2)/(Q - U)^2 + 
        Log[Q]*((2*(mD - 2*mu2 + 3*Q) + ((mD + mu2)*(-Q + U))/(mD - mu2))/
           (Q - U)^2 + (2*(mD - 2*mu2 + Q + U)*Log[q2])/(Q - U)^2) + 
        ((-2*(mD - 2*mu2 + 3*Q) + ((5*mD - 3*mu2)*(Q - U))/(mD - mu2))*
          Log[U])/(Q - U)^2 - (2*(mu2 - U)*Log[U]^2)/(Q - U)^2 + 
        Log[mD]*((2*mD)/(Q^2 - Q*U) + (2*mD*(-1 + ((mD - 2*mu2)*(-Q + U))/
              (mD - mu2)^2)*Log[Q])/(Q - U)^2 + 
          (2*(mD + (mD*(mD - 2*mu2)*(Q - U))/(mD - mu2)^2)*Log[U])/
           (Q - U)^2) + Log[mu2]*((-4*mu2)/(Q^2 - Q*U) + 
          (2*mu2^2*Log[Q])/((mD - mu2)^2*(-Q + U)) + (2*mu2^2*Log[U])/
           ((mD - mu2)^2*(Q - U))) + Log[q2]*((2*(mD - 2*mu2 + 2*Q))/
           (Q*(-Q + U)) - (2*(mD - 2*mu2 + Q + U)*Log[U])/(Q - U)^2) + 
        (4*(mu2 - U)*PolyLog[2, 1 - mu2/Q])/(Q - U)^2 - 
        (4*(mu2 - U)*PolyLog[2, 1 - mu2/U])/(Q - U)^2) + 
      Nc*Xt^4*((mD^2*(5*Q + U) + 2*mu2*(-3*Q*(Q + U) + mu2*(8*Q + U)) + 
          mD*(-3*mu2*(7*Q + U) + 4*Q*(Q + 2*U)))/((mD - mu2)*Q*(Q - U)^3) + 
        ((mu2 - U)*(3*mu2 - 2*Q - U)*Log[Q]^2)/(Q - U)^4 + 
        Log[Q]*((Q*(-4*mD^2 + 3*mD*(8*mu2 - Q) + 5*mu2*(-4*mu2 + Q)) - 
            4*(mD - mu2)*(2*mD - 4*mu2 + 5*Q)*U - (mD + mu2)*U^2)/
           (2*(mD - mu2)*(Q - U)^4) - ((Q*(2*mD - 4*mu2 + Q) + 
             4*(mD - 2*mu2 + Q)*U + U^2)*Log[q2])/(Q - U)^4) + 
        ((4*mD^2*(Q + 2*U) + mu2*(-3*(Q + U)*(Q + 3*U) + 4*mu2*(2*Q + 7*U)) + 
           mD*(-12*mu2*(Q + 3*U) + (Q + U)*(Q + 11*U)))*Log[U])/
         (2*(mD - mu2)*(Q - U)^4) + ((3*mu2 - 2*Q - U)*(-mu2 + U)*Log[U]^2)/
         (Q - U)^4 + Log[mu2]*((2*mu2*(-3 + (((mD - mu2)^2 + mu2*Q)*(Q - U))/
              ((mD - mu2)^2*Q)))/(-Q + U)^3 + 
          (mu2^2*(-6*(mD - mu2)^2 + (Q - U)*(Q + U))*Log[Q])/
           ((mD - mu2)^2*(Q - U)^4) + (mu2^2*(6*(mD - mu2)^2 - Q^2 + U^2)*
            Log[U])/((mD - mu2)^2*(Q - U)^4)) + 
        Log[q2]*((3*Q*(Q + U) + mD*(5*Q + U) - 2*mu2*(5*Q + U))/
           (Q*(Q - U)^3) + ((Q*(2*mD - 4*mu2 + Q) + 4*(mD - 2*mu2 + Q)*U + 
             U^2)*Log[U])/(Q - U)^4) + Log[mD]*
         ((mD*(6 - ((mD^2 - 2*mD*(mu2 + Q) + mu2*(mu2 + 4*Q))*(Q - U))/
              ((mD - mu2)^2*Q)))/(-Q + U)^3 + 
          (mD*(Q*(2*(mD - mu2)^2 + (mD - 2*mu2)*Q) + 4*(mD - mu2)^2*U - 
             (mD - 2*mu2)*U^2)*Log[Q])/((mD - mu2)^2*(Q - U)^4) + 
          (mD*(Q*(-2*(mD - mu2)^2 - (mD - 2*mu2)*Q) - 4*(mD - mu2)^2*U + 
             (mD - 2*mu2)*U^2)*Log[U])/((mD - mu2)^2*(Q - U)^4)) + 
        (2*(mu2 - U)*(3*mu2 - 2*Q - U)*PolyLog[2, 1 - mu2/Q])/(Q - U)^4 - 
        (2*(-mu2 + U)*(-3*mu2 + 2*Q + U)*PolyLog[2, 1 - mu2/U])/(Q - U)^4))/
     cbe^2 + (cbe^2*(Nc*Xb^2*(4/(mD - Q) - ((-A0 + mD + Q)*Log[Q]^2)/
           (mD - Q)^2 + Log[A0]*(((A0 + mD - Q)*Log[mD])/(mD - Q)^2 - 
            ((A0 + mD - Q)*Log[Q])/(mD - Q)^2) + (2*Log[q2])/(mD - Q) + 
          Log[mD]*((-4*mD)/(mD - Q)^2 + ((-A0 + mD + Q)*Log[Q])/(mD - Q)^2 - 
            (2*mD*Log[q2])/(mD - Q)^2) + Log[Q]*((2*(mD + Q))/(mD - Q)^2 + 
            (2*mD*Log[q2])/(mD - Q)^2) + (delta[A0, Q, mD]*phi[A0, Q, mD])/
           (mD*(mD - Q)^2) + ((A0*(A0 - 5*Q)*(-mD + Q) + (mD - 2*Q)*
              delta[A0, Q, Q])*phi[A0, Q, Q])/((mD - Q)^2*Q^2)) + 
        Nc*(Log[Q]^2 + Log[Q]*(-2 - 2*Log[q2]) + 2*Log[q2]^2 - 2*Log[U] + 
          Log[mD]*Log[U] - 2*Log[q2]*Log[U] + Log[A0]*(4 - Log[mD] + 
            Log[U]) + ((A0^2 - 5*A0*Q - delta[A0, Q, Q])*phi[A0, Q, Q])/Q^2 - 
          ((A0 + mD - U)*phi[A0, U, mD])/mD) + 
        Nc*Xb*Yt*(Log[A0]*((2*Log[mD])/(-mD + Q) + (2*Log[Q])/(mD - Q)) + 
          (2*Log[mD]*Log[U])/(mD - Q) + (2*Log[Q]*Log[U])/(-mD + Q) - 
          (2*(A0 + mD - U)*phi[A0, U, mD])/(mD*(mD - Q)) + 
          (2*(A0 + Q - U)*phi[A0, U, Q])/((mD - Q)*Q)) + 
        Yt^2*(Nc*((((A0 + Q - U)^2 - delta[A0, U, Q])*Log[A0])/
             (Q*delta[A0, U, Q]) - (2*(A0 + Q - U)*Log[Q])/delta[A0, U, Q] + 
            ((-A0^2 + Q^2 + 2*A0*U - U^2 + delta[A0, U, Q])*Log[U])/
             (Q*delta[A0, U, Q]) + ((A0 - U + ((A0 + Q - U)*(-A0^2 + 
                  (Q - U)*U + A0*(Q + 2*U)))/delta[A0, U, Q])*phi[A0, U, Q])/
             Q^2) + Nc*Xb^2*(Log[A0]*((-(A0 + Q - U)^2 + delta[A0, U, Q])/(
                (mD - Q)*Q*delta[A0, U, Q]) - Log[mD]/(mD - Q)^2 + 
              Log[Q]/(mD - Q)^2) + ((-A0^2 + Q^2 + 2*A0*U - U^2 + delta[A0, 
                U, Q])*Log[U])/(Q*(-mD + Q)*delta[A0, U, Q]) + 
            (Log[mD]*Log[U])/(mD - Q)^2 + Log[Q]*((2*(A0 + Q - U))/((mD - Q)*
                delta[A0, U, Q]) - Log[U]/(mD - Q)^2) - 
            ((A0 + mD - U)*phi[A0, U, mD])/(mD*(mD - Q)^2) + 
            (((mD - Q)*(A0 + Q - U)*(A0^2 + U*(-Q + U) - A0*(Q + 2*U)) + 
               (-(A0*mD) + 2*A0*Q + Q^2 + mD*U - 2*Q*U)*delta[A0, U, Q])*
              phi[A0, U, Q])/((mD - Q)^2*Q^2*delta[A0, U, Q]))) + 
        Xt*(Nc*Xb*((2*(A0 - 2*Q)*Log[Q]^2)/((mD - Q)*(Q - U)) + 
            (2*(-A0 + Q + U)*Log[Q]*Log[U])/((mD - Q)*(Q - U)) + 
            Log[A0]*((2*Log[mD])/(-mD + Q) + 2*((mD - Q)^(-1) + 
                (-Q + U)^(-1))*Log[Q] + (2*Log[U])/(Q - U)) + 
            Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/((mD - Q)*(Q - U)) - 
              (2*(-A0 + mD + U)*Log[U])/((mD - Q)*(Q - U))) + 
            (2*delta[A0, Q, mD]*phi[A0, Q, mD])/((mD^2 - mD*Q)*(Q - U)) - 
            (2*delta[A0, Q, Q]*phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)) - 
            (2*delta[A0, U, mD]*phi[A0, U, mD])/((mD^2 - mD*Q)*(Q - U)) + 
            (2*delta[A0, U, Q]*phi[A0, U, Q])/((mD - Q)*Q*(Q - U))) + 
          Yt*(Nc*((2*Log[Q]^2)/(Q - U) + (2*Log[Q]*Log[U])/(-Q + U) + 
              Log[A0]*((2*Log[Q])/(-Q + U) + (2*Log[U])/(Q - U)) - 
              (2*(-(A0*(A0 - 5*Q)) + delta[A0, Q, Q])*phi[A0, Q, Q])/(Q^2*
                (Q - U)) - (2*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2 - 
                 delta[A0, U, Q])*phi[A0, U, Q])/(Q^2*(Q - U))) + 
            Nc*Xb^2*((2*(A0 - mD - Q)*Log[Q]^2)/((mD - Q)^2*(Q - U)) + 
              Log[A0]*((-2*Log[mD])/(mD - Q)^2 + (2*Log[Q])/(mD - Q)^2) + 
              (2*(-A0 + mD + U)*Log[Q]*Log[U])/((mD - Q)^2*(Q - U)) + 
              Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/((mD - Q)^2*(Q - U)) - 
                (2*(-A0 + mD + U)*Log[U])/((mD - Q)^2*(Q - U))) + 
              (2*delta[A0, Q, mD]*phi[A0, Q, mD])/(mD*(mD - Q)^2*(Q - U)) + 
              (2*(A0*(A0 - 5*Q)*(-mD + Q) + (mD - 2*Q)*delta[A0, Q, Q])*
                phi[A0, Q, Q])/((mD - Q)^2*Q^2*(Q - U)) - (2*delta[A0, U, mD]*
                phi[A0, U, mD])/(mD*(mD - Q)^2*(Q - U)) + 
              (2*((mD - Q)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2) - 
                 (mD - 2*Q)*delta[A0, U, Q])*phi[A0, U, Q])/((mD - Q)^2*Q^2*
                (Q - U))))) + 
        Xt^2*(Nc*(4/(Q - U) + Log[Q]*((-4*Q)/(Q - U)^2 - (2*Q*Log[q2])/
               (Q - U)^2) + (2*(Q + U)*Log[U])/(Q - U)^2 + 
            Log[q2]*(2/(Q - U) + (2*Q*Log[U])/(Q - U)^2) + 
            Log[A0]*(((A0 - mD + Q)*Log[Q])/(Q - U)^2 - ((A0 - mD + Q)*
                Log[U])/(Q - U)^2) + Log[mD]*(((-A0 + mD + Q)*Log[Q])/(Q - U)^
                2 - ((-A0 + mD + Q)*Log[U])/(Q - U)^2) + 
            (delta[A0, Q, mD]*phi[A0, Q, mD])/(mD*(Q - U)^2) + 
            (((A0 + mD - U)*(Q - U) - delta[A0, U, mD])*phi[A0, U, mD])/
             (mD*(Q - U)^2)) + Nc*Xb*Yt*((2*(A0 - 2*Q)*Log[Q]^2)/
             ((mD - Q)*(Q - U)^2) + (2*(A0 - 2*Q)*Log[Q]*Log[U])/
             ((-mD + Q)*(Q - U)^2) + Log[A0]*((-2*Log[Q])/(Q - U)^2 + 
              (2*Log[U])/(Q - U)^2) + Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/(
                (mD - Q)*(Q - U)^2) + (2*(A0 - mD - Q)*Log[U])/((mD - Q)*
                (Q - U)^2)) + (2*delta[A0, Q, mD]*phi[A0, Q, mD])/
             (mD*(mD - Q)*(Q - U)^2) - (2*delta[A0, Q, Q]*phi[A0, Q, Q])/
             ((mD - Q)*Q*(Q - U)^2) + (2*((A0 + mD - U)*(Q - U) - delta[A0, 
                U, mD])*phi[A0, U, mD])/(mD*(mD - Q)*(Q - U)^2) + 
            ((-2*(Q - U)*(A0 + Q - U) + 2*delta[A0, U, Q])*phi[A0, U, Q])/
             ((mD - Q)*Q*(Q - U)^2)) + 
          Yt^2*(Nc*(Log[Q]^2/(Q - U)^2 + ((A0^2 - Q^2 - 2*A0*U + U^2 - 
                 delta[A0, U, Q])*Log[U])/(Q*(Q - U)*delta[A0, U, Q]) + 
              Log[Q]*((2*(A0 + Q - U))/((Q - U)*delta[A0, U, Q]) - 
                Log[U]/(Q - U)^2) + Log[A0]*((-(A0 + Q - U)^2 + delta[A0, U, 
                   Q])/(Q*(Q - U)*delta[A0, U, Q]) - Log[Q]/(Q - U)^2 + 
                Log[U]/(Q - U)^2) + ((A0*(A0 - 5*Q) - delta[A0, Q, Q])*
                phi[A0, Q, Q])/(Q^2*(Q - U)^2) + (((Q - U)*(A0 + Q - U)*
                  (A0^2 + U*(-Q + U) - A0*(Q + 2*U)) + delta[A0, U, Q]*
                  (-A0^2 - 2*(Q - U)^2 + A0*(2*Q + 3*U) + delta[A0, U, Q]))*
                phi[A0, U, Q])/(Q^2*(Q - U)^2*delta[A0, U, Q])) + 
            Nc*Xb^2*(((-(A0 + Q - U)^2 + delta[A0, U, Q])*Log[A0])/(Q*
                (-mD + Q)*(Q - U)*delta[A0, U, Q]) - ((-A0 + mD + Q)*
                Log[Q]^2)/((mD - Q)^2*(Q - U)^2) + ((A0^2 - Q^2 - 2*A0*U + 
                 U^2 - delta[A0, U, Q])*Log[U])/(Q*(-mD + Q)*(Q - U)*
                delta[A0, U, Q]) + Log[mD]*(((-A0 + mD + Q)*Log[Q])/
                 ((mD - Q)^2*(Q - U)^2) - ((-A0 + mD + Q)*Log[U])/
                 ((mD - Q)^2*(Q - U)^2)) + Log[Q]*((2*(A0 + Q - U))/
                 ((-mD + Q)*(Q - U)*delta[A0, U, Q]) + ((-A0 + mD + Q)*
                  Log[U])/((mD - Q)^2*(Q - U)^2)) + (delta[A0, Q, mD]*
                phi[A0, Q, mD])/(mD*(mD - Q)^2*(Q - U)^2) + 
              ((A0*(A0 - 5*Q)*(-mD + Q) + (mD - 2*Q)*delta[A0, Q, Q])*
                phi[A0, Q, Q])/((mD - Q)^2*Q^2*(Q - U)^2) + 
              (((A0 + mD - U)*(Q - U) - delta[A0, U, mD])*phi[A0, U, mD])/(mD*
                (mD - Q)^2*(Q - U)^2) + (((mD - Q)*(Q - U)*(A0 + Q - U)*
                  (-A0^2 + (Q - U)*U + A0*(Q + 2*U)) + delta[A0, U, Q]*
                  (A0^2*(mD - Q) + (2*mD - 3*Q)*(Q - U)^2 + A0*(-2*mD*Q + 
                     Q^2 - 3*mD*U + 4*Q*U) - (mD - 2*Q)*delta[A0, U, Q]))*
                phi[A0, U, Q])/((mD - Q)^2*Q^2*(Q - U)^2*delta[A0, U, 
                 Q]))))) + Nc*Xb^2*(4/(-mD + Q) - (2*mD*Log[mD]^2)/
         (mD - Q)^2 + (4*mu2*Log[mu2])/((mD - Q)*(mu2 - Q)) + 
        (2*mD*Log[Q]^2)/(mD - Q)^2 + (2*Log[q2])/(-mD + Q) + 
        (2*Log[U])/(-mD + Q) + Log[Q]*((4*(mD*mu2 - Q^2))/
           ((mD - Q)^2*(-mu2 + Q)) - (2*mD*Log[q2])/(mD - Q)^2 - 
          (2*mD*Log[U])/(mD - Q)^2) + Log[mD]*((4*mD)/(mD - Q)^2 + 
          (2*mD*Log[q2])/(mD - Q)^2 + (2*mD*Log[U])/(mD - Q)^2) - 
        (4*mD*PolyLog[2, 1 - mu2/mD])/(mD - Q)^2 + 
        (4*mD*PolyLog[2, 1 - mu2/Q])/(mD - Q)^2) + 
      Nc*(4*mu2*((-mu2 + Q)^(-1) + (-mu2 + U)^(-1))*Log[mu2] - 2*Log[Q]^2 - 
        2*Log[q2]^2 + (4*mu2*Log[U])/(mu2 - U) + 2*Log[q2]*Log[U] - 
        2*Log[U]^2 + Log[Q]*((4*mu2)/(mu2 - Q) + 2*Log[q2] + 2*Log[U]) - 
        4*PolyLog[2, 1 - mu2/Q] - 4*PolyLog[2, 1 - mu2/U]) + 
      Xt^2*(Nc*Xb^2*(2/((mD - Q)*(Q - U)) - (2*mD*Q*Log[Q]^2)/
           ((mD - Q)^2*(Q - U)^2) + (2*Q*Log[U])/((mD - Q)*(Q - U)^2) + 
          Log[mD]*((2*mD)/((mD - Q)^2*(-Q + U)) + (2*mD*Q*Log[Q])/
             ((mD - Q)^2*(Q - U)^2) - (2*mD*Q*Log[U])/((mD - Q)^2*
              (Q - U)^2)) + Log[Q]*((2*(Q^2 - mD*U))/((mD - Q)^2*(Q - U)^2) + 
            (2*mD*Q*Log[U])/((mD - Q)^2*(Q - U)^2))) + 
        Nc*(4/(-Q + U) + (4*mu2*Log[mu2])/((mu2 - U)*(Q - U)) + 
          ((-4*mu2*Q + 4*U^2)*Log[U])/((mu2 - U)*(Q - U)^2) + 
          (2*Q*Log[U]^2)/(Q - U)^2 + Log[q2]*(2/(-Q + U) - 
            (2*Q*Log[U])/(Q - U)^2) + Log[Q]*((2*(Q + U))/(Q - U)^2 + 
            (2*Q*Log[q2])/(Q - U)^2 - (2*Q*Log[U])/(Q - U)^2) - 
          (4*Q*PolyLog[2, 1 - mu2/Q])/(Q - U)^2 + (4*Q*PolyLog[2, 1 - mu2/U])/
           (Q - U)^2)))/sbe^2 + 
    Xt^4*(Nc*(24/(Q - U)^2 + ((7*Q + U)*Log[Q]^2)/(Q - U)^3 + 
        (6*(Q + 3*U)*Log[U])/(Q - U)^3 + (3*Log[U]^2)/(Q - U)^2 + 
        Log[Q]*((-2*(7*Q + 5*U))/(Q - U)^3 - (8*(Q + U)*Log[q2])/(Q - U)^3 + 
          (2*(-5*Q + U)*Log[U])/(Q - U)^3) + Log[A0]*(-8/(Q - U)^2 + 
          (4*(Q + U)*Log[Q])/(Q - U)^3 - (4*(Q + U)*Log[U])/(Q - U)^3) + 
        Log[q2]*(16/(Q - U)^2 + (8*(Q + U)*Log[U])/(Q - U)^3) + 
        (6*PolyLog[2, 1 - Q/U])/(Q - U)^2) + 
      Nc*Xb^2*((11*Q + U)/(Q*(Q - U)^3) + 
        ((-3*(mD + Q - 2*U) + (2*(Q - U)*(4*Q^2 + mD*(-5*Q + U)))/(mD - Q)^2)*
          Log[Q]^2)/(Q - U)^4 + (2*(U*(-7*Q + U) + mD*(Q + 5*U))*Log[U])/
         ((mD - Q)*(Q - U)^4) + (3*(mD + Q - 2*U)*Log[U]^2)/(Q - U)^4 + 
        Log[mD]*(mD/((mD - Q)*Q*(Q - U)^2) + (2*mD*(mD - U)*(3*mD - 2*Q - U)*
            Log[Q])/((mD - Q)^2*(Q - U)^4) - (2*mD*(mD - U)*(3*mD - 2*Q - U)*
            Log[U])/((mD - Q)^2*(Q - U)^4)) + 
        Log[q2]*((5*Q + U)/(Q*(Q - U)^3) + (2*(Q + 2*U)*Log[U])/(Q - U)^4) + 
        Log[Q]*((5*Q^2 + 10*Q*U - 3*U^2 - 4*mD*(2*Q + U))/
           ((mD - Q)*(Q - U)^4) - (2*(Q + 2*U)*Log[q2])/(Q - U)^4 + 
          (2*(5*mD*Q - 4*Q^2 - mD*U)*Log[U])/((mD - Q)^2*(Q - U)^3)) - 
        (6*(mD - U)^2*PolyLog[2, 1 - mD/Q])/((mD - Q)*(Q - U)^4) + 
        (6*(mD - U)^2*PolyLog[2, 1 - mD/U])/((mD - Q)*(Q - U)^4) - 
        (6*PolyLog[2, 1 - Q/U])/((mD - Q)*(Q - U)^2))) + 
    Xt^2*(Nc*Xb*Yb*((2*(-A0 + Q + U)*Log[Q]^2)/((mD - Q)*(Q - U)^2) - 
        (2*(-A0 + Q + U)*Log[Q]*Log[U])/((mD - Q)*(Q - U)^2) + 
        Log[A0]*((2*Log[Q])/(Q - U)^2 - (2*Log[U])/(Q - U)^2) + 
        Log[mD]*((-2*(-A0 + mD + U)*Log[Q])/((mD - Q)*(Q - U)^2) + 
          (2*(-A0 + mD + U)*Log[U])/((mD - Q)*(Q - U)^2)) - 
        (2*((A0 + mD - Q)*(Q - U) + delta[A0, Q, mD])*phi[A0, Q, mD])/
         (mD*(mD - Q)*(Q - U)^2) + (2*(A0*(Q - U) + delta[A0, Q, Q])*
          phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)^2) + 
        (2*delta[A0, U, mD]*phi[A0, U, mD])/(mD*(mD - Q)*(Q - U)^2) - 
        (2*delta[A0, U, Q]*phi[A0, U, Q])/((mD - Q)*Q*(Q - U)^2)) + 
      Nc*Xb*Yt*((2*(A0 - 2*Q)*Log[Q]^2)/((mD - Q)*(Q - U)^2) + 
        (2*(A0 - 2*Q)*Log[Q]*Log[U])/((-mD + Q)*(Q - U)^2) + 
        Log[A0]*((-2*Log[Q])/(Q - U)^2 + (2*Log[U])/(Q - U)^2) + 
        Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/((mD - Q)*(Q - U)^2) + 
          (2*(A0 - mD - Q)*Log[U])/((mD - Q)*(Q - U)^2)) + 
        (2*delta[A0, Q, mD]*phi[A0, Q, mD])/(mD*(mD - Q)*(Q - U)^2) - 
        (2*delta[A0, Q, Q]*phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)^2) + 
        (2*((A0 + mD - U)*(Q - U) - delta[A0, U, mD])*phi[A0, U, mD])/
         (mD*(mD - Q)*(Q - U)^2) + ((-2*(Q - U)*(A0 + Q - U) + 
           2*delta[A0, U, Q])*phi[A0, U, Q])/((mD - Q)*Q*(Q - U)^2)) + 
      Nc*(8/(Q - U) + (2*A0*Log[Q]^2)/(Q - U)^2 + ((-6*Q + 14*U)*Log[U])/
         (Q - U)^2 + (6*Log[U]^2)/(Q - U) + 
        Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/(Q - U)^2 + 
          (2*(A0 - mD - Q)*Log[U])/(Q - U)^2) + 
        Log[A0]*((-2*(mD + 2*Q - 3*U)*Log[Q])/(Q - U)^2 + 
          (2*(mD + 2*Q - 3*U)*Log[U])/(Q - U)^2) + 
        Log[Q]*((2*(Q - 5*U))/(Q - U)^2 + (4*(2*Q - 3*U)*Log[q2])/(Q - U)^2 - 
          (2*(A0 + 3*Q - 3*U)*Log[U])/(Q - U)^2) + 
        Log[q2]*(4/(Q - U) + ((-8*Q + 12*U)*Log[U])/(Q - U)^2) + 
        (2*delta[A0, Q, mD]*phi[A0, Q, mD])/(mD*(Q - U)^2) - 
        (2*(A0*(Q - U) + delta[A0, Q, Q])*phi[A0, Q, Q])/(Q*(Q - U)^2) + 
        (2*((A0 + mD - U)*(Q - U) - delta[A0, U, mD])*phi[A0, U, mD])/
         (mD*(Q - U)^2) + (2*delta[A0, U, Q]*phi[A0, U, Q])/(Q*(Q - U)^2) + 
        (12*PolyLog[2, 1 - Q/U])/(Q - U)) + 
      Nc*Xb^2*(-2/(Q^2 - Q*U) + ((-3*mD^2 - 2*mD*Q + Q^2 + 4*mD*U)*Log[Q]^2)/
         ((mD - Q)^2*(Q - U)^2) - (2*(mD - U)*Log[U])/((mD - Q)*(Q - U)^2) + 
        (3*Log[U]^2)/(Q - U)^2 + Log[q2]*(-2/(Q^2 - Q*U) - 
          (2*Log[U])/(Q - U)^2) + Log[mD]*((2*mD)/((mD - Q)*Q*(Q - U)) + 
          (4*mD*(mD - U)*Log[Q])/((mD - Q)^2*(Q - U)^2) + 
          (4*mD*(-mD + U)*Log[U])/((mD - Q)^2*(Q - U)^2)) + 
        Log[Q]*(2/(Q - U)^2 + (2*Log[q2])/(Q - U)^2 + 
          ((8*mD*Q - 4*Q^2 - 4*mD*U)*Log[U])/((mD - Q)^2*(Q - U)^2)) + 
        (6*(mD + Q - 2*U)*PolyLog[2, 1 - mD/Q])/((-mD + Q)*(Q - U)^2) + 
        (6*(mD - U)^2*PolyLog[2, 1 - mD/U])/((mD - Q)^2*(Q - U)^2) + 
        (6*(-2*mD + Q + U)*PolyLog[2, 1 - Q/U])/((mD - Q)^2*(Q - U))))) + 
  ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^4*gt^2*(Nc*Xb*Yb*((2*Log[Q]^2)/(-mD + Q) + 
      Log[A0]*((2*Log[mD])/(mD - Q) + (2*Log[Q])/(-mD + Q)) + 
      Log[Q]*(4/(mD - Q) + (4*Log[q2])/(mD - Q)) + 
      Log[mD]*(4/(-mD + Q) + (2*Log[Q])/(mD - Q) + (4*Log[q2])/(-mD + Q)) + 
      (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2 - delta[A0, Q, mD])*
        phi[A0, Q, mD])/(mD^2*(mD - Q)) + 
      (2*(-(A0*(A0 - 5*Q)) + delta[A0, Q, Q])*phi[A0, Q, Q])/
       ((mD - Q)*Q^2)) + Nc*Xb^3*Yb*(-16/(mD - Q)^2 + 
      (2*(-2*A0 + mD + 3*Q)*Log[Q]^2)/(mD - Q)^3 + 
      Log[A0]*((-2*(2*A0 + mD - Q)*Log[mD])/(mD - Q)^3 + 
        (2*(2*A0 + mD - Q)*Log[Q])/(mD - Q)^3) - (8*Log[q2])/(mD - Q)^2 + 
      Log[Q]*((-4*(mD + 3*Q))/(mD - Q)^3 - (4*(mD + Q)*Log[q2])/(mD - Q)^3) + 
      Log[mD]*((4*(3*mD + Q))/(mD - Q)^3 - (2*(-2*A0 + mD + 3*Q)*Log[Q])/
         (mD - Q)^3 + (4*(mD + Q)*Log[q2])/(mD - Q)^3) + 
      (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2) + 
         (-3*mD + Q)*delta[A0, Q, mD])*phi[A0, Q, mD])/(mD^2*(mD - Q)^3) + 
      (2*(A0*(A0 - 5*Q)*(mD - Q) - (mD - 3*Q)*delta[A0, Q, Q])*phi[A0, Q, Q])/
       ((mD - Q)^3*Q^2)) + Nc*(-1 + (4*Pi^2)/3 + 4*Log[A0]^2 + Log[Q]^2 + 
      Log[Q]*(-2 - 2*Log[q2]) + 6*Log[q2] + 7*Log[q2]^2 + 
      Log[A0]*(2*Log[mD] - 8*Log[q2] - 2*Log[U]) + 
      Log[mD]*(-4 - 4*Log[q2] + 2*Log[U]) - 
      (2*(-(A0*(A0 - 5*Q)) + delta[A0, Q, Q])*phi[A0, Q, Q])/Q^2 + 
      (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2 - delta[A0, U, mD])*
        phi[A0, U, mD])/mD^2) + 
    Yt*(Nc*Xb*(Log[A0]*((2*Log[mD])/(mD - Q) + (2*Log[Q])/(-mD + Q)) + 
        Log[mD]*(4/(-mD + Q) + (4*Log[q2])/(-mD + Q) + (2*Log[U])/(mD - Q)) + 
        Log[Q]*(4/(mD - Q) + (4*Log[q2])/(mD - Q) + (2*Log[U])/(-mD + Q)) + 
        (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2 - 
           delta[A0, U, mD])*phi[A0, U, mD])/(mD^2*(mD - Q)) - 
        (2*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2 - delta[A0, U, Q])*
          phi[A0, U, Q])/((mD - Q)*Q^2)) + 
      Nc*Xb^3*(-16/(mD - Q)^2 + Log[A0]*((-2*(2*A0 + mD + Q - 2*U)*Log[mD])/
           (mD - Q)^3 + (2*(2*A0 + mD + Q - 2*U)*Log[Q])/(mD - Q)^3) - 
        (8*Log[q2])/(mD - Q)^2 + Log[mD]*((4*(3*mD + Q))/(mD - Q)^3 + 
          (4*(mD + Q)*Log[q2])/(mD - Q)^3 - (2*(-2*A0 + mD + Q + 2*U)*Log[U])/
           (mD - Q)^3) + Log[Q]*((-4*(mD + 3*Q))/(mD - Q)^3 - 
          (4*(mD + Q)*Log[q2])/(mD - Q)^3 + (2*(-2*A0 + mD + Q + 2*U)*Log[U])/
           (mD - Q)^3) + (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 
             3*mD*U + U^2) + (-3*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
         (mD^2*(mD - Q)^3) + (2*((mD - Q)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 
             3*Q*U + U^2) - (mD - 3*Q)*delta[A0, U, Q])*phi[A0, U, Q])/
         ((mD - Q)^3*Q^2))) + 
    (sbe^2*(Nc*Yb^2*((((A0 + mD - Q)^2 - delta[A0, Q, mD])*Log[A0])/
          (mD*delta[A0, Q, mD]) - (2*(A0 + mD - Q)*Log[mD])/
          delta[A0, Q, mD] + ((-A0^2 + mD^2 + 2*A0*Q - Q^2 + 
            delta[A0, Q, mD])*Log[Q])/(mD*delta[A0, Q, mD]) + 
         ((A0 - Q + ((A0 + mD - Q)*(-A0^2 + (mD - Q)*Q + A0*(mD + 2*Q)))/
             delta[A0, Q, mD])*phi[A0, Q, mD])/mD^2) + 
       Nc*Xb*Yb*((2*Log[mD]*Log[Q])/(mD - Q) + (2*Log[Q]^2)/(-mD + Q) + 
         Log[A0]*((2*Log[mD])/(-mD + Q) + (2*Log[Q])/(mD - Q)) - 
         (2*(A0 + mD - Q)*phi[A0, Q, mD])/(mD*(mD - Q)) + 
         (2*A0*phi[A0, Q, Q])/((mD - Q)*Q)) + 
       Nc*(Log[Q]^2 + Log[Q]*(-2 - 2*Log[q2]) + 2*Log[q2]^2 + 
         Log[A0]*(4 + Log[mD] - Log[U]) + Log[mD]*(-2 - 2*Log[q2] + Log[U]) - 
         (A0*phi[A0, Q, Q])/Q + ((A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + 
            U^2 - delta[A0, U, mD])*phi[A0, U, mD])/mD^2) + 
       Xb^2*(Nc*Yb^2*(((-A0^2 + mD^2 + 2*A0*Q - Q^2 + delta[A0, Q, mD])*
             Log[Q])/(mD*(mD - Q)*delta[A0, Q, mD]) + Log[Q]^2/(mD - Q)^2 + 
           Log[mD]*((-2*(A0 + mD - Q))/((mD - Q)*delta[A0, Q, mD]) - 
             Log[Q]/(mD - Q)^2) + Log[A0]*(((A0 + mD - Q)^2 - delta[A0, Q, 
                mD])/(mD*(mD - Q)*delta[A0, Q, mD]) + Log[mD]/(mD - Q)^2 - 
             Log[Q]/(mD - Q)^2) + (((mD - Q)*(A0 + mD - Q)*(-A0^2 + 
                (mD - Q)*Q + A0*(mD + 2*Q)) + ((mD - Q)^2 + A0*(2*mD - Q))*
               delta[A0, Q, mD])*phi[A0, Q, mD])/(mD^2*(mD - Q)^2*
             delta[A0, Q, mD]) - (A0*phi[A0, Q, Q])/((mD - Q)^2*Q)) + 
         Nc*(4/(-mD + Q) + Log[A0]*(-(((A0 + Q - U)*Log[mD])/(mD - Q)^2) + 
             ((A0 + Q - U)*Log[Q])/(mD - Q)^2) + (2*Log[q2])/(-mD + Q) + 
           Log[mD]*((2*(mD + Q))/(mD - Q)^2 + (2*Q*Log[q2])/(mD - Q)^2 - 
             ((-A0 + Q + U)*Log[U])/(mD - Q)^2) + Log[Q]*((-4*Q)/(mD - Q)^2 - 
             (2*Q*Log[q2])/(mD - Q)^2 + ((-A0 + Q + U)*Log[U])/(mD - Q)^2) + 
           (((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2) + 
              (-2*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
            (mD^2*(mD - Q)^2) + (delta[A0, U, Q]*phi[A0, U, Q])/
            ((mD - Q)^2*Q))) + 
       Xt^2*(Nc*Yb^2*(((-A0^2 + mD^2 + 2*A0*Q - Q^2 + delta[A0, Q, mD])*
             Log[Q])/(mD*(Q - U)*delta[A0, Q, mD]) + 
           Log[A0]*(((A0 + mD - Q)^2 - delta[A0, Q, mD])/(mD*(Q - U)*delta[
                A0, Q, mD]) + Log[Q]/(Q - U)^2 - Log[U]/(Q - U)^2) + 
           Log[mD]*((-2*(A0 + mD - Q))/((Q - U)*delta[A0, Q, mD]) - 
             Log[Q]/(Q - U)^2 + Log[U]/(Q - U)^2) + 
           ((-((A0 + mD - Q)*(A0^2 + Q*(-mD + Q) - A0*(mD + 2*Q))*(Q - U)) + 
              delta[A0, Q, mD]*(-A0^2 - 2*mD^2 + 3*mD*Q + 3*A0*(mD + Q) - 
                A0*U + Q*(-2*Q + U) + delta[A0, Q, mD]))*phi[A0, Q, mD])/
            (mD^2*(Q - U)^2*delta[A0, Q, mD]) + 
           ((A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2 - 
              delta[A0, U, mD])*phi[A0, U, mD])/(mD^2*(Q - U)^2)) + 
         Nc*(4/(-Q + U) - ((-A0 + Q + U)*Log[Q]^2)/(Q - U)^2 - 
           (4*U*Log[U])/(Q - U)^2 + Log[q2]*(2/(-Q + U) - (2*U*Log[U])/
              (Q - U)^2) + Log[A0]*(-(((A0 - Q + U)*Log[Q])/(Q - U)^2) + 
             ((A0 - Q + U)*Log[U])/(Q - U)^2) + 
           Log[Q]*((2*(Q + U))/(Q - U)^2 + (2*U*Log[q2])/(Q - U)^2 + 
             ((-A0 + Q + U)*Log[U])/(Q - U)^2) - 
           ((A0*(Q - U) + delta[A0, Q, Q])*phi[A0, Q, Q])/(Q*(Q - U)^2) + 
           (delta[A0, U, Q]*phi[A0, U, Q])/(Q*(Q - U)^2)) + 
         Nc*Xb^2*Yb^2*((((A0 + mD - Q)^2 - delta[A0, Q, mD])*Log[A0])/
            (mD*(mD - Q)*(Q - U)*delta[A0, Q, mD]) - ((-A0 + Q + U)*Log[Q]^2)/
            ((mD - Q)^2*(Q - U)^2) + Log[mD]*((-2*(A0 + mD - Q))/
              ((mD - Q)*(Q - U)*delta[A0, Q, mD]) + ((-A0 + Q + U)*Log[Q])/
              ((mD - Q)^2*(Q - U)^2) - ((-A0 + Q + U)*Log[U])/
              ((mD - Q)^2*(Q - U)^2)) + Log[Q]*((-A0^2 + mD^2 + 2*A0*Q - Q^
                2 + delta[A0, Q, mD])/((mD^2 - mD*Q)*(Q - U)*delta[A0, Q, 
                mD]) + ((-A0 + Q + U)*Log[U])/((mD - Q)^2*(Q - U)^2)) + 
           (((mD - Q)*(A0 + mD - Q)*(-A0^2 + (mD - Q)*Q + A0*(mD + 2*Q))*(Q - 
                U) + delta[A0, Q, mD]*(A0^2*(-mD + Q) - (mD - Q)^2*
                 (2*mD - 2*Q + U) + A0*(3*mD^2 + mD*(Q - 2*U) + 
                  Q*(-3*Q + U)) + (2*mD - Q)*delta[A0, Q, mD]))*
             phi[A0, Q, mD])/(mD^2*(mD - Q)^2*(Q - U)^2*delta[A0, Q, mD]) - 
           ((A0*(Q - U) + delta[A0, Q, Q])*phi[A0, Q, Q])/((mD - Q)^2*Q*
             (Q - U)^2) + (((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 
                3*mD*U + U^2) + (-2*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
            (mD^2*(mD - Q)^2*(Q - U)^2) + (delta[A0, U, Q]*phi[A0, U, Q])/
            ((mD - Q)^2*Q*(Q - U)^2)) + Nc*Xb*Yb*((2*(-A0 + Q + U)*Log[Q]^2)/
            ((mD - Q)*(Q - U)^2) - (2*(-A0 + Q + U)*Log[Q]*Log[U])/
            ((mD - Q)*(Q - U)^2) + Log[A0]*((2*Log[Q])/(Q - U)^2 - 
             (2*Log[U])/(Q - U)^2) + Log[mD]*((-2*(-A0 + mD + U)*Log[Q])/
              ((mD - Q)*(Q - U)^2) + (2*(-A0 + mD + U)*Log[U])/
              ((mD - Q)*(Q - U)^2)) - (2*((A0 + mD - Q)*(Q - U) + 
              delta[A0, Q, mD])*phi[A0, Q, mD])/(mD*(mD - Q)*(Q - U)^2) + 
           (2*(A0*(Q - U) + delta[A0, Q, Q])*phi[A0, Q, Q])/
            ((mD - Q)*Q*(Q - U)^2) + (2*delta[A0, U, mD]*phi[A0, U, mD])/
            (mD*(mD - Q)*(Q - U)^2) - (2*delta[A0, U, Q]*phi[A0, U, Q])/
            ((mD - Q)*Q*(Q - U)^2))) + 
       Xt*(Nc*Yb*(Log[A0]*((2*Log[Q])/(-Q + U) + (2*Log[U])/(Q - U)) + 
           Log[mD]*((2*Log[Q])/(Q - U) + (2*Log[U])/(-Q + U)) + 
           (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2 - 
              delta[A0, Q, mD])*phi[A0, Q, mD])/(mD^2*(Q - U)) - 
           (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2 - 
              delta[A0, U, mD])*phi[A0, U, mD])/(mD^2*(Q - U))) + 
         Nc*Xb^2*Yb*(((-2*A0 + 4*Q)*Log[Q]^2)/((mD - Q)^2*(Q - U)) + 
           Log[A0]*((2*Log[mD])/(mD - Q)^2 - (2*Log[Q])/(mD - Q)^2) - 
           (2*(-A0 + Q + U)*Log[Q]*Log[U])/((mD - Q)^2*(Q - U)) + 
           Log[mD]*((2*(A0 - 2*Q)*Log[Q])/((mD - Q)^2*(Q - U)) + 
             (2*(-A0 + Q + U)*Log[U])/((mD - Q)^2*(Q - U))) + 
           (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2) + 
              (-2*mD + Q)*delta[A0, Q, mD])*phi[A0, Q, mD])/
            (mD^2*(mD - Q)^2*(Q - U)) + (2*delta[A0, Q, Q]*phi[A0, Q, Q])/
            ((mD - Q)^2*Q*(Q - U)) - (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 
                2*A0*U - 3*mD*U + U^2) + (-2*mD + Q)*delta[A0, U, mD])*
             phi[A0, U, mD])/(mD^2*(mD - Q)^2*(Q - U)) - 
           (2*delta[A0, U, Q]*phi[A0, U, Q])/((mD - Q)^2*Q*(Q - U))) + 
         Nc*Xb*((2*(A0 - 2*Q)*Log[Q]^2)/((mD - Q)*(Q - U)) + 
           (2*(-A0 + Q + U)*Log[Q]*Log[U])/((mD - Q)*(Q - U)) + 
           Log[A0]*((2*Log[mD])/(-mD + Q) + 2*((mD - Q)^(-1) + (-Q + U)^(-1))*
              Log[Q] + (2*Log[U])/(Q - U)) + Log[mD]*
            ((2*(-A0 + mD + Q)*Log[Q])/((mD - Q)*(Q - U)) - 
             (2*(-A0 + mD + U)*Log[U])/((mD - Q)*(Q - U))) + 
           (2*delta[A0, Q, mD]*phi[A0, Q, mD])/((mD^2 - mD*Q)*(Q - U)) - 
           (2*delta[A0, Q, Q]*phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)) - 
           (2*delta[A0, U, mD]*phi[A0, U, mD])/((mD^2 - mD*Q)*(Q - U)) + 
           (2*delta[A0, U, Q]*phi[A0, U, Q])/((mD - Q)*Q*(Q - U))))))/cbe^2 + 
    Xt*(Nc*Yb*(Log[A0]*((2*Log[Q])/(-Q + U) + (2*Log[U])/(Q - U)) + 
        Log[mD]*((2*Log[Q])/(Q - U) + (2*Log[U])/(-Q + U)) + 
        (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2 - 
           delta[A0, Q, mD])*phi[A0, Q, mD])/(mD^2*(Q - U)) - 
        (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2 - 
           delta[A0, U, mD])*phi[A0, U, mD])/(mD^2*(Q - U))) + 
      Nc*Xb^2*Yb*(((-2*A0 + 4*Q)*Log[Q]^2)/((mD - Q)^2*(Q - U)) + 
        Log[A0]*((2*Log[mD])/(mD - Q)^2 - (2*Log[Q])/(mD - Q)^2) - 
        (2*(-A0 + Q + U)*Log[Q]*Log[U])/((mD - Q)^2*(Q - U)) + 
        Log[mD]*((2*(A0 - 2*Q)*Log[Q])/((mD - Q)^2*(Q - U)) + 
          (2*(-A0 + Q + U)*Log[U])/((mD - Q)^2*(Q - U))) + 
        (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2) + 
           (-2*mD + Q)*delta[A0, Q, mD])*phi[A0, Q, mD])/
         (mD^2*(mD - Q)^2*(Q - U)) + (2*delta[A0, Q, Q]*phi[A0, Q, Q])/
         ((mD - Q)^2*Q*(Q - U)) - 
        (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2) + 
           (-2*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
         (mD^2*(mD - Q)^2*(Q - U)) - (2*delta[A0, U, Q]*phi[A0, U, Q])/
         ((mD - Q)^2*Q*(Q - U))) + 
      Nc*Xb*((2*(A0 - 2*Q)*Log[Q]^2)/((mD - Q)*(Q - U)) + 
        (2*(-A0 + Q + U)*Log[Q]*Log[U])/((mD - Q)*(Q - U)) + 
        Log[A0]*((2*Log[mD])/(-mD + Q) + 2*((mD - Q)^(-1) + (-Q + U)^(-1))*
           Log[Q] + (2*Log[U])/(Q - U)) + 
        Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/((mD - Q)*(Q - U)) - 
          (2*(-A0 + mD + U)*Log[U])/((mD - Q)*(Q - U))) + 
        (2*delta[A0, Q, mD]*phi[A0, Q, mD])/((mD^2 - mD*Q)*(Q - U)) - 
        (2*delta[A0, Q, Q]*phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)) - 
        (2*delta[A0, U, mD]*phi[A0, U, mD])/((mD^2 - mD*Q)*(Q - U)) + 
        (2*delta[A0, U, Q]*phi[A0, U, Q])/((mD - Q)*Q*(Q - U))) + 
      Yt*(Nc*((2*Log[Q]^2)/(Q - U) + (2*Log[Q]*Log[U])/(-Q + U) + 
          Log[A0]*((2*Log[Q])/(-Q + U) + (2*Log[U])/(Q - U)) - 
          (2*(-(A0*(A0 - 5*Q)) + delta[A0, Q, Q])*phi[A0, Q, Q])/
           (Q^2*(Q - U)) - (2*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2 - 
             delta[A0, U, Q])*phi[A0, U, Q])/(Q^2*(Q - U))) + 
        Nc*Xb*Yb*((2*Log[Q]^2)/((mD - Q)*(-Q + U)) + (2*Log[Q]*Log[U])/
           ((mD - Q)*(Q - U)) + Log[mD]*((2*Log[Q])/((mD - Q)*(Q - U)) + 
            (2*Log[U])/((mD - Q)*(-Q + U))) + 
          (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2 - 
             delta[A0, Q, mD])*phi[A0, Q, mD])/(mD^2*(mD - Q)*(Q - U)) + 
          (2*(-(A0*(A0 - 5*Q)) + delta[A0, Q, Q])*phi[A0, Q, Q])/
           ((mD - Q)*Q^2*(Q - U)) - (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 
             3*mD*U + U^2 - delta[A0, U, mD])*phi[A0, U, mD])/
           (mD^2*(mD - Q)*(Q - U)) + (2*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 
             3*Q*U + U^2 - delta[A0, U, Q])*phi[A0, U, Q])/
           ((mD - Q)*Q^2*(Q - U))) + Nc*Xb^3*Yb*
         ((2*(-2*A0 + mD + 3*Q)*Log[Q]^2)/((mD - Q)^3*(Q - U)) + 
          Log[A0]*((4*Log[mD])/(mD - Q)^3 + (4*Log[Q])/(-mD + Q)^3) - 
          (2*(-2*A0 + mD + Q + 2*U)*Log[Q]*Log[U])/((mD - Q)^3*(Q - U)) + 
          Log[mD]*((2*(-2*A0 + mD + 3*Q)*Log[Q])/((mD - Q)^3*(-Q + U)) + 
            (2*(-2*A0 + mD + Q + 2*U)*Log[U])/((mD - Q)^3*(Q - U))) + 
          (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2) + 
             (-3*mD + Q)*delta[A0, Q, mD])*phi[A0, Q, mD])/
           (mD^2*(mD - Q)^3*(Q - U)) + (2*(A0*(A0 - 5*Q)*(mD - Q) - 
             (mD - 3*Q)*delta[A0, Q, Q])*phi[A0, Q, Q])/((mD - Q)^3*Q^2*
            (Q - U)) - (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*
                U + U^2) + (-3*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
           (mD^2*(mD - Q)^3*(Q - U)) + 
          (2*(-((mD - Q)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2)) + 
             (mD - 3*Q)*delta[A0, U, Q])*phi[A0, U, Q])/((mD - Q)^3*Q^2*
            (Q - U))) + Nc*Xb^2*((2*(A0 - mD - Q)*Log[Q]^2)/
           ((mD - Q)^2*(Q - U)) + Log[A0]*((-2*Log[mD])/(mD - Q)^2 + 
            (2*Log[Q])/(mD - Q)^2) + (2*(-A0 + mD + U)*Log[Q]*Log[U])/
           ((mD - Q)^2*(Q - U)) + Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/
             ((mD - Q)^2*(Q - U)) - (2*(-A0 + mD + U)*Log[U])/
             ((mD - Q)^2*(Q - U))) + (2*delta[A0, Q, mD]*phi[A0, Q, mD])/
           (mD*(mD - Q)^2*(Q - U)) + (2*(A0*(A0 - 5*Q)*(-mD + Q) + 
             (mD - 2*Q)*delta[A0, Q, Q])*phi[A0, Q, Q])/((mD - Q)^2*Q^2*
            (Q - U)) - (2*delta[A0, U, mD]*phi[A0, U, mD])/
           (mD*(mD - Q)^2*(Q - U)) + 
          (2*((mD - Q)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2) - 
             (mD - 2*Q)*delta[A0, U, Q])*phi[A0, U, Q])/((mD - Q)^2*Q^2*
            (Q - U))))) + Nc*Xb^4*(16/(mD - Q)^2 - (4*(mD + Q)*Log[Q]^2)/
       (mD - Q)^3 + (8*Log[q2])/(mD - Q)^2 + 
      Log[mD]*((-2*(7*mD + Q))/(mD - Q)^3 + (4*(mD + Q)*Log[Q])/(mD - Q)^3 - 
        (4*(mD + Q)*Log[q2])/(mD - Q)^3) + 
      Log[Q]*((2*(3*mD + 5*Q))/(mD - Q)^3 + (4*(mD + Q)*Log[q2])/
         (mD - Q)^3) - (6*PolyLog[2, 1 - mD/Q])/(mD - Q)^2) + 
    Nc*Xb^2*(8/(-mD + Q) + (2*(A0 + 3*mD - 3*Q)*Log[Q]^2)/(mD - Q)^2 + 
      Log[A0]*((2*(mD - 2*Q + U)*Log[mD])/(mD - Q)^2 - 
        (2*(mD - 2*Q + U)*Log[Q])/(mD - Q)^2) + (4*Log[q2])/(-mD + Q) + 
      Log[mD]*((2*(3*mD + Q))/(mD - Q)^2 - (2*(A0 + 3*mD - 3*Q)*Log[Q])/
         (mD - Q)^2 + (4*mD*Log[q2])/(mD - Q)^2 - (2*(-A0 + Q + U)*Log[U])/
         (mD - Q)^2) + Log[Q]*((-2*(mD + 3*Q))/(mD - Q)^2 - 
        (4*mD*Log[q2])/(mD - Q)^2 + (2*(-A0 + Q + U)*Log[U])/(mD - Q)^2) + 
      (2*delta[A0, Q, mD]*phi[A0, Q, mD])/(mD*(mD - Q)^2) + 
      (2*(A0*(A0 - 5*Q)*(-mD + Q) + (mD - 2*Q)*delta[A0, Q, Q])*
        phi[A0, Q, Q])/((mD - Q)^2*Q^2) + 
      (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2) + 
         (-2*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/(mD^2*(mD - Q)^2) + 
      (2*delta[A0, U, Q]*phi[A0, U, Q])/((mD - Q)^2*Q) + 
      (12*PolyLog[2, 1 - mD/Q])/(mD - Q)) + 
    (cbe^2*(Nc*Xb^4*((A0*(mD + 5*Q))/(Q*(-mD + Q)^3) + 
          Log[A0]*((A0*(mD + 5*Q))/((mD - Q)^3*Q) - (2*A0*(2*mD + Q)*Log[mD])/
             (mD - Q)^4 + (2*A0*(2*mD + Q)*Log[Q])/(mD - Q)^4) + 
          (A0*(mD + 5*Q)*Log[q2])/(Q*(-mD + Q)^3) + 
          Log[Q]*((-2*A0*(2*mD + Q))/(mD - Q)^4 - (2*A0*(2*mD + Q)*Log[q2])/
             (mD - Q)^4) + Log[mD]*((2*A0*(2*mD + Q))/(mD - Q)^4 + 
            (2*A0*(2*mD + Q)*Log[q2])/(mD - Q)^4)) + 
        Nc*(1/2 + Pi^2/3 - A0/Q + Log[A0]^2 + (1 - A0/Q)*Log[q2] + 
          2*Log[q2]^2 + Log[A0]*((A0 + Q)/Q + Log[mD] - 2*Log[q2] - Log[U]) + 
          Log[mD]*(-2 - 2*Log[q2] + Log[U]) + 
          ((A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2 - 
             delta[A0, U, mD])*phi[A0, U, mD])/mD^2) + 
        Nc*Xb^2*((2*(A0 - 2*Q))/((mD - Q)*Q) + Log[A0]*
           ((2*A0)/(Q*(-mD + Q)) + ((A0 - Q + U)*Log[mD])/(mD - Q)^2 - 
            ((A0 - Q + U)*Log[Q])/(mD - Q)^2) + (2*(A0 - Q)*Log[q2])/
           ((mD - Q)*Q) + Log[mD]*((2*(-A0 + mD + Q))/(mD - Q)^2 - 
            (2*(A0 - Q)*Log[q2])/(mD - Q)^2 - ((-A0 + Q + U)*Log[U])/
             (mD - Q)^2) + Log[Q]*((2*(A0 - 2*Q))/(mD - Q)^2 + 
            (2*(A0 - Q)*Log[q2])/(mD - Q)^2 + ((-A0 + Q + U)*Log[U])/
             (mD - Q)^2) + (((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*
                mD*U + U^2) + (-2*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
           (mD^2*(mD - Q)^2) + (delta[A0, U, Q]*phi[A0, U, Q])/
           ((mD - Q)^2*Q)) + 
        Yt*(Nc*Xb*(Log[A0]*((2*Log[mD])/(mD - Q) + (2*Log[Q])/(-mD + Q)) + 
            Log[mD]*(4/(-mD + Q) + (4*Log[q2])/(-mD + Q) + (2*Log[U])/(mD - 
                Q)) + Log[Q]*(4/(mD - Q) + (4*Log[q2])/(mD - Q) + 
              (2*Log[U])/(-mD + Q)) + (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 
               3*mD*U + U^2 - delta[A0, U, mD])*phi[A0, U, mD])/
             (mD^2*(mD - Q)) - (2*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^
                2 - delta[A0, U, Q])*phi[A0, U, Q])/((mD - Q)*Q^2)) + 
          Nc*Xb^3*(-16/(mD - Q)^2 + Log[A0]*((-2*(2*A0 + mD + Q - 2*U)*
                Log[mD])/(mD - Q)^3 + (2*(2*A0 + mD + Q - 2*U)*Log[Q])/
               (mD - Q)^3) - (8*Log[q2])/(mD - Q)^2 + 
            Log[mD]*((4*(3*mD + Q))/(mD - Q)^3 + (4*(mD + Q)*Log[q2])/
               (mD - Q)^3 - (2*(-2*A0 + mD + Q + 2*U)*Log[U])/(mD - Q)^3) + 
            Log[Q]*((-4*(mD + 3*Q))/(mD - Q)^3 - (4*(mD + Q)*Log[q2])/
               (mD - Q)^3 + (2*(-2*A0 + mD + Q + 2*U)*Log[U])/(mD - Q)^3) + 
            (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2) + 
               (-3*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
             (mD^2*(mD - Q)^3) + (2*(-((mD - Q)*(A0^2 - 3*A0*Q + 2*Q^2 - 
                  2*A0*U - 3*Q*U + U^2)) + (mD - 3*Q)*delta[A0, U, Q])*
              phi[A0, U, Q])/(Q^2*(-mD + Q)^3))) + 
        Yt^2*(Nc*Xb^4*((mD + 11*Q)/(Q*(-mD + Q)^3) + 
            Log[A0]*((-((A0 + Q - U)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 
                   3*Q*U + U^2)) + (A0 + 2*Q - U)*delta[A0, U, Q])/(2*
                (mD - Q)^2*Q^2*delta[A0, U, Q]) - ((3*A0 + 2*mD + Q - 3*U)*
                Log[mD])/(mD - Q)^4 + ((3*A0 + 2*mD + Q - 3*U)*Log[Q])/
               (mD - Q)^4) + ((mD + 5*Q)*Log[q2])/(Q*(-mD + Q)^3) + 
            (((A0 - Q - U)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2) + 
               (-A0 + 2*Q + U)*delta[A0, U, Q])*Log[U])/(2*(mD - Q)^2*Q^2*
              delta[A0, U, Q]) + Log[mD]*((2*(5*mD + Q))/(mD - Q)^4 + 
              (2*(2*mD + Q)*Log[q2])/(mD - Q)^4 - ((-3*A0 + 2*mD + Q + 3*U)*
                Log[U])/(mD - Q)^4) + Log[Q]*(-((mD^2 + 2*mD*Q + 9*Q^2 - 
                 ((mD - Q)^2*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2))/
                  delta[A0, U, Q])/((mD - Q)^4*Q)) - (2*(2*mD + Q)*Log[q2])/
               (mD - Q)^4 + ((-3*A0 + 2*mD + Q + 3*U)*Log[U])/(mD - Q)^4) + 
            (((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2) + 
               (-4*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
             (mD^2*(mD - Q)^4) + ((-((mD - Q)*(A0^2*(3*mD - 7*Q) + 
                  2*(mD - 5*Q)*Q^2 - 6*(mD - 3*Q)*Q*U + (3*mD - 7*Q)*U^2 + 
                  2*A0*(-3*mD*(Q + U) + Q*(9*Q + 7*U)))) + ((mD - Q)^2*
                 (A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2)*(A0^2 + 
                  U*(-Q + U) - A0*(Q + 2*U)))/delta[A0, U, Q] + 2*
                (mD^2 - 4*mD*Q + 6*Q^2)*delta[A0, U, Q])*phi[A0, U, Q])/
             (2*(mD - Q)^4*Q^3)) + Nc*(-Q^(-1) + 
            ((-((A0 + Q - U)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + 
                  U^2)) + (A0 + 2*Q - U)*delta[A0, U, Q])*Log[A0])/
             (2*Q^2*delta[A0, U, Q]) + ((A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*
                U + U^2 - delta[A0, U, Q])*Log[Q])/(Q*delta[A0, U, Q]) - 
            Log[q2]/Q + (((A0 - Q - U)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 
                 3*Q*U + U^2) + (-A0 + 2*Q + U)*delta[A0, U, Q])*Log[U])/
             (2*Q^2*delta[A0, U, Q]) + (((A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 
                 3*Q*U + U^2)*(A0^2 + U*(-Q + U) - A0*(Q + 2*U)) + 
               delta[A0, U, Q]*(-3*A0^2 - 2*Q^2 + 6*Q*U - 3*U^2 + 
                 6*A0*(Q + U) + 2*delta[A0, U, Q]))*phi[A0, U, Q])/
             (2*Q^3*delta[A0, U, Q])) + Nc*Xb^2*(2/((mD - Q)*Q) + 
            Log[A0]*(((A0 + Q - U)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + 
                  U^2) + (-A0 - 2*Q + U)*delta[A0, U, Q])/((mD - Q)*Q^2*
                delta[A0, U, Q]) + Log[mD]/(mD - Q)^2 - Log[Q]/(mD - Q)^2) + 
            (2*Log[q2])/(mD*Q - Q^2) + ((-((A0 - Q - U)*(A0^2 - 3*A0*Q + 
                  2*Q^2 - 2*A0*U - 3*Q*U + U^2)) + (A0 - 2*Q - U)*
                delta[A0, U, Q])*Log[U])/((mD - Q)*Q^2*delta[A0, U, Q]) + 
            Log[Q]*((2*(-((mD - Q)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + 
                    U^2)) + mD*delta[A0, U, Q]))/((mD - Q)^2*Q*delta[A0, U, 
                 Q]) + (2*Log[q2])/(mD - Q)^2 - Log[U]/(mD - Q)^2) + 
            Log[mD]*(-2/(mD - Q)^2 - (2*Log[q2])/(mD - Q)^2 + 
              Log[U]/(mD - Q)^2) + ((A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*
                U + U^2 - delta[A0, U, mD])*phi[A0, U, mD])/
             (mD^2*(mD - Q)^2) + (((mD - Q)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 
                 3*Q*U + U^2)*(-A0^2 + (Q - U)*U + A0*(Q + 2*U)) + 
               delta[A0, U, Q]*(A0^2*(3*mD - 4*Q) + 2*(mD - 2*Q)*Q^2 + 
                 3*Q*(-2*mD + 3*Q)*U + (3*mD - 4*Q)*U^2 + A0*(-6*mD*(Q + U) + 
                   Q*(9*Q + 8*U)) + (-2*mD + 3*Q)*delta[A0, U, Q]))*
              phi[A0, U, Q])/((mD - Q)^2*Q^3*delta[A0, U, Q])))) + 
      Nc*Xb^2*((2*(-2*mu2 + 3*Q + U))/((mD - Q)*Q) + (2*(mD - mu2)*Log[mD]^2)/
         (mD - Q)^2 + (4*mu2*Log[mu2])/(mD*Q - Q^2) - (2*(mD - mu2)*Log[Q]^2)/
         (mD - Q)^2 + (2*(-2*mu2 + 2*Q + U)*Log[q2])/((mD - Q)*Q) + 
        (2*U*Log[U])/(Q*(-mD + Q)) + Log[Q]*((2*(mD - 2*mu2 + 2*Q + U))/
           (mD - Q)^2 + (2*(2*mD - 2*mu2 + U)*Log[q2])/(mD - Q)^2 - 
          (2*U*Log[U])/(mD - Q)^2) + Log[mD]*((-2*(3*mD - 2*mu2 + U))/
           (mD - Q)^2 - (2*(2*mD - 2*mu2 + U)*Log[q2])/(mD - Q)^2 + 
          (2*U*Log[U])/(mD - Q)^2) + (4*(mD - mu2)*PolyLog[2, 1 - mu2/mD])/
         (mD - Q)^2 - (4*(mD - mu2)*PolyLog[2, 1 - mu2/Q])/(mD - Q)^2) + 
      Nc*Xb^4*((mD*(-2*mu2 + 9*Q + U) + Q*(-16*mu2 + 3*Q + 5*U))/
         (Q*(-mD + Q)^3) + ((-mD + mu2)*(mD - 3*mu2 + 2*Q)*Log[mD]^2)/
         (mD - Q)^4 + ((mD - mu2)*(mD - 3*mu2 + 2*Q)*Log[Q]^2)/(mD - Q)^4 + 
        Log[mu2]*((2*mu2*(mD + 2*Q))/(Q*(-mD + Q)^3) - (6*mu2^2*Log[Q])/
           (mD - Q)^4) + ((mD*(-2*mu2 + 5*Q + U) + Q*(-10*mu2 + Q + 5*U))*
          Log[q2])/(Q*(-mD + Q)^3) + ((mD + 5*Q)*U*Log[U])/((mD - Q)^3*Q) + 
        Log[mD]*((2*(3*mD^2 + Q*(-2*mu2 + U) + mD*(-7*mu2 + 3*Q + 2*U)))/
           (mD - Q)^4 + (6*mu2^2*Log[mu2])/(mD - Q)^4 + 
          (2*(mD^2 + Q*(-2*mu2 + U) + 2*mD*(-2*mu2 + Q + U))*Log[q2])/
           (mD - Q)^4 - (2*(2*mD + Q)*U*Log[U])/(mD - Q)^4) + 
        Log[Q]*(-((mD^2 + Q*(-10*mu2 + Q + 2*U) + mD*(-8*mu2 + 10*Q + 4*U))/
            (mD - Q)^4) - (2*(mD^2 + Q*(-2*mu2 + U) + 2*mD*(-2*mu2 + Q + U))*
            Log[q2])/(mD - Q)^4 + (2*(2*mD + Q)*U*Log[U])/(mD - Q)^4) - 
        (2*(mD - mu2)*(mD - 3*mu2 + 2*Q)*PolyLog[2, 1 - mu2/mD])/(mD - Q)^4 + 
        (2*(mD - mu2)*(mD - 3*mu2 + 2*Q)*PolyLog[2, 1 - mu2/Q])/(mD - Q)^4) + 
      Nc*(-1/2 - U/Q + mu2*(2/Q + (-mu2 + U)^(-1)) + Log[Q]^2 + 
        Log[Q]*(-1 - 2*Log[q2]) + 2*Log[q2]^2 + 
        (U*(mu2^2 + 2*mu2*(Q - U) + U*(Q + U))*Log[U])/(Q*(mu2 - U)^2) + 
        (1 - (2*mu2^2)/(mu2 - U)^2)*Log[U]^2 + 
        Log[q2]*((2*mu2*(mu2 + Q) - 3*mu2*U + U^2)/(Q*(mu2 - U)) + 
          (-2 + (2*mu2^2)/(mu2 - U)^2)*Log[U]) + 
        Log[mu2]*(-((mu2*(2*mu2^2 + mu2*(Q - 4*U) + 2*U*(Q + U)))/
            (Q*(mu2 - U)^2)) - (2*mu2^2*Log[q2])/(mu2 - U)^2 + 
          (2*mu2^2*Log[U])/(mu2 - U)^2) + 2*PolyLog[2, 1 - mu2/Q] + 
        (2 - (4*mu2^2)/(mu2 - U)^2)*PolyLog[2, 1 - mu2/U]))/sbe^2 + 
    (Nc*Xb*((4*Sqrt[mu2]*sgn*Log[mD]^2)/(-mD + Q) + 
        (4*Sqrt[mu2]*sgn*Log[Q]^2)/(mD - Q) + 
        Log[mD]*((8*Sqrt[mu2]*sgn)/(mD - Q) + (8*Sqrt[mu2]*sgn*Log[q2])/
           (mD - Q)) + Log[Q]*((8*Sqrt[mu2]*sgn)/(-mD + Q) + 
          (8*Sqrt[mu2]*sgn*Log[q2])/(-mD + Q)) + 
        (8*Sqrt[mu2]*sgn*PolyLog[2, 1 - mu2/mD])/(-mD + Q) + 
        (8*Sqrt[mu2]*sgn*PolyLog[2, 1 - mu2/Q])/(mD - Q)) + 
      Nc*Xb^3*((32*Sqrt[mu2]*sgn)/(mD - Q)^2 + 
        (4*Sqrt[mu2]*(mD - 2*mu2 + Q)*sgn*Log[mD]^2)/(mD - Q)^3 + 
        (16*mu2^(3/2)*sgn*Log[mu2]*Log[Q])/(-mD + Q)^3 - 
        (4*Sqrt[mu2]*(mD - 2*mu2 + Q)*sgn*Log[Q]^2)/(mD - Q)^3 + 
        (16*Sqrt[mu2]*sgn*Log[q2])/(mD - Q)^2 + 
        Log[mD]*((8*Sqrt[mu2]*(3*mD + Q)*sgn)/(-mD + Q)^3 + 
          (16*mu2^(3/2)*sgn*Log[mu2])/(mD - Q)^3 - 
          (8*Sqrt[mu2]*(mD + Q)*sgn*Log[q2])/(mD - Q)^3) + 
        Log[Q]*((8*Sqrt[mu2]*(mD + 3*Q)*sgn)/(mD - Q)^3 + 
          (8*Sqrt[mu2]*(mD + Q)*sgn*Log[q2])/(mD - Q)^3) + 
        (8*Sqrt[mu2]*(mD - 2*mu2 + Q)*sgn*PolyLog[2, 1 - mu2/mD])/
         (mD - Q)^3 - (8*Sqrt[mu2]*(mD - 2*mu2 + Q)*sgn*
          PolyLog[2, 1 - mu2/Q])/(mD - Q)^3) + 
      Nc*Xt*((4*Sqrt[mu2]*(mu2 + Q)*sgn*Log[Q]^2)/((mu2 - Q)*(-Q + U)) + 
        (8*Sqrt[mu2]*Q*sgn*Log[Q]*Log[q2])/((mu2 - Q)*(Q - U)) + 
        (8*Sqrt[mu2]*sgn*U*Log[q2]*Log[U])/((mu2 - U)*(-Q + U)) - 
        (4*Sqrt[mu2]*sgn*(mu2 + U)*Log[U]^2)/((mu2 - U)*(-Q + U)) + 
        Log[mu2]*((8*mu2^(3/2)*sgn*Log[Q])/((mu2 - Q)*(Q - U)) + 
          (8*mu2^(3/2)*sgn*Log[q2])/((mu2 - Q)*(-mu2 + U)) + 
          (8*mu2^(3/2)*sgn*Log[U])/((mu2 - U)*(-Q + U))) + 
        (8*Sqrt[mu2]*(mu2 + Q)*sgn*PolyLog[2, 1 - mu2/Q])/
         ((mu2 - Q)*(-Q + U)) - (8*Sqrt[mu2]*sgn*(mu2 + U)*
          PolyLog[2, 1 - mu2/U])/((mu2 - U)*(-Q + U))))/(cbe*sbe) + 
    (Nc*(-2*Log[mD]^2 + 4*mu2*((mD - mu2)^(-1) + (-mu2 + Q)^(-1))*Log[mu2] - 
        2*Log[Q]^2 - 2*Log[q2]^2 + Log[Q]*((4*mu2)/(mu2 - Q) + 2*Log[q2]) + 
        Log[mD]*((4*mu2)/(-mD + mu2) + 2*Log[Q] + 2*Log[q2]) - 
        4*PolyLog[2, 1 - mu2/mD] - 4*PolyLog[2, 1 - mu2/Q]) + 
      Nc*Xb^2*(4/(mD - Q) + (2*Q*Log[mD]^2)/(mD - Q)^2 + 
        (4*mu2*Log[mu2])/((mD - mu2)*(mD - Q)) + (2*Log[q2])/(mD - Q) + 
        Log[mD]*((-4*(mD^2 - mu2*Q))/((mD - mu2)*(mD - Q)^2) - 
          (2*Q*Log[Q])/(mD - Q)^2 - (2*Q*Log[q2])/(mD - Q)^2) + 
        Log[Q]*((2*(mD + Q))/(mD - Q)^2 + (2*Q*Log[q2])/(mD - Q)^2) + 
        (4*Q*PolyLog[2, 1 - mu2/mD])/(mD - Q)^2 - (4*Q*PolyLog[2, 1 - mu2/Q])/
         (mD - Q)^2) + Xt^2*(Nc*Xb^2*(2/((mD - Q)*(Q - U)) - 
          (2*Q*U*Log[Q]^2)/((mD - Q)^2*(Q - U)^2) + (2*U*Log[U])/
           ((mD - Q)*(Q - U)^2) + Log[mD]*((2*Q)/((mD - Q)^2*(-Q + U)) + 
            (2*Q*U*Log[Q])/((mD - Q)^2*(Q - U)^2) - (2*Q*U*Log[U])/
             ((mD - Q)^2*(Q - U)^2)) + Log[Q]*((2*(Q^2 - mD*U))/
             ((mD - Q)^2*(Q - U)^2) + (2*Q*U*Log[U])/((mD - Q)^2*
              (Q - U)^2))) + Nc*(4/(Q - U) + (4*mu2*Log[mu2])/
           ((mu2 - Q)*(-Q + U)) + (2*U*Log[Q]^2)/(Q - U)^2 + 
          Log[Q]*((4*(Q^2 - mu2*U))/((mu2 - Q)*(Q - U)^2) - 
            (2*U*Log[q2])/(Q - U)^2) + (4*U*Log[U])/(Q - U)^2 - 
          (2*U*Log[U]^2)/(Q - U)^2 + Log[q2]*(2/(Q - U) + 
            (2*U*Log[U])/(Q - U)^2) + Log[mD]*(2/(Q - U) - 
            (2*U*Log[Q])/(Q - U)^2 + (2*U*Log[U])/(Q - U)^2) + 
          (4*U*PolyLog[2, 1 - mu2/Q])/(Q - U)^2 - (4*U*PolyLog[2, 1 - mu2/U])/
           (Q - U)^2)))/cbe^2 + 
    Xt^2*(Nc^2*((2*Q*U*Log[Q]^2)/(Q - U)^3 + Log[mD]*(-((Q + U)/(Q - U)^2) + 
          (2*Q*U*Log[Q])/(Q - U)^3 - (2*Q*U*Log[U])/(Q - U)^3) + 
        Log[Q]*(-((Q + U)/(Q - U)^2) - (4*Q*U*Log[q2])/(Q - U)^3 - 
          (2*Q*U*Log[U])/(Q - U)^3) + Log[q2]*((2*(Q + U))/(Q - U)^2 + 
          (4*Q*U*Log[U])/(Q - U)^3)) + Nc*((-5*Q + U)/(Q^2 - Q*U) + 
        ((-3*Q + U)*Log[Q]^2)/(Q - U)^2 + (U*(-5*Q + U)*Log[U])/
         (Q*(Q - U)^2) + (3*Log[U]^2)/(-Q + U) + 
        Log[Q]*((3*Q + U)/(Q - U)^2 + (2*U*Log[q2])/(Q - U)^2 + 
          ((6*Q - 4*U)*Log[U])/(Q - U)^2) + Log[q2]*((-3*Q + U)/(Q^2 - Q*U) - 
          (2*U*Log[U])/(Q - U)^2) + (6*PolyLog[2, 1 - Q/U])/(-Q + U)) + 
      Xb^4*(Nc^2*((-2*(Q + U))/((mD - Q)^2*(Q - U)^2) + 
          (2*Q*(mD + Q)*U*Log[Q]^2)/((mD - Q)^3*(Q - U)^3) - 
          (4*Q*U*Log[U])/((mD - Q)^2*(Q - U)^3) + 
          Log[mD]*(((mD + Q)*(Q + U))/((mD - Q)^3*(Q - U)^2) + 
            (2*Q*(mD + Q)*U*Log[Q])/((-mD + Q)^3*(Q - U)^3) + 
            (2*Q*(mD + Q)*U*Log[U])/((mD - Q)^3*(Q - U)^3)) + 
          Log[Q]*((-(Q^2*(mD + Q)) + 4*(mD - Q)*Q*U + (mD + Q)*U^2)/
             ((mD - Q)^3*(Q - U)^3) + (2*Q*(mD + Q)*U*Log[U])/
             ((-mD + Q)^3*(Q - U)^3))) + Nc*((mD + 11*Q)/(Q*(-mD + Q)^3) + 
          ((-2*Q*(3*mD + Q) + 4*(mD + Q)*U)*Log[Q]^2)/((mD - Q)^3*
            (Q - U)^2) + ((mD + 5*Q)*Log[q2])/(Q*(-mD + Q)^3) + 
          (U*Log[U])/((mD - Q)^2*Q*(-Q + U)) + 
          Log[Q]*((-4*(mD + 2*Q))/(mD - Q)^4 + 3/((mD - Q)^2*(Q - U)) - 
            (2*(2*mD + Q)*Log[q2])/(mD - Q)^4 - 
            (2*(mD - U)*(-3*mD*Q + 2*mD*U + Q*U)*Log[U])/((mD - Q)^4*
              (Q - U)^2)) + Log[mD]*((-2*(mD^2 - 7*mD*Q + 5*mD*U + Q*U))/
             ((mD - Q)^4*(Q - U)) + ((2*Q*(3*mD + Q) - 4*(mD + Q)*U)*Log[Q])/
             ((mD - Q)^3*(Q - U)^2) + (2*(2*mD + Q)*Log[q2])/(mD - Q)^4 + 
            (2*(mD - U)*(-3*mD*Q + 2*mD*U + Q*U)*Log[U])/((mD - Q)^4*
              (Q - U)^2)) + (6*PolyLog[2, 1 - mD/Q])/((mD - Q)^2*(-Q + U)) + 
          (6*(mD - U)^2*PolyLog[2, 1 - mD/U])/((mD - Q)^4*(Q - U)) - 
          (6*(mD - U)^2*PolyLog[2, 1 - Q/U])/((mD - Q)^4*(Q - U)))) + 
      Xb^2*(Nc^2*((4*Q*U*Log[Q]^2)/((-mD + Q)*(Q - U)^3) + 
          Log[Q]*((2*(Q + U))/((mD - Q)*(Q - U)^2) + (4*Q*U*Log[U])/
             ((mD - Q)*(Q - U)^3)) + Log[mD]*
           ((2*(Q + U))/((-mD + Q)*(Q - U)^2) + (4*Q*U*Log[Q])/
             ((mD - Q)*(Q - U)^3) + (4*Q*U*Log[U])/((-mD + Q)*(Q - U)^3))) + 
        Nc*(2/(mD*Q - Q^2) + ((3*mD^2 + 6*mD*Q - 5*Q^2 - 8*mD*U + 4*Q*U)*
            Log[Q]^2)/((mD - Q)^2*(Q - U)^2) + (2*Log[q2])/(mD*Q - Q^2) + 
          (2*U*Log[U])/((mD - Q)*Q*(Q - U)) - (3*Log[U]^2)/(Q - U)^2 + 
          Log[mD]*((2*(mD - U))/((mD - Q)^2*(-Q + U)) + 
            ((-6*mD^2 + 2*Q*(Q - 2*U) + 8*mD*U)*Log[Q])/((mD - Q)^2*
              (Q - U)^2) - (2*Log[q2])/(mD - Q)^2 + 
            (2*(3*mD^2 - 4*mD*U + U^2)*Log[U])/((mD - Q)^2*(Q - U)^2)) + 
          Log[Q]*(2/(mD - Q)^2 + (2*Log[q2])/(mD - Q)^2 - 
            (2*(6*mD*Q - 3*Q^2 - 4*mD*U + U^2)*Log[U])/((mD - Q)^2*
              (Q - U)^2)) + (6*(mD + Q - 2*U)*PolyLog[2, 1 - mD/Q])/
           ((mD - Q)*(Q - U)^2) - (6*(mD - U)^2*PolyLog[2, 1 - mD/U])/
           ((mD - Q)^2*(Q - U)^2) + ((12*mD - 6*(Q + U))*PolyLog[2, 1 - Q/U])/
           ((mD - Q)^2*(Q - U))))))) /. { delta[x_,y_,z_] :> TDelta[x,y,z], phi[x_,y_,z_] :> TPhi[x,y,z] }
    ];


(* Tau Yukawa lambda 2-loop threshold correction, computation with general masses *)
lambda2LHSSTau = With[{
    k = 1/(4*Pi)^2,
    sbe = Sqrt[TanBeta^2/(1+TanBeta^2)],
    cbe = Sqrt[1/(1+TanBeta^2)],
    Nc = 3, (* number of colors *)
    Q  = msq2[3,3]*(1+0.03),
    mD = msd2[3,3]*(1-0.03),
    L  = msl2[3,3]*(1-0.01),
    R  = mse2[3,3]*(1+0.01),
    mu2 = MuInput^2*(1-0.05),
    q2 = SCALE^2, (* renormalization/matching scale *)
    gt = Yu[3,3], (* SM Yukawa coupling *)
    gb = Yd[3,3],
    gtau = Ye[3,3],
    A0 = mAInput^2*(1-0.015),
    Xt = xt,
    Yt = yt,
    Xb = xb,
    Yb = yb,
    Xtau = xtau,
    Ytau = ytau,
    mQ3 = Sqrt[msq2[3,3]],
    mU3 = Sqrt[msu2[3,3]],
    mD3 = Sqrt[msd2[3,3]],
    M3 = M3Input,
    Mu = MuInput,
    Q2 = SCALE^2,
    Xtildet = xtt,
    Xtildeb = xbb,
    cosb = Cos[ArcTan[TanBeta]],
    sinb = Sin[ArcTan[TanBeta]],
    MA = mAInput,
    CF = 4/3
    },
    ybMSSM[mQ3_,mU3_,mD3_,M3_,Mu_,TanBeta_,Xt_,Xb_] := Module[{deltagsb, deltagbyL1, deltagbyL2, deltagbyL3, deltagbyL4},
          deltagsb = - g3^2*CF*k*(1+Log[M3^2/Q2]+TCF[6][mQ3/M3]+TCF[6][mD3/M3]-Xb/M3*TCF[9][mQ3/M3,mD3/M3]);
          deltagbyL1 = - gb^2/cosb^2*k*(3/4*Log[Mu^2/Q2]+3/8*sinb^2*(2*Log[MA^2/Q2]-1)+TCF[6][mQ3/Mu]+1/2*TCF[6][mD3/Mu]);
          deltagbyL2 = - gt^2/sinb^2*k*(1/4*Log[Mu^2/Q2]+1/8*cosb^2*(2*Log[MA^2/Q2]-1)+sinb^2*(Log[MA^2/Q2]-1));
          deltagbyL3 = - gt^2/sinb^2*k*(1/2*TCF[6][mU3/Mu]+(Xt*TanBeta)/(2*Mu)*TCF[9][mQ3/Mu,mU3/Mu] );
          deltagbyL4 = - 1/2*(-gt^2*Nc*k*Xtildet/6*TCF[5][xQU]-gb^2*Nc*k*Xtildeb/6*TCF[5][xQD]);
          (gb/(1-deltagsb-(deltagbyL1+deltagbyL2+deltagbyL3+deltagbyL4)))
    ];
    k^2*Nc*((ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^4*gtau^2*Xtau*(Xb^3*(-8/(mD - Q)^2 +
        Log[L]*((8*L)/((mD - Q)^2*(L - R)) + (4*L*(mD + Q)*Log[Q])/
           ((mD - Q)^3*(L - R))) + (-8/(mD - Q)^2 - (4*(mD + Q)*Log[Q])/
           (mD - Q)^3)*Log[q2] + (8*R*Log[R])/((mD - Q)^2*(-L + R)) +
        Log[Q]*((-4*(mD + Q))/(mD - Q)^3 - (4*(mD + Q)*R*Log[R])/
           ((mD - Q)^3*(L - R))) + Log[mD]*((4*(mD + Q))/(mD - Q)^3 +
          (4*L*(mD + Q)*Log[L])/((mD - Q)^3*(-L + R)) + (4*(mD + Q)*Log[q2])/
           (mD - Q)^3 + (4*(mD + Q)*R*Log[R])/((mD - Q)^3*(L - R)))) +
      Xb*((4*L*Log[L]*Log[Q])/((mD - Q)*(-L + R)) + (4*Log[Q]*Log[q2])/
         (mD - Q) + Log[Q]*(4/(mD - Q) + (4*R*Log[R])/((mD - Q)*(L - R))) +
        Log[mD]*(4/(-mD + Q) + (4*L*Log[L])/((mD - Q)*(L - R)) +
          (4*Log[q2])/(-mD + Q) + (4*R*Log[R])/((mD - Q)*(-L + R))))))/(cbe*sbe) +
   (ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2*gtau^4*(Xb*Xtau*(Log[L]*(4/(-L + R) + (4*Q*Log[Q])/
           ((-mD + Q)*(L - R)) + (4*Log[q2])/(-L + R)) + (4*Log[R])/(L - R) +
        (4*Q*Log[Q]*Log[R])/((mD - Q)*(L - R)) + (4*Log[q2]*Log[R])/(L - R) +
        Log[mD]*((4*mD*Log[L])/((mD - Q)*(L - R)) + (4*mD*Log[R])/
           ((mD - Q)*(-L + R)))) + Xb*Xtau^3*(-8/(L - R)^2 +
        Log[L]*((4*(L + R))/(L - R)^3 + (4*Q*(L + R)*Log[Q])/
           ((mD - Q)*(L - R)^3) + (4*(L + R)*Log[q2])/(L - R)^3) -
        (4*(L + R)*Log[R])/(L - R)^3 + Log[q2]*(-8/(L - R)^2 -
          (4*(L + R)*Log[R])/(L - R)^3) + Log[mD]*((8*mD)/((mD - Q)*(L - R)^2) -
          (4*mD*(L + R)*Log[L])/((mD - Q)*(L - R)^3) + (4*mD*(L + R)*Log[R])/
           ((mD - Q)*(L - R)^3)) + Log[Q]*((8*Q)/((-mD + Q)*(L - R)^2) -
          (4*Q*(L + R)*Log[R])/((mD - Q)*(L - R)^3)))))/(cbe*sbe)) +
 gtau^6*k^2*(5 + 3*Log[L]^2 + Log[L]*(-6 - 6*Log[q2]) + 5*Log[q2]^2 +
   Log[q2]*(10 - 4*Log[R]) - 4*Log[R] + 2*Log[R]^2 +
   (sbe^2*(3/2 - A0/L + Pi^2 - (2*A0)/R + 3*Log[A0]^2 + 3*Log[L]^2 +
      Log[L]*(-9/2 - 3*Log[q2]) + 2*Log[q2]^2 +
      Log[A0]*(7 + A0*(L^(-1) + 2/R) - 3*Log[L] - 3*Log[R]) +
      Log[q2]*(-((A0*(2*L + R))/(L*R)) - Log[R]) - (5*Log[R])/2 +
      2*Log[R]^2 + ((A0^2 - 7*A0*L - delta[A0, L, L])*phi[A0, L, L])/L^2 +
      Ytau^2*(-L^(-1) - 2/R +
        ((((L - R)^4 - A0^3*(L + R) - A0*(L - R)^2*(3*L + 5*R) +
             A0^2*(3*L^2 + 4*L*R + 5*R^2))*delta[A0, R, L] +
           delta[A0, L, R]*(-(L*(-A0^2 + (L - R)^2)*R) +
             (-L^2 + 3*L*R - R^2 + A0*(L + R))*delta[A0, R, L]))*Log[A0])/
         (2*L*R^2*delta[A0, L, R]*delta[A0, R, L]) +
        (((A0^3 - (L - R)^3 - 3*A0^2*(L + R) + 3*A0*(L^2 - R^2))*
            delta[A0, R, L] + delta[A0, L, R]*(-(R*(A0^2 - L^2 - 2*A0*R +
                R^2)) + (-A0 + L + 2*R)*delta[A0, R, L]))*Log[L])/
         (2*R^2*delta[A0, L, R]*delta[A0, R, L]) + (-L^(-1) - 2/R)*Log[q2] +
        (((A0^3 + (L - R)^3 - A0^2*(L + 5*R) - A0*(L^2 + 4*L*R - 5*R^2))*
            delta[A0, R, L] - delta[A0, L, R]*(2*L*(A0 + L - R)*R +
             (A0 + L - 3*R)*delta[A0, R, L]))*Log[R])/(2*L*R*delta[A0, L, R]*
          delta[A0, R, L]) + ((A0^4 + (L - R)^4 - 4*A0^3*(L + R) -
           4*A0*(L - R)^2*(L + R) + 2*A0^2*(3*L^2 + 2*L*R + R^2) -
           3*(A0^2 + (L - R)^2 - 2*A0*(L + R))*delta[A0, L, R] +
           2*delta[A0, L, R]^2)*phi[A0, L, R])/(2*R^3*delta[A0, L, R]) +
        ((-(A0 + L - R)^2 + delta[A0, R, L])*phi[A0, R, L])/
         (2*L*delta[A0, R, L])) + ((A0^2 - 6*A0*R - delta[A0, R, R])*
        phi[A0, R, R])/R^2 + Xtau*Ytau*((6*Log[L]^2)/(L - R) +
        (12*Log[R])/(L - R) + (12*Log[q2]*Log[R])/(L - R) -
        (4*Log[R]^2)/(L - R) + Log[A0]*((2*Log[L])/(L - R) -
          (2*Log[R])/(L - R)) + Log[L]*(-12/(L - R) - (12*Log[q2])/(L - R) -
          (2*Log[R])/(L - R)) - (2*(-(A0*(A0 - 7*L)) + delta[A0, L, L])*
          phi[A0, L, L])/(L^2*(L - R)) + (2*(A0 + L - R)*phi[A0, R, L])/
         (L*(L - R)) + (2*(-(A0*(A0 - 6*R)) + delta[A0, R, R])*phi[A0, R, R])/
         ((L - R)*R^2)) + Xtau^3*Ytau*(-48/(L - R)^2 +
        (2*(4*A0 - 5*L - 3*R)*Log[L]^2)/(L - R)^3 - (12*(L + 3*R)*Log[R])/
         (L - R)^3 + (4*(-A0 + L + R)*Log[R]^2)/(L - R)^3 +
        Log[A0]*((2*(-6*A0 + L - R)*Log[L])/(L - R)^3 +
          (2*(6*A0 - L + R)*Log[R])/(L - R)^3) +
        Log[q2]*(-24/(L - R)^2 - (12*(L + R)*Log[R])/(L - R)^3) +
        Log[L]*((12*(3*L + R))/(L - R)^3 + (12*(L + R)*Log[q2])/(L - R)^3 +
          (2*(-2*A0 + 3*L + R)*Log[R])/(L - R)^3) +
        (2*(A0*(A0 - 7*L)*(L - R) + (-5*L + R)*delta[A0, L, L])*
          phi[A0, L, L])/(L^2*(L - R)^3) -
        (2*((L - R)*(A0 + L - R) - 2*delta[A0, R, L])*phi[A0, R, L])/
         (L*(L - R)^3) + (2*(A0*(A0 - 6*R)*(L - R) - (L - 3*R)*
            delta[A0, R, R])*phi[A0, R, R])/((L - R)^3*R^2)) +
      Xtau^2*((4*A0*L - 2*A0*R - 4*L*R)/(L^2*R - L*R^2) +
        ((A0 + L - 3*R)*Log[L]^2)/(L - R)^2 + ((2*A0 + L - 5*R)*Log[R])/
         (L - R)^2 - (2*Log[R]^2)/(L - R) +
        Log[q2]*((4*A0*L - 2*A0*R - 2*L*R)/(L^2*R - L*R^2) +
          (2*(A0 - L)*Log[R])/(L - R)^2) +
        Log[L]*((-2*A0 + L + 3*R)/(L - R)^2 - (2*(A0 - L)*Log[q2])/
           (L - R)^2 + ((-A0 + L + R)*Log[R])/(L - R)^2) +
        Log[A0]*((2*A0*(-2*L + R))/(L*(L - R)*R) + ((A0 - 5*L + 5*R)*Log[L])/
           (L - R)^2 - ((A0 - 5*L + 5*R)*Log[R])/(L - R)^2) +
        ((A0*(A0 - 7*L)*(L - R) + (-2*L + R)*delta[A0, L, L])*phi[A0, L, L])/
         (L^2*(L - R)^2) + (delta[A0, R, L]*phi[A0, R, L])/(L*(L - R)^2) +
        ((-(A0*(A0 - 6*R)) + delta[A0, R, R])*phi[A0, R, R])/((L - R)*R^2) +
        Ytau^2*((4*L - 2*R)/(L^2*R - L*R^2) + (3*Log[L]^2)/(L - R)^2 +
          ((-((L - R)*((L - R)^2*(2*L - R) + A0^2*(2*L + R) - 4*A0*L*
                 (L + 2*R))*delta[A0, R, L]) + delta[A0, L, R]*
              (2*L*(L - R)*(A0 + L - R)*R + (2*L^2 + L*R - R^2)*delta[A0, R,
                 L]))*Log[R])/(L*(L - R)^2*R*delta[A0, L, R]*
            delta[A0, R, L]) + (2*Log[R]^2)/(L - R)^2 +
          Log[L]*(((L - R)*(-A0^3 + L*(L - R)^2 + A0^2*(3*L + 4*R) -
                A0*(3*L^2 + 2*L*R + 7*R^2))*delta[A0, R, L] +
              delta[A0, L, R]*(-((L - R)*R*(-A0^2 + L^2 + 2*A0*R - R^2)) +
                (-L^2 + A0*(L - R) - 2*L*R + R^2)*delta[A0, R, L]))/
             ((L - R)^2*R^2*delta[A0, L, R]*delta[A0, R, L]) -
            (2*Log[q2])/(L - R)^2 - (5*Log[R])/(L - R)^2) +
          Log[A0]*(((A0^3*L - (L - R)^4 + A0*L*(3*L^2 - 2*L*R - R^2) +
                A0^2*(-3*L^2 - 2*L*R + R^2))*delta[A0, R, L] +
              delta[A0, L, R]*(L*(-A0^2 + (L - R)^2)*R + (-(A0*L) + L^2 -
                  3*L*R + R^2)*delta[A0, R, L]))/(L*(L - R)*R^2*
              delta[A0, L, R]*delta[A0, R, L]) + Log[L]/(L - R)^2 -
            Log[R]/(L - R)^2) + Log[q2]*((4*L - 2*R)/(L^2*R - L*R^2) +
            (2*Log[R])/(L - R)^2) + ((A0*(A0 - 7*L) - delta[A0, L, L])*
            phi[A0, L, L])/(L^2*(L - R)^2) +
          ((-((L - R)*(A0^4 - 4*A0*L^2*(L - R) + (L - R)^4 - 4*A0^3*(L + R) +
                A0^2*(6*L^2 + 4*L*R + 6*R^2))) + (A0^2*(3*L - 5*R) +
               (3*L - 5*R)*(L - R)^2 + A0*(-6*L^2 + 4*L*R + 14*R^2))*
              delta[A0, L, R] - 2*(L - 2*R)*delta[A0, L, R]^2)*phi[A0, L, R])/
           ((L - R)^2*R^3*delta[A0, L, R]) +
          (((L - R)*(A0 + L - R)^2 + A0*delta[A0, R, L])*phi[A0, R, L])/
           (L*(L - R)^2*delta[A0, R, L]) + ((A0*(A0 - 6*R) - delta[A0, R, R])*
            phi[A0, R, R])/((L - R)^2*R^2))) +
      Xtau^4*((3*L*(L - R)*R + A0*(-2*L^2 - 5*L*R + R^2))/(L*(L - R)^3*R) +
        Log[L]*((3*(4*A0*L - L^2 + R^2))/(2*(L - R)^4) +
          (3*(2*A0*L - L^2 + R^2)*Log[q2])/(L - R)^4) +
        (3*(-4*A0*L + L^2 - R^2)*Log[R])/(2*(L - R)^4) +
        Log[q2]*((6*L*(L - R)*R + A0*(-2*L^2 - 5*L*R + R^2))/
           (L*(L - R)^3*R) + (3*(-2*A0*L + L^2 - R^2)*Log[R])/(L - R)^4) +
        Log[A0]*((6*L*R*(-L + R) + A0*(2*L^2 + 5*L*R - R^2))/
           (L*(L - R)^3*R) + (3*(-2*A0*L + L^2 - R^2)*Log[L])/(L - R)^4 +
          (3*(2*A0*L - L^2 + R^2)*Log[R])/(L - R)^4) +
        Ytau^2*((-2*L^2 - 11*L*R + R^2)/(L*(L - R)^3*R) +
          ((7*A0 - 11*L - 3*R)*Log[L]^2)/(L - R)^4 -
          (((L - R)^2*(A0^3 - 3*L^3 + 11*L^2*R - 5*L*R^2 - 3*R^3 - A0^2*
                (5*L + 3*R) + A0*(7*L^2 + 4*L*R + 5*R^2))*delta[A0, R, L] +
             delta[A0, L, R]*(2*L*(L - R)^2*(A0 + L - R)*R + (3*L^3 -
                 A0*(L - R)^2 + 7*L^2*R + 13*L*R^2 + R^3)*delta[A0, R, L]))*
            Log[R])/(2*L*(L - R)^4*R*delta[A0, L, R]*delta[A0, R, L]) -
          (2*(-2*A0 + L + 3*R)*Log[R]^2)/(L - R)^4 +
          Log[q2]*((-2*L^2 - 5*L*R + R^2)/(L*(L - R)^3*R) -
            (6*L*Log[R])/(L - R)^4) + Log[A0]*
           (((L - R)*(-A0^3 + L^3 - 3*L^2*R - L*R^2 + 3*R^3 + 3*A0^2*
                 (L + R) - A0*(3*L^2 + 5*R^2))*delta[A0, R, L] +
              delta[A0, L, R]*(-(L*(-A0^2 + (L - R)^2)*R) +
                (-L^2 + A0*(L - R) + 3*L*R + 3*R^2)*delta[A0, R, L]))/
             (2*L*(L - R)^2*R^2*delta[A0, L, R]*delta[A0, R, L]) -
            (3*(A0 - L + R)*Log[L])/(L - R)^4 + (3*(A0 - L + R)*Log[R])/
             (L - R)^4) + Log[L]*(((L - R)^2*(A0^3 - L^3 + L^2*R + 9*L*R^2 -
                9*R^3 - A0^2*(3*L + 5*R) + A0*(3*L^2 + 4*L*R + 9*R^2))*delta[
                A0, R, L] + delta[A0, L, R]*((L - R)^2*R*(-A0^2 + L^2 +
                  2*A0*R - R^2) + (L^3 - A0*(L - R)^2 + 2*L^2*R + 17*L*R^2 +
                  4*R^3)*delta[A0, R, L]))/(2*(L - R)^4*R^2*delta[A0, L, R]*
              delta[A0, R, L]) + (6*L*Log[q2])/(L - R)^4 +
            ((-11*A0 + 13*L + 9*R)*Log[R])/(L - R)^4) +
          ((A0*(A0 - 7*L)*(L - R) + (-8*L + R)*delta[A0, L, L])*
            phi[A0, L, L])/(L^2*(L - R)^4) +
          (((L - R)^2*(A0^4 - 4*A0^3*(L + R) - 4*A0*(L - R)^2*(L + R) +
               (L - R)^2*(L^2 - 2*L*R - 3*R^2) + A0^2*(6*L^2 + 4*L*R +
                 6*R^2)) - (L - R)*(3*L^3 + 3*A0^2*(L - 3*R) - 15*L^2*R + 29*
                L*R^2 - 17*R^3 - 6*A0*(L^2 - 2*L*R - 3*R^2))*
              delta[A0, L, R] + 2*(L^2 - 5*L*R + 12*R^2)*delta[A0, L, R]^2)*
            phi[A0, L, R])/(2*(L - R)^4*R^3*delta[A0, L, R]) -
          (((L - R)^2*(A0 + L - R)^2 + (4*A0 + 3*L - 3*R)*(L - R)*
              delta[A0, R, L] - 6*delta[A0, R, L]^2)*phi[A0, R, L])/
           (2*L*(L - R)^4*delta[A0, R, L]) +
          ((A0*(A0 - 6*R)*(-L + R) + (L - 5*R)*delta[A0, R, R])*
            phi[A0, R, R])/((L - R)^4*R^2)))))/cbe^2 +
   Xtau^4*((2*(2*L^2 - 27*L*R + R^2))/(L*(L - R)^2*R) +
     ((-7*L - 9*R)*Log[L]^2)/(L - R)^3 + ((-6*L^2 - 44*L*R + 2*R^2)*Log[R])/
      (L*(L - R)^3) + ((3*L + 5*R)*Log[R]^2)/(L - R)^3 +
     Log[L]*((-4*(L^2 - 10*L*R - 3*R^2))/((L - R)^3*R) +
       (2*(5*L + 7*R)*Log[q2])/(L - R)^3 + (4*(L + R)*Log[R])/(L - R)^3) +
     Log[q2]*((2*(2*L^2 - 15*L*R + R^2))/(L*(L - R)^2*R) -
       (2*(5*L + 7*R)*Log[R])/(L - R)^3) - (6*PolyLog[2, 1 - L/R])/
      (L - R)^2) + Xtau^2*((-2*L^2 - 3*L*R + R^2)/(L*(L - R)*R) +
     ((7*L - 9*R)*Log[L]^2)/(L - R)^2 + ((16*L^2 - 21*L*R + R^2)*Log[R])/
      (L*(L - R)^2) - (9*Log[R]^2)/(L - R) +
     Log[L]*((2*L^2 - 17*L*R + 19*R^2)/((L - R)^2*R) -
       (2*(8*L - 9*R)*Log[q2])/(L - R)^2 + (2*L*Log[R])/(L - R)^2) +
     Log[q2]*((-2*L^2 - L*R + R^2)/(L*(L - R)*R) + (2*(8*L - 9*R)*Log[R])/
        (L - R)^2) - (6*PolyLog[2, 1 - L/R])/(L - R)) +
   Xtau^6*((-2*L^2 - 11*L*R + R^2)/(L*(L - R)^3*R) +
     ((-9*L - 5*R)*Log[L]^2)/(L - R)^4 + ((-10*L^2 - 3*L*R + R^2)*Log[R])/
      (L*(L - R)^4) + ((-3*L - 5*R)*Log[R]^2)/(L - R)^4 +
     Log[q2]*((-2*L^2 - 5*L*R + R^2)/(L*(L - R)^3*R) -
       (6*L*Log[R])/(L - R)^4) + Log[L]*((2*L^2 + 13*L*R - 3*R^2)/
        ((L - R)^4*R) + (6*L*Log[q2])/(L - R)^4 + (2*(6*L + 5*R)*Log[R])/
        (L - R)^4) - (2*PolyLog[2, 1 - L/R])/(L - R)^3 +
     (4*PolyLog[2, (L - R)/L])/(L - R)^3) +
   ((4*mu2^3*(2*L + R) - L*R*(4*L^2 + 3*L*R + 2*R^2) -
       3*mu2^2*(4*L^2 + 7*L*R + 2*R^2) + mu2*(4*L^3 + 17*L^2*R + 13*L*R^2 +
         2*R^3))/(2*L*(-L + mu2)*(mu2 - R)*R) - (2*L*(L - 2*mu2)*Log[L]^2)/
      (L - mu2)^2 + ((4*mu2^4 - 4*mu2^3*R - 2*L^2*R^2 +
        2*mu2^2*R*(-4*L + R) + 4*L*mu2*R*(L + R))*Log[mu2]^2)/
      ((L - mu2)^2*(mu2 - R)^2) - 2*Log[q2]^2 +
     ((-2*mu2*(mu2 - R)^2*R + L^2*(9*mu2^2 - 4*mu2*R + R^2) +
        L*(-13*mu2^3 + 14*mu2^2*R - 9*mu2*R^2 + 2*R^3))*Log[R])/
      (2*L*(L - mu2)*(mu2 - R)^2) + ((2*mu2 - R)*R*Log[R]^2)/(mu2 - R)^2 +
     Log[q2]*(-(((2*L + R)*(L - 2*mu2 + R))/(L*R)) + Log[R]) +
     Log[mu2]*((mu2*(-2*L^2*R^2*(L + R) - 2*mu2^4*(2*L + R) +
          L*mu2*R*(3*L^2 + 8*L*R - 2*R^2) + mu2^3*(8*L^2 + L*R + 4*R^2) -
          2*mu2^2*(2*L^3 + 4*L^2*R - L*R^2 + R^3)))/(L*(L - mu2)^2*
         (mu2 - R)^2*R) + ((2*L*mu2^3 - 3*mu2^4 - L*mu2^2*(L - 8*R) +
          2*L^2*R^2 - 4*L*mu2*R*(L + R))*Log[R])/((L - mu2)^2*(mu2 - R)^2)) +
     Log[L]*((4*L^3*(mu2 - R) + mu2^2*(15*mu2 - 13*R)*R +
         L^2*(-8*mu2^2 + 15*mu2*R - 5*R^2) + 2*L*mu2*(2*mu2^2 - 7*mu2*R +
           3*R^2))/(2*(L - mu2)^2*(mu2 - R)*R) +
       ((-5*mu2^4 - 2*mu2^3*(L - 4*R) + 2*L^2*R^2 - 4*L*mu2*R*(L + R) +
          mu2^2*(L^2 + 8*L*R - 4*R^2))*Log[mu2])/((L - mu2)^2*(mu2 - R)^2) +
       3*Log[q2] + ((2*mu2^4 - 2*mu2^3*R - L^2*R^2 + mu2^2*R*(-4*L + R) +
          2*L*mu2*R*(L + R))*Log[R])/((L - mu2)^2*(mu2 - R)^2)) +
     (2*(-L^2 + 2*L*mu2 + mu2^2)*PolyLog[2, 1 - L/mu2])/(L - mu2)^2 +
     (4*mu2^2*PolyLog[2, (-L + mu2)/mu2])/(L - mu2)^2 +
     (2*(mu2^2 + 2*mu2*R - R^2)*PolyLog[2, 1 - R/mu2])/(mu2 - R)^2 +
     Xtau^2*((4*L^2 - 8*L*mu2 + 6*L*R + 4*mu2*R - 2*R^2)/(L^2*R - L*R^2) +
       (2*L*(3*L^2 - 6*L*mu2 + 5*mu2^2 + 2*L*R - 4*mu2*R)*Log[L]^2)/
        ((L - mu2)^2*(L - R)^2) +
       ((L^3*(5*mu2 - 11*R) + 2*mu2*(mu2 - R)*R^2 +
          L^2*(-5*mu2^2 + 16*mu2*R + R^2) + L*(4*mu2^3 - 13*mu2^2*R +
            mu2*R^2 + 2*R^3))*Log[R])/(L*(L - mu2)*(L - R)^2*(mu2 - R)) +
       (2*R*(4*mu2^2 - 6*mu2*R + 3*R^2 + L*(-2*mu2 + R))*Log[R]^2)/
        ((L - R)^2*(mu2 - R)^2) +
       Log[L]*((-4*L^3*(mu2 - R) + L*R*(5*mu2^2 - 2*mu2*R - 7*R^2) -
           mu2*R*(4*mu2^2 - 7*mu2*R + R^2) + L^2*(4*mu2^2 - 13*mu2*R +
             11*R^2))/((L - mu2)*(L - R)^2*(mu2 - R)*R) +
         (2*(-2*mu2^5 + 2*L^3*R^2 - 2*L^2*mu2*R*(2*L + 3*R) +
            mu2^4*(3*L + 7*R) - 2*mu2^3*(2*L^2 + 5*L*R + 3*R^2) +
            mu2^2*(L^3 + 13*L^2*R + 4*L*R^2 + 2*R^3))*Log[mu2])/
          ((L - mu2)^2*(L - R)^2*(mu2 - R)^2) +
         ((-8*L + 4*mu2 + 2*R)*Log[q2])/(L - R)^2 +
         ((-6*L)/(L - R)^2 + (4*L*(L - 2*mu2))/((L - mu2)^2*(L - R)) -
           (4*R)/(L - R)^2 - (4*(L + R))/(L - R)^2 + (2*R*(-2*mu2 + R))/
            ((mu2 - R)^2*(-L + R)))*Log[R]) +
       Log[q2]*((4*L^2 - 8*L*mu2 + 4*L*R + 4*mu2*R - 2*R^2)/(L^2*R - L*R^2) +
         ((8*L - 2*(2*mu2 + R))*Log[R])/(L - R)^2) +
       Log[mu2]*((4*mu2*(L^2*(2*mu2 - R) + mu2*(mu2 - R)*R +
            L*mu2*(-2*mu2 + R)))/(L*(L - mu2)*(L - R)*(mu2 - R)*R) +
         (2*(2*mu2^5 - 2*L^3*R^2 + 2*L^2*mu2*R*(2*L + 3*R) -
            mu2^4*(3*L + 7*R) + 2*mu2^3*(2*L^2 + 5*L*R + 3*R^2) -
            mu2^2*(L^3 + 13*L^2*R + 4*L*R^2 + 2*R^3))*Log[R])/
          ((L - mu2)^2*(L - R)^2*(mu2 - R)^2)) +
       (4*(-L + mu2)*PolyLog[2, (-L + mu2)/mu2])/(L - R)^2 +
       (4*(L - mu2)*PolyLog[2, 1 - R/mu2])/(L - R)^2) +
     Xtau^4*((-2*L^4*(mu2 - R) + mu2*R^2*(2*mu2^2 - 3*mu2*R + R^2) +
         L^3*(6*mu2^2 - 36*mu2*R + 32*R^2) - L*R*(16*mu2^3 - 5*mu2^2*R -
           14*mu2*R^2 + R^3) - L^2*(4*mu2^3 - 46*mu2^2*R + 31*mu2*R^2 +
           15*R^3))/(L*(L - mu2)*(L - R)^3*(mu2 - R)*R) -
       (2*L*(5*L^3 + mu2*(5*mu2 - 2*R)*R + 5*L^2*(-2*mu2 + R) +
          L*(6*mu2^2 - 10*mu2*R + R^2))*Log[L]^2)/((L - mu2)^2*(L - R)^4) +
       ((-2*mu2*(mu2 - R)^2*R^3 + L^4*(-27*mu2^2 + 64*mu2*R - 33*R^2) +
          L^3*(47*mu2^3 - 162*mu2^2*R + 153*mu2*R^2 - 50*R^3) +
          L*R*(-12*mu2^4 - 19*mu2^3*R + 70*mu2^2*R^2 - 45*mu2*R^3 + 2*R^4) +
          L^2*(-24*mu2^4 + 118*mu2^3*R - 101*mu2^2*R^2 - 26*mu2*R^3 +
            45*R^4))*Log[R])/(2*L*(L - mu2)*(L - R)^4*(mu2 - R)^2) -
       (R*(L + R)*(8*mu2^2 - 14*mu2*R + 7*R^2 + L*(-2*mu2 + R))*Log[R]^2)/
        ((L - R)^4*(mu2 - R)^2) + Log[q2]*
        ((-2*L^3 + L^2*(4*mu2 - 21*R) + R^2*(-2*mu2 + R) + 10*L*R*(mu2 + R))/
          (L*(L - R)^3*R) + ((-13*L^2 + 12*L*mu2 - 6*L*R + 7*R^2)*Log[R])/
          (L - R)^4) + Log[L]*((4*L^5*(mu2 - R) + mu2^2*R^3*
            (-13*mu2 + 15*R) + L^4*(-8*mu2^2 + 77*mu2*R - 71*R^2) +
           2*L*mu2*R*(-18*mu2^3 + 4*mu2^2*R + 37*mu2*R^2 - 25*R^3) +
           2*L^3*(2*mu2^3 - 87*mu2^2*R + 81*mu2*R^2 + 6*R^3) +
           L^2*R*(145*mu2^3 - 123*mu2^2*R - 49*mu2*R^2 + 27*R^3))/
          (2*(L - mu2)^2*(L - R)^4*(mu2 - R)*R) +
         ((12*L*mu2^5 + 2*L^2*R^2*(-2*L^2 - 2*L*R + R^2) -
            mu2^4*(25*L^2 + 28*L*R + R^2) + 6*L*mu2^3*(3*L^2 + 10*L*R +
              3*R^2) - L*mu2^2*(3*L^3 + 44*L^2*R + 41*L*R^2 - 4*R^3) +
            4*L*mu2*R*(2*L^3 + 7*L^2*R + L*R^2 - R^3))*Log[mu2])/
          ((L - mu2)^2*(L - R)^4*(mu2 - R)^2) +
         ((13*L^2 - 7*R^2 + 6*L*(-2*mu2 + R))*Log[q2])/(L - R)^4 +
         ((10*L)/(L - R)^3 - 2/(L - R)^2 + (22*L*R)/(L - R)^4 +
           (6*R)/(-L + R)^3 - (2*L*(L - 2*mu2)*(L + R))/((L - mu2)^2*
             (L - R)^3) + ((2*mu2 - R)*R*(L + R))/((mu2 - R)^2*(-L + R)^3) +
           (4*(L + R)^2)/(L - R)^4)*Log[R]) +
       Log[mu2]*((-2*mu2*(-(mu2^2*(mu2 - R)^2*R^2) + L*mu2^2*R*
             (2*mu2^2 - 5*mu2*R + 2*R^2) + L^4*(2*mu2^2 - 3*mu2*R + 2*R^2) +
            L^3*(-4*mu2^3 + 8*mu2^2*R - 9*mu2*R^2 + 2*R^3) +
            L^2*(2*mu2^4 - 5*mu2^3*R + 7*mu2^2*R^2 - R^4)))/
          (L*(L - mu2)^2*(L - R)^3*(mu2 - R)^2*R) +
         ((-12*L*mu2^5 + 2*L^2*R^2*(2*L^2 + 2*L*R - R^2) +
            mu2^4*(25*L^2 + 28*L*R + R^2) - 6*L*mu2^3*(3*L^2 + 10*L*R +
              3*R^2) + L*mu2^2*(3*L^3 + 44*L^2*R + 41*L*R^2 - 4*R^3) -
            4*L*mu2*R*(2*L^3 + 7*L^2*R + L*R^2 - R^3))*Log[R])/
          ((L - mu2)^2*(L - R)^4*(mu2 - R)^2)) +
       ((4*L^2 + 6*mu2^2 - 2*R^2 + 4*L*(-3*mu2 + R))*
         PolyLog[2, (-L + mu2)/mu2])/(L - R)^4 +
       (2*(-2*L^2 + 6*L*mu2 - 3*mu2^2 - 2*L*R + R^2)*PolyLog[2, 1 - R/mu2])/
          (L - R)^4))/cbe^2) /. { delta[x_,y_,z_] :> TDelta[x,y,z], phi[x_,y_,z_] :> TPhi[x,y,z] }
          ];
(* Two loop alphas-alphab corrections, computation for general masses *)
lambda2LPhiHSSAlphaBAlphaSFull = With[{
    k = 1/(4*Pi)^2,
    sbe = Sqrt[TanBeta^2/(1+TanBeta^2)],
    cbe = Sqrt[1/(1+TanBeta^2)],
    Nc = 3, (* number of colors *)
    Q  = msq2[3,3]*(1+0.02),
    U  = msu2[3,3]*(1-0.02),
    mD  = msd2[3,3]*(1-0.015),
    mu2 = MuInput^2,
    q2 = SCALE^2, (* renormalization/matching scale *)
    g = M3Input^2,
    sgn = M3Input/Abs[M3Input],
    gt = Yu[3,3], (* SM Yukawa coupling *)
    gb = Yd[3,3],
    Xt = xt,
    Yt = yt,
    Xb = xb,
    Yb = yb,
    mQ3 = Sqrt[msq2[3,3]],
    mU3 = Sqrt[msu2[3,3]],
    mD3 = Sqrt[msd2[3,3]],
    M3 = M3Input,
    Mu = MuInput,
    Q2 = SCALE^2,
    Xtildet = xtt,
    Xtildeb = xbb,
    cosb = Cos[ArcTan[TanBeta]],
    sinb = Sin[ArcTan[TanBeta]],
    MA = mAInput,
    CF = 4/3
    },
    ybMSSM[mQ3_,mU3_,mD3_,M3_,Mu_,TanBeta_,Xt_,Xb_] := Module[{deltagsb, deltagbyL1, deltagbyL2, deltagbyL3, deltagbyL4},
          deltagsb = - g3^2*CF*k*(1+Log[M3^2/Q2]+TCF[6][mQ3/M3]+TCF[6][mD3/M3]-Xb/M3*TCF[9][mQ3/M3,mD3/M3]);
          deltagbyL1 = - gb^2/cosb^2*k*(3/4*Log[Mu^2/Q2]+3/8*sinb^2*(2*Log[MA^2/Q2]-1)+TCF[6][mQ3/Mu]+1/2*TCF[6][mD3/Mu]);
          deltagbyL2 = - gt^2/sinb^2*k*(1/4*Log[Mu^2/Q2]+1/8*cosb^2*(2*Log[MA^2/Q2]-1)+sinb^2*(Log[MA^2/Q2]-1));
          deltagbyL3 = - gt^2/sinb^2*k*(1/2*TCF[6][mU3/Mu]+(Xt*TanBeta)/(2*Mu)*TCF[9][mQ3/Mu,mU3/Mu] );
          deltagbyL4 = - 1/2*(-gt^2*Nc*k*Xtildet/6*TCF[5][xQU]-gb^2*Nc*k*Xtildeb/6*TCF[5][xQD]);
          (gb/(1-deltagsb-(deltagbyL1+deltagbyL2+deltagbyL3+deltagbyL4)))
    ];

    (
       g3^2*ybMSSM[mQ3, mU3, mD3, M3, Mu, TanBeta, Xt, Xb]^4*k^2*((8*(-3*mD^2*Q^2 + 2*g^3*(mD + Q) + 6*g*mD*Q*(mD + Q) - 
            g^2*(2*mD^2 + 9*mD*Q + 2*Q^2)))/((g - mD)*mD*(g - Q)*Q) + 
         (24*(2*g^2 - 2*g*mD + mD^2)*Log[mD])/(g - mD)^2 - 
         (8*(3*g^2 - 2*g*mD + mD^2)*Log[mD]^2)/(g - mD)^2 - 
         (8*(3*g^2 - 2*g*Q + Q^2)*Log[Q]^2)/(g - Q)^2 + 
         ((16*(-3*mD^2*Q^2 + g^3*(mD + Q) + 3*g*mD*Q*(mD + Q) - 
              g^2*(mD^2 + 3*mD*Q + Q^2)))/((g - mD)*mD*(g - Q)*Q) + 
           (16*(2*g^2 - 2*g*mD + mD^2)*Log[mD])/(g - mD)^2)*Log[q2] - 16*Log[q2]^2 + 
         Log[Q]*((24*(2*g^2 - 2*g*Q + Q^2))/(g - Q)^2 + 
           (16*(2*g^2 - 2*g*Q + Q^2)*Log[q2])/(g - Q)^2) + 
         Log[g]*((-8*g^2*(2*g^3*(mD + Q) + 2*g*(-mD + Q)^2*(mD + Q) + 
              g^2*(-4*mD^2 + 2*mD*Q - 4*Q^2) + mD*Q*(mD^2 + Q^2)))/
            ((g - mD)^2*mD*(g - Q)^2*Q) + (16*g^2*Log[mD])/(g - mD)^2 + 
           (16*g^2*Log[Q])/(g - Q)^2 - (16*g^2*(2*g^2 + mD^2 + Q^2 - 2*g*(mD + Q))*
             Log[q2])/((g - mD)^2*(g - Q)^2)) + 
         Xb^2*((-32*g)/(mD*Q) + (32*g^2*(g - mD - Q)*Log[g])/
            ((g - mD)*mD*(g - Q)*Q) + (32*(3*g - 2*mD)*Log[mD])/
            ((g - mD)*(mD - Q)) + (16*(-3*mD + Q)*Log[mD]^2)/(-mD + Q)^2 + 
           (16*(mD - 3*Q)*Log[Q]^2)/(-mD + Q)^2 + 
           ((-32*g)/(mD*Q) - (64*Log[mD])/(-mD + Q))*Log[q2] + 
           Log[Q]*((32*(3*g - 2*Q))/((g - Q)*(-mD + Q)) + (32*(mD + Q)*Log[mD])/
              (-mD + Q)^2 + (64*Log[q2])/(-mD + Q))) - 
         (32*g^2*PolyLog[2, 1 - g/mD])/(g - mD)^2 - (32*g^2*PolyLog[2, 1 - g/Q])/
          (g - Q)^2 + Xb^4*((16*(6*mD*Q + g*(mD + Q)))/(mD*Q*(-mD + Q)^2) - 
           (16*g*(mD + Q)*Log[g])/(mD*Q*(-mD + Q)^2) + 
           (16*(2*g + 7*mD + 3*Q)*Log[mD])/(-mD + Q)^3 + 
           (16*(g*(mD - Q) + 2*mD*(mD + Q))*Log[mD]^2)/(-mD + Q)^4 + 
           (16*(g*(-mD + Q) + 2*Q*(mD + Q))*Log[Q]^2)/(-mD + Q)^4 + 
           ((16*(4*mD*Q + g*(mD + Q)))/(mD*Q*(-mD + Q)^2) + 
             (32*(g + mD + Q)*Log[mD])/(-mD + Q)^3)*Log[q2] + 
           Log[Q]*((-16*(2*g + 3*mD + 7*Q))/(-mD + Q)^3 - (32*(mD + Q)^2*Log[mD])/
              (-mD + Q)^4 - (32*(g + mD + Q)*Log[q2])/(-mD + Q)^3) + 
           (16*(-2*g + mD + Q)*PolyLog[2, 1 - g/mD])/(-mD + Q)^3 + 
           (16*(2*g - mD - Q)*PolyLog[2, 1 - g/Q])/(-mD + Q)^3) + 
         Xb^3*((-256*Sqrt[g]*sgn)/(-mD + Q)^2 - (64*Sqrt[g]*(3*mD + Q)*sgn*Log[mD])/
            (-mD + Q)^3 + (32*Sqrt[g]*(2*g - mD - Q)*sgn*Log[mD]^2)/(mD - Q)^3 + 
           (32*Sqrt[g]*(2*g - mD - Q)*sgn*Log[Q]^2)/(-mD + Q)^3 + 
           Log[g]*((128*g^(3/2)*sgn*Log[mD])/(-mD + Q)^3 - (128*g^(3/2)*sgn*Log[Q])/
              (-mD + Q)^3) + ((-128*Sqrt[g]*sgn)/(-mD + Q)^2 - 
             (64*Sqrt[g]*(mD + Q)*sgn*Log[mD])/(-mD + Q)^3)*Log[q2] + 
           Log[Q]*((64*Sqrt[g]*(mD + 3*Q)*sgn)/(-mD + Q)^3 + 
             (64*Sqrt[g]*(mD + Q)*sgn*Log[q2])/(-mD + Q)^3) + 
           (64*Sqrt[g]*(2*g - mD - Q)*sgn*PolyLog[2, 1 - g/mD])/(mD - Q)^3 + 
           (64*Sqrt[g]*(2*g - mD - Q)*sgn*PolyLog[2, 1 - g/Q])/(-mD + Q)^3) + 
         Xb*((64*Sqrt[g]*sgn*Log[mD])/(-mD + Q) + (64*g^(3/2)*sgn*Log[mD]^2)/
            ((g - mD)*(mD - Q)) + (64*g^(3/2)*sgn*Log[Q]^2)/((g - Q)*(-mD + Q)) - 
           (64*g^(3/2)*sgn*Log[mD]*Log[q2])/((g - mD)*(mD - Q)) + 
           Log[g]*((-64*g^(3/2)*sgn*Log[mD])/((g - mD)*(mD - Q)) - 
             (64*g^(3/2)*sgn*Log[Q])/((g - Q)*(-mD + Q)) + (64*g^(3/2)*sgn*Log[q2])/
              ((g - mD)*(g - Q))) + Log[Q]*((-64*Sqrt[g]*sgn)/(-mD + Q) - 
             (64*g^(3/2)*sgn*Log[q2])/((g - Q)*(-mD + Q))) + 
           (128*g^(3/2)*sgn*PolyLog[2, 1 - g/mD])/((g - mD)*(mD - Q)) + 
           (128*g^(3/2)*sgn*PolyLog[2, 1 - g/Q])/((g - Q)*(-mD + Q))))
      ) /. { delta[x_,y_,z_] :> TDelta[x,y,z], phi[x_,y_,z_] :> TPhi[x,y,z] }
      ];

(* Two loop alphas-alphab corrections, computation for degenerate masses *)
lambda2LPhiHSSAlphaBAlphaSDegenerate = With[{
    k = 1/(4*Pi)^2,
    sbe = Sqrt[TanBeta^2/(1+TanBeta^2)],
    cbe = Sqrt[1/(1+TanBeta^2)],
    Nc = 3, (* number of colors *)
    q2 = SCALE^2, (* renormalization/matching scale *)
    MS2 = msq2[3,3],
    gb = Yd[3,3],
    Xb = xb,
    mQ3 = Sqrt[msq2[3,3]],
    mU3 = Sqrt[msu2[3,3]],
    mD3 = Sqrt[msd2[3,3]],
    M3 = M3Input,
    Mu = MuInput,
    Q2 = SCALE^2,
    gt = Yu[3,3],
    Xt = xt,
    Xtildet = xtt,
    Xtildeb = xbb,
    cosb = Cos[ArcTan[TanBeta]],
    sinb = Sin[ArcTan[TanBeta]],
    MA = mAInput,
    CF = 4/3
  },
  ybMSSM[mQ3_,mU3_,mD3_,M3_,Mu_,TanBeta_,Xt_,Xb_] := Module[{deltagsb, deltagbyL1, deltagbyL2, deltagbyL3, deltagbyL4},
          deltagsb = - g3^2*CF*k*(1+Log[M3^2/Q2]+TCF[6][mQ3/M3]+TCF[6][mD3/M3]-Xb/M3*TCF[9][mQ3/M3,mD3/M3]);
          deltagbyL1 = - gb^2/cosb^2*k*(3/4*Log[Mu^2/Q2]+3/8*sinb^2*(2*Log[MA^2/Q2]-1)+TCF[6][mQ3/Mu]+1/2*TCF[6][mD3/Mu]);
          deltagbyL2 = - gt^2/sinb^2*k*(1/4*Log[Mu^2/Q2]+1/8*cosb^2*(2*Log[MA^2/Q2]-1)+sinb^2*(Log[MA^2/Q2]-1));
          deltagbyL3 = - gt^2/sinb^2*k*(1/2*TCF[6][mU3/Mu]+(Xt*TanBeta)/(2*Mu)*TCF[9][mQ3/Mu,mU3/Mu] );
          deltagbyL4 = - 1/2*(-gt^2*Nc*k*Xtildet/6*TCF[5][xQU]-gb^2*Nc*k*Xtildeb/6*TCF[5][xQD]);
          (gb/(1-deltagsb-(deltagbyL1+deltagbyL2+deltagbyL3+deltagbyL4)))
    ];
  g3^2*(ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb])^4*k^2*((-32*Xb)/Sqrt[MS2] + (16*Xb^2)/MS2 + (16*Xb^3)/(3*MS2^(3/2)) -
  (4*Xb^4)/(3*MS2^2) + 32*Log[MS2] + (32*Xb*Log[MS2])/Sqrt[MS2] -
  (32*Xb^2*Log[MS2])/MS2 - (32*Xb^3*Log[MS2])/(3*MS2^(3/2)) - 16*Log[MS2]^2 -
  32*Log[q2] - (32*Xb*Log[q2])/Sqrt[MS2] + (32*Xb^2*Log[q2])/MS2 +
  (32*Xb^3*Log[q2])/(3*MS2^(3/2)) + 32*Log[MS2]*Log[q2] - 16*Log[q2]^2)
];

(* Two loop alphas-alphab corrections, computation for degenerate squark *)
lambda2LPhiHSSAlphaBAlphaSDegenerateSquark = With[{
    k = 1/(4*Pi)^2,
    sbe = Sqrt[TanBeta^2/(1+TanBeta^2)],
    cbe = Sqrt[1/(1+TanBeta^2)],
    Nc = 3, (* number of colors *)
    q2 = SCALE^2, (* renormalization/matching scale *)
    MS2 = msq2[3,3],
    gb = Yd[3,3],
    Xb = xb,      
    mQ3 = Sqrt[msq2[3,3]],
    mU3 = Sqrt[msu2[3,3]],
    mD3 = Sqrt[msd2[3,3]],
    M3 = M3Input,
    g = M3Input^2,
    sgn = M3Input/Abs[M3Input],  
    Mu = MuInput,
    Q2 = SCALE^2,
    gt = Yu[3,3],
    Xt = xt,
    Xtildet = xtt,
    Xtildeb = xbb,
    cosb = Cos[ArcTan[TanBeta]],
    sinb = Sin[ArcTan[TanBeta]],
    MA = mAInput,
    CF = 4/3
  },
  ybMSSM[mQ3_,mU3_,mD3_,M3_,Mu_,TanBeta_,Xt_,Xb_] := Module[{deltagsb, deltagbyL1, deltagbyL2, deltagbyL3, deltagbyL4},
          deltagsb = - g3^2*CF*k*(1+Log[M3^2/Q2]+TCF[6][mQ3/M3]+TCF[6][mD3/M3]-Xb/M3*TCF[9][mQ3/M3,mD3/M3]);
          deltagbyL1 = - gb^2/cosb^2*k*(3/4*Log[Mu^2/Q2]+3/8*sinb^2*(2*Log[MA^2/Q2]-1)+TCF[6][mQ3/Mu]+1/2*TCF[6][mD3/Mu]);
          deltagbyL2 = - gt^2/sinb^2*k*(1/4*Log[Mu^2/Q2]+1/8*cosb^2*(2*Log[MA^2/Q2]-1)+sinb^2*(Log[MA^2/Q2]-1));
          deltagbyL3 = - gt^2/sinb^2*k*(1/2*TCF[6][mU3/Mu]+(Xt*TanBeta)/(2*Mu)*TCF[9][mQ3/Mu,mU3/Mu] );
          deltagbyL4 = - 1/2*(-gt^2*Nc*k*Xtildet/6*TCF[5][xQU]-gb^2*Nc*k*Xtildeb/6*TCF[5][xQD]);
          (gb/(1-deltagsb-(deltagbyL1+deltagbyL2+deltagbyL3+deltagbyL4)))
      ];
(-8*g3^2*ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^4*k^2*(6*MS2^3*(3*g^2 - 2*g*MS2 + MS2^2 - 4*g^(3/2)*sgn*Xb)*
    Log[MS2]^2 + g^(3/2)*Log[g]*(6*Sqrt[g]*MS2^2*(2*g + MS2) - 
     24*MS2^2*(g + MS2)*sgn*Xb - 12*Sqrt[g]*(g - 2*MS2)*MS2*Xb^2 + 
     4*(g - 2*MS2)*MS2*sgn*Xb^3 + Sqrt[g]*(2*g - 3*MS2)*Xb^4 - 
     12*MS2^3*(Sqrt[g] - 2*sgn*Xb)*(Log[MS2] - Log[q2])) + 
   MS2*Log[MS2]*(-18*MS2^2*(2*g^2 - 2*g*MS2 + MS2^2) + 
     48*g^(3/2)*MS2^2*sgn*Xb + 12*(g - 2*MS2)*(2*g - MS2)*MS2*Xb^2 + 
     4*Sqrt[g]*MS2^2*sgn*Xb^3 + (-3*g^2 + 6*g*MS2 - 2*MS2^2)*Xb^4 - 
     12*MS2^2*(2*g^2 - 2*g*MS2 + MS2^2 - 2*g^(3/2)*sgn*Xb)*Log[q2]) - 
   (g - MS2)*(3*MS2^2*(4*g^2 - 9*g*MS2 + 3*MS2^2) + 
     24*Sqrt[g]*MS2^2*(-g + MS2)*sgn*Xb - 12*MS2*(g^2 - 3*g*MS2 + MS2^2)*
      Xb^2 + 4*Sqrt[g]*(g - 2*MS2)*MS2*sgn*Xb^3 + (2*g^2 - 4*g*MS2 + MS2^2)*
      Xb^4 + 2*Log[q2]*(3*MS2^2*(2*g^2 - 3*g*MS2 + 3*MS2^2) - 
       12*g^(3/2)*MS2^2*sgn*Xb - 6*(g - 2*MS2)*(g - MS2)*MS2*Xb^2 + 
       2*Sqrt[g]*(g - MS2)*MS2*sgn*Xb^3 + (g - MS2)^2*Xb^4 + 
       3*MS2^3*(-g + MS2)*Log[q2])) + 24*g^(3/2)*MS2^3*(Sqrt[g] - 2*sgn*Xb)*
    PolyLog[2, 1 - g/MS2]))/(3*(g - MS2)^2*MS2^3)
];

lambda3LATASASDegenerate = With[{
    k = 1/(4*Pi)^2,
    gt = Yu[3,3],
    MR = SCALE,
    MS = Sqrt[Sqrt[Abs[msq2[3,3] msu2[3,3]]]],
    Xtt = xt/Sqrt[Sqrt[Abs[msq2[3,3] msu2[3,3]]]]
    },
    (1/2)*k^3*(gt^4 g3^4) (
    (* 16/27 Xtt^3 (3568 - 1664 Log[MS^2/MR^2] + 444 Log[MS^2/MR^2]^2 - 2259 Zeta[3]) +  *)
     8/9 Xtt^2 (149 - 448 Log[MS^2/MR^2] + 912 Log[MS^2/MR^2]^2 - 990 Zeta[3]) - 
     64/27 Xtt (823 - 200 Log[MS^2/MR^2] + 954 Log[MS^2/MR^2]^2 - 477 Zeta[3]) - 
     4/135 (176 \[Pi]^4 + 960 \[Pi]^2 Log[2]^2 + 9900 Log[MS^2/MR^2]^2 - 24840 Log[MS^2/MR^2]^3 - 
        5 (4577 + 192 Log[2]^4 + 4608 PolyLog[4, 1/2] - 9864 Zeta[3]) + 
        10 Log[MS^2/MR^2] (-877 + 432 Zeta[3]))
    ) /. { a:PolyLog[__] :> N[a], a:Zeta[__] :> N[a] }
];

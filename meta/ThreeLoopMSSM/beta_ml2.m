{(-6*g1^2*M1^2)/5 - 6*g2^2*M2^2 + 2*MatMul[Adj[he], he] + 
  2*mh1*MatMul[Adj[Ye], Ye] + MatMul[ml, Adj[Ye], Ye] + 
  2*MatMul[Adj[Ye], me, Ye] + MatMul[Adj[Ye], Ye, ml], 
 (621*g1^4*M1^2)/25 + (18*g1^2*g2^2*M1^2)/5 + (18*g1^2*g2^2*M1*M2)/5 + 
  (18*g1^2*g2^2*M2^2)/5 + 33*g2^4*M2^2 + (12*g1^2*MatMul[Adj[he], he])/5 - 
  (12*g1^2*M1*MatMul[Adj[he], Ye])/5 - (12*g1^2*M1*MatMul[Adj[Ye], he])/5 + 
  (24*g1^2*M1^2*MatMul[Adj[Ye], Ye])/5 + (12*g1^2*mh1*MatMul[Adj[Ye], Ye])/
   5 + (6*g1^2*MatMul[ml, Adj[Ye], Ye])/5 + (12*g1^2*MatMul[Adj[Ye], me, Ye])/
   5 + (6*g1^2*MatMul[Adj[Ye], Ye, ml])/5 - 
  4*MatMul[Adj[he], he, Adj[Ye], Ye] - 4*MatMul[Adj[he], Ye, Adj[Ye], he] - 
  4*MatMul[Adj[Ye], he, Adj[he], Ye] - 4*MatMul[Adj[Ye], Ye, Adj[he], he] - 
  8*mh1*MatMul[Adj[Ye], Ye, Adj[Ye], Ye] - 
  2*MatMul[ml, Adj[Ye], Ye, Adj[Ye], Ye] - 
  4*MatMul[Adj[Ye], me, Ye, Adj[Ye], Ye] - 
  4*MatMul[Adj[Ye], Ye, ml, Adj[Ye], Ye] - 
  4*MatMul[Adj[Ye], Ye, Adj[Ye], me, Ye] - 
  2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, ml] + 
  g2^4*(3*mh1 + 3*mh2 + 3*trace[ml] + 9*trace[mq]) + 
  g1^4*((9*mh1)/25 + (9*mh2)/25 + (6*trace[mb])/25 + (18*trace[me])/25 + 
    (9*trace[ml])/25 + (3*trace[mq])/25 + (24*trace[mt])/25) + 
  MatMul[Adj[he], Ye]*(-6*trace[hb, Adj[Yb]] - 2*trace[he, Adj[Ye]]) + 
  MatMul[Adj[Ye], he]*(-6*trace[Adj[hb], Yb] - 2*trace[Adj[he], Ye]) + 
  MatMul[Adj[he], he]*(-6*trace[Adj[Yb], Yb] - 2*trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Ye], me, Ye]*(-6*trace[Adj[Yb], Yb] - 2*trace[Adj[Ye], Ye]) + 
  MatMul[ml, Adj[Ye], Ye]*(-3*trace[Adj[Yb], Yb] - trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Ye], Ye, ml]*(-3*trace[Adj[Yb], Yb] - trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Ye], Ye]*(-6*trace[hb, Adj[hb]] - 2*trace[he, Adj[he]] - 
    12*mh1*trace[Adj[Yb], Yb] - 4*mh1*trace[Adj[Ye], Ye] - 
    6*trace[Yb, Adj[Yb], mb] - 2*trace[Ye, Adj[Ye], me] - 
    6*trace[Adj[Yb], Yb, mq] - 2*trace[Adj[Ye], Ye, ml]), 
 g2^6*(-3*mh1 - 3*mh2 - 3*trace[ml] - 9*trace[mq]) + 
  g1^2*g2^4*((-9*mh1)/5 - (9*mh2)/5 - (9*trace[ml])/5 - (27*trace[mq])/5) + 
  g1^6*((-621*mh1)/125 - (621*mh2)/125 - (414*trace[mb])/125 - 
    (1242*trace[me])/125 - (621*trace[ml])/125 - (207*trace[mq])/125 - 
    (1656*trace[mt])/125) + g1^4*g2^2*((-27*mh1)/25 - (27*mh2)/25 - 
    (18*trace[mb])/25 - (54*trace[me])/25 - (27*trace[ml])/25 - 
    (9*trace[mq])/25 - (72*trace[mt])/25) + MatMul[Adj[he], Ye, Adj[Ye], Ye]*
   (12*trace[hb, Adj[Yb]] + 4*trace[he, Adj[Ye]]) + 
  MatMul[Adj[Ye], Ye, Adj[he], Ye]*(12*trace[hb, Adj[Yb]] + 
    4*trace[he, Adj[Ye]]) + g2^2*MatMul[Adj[he], Ye]*
   (36*trace[hb, Adj[Yb]] + 12*trace[he, Adj[Ye]]) + 
  g2^2*M2*MatMul[Adj[Ye], Ye]*(-36*trace[hb, Adj[Yb]] - 
    12*trace[he, Adj[Ye]] - 36*trace[Adj[hb], Yb] - 12*trace[Adj[he], Ye]) + 
  MatMul[Adj[Ye], he, Adj[Ye], Ye]*(12*trace[Adj[hb], Yb] + 
    4*trace[Adj[he], Ye]) + MatMul[Adj[Ye], Ye, Adj[Ye], he]*
   (12*trace[Adj[hb], Yb] + 4*trace[Adj[he], Ye]) + 
  g2^2*MatMul[Adj[Ye], he]*(36*trace[Adj[hb], Yb] + 12*trace[Adj[he], Ye]) + 
  g1^4*M1*((42*trace[hb, Adj[Yb]])/5 + (54*trace[he, Adj[Ye]])/5 + 
    (78*trace[ht, Adj[Yt]])/5 + (42*trace[Adj[hb], Yb])/5 + 
    (54*trace[Adj[he], Ye])/5 + (78*trace[Adj[ht], Yt])/5) + 
  g2^4*M2*(90*trace[hb, Adj[Yb]] + 30*trace[he, Adj[Ye]] + 
    90*trace[ht, Adj[Yt]] + 90*trace[Adj[hb], Yb] + 30*trace[Adj[he], Ye] + 
    90*trace[Adj[ht], Yt]) + g2^2*M2*MatMul[Adj[he], Ye]*
   (-36*trace[Adj[Yb], Yb] - 12*trace[Adj[Ye], Ye]) + 
  g2^2*M2*MatMul[Adj[Ye], he]*(-36*trace[Adj[Yb], Yb] - 
    12*trace[Adj[Ye], Ye]) + MatMul[ml, Adj[Ye], Ye, Adj[Ye], Ye]*
   (6*trace[Adj[Yb], Yb] + 2*trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Ye], Ye, Adj[Ye], Ye, ml]*(6*trace[Adj[Yb], Yb] + 
    2*trace[Adj[Ye], Ye]) + MatMul[Adj[he], he, Adj[Ye], Ye]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[Adj[he], Ye, Adj[Ye], he]*(12*trace[Adj[Yb], Yb] + 
    4*trace[Adj[Ye], Ye]) + MatMul[Adj[Ye], he, Adj[he], Ye]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Ye], Ye, Adj[he], he]*(12*trace[Adj[Yb], Yb] + 
    4*trace[Adj[Ye], Ye]) + MatMul[Adj[Ye], me, Ye, Adj[Ye], Ye]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Ye], Ye, ml, Adj[Ye], Ye]*(12*trace[Adj[Yb], Yb] + 
    4*trace[Adj[Ye], Ye]) + MatMul[Adj[Ye], Ye, Adj[Ye], me, Ye]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  g2^2*MatMul[ml, Adj[Ye], Ye]*(18*trace[Adj[Yb], Yb] + 
    6*trace[Adj[Ye], Ye]) + g2^2*MatMul[Adj[Ye], Ye, ml]*
   (18*trace[Adj[Yb], Yb] + 6*trace[Adj[Ye], Ye]) + 
  g2^2*MatMul[Adj[he], he]*(36*trace[Adj[Yb], Yb] + 12*trace[Adj[Ye], Ye]) + 
  g2^2*MatMul[Adj[Ye], me, Ye]*(36*trace[Adj[Yb], Yb] + 
    12*trace[Adj[Ye], Ye]) + g2^2*M2^2*MatMul[Adj[Ye], Ye]*
   (72*trace[Adj[Yb], Yb] + 24*trace[Adj[Ye], Ye]) + 
  g2^4*M2^2*(-270*trace[Adj[Yb], Yb] - 90*trace[Adj[Ye], Ye] - 
    270*trace[Adj[Yt], Yt]) + g1^4*M1^2*((-126*trace[Adj[Yb], Yb])/5 - 
    (162*trace[Adj[Ye], Ye])/5 - (234*trace[Adj[Yt], Yt])/5) + 
  MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*(12*trace[hb, Adj[hb]] + 
    4*trace[he, Adj[he]] + 36*mh1*trace[Adj[Yb], Yb] + 
    12*mh1*trace[Adj[Ye], Ye] + 12*trace[Yb, Adj[Yb], mb] + 
    4*trace[Ye, Adj[Ye], me] + 12*trace[Adj[Yb], Yb, mq] + 
    4*trace[Adj[Ye], Ye, ml]) + g2^2*MatMul[Adj[Ye], Ye]*
   (36*trace[hb, Adj[hb]] + 12*trace[he, Adj[he]] + 
    72*mh1*trace[Adj[Yb], Yb] + 24*mh1*trace[Adj[Ye], Ye] + 
    36*trace[Yb, Adj[Yb], mb] + 12*trace[Ye, Adj[Ye], me] + 
    36*trace[Adj[Yb], Yb, mq] + 12*trace[Adj[Ye], Ye, ml]) + 
  g2^4*(-54*trace[hb, Adj[hb]] - 18*trace[he, Adj[he]] - 
    54*trace[ht, Adj[ht]] - 54*mh1*trace[Adj[Yb], Yb] - 
    18*mh1*trace[Adj[Ye], Ye] - 54*mh2*trace[Adj[Yt], Yt] - 
    54*trace[Yb, Adj[Yb], mb] - 18*trace[Ye, Adj[Ye], me] - 
    54*trace[Yt, Adj[Yt], mt] - 54*trace[Adj[Yb], Yb, mq] - 
    18*trace[Adj[Ye], Ye, ml] - 54*trace[Adj[Yt], Yt, mq]) + 
  g1^4*((-126*trace[hb, Adj[hb]])/25 - (162*trace[he, Adj[he]])/25 - 
    (234*trace[ht, Adj[ht]])/25 - (126*mh1*trace[Adj[Yb], Yb])/25 - 
    (162*mh1*trace[Adj[Ye], Ye])/25 - (234*mh2*trace[Adj[Yt], Yt])/25 - 
    (126*trace[Yb, Adj[Yb], mb])/25 - (162*trace[Ye, Adj[Ye], me])/25 - 
    (234*trace[Yt, Adj[Yt], mt])/25 - (126*trace[Adj[Yb], Yb, mq])/25 - 
    (162*trace[Adj[Ye], Ye, ml])/25 - (234*trace[Adj[Yt], Yt, mq])/25) + 
  MatMul[Adj[he], Ye]*(-36*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
    12*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb] - 12*trace[hb, Adj[Yb]]*
     trace[Adj[Ye], Ye] - 4*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] + 
    12*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 72*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
    24*trace[Adj[Ye], he, Adj[Ye], Ye] + 
    12*trace[Adj[Yt], ht, Adj[Yb], Yb]) + MatMul[Adj[Ye], he]*
   (-36*trace[Adj[hb], Yb]*trace[Adj[Yb], Yb] - 12*trace[Adj[he], Ye]*
     trace[Adj[Yb], Yb] - 12*trace[Adj[hb], Yb]*trace[Adj[Ye], Ye] - 
    4*trace[Adj[he], Ye]*trace[Adj[Ye], Ye] + 
    12*trace[Adj[ht], Yt, Adj[Yb], Yb] + 72*trace[Adj[Yb], Yb, Adj[hb], Yb] + 
    24*trace[Adj[Ye], Ye, Adj[he], Ye] + 
    12*trace[Adj[Yt], Yt, Adj[hb], Yb]) + MatMul[ml, Adj[Ye], Ye]*
   (-9*trace[Adj[Yb], Yb]^2 - 6*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
    trace[Adj[Ye], Ye]^2 + 18*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    6*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 6*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + 
  MatMul[Adj[Ye], Ye, ml]*(-9*trace[Adj[Yb], Yb]^2 - 
    6*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - trace[Adj[Ye], Ye]^2 + 
    18*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 6*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    6*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + MatMul[Adj[he], he]*
   (-18*trace[Adj[Yb], Yb]^2 - 12*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
    2*trace[Adj[Ye], Ye]^2 + 36*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    12*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    12*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + MatMul[Adj[Ye], me, Ye]*
   (-18*trace[Adj[Yb], Yb]^2 - 12*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
    2*trace[Adj[Ye], Ye]^2 + 36*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    12*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    12*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + MatMul[Adj[Ye], Ye]*
   (-36*trace[hb, Adj[Yb]]*trace[Adj[hb], Yb] - 12*trace[he, Adj[Ye]]*
     trace[Adj[hb], Yb] - 12*trace[hb, Adj[Yb]]*trace[Adj[he], Ye] - 
    4*trace[he, Adj[Ye]]*trace[Adj[he], Ye] - 36*trace[hb, Adj[hb]]*
     trace[Adj[Yb], Yb] - 12*trace[he, Adj[he]]*trace[Adj[Yb], Yb] - 
    54*mh1*trace[Adj[Yb], Yb]^2 - 12*trace[hb, Adj[hb]]*trace[Adj[Ye], Ye] - 
    4*trace[he, Adj[he]]*trace[Adj[Ye], Ye] - 36*mh1*trace[Adj[Yb], Yb]*
     trace[Adj[Ye], Ye] - 6*mh1*trace[Adj[Ye], Ye]^2 - 
    36*trace[Adj[Yb], Yb]*trace[Yb, Adj[Yb], mb] - 
    12*trace[Adj[Ye], Ye]*trace[Yb, Adj[Yb], mb] - 
    12*trace[Adj[Yb], Yb]*trace[Ye, Adj[Ye], me] - 
    4*trace[Adj[Ye], Ye]*trace[Ye, Adj[Ye], me] - 
    36*trace[Adj[Yb], Yb]*trace[Adj[Yb], Yb, mq] - 
    12*trace[Adj[Ye], Ye]*trace[Adj[Yb], Yb, mq] - 
    12*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye, ml] - 
    4*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, ml] + 
    12*trace[hb, Adj[ht], Yt, Adj[Yb]] + 72*trace[Adj[hb], hb, Adj[Yb], Yb] + 
    12*trace[Adj[hb], hb, Adj[Yt], Yt] + 24*trace[Adj[he], he, Adj[Ye], Ye] + 
    12*trace[Adj[ht], ht, Adj[Yb], Yb] + 72*trace[Adj[Yb], hb, Adj[hb], Yb] + 
    108*mh1*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    24*trace[Adj[Ye], he, Adj[he], Ye] + 
    36*mh1*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    12*trace[Adj[Yt], ht, Adj[hb], Yb] + 
    24*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    12*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    72*trace[Yb, Adj[Yb], Yb, Adj[Yb], mb] + 
    12*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
    24*trace[Ye, Adj[Ye], Ye, Adj[Ye], me] + 
    12*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
    72*trace[Adj[Yb], Yb, Adj[Yb], Yb, mq] + 
    12*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
    24*trace[Adj[Ye], Ye, Adj[Ye], Ye, ml] + 
    12*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq]) + 
  g2^4*g3^2*M2^2*(1080 - 1296*Zeta[3]) + g2^4*g3^2*M2*M3*
   (720 - 864*Zeta[3]) + g2^4*g3^2*M3^2*(432 - 432*Zeta[3]) + 
  g2^4*M2^2*MatMul[Adj[Ye], Ye]*(-135 - 378*Zeta[3]) + 
  g1^4*g3^2*M1^2*(792/5 - (4752*Zeta[3])/25) + 
  g1^6*M1^2*(55767/125 - (21492*Zeta[3])/125) + 
  g1^4*g3^2*M1*M3*(528/5 - (3168*Zeta[3])/25) + 
  g1^2*g2^4*M2^2*(171/5 - (486*Zeta[3])/5) + 
  g1^2*g2^4*M1*M2*(18 - (324*Zeta[3])/5) + 
  g1^4*g3^2*M3^2*(1584/25 - (1584*Zeta[3])/25) + 
  g2^4*MatMul[Adj[he], he]*(-45/2 - 63*Zeta[3]) + 
  g2^4*MatMul[Adj[Ye], me, Ye]*(-45/2 - 63*Zeta[3]) + 
  g1^4*g2^2*M1^2*(81/25 - (1458*Zeta[3])/25) + 
  g1^2*M1^2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*(36 - (216*Zeta[3])/5) + 
  g1^4*g2^2*M1*M2*(54/25 - (972*Zeta[3])/25) + 
  g2^2*M2*MatMul[Adj[he], Ye, Adj[Ye], Ye]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Ye], he, Adj[Ye], Ye]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Ye], Ye, Adj[he], Ye]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Ye], Ye, Adj[Ye], he]*(6 - 36*Zeta[3]) + 
  g1^2*g2^4*M1^2*(72/5 - (162*Zeta[3])/5) + g1^2*g2^2*M1*MatMul[Adj[he], Ye]*
   (81/5 - (162*Zeta[3])/5) + g1^2*g2^2*M2*MatMul[Adj[he], Ye]*
   (81/5 - (162*Zeta[3])/5) + g1^2*g2^2*M1*MatMul[Adj[Ye], he]*
   (81/5 - (162*Zeta[3])/5) + g1^2*g2^2*M2*MatMul[Adj[Ye], he]*
   (81/5 - (162*Zeta[3])/5) + g2^4*MatMul[ml, Adj[Ye], Ye]*
   (-45/4 - (63*Zeta[3])/2) + g2^4*MatMul[Adj[Ye], Ye, ml]*
   (-45/4 - (63*Zeta[3])/2) + g1^2*MatMul[Adj[he], he, Adj[Ye], Ye]*
   (18 - (108*Zeta[3])/5) + g1^2*MatMul[Adj[he], Ye, Adj[Ye], he]*
   (18 - (108*Zeta[3])/5) + g1^2*MatMul[Adj[Ye], he, Adj[he], Ye]*
   (18 - (108*Zeta[3])/5) + g1^2*MatMul[Adj[Ye], Ye, Adj[he], he]*
   (18 - (108*Zeta[3])/5) + g1^2*MatMul[Adj[Ye], me, Ye, Adj[Ye], Ye]*
   (18 - (108*Zeta[3])/5) + g1^2*MatMul[Adj[Ye], Ye, ml, Adj[Ye], Ye]*
   (18 - (108*Zeta[3])/5) + g1^2*MatMul[Adj[Ye], Ye, Adj[Ye], me, Ye]*
   (18 - (108*Zeta[3])/5) + g1^4*g2^2*M2^2*(108/25 - (486*Zeta[3])/25) + 
  g1^2*MatMul[ml, Adj[Ye], Ye, Adj[Ye], Ye]*(9 - (54*Zeta[3])/5) + 
  g1^2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, ml]*(9 - (54*Zeta[3])/5) + 
  g1^4*M1*MatMul[Adj[he], Ye]*(549/5 - (162*Zeta[3])/25) + 
  g1^4*M1*MatMul[Adj[Ye], he]*(549/5 - (162*Zeta[3])/25) + 
  12*MatMul[Adj[he], he, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
  12*MatMul[Adj[he], Ye, Adj[Ye], he, Adj[Ye], Ye]*Zeta[3] + 
  12*MatMul[Adj[he], Ye, Adj[Ye], Ye, Adj[Ye], he]*Zeta[3] + 
  12*MatMul[Adj[Ye], he, Adj[he], Ye, Adj[Ye], Ye]*Zeta[3] + 
  12*MatMul[Adj[Ye], he, Adj[Ye], Ye, Adj[he], Ye]*Zeta[3] + 
  12*MatMul[Adj[Ye], Ye, Adj[he], he, Adj[Ye], Ye]*Zeta[3] + 
  12*MatMul[Adj[Ye], Ye, Adj[he], Ye, Adj[Ye], he]*Zeta[3] + 
  12*MatMul[Adj[Ye], Ye, Adj[Ye], he, Adj[he], Ye]*Zeta[3] + 
  12*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, Adj[he], he]*Zeta[3] + 
  36*mh1*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
  6*MatMul[ml, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
  12*MatMul[Adj[Ye], me, Ye, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
  12*MatMul[Adj[Ye], Ye, ml, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
  12*MatMul[Adj[Ye], Ye, Adj[Ye], me, Ye, Adj[Ye], Ye]*Zeta[3] + 
  12*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, ml, Adj[Ye], Ye]*Zeta[3] + 
  12*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], me, Ye]*Zeta[3] + 
  6*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye, ml]*Zeta[3] + 
  g1^4*MatMul[ml, Adj[Ye], Ye]*(-549/20 + (81*Zeta[3])/50) + 
  g1^4*MatMul[Adj[Ye], Ye, ml]*(-549/20 + (81*Zeta[3])/50) + 
  g1^4*MatMul[Adj[he], he]*(-549/10 + (81*Zeta[3])/25) + 
  g1^4*MatMul[Adj[Ye], me, Ye]*(-549/10 + (81*Zeta[3])/25) + 
  g1^2*g2^2*MatMul[ml, Adj[Ye], Ye]*(-81/10 + (81*Zeta[3])/5) + 
  g1^2*g2^2*MatMul[Adj[Ye], Ye, ml]*(-81/10 + (81*Zeta[3])/5) + 
  g2^2*MatMul[ml, Adj[Ye], Ye, Adj[Ye], Ye]*(-3 + 18*Zeta[3]) + 
  g2^2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, ml]*(-3 + 18*Zeta[3]) + 
  g1^4*M1^2*MatMul[Adj[Ye], Ye]*(-1647/5 + (486*Zeta[3])/25) + 
  g1^2*M1*MatMul[Adj[he], Ye, Adj[Ye], Ye]*(-18 + (108*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Ye], he, Adj[Ye], Ye]*(-18 + (108*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Ye], Ye, Adj[he], Ye]*(-18 + (108*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Ye], Ye, Adj[Ye], he]*(-18 + (108*Zeta[3])/5) + 
  g1^2*g2^2*MatMul[Adj[he], he]*(-81/5 + (162*Zeta[3])/5) + 
  g1^2*g2^2*MatMul[Adj[Ye], me, Ye]*(-81/5 + (162*Zeta[3])/5) + 
  g2^2*MatMul[Adj[he], he, Adj[Ye], Ye]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[he], Ye, Adj[Ye], he]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Ye], he, Adj[he], Ye]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Ye], Ye, Adj[he], he]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Ye], me, Ye, Adj[Ye], Ye]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Ye], Ye, ml, Adj[Ye], Ye]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Ye], Ye, Adj[Ye], me, Ye]*(-6 + 36*Zeta[3]) + 
  g1^2*g2^2*M1^2*MatMul[Adj[Ye], Ye]*(-162/5 + (324*Zeta[3])/5) + 
  g1^2*g2^2*M1*M2*MatMul[Adj[Ye], Ye]*(-162/5 + (324*Zeta[3])/5) + 
  g1^2*g2^2*M2^2*MatMul[Adj[Ye], Ye]*(-162/5 + (324*Zeta[3])/5) + 
  g2^2*M2^2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*(-12 + 72*Zeta[3]) + 
  g2^4*M2*MatMul[Adj[he], Ye]*(45 + 126*Zeta[3]) + 
  g2^4*M2*MatMul[Adj[Ye], he]*(45 + 126*Zeta[3]) + 
  g2^6*M2^2*(2157 + 3780*Zeta[3]) + g2^4*MatMul[Adj[Ye], Ye]*
   ((-45*mh1)/2 - 63*mh1*Zeta[3]) + g1^2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*
   (36*mh1 - (216*mh1*Zeta[3])/5) + g1^4*MatMul[Adj[Ye], Ye]*
   ((-2817*mh1)/50 - (36*mh2)/25 - (24*trace[mb])/25 - (72*trace[me])/25 - 
    (36*trace[ml])/25 - (12*trace[mq])/25 - (96*trace[mt])/25 + 
    (81*mh1*Zeta[3])/25) + g1^2*g2^2*MatMul[Adj[Ye], Ye]*
   ((-81*mh1)/5 + (162*mh1*Zeta[3])/5) + 
  g2^2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*(-12*mh1 + 72*mh1*Zeta[3]) + 
  g1^2*MatMul[Adj[he], Ye]*(16*trace[hb, Adj[Yb]] - 
    24*trace[hb, Adj[Yb]]*Zeta[3]) + g3^2*MatMul[Adj[he], Ye]*
   (-64*trace[hb, Adj[Yb]] + 96*trace[hb, Adj[Yb]]*Zeta[3]) + 
  g3^2*M3*MatMul[Adj[Ye], Ye]*(64*trace[hb, Adj[Yb]] + 
    64*trace[Adj[hb], Yb] - 96*trace[hb, Adj[Yb]]*Zeta[3] - 
    96*trace[Adj[hb], Yb]*Zeta[3]) + g1^2*MatMul[Adj[Ye], he]*
   (16*trace[Adj[hb], Yb] - 24*trace[Adj[hb], Yb]*Zeta[3]) + 
  g1^2*M1*MatMul[Adj[Ye], Ye]*(-16*trace[hb, Adj[Yb]] - 
    16*trace[Adj[hb], Yb] + 24*trace[hb, Adj[Yb]]*Zeta[3] + 
    24*trace[Adj[hb], Yb]*Zeta[3]) + g3^2*MatMul[Adj[Ye], he]*
   (-64*trace[Adj[hb], Yb] + 96*trace[Adj[hb], Yb]*Zeta[3]) + 
  g3^2*M3*MatMul[Adj[he], Ye]*(64*trace[Adj[Yb], Yb] - 
    96*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*M3*MatMul[Adj[Ye], he]*
   (64*trace[Adj[Yb], Yb] - 96*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g1^2*M1^2*MatMul[Adj[Ye], Ye]*(32*trace[Adj[Yb], Yb] - 
    48*trace[Adj[Yb], Yb]*Zeta[3]) + g1^2*MatMul[Adj[he], he]*
   (16*trace[Adj[Yb], Yb] - 24*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g1^2*MatMul[Adj[Ye], me, Ye]*(16*trace[Adj[Yb], Yb] - 
    24*trace[Adj[Yb], Yb]*Zeta[3]) + g1^2*MatMul[ml, Adj[Ye], Ye]*
   (8*trace[Adj[Yb], Yb] - 12*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g1^2*MatMul[Adj[Ye], Ye, ml]*(8*trace[Adj[Yb], Yb] - 
    12*trace[Adj[Yb], Yb]*Zeta[3]) + g1^2*M1*MatMul[Adj[he], Ye]*
   (-16*trace[Adj[Yb], Yb] + 24*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g1^2*M1*MatMul[Adj[Ye], he]*(-16*trace[Adj[Yb], Yb] + 
    24*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*MatMul[ml, Adj[Ye], Ye]*
   (-32*trace[Adj[Yb], Yb] + 48*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g3^2*MatMul[Adj[Ye], Ye, ml]*(-32*trace[Adj[Yb], Yb] + 
    48*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*MatMul[Adj[he], he]*
   (-64*trace[Adj[Yb], Yb] + 96*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g3^2*MatMul[Adj[Ye], me, Ye]*(-64*trace[Adj[Yb], Yb] + 
    96*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*M3^2*MatMul[Adj[Ye], Ye]*
   (-128*trace[Adj[Yb], Yb] + 192*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g1^2*MatMul[Adj[Ye], Ye]*(16*trace[hb, Adj[hb]] + 
    32*mh1*trace[Adj[Yb], Yb] + 16*trace[Yb, Adj[Yb], mb] + 
    16*trace[Adj[Yb], Yb, mq] - 24*trace[hb, Adj[hb]]*Zeta[3] - 
    48*mh1*trace[Adj[Yb], Yb]*Zeta[3] - 24*trace[Yb, Adj[Yb], mb]*Zeta[3] - 
    24*trace[Adj[Yb], Yb, mq]*Zeta[3]) + g3^2*MatMul[Adj[Ye], Ye]*
   (-64*trace[hb, Adj[hb]] - 128*mh1*trace[Adj[Yb], Yb] - 
    64*trace[Yb, Adj[Yb], mb] - 64*trace[Adj[Yb], Yb, mq] + 
    96*trace[hb, Adj[hb]]*Zeta[3] + 192*mh1*trace[Adj[Yb], Yb]*Zeta[3] + 
    96*trace[Yb, Adj[Yb], mb]*Zeta[3] + 96*trace[Adj[Yb], Yb, mq]*Zeta[3])}

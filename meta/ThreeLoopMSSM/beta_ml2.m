{(-6*g1^2*M1^2)/5 - 6*g2^2*M2^2 + 2*me*Ye*Adj[Ye] + 2*MatMul[hec, he] + 
  2*mh1*MatMul[Adj[Ye], Ye] + 2*ml*MatMul[Adj[Ye], Ye], 
 (621*g1^4*M1^2)/25 + (18*g1^2*g2^2*M1^2)/5 + (18*g1^2*g2^2*M1*M2)/5 + 
  (18*g1^2*g2^2*M2^2)/5 + 33*g2^4*M2^2 + (12*g1^2*me*Ye*Adj[Ye])/5 + 
  (12*g1^2*MatMul[hec, he])/5 - (12*g1^2*M1*MatMul[hec, Ye])/5 - 
  (12*g1^2*M1*MatMul[Adj[Ye], he])/5 + (24*g1^2*M1^2*MatMul[Adj[Ye], Ye])/5 + 
  (12*g1^2*mh1*MatMul[Adj[Ye], Ye])/5 + (12*g1^2*ml*MatMul[Adj[Ye], Ye])/5 - 
  4*ml*MatMul[Adj[Ye], Ye]^2 - 4*me*Adj[Ye]*MatMul[Ye, Adj[Ye], Ye] - 
  4*me*Ye*MatMul[Adj[Ye], Ye, Adj[Ye]] - 4*MatMul[hec, he, Adj[Ye], Ye] - 
  4*MatMul[hec, Ye, Adj[Ye], he] - 4*MatMul[Adj[Ye], he, hec, Ye] - 
  4*MatMul[Adj[Ye], Ye, hec, he] - 8*mh1*MatMul[Adj[Ye], Ye, Adj[Ye], Ye] - 
  4*ml*MatMul[Adj[Ye], Ye, Adj[Ye], Ye] + 
  g2^4*(3*mh1 + 3*mh2 + 3*trace[ml] + 9*trace[mq]) + 
  g1^4*((9*mh1)/25 + (9*mh2)/25 + (6*trace[mb])/25 + (18*trace[me])/25 + 
    (9*trace[ml])/25 + (3*trace[mq])/25 + (24*trace[mt])/25) + 
  MatMul[hec, Ye]*(-6*trace[hb, Ybc] - 2*trace[he, Adj[Ye]]) + 
  MatMul[Adj[Ye], he]*(-6*trace[hbc, Yb] - 2*trace[hec, Ye]) + 
  me*Ye*Adj[Ye]*(-6*trace[Ybc, Yb] - 2*trace[Adj[Ye], Ye]) + 
  MatMul[hec, he]*(-6*trace[Ybc, Yb] - 2*trace[Adj[Ye], Ye]) + 
  2*ml*MatMul[Adj[Ye], Ye]*(-3*trace[Ybc, Yb] - trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Ye], Ye]*(-6*trace[hb, hbc] - 2*trace[he, hec] - 
    12*mh1*trace[Ybc, Yb] - 4*mh1*trace[Adj[Ye], Ye] - 6*trace[Yb, Ybc, mb] - 
    6*trace[Ybc, Yb, mq] - 2*trace[Ye, Adj[Ye], me] - 
    2*trace[Adj[Ye], Ye, ml]), 
 g2^6*(-3*mh1 - 3*mh2 - 3*trace[ml] - 9*trace[mq]) + 
  g1^2*g2^4*((-9*mh1)/5 - (9*mh2)/5 - (9*trace[ml])/5 - (27*trace[mq])/5) + 
  g1^6*((-621*mh1)/125 - (621*mh2)/125 - (414*trace[mb])/125 - 
    (1242*trace[me])/125 - (621*trace[ml])/125 - (207*trace[mq])/125 - 
    (1656*trace[mt])/125) + g1^4*g2^2*((-27*mh1)/25 - (27*mh2)/25 - 
    (18*trace[mb])/25 - (54*trace[me])/25 - (27*trace[ml])/25 - 
    (9*trace[mq])/25 - (72*trace[mt])/25) + MatMul[hec, Ye, Adj[Ye], Ye]*
   (12*trace[hb, Ybc] + 4*trace[he, Adj[Ye]]) + 
  MatMul[Adj[Ye], Ye, hec, Ye]*(12*trace[hb, Ybc] + 4*trace[he, Adj[Ye]]) + 
  g2^2*MatMul[hec, Ye]*(36*trace[hb, Ybc] + 12*trace[he, Adj[Ye]]) + 
  g2^2*M2*MatMul[Adj[Ye], Ye]*(-36*trace[hb, Ybc] - 36*trace[hbc, Yb] - 
    12*trace[he, Adj[Ye]] - 12*trace[hec, Ye]) + 
  MatMul[Adj[Ye], he, Adj[Ye], Ye]*(12*trace[hbc, Yb] + 4*trace[hec, Ye]) + 
  MatMul[Adj[Ye], Ye, Adj[Ye], he]*(12*trace[hbc, Yb] + 4*trace[hec, Ye]) + 
  g2^2*MatMul[Adj[Ye], he]*(36*trace[hbc, Yb] + 12*trace[hec, Ye]) + 
  g1^4*M1*((42*trace[hb, Ybc])/5 + (42*trace[hbc, Yb])/5 + 
    (54*trace[he, Adj[Ye]])/5 + (54*trace[hec, Ye])/5 + 
    (78*trace[ht, Ytc])/5 + (78*trace[htc, Yt])/5) + 
  g2^4*M2*(90*trace[hb, Ybc] + 90*trace[hbc, Yb] + 30*trace[he, Adj[Ye]] + 
    30*trace[hec, Ye] + 90*trace[ht, Ytc] + 90*trace[htc, Yt]) + 
  g2^4*M2^2*(-270*trace[Ybc, Yb] - 270*trace[Ytc, Yt] - 
    90*trace[Adj[Ye], Ye]) + g1^4*M1^2*((-126*trace[Ybc, Yb])/5 - 
    (234*trace[Ytc, Yt])/5 - (162*trace[Adj[Ye], Ye])/5) + 
  g2^2*M2*MatMul[hec, Ye]*(-36*trace[Ybc, Yb] - 12*trace[Adj[Ye], Ye]) + 
  g2^2*M2*MatMul[Adj[Ye], he]*(-36*trace[Ybc, Yb] - 12*trace[Adj[Ye], Ye]) + 
  2*ml*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*(6*trace[Ybc, Yb] + 
    2*trace[Adj[Ye], Ye]) + ml*MatMul[Adj[Ye], Ye]^2*
   (12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  me*Adj[Ye]*MatMul[Ye, Adj[Ye], Ye]*(12*trace[Ybc, Yb] + 
    4*trace[Adj[Ye], Ye]) + me*Ye*MatMul[Adj[Ye], Ye, Adj[Ye]]*
   (12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[hec, he, Adj[Ye], Ye]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[hec, Ye, Adj[Ye], he]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Ye], he, hec, Ye]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Ye], Ye, hec, he]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  2*g2^2*ml*MatMul[Adj[Ye], Ye]*(18*trace[Ybc, Yb] + 6*trace[Adj[Ye], Ye]) + 
  g2^2*me*Ye*Adj[Ye]*(36*trace[Ybc, Yb] + 12*trace[Adj[Ye], Ye]) + 
  g2^2*MatMul[hec, he]*(36*trace[Ybc, Yb] + 12*trace[Adj[Ye], Ye]) + 
  g2^2*M2^2*MatMul[Adj[Ye], Ye]*(72*trace[Ybc, Yb] + 24*trace[Adj[Ye], Ye]) + 
  g2^4*(-54*trace[hb, hbc] - 18*trace[he, hec] - 54*trace[ht, htc] - 
    54*mh1*trace[Ybc, Yb] - 54*mh2*trace[Ytc, Yt] - 
    18*mh1*trace[Adj[Ye], Ye] - 54*trace[Yb, Ybc, mb] - 
    54*trace[Ybc, Yb, mq] - 18*trace[Ye, Adj[Ye], me] - 
    54*trace[Yt, Ytc, mt] - 54*trace[Ytc, Yt, mq] - 
    18*trace[Adj[Ye], Ye, ml]) + g1^4*((-126*trace[hb, hbc])/25 - 
    (162*trace[he, hec])/25 - (234*trace[ht, htc])/25 - 
    (126*mh1*trace[Ybc, Yb])/25 - (234*mh2*trace[Ytc, Yt])/25 - 
    (162*mh1*trace[Adj[Ye], Ye])/25 - (126*trace[Yb, Ybc, mb])/25 - 
    (126*trace[Ybc, Yb, mq])/25 - (162*trace[Ye, Adj[Ye], me])/25 - 
    (234*trace[Yt, Ytc, mt])/25 - (234*trace[Ytc, Yt, mq])/25 - 
    (162*trace[Adj[Ye], Ye, ml])/25) + MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*
   (12*trace[hb, hbc] + 4*trace[he, hec] + 36*mh1*trace[Ybc, Yb] + 
    12*mh1*trace[Adj[Ye], Ye] + 12*trace[Yb, Ybc, mb] + 
    12*trace[Ybc, Yb, mq] + 4*trace[Ye, Adj[Ye], me] + 
    4*trace[Adj[Ye], Ye, ml]) + g2^2*MatMul[Adj[Ye], Ye]*
   (36*trace[hb, hbc] + 12*trace[he, hec] + 72*mh1*trace[Ybc, Yb] + 
    24*mh1*trace[Adj[Ye], Ye] + 36*trace[Yb, Ybc, mb] + 
    36*trace[Ybc, Yb, mq] + 12*trace[Ye, Adj[Ye], me] + 
    12*trace[Adj[Ye], Ye, ml]) + MatMul[hec, Ye]*
   (-36*trace[hb, Ybc]*trace[Ybc, Yb] - 12*trace[he, Adj[Ye]]*
     trace[Ybc, Yb] - 12*trace[hb, Ybc]*trace[Adj[Ye], Ye] - 
    4*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] + 12*trace[hb, Ytc, Yt, Ybc] + 
    72*trace[Ybc, hb, Ybc, Yb] + 12*trace[Ytc, ht, Ybc, Yb] + 
    24*trace[Adj[Ye], he, Adj[Ye], Ye]) + MatMul[Adj[Ye], he]*
   (-36*trace[hbc, Yb]*trace[Ybc, Yb] - 12*trace[hec, Ye]*trace[Ybc, Yb] - 
    12*trace[hbc, Yb]*trace[Adj[Ye], Ye] - 4*trace[hec, Ye]*
     trace[Adj[Ye], Ye] + 12*trace[htc, Yt, Ybc, Yb] + 
    72*trace[Ybc, Yb, hbc, Yb] + 12*trace[Ytc, Yt, hbc, Yb] + 
    24*trace[Adj[Ye], Ye, hec, Ye]) + 2*ml*MatMul[Adj[Ye], Ye]*
   (-9*trace[Ybc, Yb]^2 - 6*trace[Ybc, Yb]*trace[Adj[Ye], Ye] - 
    trace[Adj[Ye], Ye]^2 + 18*trace[Ybc, Yb, Ybc, Yb] + 
    6*trace[Ytc, Yt, Ybc, Yb] + 6*trace[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  me*Ye*Adj[Ye]*(-18*trace[Ybc, Yb]^2 - 12*trace[Ybc, Yb]*
     trace[Adj[Ye], Ye] - 2*trace[Adj[Ye], Ye]^2 + 
    36*trace[Ybc, Yb, Ybc, Yb] + 12*trace[Ytc, Yt, Ybc, Yb] + 
    12*trace[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  MatMul[hec, he]*(-18*trace[Ybc, Yb]^2 - 12*trace[Ybc, Yb]*
     trace[Adj[Ye], Ye] - 2*trace[Adj[Ye], Ye]^2 + 
    36*trace[Ybc, Yb, Ybc, Yb] + 12*trace[Ytc, Yt, Ybc, Yb] + 
    12*trace[Adj[Ye], Ye, Adj[Ye], Ye]) + MatMul[Adj[Ye], Ye]*
   (-36*trace[hb, Ybc]*trace[hbc, Yb] - 12*trace[hbc, Yb]*
     trace[he, Adj[Ye]] - 12*trace[hb, Ybc]*trace[hec, Ye] - 
    4*trace[he, Adj[Ye]]*trace[hec, Ye] - 36*trace[hb, hbc]*trace[Ybc, Yb] - 
    12*trace[he, hec]*trace[Ybc, Yb] - 54*mh1*trace[Ybc, Yb]^2 - 
    12*trace[hb, hbc]*trace[Adj[Ye], Ye] - 4*trace[he, hec]*
     trace[Adj[Ye], Ye] - 36*mh1*trace[Ybc, Yb]*trace[Adj[Ye], Ye] - 
    6*mh1*trace[Adj[Ye], Ye]^2 - 36*trace[Ybc, Yb]*trace[Yb, Ybc, mb] - 
    12*trace[Adj[Ye], Ye]*trace[Yb, Ybc, mb] - 36*trace[Ybc, Yb]*
     trace[Ybc, Yb, mq] - 12*trace[Adj[Ye], Ye]*trace[Ybc, Yb, mq] - 
    12*trace[Ybc, Yb]*trace[Ye, Adj[Ye], me] - 4*trace[Adj[Ye], Ye]*
     trace[Ye, Adj[Ye], me] - 12*trace[Ybc, Yb]*trace[Adj[Ye], Ye, ml] - 
    4*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, ml] + 
    12*trace[hb, htc, Yt, Ybc] + 72*trace[hbc, hb, Ybc, Yb] + 
    12*trace[hbc, hb, Ytc, Yt] + 24*trace[hec, he, Adj[Ye], Ye] + 
    12*trace[htc, ht, Ybc, Yb] + 72*trace[Ybc, hb, hbc, Yb] + 
    108*mh1*trace[Ybc, Yb, Ybc, Yb] + 12*trace[Ytc, ht, hbc, Yb] + 
    24*mh1*trace[Ytc, Yt, Ybc, Yb] + 12*mh2*trace[Ytc, Yt, Ybc, Yb] + 
    24*trace[Adj[Ye], he, hec, Ye] + 36*mh1*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    72*trace[Yb, Ybc, Yb, Ybc, mb] + 12*trace[Yb, Ytc, Yt, Ybc, mb] + 
    72*trace[Ybc, Yb, Ybc, Yb, mq] + 12*trace[Ybc, Yb, Ytc, Yt, mq] + 
    24*trace[Ye, Adj[Ye], Ye, Adj[Ye], me] + 12*trace[Yt, Ybc, Yb, Ytc, mt] + 
    12*trace[Ytc, Yt, Ybc, Yb, mq] + 24*trace[Adj[Ye], Ye, Adj[Ye], Ye, 
      ml]) + g2^4*g3^2*M2^2*(1080 - 1296*Zeta[3]) + 
  g2^4*g3^2*M2*M3*(720 - 864*Zeta[3]) + g2^4*g3^2*M3^2*(432 - 432*Zeta[3]) + 
  g2^4*M2^2*MatMul[Adj[Ye], Ye]*(-135 - 378*Zeta[3]) + 
  g1^4*g3^2*M1^2*(792/5 - (4752*Zeta[3])/25) + 
  g1^6*M1^2*(55767/125 - (21492*Zeta[3])/125) + 
  g1^4*g3^2*M1*M3*(528/5 - (3168*Zeta[3])/25) + 
  g1^2*g2^4*M2^2*(171/5 - (486*Zeta[3])/5) + 
  g1^2*g2^4*M1*M2*(18 - (324*Zeta[3])/5) + 
  g1^4*g3^2*M3^2*(1584/25 - (1584*Zeta[3])/25) + 
  g2^4*me*Ye*Adj[Ye]*(-45/2 - 63*Zeta[3]) + g2^4*MatMul[hec, he]*
   (-45/2 - 63*Zeta[3]) + g1^4*g2^2*M1^2*(81/25 - (1458*Zeta[3])/25) + 
  g1^2*M1^2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*(36 - (216*Zeta[3])/5) + 
  g1^4*g2^2*M1*M2*(54/25 - (972*Zeta[3])/25) + 
  g2^2*M2*MatMul[hec, Ye, Adj[Ye], Ye]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Ye], he, Adj[Ye], Ye]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Ye], Ye, hec, Ye]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Ye], Ye, Adj[Ye], he]*(6 - 36*Zeta[3]) + 
  g1^2*g2^4*M1^2*(72/5 - (162*Zeta[3])/5) + g1^2*g2^2*M1*MatMul[hec, Ye]*
   (81/5 - (162*Zeta[3])/5) + g1^2*g2^2*M2*MatMul[hec, Ye]*
   (81/5 - (162*Zeta[3])/5) + g1^2*g2^2*M1*MatMul[Adj[Ye], he]*
   (81/5 - (162*Zeta[3])/5) + g1^2*g2^2*M2*MatMul[Adj[Ye], he]*
   (81/5 - (162*Zeta[3])/5) + 2*g2^4*ml*MatMul[Adj[Ye], Ye]*
   (-45/4 - (63*Zeta[3])/2) + g1^2*ml*MatMul[Adj[Ye], Ye]^2*
   (18 - (108*Zeta[3])/5) + g1^2*me*Adj[Ye]*MatMul[Ye, Adj[Ye], Ye]*
   (18 - (108*Zeta[3])/5) + g1^2*me*Ye*MatMul[Adj[Ye], Ye, Adj[Ye]]*
   (18 - (108*Zeta[3])/5) + g1^2*MatMul[hec, he, Adj[Ye], Ye]*
   (18 - (108*Zeta[3])/5) + g1^2*MatMul[hec, Ye, Adj[Ye], he]*
   (18 - (108*Zeta[3])/5) + g1^2*MatMul[Adj[Ye], he, hec, Ye]*
   (18 - (108*Zeta[3])/5) + g1^2*MatMul[Adj[Ye], Ye, hec, he]*
   (18 - (108*Zeta[3])/5) + g1^4*g2^2*M2^2*(108/25 - (486*Zeta[3])/25) + 
  2*g1^2*ml*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*(9 - (54*Zeta[3])/5) + 
  g1^4*M1*MatMul[hec, Ye]*(549/5 - (162*Zeta[3])/25) + 
  g1^4*M1*MatMul[Adj[Ye], he]*(549/5 - (162*Zeta[3])/25) + 
  12*me*MatMul[Ye, Adj[Ye], Ye]*MatMul[Adj[Ye], Ye, Adj[Ye]]*Zeta[3] + 
  24*ml*MatMul[Adj[Ye], Ye]*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
  12*me*Adj[Ye]*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
  12*me*Ye*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye]]*Zeta[3] + 
  12*MatMul[hec, he, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
  12*MatMul[hec, Ye, Adj[Ye], he, Adj[Ye], Ye]*Zeta[3] + 
  12*MatMul[hec, Ye, Adj[Ye], Ye, Adj[Ye], he]*Zeta[3] + 
  12*MatMul[Adj[Ye], he, hec, Ye, Adj[Ye], Ye]*Zeta[3] + 
  12*MatMul[Adj[Ye], he, Adj[Ye], Ye, hec, Ye]*Zeta[3] + 
  12*MatMul[Adj[Ye], Ye, hec, he, Adj[Ye], Ye]*Zeta[3] + 
  12*MatMul[Adj[Ye], Ye, hec, Ye, Adj[Ye], he]*Zeta[3] + 
  12*MatMul[Adj[Ye], Ye, Adj[Ye], he, hec, Ye]*Zeta[3] + 
  12*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, hec, he]*Zeta[3] + 
  36*mh1*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
  12*ml*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
  2*g1^4*ml*MatMul[Adj[Ye], Ye]*(-549/20 + (81*Zeta[3])/50) + 
  g1^4*me*Ye*Adj[Ye]*(-549/10 + (81*Zeta[3])/25) + 
  g1^4*MatMul[hec, he]*(-549/10 + (81*Zeta[3])/25) + 
  2*g1^2*g2^2*ml*MatMul[Adj[Ye], Ye]*(-81/10 + (81*Zeta[3])/5) + 
  2*g2^2*ml*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*(-3 + 18*Zeta[3]) + 
  g1^4*M1^2*MatMul[Adj[Ye], Ye]*(-1647/5 + (486*Zeta[3])/25) + 
  g1^2*M1*MatMul[hec, Ye, Adj[Ye], Ye]*(-18 + (108*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Ye], he, Adj[Ye], Ye]*(-18 + (108*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Ye], Ye, hec, Ye]*(-18 + (108*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Ye], Ye, Adj[Ye], he]*(-18 + (108*Zeta[3])/5) + 
  g1^2*g2^2*me*Ye*Adj[Ye]*(-81/5 + (162*Zeta[3])/5) + 
  g1^2*g2^2*MatMul[hec, he]*(-81/5 + (162*Zeta[3])/5) + 
  g2^2*ml*MatMul[Adj[Ye], Ye]^2*(-6 + 36*Zeta[3]) + 
  g2^2*me*Adj[Ye]*MatMul[Ye, Adj[Ye], Ye]*(-6 + 36*Zeta[3]) + 
  g2^2*me*Ye*MatMul[Adj[Ye], Ye, Adj[Ye]]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[hec, he, Adj[Ye], Ye]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[hec, Ye, Adj[Ye], he]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Ye], he, hec, Ye]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Ye], Ye, hec, he]*(-6 + 36*Zeta[3]) + 
  g1^2*g2^2*M1^2*MatMul[Adj[Ye], Ye]*(-162/5 + (324*Zeta[3])/5) + 
  g1^2*g2^2*M1*M2*MatMul[Adj[Ye], Ye]*(-162/5 + (324*Zeta[3])/5) + 
  g1^2*g2^2*M2^2*MatMul[Adj[Ye], Ye]*(-162/5 + (324*Zeta[3])/5) + 
  g2^2*M2^2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*(-12 + 72*Zeta[3]) + 
  g2^4*M2*MatMul[hec, Ye]*(45 + 126*Zeta[3]) + g2^4*M2*MatMul[Adj[Ye], he]*
   (45 + 126*Zeta[3]) + g2^6*M2^2*(2157 + 3780*Zeta[3]) + 
  g2^4*MatMul[Adj[Ye], Ye]*((-45*mh1)/2 - 63*mh1*Zeta[3]) + 
  g1^2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*(36*mh1 - (216*mh1*Zeta[3])/5) + 
  g1^4*MatMul[Adj[Ye], Ye]*((-2817*mh1)/50 - (36*mh2)/25 - 
    (24*trace[mb])/25 - (72*trace[me])/25 - (36*trace[ml])/25 - 
    (12*trace[mq])/25 - (96*trace[mt])/25 + (81*mh1*Zeta[3])/25) + 
  g1^2*g2^2*MatMul[Adj[Ye], Ye]*((-81*mh1)/5 + (162*mh1*Zeta[3])/5) + 
  g2^2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*(-12*mh1 + 72*mh1*Zeta[3]) + 
  g1^2*MatMul[hec, Ye]*(16*trace[hb, Ybc] - 24*trace[hb, Ybc]*Zeta[3]) + 
  g3^2*MatMul[hec, Ye]*(-64*trace[hb, Ybc] + 96*trace[hb, Ybc]*Zeta[3]) + 
  g3^2*M3*MatMul[Adj[Ye], Ye]*(64*trace[hb, Ybc] + 64*trace[hbc, Yb] - 
    96*trace[hb, Ybc]*Zeta[3] - 96*trace[hbc, Yb]*Zeta[3]) + 
  g1^2*MatMul[Adj[Ye], he]*(16*trace[hbc, Yb] - 24*trace[hbc, Yb]*Zeta[3]) + 
  g1^2*M1*MatMul[Adj[Ye], Ye]*(-16*trace[hb, Ybc] - 16*trace[hbc, Yb] + 
    24*trace[hb, Ybc]*Zeta[3] + 24*trace[hbc, Yb]*Zeta[3]) + 
  g3^2*MatMul[Adj[Ye], he]*(-64*trace[hbc, Yb] + 96*trace[hbc, Yb]*Zeta[3]) + 
  g3^2*M3*MatMul[hec, Ye]*(64*trace[Ybc, Yb] - 96*trace[Ybc, Yb]*Zeta[3]) + 
  g3^2*M3*MatMul[Adj[Ye], he]*(64*trace[Ybc, Yb] - 
    96*trace[Ybc, Yb]*Zeta[3]) + g1^2*M1^2*MatMul[Adj[Ye], Ye]*
   (32*trace[Ybc, Yb] - 48*trace[Ybc, Yb]*Zeta[3]) + 
  g1^2*me*Ye*Adj[Ye]*(16*trace[Ybc, Yb] - 24*trace[Ybc, Yb]*Zeta[3]) + 
  g1^2*MatMul[hec, he]*(16*trace[Ybc, Yb] - 24*trace[Ybc, Yb]*Zeta[3]) + 
  2*g1^2*ml*MatMul[Adj[Ye], Ye]*(8*trace[Ybc, Yb] - 
    12*trace[Ybc, Yb]*Zeta[3]) + g1^2*M1*MatMul[hec, Ye]*
   (-16*trace[Ybc, Yb] + 24*trace[Ybc, Yb]*Zeta[3]) + 
  g1^2*M1*MatMul[Adj[Ye], he]*(-16*trace[Ybc, Yb] + 
    24*trace[Ybc, Yb]*Zeta[3]) + 2*g3^2*ml*MatMul[Adj[Ye], Ye]*
   (-32*trace[Ybc, Yb] + 48*trace[Ybc, Yb]*Zeta[3]) + 
  g3^2*me*Ye*Adj[Ye]*(-64*trace[Ybc, Yb] + 96*trace[Ybc, Yb]*Zeta[3]) + 
  g3^2*MatMul[hec, he]*(-64*trace[Ybc, Yb] + 96*trace[Ybc, Yb]*Zeta[3]) + 
  g3^2*M3^2*MatMul[Adj[Ye], Ye]*(-128*trace[Ybc, Yb] + 
    192*trace[Ybc, Yb]*Zeta[3]) + g1^2*MatMul[Adj[Ye], Ye]*
   (16*trace[hb, hbc] + 32*mh1*trace[Ybc, Yb] + 16*trace[Yb, Ybc, mb] + 
    16*trace[Ybc, Yb, mq] - 24*trace[hb, hbc]*Zeta[3] - 
    48*mh1*trace[Ybc, Yb]*Zeta[3] - 24*trace[Yb, Ybc, mb]*Zeta[3] - 
    24*trace[Ybc, Yb, mq]*Zeta[3]) + g3^2*MatMul[Adj[Ye], Ye]*
   (-64*trace[hb, hbc] - 128*mh1*trace[Ybc, Yb] - 64*trace[Yb, Ybc, mb] - 
    64*trace[Ybc, Yb, mq] + 96*trace[hb, hbc]*Zeta[3] + 
    192*mh1*trace[Ybc, Yb]*Zeta[3] + 96*trace[Yb, Ybc, mb]*Zeta[3] + 
    96*trace[Ybc, Yb, mq]*Zeta[3])}

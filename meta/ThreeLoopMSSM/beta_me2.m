{4*he*hec - (24*g1^2*M1^2)/5 + 4*ml*Ye*Adj[Ye] + 4*me*MatMul[Ye, Adj[Ye]] + 
  4*mh1*MatMul[Ye, Adj[Ye]], (-12*g1^2*he*hec)/5 + 12*g2^2*he*hec + 
  (2808*g1^4*M1^2)/25 + (12*g1^2*hec*M1*Ye)/5 - 12*g2^2*hec*M2*Ye - 
  (12*g1^2*ml*Ye*Adj[Ye])/5 + 12*g2^2*ml*Ye*Adj[Ye] + 
  (12*g1^2*M1*MatMul[he, Adj[Ye]])/5 - 12*g2^2*M2*MatMul[he, Adj[Ye]] - 
  (24*g1^2*M1^2*MatMul[Ye, Adj[Ye]])/5 + 24*g2^2*M2^2*MatMul[Ye, Adj[Ye]] - 
  (12*g1^2*me*MatMul[Ye, Adj[Ye]])/5 + 12*g2^2*me*MatMul[Ye, Adj[Ye]] - 
  (12*g1^2*mh1*MatMul[Ye, Adj[Ye]])/5 + 12*g2^2*mh1*MatMul[Ye, Adj[Ye]] - 
  4*me*MatMul[Ye, Adj[Ye]]^2 - 4*hec*MatMul[he, Adj[Ye], Ye] - 
  4*hec*MatMul[Ye, Adj[Ye], he] - 4*ml*Adj[Ye]*MatMul[Ye, Adj[Ye], Ye] - 
  4*ml*Ye*MatMul[Adj[Ye], Ye, Adj[Ye]] - 4*MatMul[he, hec, Ye, Adj[Ye]] - 
  4*MatMul[Ye, hec, he, Adj[Ye]] - 4*me*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]] - 
  8*mh1*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]] + 
  g1^4*((36*mh1)/25 + (36*mh2)/25 + (24*trace[mb])/25 + (72*trace[me])/25 + 
    (36*trace[ml])/25 + (12*trace[mq])/25 + (96*trace[mt])/25) + 
  Ye*(-12*hec*trace[hb, Adj[Yb]] - 4*hec*trace[he, Adj[Ye]]) + 
  MatMul[he, Adj[Ye]]*(-12*trace[hbc, Yb] - 4*trace[hec, Ye]) + 
  ml*Ye*Adj[Ye]*(-12*trace[Adj[Yb], Yb] - 4*trace[Adj[Ye], Ye]) + 
  2*me*MatMul[Ye, Adj[Ye]]*(-6*trace[Adj[Yb], Yb] - 2*trace[Adj[Ye], Ye]) + 
  he*(-12*hec*trace[Adj[Yb], Yb] - 4*hec*trace[Adj[Ye], Ye]) + 
  MatMul[Ye, Adj[Ye]]*(-12*trace[hb, hbc] - 4*trace[he, hec] - 
    24*mh1*trace[Adj[Yb], Yb] - 8*mh1*trace[Adj[Ye], Ye] - 
    12*trace[Yb, Adj[Yb], mb] - 4*trace[Ye, Adj[Ye], me] - 
    12*trace[Adj[Yb], Yb, mq] - 4*trace[Adj[Ye], Ye, ml]), 
 g1^6*((-2808*mh1)/125 - (2808*mh2)/125 - (1872*trace[mb])/125 - 
    (5616*trace[me])/125 - (2808*trace[ml])/125 - (936*trace[mq])/125 - 
    (7488*trace[mt])/125) + MatMul[Ye, hec, Ye, Adj[Ye]]*
   (12*trace[hb, Adj[Yb]] + 4*trace[he, Adj[Ye]]) + 
  MatMul[Ye, Adj[Ye], Ye]*(12*hec*trace[hb, Adj[Yb]] + 
    4*hec*trace[he, Adj[Ye]]) + MatMul[he, Adj[Ye], Ye, Adj[Ye]]*
   (12*trace[hbc, Yb] + 4*trace[hec, Ye]) + MatMul[Ye, Adj[Ye], he, Adj[Ye]]*
   (12*trace[hbc, Yb] + 4*trace[hec, Ye]) + 
  g1^4*M1*((168*trace[hb, Adj[Yb]])/5 + (168*trace[hbc, Yb])/5 + 
    (216*trace[he, Adj[Ye]])/5 + (216*trace[hec, Ye])/5 + 
    (312*trace[ht, Adj[Yt]])/5 + (312*trace[htc, Yt])/5) + 
  2*me*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*(6*trace[Adj[Yb], Yb] + 
    2*trace[Adj[Ye], Ye]) + me*MatMul[Ye, Adj[Ye]]^2*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  ml*Adj[Ye]*MatMul[Ye, Adj[Ye], Ye]*(12*trace[Adj[Yb], Yb] + 
    4*trace[Adj[Ye], Ye]) + ml*Ye*MatMul[Adj[Ye], Ye, Adj[Ye]]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[he, hec, Ye, Adj[Ye]]*(12*trace[Adj[Yb], Yb] + 
    4*trace[Adj[Ye], Ye]) + MatMul[Ye, hec, he, Adj[Ye]]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[he, Adj[Ye], Ye]*(12*hec*trace[Adj[Yb], Yb] + 
    4*hec*trace[Adj[Ye], Ye]) + MatMul[Ye, Adj[Ye], he]*
   (12*hec*trace[Adj[Yb], Yb] + 4*hec*trace[Adj[Ye], Ye]) + 
  g1^4*M1^2*((-504*trace[Adj[Yb], Yb])/5 - (648*trace[Adj[Ye], Ye])/5 - 
    (936*trace[Adj[Yt], Yt])/5) + MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*
   (12*trace[hb, hbc] + 4*trace[he, hec] + 36*mh1*trace[Adj[Yb], Yb] + 
    12*mh1*trace[Adj[Ye], Ye] + 12*trace[Yb, Adj[Yb], mb] + 
    4*trace[Ye, Adj[Ye], me] + 12*trace[Adj[Yb], Yb, mq] + 
    4*trace[Adj[Ye], Ye, ml]) + g1^4*((-504*trace[hb, hbc])/25 - 
    (648*trace[he, hec])/25 - (936*trace[ht, htc])/25 - 
    (504*mh1*trace[Adj[Yb], Yb])/25 - (648*mh1*trace[Adj[Ye], Ye])/25 - 
    (936*mh2*trace[Adj[Yt], Yt])/25 - (504*trace[Yb, Adj[Yb], mb])/25 - 
    (648*trace[Ye, Adj[Ye], me])/25 - (936*trace[Yt, Adj[Yt], mt])/25 - 
    (504*trace[Adj[Yb], Yb, mq])/25 - (648*trace[Adj[Ye], Ye, ml])/25 - 
    (936*trace[Adj[Yt], Yt, mq])/25) + 
  Ye*(-72*hec*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
    24*hec*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb] - 
    24*hec*trace[hb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
    8*hec*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] + 
    24*hec*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
    144*hec*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
    48*hec*trace[Adj[Ye], he, Adj[Ye], Ye] + 
    24*hec*trace[Adj[Yt], ht, Adj[Yb], Yb]) + 
  MatMul[he, Adj[Ye]]*(-72*trace[hbc, Yb]*trace[Adj[Yb], Yb] - 
    24*trace[hec, Ye]*trace[Adj[Yb], Yb] - 24*trace[hbc, Yb]*
     trace[Adj[Ye], Ye] - 8*trace[hec, Ye]*trace[Adj[Ye], Ye] + 
    24*trace[htc, Yt, Adj[Yb], Yb] + 144*trace[Adj[Yb], Yb, hbc, Yb] + 
    48*trace[Adj[Ye], Ye, hec, Ye] + 24*trace[Adj[Yt], Yt, hbc, Yb]) + 
  2*me*MatMul[Ye, Adj[Ye]]*(-18*trace[Adj[Yb], Yb]^2 - 
    12*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 2*trace[Adj[Ye], Ye]^2 + 
    36*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 12*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    12*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + 
  ml*Ye*Adj[Ye]*(-36*trace[Adj[Yb], Yb]^2 - 24*trace[Adj[Yb], Yb]*
     trace[Adj[Ye], Ye] - 4*trace[Adj[Ye], Ye]^2 + 
    72*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 24*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    24*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + 
  he*(-36*hec*trace[Adj[Yb], Yb]^2 - 24*hec*trace[Adj[Yb], Yb]*
     trace[Adj[Ye], Ye] - 4*hec*trace[Adj[Ye], Ye]^2 + 
    72*hec*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    24*hec*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    24*hec*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + 
  MatMul[Ye, Adj[Ye]]*(-72*trace[hb, Adj[Yb]]*trace[hbc, Yb] - 
    24*trace[hbc, Yb]*trace[he, Adj[Ye]] - 24*trace[hb, Adj[Yb]]*
     trace[hec, Ye] - 8*trace[he, Adj[Ye]]*trace[hec, Ye] - 
    72*trace[hb, hbc]*trace[Adj[Yb], Yb] - 24*trace[he, hec]*
     trace[Adj[Yb], Yb] - 108*mh1*trace[Adj[Yb], Yb]^2 - 
    24*trace[hb, hbc]*trace[Adj[Ye], Ye] - 8*trace[he, hec]*
     trace[Adj[Ye], Ye] - 72*mh1*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
    12*mh1*trace[Adj[Ye], Ye]^2 - 72*trace[Adj[Yb], Yb]*
     trace[Yb, Adj[Yb], mb] - 24*trace[Adj[Ye], Ye]*trace[Yb, Adj[Yb], mb] - 
    24*trace[Adj[Yb], Yb]*trace[Ye, Adj[Ye], me] - 
    8*trace[Adj[Ye], Ye]*trace[Ye, Adj[Ye], me] - 
    72*trace[Adj[Yb], Yb]*trace[Adj[Yb], Yb, mq] - 
    24*trace[Adj[Ye], Ye]*trace[Adj[Yb], Yb, mq] - 
    24*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye, ml] - 
    8*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, ml] + 
    24*trace[hb, htc, Yt, Adj[Yb]] + 144*trace[hbc, hb, Adj[Yb], Yb] + 
    24*trace[hbc, hb, Adj[Yt], Yt] + 48*trace[hec, he, Adj[Ye], Ye] + 
    24*trace[htc, ht, Adj[Yb], Yb] + 144*trace[Adj[Yb], hb, hbc, Yb] + 
    216*mh1*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    48*trace[Adj[Ye], he, hec, Ye] + 72*mh1*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    24*trace[Adj[Yt], ht, hbc, Yb] + 48*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    24*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    144*trace[Yb, Adj[Yb], Yb, Adj[Yb], mb] + 
    24*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
    48*trace[Ye, Adj[Ye], Ye, Adj[Ye], me] + 
    24*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
    144*trace[Adj[Yb], Yb, Adj[Yb], Yb, mq] + 
    24*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
    48*trace[Adj[Ye], Ye, Adj[Ye], Ye, ml] + 
    24*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq]) + 
  g1^4*g3^2*M1^2*(3168/5 - (19008*Zeta[3])/25) + 
  g1^6*M1^2*(191964/125 - (85968*Zeta[3])/125) + 
  g1^4*g3^2*M1*M3*(2112/5 - (12672*Zeta[3])/25) + 
  g1^4*g3^2*M3^2*(6336/25 - (6336*Zeta[3])/25) + 
  g1^4*g2^2*M1^2*(972/5 - (5832*Zeta[3])/25) + g1^4*M1^2*MatMul[Ye, Adj[Ye]]*
   (-9018/25 - (972*Zeta[3])/5) + g1^4*g2^2*M1*M2*
   (648/5 - (3888*Zeta[3])/25) + g2^4*M2^2*MatMul[Ye, Adj[Ye]]*
   (-474 - 108*Zeta[3]) + g1^4*g2^2*M2^2*(1944/25 - (1944*Zeta[3])/25) + 
  g2^2*M2^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*(36 - 72*Zeta[3]) + 
  g1^2*g2^2*M1*MatMul[he, Adj[Ye]]*(54 - (324*Zeta[3])/5) + 
  g1^2*g2^2*M2*MatMul[he, Adj[Ye]]*(54 - (324*Zeta[3])/5) + 
  g2^2*me*MatMul[Ye, Adj[Ye]]^2*(18 - 36*Zeta[3]) + 
  g2^2*ml*Adj[Ye]*MatMul[Ye, Adj[Ye], Ye]*(18 - 36*Zeta[3]) + 
  g2^2*ml*Ye*MatMul[Adj[Ye], Ye, Adj[Ye]]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[he, hec, Ye, Adj[Ye]]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[Ye, hec, he, Adj[Ye]]*(18 - 36*Zeta[3]) + 
  g1^4*ml*Ye*Adj[Ye]*(-1503/25 - (162*Zeta[3])/5) + 
  g1^2*M1*MatMul[he, Adj[Ye], Ye, Adj[Ye]]*(-18/5 - (108*Zeta[3])/5) + 
  g1^2*M1*MatMul[Ye, hec, Ye, Adj[Ye]]*(-18/5 - (108*Zeta[3])/5) + 
  g1^2*M1*MatMul[Ye, Adj[Ye], he, Adj[Ye]]*(-18/5 - (108*Zeta[3])/5) + 
  g2^4*ml*Ye*Adj[Ye]*(-87 - 18*Zeta[3]) + 
  2*g2^2*me*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*(9 - 18*Zeta[3]) + 
  2*g1^4*me*MatMul[Ye, Adj[Ye]]*(-1503/50 - (81*Zeta[3])/5) + 
  2*g2^4*me*MatMul[Ye, Adj[Ye]]*(-87/2 - 9*Zeta[3]) + 
  2*g1^2*me*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*(9/5 + (54*Zeta[3])/5) + 
  2*me*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye]]*(6 + 12*Zeta[3]) + 
  g1^2*me*MatMul[Ye, Adj[Ye]]^2*(18/5 + (108*Zeta[3])/5) + 
  g1^2*ml*Adj[Ye]*MatMul[Ye, Adj[Ye], Ye]*(18/5 + (108*Zeta[3])/5) + 
  g1^2*ml*Ye*MatMul[Adj[Ye], Ye, Adj[Ye]]*(18/5 + (108*Zeta[3])/5) + 
  g1^2*MatMul[he, hec, Ye, Adj[Ye]]*(18/5 + (108*Zeta[3])/5) + 
  g1^2*MatMul[Ye, hec, he, Adj[Ye]]*(18/5 + (108*Zeta[3])/5) + 
  ml*MatMul[Ye, Adj[Ye], Ye]*MatMul[Adj[Ye], Ye, Adj[Ye]]*(12 + 24*Zeta[3]) + 
  2*me*MatMul[Ye, Adj[Ye]]*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*
   (12 + 24*Zeta[3]) + ml*Adj[Ye]*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye]*
   (12 + 24*Zeta[3]) + ml*Ye*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye]]*
   (12 + 24*Zeta[3]) + MatMul[he, hec, Ye, Adj[Ye], Ye, Adj[Ye]]*
   (12 + 24*Zeta[3]) + MatMul[he, Adj[Ye], Ye, hec, Ye, Adj[Ye]]*
   (12 + 24*Zeta[3]) + MatMul[Ye, hec, he, Adj[Ye], Ye, Adj[Ye]]*
   (12 + 24*Zeta[3]) + MatMul[Ye, hec, Ye, Adj[Ye], he, Adj[Ye]]*
   (12 + 24*Zeta[3]) + MatMul[Ye, Adj[Ye], he, hec, Ye, Adj[Ye]]*
   (12 + 24*Zeta[3]) + MatMul[Ye, Adj[Ye], Ye, hec, he, Adj[Ye]]*
   (12 + 24*Zeta[3]) + 2*g1^2*g2^2*me*MatMul[Ye, Adj[Ye]]*
   (-27 + (162*Zeta[3])/5) + g2^2*M2*MatMul[he, Adj[Ye], Ye, Adj[Ye]]*
   (-18 + 36*Zeta[3]) + g2^2*M2*MatMul[Ye, hec, Ye, Adj[Ye]]*
   (-18 + 36*Zeta[3]) + g2^2*M2*MatMul[Ye, Adj[Ye], he, Adj[Ye]]*
   (-18 + 36*Zeta[3]) + g2^4*M2*MatMul[he, Adj[Ye]]*(174 + 36*Zeta[3]) + 
  g1^2*M1^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*(36/5 + (216*Zeta[3])/5) + 
  g1^2*g2^2*ml*Ye*Adj[Ye]*(-54 + (324*Zeta[3])/5) + 
  g1^4*M1*MatMul[he, Adj[Ye]]*(3006/25 + (324*Zeta[3])/5) + 
  g1^2*g2^2*M1^2*MatMul[Ye, Adj[Ye]]*(-108 + (648*Zeta[3])/5) + 
  g1^2*g2^2*M1*M2*MatMul[Ye, Adj[Ye]]*(-108 + (648*Zeta[3])/5) + 
  g1^2*g2^2*M2^2*MatMul[Ye, Adj[Ye]]*(-108 + (648*Zeta[3])/5) + 
  g1^2*g2^2*M1*Ye*(54*hec - (324*hec*Zeta[3])/5) + 
  g1^2*g2^2*M2*Ye*(54*hec - (324*hec*Zeta[3])/5) + 
  g2^2*MatMul[he, Adj[Ye], Ye]*(18*hec - 36*hec*Zeta[3]) + 
  g2^2*MatMul[Ye, Adj[Ye], he]*(18*hec - 36*hec*Zeta[3]) + 
  g1^4*he*((-1503*hec)/25 - (162*hec*Zeta[3])/5) + 
  g1^2*M1*MatMul[Ye, Adj[Ye], Ye]*((-18*hec)/5 - (108*hec*Zeta[3])/5) + 
  g2^4*he*(-87*hec - 18*hec*Zeta[3]) + g1^2*MatMul[he, Adj[Ye], Ye]*
   ((18*hec)/5 + (108*hec*Zeta[3])/5) + g1^2*MatMul[Ye, Adj[Ye], he]*
   ((18*hec)/5 + (108*hec*Zeta[3])/5) + MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye]*
   (12*hec + 24*hec*Zeta[3]) + MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye]*
   (12*hec + 24*hec*Zeta[3]) + MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he]*
   (12*hec + 24*hec*Zeta[3]) + g2^2*M2*MatMul[Ye, Adj[Ye], Ye]*
   (-18*hec + 36*hec*Zeta[3]) + g2^4*M2*Ye*(174*hec + 36*hec*Zeta[3]) + 
  g1^2*g2^2*he*(-54*hec + (324*hec*Zeta[3])/5) + 
  g1^4*M1*Ye*((3006*hec)/25 + (324*hec*Zeta[3])/5) + 
  g2^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*(36*mh1 - 72*mh1*Zeta[3]) + 
  g1^4*MatMul[Ye, Adj[Ye]]*((-1467*mh1)/25 + (36*mh2)/25 + 
    (24*trace[mb])/25 + (72*trace[me])/25 + (36*trace[ml])/25 + 
    (12*trace[mq])/25 + (96*trace[mt])/25 - (162*mh1*Zeta[3])/5) + 
  g2^4*MatMul[Ye, Adj[Ye]]*(-99*mh1 - 12*mh2 - 12*trace[ml] - 36*trace[mq] - 
    18*mh1*Zeta[3]) + g1^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*
   ((36*mh1)/5 + (216*mh1*Zeta[3])/5) + g1^2*g2^2*MatMul[Ye, Adj[Ye]]*
   (-54*mh1 + (324*mh1*Zeta[3])/5) + 
  MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye]]*(36*mh1 + 72*mh1*Zeta[3]) + 
  g3^2*Ye*(-128*hec*trace[hb, Adj[Yb]] + 192*hec*trace[hb, Adj[Yb]]*
     Zeta[3]) + g3^2*M3*MatMul[Ye, Adj[Ye]]*(128*trace[hb, Adj[Yb]] + 
    128*trace[hbc, Yb] - 192*trace[hb, Adj[Yb]]*Zeta[3] - 
    192*trace[hbc, Yb]*Zeta[3]) + g3^2*MatMul[he, Adj[Ye]]*
   (-128*trace[hbc, Yb] + 192*trace[hbc, Yb]*Zeta[3]) + 
  g2^2*Ye*(54*hec*trace[hb, Adj[Yb]] + 18*hec*trace[he, Adj[Ye]] - 
    108*hec*trace[hb, Adj[Yb]]*Zeta[3] - 36*hec*trace[he, Adj[Ye]]*Zeta[3]) + 
  g1^2*Ye*((214*hec*trace[hb, Adj[Yb]])/5 + (18*hec*trace[he, Adj[Ye]])/5 + 
    (84*hec*trace[hb, Adj[Yb]]*Zeta[3])/5 + 
    (108*hec*trace[he, Adj[Ye]]*Zeta[3])/5) + g2^2*MatMul[he, Adj[Ye]]*
   (54*trace[hbc, Yb] + 18*trace[hec, Ye] - 108*trace[hbc, Yb]*Zeta[3] - 
    36*trace[hec, Ye]*Zeta[3]) + g1^2*M1*MatMul[Ye, Adj[Ye]]*
   ((-214*trace[hb, Adj[Yb]])/5 - (214*trace[hbc, Yb])/5 - 
    (18*trace[he, Adj[Ye]])/5 - (18*trace[hec, Ye])/5 - 
    (84*trace[hb, Adj[Yb]]*Zeta[3])/5 - (84*trace[hbc, Yb]*Zeta[3])/5 - 
    (108*trace[he, Adj[Ye]]*Zeta[3])/5 - (108*trace[hec, Ye]*Zeta[3])/5) + 
  g1^2*MatMul[he, Adj[Ye]]*((214*trace[hbc, Yb])/5 + (18*trace[hec, Ye])/5 + 
    (84*trace[hbc, Yb]*Zeta[3])/5 + (108*trace[hec, Ye]*Zeta[3])/5) + 
  g2^2*M2*MatMul[Ye, Adj[Ye]]*(-54*trace[hb, Adj[Yb]] - 54*trace[hbc, Yb] - 
    18*trace[he, Adj[Ye]] - 18*trace[hec, Ye] + 108*trace[hb, Adj[Yb]]*
     Zeta[3] + 108*trace[hbc, Yb]*Zeta[3] + 36*trace[he, Adj[Ye]]*Zeta[3] + 
    36*trace[hec, Ye]*Zeta[3]) + g3^2*M3*MatMul[he, Adj[Ye]]*
   (128*trace[Adj[Yb], Yb] - 192*trace[Adj[Yb], Yb]*Zeta[3]) + 
  2*g3^2*me*MatMul[Ye, Adj[Ye]]*(-64*trace[Adj[Yb], Yb] + 
    96*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*ml*Ye*Adj[Ye]*
   (-128*trace[Adj[Yb], Yb] + 192*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g3^2*M3^2*MatMul[Ye, Adj[Ye]]*(-256*trace[Adj[Yb], Yb] + 
    384*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g3^2*M3*Ye*(128*hec*trace[Adj[Yb], Yb] - 192*hec*trace[Adj[Yb], Yb]*
     Zeta[3]) + g3^2*he*(-128*hec*trace[Adj[Yb], Yb] + 
    192*hec*trace[Adj[Yb], Yb]*Zeta[3]) + g2^2*M2^2*MatMul[Ye, Adj[Ye]]*
   (108*trace[Adj[Yb], Yb] + 36*trace[Adj[Ye], Ye] - 
    216*trace[Adj[Yb], Yb]*Zeta[3] - 72*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g2^2*ml*Ye*Adj[Ye]*(54*trace[Adj[Yb], Yb] + 18*trace[Adj[Ye], Ye] - 
    108*trace[Adj[Yb], Yb]*Zeta[3] - 36*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g1^2*M1*MatMul[he, Adj[Ye]]*((-214*trace[Adj[Yb], Yb])/5 - 
    (18*trace[Adj[Ye], Ye])/5 - (84*trace[Adj[Yb], Yb]*Zeta[3])/5 - 
    (108*trace[Adj[Ye], Ye]*Zeta[3])/5) + 2*g2^2*me*MatMul[Ye, Adj[Ye]]*
   (27*trace[Adj[Yb], Yb] + 9*trace[Adj[Ye], Ye] - 
    54*trace[Adj[Yb], Yb]*Zeta[3] - 18*trace[Adj[Ye], Ye]*Zeta[3]) + 
  2*g1^2*me*MatMul[Ye, Adj[Ye]]*((107*trace[Adj[Yb], Yb])/5 + 
    (9*trace[Adj[Ye], Ye])/5 + (42*trace[Adj[Yb], Yb]*Zeta[3])/5 + 
    (54*trace[Adj[Ye], Ye]*Zeta[3])/5) + g1^2*ml*Ye*Adj[Ye]*
   ((214*trace[Adj[Yb], Yb])/5 + (18*trace[Adj[Ye], Ye])/5 + 
    (84*trace[Adj[Yb], Yb]*Zeta[3])/5 + (108*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g2^2*M2*MatMul[he, Adj[Ye]]*(-54*trace[Adj[Yb], Yb] - 
    18*trace[Adj[Ye], Ye] + 108*trace[Adj[Yb], Yb]*Zeta[3] + 
    36*trace[Adj[Ye], Ye]*Zeta[3]) + g1^2*M1^2*MatMul[Ye, Adj[Ye]]*
   ((428*trace[Adj[Yb], Yb])/5 + (36*trace[Adj[Ye], Ye])/5 + 
    (168*trace[Adj[Yb], Yb]*Zeta[3])/5 + (216*trace[Adj[Ye], Ye]*Zeta[3])/
     5) + g2^2*he*(54*hec*trace[Adj[Yb], Yb] + 18*hec*trace[Adj[Ye], Ye] - 
    108*hec*trace[Adj[Yb], Yb]*Zeta[3] - 36*hec*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g1^2*M1*Ye*((-214*hec*trace[Adj[Yb], Yb])/5 - (18*hec*trace[Adj[Ye], Ye])/
     5 - (84*hec*trace[Adj[Yb], Yb]*Zeta[3])/5 - 
    (108*hec*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*he*((214*hec*trace[Adj[Yb], Yb])/5 + (18*hec*trace[Adj[Ye], Ye])/5 + 
    (84*hec*trace[Adj[Yb], Yb]*Zeta[3])/5 + 
    (108*hec*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g2^2*M2*Ye*(-54*hec*trace[Adj[Yb], Yb] - 18*hec*trace[Adj[Ye], Ye] + 
    108*hec*trace[Adj[Yb], Yb]*Zeta[3] + 36*hec*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g3^2*MatMul[Ye, Adj[Ye]]*(-128*trace[hb, hbc] - 
    256*mh1*trace[Adj[Yb], Yb] - 128*trace[Yb, Adj[Yb], mb] - 
    128*trace[Adj[Yb], Yb, mq] + 192*trace[hb, hbc]*Zeta[3] + 
    384*mh1*trace[Adj[Yb], Yb]*Zeta[3] + 192*trace[Yb, Adj[Yb], mb]*Zeta[3] + 
    192*trace[Adj[Yb], Yb, mq]*Zeta[3]) + g2^2*MatMul[Ye, Adj[Ye]]*
   (54*trace[hb, hbc] + 18*trace[he, hec] + 108*mh1*trace[Adj[Yb], Yb] + 
    36*mh1*trace[Adj[Ye], Ye] + 54*trace[Yb, Adj[Yb], mb] + 
    18*trace[Ye, Adj[Ye], me] + 54*trace[Adj[Yb], Yb, mq] + 
    18*trace[Adj[Ye], Ye, ml] - 108*trace[hb, hbc]*Zeta[3] - 
    36*trace[he, hec]*Zeta[3] - 216*mh1*trace[Adj[Yb], Yb]*Zeta[3] - 
    72*mh1*trace[Adj[Ye], Ye]*Zeta[3] - 108*trace[Yb, Adj[Yb], mb]*Zeta[3] - 
    36*trace[Ye, Adj[Ye], me]*Zeta[3] - 108*trace[Adj[Yb], Yb, mq]*Zeta[3] - 
    36*trace[Adj[Ye], Ye, ml]*Zeta[3]) + g1^2*MatMul[Ye, Adj[Ye]]*
   ((214*trace[hb, hbc])/5 + (18*trace[he, hec])/5 + 
    (428*mh1*trace[Adj[Yb], Yb])/5 + (36*mh1*trace[Adj[Ye], Ye])/5 + 
    (214*trace[Yb, Adj[Yb], mb])/5 + (18*trace[Ye, Adj[Ye], me])/5 + 
    (214*trace[Adj[Yb], Yb, mq])/5 + (18*trace[Adj[Ye], Ye, ml])/5 + 
    (84*trace[hb, hbc]*Zeta[3])/5 + (108*trace[he, hec]*Zeta[3])/5 + 
    (168*mh1*trace[Adj[Yb], Yb]*Zeta[3])/5 + 
    (216*mh1*trace[Adj[Ye], Ye]*Zeta[3])/5 + 
    (84*trace[Yb, Adj[Yb], mb]*Zeta[3])/5 + 
    (108*trace[Ye, Adj[Ye], me]*Zeta[3])/5 + 
    (84*trace[Adj[Yb], Yb, mq]*Zeta[3])/5 + 
    (108*trace[Adj[Ye], Ye, ml]*Zeta[3])/5)}

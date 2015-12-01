{(-2*g1^2*M1^2)/15 - 6*g2^2*M2^2 - (32*g3^2*M3^2)/3 + 2*mb*Yb*Adj[Yb] + 
  2*mt*Yt*Adj[Yt] + 2*MatMul[hbc, hb] + 2*MatMul[htc, ht] + 
  2*mh1*MatMul[Adj[Yb], Yb] + 2*mq*MatMul[Adj[Yb], Yb] + 
  2*mh2*MatMul[Adj[Yt], Yt] + 2*mq*MatMul[Adj[Yt], Yt], 
 (199*g1^4*M1^2)/75 + (2*g1^2*g2^2*M1^2)/5 + (32*g1^2*g3^2*M1^2)/45 + 
  (2*g1^2*g2^2*M1*M2)/5 + (2*g1^2*g2^2*M2^2)/5 + 33*g2^4*M2^2 + 
  32*g2^2*g3^2*M2^2 + (32*g1^2*g3^2*M1*M3)/45 + 32*g2^2*g3^2*M2*M3 + 
  (32*g1^2*g3^2*M3^2)/45 + 32*g2^2*g3^2*M3^2 - (128*g3^4*M3^2)/3 + 
  (4*g1^2*mb*Yb*Adj[Yb])/5 + (8*g1^2*mt*Yt*Adj[Yt])/5 + 
  (4*g1^2*MatMul[hbc, hb])/5 - (4*g1^2*M1*MatMul[hbc, Yb])/5 + 
  (8*g1^2*MatMul[htc, ht])/5 - (8*g1^2*M1*MatMul[htc, Yt])/5 - 
  (4*g1^2*M1*MatMul[Adj[Yb], hb])/5 + (8*g1^2*M1^2*MatMul[Adj[Yb], Yb])/5 + 
  (4*g1^2*mh1*MatMul[Adj[Yb], Yb])/5 + (4*g1^2*mq*MatMul[Adj[Yb], Yb])/5 - 
  4*mq*MatMul[Adj[Yb], Yb]^2 - (8*g1^2*M1*MatMul[Adj[Yt], ht])/5 + 
  (16*g1^2*M1^2*MatMul[Adj[Yt], Yt])/5 + (8*g1^2*mh2*MatMul[Adj[Yt], Yt])/5 + 
  (8*g1^2*mq*MatMul[Adj[Yt], Yt])/5 - 4*mq*MatMul[Adj[Yt], Yt]^2 - 
  4*mb*Adj[Yb]*MatMul[Yb, Adj[Yb], Yb] - 
  4*mt*Adj[Yt]*MatMul[Yt, Adj[Yt], Yt] - 
  4*mb*Yb*MatMul[Adj[Yb], Yb, Adj[Yb]] - 
  4*mt*Yt*MatMul[Adj[Yt], Yt, Adj[Yt]] - 4*MatMul[hbc, hb, Adj[Yb], Yb] - 
  4*MatMul[hbc, Yb, Adj[Yb], hb] - 4*MatMul[htc, ht, Adj[Yt], Yt] - 
  4*MatMul[htc, Yt, Adj[Yt], ht] - 4*MatMul[Adj[Yb], hb, hbc, Yb] - 
  4*MatMul[Adj[Yb], Yb, hbc, hb] - 8*mh1*MatMul[Adj[Yb], Yb, Adj[Yb], Yb] - 
  4*mq*MatMul[Adj[Yb], Yb, Adj[Yb], Yb] - 4*MatMul[Adj[Yt], ht, htc, Yt] - 
  4*MatMul[Adj[Yt], Yt, htc, ht] - 8*mh2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt] - 
  4*mq*MatMul[Adj[Yt], Yt, Adj[Yt], Yt] + 
  g2^4*(3*mh1 + 3*mh2 + 3*trace[ml] + 9*trace[mq]) + 
  g1^4*(mh1/25 + mh2/25 + (2*trace[mb])/75 + (2*trace[me])/25 + 
    trace[ml]/25 + trace[mq]/75 + (8*trace[mt])/75) + 
  g3^4*((16*trace[mb])/3 + (32*trace[mq])/3 + (16*trace[mt])/3) + 
  MatMul[hbc, Yb]*(-6*trace[hb, Adj[Yb]] - 2*trace[he, Adj[Ye]]) + 
  MatMul[Adj[Yb], hb]*(-6*trace[hbc, Yb] - 2*trace[hec, Ye]) - 
  6*MatMul[htc, Yt]*trace[ht, Adj[Yt]] - 6*MatMul[Adj[Yt], ht]*
   trace[htc, Yt] + mb*Yb*Adj[Yb]*(-6*trace[Adj[Yb], Yb] - 
    2*trace[Adj[Ye], Ye]) + MatMul[hbc, hb]*(-6*trace[Adj[Yb], Yb] - 
    2*trace[Adj[Ye], Ye]) + 2*mq*MatMul[Adj[Yb], Yb]*
   (-3*trace[Adj[Yb], Yb] - trace[Adj[Ye], Ye]) - 
  6*mt*Yt*Adj[Yt]*trace[Adj[Yt], Yt] - 6*MatMul[htc, ht]*trace[Adj[Yt], Yt] - 
  6*mq*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  MatMul[Adj[Yb], Yb]*(-6*trace[hb, hbc] - 2*trace[he, hec] - 
    12*mh1*trace[Adj[Yb], Yb] - 4*mh1*trace[Adj[Ye], Ye] - 
    6*trace[Yb, Adj[Yb], mb] - 2*trace[Ye, Adj[Ye], me] - 
    6*trace[Adj[Yb], Yb, mq] - 2*trace[Adj[Ye], Ye, ml]) + 
  MatMul[Adj[Yt], Yt]*(-6*trace[ht, htc] - 12*mh2*trace[Adj[Yt], Yt] - 
    6*trace[Yt, Adj[Yt], mt] - 6*trace[Adj[Yt], Yt, mq]), 
 (-32*g1^2*g2^2*g3^2*M1^2)/5 - (32*g1^2*g2^2*g3^2*M1*M2)/5 - 
  (32*g1^2*g2^2*g3^2*M2^2)/5 - (32*g1^2*g2^2*g3^2*M1*M3)/5 - 
  (32*g1^2*g2^2*g3^2*M2*M3)/5 - (32*g1^2*g2^2*g3^2*M3^2)/5 - 
  8*g2^2*g3^2*mb*Yb*Adj[Yb] - 8*g2^2*g3^2*mt*Yt*Adj[Yt] - 
  8*g2^2*g3^2*MatMul[hbc, hb] + 8*g2^2*g3^2*M2*MatMul[hbc, Yb] + 
  8*g2^2*g3^2*M3*MatMul[hbc, Yb] - 8*g2^2*g3^2*MatMul[htc, ht] + 
  8*g2^2*g3^2*M2*MatMul[htc, Yt] + 8*g2^2*g3^2*M3*MatMul[htc, Yt] + 
  8*g2^2*g3^2*M2*MatMul[Adj[Yb], hb] + 8*g2^2*g3^2*M3*MatMul[Adj[Yb], hb] - 
  16*g2^2*g3^2*M2^2*MatMul[Adj[Yb], Yb] - 16*g2^2*g3^2*M2*M3*
   MatMul[Adj[Yb], Yb] - 16*g2^2*g3^2*M3^2*MatMul[Adj[Yb], Yb] - 
  8*g2^2*g3^2*mh1*MatMul[Adj[Yb], Yb] - 8*g2^2*g3^2*mq*MatMul[Adj[Yb], Yb] + 
  (128*g3^2*mq*MatMul[Adj[Yb], Yb]^2)/3 + 
  8*g2^2*g3^2*M2*MatMul[Adj[Yt], ht] + 8*g2^2*g3^2*M3*MatMul[Adj[Yt], ht] - 
  16*g2^2*g3^2*M2^2*MatMul[Adj[Yt], Yt] - 16*g2^2*g3^2*M2*M3*
   MatMul[Adj[Yt], Yt] - 16*g2^2*g3^2*M3^2*MatMul[Adj[Yt], Yt] - 
  8*g2^2*g3^2*mh2*MatMul[Adj[Yt], Yt] - 8*g2^2*g3^2*mq*MatMul[Adj[Yt], Yt] + 
  (128*g3^2*mq*MatMul[Adj[Yt], Yt]^2)/3 + 
  (128*g3^2*mb*Adj[Yb]*MatMul[Yb, Adj[Yb], Yb])/3 + 
  (128*g3^2*mt*Adj[Yt]*MatMul[Yt, Adj[Yt], Yt])/3 + 
  (128*g3^2*mb*Yb*MatMul[Adj[Yb], Yb, Adj[Yb]])/3 + 
  8*mt*MatMul[Yt, Adj[Yb], Yb]*MatMul[Adj[Yb], Yb, Adj[Yt]] + 
  8*mb*MatMul[Yb, Adj[Yt], Yt]*MatMul[Adj[Yt], Yt, Adj[Yb]] + 
  (128*g3^2*mt*Yt*MatMul[Adj[Yt], Yt, Adj[Yt]])/3 + 
  (128*g3^2*MatMul[hbc, hb, Adj[Yb], Yb])/3 + 
  (128*g3^2*MatMul[hbc, Yb, Adj[Yb], hb])/3 - 
  (128*g3^2*M3*MatMul[hbc, Yb, Adj[Yb], Yb])/3 + 
  (128*g3^2*MatMul[htc, ht, Adj[Yt], Yt])/3 + 
  (128*g3^2*MatMul[htc, Yt, Adj[Yt], ht])/3 - 
  (128*g3^2*M3*MatMul[htc, Yt, Adj[Yt], Yt])/3 + 
  (128*g3^2*MatMul[Adj[Yb], hb, hbc, Yb])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yb], hb, Adj[Yb], Yb])/3 + 
  (128*g3^2*MatMul[Adj[Yb], Yb, hbc, hb])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yb], Yb, hbc, Yb])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yb], Yb, Adj[Yb], hb])/3 + 
  (256*g3^2*M3^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb])/3 + 
  (256*g3^2*mh1*MatMul[Adj[Yb], Yb, Adj[Yb], Yb])/3 + 
  (128*g3^2*mq*MatMul[Adj[Yb], Yb, Adj[Yb], Yb])/3 + 
  8*mq*MatMul[Adj[Yb], Yb]*MatMul[Adj[Yb], Yb, Adj[Yt], Yt] + 
  8*mq*MatMul[Adj[Yt], Yt]*MatMul[Adj[Yb], Yb, Adj[Yt], Yt] + 
  (128*g3^2*MatMul[Adj[Yt], ht, htc, Yt])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yt], ht, Adj[Yt], Yt])/3 + 
  (128*g3^2*MatMul[Adj[Yt], Yt, htc, ht])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yt], Yt, htc, Yt])/3 + 
  8*mq*MatMul[Adj[Yb], Yb]*MatMul[Adj[Yt], Yt, Adj[Yb], Yb] + 
  8*mq*MatMul[Adj[Yt], Yt]*MatMul[Adj[Yt], Yt, Adj[Yb], Yb] - 
  (128*g3^2*M3*MatMul[Adj[Yt], Yt, Adj[Yt], ht])/3 + 
  (256*g3^2*M3^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt])/3 + 
  (256*g3^2*mh2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt])/3 + 
  (128*g3^2*mq*MatMul[Adj[Yt], Yt, Adj[Yt], Yt])/3 + 
  8*mb*Adj[Yb]*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  8*mt*Adj[Yt]*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt] + 
  8*mb*Yb*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb]] + 
  8*mt*Yt*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt]] + 
  8*MatMul[hbc, hb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  8*MatMul[hbc, Yb, Adj[Yt], ht, Adj[Yb], Yb] + 
  8*MatMul[hbc, Yb, Adj[Yt], Yt, Adj[Yb], hb] + 
  8*MatMul[htc, ht, Adj[Yb], Yb, Adj[Yt], Yt] + 
  8*MatMul[htc, Yt, Adj[Yb], hb, Adj[Yt], Yt] + 
  8*MatMul[htc, Yt, Adj[Yb], Yb, Adj[Yt], ht] + 
  8*MatMul[Adj[Yb], hb, htc, Yt, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yb], hb, Adj[Yt], Yt, hbc, Yb] + 
  8*MatMul[Adj[Yb], Yb, htc, ht, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yb], Yb, htc, Yt, Adj[Yb], hb] + 
  8*MatMul[Adj[Yb], Yb, Adj[Yt], ht, hbc, Yb] + 
  8*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, hbc, hb] + 
  (16*mh1 + 8*mh2)*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  8*mq*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yt], ht, hbc, Yb, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yt], ht, Adj[Yb], Yb, htc, Yt] + 
  8*MatMul[Adj[Yt], Yt, hbc, hb, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yt], Yt, hbc, Yb, Adj[Yt], ht] + 
  8*MatMul[Adj[Yt], Yt, Adj[Yb], hb, htc, Yt] + 
  8*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, htc, ht] + 
  (8*mh1 + 16*mh2)*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], Yt] + 
  8*mq*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], Yt] + 
  g2^4*g3^2*(-16*mh1 - 16*mh2 - 16*trace[ml] - 48*trace[mq]) + 
  g2^6*(-3*mh1 - 3*mh2 - 3*trace[ml] - 9*trace[mq]) + 
  g1^2*g2^4*(-mh1/5 - mh2/5 - trace[ml]/5 - (3*trace[mq])/5) + 
  g2^2*g3^4*(-16*trace[mb] - 32*trace[mq] - 16*trace[mt]) + 
  g1^6*((-199*mh1)/375 - (199*mh2)/375 - (398*trace[mb])/1125 - 
    (398*trace[me])/375 - (199*trace[ml])/375 - (199*trace[mq])/1125 - 
    (1592*trace[mt])/1125) + g1^4*g3^2*((-16*mh1)/75 - (16*mh2)/75 - 
    (32*trace[mb])/225 - (32*trace[me])/75 - (16*trace[ml])/75 - 
    (16*trace[mq])/225 - (128*trace[mt])/225) + 
  g1^2*g3^4*((-16*trace[mb])/45 - (32*trace[mq])/45 - (16*trace[mt])/45) + 
  g1^4*g2^2*((-3*mh1)/25 - (3*mh2)/25 - (2*trace[mb])/25 - (6*trace[me])/25 - 
    (3*trace[ml])/25 - trace[mq]/25 - (8*trace[mt])/25) + 
  g3^6*((320*trace[mb])/9 + (640*trace[mq])/9 + (320*trace[mt])/9) + 
  MatMul[hbc, Yb, Adj[Yb], Yb]*(12*trace[hb, Adj[Yb]] + 
    4*trace[he, Adj[Ye]]) + MatMul[Adj[Yb], Yb, hbc, Yb]*
   (12*trace[hb, Adj[Yb]] + 4*trace[he, Adj[Ye]]) + 
  g2^2*MatMul[hbc, Yb]*(36*trace[hb, Adj[Yb]] + 12*trace[he, Adj[Ye]]) + 
  g2^2*M2*MatMul[Adj[Yb], Yb]*(-36*trace[hb, Adj[Yb]] - 36*trace[hbc, Yb] - 
    12*trace[he, Adj[Ye]] - 12*trace[hec, Ye]) + 
  MatMul[Adj[Yb], hb, Adj[Yb], Yb]*(12*trace[hbc, Yb] + 4*trace[hec, Ye]) + 
  MatMul[Adj[Yb], Yb, Adj[Yb], hb]*(12*trace[hbc, Yb] + 4*trace[hec, Ye]) + 
  g2^2*MatMul[Adj[Yb], hb]*(36*trace[hbc, Yb] + 12*trace[hec, Ye]) + 
  36*g2^2*MatMul[htc, Yt]*trace[ht, Adj[Yt]] + 
  12*MatMul[htc, Yt, Adj[Yt], Yt]*trace[ht, Adj[Yt]] + 
  12*MatMul[Adj[Yt], Yt, htc, Yt]*trace[ht, Adj[Yt]] + 
  g2^2*M2*MatMul[Adj[Yt], Yt]*(-36*trace[ht, Adj[Yt]] - 36*trace[htc, Yt]) + 
  36*g2^2*MatMul[Adj[Yt], ht]*trace[htc, Yt] + 
  12*MatMul[Adj[Yt], ht, Adj[Yt], Yt]*trace[htc, Yt] + 
  12*MatMul[Adj[Yt], Yt, Adj[Yt], ht]*trace[htc, Yt] + 
  g1^4*M1*((14*trace[hb, Adj[Yb]])/15 + (14*trace[hbc, Yb])/15 + 
    (6*trace[he, Adj[Ye]])/5 + (6*trace[hec, Ye])/5 + 
    (26*trace[ht, Adj[Yt]])/15 + (26*trace[htc, Yt])/15) + 
  g2^4*M2*(90*trace[hb, Adj[Yb]] + 90*trace[hbc, Yb] + 
    30*trace[he, Adj[Ye]] + 30*trace[hec, Ye] + 90*trace[ht, Adj[Yt]] + 
    90*trace[htc, Yt]) + g3^4*M3*((320*trace[hb, Adj[Yb]])/3 + 
    (320*trace[hbc, Yb])/3 + (320*trace[ht, Adj[Yt]])/3 + 
    (320*trace[htc, Yt])/3) + g2^2*M2*MatMul[hbc, Yb]*
   (-36*trace[Adj[Yb], Yb] - 12*trace[Adj[Ye], Ye]) + 
  g2^2*M2*MatMul[Adj[Yb], hb]*(-36*trace[Adj[Yb], Yb] - 
    12*trace[Adj[Ye], Ye]) + 2*mq*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*
   (6*trace[Adj[Yb], Yb] + 2*trace[Adj[Ye], Ye]) + 
  mq*MatMul[Adj[Yb], Yb]^2*(12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  mb*Adj[Yb]*MatMul[Yb, Adj[Yb], Yb]*(12*trace[Adj[Yb], Yb] + 
    4*trace[Adj[Ye], Ye]) + mb*Yb*MatMul[Adj[Yb], Yb, Adj[Yb]]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[hbc, hb, Adj[Yb], Yb]*(12*trace[Adj[Yb], Yb] + 
    4*trace[Adj[Ye], Ye]) + MatMul[hbc, Yb, Adj[Yb], hb]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Yb], hb, hbc, Yb]*(12*trace[Adj[Yb], Yb] + 
    4*trace[Adj[Ye], Ye]) + MatMul[Adj[Yb], Yb, hbc, hb]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  2*g2^2*mq*MatMul[Adj[Yb], Yb]*(18*trace[Adj[Yb], Yb] + 
    6*trace[Adj[Ye], Ye]) + g2^2*mb*Yb*Adj[Yb]*(36*trace[Adj[Yb], Yb] + 
    12*trace[Adj[Ye], Ye]) + g2^2*MatMul[hbc, hb]*
   (36*trace[Adj[Yb], Yb] + 12*trace[Adj[Ye], Ye]) + 
  g2^2*M2^2*MatMul[Adj[Yb], Yb]*(72*trace[Adj[Yb], Yb] + 
    24*trace[Adj[Ye], Ye]) + g3^4*M3^2*(-320*trace[Adj[Yb], Yb] - 
    320*trace[Adj[Yt], Yt]) + g2^4*M2^2*(-270*trace[Adj[Yb], Yb] - 
    90*trace[Adj[Ye], Ye] - 270*trace[Adj[Yt], Yt]) + 
  g1^4*M1^2*((-14*trace[Adj[Yb], Yb])/5 - (18*trace[Adj[Ye], Ye])/5 - 
    (26*trace[Adj[Yt], Yt])/5) + 36*g2^2*mt*Yt*Adj[Yt]*trace[Adj[Yt], Yt] + 
  36*g2^2*MatMul[htc, ht]*trace[Adj[Yt], Yt] - 36*g2^2*M2*MatMul[htc, Yt]*
   trace[Adj[Yt], Yt] - 36*g2^2*M2*MatMul[Adj[Yt], ht]*trace[Adj[Yt], Yt] + 
  72*g2^2*M2^2*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  36*g2^2*mq*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*mq*MatMul[Adj[Yt], Yt]^2*trace[Adj[Yt], Yt] + 
  12*mt*Adj[Yt]*MatMul[Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*mt*Yt*MatMul[Adj[Yt], Yt, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  12*MatMul[htc, ht, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*MatMul[htc, Yt, Adj[Yt], ht]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[Yt], ht, htc, Yt]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[Yt], Yt, htc, ht]*trace[Adj[Yt], Yt] + 
  12*mq*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*(12*trace[hb, hbc] + 4*trace[he, hec] + 
    36*mh1*trace[Adj[Yb], Yb] + 12*mh1*trace[Adj[Ye], Ye] + 
    12*trace[Yb, Adj[Yb], mb] + 4*trace[Ye, Adj[Ye], me] + 
    12*trace[Adj[Yb], Yb, mq] + 4*trace[Adj[Ye], Ye, ml]) + 
  g2^2*MatMul[Adj[Yb], Yb]*(36*trace[hb, hbc] + 12*trace[he, hec] + 
    72*mh1*trace[Adj[Yb], Yb] + 24*mh1*trace[Adj[Ye], Ye] + 
    36*trace[Yb, Adj[Yb], mb] + 12*trace[Ye, Adj[Ye], me] + 
    36*trace[Adj[Yb], Yb, mq] + 12*trace[Adj[Ye], Ye, ml]) + 
  g3^4*(-64*trace[hb, hbc] - 64*trace[ht, htc] - 64*mh1*trace[Adj[Yb], Yb] - 
    64*mh2*trace[Adj[Yt], Yt] - 64*trace[Yb, Adj[Yb], mb] - 
    64*trace[Yt, Adj[Yt], mt] - 64*trace[Adj[Yb], Yb, mq] - 
    64*trace[Adj[Yt], Yt, mq]) + 
  g2^4*(-54*trace[hb, hbc] - 18*trace[he, hec] - 54*trace[ht, htc] - 
    54*mh1*trace[Adj[Yb], Yb] - 18*mh1*trace[Adj[Ye], Ye] - 
    54*mh2*trace[Adj[Yt], Yt] - 54*trace[Yb, Adj[Yb], mb] - 
    18*trace[Ye, Adj[Ye], me] - 54*trace[Yt, Adj[Yt], mt] - 
    54*trace[Adj[Yb], Yb, mq] - 18*trace[Adj[Ye], Ye, ml] - 
    54*trace[Adj[Yt], Yt, mq]) + g1^4*((-14*trace[hb, hbc])/25 - 
    (18*trace[he, hec])/25 - (26*trace[ht, htc])/25 - 
    (14*mh1*trace[Adj[Yb], Yb])/25 - (18*mh1*trace[Adj[Ye], Ye])/25 - 
    (26*mh2*trace[Adj[Yt], Yt])/25 - (14*trace[Yb, Adj[Yb], mb])/25 - 
    (18*trace[Ye, Adj[Ye], me])/25 - (26*trace[Yt, Adj[Yt], mt])/25 - 
    (14*trace[Adj[Yb], Yb, mq])/25 - (18*trace[Adj[Ye], Ye, ml])/25 - 
    (26*trace[Adj[Yt], Yt, mq])/25) + MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*
   (12*trace[ht, htc] + 36*mh2*trace[Adj[Yt], Yt] + 
    12*trace[Yt, Adj[Yt], mt] + 12*trace[Adj[Yt], Yt, mq]) + 
  g2^2*MatMul[Adj[Yt], Yt]*(36*trace[ht, htc] + 72*mh2*trace[Adj[Yt], Yt] + 
    36*trace[Yt, Adj[Yt], mt] + 36*trace[Adj[Yt], Yt, mq]) + 
  MatMul[hbc, Yb]*(-36*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
    12*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb] - 12*trace[hb, Adj[Yb]]*
     trace[Adj[Ye], Ye] - 4*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] + 
    12*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 72*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
    24*trace[Adj[Ye], he, Adj[Ye], Ye] + 
    12*trace[Adj[Yt], ht, Adj[Yb], Yb]) + 
  MatMul[htc, Yt]*(-36*trace[ht, Adj[Yt]]*trace[Adj[Yt], Yt] + 
    12*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 12*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
    72*trace[Adj[Yt], ht, Adj[Yt], Yt]) + MatMul[Adj[Yb], hb]*
   (-36*trace[hbc, Yb]*trace[Adj[Yb], Yb] - 12*trace[hec, Ye]*
     trace[Adj[Yb], Yb] - 12*trace[hbc, Yb]*trace[Adj[Ye], Ye] - 
    4*trace[hec, Ye]*trace[Adj[Ye], Ye] + 12*trace[htc, Yt, Adj[Yb], Yb] + 
    72*trace[Adj[Yb], Yb, hbc, Yb] + 24*trace[Adj[Ye], Ye, hec, Ye] + 
    12*trace[Adj[Yt], Yt, hbc, Yb]) + MatMul[Adj[Yt], ht]*
   (-36*trace[htc, Yt]*trace[Adj[Yt], Yt] + 12*trace[htc, Yt, Adj[Yb], Yb] + 
    12*trace[Adj[Yt], Yt, hbc, Yb] + 72*trace[Adj[Yt], Yt, htc, Yt]) + 
  2*mq*MatMul[Adj[Yb], Yb]*(-9*trace[Adj[Yb], Yb]^2 - 
    6*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - trace[Adj[Ye], Ye]^2 + 
    18*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 6*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    6*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + 
  mb*Yb*Adj[Yb]*(-18*trace[Adj[Yb], Yb]^2 - 12*trace[Adj[Yb], Yb]*
     trace[Adj[Ye], Ye] - 2*trace[Adj[Ye], Ye]^2 + 
    36*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 12*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    12*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + 
  MatMul[hbc, hb]*(-18*trace[Adj[Yb], Yb]^2 - 12*trace[Adj[Yb], Yb]*
     trace[Adj[Ye], Ye] - 2*trace[Adj[Ye], Ye]^2 + 
    36*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 12*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    12*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + 2*mq*MatMul[Adj[Yt], Yt]*
   (-9*trace[Adj[Yt], Yt]^2 + 6*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    18*trace[Adj[Yt], Yt, Adj[Yt], Yt]) + 
  mt*Yt*Adj[Yt]*(-18*trace[Adj[Yt], Yt]^2 + 
    12*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    36*trace[Adj[Yt], Yt, Adj[Yt], Yt]) + 
  MatMul[htc, ht]*(-18*trace[Adj[Yt], Yt]^2 + 
    12*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    36*trace[Adj[Yt], Yt, Adj[Yt], Yt]) + MatMul[Adj[Yb], Yb]*
   (-36*trace[hb, Adj[Yb]]*trace[hbc, Yb] - 12*trace[hbc, Yb]*
     trace[he, Adj[Ye]] - 12*trace[hb, Adj[Yb]]*trace[hec, Ye] - 
    4*trace[he, Adj[Ye]]*trace[hec, Ye] - 36*trace[hb, hbc]*
     trace[Adj[Yb], Yb] - 12*trace[he, hec]*trace[Adj[Yb], Yb] - 
    54*mh1*trace[Adj[Yb], Yb]^2 - 12*trace[hb, hbc]*trace[Adj[Ye], Ye] - 
    4*trace[he, hec]*trace[Adj[Ye], Ye] - 36*mh1*trace[Adj[Yb], Yb]*
     trace[Adj[Ye], Ye] - 6*mh1*trace[Adj[Ye], Ye]^2 - 
    36*trace[Adj[Yb], Yb]*trace[Yb, Adj[Yb], mb] - 
    12*trace[Adj[Ye], Ye]*trace[Yb, Adj[Yb], mb] - 
    12*trace[Adj[Yb], Yb]*trace[Ye, Adj[Ye], me] - 
    4*trace[Adj[Ye], Ye]*trace[Ye, Adj[Ye], me] - 
    36*trace[Adj[Yb], Yb]*trace[Adj[Yb], Yb, mq] - 
    12*trace[Adj[Ye], Ye]*trace[Adj[Yb], Yb, mq] - 
    12*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye, ml] - 
    4*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, ml] + 
    12*trace[hb, htc, Yt, Adj[Yb]] + 72*trace[hbc, hb, Adj[Yb], Yb] + 
    12*trace[hbc, hb, Adj[Yt], Yt] + 24*trace[hec, he, Adj[Ye], Ye] + 
    12*trace[htc, ht, Adj[Yb], Yb] + 72*trace[Adj[Yb], hb, hbc, Yb] + 
    108*mh1*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    24*trace[Adj[Ye], he, hec, Ye] + 36*mh1*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    12*trace[Adj[Yt], ht, hbc, Yb] + 24*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    12*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    72*trace[Yb, Adj[Yb], Yb, Adj[Yb], mb] + 
    12*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
    24*trace[Ye, Adj[Ye], Ye, Adj[Ye], me] + 
    12*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
    72*trace[Adj[Yb], Yb, Adj[Yb], Yb, mq] + 
    12*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
    24*trace[Adj[Ye], Ye, Adj[Ye], Ye, ml] + 
    12*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq]) + 
  MatMul[Adj[Yt], Yt]*(-36*trace[ht, Adj[Yt]]*trace[htc, Yt] - 
    36*trace[ht, htc]*trace[Adj[Yt], Yt] - 54*mh2*trace[Adj[Yt], Yt]^2 - 
    36*trace[Adj[Yt], Yt]*trace[Yt, Adj[Yt], mt] - 
    36*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, mq] + 
    12*trace[hb, htc, Yt, Adj[Yb]] + 12*trace[hbc, hb, Adj[Yt], Yt] + 
    12*trace[htc, ht, Adj[Yb], Yb] + 72*trace[htc, ht, Adj[Yt], Yt] + 
    12*trace[Adj[Yt], ht, hbc, Yb] + 72*trace[Adj[Yt], ht, htc, Yt] + 
    12*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    24*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    108*mh2*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
    12*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
    12*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
    72*trace[Yt, Adj[Yt], Yt, Adj[Yt], mt] + 
    12*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
    12*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
    72*trace[Adj[Yt], Yt, Adj[Yt], Yt, mq]) + 
  g2^4*g3^2*M2^2*(664 - 1296*Zeta[3]) + g3^4*M3^2*MatMul[Adj[Yb], Yb]*
   (32 - 1088*Zeta[3]) + g3^4*M3^2*MatMul[Adj[Yt], Yt]*(32 - 1088*Zeta[3]) + 
  g2^2*g3^4*M3^2*(192 - 864*Zeta[3]) + g2^4*g3^2*M2*M3*(400 - 864*Zeta[3]) + 
  g2^2*g3^4*M2*M3*(64 - 576*Zeta[3]) + g2^4*g3^2*M3^2*(272 - 432*Zeta[3]) + 
  g2^4*M2^2*MatMul[Adj[Yb], Yb]*(-135 - 378*Zeta[3]) + 
  g2^4*M2^2*MatMul[Adj[Yt], Yt]*(-135 - 378*Zeta[3]) + 
  g2^2*g3^4*M2^2*(80 - 288*Zeta[3]) + g1^2*g3^4*M3^2*
   (2464/15 - (1056*Zeta[3])/5) + g3^4*mb*Yb*Adj[Yb]*
   (16/3 - (544*Zeta[3])/3) + g3^4*mt*Yt*Adj[Yt]*(16/3 - (544*Zeta[3])/3) + 
  g3^4*MatMul[hbc, hb]*(16/3 - (544*Zeta[3])/3) + 
  g3^4*MatMul[htc, ht]*(16/3 - (544*Zeta[3])/3) + 
  g1^2*g3^4*M1*M3*(4864/45 - (704*Zeta[3])/5) + 
  g1^2*g2^4*M2^2*(379/5 - (486*Zeta[3])/5) + 2*g3^4*mq*MatMul[Adj[Yb], Yb]*
   (8/3 - (272*Zeta[3])/3) + 2*g3^4*mq*MatMul[Adj[Yt], Yt]*
   (8/3 - (272*Zeta[3])/3) + g1^2*g3^4*M1^2*(592/9 - (352*Zeta[3])/5) + 
  g1^2*g2^4*M1*M2*(50 - (324*Zeta[3])/5) + g2^4*mb*Yb*Adj[Yb]*
   (-45/2 - 63*Zeta[3]) + g2^4*mt*Yt*Adj[Yt]*(-45/2 - 63*Zeta[3]) + 
  g2^4*MatMul[hbc, hb]*(-45/2 - 63*Zeta[3]) + 
  g2^4*MatMul[htc, ht]*(-45/2 - 63*Zeta[3]) + 
  g2^2*M2*MatMul[hbc, Yb, Adj[Yb], Yb]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[htc, Yt, Adj[Yt], Yt]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Yb], hb, Adj[Yb], Yb]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Yb], Yb, hbc, Yb]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Yb], Yb, Adj[Yb], hb]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Yt], ht, Adj[Yt], Yt]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Yt], Yt, htc, Yt]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Yt], Yt, Adj[Yt], ht]*(6 - 36*Zeta[3]) + 
  g1^2*g2^4*M1^2*(152/5 - (162*Zeta[3])/5) + 2*g2^4*mq*MatMul[Adj[Yb], Yb]*
   (-45/4 - (63*Zeta[3])/2) + 2*g2^4*mq*MatMul[Adj[Yt], Yt]*
   (-45/4 - (63*Zeta[3])/2) + g1^2*M1^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*
   (44/3 - 24*Zeta[3]) + g1^4*g3^2*M1^2*(776/75 - (528*Zeta[3])/25) + 
  g1^6*M1^2*(57511/1125 - (2388*Zeta[3])/125) + 
  g1^2*g2^2*M1*MatMul[htc, Yt]*(59/5 - 18*Zeta[3]) + 
  g1^2*g2^2*M2*MatMul[htc, Yt]*(59/5 - 18*Zeta[3]) + 
  g1^2*g2^2*M1*MatMul[Adj[Yt], ht]*(59/5 - 18*Zeta[3]) + 
  g1^2*g2^2*M2*MatMul[Adj[Yt], ht]*(59/5 - 18*Zeta[3]) + 
  g1^2*g3^2*M1*MatMul[hbc, Yb]*(152/15 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[hbc, Yb]*(152/15 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M1*MatMul[Adj[Yb], hb]*(152/15 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[Adj[Yb], hb]*(152/15 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M1*MatMul[htc, Yt]*(136/5 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[htc, Yt]*(136/5 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M1*MatMul[Adj[Yt], ht]*(136/5 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[Adj[Yt], ht]*(136/5 - (256*Zeta[3])/15) + 
  g1^4*g3^2*M1*M3*(1552/225 - (352*Zeta[3])/25) + 
  g1^2*mq*MatMul[Adj[Yt], Yt]^2*(22/3 - 12*Zeta[3]) + 
  g1^2*mt*Adj[Yt]*MatMul[Yt, Adj[Yt], Yt]*(22/3 - 12*Zeta[3]) + 
  g1^2*mt*Yt*MatMul[Adj[Yt], Yt, Adj[Yt]]*(22/3 - 12*Zeta[3]) + 
  g1^2*MatMul[htc, ht, Adj[Yt], Yt]*(22/3 - 12*Zeta[3]) + 
  g1^2*MatMul[htc, Yt, Adj[Yt], ht]*(22/3 - 12*Zeta[3]) + 
  g1^2*MatMul[Adj[Yt], ht, htc, Yt]*(22/3 - 12*Zeta[3]) + 
  g1^2*MatMul[Adj[Yt], Yt, htc, ht]*(22/3 - 12*Zeta[3]) + 
  g1^4*g3^2*M3^2*(208/45 - (176*Zeta[3])/25) + 
  g1^4*g2^2*M1^2*(33/25 - (162*Zeta[3])/25) + 
  2*g1^2*mq*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*(11/3 - 6*Zeta[3]) + 
  g1^2*M1^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*(28/15 - (24*Zeta[3])/5) + 
  g1^4*g2^2*M1*M2*(22/25 - (108*Zeta[3])/25) + 
  g1^4*M1*MatMul[htc, Yt]*(3767/75 - (286*Zeta[3])/75) + 
  g1^4*M1*MatMul[Adj[Yt], ht]*(3767/75 - (286*Zeta[3])/75) + 
  g1^2*g2^2*M1*MatMul[hbc, Yb]*(41/5 - (18*Zeta[3])/5) + 
  g1^2*g2^2*M2*MatMul[hbc, Yb]*(41/5 - (18*Zeta[3])/5) + 
  g1^2*g2^2*M1*MatMul[Adj[Yb], hb]*(41/5 - (18*Zeta[3])/5) + 
  g1^2*g2^2*M2*MatMul[Adj[Yb], hb]*(41/5 - (18*Zeta[3])/5) + 
  g1^2*mq*MatMul[Adj[Yb], Yb]^2*(14/15 - (12*Zeta[3])/5) + 
  g1^2*mb*Adj[Yb]*MatMul[Yb, Adj[Yb], Yb]*(14/15 - (12*Zeta[3])/5) + 
  g1^2*mb*Yb*MatMul[Adj[Yb], Yb, Adj[Yb]]*(14/15 - (12*Zeta[3])/5) + 
  g1^2*MatMul[hbc, hb, Adj[Yb], Yb]*(14/15 - (12*Zeta[3])/5) + 
  g1^2*MatMul[hbc, Yb, Adj[Yb], hb]*(14/15 - (12*Zeta[3])/5) + 
  g1^2*MatMul[Adj[Yb], hb, hbc, Yb]*(14/15 - (12*Zeta[3])/5) + 
  g1^2*MatMul[Adj[Yb], Yb, hbc, hb]*(14/15 - (12*Zeta[3])/5) + 
  g1^4*g2^2*M2^2*(4/5 - (54*Zeta[3])/25) + 
  2*g1^2*mq*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*(7/15 - (6*Zeta[3])/5) + 
  g1^4*M1*MatMul[hbc, Yb]*(633/25 - (14*Zeta[3])/15) + 
  g1^4*M1*MatMul[Adj[Yb], hb]*(633/25 - (14*Zeta[3])/15) + 
  2*g1^4*mq*MatMul[Adj[Yb], Yb]*(-633/100 + (7*Zeta[3])/30) + 
  g1^4*mb*Yb*Adj[Yb]*(-633/50 + (7*Zeta[3])/15) + 
  g1^4*MatMul[hbc, hb]*(-633/50 + (7*Zeta[3])/15) + 
  2*g1^4*mq*MatMul[Adj[Yt], Yt]*(-3767/300 + (143*Zeta[3])/150) + 
  12*mb*MatMul[Yb, Adj[Yb], Yb]*MatMul[Adj[Yb], Yb, Adj[Yb]]*Zeta[3] + 
  12*mt*MatMul[Yt, Adj[Yt], Yt]*MatMul[Adj[Yt], Yt, Adj[Yt]]*Zeta[3] + 
  24*mq*MatMul[Adj[Yb], Yb]*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
  24*mq*MatMul[Adj[Yt], Yt]*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] + 
  12*mb*Adj[Yb]*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
  12*mt*Adj[Yt]*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] + 
  12*mb*Yb*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb]]*Zeta[3] + 
  12*mt*Yt*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt]]*Zeta[3] + 
  12*MatMul[hbc, hb, Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
  12*MatMul[hbc, Yb, Adj[Yb], hb, Adj[Yb], Yb]*Zeta[3] + 
  12*MatMul[hbc, Yb, Adj[Yb], Yb, Adj[Yb], hb]*Zeta[3] + 
  12*MatMul[htc, ht, Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] + 
  12*MatMul[htc, Yt, Adj[Yt], ht, Adj[Yt], Yt]*Zeta[3] + 
  12*MatMul[htc, Yt, Adj[Yt], Yt, Adj[Yt], ht]*Zeta[3] + 
  12*MatMul[Adj[Yb], hb, hbc, Yb, Adj[Yb], Yb]*Zeta[3] + 
  12*MatMul[Adj[Yb], hb, Adj[Yb], Yb, hbc, Yb]*Zeta[3] + 
  12*MatMul[Adj[Yb], Yb, hbc, hb, Adj[Yb], Yb]*Zeta[3] + 
  12*MatMul[Adj[Yb], Yb, hbc, Yb, Adj[Yb], hb]*Zeta[3] + 
  12*MatMul[Adj[Yb], Yb, Adj[Yb], hb, hbc, Yb]*Zeta[3] + 
  12*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, hbc, hb]*Zeta[3] + 
  36*mh1*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
  12*mq*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
  12*MatMul[Adj[Yt], ht, htc, Yt, Adj[Yt], Yt]*Zeta[3] + 
  12*MatMul[Adj[Yt], ht, Adj[Yt], Yt, htc, Yt]*Zeta[3] + 
  12*MatMul[Adj[Yt], Yt, htc, ht, Adj[Yt], Yt]*Zeta[3] + 
  12*MatMul[Adj[Yt], Yt, htc, Yt, Adj[Yt], ht]*Zeta[3] + 
  12*MatMul[Adj[Yt], Yt, Adj[Yt], ht, htc, Yt]*Zeta[3] + 
  12*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, htc, ht]*Zeta[3] + 
  36*mh2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] + 
  12*mq*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] + 
  2*g1^2*g2^2*mq*MatMul[Adj[Yb], Yb]*(-41/10 + (9*Zeta[3])/5) + 
  g1^4*mt*Yt*Adj[Yt]*(-3767/150 + (143*Zeta[3])/75) + 
  g1^4*MatMul[htc, ht]*(-3767/150 + (143*Zeta[3])/75) + 
  g1^2*M1*MatMul[hbc, Yb, Adj[Yb], Yb]*(-14/15 + (12*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Yb], hb, Adj[Yb], Yb]*(-14/15 + (12*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Yb], Yb, hbc, Yb]*(-14/15 + (12*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Yb], Yb, Adj[Yb], hb]*(-14/15 + (12*Zeta[3])/5) + 
  g1^4*M1^2*MatMul[Adj[Yb], Yb]*(-1899/25 + (14*Zeta[3])/5) + 
  g1^2*g2^2*mb*Yb*Adj[Yb]*(-41/5 + (18*Zeta[3])/5) + 
  g1^2*g2^2*MatMul[hbc, hb]*(-41/5 + (18*Zeta[3])/5) + 
  g1^2*g2^2*M1^2*MatMul[Adj[Yb], Yb]*(-82/5 + (36*Zeta[3])/5) + 
  g1^2*g2^2*M1*M2*MatMul[Adj[Yb], Yb]*(-82/5 + (36*Zeta[3])/5) + 
  g1^2*g2^2*M2^2*MatMul[Adj[Yb], Yb]*(-82/5 + (36*Zeta[3])/5) + 
  2*g1^2*g3^2*mq*MatMul[Adj[Yt], Yt]*(-68/5 + (128*Zeta[3])/15) + 
  2*g1^2*g3^2*mq*MatMul[Adj[Yb], Yb]*(-76/15 + (128*Zeta[3])/15) + 
  2*g1^2*g2^2*mq*MatMul[Adj[Yt], Yt]*(-59/10 + 9*Zeta[3]) + 
  g1^4*M1^2*MatMul[Adj[Yt], Yt]*(-3767/25 + (286*Zeta[3])/25) + 
  g1^2*M1*MatMul[htc, Yt, Adj[Yt], Yt]*(-22/3 + 12*Zeta[3]) + 
  g1^2*M1*MatMul[Adj[Yt], ht, Adj[Yt], Yt]*(-22/3 + 12*Zeta[3]) + 
  g1^2*M1*MatMul[Adj[Yt], Yt, htc, Yt]*(-22/3 + 12*Zeta[3]) + 
  g1^2*M1*MatMul[Adj[Yt], Yt, Adj[Yt], ht]*(-22/3 + 12*Zeta[3]) + 
  g1^2*g3^2*mt*Yt*Adj[Yt]*(-136/5 + (256*Zeta[3])/15) + 
  g1^2*g3^2*MatMul[htc, ht]*(-136/5 + (256*Zeta[3])/15) + 
  g1^2*g3^2*mb*Yb*Adj[Yb]*(-152/15 + (256*Zeta[3])/15) + 
  g1^2*g3^2*MatMul[hbc, hb]*(-152/15 + (256*Zeta[3])/15) + 
  g1^2*g2^2*mt*Yt*Adj[Yt]*(-59/5 + 18*Zeta[3]) + 
  g1^2*g2^2*MatMul[htc, ht]*(-59/5 + 18*Zeta[3]) + 
  2*g2^2*mq*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*(-3 + 18*Zeta[3]) + 
  2*g2^2*mq*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*(-3 + 18*Zeta[3]) + 
  g1^2*g3^2*M1^2*MatMul[Adj[Yt], Yt]*(-272/5 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M1*M3*MatMul[Adj[Yt], Yt]*(-272/5 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M3^2*MatMul[Adj[Yt], Yt]*(-272/5 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M1^2*MatMul[Adj[Yb], Yb]*(-304/15 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M1*M3*MatMul[Adj[Yb], Yb]*(-304/15 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M3^2*MatMul[Adj[Yb], Yb]*(-304/15 + (512*Zeta[3])/15) + 
  g1^2*g2^2*M1^2*MatMul[Adj[Yt], Yt]*(-118/5 + 36*Zeta[3]) + 
  g1^2*g2^2*M1*M2*MatMul[Adj[Yt], Yt]*(-118/5 + 36*Zeta[3]) + 
  g1^2*g2^2*M2^2*MatMul[Adj[Yt], Yt]*(-118/5 + 36*Zeta[3]) + 
  g2^2*mq*MatMul[Adj[Yb], Yb]^2*(-6 + 36*Zeta[3]) + 
  g2^2*mq*MatMul[Adj[Yt], Yt]^2*(-6 + 36*Zeta[3]) + 
  g2^2*mb*Adj[Yb]*MatMul[Yb, Adj[Yb], Yb]*(-6 + 36*Zeta[3]) + 
  g2^2*mt*Adj[Yt]*MatMul[Yt, Adj[Yt], Yt]*(-6 + 36*Zeta[3]) + 
  g2^2*mb*Yb*MatMul[Adj[Yb], Yb, Adj[Yb]]*(-6 + 36*Zeta[3]) + 
  g2^2*mt*Yt*MatMul[Adj[Yt], Yt, Adj[Yt]]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[hbc, hb, Adj[Yb], Yb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[hbc, Yb, Adj[Yb], hb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[htc, ht, Adj[Yt], Yt]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[htc, Yt, Adj[Yt], ht]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Yb], hb, hbc, Yb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Yb], Yb, hbc, hb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Yt], ht, htc, Yt]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Yt], Yt, htc, ht]*(-6 + 36*Zeta[3]) + 
  g2^2*M2^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*(-12 + 72*Zeta[3]) + 
  g2^2*M2^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*(-12 + 72*Zeta[3]) + 
  g2^4*M2*MatMul[hbc, Yb]*(45 + 126*Zeta[3]) + 
  g2^4*M2*MatMul[htc, Yt]*(45 + 126*Zeta[3]) + g2^4*M2*MatMul[Adj[Yb], hb]*
   (45 + 126*Zeta[3]) + g2^4*M2*MatMul[Adj[Yt], ht]*(45 + 126*Zeta[3]) + 
  g3^4*M3*MatMul[hbc, Yb]*(-32/3 + (1088*Zeta[3])/3) + 
  g3^4*M3*MatMul[htc, Yt]*(-32/3 + (1088*Zeta[3])/3) + 
  g3^4*M3*MatMul[Adj[Yb], hb]*(-32/3 + (1088*Zeta[3])/3) + 
  g3^4*M3*MatMul[Adj[Yt], ht]*(-32/3 + (1088*Zeta[3])/3) + 
  g2^6*M2^2*(2157 + 3780*Zeta[3]) + g3^6*M3^2*(20512/9 + 7680*Zeta[3]) + 
  g3^4*MatMul[Adj[Yb], Yb]*((16*mh1)/3 - (544*mh1*Zeta[3])/3) + 
  g2^4*MatMul[Adj[Yb], Yb]*((-45*mh1)/2 - 63*mh1*Zeta[3]) + 
  g1^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*((28*mh1)/15 - (24*mh1*Zeta[3])/5) + 
  g1^4*MatMul[Adj[Yb], Yb]*((-657*mh1)/50 - (12*mh2)/25 - (8*trace[mb])/25 - 
    (24*trace[me])/25 - (12*trace[ml])/25 - (4*trace[mq])/25 - 
    (32*trace[mt])/25 + (7*mh1*Zeta[3])/15) + g1^2*g2^2*MatMul[Adj[Yb], Yb]*
   ((-41*mh1)/5 + (18*mh1*Zeta[3])/5) + g1^2*g3^2*MatMul[Adj[Yb], Yb]*
   ((-152*mh1)/15 + (256*mh1*Zeta[3])/15) + 
  g2^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*(-12*mh1 + 72*mh1*Zeta[3]) + 
  g3^4*MatMul[Adj[Yt], Yt]*((16*mh2)/3 - (544*mh2*Zeta[3])/3) + 
  g2^4*MatMul[Adj[Yt], Yt]*((-45*mh2)/2 - 63*mh2*Zeta[3]) + 
  g1^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*((44*mh2)/3 - 24*mh2*Zeta[3]) + 
  g1^4*MatMul[Adj[Yt], Yt]*((-24*mh1)/25 - (3911*mh2)/150 - 
    (16*trace[mb])/25 - (48*trace[me])/25 - (24*trace[ml])/25 - 
    (8*trace[mq])/25 - (64*trace[mt])/25 + (143*mh2*Zeta[3])/75) + 
  g1^2*g3^2*MatMul[Adj[Yt], Yt]*((-136*mh2)/5 + (256*mh2*Zeta[3])/15) + 
  g1^2*g2^2*MatMul[Adj[Yt], Yt]*((-59*mh2)/5 + 18*mh2*Zeta[3]) + 
  g2^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*(-12*mh2 + 72*mh2*Zeta[3]) + 
  g3^2*MatMul[hbc, Yb]*(-16*trace[hb, Adj[Yb]] + 16*trace[he, Adj[Ye]] + 
    96*trace[hb, Adj[Yb]]*Zeta[3]) + g3^2*M3*MatMul[Adj[Yb], Yb]*
   (16*trace[hb, Adj[Yb]] + 16*trace[hbc, Yb] - 16*trace[he, Adj[Ye]] - 
    16*trace[hec, Ye] - 96*trace[hb, Adj[Yb]]*Zeta[3] - 
    96*trace[hbc, Yb]*Zeta[3]) + g3^2*MatMul[Adj[Yb], hb]*
   (-16*trace[hbc, Yb] + 16*trace[hec, Ye] + 96*trace[hbc, Yb]*Zeta[3]) + 
  g1^2*MatMul[hbc, Yb]*((32*trace[hb, Adj[Yb]])/5 - 
    (16*trace[he, Adj[Ye]])/5 - (48*trace[hb, Adj[Yb]]*Zeta[3])/5 + 
    (24*trace[he, Adj[Ye]]*Zeta[3])/5) + g1^2*M1*MatMul[Adj[Yb], Yb]*
   ((-32*trace[hb, Adj[Yb]])/5 - (32*trace[hbc, Yb])/5 + 
    (16*trace[he, Adj[Ye]])/5 + (16*trace[hec, Ye])/5 + 
    (48*trace[hb, Adj[Yb]]*Zeta[3])/5 + (48*trace[hbc, Yb]*Zeta[3])/5 - 
    (24*trace[he, Adj[Ye]]*Zeta[3])/5 - (24*trace[hec, Ye]*Zeta[3])/5) + 
  g1^2*MatMul[Adj[Yb], hb]*((32*trace[hbc, Yb])/5 - (16*trace[hec, Ye])/5 - 
    (48*trace[hbc, Yb]*Zeta[3])/5 + (24*trace[hec, Ye]*Zeta[3])/5) + 
  g1^2*MatMul[htc, Yt]*(4*trace[ht, Adj[Yt]] - 
    (48*trace[ht, Adj[Yt]]*Zeta[3])/5) + g3^2*MatMul[htc, Yt]*
   (-16*trace[ht, Adj[Yt]] + 96*trace[ht, Adj[Yt]]*Zeta[3]) + 
  g3^2*M3*MatMul[Adj[Yt], Yt]*(16*trace[ht, Adj[Yt]] + 16*trace[htc, Yt] - 
    96*trace[ht, Adj[Yt]]*Zeta[3] - 96*trace[htc, Yt]*Zeta[3]) + 
  g1^2*MatMul[Adj[Yt], ht]*(4*trace[htc, Yt] - (48*trace[htc, Yt]*Zeta[3])/
     5) + g1^2*M1*MatMul[Adj[Yt], Yt]*(-4*trace[ht, Adj[Yt]] - 
    4*trace[htc, Yt] + (48*trace[ht, Adj[Yt]]*Zeta[3])/5 + 
    (48*trace[htc, Yt]*Zeta[3])/5) + g3^2*MatMul[Adj[Yt], ht]*
   (-16*trace[htc, Yt] + 96*trace[htc, Yt]*Zeta[3]) + 
  g3^2*M3*MatMul[hbc, Yb]*(16*trace[Adj[Yb], Yb] - 16*trace[Adj[Ye], Ye] - 
    96*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*M3*MatMul[Adj[Yb], hb]*
   (16*trace[Adj[Yb], Yb] - 16*trace[Adj[Ye], Ye] - 
    96*trace[Adj[Yb], Yb]*Zeta[3]) + 2*g3^2*mq*MatMul[Adj[Yb], Yb]*
   (-8*trace[Adj[Yb], Yb] + 8*trace[Adj[Ye], Ye] + 
    48*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*mb*Yb*Adj[Yb]*
   (-16*trace[Adj[Yb], Yb] + 16*trace[Adj[Ye], Ye] + 
    96*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*MatMul[hbc, hb]*
   (-16*trace[Adj[Yb], Yb] + 16*trace[Adj[Ye], Ye] + 
    96*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*M3^2*MatMul[Adj[Yb], Yb]*
   (-32*trace[Adj[Yb], Yb] + 32*trace[Adj[Ye], Ye] + 
    192*trace[Adj[Yb], Yb]*Zeta[3]) + g1^2*M1*MatMul[hbc, Yb]*
   ((-32*trace[Adj[Yb], Yb])/5 + (16*trace[Adj[Ye], Ye])/5 + 
    (48*trace[Adj[Yb], Yb]*Zeta[3])/5 - (24*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Yb], hb]*((-32*trace[Adj[Yb], Yb])/5 + 
    (16*trace[Adj[Ye], Ye])/5 + (48*trace[Adj[Yb], Yb]*Zeta[3])/5 - 
    (24*trace[Adj[Ye], Ye]*Zeta[3])/5) + 2*g1^2*mq*MatMul[Adj[Yb], Yb]*
   ((16*trace[Adj[Yb], Yb])/5 - (8*trace[Adj[Ye], Ye])/5 - 
    (24*trace[Adj[Yb], Yb]*Zeta[3])/5 + (12*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*mb*Yb*Adj[Yb]*((32*trace[Adj[Yb], Yb])/5 - (16*trace[Adj[Ye], Ye])/5 - 
    (48*trace[Adj[Yb], Yb]*Zeta[3])/5 + (24*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*MatMul[hbc, hb]*((32*trace[Adj[Yb], Yb])/5 - 
    (16*trace[Adj[Ye], Ye])/5 - (48*trace[Adj[Yb], Yb]*Zeta[3])/5 + 
    (24*trace[Adj[Ye], Ye]*Zeta[3])/5) + g1^2*M1^2*MatMul[Adj[Yb], Yb]*
   ((64*trace[Adj[Yb], Yb])/5 - (32*trace[Adj[Ye], Ye])/5 - 
    (96*trace[Adj[Yb], Yb]*Zeta[3])/5 + (48*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g3^2*M3*MatMul[htc, Yt]*(16*trace[Adj[Yt], Yt] - 
    96*trace[Adj[Yt], Yt]*Zeta[3]) + g3^2*M3*MatMul[Adj[Yt], ht]*
   (16*trace[Adj[Yt], Yt] - 96*trace[Adj[Yt], Yt]*Zeta[3]) + 
  g1^2*M1^2*MatMul[Adj[Yt], Yt]*(8*trace[Adj[Yt], Yt] - 
    (96*trace[Adj[Yt], Yt]*Zeta[3])/5) + g1^2*mt*Yt*Adj[Yt]*
   (4*trace[Adj[Yt], Yt] - (48*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  g1^2*MatMul[htc, ht]*(4*trace[Adj[Yt], Yt] - 
    (48*trace[Adj[Yt], Yt]*Zeta[3])/5) + 2*g1^2*mq*MatMul[Adj[Yt], Yt]*
   (2*trace[Adj[Yt], Yt] - (24*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  g1^2*M1*MatMul[htc, Yt]*(-4*trace[Adj[Yt], Yt] + 
    (48*trace[Adj[Yt], Yt]*Zeta[3])/5) + g1^2*M1*MatMul[Adj[Yt], ht]*
   (-4*trace[Adj[Yt], Yt] + (48*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  2*g3^2*mq*MatMul[Adj[Yt], Yt]*(-8*trace[Adj[Yt], Yt] + 
    48*trace[Adj[Yt], Yt]*Zeta[3]) + g3^2*mt*Yt*Adj[Yt]*
   (-16*trace[Adj[Yt], Yt] + 96*trace[Adj[Yt], Yt]*Zeta[3]) + 
  g3^2*MatMul[htc, ht]*(-16*trace[Adj[Yt], Yt] + 96*trace[Adj[Yt], Yt]*
     Zeta[3]) + g3^2*M3^2*MatMul[Adj[Yt], Yt]*(-32*trace[Adj[Yt], Yt] + 
    192*trace[Adj[Yt], Yt]*Zeta[3]) + g3^2*MatMul[Adj[Yb], Yb]*
   (-16*trace[hb, hbc] + 16*trace[he, hec] - 32*mh1*trace[Adj[Yb], Yb] + 
    32*mh1*trace[Adj[Ye], Ye] - 16*trace[Yb, Adj[Yb], mb] + 
    16*trace[Ye, Adj[Ye], me] - 16*trace[Adj[Yb], Yb, mq] + 
    16*trace[Adj[Ye], Ye, ml] + 96*trace[hb, hbc]*Zeta[3] + 
    192*mh1*trace[Adj[Yb], Yb]*Zeta[3] + 96*trace[Yb, Adj[Yb], mb]*Zeta[3] + 
    96*trace[Adj[Yb], Yb, mq]*Zeta[3]) + g1^2*MatMul[Adj[Yb], Yb]*
   ((32*trace[hb, hbc])/5 - (16*trace[he, hec])/5 + 
    (64*mh1*trace[Adj[Yb], Yb])/5 - (32*mh1*trace[Adj[Ye], Ye])/5 + 
    (32*trace[Yb, Adj[Yb], mb])/5 - (16*trace[Ye, Adj[Ye], me])/5 + 
    (32*trace[Adj[Yb], Yb, mq])/5 - (16*trace[Adj[Ye], Ye, ml])/5 - 
    (48*trace[hb, hbc]*Zeta[3])/5 + (24*trace[he, hec]*Zeta[3])/5 - 
    (96*mh1*trace[Adj[Yb], Yb]*Zeta[3])/5 + 
    (48*mh1*trace[Adj[Ye], Ye]*Zeta[3])/5 - 
    (48*trace[Yb, Adj[Yb], mb]*Zeta[3])/5 + 
    (24*trace[Ye, Adj[Ye], me]*Zeta[3])/5 - 
    (48*trace[Adj[Yb], Yb, mq]*Zeta[3])/5 + 
    (24*trace[Adj[Ye], Ye, ml]*Zeta[3])/5) + g1^2*MatMul[Adj[Yt], Yt]*
   (4*trace[ht, htc] + 8*mh2*trace[Adj[Yt], Yt] + 4*trace[Yt, Adj[Yt], mt] + 
    4*trace[Adj[Yt], Yt, mq] - (48*trace[ht, htc]*Zeta[3])/5 - 
    (96*mh2*trace[Adj[Yt], Yt]*Zeta[3])/5 - 
    (48*trace[Yt, Adj[Yt], mt]*Zeta[3])/5 - 
    (48*trace[Adj[Yt], Yt, mq]*Zeta[3])/5) + g3^2*MatMul[Adj[Yt], Yt]*
   (-16*trace[ht, htc] - 32*mh2*trace[Adj[Yt], Yt] - 
    16*trace[Yt, Adj[Yt], mt] - 16*trace[Adj[Yt], Yt, mq] + 
    96*trace[ht, htc]*Zeta[3] + 192*mh2*trace[Adj[Yt], Yt]*Zeta[3] + 
    96*trace[Yt, Adj[Yt], mt]*Zeta[3] + 96*trace[Adj[Yt], Yt, mq]*Zeta[3])}

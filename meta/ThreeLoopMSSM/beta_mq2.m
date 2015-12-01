{(-2*g1^2*M1^2)/15 - 6*g2^2*M2^2 - (32*g3^2*M3^2)/3 + 2*MatMul[Adj[hb], hb] + 
  2*MatMul[Adj[ht], ht] + 2*mh1*MatMul[Adj[Yb], Yb] + 
  2*mh2*MatMul[Adj[Yt], Yt] + MatMul[mq, Adj[Yb], Yb] + 
  MatMul[mq, Adj[Yt], Yt] + 2*MatMul[Adj[Yb], mb, Yb] + 
  MatMul[Adj[Yb], Yb, mq] + 2*MatMul[Adj[Yt], mt, Yt] + 
  MatMul[Adj[Yt], Yt, mq], (199*g1^4*M1^2)/75 + (2*g1^2*g2^2*M1^2)/5 + 
  (32*g1^2*g3^2*M1^2)/45 + (2*g1^2*g2^2*M1*M2)/5 + (2*g1^2*g2^2*M2^2)/5 + 
  33*g2^4*M2^2 + 32*g2^2*g3^2*M2^2 + (32*g1^2*g3^2*M1*M3)/45 + 
  32*g2^2*g3^2*M2*M3 + (32*g1^2*g3^2*M3^2)/45 + 32*g2^2*g3^2*M3^2 - 
  (128*g3^4*M3^2)/3 + (4*g1^2*MatMul[Adj[hb], hb])/5 - 
  (4*g1^2*M1*MatMul[Adj[hb], Yb])/5 + (8*g1^2*MatMul[Adj[ht], ht])/5 - 
  (8*g1^2*M1*MatMul[Adj[ht], Yt])/5 - (4*g1^2*M1*MatMul[Adj[Yb], hb])/5 + 
  (8*g1^2*M1^2*MatMul[Adj[Yb], Yb])/5 + (4*g1^2*mh1*MatMul[Adj[Yb], Yb])/5 - 
  (8*g1^2*M1*MatMul[Adj[Yt], ht])/5 + (16*g1^2*M1^2*MatMul[Adj[Yt], Yt])/5 + 
  (8*g1^2*mh2*MatMul[Adj[Yt], Yt])/5 + (2*g1^2*MatMul[mq, Adj[Yb], Yb])/5 + 
  (4*g1^2*MatMul[mq, Adj[Yt], Yt])/5 + (4*g1^2*MatMul[Adj[Yb], mb, Yb])/5 + 
  (2*g1^2*MatMul[Adj[Yb], Yb, mq])/5 + (8*g1^2*MatMul[Adj[Yt], mt, Yt])/5 + 
  (4*g1^2*MatMul[Adj[Yt], Yt, mq])/5 - 4*MatMul[Adj[hb], hb, Adj[Yb], Yb] - 
  4*MatMul[Adj[hb], Yb, Adj[Yb], hb] - 4*MatMul[Adj[ht], ht, Adj[Yt], Yt] - 
  4*MatMul[Adj[ht], Yt, Adj[Yt], ht] - 4*MatMul[Adj[Yb], hb, Adj[hb], Yb] - 
  4*MatMul[Adj[Yb], Yb, Adj[hb], hb] - 
  8*mh1*MatMul[Adj[Yb], Yb, Adj[Yb], Yb] - 
  4*MatMul[Adj[Yt], ht, Adj[ht], Yt] - 4*MatMul[Adj[Yt], Yt, Adj[ht], ht] - 
  8*mh2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt] - 
  2*MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb] - 
  2*MatMul[mq, Adj[Yt], Yt, Adj[Yt], Yt] - 
  4*MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb] - 
  4*MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb] - 
  4*MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb] - 
  2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq] - 
  4*MatMul[Adj[Yt], mt, Yt, Adj[Yt], Yt] - 
  4*MatMul[Adj[Yt], Yt, mq, Adj[Yt], Yt] - 
  4*MatMul[Adj[Yt], Yt, Adj[Yt], mt, Yt] - 
  2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, mq] + 
  g2^4*(3*mh1 + 3*mh2 + 3*trace[ml] + 9*trace[mq]) + 
  g1^4*(mh1/25 + mh2/25 + (2*trace[mb])/75 + (2*trace[me])/25 + 
    trace[ml]/25 + trace[mq]/75 + (8*trace[mt])/75) + 
  g3^4*((16*trace[mb])/3 + (32*trace[mq])/3 + (16*trace[mt])/3) + 
  MatMul[Adj[hb], Yb]*(-6*trace[hb, Adj[Yb]] - 2*trace[he, Adj[Ye]]) - 
  6*MatMul[Adj[ht], Yt]*trace[ht, Adj[Yt]] + MatMul[Adj[Yb], hb]*
   (-6*trace[Adj[hb], Yb] - 2*trace[Adj[he], Ye]) - 
  6*MatMul[Adj[Yt], ht]*trace[Adj[ht], Yt] + MatMul[Adj[hb], hb]*
   (-6*trace[Adj[Yb], Yb] - 2*trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Yb], mb, Yb]*(-6*trace[Adj[Yb], Yb] - 2*trace[Adj[Ye], Ye]) + 
  MatMul[mq, Adj[Yb], Yb]*(-3*trace[Adj[Yb], Yb] - trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Yb], Yb, mq]*(-3*trace[Adj[Yb], Yb] - trace[Adj[Ye], Ye]) - 
  6*MatMul[Adj[ht], ht]*trace[Adj[Yt], Yt] - 3*MatMul[mq, Adj[Yt], Yt]*
   trace[Adj[Yt], Yt] - 6*MatMul[Adj[Yt], mt, Yt]*trace[Adj[Yt], Yt] - 
  3*MatMul[Adj[Yt], Yt, mq]*trace[Adj[Yt], Yt] + 
  MatMul[Adj[Yb], Yb]*(-6*trace[hb, Adj[hb]] - 2*trace[he, Adj[he]] - 
    12*mh1*trace[Adj[Yb], Yb] - 4*mh1*trace[Adj[Ye], Ye] - 
    6*trace[Yb, Adj[Yb], mb] - 2*trace[Ye, Adj[Ye], me] - 
    6*trace[Adj[Yb], Yb, mq] - 2*trace[Adj[Ye], Ye, ml]) + 
  MatMul[Adj[Yt], Yt]*(-6*trace[ht, Adj[ht]] - 12*mh2*trace[Adj[Yt], Yt] - 
    6*trace[Yt, Adj[Yt], mt] - 6*trace[Adj[Yt], Yt, mq]), 
 (-32*g1^2*g2^2*g3^2*M1^2)/5 - (32*g1^2*g2^2*g3^2*M1*M2)/5 - 
  (32*g1^2*g2^2*g3^2*M2^2)/5 - (32*g1^2*g2^2*g3^2*M1*M3)/5 - 
  (32*g1^2*g2^2*g3^2*M2*M3)/5 - (32*g1^2*g2^2*g3^2*M3^2)/5 - 
  8*g2^2*g3^2*MatMul[Adj[hb], hb] + 8*g2^2*g3^2*M2*MatMul[Adj[hb], Yb] + 
  8*g2^2*g3^2*M3*MatMul[Adj[hb], Yb] - 8*g2^2*g3^2*MatMul[Adj[ht], ht] + 
  8*g2^2*g3^2*M2*MatMul[Adj[ht], Yt] + 8*g2^2*g3^2*M3*MatMul[Adj[ht], Yt] + 
  8*g2^2*g3^2*M2*MatMul[Adj[Yb], hb] + 8*g2^2*g3^2*M3*MatMul[Adj[Yb], hb] - 
  16*g2^2*g3^2*M2^2*MatMul[Adj[Yb], Yb] - 16*g2^2*g3^2*M2*M3*
   MatMul[Adj[Yb], Yb] - 16*g2^2*g3^2*M3^2*MatMul[Adj[Yb], Yb] - 
  8*g2^2*g3^2*mh1*MatMul[Adj[Yb], Yb] + 8*g2^2*g3^2*M2*MatMul[Adj[Yt], ht] + 
  8*g2^2*g3^2*M3*MatMul[Adj[Yt], ht] - 16*g2^2*g3^2*M2^2*
   MatMul[Adj[Yt], Yt] - 16*g2^2*g3^2*M2*M3*MatMul[Adj[Yt], Yt] - 
  16*g2^2*g3^2*M3^2*MatMul[Adj[Yt], Yt] - 
  8*g2^2*g3^2*mh2*MatMul[Adj[Yt], Yt] - 4*g2^2*g3^2*MatMul[mq, Adj[Yb], Yb] - 
  4*g2^2*g3^2*MatMul[mq, Adj[Yt], Yt] - 8*g2^2*g3^2*MatMul[Adj[Yb], mb, Yb] - 
  4*g2^2*g3^2*MatMul[Adj[Yb], Yb, mq] - 8*g2^2*g3^2*MatMul[Adj[Yt], mt, Yt] - 
  4*g2^2*g3^2*MatMul[Adj[Yt], Yt, mq] + 
  (128*g3^2*MatMul[Adj[hb], hb, Adj[Yb], Yb])/3 + 
  (128*g3^2*MatMul[Adj[hb], Yb, Adj[Yb], hb])/3 - 
  (128*g3^2*M3*MatMul[Adj[hb], Yb, Adj[Yb], Yb])/3 + 
  (128*g3^2*MatMul[Adj[ht], ht, Adj[Yt], Yt])/3 + 
  (128*g3^2*MatMul[Adj[ht], Yt, Adj[Yt], ht])/3 - 
  (128*g3^2*M3*MatMul[Adj[ht], Yt, Adj[Yt], Yt])/3 + 
  (128*g3^2*MatMul[Adj[Yb], hb, Adj[hb], Yb])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yb], hb, Adj[Yb], Yb])/3 + 
  (128*g3^2*MatMul[Adj[Yb], Yb, Adj[hb], hb])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yb], Yb, Adj[hb], Yb])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yb], Yb, Adj[Yb], hb])/3 + 
  (256*g3^2*M3^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb])/3 + 
  (256*g3^2*mh1*MatMul[Adj[Yb], Yb, Adj[Yb], Yb])/3 + 
  (128*g3^2*MatMul[Adj[Yt], ht, Adj[ht], Yt])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yt], ht, Adj[Yt], Yt])/3 + 
  (128*g3^2*MatMul[Adj[Yt], Yt, Adj[ht], ht])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yt], Yt, Adj[ht], Yt])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yt], Yt, Adj[Yt], ht])/3 + 
  (256*g3^2*M3^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt])/3 + 
  (256*g3^2*mh2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt])/3 + 
  (64*g3^2*MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb])/3 + 
  (64*g3^2*MatMul[mq, Adj[Yt], Yt, Adj[Yt], Yt])/3 + 
  (128*g3^2*MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb])/3 + 
  (128*g3^2*MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb])/3 + 
  (128*g3^2*MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb])/3 + 
  (64*g3^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq])/3 + 
  (128*g3^2*MatMul[Adj[Yt], mt, Yt, Adj[Yt], Yt])/3 + 
  (128*g3^2*MatMul[Adj[Yt], Yt, mq, Adj[Yt], Yt])/3 + 
  (128*g3^2*MatMul[Adj[Yt], Yt, Adj[Yt], mt, Yt])/3 + 
  (64*g3^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, mq])/3 + 
  8*MatMul[Adj[hb], hb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  8*MatMul[Adj[hb], Yb, Adj[Yt], ht, Adj[Yb], Yb] + 
  8*MatMul[Adj[hb], Yb, Adj[Yt], Yt, Adj[Yb], hb] + 
  8*MatMul[Adj[ht], ht, Adj[Yb], Yb, Adj[Yt], Yt] + 
  8*MatMul[Adj[ht], Yt, Adj[Yb], hb, Adj[Yt], Yt] + 
  8*MatMul[Adj[ht], Yt, Adj[Yb], Yb, Adj[Yt], ht] + 
  8*MatMul[Adj[Yb], hb, Adj[ht], Yt, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yb], hb, Adj[Yt], Yt, Adj[hb], Yb] + 
  8*MatMul[Adj[Yb], Yb, Adj[ht], ht, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yb], Yb, Adj[ht], Yt, Adj[Yb], hb] + 
  8*MatMul[Adj[Yb], Yb, Adj[Yt], ht, Adj[hb], Yb] + 
  8*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[hb], hb] + 
  (16*mh1 + 8*mh2)*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yt], ht, Adj[hb], Yb, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yt], ht, Adj[Yb], Yb, Adj[ht], Yt] + 
  8*MatMul[Adj[Yt], Yt, Adj[hb], hb, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yt], Yt, Adj[hb], Yb, Adj[Yt], ht] + 
  8*MatMul[Adj[Yt], Yt, Adj[Yb], hb, Adj[ht], Yt] + 
  8*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[ht], ht] + 
  (8*mh1 + 16*mh2)*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], Yt] + 
  4*MatMul[mq, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  4*MatMul[mq, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yb], mb, Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yb], Yb, mq, Adj[Yt], Yt, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yb], Yb, Adj[Yt], mt, Yt, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, mq, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], mb, Yb] + 
  4*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
  8*MatMul[Adj[Yt], mt, Yt, Adj[Yb], Yb, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yt], Yt, mq, Adj[Yb], Yb, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yt], Yt, Adj[Yb], mb, Yb, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, mq, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], mt, Yt] + 
  4*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
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
  MatMul[Adj[hb], Yb, Adj[Yb], Yb]*(12*trace[hb, Adj[Yb]] + 
    4*trace[he, Adj[Ye]]) + MatMul[Adj[Yb], Yb, Adj[hb], Yb]*
   (12*trace[hb, Adj[Yb]] + 4*trace[he, Adj[Ye]]) + 
  g2^2*MatMul[Adj[hb], Yb]*(36*trace[hb, Adj[Yb]] + 12*trace[he, Adj[Ye]]) + 
  36*g2^2*MatMul[Adj[ht], Yt]*trace[ht, Adj[Yt]] + 
  12*MatMul[Adj[ht], Yt, Adj[Yt], Yt]*trace[ht, Adj[Yt]] + 
  12*MatMul[Adj[Yt], Yt, Adj[ht], Yt]*trace[ht, Adj[Yt]] + 
  g2^2*M2*MatMul[Adj[Yb], Yb]*(-36*trace[hb, Adj[Yb]] - 
    12*trace[he, Adj[Ye]] - 36*trace[Adj[hb], Yb] - 12*trace[Adj[he], Ye]) + 
  MatMul[Adj[Yb], hb, Adj[Yb], Yb]*(12*trace[Adj[hb], Yb] + 
    4*trace[Adj[he], Ye]) + MatMul[Adj[Yb], Yb, Adj[Yb], hb]*
   (12*trace[Adj[hb], Yb] + 4*trace[Adj[he], Ye]) + 
  g2^2*MatMul[Adj[Yb], hb]*(36*trace[Adj[hb], Yb] + 12*trace[Adj[he], Ye]) + 
  g2^2*M2*MatMul[Adj[Yt], Yt]*(-36*trace[ht, Adj[Yt]] - 
    36*trace[Adj[ht], Yt]) + 36*g2^2*MatMul[Adj[Yt], ht]*trace[Adj[ht], Yt] + 
  12*MatMul[Adj[Yt], ht, Adj[Yt], Yt]*trace[Adj[ht], Yt] + 
  12*MatMul[Adj[Yt], Yt, Adj[Yt], ht]*trace[Adj[ht], Yt] + 
  g1^4*M1*((14*trace[hb, Adj[Yb]])/15 + (6*trace[he, Adj[Ye]])/5 + 
    (26*trace[ht, Adj[Yt]])/15 + (14*trace[Adj[hb], Yb])/15 + 
    (6*trace[Adj[he], Ye])/5 + (26*trace[Adj[ht], Yt])/15) + 
  g2^4*M2*(90*trace[hb, Adj[Yb]] + 30*trace[he, Adj[Ye]] + 
    90*trace[ht, Adj[Yt]] + 90*trace[Adj[hb], Yb] + 30*trace[Adj[he], Ye] + 
    90*trace[Adj[ht], Yt]) + g3^4*M3*((320*trace[hb, Adj[Yb]])/3 + 
    (320*trace[ht, Adj[Yt]])/3 + (320*trace[Adj[hb], Yb])/3 + 
    (320*trace[Adj[ht], Yt])/3) + g2^2*M2*MatMul[Adj[hb], Yb]*
   (-36*trace[Adj[Yb], Yb] - 12*trace[Adj[Ye], Ye]) + 
  g2^2*M2*MatMul[Adj[Yb], hb]*(-36*trace[Adj[Yb], Yb] - 
    12*trace[Adj[Ye], Ye]) + MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb]*
   (6*trace[Adj[Yb], Yb] + 2*trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq]*(6*trace[Adj[Yb], Yb] + 
    2*trace[Adj[Ye], Ye]) + MatMul[Adj[hb], hb, Adj[Yb], Yb]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[Adj[hb], Yb, Adj[Yb], hb]*(12*trace[Adj[Yb], Yb] + 
    4*trace[Adj[Ye], Ye]) + MatMul[Adj[Yb], hb, Adj[hb], Yb]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Yb], Yb, Adj[hb], hb]*(12*trace[Adj[Yb], Yb] + 
    4*trace[Adj[Ye], Ye]) + MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb]*(12*trace[Adj[Yb], Yb] + 
    4*trace[Adj[Ye], Ye]) + MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  g2^2*MatMul[mq, Adj[Yb], Yb]*(18*trace[Adj[Yb], Yb] + 
    6*trace[Adj[Ye], Ye]) + g2^2*MatMul[Adj[Yb], Yb, mq]*
   (18*trace[Adj[Yb], Yb] + 6*trace[Adj[Ye], Ye]) + 
  g2^2*MatMul[Adj[hb], hb]*(36*trace[Adj[Yb], Yb] + 12*trace[Adj[Ye], Ye]) + 
  g2^2*MatMul[Adj[Yb], mb, Yb]*(36*trace[Adj[Yb], Yb] + 
    12*trace[Adj[Ye], Ye]) + g2^2*M2^2*MatMul[Adj[Yb], Yb]*
   (72*trace[Adj[Yb], Yb] + 24*trace[Adj[Ye], Ye]) + 
  g3^4*M3^2*(-320*trace[Adj[Yb], Yb] - 320*trace[Adj[Yt], Yt]) + 
  g2^4*M2^2*(-270*trace[Adj[Yb], Yb] - 90*trace[Adj[Ye], Ye] - 
    270*trace[Adj[Yt], Yt]) + g1^4*M1^2*((-14*trace[Adj[Yb], Yb])/5 - 
    (18*trace[Adj[Ye], Ye])/5 - (26*trace[Adj[Yt], Yt])/5) + 
  36*g2^2*MatMul[Adj[ht], ht]*trace[Adj[Yt], Yt] - 
  36*g2^2*M2*MatMul[Adj[ht], Yt]*trace[Adj[Yt], Yt] - 
  36*g2^2*M2*MatMul[Adj[Yt], ht]*trace[Adj[Yt], Yt] + 
  72*g2^2*M2^2*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  18*g2^2*MatMul[mq, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  36*g2^2*MatMul[Adj[Yt], mt, Yt]*trace[Adj[Yt], Yt] + 
  18*g2^2*MatMul[Adj[Yt], Yt, mq]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[ht], ht, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[ht], Yt, Adj[Yt], ht]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[Yt], ht, Adj[ht], Yt]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[Yt], Yt, Adj[ht], ht]*trace[Adj[Yt], Yt] + 
  6*MatMul[mq, Adj[Yt], Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[Yt], mt, Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[Yt], Yt, mq, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[Yt], Yt, Adj[Yt], mt, Yt]*trace[Adj[Yt], Yt] + 
  6*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, mq]*trace[Adj[Yt], Yt] + 
  MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*(12*trace[hb, Adj[hb]] + 
    4*trace[he, Adj[he]] + 36*mh1*trace[Adj[Yb], Yb] + 
    12*mh1*trace[Adj[Ye], Ye] + 12*trace[Yb, Adj[Yb], mb] + 
    4*trace[Ye, Adj[Ye], me] + 12*trace[Adj[Yb], Yb, mq] + 
    4*trace[Adj[Ye], Ye, ml]) + g2^2*MatMul[Adj[Yb], Yb]*
   (36*trace[hb, Adj[hb]] + 12*trace[he, Adj[he]] + 
    72*mh1*trace[Adj[Yb], Yb] + 24*mh1*trace[Adj[Ye], Ye] + 
    36*trace[Yb, Adj[Yb], mb] + 12*trace[Ye, Adj[Ye], me] + 
    36*trace[Adj[Yb], Yb, mq] + 12*trace[Adj[Ye], Ye, ml]) + 
  g3^4*(-64*trace[hb, Adj[hb]] - 64*trace[ht, Adj[ht]] - 
    64*mh1*trace[Adj[Yb], Yb] - 64*mh2*trace[Adj[Yt], Yt] - 
    64*trace[Yb, Adj[Yb], mb] - 64*trace[Yt, Adj[Yt], mt] - 
    64*trace[Adj[Yb], Yb, mq] - 64*trace[Adj[Yt], Yt, mq]) + 
  g2^4*(-54*trace[hb, Adj[hb]] - 18*trace[he, Adj[he]] - 
    54*trace[ht, Adj[ht]] - 54*mh1*trace[Adj[Yb], Yb] - 
    18*mh1*trace[Adj[Ye], Ye] - 54*mh2*trace[Adj[Yt], Yt] - 
    54*trace[Yb, Adj[Yb], mb] - 18*trace[Ye, Adj[Ye], me] - 
    54*trace[Yt, Adj[Yt], mt] - 54*trace[Adj[Yb], Yb, mq] - 
    18*trace[Adj[Ye], Ye, ml] - 54*trace[Adj[Yt], Yt, mq]) + 
  g1^4*((-14*trace[hb, Adj[hb]])/25 - (18*trace[he, Adj[he]])/25 - 
    (26*trace[ht, Adj[ht]])/25 - (14*mh1*trace[Adj[Yb], Yb])/25 - 
    (18*mh1*trace[Adj[Ye], Ye])/25 - (26*mh2*trace[Adj[Yt], Yt])/25 - 
    (14*trace[Yb, Adj[Yb], mb])/25 - (18*trace[Ye, Adj[Ye], me])/25 - 
    (26*trace[Yt, Adj[Yt], mt])/25 - (14*trace[Adj[Yb], Yb, mq])/25 - 
    (18*trace[Adj[Ye], Ye, ml])/25 - (26*trace[Adj[Yt], Yt, mq])/25) + 
  MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*(12*trace[ht, Adj[ht]] + 
    36*mh2*trace[Adj[Yt], Yt] + 12*trace[Yt, Adj[Yt], mt] + 
    12*trace[Adj[Yt], Yt, mq]) + g2^2*MatMul[Adj[Yt], Yt]*
   (36*trace[ht, Adj[ht]] + 72*mh2*trace[Adj[Yt], Yt] + 
    36*trace[Yt, Adj[Yt], mt] + 36*trace[Adj[Yt], Yt, mq]) + 
  MatMul[Adj[hb], Yb]*(-36*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
    12*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb] - 12*trace[hb, Adj[Yb]]*
     trace[Adj[Ye], Ye] - 4*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] + 
    12*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 72*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
    24*trace[Adj[Ye], he, Adj[Ye], Ye] + 
    12*trace[Adj[Yt], ht, Adj[Yb], Yb]) + MatMul[Adj[ht], Yt]*
   (-36*trace[ht, Adj[Yt]]*trace[Adj[Yt], Yt] + 
    12*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 12*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
    72*trace[Adj[Yt], ht, Adj[Yt], Yt]) + MatMul[Adj[Yb], hb]*
   (-36*trace[Adj[hb], Yb]*trace[Adj[Yb], Yb] - 12*trace[Adj[he], Ye]*
     trace[Adj[Yb], Yb] - 12*trace[Adj[hb], Yb]*trace[Adj[Ye], Ye] - 
    4*trace[Adj[he], Ye]*trace[Adj[Ye], Ye] + 
    12*trace[Adj[ht], Yt, Adj[Yb], Yb] + 72*trace[Adj[Yb], Yb, Adj[hb], Yb] + 
    24*trace[Adj[Ye], Ye, Adj[he], Ye] + 
    12*trace[Adj[Yt], Yt, Adj[hb], Yb]) + MatMul[Adj[Yt], ht]*
   (-36*trace[Adj[ht], Yt]*trace[Adj[Yt], Yt] + 
    12*trace[Adj[ht], Yt, Adj[Yb], Yb] + 12*trace[Adj[Yt], Yt, Adj[hb], Yb] + 
    72*trace[Adj[Yt], Yt, Adj[ht], Yt]) + MatMul[mq, Adj[Yb], Yb]*
   (-9*trace[Adj[Yb], Yb]^2 - 6*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
    trace[Adj[Ye], Ye]^2 + 18*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    6*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 6*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + 
  MatMul[Adj[Yb], Yb, mq]*(-9*trace[Adj[Yb], Yb]^2 - 
    6*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - trace[Adj[Ye], Ye]^2 + 
    18*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 6*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    6*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + MatMul[Adj[hb], hb]*
   (-18*trace[Adj[Yb], Yb]^2 - 12*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
    2*trace[Adj[Ye], Ye]^2 + 36*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    12*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    12*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + MatMul[Adj[Yb], mb, Yb]*
   (-18*trace[Adj[Yb], Yb]^2 - 12*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
    2*trace[Adj[Ye], Ye]^2 + 36*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    12*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    12*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + MatMul[mq, Adj[Yt], Yt]*
   (-9*trace[Adj[Yt], Yt]^2 + 6*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    18*trace[Adj[Yt], Yt, Adj[Yt], Yt]) + MatMul[Adj[Yt], Yt, mq]*
   (-9*trace[Adj[Yt], Yt]^2 + 6*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    18*trace[Adj[Yt], Yt, Adj[Yt], Yt]) + MatMul[Adj[ht], ht]*
   (-18*trace[Adj[Yt], Yt]^2 + 12*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    36*trace[Adj[Yt], Yt, Adj[Yt], Yt]) + MatMul[Adj[Yt], mt, Yt]*
   (-18*trace[Adj[Yt], Yt]^2 + 12*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    36*trace[Adj[Yt], Yt, Adj[Yt], Yt]) + MatMul[Adj[Yb], Yb]*
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
  MatMul[Adj[Yt], Yt]*(-36*trace[ht, Adj[Yt]]*trace[Adj[ht], Yt] - 
    36*trace[ht, Adj[ht]]*trace[Adj[Yt], Yt] - 54*mh2*trace[Adj[Yt], Yt]^2 - 
    36*trace[Adj[Yt], Yt]*trace[Yt, Adj[Yt], mt] - 
    36*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, mq] + 
    12*trace[hb, Adj[ht], Yt, Adj[Yb]] + 12*trace[Adj[hb], hb, Adj[Yt], Yt] + 
    12*trace[Adj[ht], ht, Adj[Yb], Yb] + 72*trace[Adj[ht], ht, Adj[Yt], Yt] + 
    12*trace[Adj[Yt], ht, Adj[hb], Yb] + 72*trace[Adj[Yt], ht, Adj[ht], Yt] + 
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
   (2464/15 - (1056*Zeta[3])/5) + g3^4*MatMul[Adj[hb], hb]*
   (16/3 - (544*Zeta[3])/3) + g3^4*MatMul[Adj[ht], ht]*
   (16/3 - (544*Zeta[3])/3) + g3^4*MatMul[Adj[Yb], mb, Yb]*
   (16/3 - (544*Zeta[3])/3) + g3^4*MatMul[Adj[Yt], mt, Yt]*
   (16/3 - (544*Zeta[3])/3) + g1^2*g3^4*M1*M3*(4864/45 - (704*Zeta[3])/5) + 
  g1^2*g2^4*M2^2*(379/5 - (486*Zeta[3])/5) + g3^4*MatMul[mq, Adj[Yb], Yb]*
   (8/3 - (272*Zeta[3])/3) + g3^4*MatMul[mq, Adj[Yt], Yt]*
   (8/3 - (272*Zeta[3])/3) + g3^4*MatMul[Adj[Yb], Yb, mq]*
   (8/3 - (272*Zeta[3])/3) + g3^4*MatMul[Adj[Yt], Yt, mq]*
   (8/3 - (272*Zeta[3])/3) + g1^2*g3^4*M1^2*(592/9 - (352*Zeta[3])/5) + 
  g1^2*g2^4*M1*M2*(50 - (324*Zeta[3])/5) + g2^4*MatMul[Adj[hb], hb]*
   (-45/2 - 63*Zeta[3]) + g2^4*MatMul[Adj[ht], ht]*(-45/2 - 63*Zeta[3]) + 
  g2^4*MatMul[Adj[Yb], mb, Yb]*(-45/2 - 63*Zeta[3]) + 
  g2^4*MatMul[Adj[Yt], mt, Yt]*(-45/2 - 63*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[hb], Yb, Adj[Yb], Yb]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[ht], Yt, Adj[Yt], Yt]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Yb], hb, Adj[Yb], Yb]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Yb], Yb, Adj[hb], Yb]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Yb], Yb, Adj[Yb], hb]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Yt], ht, Adj[Yt], Yt]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Yt], Yt, Adj[ht], Yt]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Adj[Yt], Yt, Adj[Yt], ht]*(6 - 36*Zeta[3]) + 
  g1^2*g2^4*M1^2*(152/5 - (162*Zeta[3])/5) + g2^4*MatMul[mq, Adj[Yb], Yb]*
   (-45/4 - (63*Zeta[3])/2) + g2^4*MatMul[mq, Adj[Yt], Yt]*
   (-45/4 - (63*Zeta[3])/2) + g2^4*MatMul[Adj[Yb], Yb, mq]*
   (-45/4 - (63*Zeta[3])/2) + g2^4*MatMul[Adj[Yt], Yt, mq]*
   (-45/4 - (63*Zeta[3])/2) + g1^2*M1^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*
   (44/3 - 24*Zeta[3]) + g1^4*g3^2*M1^2*(776/75 - (528*Zeta[3])/25) + 
  g1^6*M1^2*(57511/1125 - (2388*Zeta[3])/125) + 
  g1^2*g2^2*M1*MatMul[Adj[ht], Yt]*(59/5 - 18*Zeta[3]) + 
  g1^2*g2^2*M2*MatMul[Adj[ht], Yt]*(59/5 - 18*Zeta[3]) + 
  g1^2*g2^2*M1*MatMul[Adj[Yt], ht]*(59/5 - 18*Zeta[3]) + 
  g1^2*g2^2*M2*MatMul[Adj[Yt], ht]*(59/5 - 18*Zeta[3]) + 
  g1^2*g3^2*M1*MatMul[Adj[hb], Yb]*(152/15 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[Adj[hb], Yb]*(152/15 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M1*MatMul[Adj[Yb], hb]*(152/15 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[Adj[Yb], hb]*(152/15 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M1*MatMul[Adj[ht], Yt]*(136/5 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[Adj[ht], Yt]*(136/5 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M1*MatMul[Adj[Yt], ht]*(136/5 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[Adj[Yt], ht]*(136/5 - (256*Zeta[3])/15) + 
  g1^4*g3^2*M1*M3*(1552/225 - (352*Zeta[3])/25) + 
  g1^2*MatMul[Adj[ht], ht, Adj[Yt], Yt]*(22/3 - 12*Zeta[3]) + 
  g1^2*MatMul[Adj[ht], Yt, Adj[Yt], ht]*(22/3 - 12*Zeta[3]) + 
  g1^2*MatMul[Adj[Yt], ht, Adj[ht], Yt]*(22/3 - 12*Zeta[3]) + 
  g1^2*MatMul[Adj[Yt], Yt, Adj[ht], ht]*(22/3 - 12*Zeta[3]) + 
  g1^2*MatMul[Adj[Yt], mt, Yt, Adj[Yt], Yt]*(22/3 - 12*Zeta[3]) + 
  g1^2*MatMul[Adj[Yt], Yt, mq, Adj[Yt], Yt]*(22/3 - 12*Zeta[3]) + 
  g1^2*MatMul[Adj[Yt], Yt, Adj[Yt], mt, Yt]*(22/3 - 12*Zeta[3]) + 
  g1^4*g3^2*M3^2*(208/45 - (176*Zeta[3])/25) + 
  g1^4*g2^2*M1^2*(33/25 - (162*Zeta[3])/25) + 
  g1^2*MatMul[mq, Adj[Yt], Yt, Adj[Yt], Yt]*(11/3 - 6*Zeta[3]) + 
  g1^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, mq]*(11/3 - 6*Zeta[3]) + 
  g1^2*M1^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*(28/15 - (24*Zeta[3])/5) + 
  g1^4*g2^2*M1*M2*(22/25 - (108*Zeta[3])/25) + g1^4*M1*MatMul[Adj[ht], Yt]*
   (3767/75 - (286*Zeta[3])/75) + g1^4*M1*MatMul[Adj[Yt], ht]*
   (3767/75 - (286*Zeta[3])/75) + g1^2*g2^2*M1*MatMul[Adj[hb], Yb]*
   (41/5 - (18*Zeta[3])/5) + g1^2*g2^2*M2*MatMul[Adj[hb], Yb]*
   (41/5 - (18*Zeta[3])/5) + g1^2*g2^2*M1*MatMul[Adj[Yb], hb]*
   (41/5 - (18*Zeta[3])/5) + g1^2*g2^2*M2*MatMul[Adj[Yb], hb]*
   (41/5 - (18*Zeta[3])/5) + g1^2*MatMul[Adj[hb], hb, Adj[Yb], Yb]*
   (14/15 - (12*Zeta[3])/5) + g1^2*MatMul[Adj[hb], Yb, Adj[Yb], hb]*
   (14/15 - (12*Zeta[3])/5) + g1^2*MatMul[Adj[Yb], hb, Adj[hb], Yb]*
   (14/15 - (12*Zeta[3])/5) + g1^2*MatMul[Adj[Yb], Yb, Adj[hb], hb]*
   (14/15 - (12*Zeta[3])/5) + g1^2*MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb]*
   (14/15 - (12*Zeta[3])/5) + g1^2*MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb]*
   (14/15 - (12*Zeta[3])/5) + g1^2*MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb]*
   (14/15 - (12*Zeta[3])/5) + g1^4*g2^2*M2^2*(4/5 - (54*Zeta[3])/25) + 
  g1^2*MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb]*(7/15 - (6*Zeta[3])/5) + 
  g1^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq]*(7/15 - (6*Zeta[3])/5) + 
  g1^4*M1*MatMul[Adj[hb], Yb]*(633/25 - (14*Zeta[3])/15) + 
  g1^4*M1*MatMul[Adj[Yb], hb]*(633/25 - (14*Zeta[3])/15) + 
  g1^4*MatMul[mq, Adj[Yb], Yb]*(-633/100 + (7*Zeta[3])/30) + 
  g1^4*MatMul[Adj[Yb], Yb, mq]*(-633/100 + (7*Zeta[3])/30) + 
  g1^4*MatMul[Adj[hb], hb]*(-633/50 + (7*Zeta[3])/15) + 
  g1^4*MatMul[Adj[Yb], mb, Yb]*(-633/50 + (7*Zeta[3])/15) + 
  g1^4*MatMul[mq, Adj[Yt], Yt]*(-3767/300 + (143*Zeta[3])/150) + 
  g1^4*MatMul[Adj[Yt], Yt, mq]*(-3767/300 + (143*Zeta[3])/150) + 
  12*MatMul[Adj[hb], hb, Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
  12*MatMul[Adj[hb], Yb, Adj[Yb], hb, Adj[Yb], Yb]*Zeta[3] + 
  12*MatMul[Adj[hb], Yb, Adj[Yb], Yb, Adj[Yb], hb]*Zeta[3] + 
  12*MatMul[Adj[ht], ht, Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] + 
  12*MatMul[Adj[ht], Yt, Adj[Yt], ht, Adj[Yt], Yt]*Zeta[3] + 
  12*MatMul[Adj[ht], Yt, Adj[Yt], Yt, Adj[Yt], ht]*Zeta[3] + 
  12*MatMul[Adj[Yb], hb, Adj[hb], Yb, Adj[Yb], Yb]*Zeta[3] + 
  12*MatMul[Adj[Yb], hb, Adj[Yb], Yb, Adj[hb], Yb]*Zeta[3] + 
  12*MatMul[Adj[Yb], Yb, Adj[hb], hb, Adj[Yb], Yb]*Zeta[3] + 
  12*MatMul[Adj[Yb], Yb, Adj[hb], Yb, Adj[Yb], hb]*Zeta[3] + 
  12*MatMul[Adj[Yb], Yb, Adj[Yb], hb, Adj[hb], Yb]*Zeta[3] + 
  12*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, Adj[hb], hb]*Zeta[3] + 
  36*mh1*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
  12*MatMul[Adj[Yt], ht, Adj[ht], Yt, Adj[Yt], Yt]*Zeta[3] + 
  12*MatMul[Adj[Yt], ht, Adj[Yt], Yt, Adj[ht], Yt]*Zeta[3] + 
  12*MatMul[Adj[Yt], Yt, Adj[ht], ht, Adj[Yt], Yt]*Zeta[3] + 
  12*MatMul[Adj[Yt], Yt, Adj[ht], Yt, Adj[Yt], ht]*Zeta[3] + 
  12*MatMul[Adj[Yt], Yt, Adj[Yt], ht, Adj[ht], Yt]*Zeta[3] + 
  12*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, Adj[ht], ht]*Zeta[3] + 
  36*mh2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] + 
  6*MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
  6*MatMul[mq, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] + 
  12*MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
  12*MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
  12*MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb, Adj[Yb], Yb]*Zeta[3] + 
  12*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq, Adj[Yb], Yb]*Zeta[3] + 
  12*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], mb, Yb]*Zeta[3] + 
  6*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb, mq]*Zeta[3] + 
  12*MatMul[Adj[Yt], mt, Yt, Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] + 
  12*MatMul[Adj[Yt], Yt, mq, Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] + 
  12*MatMul[Adj[Yt], Yt, Adj[Yt], mt, Yt, Adj[Yt], Yt]*Zeta[3] + 
  12*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, mq, Adj[Yt], Yt]*Zeta[3] + 
  12*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], mt, Yt]*Zeta[3] + 
  6*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt, mq]*Zeta[3] + 
  g1^2*g2^2*MatMul[mq, Adj[Yb], Yb]*(-41/10 + (9*Zeta[3])/5) + 
  g1^2*g2^2*MatMul[Adj[Yb], Yb, mq]*(-41/10 + (9*Zeta[3])/5) + 
  g1^4*MatMul[Adj[ht], ht]*(-3767/150 + (143*Zeta[3])/75) + 
  g1^4*MatMul[Adj[Yt], mt, Yt]*(-3767/150 + (143*Zeta[3])/75) + 
  g1^2*M1*MatMul[Adj[hb], Yb, Adj[Yb], Yb]*(-14/15 + (12*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Yb], hb, Adj[Yb], Yb]*(-14/15 + (12*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Yb], Yb, Adj[hb], Yb]*(-14/15 + (12*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Yb], Yb, Adj[Yb], hb]*(-14/15 + (12*Zeta[3])/5) + 
  g1^4*M1^2*MatMul[Adj[Yb], Yb]*(-1899/25 + (14*Zeta[3])/5) + 
  g1^2*g2^2*MatMul[Adj[hb], hb]*(-41/5 + (18*Zeta[3])/5) + 
  g1^2*g2^2*MatMul[Adj[Yb], mb, Yb]*(-41/5 + (18*Zeta[3])/5) + 
  g1^2*g2^2*M1^2*MatMul[Adj[Yb], Yb]*(-82/5 + (36*Zeta[3])/5) + 
  g1^2*g2^2*M1*M2*MatMul[Adj[Yb], Yb]*(-82/5 + (36*Zeta[3])/5) + 
  g1^2*g2^2*M2^2*MatMul[Adj[Yb], Yb]*(-82/5 + (36*Zeta[3])/5) + 
  g1^2*g3^2*MatMul[mq, Adj[Yt], Yt]*(-68/5 + (128*Zeta[3])/15) + 
  g1^2*g3^2*MatMul[Adj[Yt], Yt, mq]*(-68/5 + (128*Zeta[3])/15) + 
  g1^2*g3^2*MatMul[mq, Adj[Yb], Yb]*(-76/15 + (128*Zeta[3])/15) + 
  g1^2*g3^2*MatMul[Adj[Yb], Yb, mq]*(-76/15 + (128*Zeta[3])/15) + 
  g1^2*g2^2*MatMul[mq, Adj[Yt], Yt]*(-59/10 + 9*Zeta[3]) + 
  g1^2*g2^2*MatMul[Adj[Yt], Yt, mq]*(-59/10 + 9*Zeta[3]) + 
  g1^4*M1^2*MatMul[Adj[Yt], Yt]*(-3767/25 + (286*Zeta[3])/25) + 
  g1^2*M1*MatMul[Adj[ht], Yt, Adj[Yt], Yt]*(-22/3 + 12*Zeta[3]) + 
  g1^2*M1*MatMul[Adj[Yt], ht, Adj[Yt], Yt]*(-22/3 + 12*Zeta[3]) + 
  g1^2*M1*MatMul[Adj[Yt], Yt, Adj[ht], Yt]*(-22/3 + 12*Zeta[3]) + 
  g1^2*M1*MatMul[Adj[Yt], Yt, Adj[Yt], ht]*(-22/3 + 12*Zeta[3]) + 
  g1^2*g3^2*MatMul[Adj[ht], ht]*(-136/5 + (256*Zeta[3])/15) + 
  g1^2*g3^2*MatMul[Adj[Yt], mt, Yt]*(-136/5 + (256*Zeta[3])/15) + 
  g1^2*g3^2*MatMul[Adj[hb], hb]*(-152/15 + (256*Zeta[3])/15) + 
  g1^2*g3^2*MatMul[Adj[Yb], mb, Yb]*(-152/15 + (256*Zeta[3])/15) + 
  g1^2*g2^2*MatMul[Adj[ht], ht]*(-59/5 + 18*Zeta[3]) + 
  g1^2*g2^2*MatMul[Adj[Yt], mt, Yt]*(-59/5 + 18*Zeta[3]) + 
  g2^2*MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb]*(-3 + 18*Zeta[3]) + 
  g2^2*MatMul[mq, Adj[Yt], Yt, Adj[Yt], Yt]*(-3 + 18*Zeta[3]) + 
  g2^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq]*(-3 + 18*Zeta[3]) + 
  g2^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, mq]*(-3 + 18*Zeta[3]) + 
  g1^2*g3^2*M1^2*MatMul[Adj[Yt], Yt]*(-272/5 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M1*M3*MatMul[Adj[Yt], Yt]*(-272/5 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M3^2*MatMul[Adj[Yt], Yt]*(-272/5 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M1^2*MatMul[Adj[Yb], Yb]*(-304/15 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M1*M3*MatMul[Adj[Yb], Yb]*(-304/15 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M3^2*MatMul[Adj[Yb], Yb]*(-304/15 + (512*Zeta[3])/15) + 
  g1^2*g2^2*M1^2*MatMul[Adj[Yt], Yt]*(-118/5 + 36*Zeta[3]) + 
  g1^2*g2^2*M1*M2*MatMul[Adj[Yt], Yt]*(-118/5 + 36*Zeta[3]) + 
  g1^2*g2^2*M2^2*MatMul[Adj[Yt], Yt]*(-118/5 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[hb], hb, Adj[Yb], Yb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[hb], Yb, Adj[Yb], hb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[ht], ht, Adj[Yt], Yt]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[ht], Yt, Adj[Yt], ht]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Yb], hb, Adj[hb], Yb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Yb], Yb, Adj[hb], hb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Yt], ht, Adj[ht], Yt]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Yt], Yt, Adj[ht], ht]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Yt], mt, Yt, Adj[Yt], Yt]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Yt], Yt, mq, Adj[Yt], Yt]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Adj[Yt], Yt, Adj[Yt], mt, Yt]*(-6 + 36*Zeta[3]) + 
  g2^2*M2^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*(-12 + 72*Zeta[3]) + 
  g2^2*M2^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*(-12 + 72*Zeta[3]) + 
  g2^4*M2*MatMul[Adj[hb], Yb]*(45 + 126*Zeta[3]) + 
  g2^4*M2*MatMul[Adj[ht], Yt]*(45 + 126*Zeta[3]) + 
  g2^4*M2*MatMul[Adj[Yb], hb]*(45 + 126*Zeta[3]) + 
  g2^4*M2*MatMul[Adj[Yt], ht]*(45 + 126*Zeta[3]) + 
  g3^4*M3*MatMul[Adj[hb], Yb]*(-32/3 + (1088*Zeta[3])/3) + 
  g3^4*M3*MatMul[Adj[ht], Yt]*(-32/3 + (1088*Zeta[3])/3) + 
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
  g3^2*MatMul[Adj[hb], Yb]*(-16*trace[hb, Adj[Yb]] + 16*trace[he, Adj[Ye]] + 
    96*trace[hb, Adj[Yb]]*Zeta[3]) + g1^2*MatMul[Adj[hb], Yb]*
   ((32*trace[hb, Adj[Yb]])/5 - (16*trace[he, Adj[Ye]])/5 - 
    (48*trace[hb, Adj[Yb]]*Zeta[3])/5 + (24*trace[he, Adj[Ye]]*Zeta[3])/5) + 
  g1^2*MatMul[Adj[ht], Yt]*(4*trace[ht, Adj[Yt]] - 
    (48*trace[ht, Adj[Yt]]*Zeta[3])/5) + g3^2*MatMul[Adj[ht], Yt]*
   (-16*trace[ht, Adj[Yt]] + 96*trace[ht, Adj[Yt]]*Zeta[3]) + 
  g3^2*M3*MatMul[Adj[Yb], Yb]*(16*trace[hb, Adj[Yb]] - 
    16*trace[he, Adj[Ye]] + 16*trace[Adj[hb], Yb] - 16*trace[Adj[he], Ye] - 
    96*trace[hb, Adj[Yb]]*Zeta[3] - 96*trace[Adj[hb], Yb]*Zeta[3]) + 
  g3^2*MatMul[Adj[Yb], hb]*(-16*trace[Adj[hb], Yb] + 16*trace[Adj[he], Ye] + 
    96*trace[Adj[hb], Yb]*Zeta[3]) + g1^2*M1*MatMul[Adj[Yb], Yb]*
   ((-32*trace[hb, Adj[Yb]])/5 + (16*trace[he, Adj[Ye]])/5 - 
    (32*trace[Adj[hb], Yb])/5 + (16*trace[Adj[he], Ye])/5 + 
    (48*trace[hb, Adj[Yb]]*Zeta[3])/5 - (24*trace[he, Adj[Ye]]*Zeta[3])/5 + 
    (48*trace[Adj[hb], Yb]*Zeta[3])/5 - (24*trace[Adj[he], Ye]*Zeta[3])/5) + 
  g1^2*MatMul[Adj[Yb], hb]*((32*trace[Adj[hb], Yb])/5 - 
    (16*trace[Adj[he], Ye])/5 - (48*trace[Adj[hb], Yb]*Zeta[3])/5 + 
    (24*trace[Adj[he], Ye]*Zeta[3])/5) + g3^2*M3*MatMul[Adj[Yt], Yt]*
   (16*trace[ht, Adj[Yt]] + 16*trace[Adj[ht], Yt] - 
    96*trace[ht, Adj[Yt]]*Zeta[3] - 96*trace[Adj[ht], Yt]*Zeta[3]) + 
  g1^2*MatMul[Adj[Yt], ht]*(4*trace[Adj[ht], Yt] - 
    (48*trace[Adj[ht], Yt]*Zeta[3])/5) + g1^2*M1*MatMul[Adj[Yt], Yt]*
   (-4*trace[ht, Adj[Yt]] - 4*trace[Adj[ht], Yt] + 
    (48*trace[ht, Adj[Yt]]*Zeta[3])/5 + (48*trace[Adj[ht], Yt]*Zeta[3])/5) + 
  g3^2*MatMul[Adj[Yt], ht]*(-16*trace[Adj[ht], Yt] + 
    96*trace[Adj[ht], Yt]*Zeta[3]) + g3^2*M3*MatMul[Adj[hb], Yb]*
   (16*trace[Adj[Yb], Yb] - 16*trace[Adj[Ye], Ye] - 
    96*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*M3*MatMul[Adj[Yb], hb]*
   (16*trace[Adj[Yb], Yb] - 16*trace[Adj[Ye], Ye] - 
    96*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*MatMul[mq, Adj[Yb], Yb]*
   (-8*trace[Adj[Yb], Yb] + 8*trace[Adj[Ye], Ye] + 
    48*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*MatMul[Adj[Yb], Yb, mq]*
   (-8*trace[Adj[Yb], Yb] + 8*trace[Adj[Ye], Ye] + 
    48*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*MatMul[Adj[hb], hb]*
   (-16*trace[Adj[Yb], Yb] + 16*trace[Adj[Ye], Ye] + 
    96*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*MatMul[Adj[Yb], mb, Yb]*
   (-16*trace[Adj[Yb], Yb] + 16*trace[Adj[Ye], Ye] + 
    96*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*M3^2*MatMul[Adj[Yb], Yb]*
   (-32*trace[Adj[Yb], Yb] + 32*trace[Adj[Ye], Ye] + 
    192*trace[Adj[Yb], Yb]*Zeta[3]) + g1^2*M1*MatMul[Adj[hb], Yb]*
   ((-32*trace[Adj[Yb], Yb])/5 + (16*trace[Adj[Ye], Ye])/5 + 
    (48*trace[Adj[Yb], Yb]*Zeta[3])/5 - (24*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Yb], hb]*((-32*trace[Adj[Yb], Yb])/5 + 
    (16*trace[Adj[Ye], Ye])/5 + (48*trace[Adj[Yb], Yb]*Zeta[3])/5 - 
    (24*trace[Adj[Ye], Ye]*Zeta[3])/5) + g1^2*MatMul[mq, Adj[Yb], Yb]*
   ((16*trace[Adj[Yb], Yb])/5 - (8*trace[Adj[Ye], Ye])/5 - 
    (24*trace[Adj[Yb], Yb]*Zeta[3])/5 + (12*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*MatMul[Adj[Yb], Yb, mq]*((16*trace[Adj[Yb], Yb])/5 - 
    (8*trace[Adj[Ye], Ye])/5 - (24*trace[Adj[Yb], Yb]*Zeta[3])/5 + 
    (12*trace[Adj[Ye], Ye]*Zeta[3])/5) + g1^2*MatMul[Adj[hb], hb]*
   ((32*trace[Adj[Yb], Yb])/5 - (16*trace[Adj[Ye], Ye])/5 - 
    (48*trace[Adj[Yb], Yb]*Zeta[3])/5 + (24*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*MatMul[Adj[Yb], mb, Yb]*((32*trace[Adj[Yb], Yb])/5 - 
    (16*trace[Adj[Ye], Ye])/5 - (48*trace[Adj[Yb], Yb]*Zeta[3])/5 + 
    (24*trace[Adj[Ye], Ye]*Zeta[3])/5) + g1^2*M1^2*MatMul[Adj[Yb], Yb]*
   ((64*trace[Adj[Yb], Yb])/5 - (32*trace[Adj[Ye], Ye])/5 - 
    (96*trace[Adj[Yb], Yb]*Zeta[3])/5 + (48*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g3^2*M3*MatMul[Adj[ht], Yt]*(16*trace[Adj[Yt], Yt] - 
    96*trace[Adj[Yt], Yt]*Zeta[3]) + g3^2*M3*MatMul[Adj[Yt], ht]*
   (16*trace[Adj[Yt], Yt] - 96*trace[Adj[Yt], Yt]*Zeta[3]) + 
  g1^2*M1^2*MatMul[Adj[Yt], Yt]*(8*trace[Adj[Yt], Yt] - 
    (96*trace[Adj[Yt], Yt]*Zeta[3])/5) + g1^2*MatMul[Adj[ht], ht]*
   (4*trace[Adj[Yt], Yt] - (48*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  g1^2*MatMul[Adj[Yt], mt, Yt]*(4*trace[Adj[Yt], Yt] - 
    (48*trace[Adj[Yt], Yt]*Zeta[3])/5) + g1^2*MatMul[mq, Adj[Yt], Yt]*
   (2*trace[Adj[Yt], Yt] - (24*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  g1^2*MatMul[Adj[Yt], Yt, mq]*(2*trace[Adj[Yt], Yt] - 
    (24*trace[Adj[Yt], Yt]*Zeta[3])/5) + g1^2*M1*MatMul[Adj[ht], Yt]*
   (-4*trace[Adj[Yt], Yt] + (48*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  g1^2*M1*MatMul[Adj[Yt], ht]*(-4*trace[Adj[Yt], Yt] + 
    (48*trace[Adj[Yt], Yt]*Zeta[3])/5) + g3^2*MatMul[mq, Adj[Yt], Yt]*
   (-8*trace[Adj[Yt], Yt] + 48*trace[Adj[Yt], Yt]*Zeta[3]) + 
  g3^2*MatMul[Adj[Yt], Yt, mq]*(-8*trace[Adj[Yt], Yt] + 
    48*trace[Adj[Yt], Yt]*Zeta[3]) + g3^2*MatMul[Adj[ht], ht]*
   (-16*trace[Adj[Yt], Yt] + 96*trace[Adj[Yt], Yt]*Zeta[3]) + 
  g3^2*MatMul[Adj[Yt], mt, Yt]*(-16*trace[Adj[Yt], Yt] + 
    96*trace[Adj[Yt], Yt]*Zeta[3]) + g3^2*M3^2*MatMul[Adj[Yt], Yt]*
   (-32*trace[Adj[Yt], Yt] + 192*trace[Adj[Yt], Yt]*Zeta[3]) + 
  g3^2*MatMul[Adj[Yb], Yb]*(-16*trace[hb, Adj[hb]] + 16*trace[he, Adj[he]] - 
    32*mh1*trace[Adj[Yb], Yb] + 32*mh1*trace[Adj[Ye], Ye] - 
    16*trace[Yb, Adj[Yb], mb] + 16*trace[Ye, Adj[Ye], me] - 
    16*trace[Adj[Yb], Yb, mq] + 16*trace[Adj[Ye], Ye, ml] + 
    96*trace[hb, Adj[hb]]*Zeta[3] + 192*mh1*trace[Adj[Yb], Yb]*Zeta[3] + 
    96*trace[Yb, Adj[Yb], mb]*Zeta[3] + 96*trace[Adj[Yb], Yb, mq]*Zeta[3]) + 
  g1^2*MatMul[Adj[Yb], Yb]*((32*trace[hb, Adj[hb]])/5 - 
    (16*trace[he, Adj[he]])/5 + (64*mh1*trace[Adj[Yb], Yb])/5 - 
    (32*mh1*trace[Adj[Ye], Ye])/5 + (32*trace[Yb, Adj[Yb], mb])/5 - 
    (16*trace[Ye, Adj[Ye], me])/5 + (32*trace[Adj[Yb], Yb, mq])/5 - 
    (16*trace[Adj[Ye], Ye, ml])/5 - (48*trace[hb, Adj[hb]]*Zeta[3])/5 + 
    (24*trace[he, Adj[he]]*Zeta[3])/5 - (96*mh1*trace[Adj[Yb], Yb]*Zeta[3])/
     5 + (48*mh1*trace[Adj[Ye], Ye]*Zeta[3])/5 - 
    (48*trace[Yb, Adj[Yb], mb]*Zeta[3])/5 + 
    (24*trace[Ye, Adj[Ye], me]*Zeta[3])/5 - 
    (48*trace[Adj[Yb], Yb, mq]*Zeta[3])/5 + 
    (24*trace[Adj[Ye], Ye, ml]*Zeta[3])/5) + g1^2*MatMul[Adj[Yt], Yt]*
   (4*trace[ht, Adj[ht]] + 8*mh2*trace[Adj[Yt], Yt] + 
    4*trace[Yt, Adj[Yt], mt] + 4*trace[Adj[Yt], Yt, mq] - 
    (48*trace[ht, Adj[ht]]*Zeta[3])/5 - (96*mh2*trace[Adj[Yt], Yt]*Zeta[3])/
     5 - (48*trace[Yt, Adj[Yt], mt]*Zeta[3])/5 - 
    (48*trace[Adj[Yt], Yt, mq]*Zeta[3])/5) + g3^2*MatMul[Adj[Yt], Yt]*
   (-16*trace[ht, Adj[ht]] - 32*mh2*trace[Adj[Yt], Yt] - 
    16*trace[Yt, Adj[Yt], mt] - 16*trace[Adj[Yt], Yt, mq] + 
    96*trace[ht, Adj[ht]]*Zeta[3] + 192*mh2*trace[Adj[Yt], Yt]*Zeta[3] + 
    96*trace[Yt, Adj[Yt], mt]*Zeta[3] + 96*trace[Adj[Yt], Yt, mq]*Zeta[3])}

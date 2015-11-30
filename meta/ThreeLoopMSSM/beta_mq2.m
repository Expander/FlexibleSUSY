{(-2*g1^2*M1^2)/15 - 6*g2^2*M2^2 - (32*g3^2*M3^2)/3 + 2*mb*Yb*Ybc + 
  2*mt*Yt*Ytc + 2*MatMul[hbc, hb] + 2*MatMul[htc, ht] + 
  2*mh1*MatMul[Ybc, Yb] + 2*mq*MatMul[Ybc, Yb] + 2*mh2*MatMul[Ytc, Yt] + 
  2*mq*MatMul[Ytc, Yt], (199*g1^4*M1^2)/75 + (2*g1^2*g2^2*M1^2)/5 + 
  (32*g1^2*g3^2*M1^2)/45 + (2*g1^2*g2^2*M1*M2)/5 + (2*g1^2*g2^2*M2^2)/5 + 
  33*g2^4*M2^2 + 32*g2^2*g3^2*M2^2 + (32*g1^2*g3^2*M1*M3)/45 + 
  32*g2^2*g3^2*M2*M3 + (32*g1^2*g3^2*M3^2)/45 + 32*g2^2*g3^2*M3^2 - 
  (128*g3^4*M3^2)/3 + (4*g1^2*mb*Yb*Ybc)/5 + (8*g1^2*mt*Yt*Ytc)/5 + 
  (4*g1^2*MatMul[hbc, hb])/5 - (4*g1^2*M1*MatMul[hbc, Yb])/5 + 
  (8*g1^2*MatMul[htc, ht])/5 - (8*g1^2*M1*MatMul[htc, Yt])/5 - 
  (4*g1^2*M1*MatMul[Ybc, hb])/5 + (8*g1^2*M1^2*MatMul[Ybc, Yb])/5 + 
  (4*g1^2*mh1*MatMul[Ybc, Yb])/5 + (4*g1^2*mq*MatMul[Ybc, Yb])/5 - 
  4*mq*MatMul[Ybc, Yb]^2 - (8*g1^2*M1*MatMul[Ytc, ht])/5 + 
  (16*g1^2*M1^2*MatMul[Ytc, Yt])/5 + (8*g1^2*mh2*MatMul[Ytc, Yt])/5 + 
  (8*g1^2*mq*MatMul[Ytc, Yt])/5 - 4*mq*MatMul[Ytc, Yt]^2 - 
  4*mb*Ybc*MatMul[Yb, Ybc, Yb] - 4*mb*Yb*MatMul[Ybc, Yb, Ybc] - 
  4*mt*Ytc*MatMul[Yt, Ytc, Yt] - 4*mt*Yt*MatMul[Ytc, Yt, Ytc] - 
  4*MatMul[hbc, hb, Ybc, Yb] - 4*MatMul[hbc, Yb, Ybc, hb] - 
  4*MatMul[htc, ht, Ytc, Yt] - 4*MatMul[htc, Yt, Ytc, ht] - 
  4*MatMul[Ybc, hb, hbc, Yb] - 4*MatMul[Ybc, Yb, hbc, hb] - 
  8*mh1*MatMul[Ybc, Yb, Ybc, Yb] - 4*mq*MatMul[Ybc, Yb, Ybc, Yb] - 
  4*MatMul[Ytc, ht, htc, Yt] - 4*MatMul[Ytc, Yt, htc, ht] - 
  8*mh2*MatMul[Ytc, Yt, Ytc, Yt] - 4*mq*MatMul[Ytc, Yt, Ytc, Yt] + 
  g2^4*(3*mh1 + 3*mh2 + 3*trace[ml] + 9*trace[mq]) + 
  g1^4*(mh1/25 + mh2/25 + (2*trace[mb])/75 + (2*trace[me])/25 + 
    trace[ml]/25 + trace[mq]/75 + (8*trace[mt])/75) + 
  g3^4*((16*trace[mb])/3 + (32*trace[mq])/3 + (16*trace[mt])/3) + 
  MatMul[hbc, Yb]*(-6*trace[hb, Ybc] - 2*trace[he, Adj[Ye]]) + 
  MatMul[Ybc, hb]*(-6*trace[hbc, Yb] - 2*trace[hec, Ye]) - 
  6*MatMul[htc, Yt]*trace[ht, Ytc] - 6*MatMul[Ytc, ht]*trace[htc, Yt] - 
  6*mt*Yt*Ytc*trace[Ytc, Yt] - 6*MatMul[htc, ht]*trace[Ytc, Yt] - 
  6*mq*MatMul[Ytc, Yt]*trace[Ytc, Yt] + 
  mb*Yb*Ybc*(-6*trace[Ybc, Yb] - 2*trace[Adj[Ye], Ye]) + 
  MatMul[hbc, hb]*(-6*trace[Ybc, Yb] - 2*trace[Adj[Ye], Ye]) + 
  2*mq*MatMul[Ybc, Yb]*(-3*trace[Ybc, Yb] - trace[Adj[Ye], Ye]) + 
  MatMul[Ytc, Yt]*(-6*trace[ht, htc] - 12*mh2*trace[Ytc, Yt] - 
    6*trace[Yt, Ytc, mt] - 6*trace[Ytc, Yt, mq]) + 
  MatMul[Ybc, Yb]*(-6*trace[hb, hbc] - 2*trace[he, hec] - 
    12*mh1*trace[Ybc, Yb] - 4*mh1*trace[Adj[Ye], Ye] - 6*trace[Yb, Ybc, mb] - 
    6*trace[Ybc, Yb, mq] - 2*trace[Ye, Adj[Ye], me] - 
    2*trace[Adj[Ye], Ye, ml]), (-32*g1^2*g2^2*g3^2*M1^2)/5 - 
  (32*g1^2*g2^2*g3^2*M1*M2)/5 - (32*g1^2*g2^2*g3^2*M2^2)/5 - 
  (32*g1^2*g2^2*g3^2*M1*M3)/5 - (32*g1^2*g2^2*g3^2*M2*M3)/5 - 
  (32*g1^2*g2^2*g3^2*M3^2)/5 - 8*g2^2*g3^2*mb*Yb*Ybc - 
  8*g2^2*g3^2*mt*Yt*Ytc - 8*g2^2*g3^2*MatMul[hbc, hb] + 
  8*g2^2*g3^2*M2*MatMul[hbc, Yb] + 8*g2^2*g3^2*M3*MatMul[hbc, Yb] - 
  8*g2^2*g3^2*MatMul[htc, ht] + 8*g2^2*g3^2*M2*MatMul[htc, Yt] + 
  8*g2^2*g3^2*M3*MatMul[htc, Yt] + 8*g2^2*g3^2*M2*MatMul[Ybc, hb] + 
  8*g2^2*g3^2*M3*MatMul[Ybc, hb] - 16*g2^2*g3^2*M2^2*MatMul[Ybc, Yb] - 
  16*g2^2*g3^2*M2*M3*MatMul[Ybc, Yb] - 16*g2^2*g3^2*M3^2*MatMul[Ybc, Yb] - 
  8*g2^2*g3^2*mh1*MatMul[Ybc, Yb] - 8*g2^2*g3^2*mq*MatMul[Ybc, Yb] + 
  (128*g3^2*mq*MatMul[Ybc, Yb]^2)/3 + 8*g2^2*g3^2*M2*MatMul[Ytc, ht] + 
  8*g2^2*g3^2*M3*MatMul[Ytc, ht] - 16*g2^2*g3^2*M2^2*MatMul[Ytc, Yt] - 
  16*g2^2*g3^2*M2*M3*MatMul[Ytc, Yt] - 16*g2^2*g3^2*M3^2*MatMul[Ytc, Yt] - 
  8*g2^2*g3^2*mh2*MatMul[Ytc, Yt] - 8*g2^2*g3^2*mq*MatMul[Ytc, Yt] + 
  (128*g3^2*mq*MatMul[Ytc, Yt]^2)/3 + (128*g3^2*mb*Ybc*MatMul[Yb, Ybc, Yb])/
   3 + (128*g3^2*mb*Yb*MatMul[Ybc, Yb, Ybc])/3 + 
  8*mt*MatMul[Ybc, Yb, Ytc]*MatMul[Yt, Ybc, Yb] + 
  (128*g3^2*mt*Ytc*MatMul[Yt, Ytc, Yt])/3 + 8*mb*MatMul[Yb, Ytc, Yt]*
   MatMul[Ytc, Yt, Ybc] + (128*g3^2*mt*Yt*MatMul[Ytc, Yt, Ytc])/3 + 
  (128*g3^2*MatMul[hbc, hb, Ybc, Yb])/3 + (128*g3^2*MatMul[hbc, Yb, Ybc, hb])/
   3 - (128*g3^2*M3*MatMul[hbc, Yb, Ybc, Yb])/3 + 
  (128*g3^2*MatMul[htc, ht, Ytc, Yt])/3 + (128*g3^2*MatMul[htc, Yt, Ytc, ht])/
   3 - (128*g3^2*M3*MatMul[htc, Yt, Ytc, Yt])/3 + 
  (128*g3^2*MatMul[Ybc, hb, hbc, Yb])/3 - 
  (128*g3^2*M3*MatMul[Ybc, hb, Ybc, Yb])/3 + 
  (128*g3^2*MatMul[Ybc, Yb, hbc, hb])/3 - 
  (128*g3^2*M3*MatMul[Ybc, Yb, hbc, Yb])/3 - 
  (128*g3^2*M3*MatMul[Ybc, Yb, Ybc, hb])/3 + 
  (256*g3^2*M3^2*MatMul[Ybc, Yb, Ybc, Yb])/3 + 
  (256*g3^2*mh1*MatMul[Ybc, Yb, Ybc, Yb])/3 + 
  (128*g3^2*mq*MatMul[Ybc, Yb, Ybc, Yb])/3 + 8*mq*MatMul[Ybc, Yb]*
   MatMul[Ybc, Yb, Ytc, Yt] + 8*mq*MatMul[Ytc, Yt]*MatMul[Ybc, Yb, Ytc, Yt] + 
  (128*g3^2*MatMul[Ytc, ht, htc, Yt])/3 - 
  (128*g3^2*M3*MatMul[Ytc, ht, Ytc, Yt])/3 + 
  (128*g3^2*MatMul[Ytc, Yt, htc, ht])/3 - 
  (128*g3^2*M3*MatMul[Ytc, Yt, htc, Yt])/3 + 8*mq*MatMul[Ybc, Yb]*
   MatMul[Ytc, Yt, Ybc, Yb] + 8*mq*MatMul[Ytc, Yt]*MatMul[Ytc, Yt, Ybc, Yb] - 
  (128*g3^2*M3*MatMul[Ytc, Yt, Ytc, ht])/3 + 
  (256*g3^2*M3^2*MatMul[Ytc, Yt, Ytc, Yt])/3 + 
  (256*g3^2*mh2*MatMul[Ytc, Yt, Ytc, Yt])/3 + 
  (128*g3^2*mq*MatMul[Ytc, Yt, Ytc, Yt])/3 + 
  8*mb*Ybc*MatMul[Yb, Ytc, Yt, Ybc, Yb] + 
  8*mb*Yb*MatMul[Ybc, Yb, Ytc, Yt, Ybc] + 
  8*mt*Ytc*MatMul[Yt, Ybc, Yb, Ytc, Yt] + 
  8*mt*Yt*MatMul[Ytc, Yt, Ybc, Yb, Ytc] + 
  8*MatMul[hbc, hb, Ytc, Yt, Ybc, Yb] + 8*MatMul[hbc, Yb, Ytc, ht, Ybc, Yb] + 
  8*MatMul[hbc, Yb, Ytc, Yt, Ybc, hb] + 8*MatMul[htc, ht, Ybc, Yb, Ytc, Yt] + 
  8*MatMul[htc, Yt, Ybc, hb, Ytc, Yt] + 8*MatMul[htc, Yt, Ybc, Yb, Ytc, ht] + 
  8*MatMul[Ybc, hb, htc, Yt, Ybc, Yb] + 8*MatMul[Ybc, hb, Ytc, Yt, hbc, Yb] + 
  8*MatMul[Ybc, Yb, htc, ht, Ybc, Yb] + 8*MatMul[Ybc, Yb, htc, Yt, Ybc, hb] + 
  8*MatMul[Ybc, Yb, Ytc, ht, hbc, Yb] + 8*MatMul[Ybc, Yb, Ytc, Yt, hbc, hb] + 
  (16*mh1 + 8*mh2)*MatMul[Ybc, Yb, Ytc, Yt, Ybc, Yb] + 
  8*mq*MatMul[Ybc, Yb, Ytc, Yt, Ybc, Yb] + 
  8*MatMul[Ytc, ht, hbc, Yb, Ytc, Yt] + 8*MatMul[Ytc, ht, Ybc, Yb, htc, Yt] + 
  8*MatMul[Ytc, Yt, hbc, hb, Ytc, Yt] + 8*MatMul[Ytc, Yt, hbc, Yb, Ytc, ht] + 
  8*MatMul[Ytc, Yt, Ybc, hb, htc, Yt] + 8*MatMul[Ytc, Yt, Ybc, Yb, htc, ht] + 
  (8*mh1 + 16*mh2)*MatMul[Ytc, Yt, Ybc, Yb, Ytc, Yt] + 
  8*mq*MatMul[Ytc, Yt, Ybc, Yb, Ytc, Yt] + 
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
  MatMul[hbc, Yb, Ybc, Yb]*(12*trace[hb, Ybc] + 4*trace[he, Adj[Ye]]) + 
  MatMul[Ybc, Yb, hbc, Yb]*(12*trace[hb, Ybc] + 4*trace[he, Adj[Ye]]) + 
  g2^2*MatMul[hbc, Yb]*(36*trace[hb, Ybc] + 12*trace[he, Adj[Ye]]) + 
  g2^2*M2*MatMul[Ybc, Yb]*(-36*trace[hb, Ybc] - 36*trace[hbc, Yb] - 
    12*trace[he, Adj[Ye]] - 12*trace[hec, Ye]) + 
  MatMul[Ybc, hb, Ybc, Yb]*(12*trace[hbc, Yb] + 4*trace[hec, Ye]) + 
  MatMul[Ybc, Yb, Ybc, hb]*(12*trace[hbc, Yb] + 4*trace[hec, Ye]) + 
  g2^2*MatMul[Ybc, hb]*(36*trace[hbc, Yb] + 12*trace[hec, Ye]) + 
  36*g2^2*MatMul[htc, Yt]*trace[ht, Ytc] + 12*MatMul[htc, Yt, Ytc, Yt]*
   trace[ht, Ytc] + 12*MatMul[Ytc, Yt, htc, Yt]*trace[ht, Ytc] + 
  g2^2*M2*MatMul[Ytc, Yt]*(-36*trace[ht, Ytc] - 36*trace[htc, Yt]) + 
  36*g2^2*MatMul[Ytc, ht]*trace[htc, Yt] + 12*MatMul[Ytc, ht, Ytc, Yt]*
   trace[htc, Yt] + 12*MatMul[Ytc, Yt, Ytc, ht]*trace[htc, Yt] + 
  g1^4*M1*((14*trace[hb, Ybc])/15 + (14*trace[hbc, Yb])/15 + 
    (6*trace[he, Adj[Ye]])/5 + (6*trace[hec, Ye])/5 + 
    (26*trace[ht, Ytc])/15 + (26*trace[htc, Yt])/15) + 
  g2^4*M2*(90*trace[hb, Ybc] + 90*trace[hbc, Yb] + 30*trace[he, Adj[Ye]] + 
    30*trace[hec, Ye] + 90*trace[ht, Ytc] + 90*trace[htc, Yt]) + 
  g3^4*M3*((320*trace[hb, Ybc])/3 + (320*trace[hbc, Yb])/3 + 
    (320*trace[ht, Ytc])/3 + (320*trace[htc, Yt])/3) + 
  g3^4*M3^2*(-320*trace[Ybc, Yb] - 320*trace[Ytc, Yt]) + 
  36*g2^2*mt*Yt*Ytc*trace[Ytc, Yt] + 36*g2^2*MatMul[htc, ht]*trace[Ytc, Yt] - 
  36*g2^2*M2*MatMul[htc, Yt]*trace[Ytc, Yt] - 36*g2^2*M2*MatMul[Ytc, ht]*
   trace[Ytc, Yt] + 72*g2^2*M2^2*MatMul[Ytc, Yt]*trace[Ytc, Yt] + 
  36*g2^2*mq*MatMul[Ytc, Yt]*trace[Ytc, Yt] + 12*mq*MatMul[Ytc, Yt]^2*
   trace[Ytc, Yt] + 12*mt*Ytc*MatMul[Yt, Ytc, Yt]*trace[Ytc, Yt] + 
  12*mt*Yt*MatMul[Ytc, Yt, Ytc]*trace[Ytc, Yt] + 
  12*MatMul[htc, ht, Ytc, Yt]*trace[Ytc, Yt] + 12*MatMul[htc, Yt, Ytc, ht]*
   trace[Ytc, Yt] + 12*MatMul[Ytc, ht, htc, Yt]*trace[Ytc, Yt] + 
  12*MatMul[Ytc, Yt, htc, ht]*trace[Ytc, Yt] + 12*mq*MatMul[Ytc, Yt, Ytc, Yt]*
   trace[Ytc, Yt] + g2^4*M2^2*(-270*trace[Ybc, Yb] - 270*trace[Ytc, Yt] - 
    90*trace[Adj[Ye], Ye]) + g2^2*M2*MatMul[hbc, Yb]*
   (-36*trace[Ybc, Yb] - 12*trace[Adj[Ye], Ye]) + 
  g2^2*M2*MatMul[Ybc, hb]*(-36*trace[Ybc, Yb] - 12*trace[Adj[Ye], Ye]) + 
  g1^4*M1^2*((-14*trace[Ybc, Yb])/5 - (26*trace[Ytc, Yt])/5 - 
    (18*trace[Adj[Ye], Ye])/5) + 2*mq*MatMul[Ybc, Yb, Ybc, Yb]*
   (6*trace[Ybc, Yb] + 2*trace[Adj[Ye], Ye]) + 
  mq*MatMul[Ybc, Yb]^2*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  mb*Ybc*MatMul[Yb, Ybc, Yb]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  mb*Yb*MatMul[Ybc, Yb, Ybc]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[hbc, hb, Ybc, Yb]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[hbc, Yb, Ybc, hb]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[Ybc, hb, hbc, Yb]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[Ybc, Yb, hbc, hb]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  2*g2^2*mq*MatMul[Ybc, Yb]*(18*trace[Ybc, Yb] + 6*trace[Adj[Ye], Ye]) + 
  g2^2*mb*Yb*Ybc*(36*trace[Ybc, Yb] + 12*trace[Adj[Ye], Ye]) + 
  g2^2*MatMul[hbc, hb]*(36*trace[Ybc, Yb] + 12*trace[Adj[Ye], Ye]) + 
  g2^2*M2^2*MatMul[Ybc, Yb]*(72*trace[Ybc, Yb] + 24*trace[Adj[Ye], Ye]) + 
  g3^4*(-64*trace[hb, hbc] - 64*trace[ht, htc] - 64*mh1*trace[Ybc, Yb] - 
    64*mh2*trace[Ytc, Yt] - 64*trace[Yb, Ybc, mb] - 64*trace[Ybc, Yb, mq] - 
    64*trace[Yt, Ytc, mt] - 64*trace[Ytc, Yt, mq]) + 
  MatMul[Ytc, Yt, Ytc, Yt]*(12*trace[ht, htc] + 36*mh2*trace[Ytc, Yt] + 
    12*trace[Yt, Ytc, mt] + 12*trace[Ytc, Yt, mq]) + 
  g2^2*MatMul[Ytc, Yt]*(36*trace[ht, htc] + 72*mh2*trace[Ytc, Yt] + 
    36*trace[Yt, Ytc, mt] + 36*trace[Ytc, Yt, mq]) + 
  g2^4*(-54*trace[hb, hbc] - 18*trace[he, hec] - 54*trace[ht, htc] - 
    54*mh1*trace[Ybc, Yb] - 54*mh2*trace[Ytc, Yt] - 
    18*mh1*trace[Adj[Ye], Ye] - 54*trace[Yb, Ybc, mb] - 
    54*trace[Ybc, Yb, mq] - 18*trace[Ye, Adj[Ye], me] - 
    54*trace[Yt, Ytc, mt] - 54*trace[Ytc, Yt, mq] - 
    18*trace[Adj[Ye], Ye, ml]) + g1^4*((-14*trace[hb, hbc])/25 - 
    (18*trace[he, hec])/25 - (26*trace[ht, htc])/25 - 
    (14*mh1*trace[Ybc, Yb])/25 - (26*mh2*trace[Ytc, Yt])/25 - 
    (18*mh1*trace[Adj[Ye], Ye])/25 - (14*trace[Yb, Ybc, mb])/25 - 
    (14*trace[Ybc, Yb, mq])/25 - (18*trace[Ye, Adj[Ye], me])/25 - 
    (26*trace[Yt, Ytc, mt])/25 - (26*trace[Ytc, Yt, mq])/25 - 
    (18*trace[Adj[Ye], Ye, ml])/25) + MatMul[Ybc, Yb, Ybc, Yb]*
   (12*trace[hb, hbc] + 4*trace[he, hec] + 36*mh1*trace[Ybc, Yb] + 
    12*mh1*trace[Adj[Ye], Ye] + 12*trace[Yb, Ybc, mb] + 
    12*trace[Ybc, Yb, mq] + 4*trace[Ye, Adj[Ye], me] + 
    4*trace[Adj[Ye], Ye, ml]) + g2^2*MatMul[Ybc, Yb]*
   (36*trace[hb, hbc] + 12*trace[he, hec] + 72*mh1*trace[Ybc, Yb] + 
    24*mh1*trace[Adj[Ye], Ye] + 36*trace[Yb, Ybc, mb] + 
    36*trace[Ybc, Yb, mq] + 12*trace[Ye, Adj[Ye], me] + 
    12*trace[Adj[Ye], Ye, ml]) + MatMul[htc, Yt]*
   (-36*trace[ht, Ytc]*trace[Ytc, Yt] + 12*trace[hb, Ytc, Yt, Ybc] + 
    12*trace[Ytc, ht, Ybc, Yb] + 72*trace[Ytc, ht, Ytc, Yt]) + 
  MatMul[Ytc, ht]*(-36*trace[htc, Yt]*trace[Ytc, Yt] + 
    12*trace[htc, Yt, Ybc, Yb] + 12*trace[Ytc, Yt, hbc, Yb] + 
    72*trace[Ytc, Yt, htc, Yt]) + 2*mq*MatMul[Ytc, Yt]*
   (-9*trace[Ytc, Yt]^2 + 6*trace[Ytc, Yt, Ybc, Yb] + 
    18*trace[Ytc, Yt, Ytc, Yt]) + mt*Yt*Ytc*(-18*trace[Ytc, Yt]^2 + 
    12*trace[Ytc, Yt, Ybc, Yb] + 36*trace[Ytc, Yt, Ytc, Yt]) + 
  MatMul[htc, ht]*(-18*trace[Ytc, Yt]^2 + 12*trace[Ytc, Yt, Ybc, Yb] + 
    36*trace[Ytc, Yt, Ytc, Yt]) + MatMul[hbc, Yb]*
   (-36*trace[hb, Ybc]*trace[Ybc, Yb] - 12*trace[he, Adj[Ye]]*
     trace[Ybc, Yb] - 12*trace[hb, Ybc]*trace[Adj[Ye], Ye] - 
    4*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] + 12*trace[hb, Ytc, Yt, Ybc] + 
    72*trace[Ybc, hb, Ybc, Yb] + 12*trace[Ytc, ht, Ybc, Yb] + 
    24*trace[Adj[Ye], he, Adj[Ye], Ye]) + 
  MatMul[Ybc, hb]*(-36*trace[hbc, Yb]*trace[Ybc, Yb] - 
    12*trace[hec, Ye]*trace[Ybc, Yb] - 12*trace[hbc, Yb]*trace[Adj[Ye], Ye] - 
    4*trace[hec, Ye]*trace[Adj[Ye], Ye] + 12*trace[htc, Yt, Ybc, Yb] + 
    72*trace[Ybc, Yb, hbc, Yb] + 12*trace[Ytc, Yt, hbc, Yb] + 
    24*trace[Adj[Ye], Ye, hec, Ye]) + 2*mq*MatMul[Ybc, Yb]*
   (-9*trace[Ybc, Yb]^2 - 6*trace[Ybc, Yb]*trace[Adj[Ye], Ye] - 
    trace[Adj[Ye], Ye]^2 + 18*trace[Ybc, Yb, Ybc, Yb] + 
    6*trace[Ytc, Yt, Ybc, Yb] + 6*trace[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  mb*Yb*Ybc*(-18*trace[Ybc, Yb]^2 - 12*trace[Ybc, Yb]*trace[Adj[Ye], Ye] - 
    2*trace[Adj[Ye], Ye]^2 + 36*trace[Ybc, Yb, Ybc, Yb] + 
    12*trace[Ytc, Yt, Ybc, Yb] + 12*trace[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  MatMul[hbc, hb]*(-18*trace[Ybc, Yb]^2 - 12*trace[Ybc, Yb]*
     trace[Adj[Ye], Ye] - 2*trace[Adj[Ye], Ye]^2 + 
    36*trace[Ybc, Yb, Ybc, Yb] + 12*trace[Ytc, Yt, Ybc, Yb] + 
    12*trace[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  MatMul[Ytc, Yt]*(-36*trace[ht, Ytc]*trace[htc, Yt] - 
    36*trace[ht, htc]*trace[Ytc, Yt] - 54*mh2*trace[Ytc, Yt]^2 - 
    36*trace[Ytc, Yt]*trace[Yt, Ytc, mt] - 36*trace[Ytc, Yt]*
     trace[Ytc, Yt, mq] + 12*trace[hb, htc, Yt, Ybc] + 
    12*trace[hbc, hb, Ytc, Yt] + 12*trace[htc, ht, Ybc, Yb] + 
    72*trace[htc, ht, Ytc, Yt] + 12*trace[Ytc, ht, hbc, Yb] + 
    72*trace[Ytc, ht, htc, Yt] + 12*mh1*trace[Ytc, Yt, Ybc, Yb] + 
    24*mh2*trace[Ytc, Yt, Ybc, Yb] + 108*mh2*trace[Ytc, Yt, Ytc, Yt] + 
    12*trace[Yb, Ytc, Yt, Ybc, mb] + 12*trace[Ybc, Yb, Ytc, Yt, mq] + 
    12*trace[Yt, Ybc, Yb, Ytc, mt] + 72*trace[Yt, Ytc, Yt, Ytc, mt] + 
    12*trace[Ytc, Yt, Ybc, Yb, mq] + 72*trace[Ytc, Yt, Ytc, Yt, mq]) + 
  MatMul[Ybc, Yb]*(-36*trace[hb, Ybc]*trace[hbc, Yb] - 
    12*trace[hbc, Yb]*trace[he, Adj[Ye]] - 12*trace[hb, Ybc]*trace[hec, Ye] - 
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
      ml]) + g2^4*g3^2*M2^2*(664 - 1296*Zeta[3]) + 
  g3^4*M3^2*MatMul[Ybc, Yb]*(32 - 1088*Zeta[3]) + 
  g3^4*M3^2*MatMul[Ytc, Yt]*(32 - 1088*Zeta[3]) + 
  g2^2*g3^4*M3^2*(192 - 864*Zeta[3]) + g2^4*g3^2*M2*M3*(400 - 864*Zeta[3]) + 
  g2^2*g3^4*M2*M3*(64 - 576*Zeta[3]) + g2^4*g3^2*M3^2*(272 - 432*Zeta[3]) + 
  g2^4*M2^2*MatMul[Ybc, Yb]*(-135 - 378*Zeta[3]) + 
  g2^4*M2^2*MatMul[Ytc, Yt]*(-135 - 378*Zeta[3]) + 
  g2^2*g3^4*M2^2*(80 - 288*Zeta[3]) + g1^2*g3^4*M3^2*
   (2464/15 - (1056*Zeta[3])/5) + g3^4*mb*Yb*Ybc*(16/3 - (544*Zeta[3])/3) + 
  g3^4*mt*Yt*Ytc*(16/3 - (544*Zeta[3])/3) + g3^4*MatMul[hbc, hb]*
   (16/3 - (544*Zeta[3])/3) + g3^4*MatMul[htc, ht]*(16/3 - (544*Zeta[3])/3) + 
  g1^2*g3^4*M1*M3*(4864/45 - (704*Zeta[3])/5) + 
  g1^2*g2^4*M2^2*(379/5 - (486*Zeta[3])/5) + 2*g3^4*mq*MatMul[Ybc, Yb]*
   (8/3 - (272*Zeta[3])/3) + 2*g3^4*mq*MatMul[Ytc, Yt]*
   (8/3 - (272*Zeta[3])/3) + g1^2*g3^4*M1^2*(592/9 - (352*Zeta[3])/5) + 
  g1^2*g2^4*M1*M2*(50 - (324*Zeta[3])/5) + 
  g2^4*mb*Yb*Ybc*(-45/2 - 63*Zeta[3]) + g2^4*mt*Yt*Ytc*(-45/2 - 63*Zeta[3]) + 
  g2^4*MatMul[hbc, hb]*(-45/2 - 63*Zeta[3]) + 
  g2^4*MatMul[htc, ht]*(-45/2 - 63*Zeta[3]) + 
  g2^2*M2*MatMul[hbc, Yb, Ybc, Yb]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[htc, Yt, Ytc, Yt]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Ybc, hb, Ybc, Yb]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Ybc, Yb, hbc, Yb]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Ybc, Yb, Ybc, hb]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Ytc, ht, Ytc, Yt]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Ytc, Yt, htc, Yt]*(6 - 36*Zeta[3]) + 
  g2^2*M2*MatMul[Ytc, Yt, Ytc, ht]*(6 - 36*Zeta[3]) + 
  g1^2*g2^4*M1^2*(152/5 - (162*Zeta[3])/5) + 2*g2^4*mq*MatMul[Ybc, Yb]*
   (-45/4 - (63*Zeta[3])/2) + 2*g2^4*mq*MatMul[Ytc, Yt]*
   (-45/4 - (63*Zeta[3])/2) + g1^2*M1^2*MatMul[Ytc, Yt, Ytc, Yt]*
   (44/3 - 24*Zeta[3]) + g1^4*g3^2*M1^2*(776/75 - (528*Zeta[3])/25) + 
  g1^6*M1^2*(57511/1125 - (2388*Zeta[3])/125) + 
  g1^2*g2^2*M1*MatMul[htc, Yt]*(59/5 - 18*Zeta[3]) + 
  g1^2*g2^2*M2*MatMul[htc, Yt]*(59/5 - 18*Zeta[3]) + 
  g1^2*g2^2*M1*MatMul[Ytc, ht]*(59/5 - 18*Zeta[3]) + 
  g1^2*g2^2*M2*MatMul[Ytc, ht]*(59/5 - 18*Zeta[3]) + 
  g1^2*g3^2*M1*MatMul[hbc, Yb]*(152/15 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[hbc, Yb]*(152/15 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M1*MatMul[Ybc, hb]*(152/15 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[Ybc, hb]*(152/15 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M1*MatMul[htc, Yt]*(136/5 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[htc, Yt]*(136/5 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M1*MatMul[Ytc, ht]*(136/5 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[Ytc, ht]*(136/5 - (256*Zeta[3])/15) + 
  g1^4*g3^2*M1*M3*(1552/225 - (352*Zeta[3])/25) + 
  g1^2*mq*MatMul[Ytc, Yt]^2*(22/3 - 12*Zeta[3]) + 
  g1^2*mt*Ytc*MatMul[Yt, Ytc, Yt]*(22/3 - 12*Zeta[3]) + 
  g1^2*mt*Yt*MatMul[Ytc, Yt, Ytc]*(22/3 - 12*Zeta[3]) + 
  g1^2*MatMul[htc, ht, Ytc, Yt]*(22/3 - 12*Zeta[3]) + 
  g1^2*MatMul[htc, Yt, Ytc, ht]*(22/3 - 12*Zeta[3]) + 
  g1^2*MatMul[Ytc, ht, htc, Yt]*(22/3 - 12*Zeta[3]) + 
  g1^2*MatMul[Ytc, Yt, htc, ht]*(22/3 - 12*Zeta[3]) + 
  g1^4*g3^2*M3^2*(208/45 - (176*Zeta[3])/25) + 
  g1^4*g2^2*M1^2*(33/25 - (162*Zeta[3])/25) + 
  2*g1^2*mq*MatMul[Ytc, Yt, Ytc, Yt]*(11/3 - 6*Zeta[3]) + 
  g1^2*M1^2*MatMul[Ybc, Yb, Ybc, Yb]*(28/15 - (24*Zeta[3])/5) + 
  g1^4*g2^2*M1*M2*(22/25 - (108*Zeta[3])/25) + 
  g1^4*M1*MatMul[htc, Yt]*(3767/75 - (286*Zeta[3])/75) + 
  g1^4*M1*MatMul[Ytc, ht]*(3767/75 - (286*Zeta[3])/75) + 
  g1^2*g2^2*M1*MatMul[hbc, Yb]*(41/5 - (18*Zeta[3])/5) + 
  g1^2*g2^2*M2*MatMul[hbc, Yb]*(41/5 - (18*Zeta[3])/5) + 
  g1^2*g2^2*M1*MatMul[Ybc, hb]*(41/5 - (18*Zeta[3])/5) + 
  g1^2*g2^2*M2*MatMul[Ybc, hb]*(41/5 - (18*Zeta[3])/5) + 
  g1^2*mq*MatMul[Ybc, Yb]^2*(14/15 - (12*Zeta[3])/5) + 
  g1^2*mb*Ybc*MatMul[Yb, Ybc, Yb]*(14/15 - (12*Zeta[3])/5) + 
  g1^2*mb*Yb*MatMul[Ybc, Yb, Ybc]*(14/15 - (12*Zeta[3])/5) + 
  g1^2*MatMul[hbc, hb, Ybc, Yb]*(14/15 - (12*Zeta[3])/5) + 
  g1^2*MatMul[hbc, Yb, Ybc, hb]*(14/15 - (12*Zeta[3])/5) + 
  g1^2*MatMul[Ybc, hb, hbc, Yb]*(14/15 - (12*Zeta[3])/5) + 
  g1^2*MatMul[Ybc, Yb, hbc, hb]*(14/15 - (12*Zeta[3])/5) + 
  g1^4*g2^2*M2^2*(4/5 - (54*Zeta[3])/25) + 2*g1^2*mq*MatMul[Ybc, Yb, Ybc, Yb]*
   (7/15 - (6*Zeta[3])/5) + g1^4*M1*MatMul[hbc, Yb]*
   (633/25 - (14*Zeta[3])/15) + g1^4*M1*MatMul[Ybc, hb]*
   (633/25 - (14*Zeta[3])/15) + 2*g1^4*mq*MatMul[Ybc, Yb]*
   (-633/100 + (7*Zeta[3])/30) + g1^4*mb*Yb*Ybc*(-633/50 + (7*Zeta[3])/15) + 
  g1^4*MatMul[hbc, hb]*(-633/50 + (7*Zeta[3])/15) + 
  2*g1^4*mq*MatMul[Ytc, Yt]*(-3767/300 + (143*Zeta[3])/150) + 
  12*mb*MatMul[Yb, Ybc, Yb]*MatMul[Ybc, Yb, Ybc]*Zeta[3] + 
  12*mt*MatMul[Yt, Ytc, Yt]*MatMul[Ytc, Yt, Ytc]*Zeta[3] + 
  24*mq*MatMul[Ybc, Yb]*MatMul[Ybc, Yb, Ybc, Yb]*Zeta[3] + 
  24*mq*MatMul[Ytc, Yt]*MatMul[Ytc, Yt, Ytc, Yt]*Zeta[3] + 
  12*mb*Ybc*MatMul[Yb, Ybc, Yb, Ybc, Yb]*Zeta[3] + 
  12*mb*Yb*MatMul[Ybc, Yb, Ybc, Yb, Ybc]*Zeta[3] + 
  12*mt*Ytc*MatMul[Yt, Ytc, Yt, Ytc, Yt]*Zeta[3] + 
  12*mt*Yt*MatMul[Ytc, Yt, Ytc, Yt, Ytc]*Zeta[3] + 
  12*MatMul[hbc, hb, Ybc, Yb, Ybc, Yb]*Zeta[3] + 
  12*MatMul[hbc, Yb, Ybc, hb, Ybc, Yb]*Zeta[3] + 
  12*MatMul[hbc, Yb, Ybc, Yb, Ybc, hb]*Zeta[3] + 
  12*MatMul[htc, ht, Ytc, Yt, Ytc, Yt]*Zeta[3] + 
  12*MatMul[htc, Yt, Ytc, ht, Ytc, Yt]*Zeta[3] + 
  12*MatMul[htc, Yt, Ytc, Yt, Ytc, ht]*Zeta[3] + 
  12*MatMul[Ybc, hb, hbc, Yb, Ybc, Yb]*Zeta[3] + 
  12*MatMul[Ybc, hb, Ybc, Yb, hbc, Yb]*Zeta[3] + 
  12*MatMul[Ybc, Yb, hbc, hb, Ybc, Yb]*Zeta[3] + 
  12*MatMul[Ybc, Yb, hbc, Yb, Ybc, hb]*Zeta[3] + 
  12*MatMul[Ybc, Yb, Ybc, hb, hbc, Yb]*Zeta[3] + 
  12*MatMul[Ybc, Yb, Ybc, Yb, hbc, hb]*Zeta[3] + 
  36*mh1*MatMul[Ybc, Yb, Ybc, Yb, Ybc, Yb]*Zeta[3] + 
  12*mq*MatMul[Ybc, Yb, Ybc, Yb, Ybc, Yb]*Zeta[3] + 
  12*MatMul[Ytc, ht, htc, Yt, Ytc, Yt]*Zeta[3] + 
  12*MatMul[Ytc, ht, Ytc, Yt, htc, Yt]*Zeta[3] + 
  12*MatMul[Ytc, Yt, htc, ht, Ytc, Yt]*Zeta[3] + 
  12*MatMul[Ytc, Yt, htc, Yt, Ytc, ht]*Zeta[3] + 
  12*MatMul[Ytc, Yt, Ytc, ht, htc, Yt]*Zeta[3] + 
  12*MatMul[Ytc, Yt, Ytc, Yt, htc, ht]*Zeta[3] + 
  36*mh2*MatMul[Ytc, Yt, Ytc, Yt, Ytc, Yt]*Zeta[3] + 
  12*mq*MatMul[Ytc, Yt, Ytc, Yt, Ytc, Yt]*Zeta[3] + 
  2*g1^2*g2^2*mq*MatMul[Ybc, Yb]*(-41/10 + (9*Zeta[3])/5) + 
  g1^4*mt*Yt*Ytc*(-3767/150 + (143*Zeta[3])/75) + 
  g1^4*MatMul[htc, ht]*(-3767/150 + (143*Zeta[3])/75) + 
  g1^2*M1*MatMul[hbc, Yb, Ybc, Yb]*(-14/15 + (12*Zeta[3])/5) + 
  g1^2*M1*MatMul[Ybc, hb, Ybc, Yb]*(-14/15 + (12*Zeta[3])/5) + 
  g1^2*M1*MatMul[Ybc, Yb, hbc, Yb]*(-14/15 + (12*Zeta[3])/5) + 
  g1^2*M1*MatMul[Ybc, Yb, Ybc, hb]*(-14/15 + (12*Zeta[3])/5) + 
  g1^4*M1^2*MatMul[Ybc, Yb]*(-1899/25 + (14*Zeta[3])/5) + 
  g1^2*g2^2*mb*Yb*Ybc*(-41/5 + (18*Zeta[3])/5) + 
  g1^2*g2^2*MatMul[hbc, hb]*(-41/5 + (18*Zeta[3])/5) + 
  g1^2*g2^2*M1^2*MatMul[Ybc, Yb]*(-82/5 + (36*Zeta[3])/5) + 
  g1^2*g2^2*M1*M2*MatMul[Ybc, Yb]*(-82/5 + (36*Zeta[3])/5) + 
  g1^2*g2^2*M2^2*MatMul[Ybc, Yb]*(-82/5 + (36*Zeta[3])/5) + 
  2*g1^2*g3^2*mq*MatMul[Ytc, Yt]*(-68/5 + (128*Zeta[3])/15) + 
  2*g1^2*g3^2*mq*MatMul[Ybc, Yb]*(-76/15 + (128*Zeta[3])/15) + 
  2*g1^2*g2^2*mq*MatMul[Ytc, Yt]*(-59/10 + 9*Zeta[3]) + 
  g1^4*M1^2*MatMul[Ytc, Yt]*(-3767/25 + (286*Zeta[3])/25) + 
  g1^2*M1*MatMul[htc, Yt, Ytc, Yt]*(-22/3 + 12*Zeta[3]) + 
  g1^2*M1*MatMul[Ytc, ht, Ytc, Yt]*(-22/3 + 12*Zeta[3]) + 
  g1^2*M1*MatMul[Ytc, Yt, htc, Yt]*(-22/3 + 12*Zeta[3]) + 
  g1^2*M1*MatMul[Ytc, Yt, Ytc, ht]*(-22/3 + 12*Zeta[3]) + 
  g1^2*g3^2*mt*Yt*Ytc*(-136/5 + (256*Zeta[3])/15) + 
  g1^2*g3^2*MatMul[htc, ht]*(-136/5 + (256*Zeta[3])/15) + 
  g1^2*g3^2*mb*Yb*Ybc*(-152/15 + (256*Zeta[3])/15) + 
  g1^2*g3^2*MatMul[hbc, hb]*(-152/15 + (256*Zeta[3])/15) + 
  g1^2*g2^2*mt*Yt*Ytc*(-59/5 + 18*Zeta[3]) + g1^2*g2^2*MatMul[htc, ht]*
   (-59/5 + 18*Zeta[3]) + 2*g2^2*mq*MatMul[Ybc, Yb, Ybc, Yb]*
   (-3 + 18*Zeta[3]) + 2*g2^2*mq*MatMul[Ytc, Yt, Ytc, Yt]*(-3 + 18*Zeta[3]) + 
  g1^2*g3^2*M1^2*MatMul[Ytc, Yt]*(-272/5 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M1*M3*MatMul[Ytc, Yt]*(-272/5 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M3^2*MatMul[Ytc, Yt]*(-272/5 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M1^2*MatMul[Ybc, Yb]*(-304/15 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M1*M3*MatMul[Ybc, Yb]*(-304/15 + (512*Zeta[3])/15) + 
  g1^2*g3^2*M3^2*MatMul[Ybc, Yb]*(-304/15 + (512*Zeta[3])/15) + 
  g1^2*g2^2*M1^2*MatMul[Ytc, Yt]*(-118/5 + 36*Zeta[3]) + 
  g1^2*g2^2*M1*M2*MatMul[Ytc, Yt]*(-118/5 + 36*Zeta[3]) + 
  g1^2*g2^2*M2^2*MatMul[Ytc, Yt]*(-118/5 + 36*Zeta[3]) + 
  g2^2*mq*MatMul[Ybc, Yb]^2*(-6 + 36*Zeta[3]) + 
  g2^2*mq*MatMul[Ytc, Yt]^2*(-6 + 36*Zeta[3]) + 
  g2^2*mb*Ybc*MatMul[Yb, Ybc, Yb]*(-6 + 36*Zeta[3]) + 
  g2^2*mb*Yb*MatMul[Ybc, Yb, Ybc]*(-6 + 36*Zeta[3]) + 
  g2^2*mt*Ytc*MatMul[Yt, Ytc, Yt]*(-6 + 36*Zeta[3]) + 
  g2^2*mt*Yt*MatMul[Ytc, Yt, Ytc]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[hbc, hb, Ybc, Yb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[hbc, Yb, Ybc, hb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[htc, ht, Ytc, Yt]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[htc, Yt, Ytc, ht]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Ybc, hb, hbc, Yb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Ybc, Yb, hbc, hb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Ytc, ht, htc, Yt]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Ytc, Yt, htc, ht]*(-6 + 36*Zeta[3]) + 
  g2^2*M2^2*MatMul[Ybc, Yb, Ybc, Yb]*(-12 + 72*Zeta[3]) + 
  g2^2*M2^2*MatMul[Ytc, Yt, Ytc, Yt]*(-12 + 72*Zeta[3]) + 
  g2^4*M2*MatMul[hbc, Yb]*(45 + 126*Zeta[3]) + 
  g2^4*M2*MatMul[htc, Yt]*(45 + 126*Zeta[3]) + 
  g2^4*M2*MatMul[Ybc, hb]*(45 + 126*Zeta[3]) + 
  g2^4*M2*MatMul[Ytc, ht]*(45 + 126*Zeta[3]) + 
  g3^4*M3*MatMul[hbc, Yb]*(-32/3 + (1088*Zeta[3])/3) + 
  g3^4*M3*MatMul[htc, Yt]*(-32/3 + (1088*Zeta[3])/3) + 
  g3^4*M3*MatMul[Ybc, hb]*(-32/3 + (1088*Zeta[3])/3) + 
  g3^4*M3*MatMul[Ytc, ht]*(-32/3 + (1088*Zeta[3])/3) + 
  g2^6*M2^2*(2157 + 3780*Zeta[3]) + g3^6*M3^2*(20512/9 + 7680*Zeta[3]) + 
  g3^4*MatMul[Ybc, Yb]*((16*mh1)/3 - (544*mh1*Zeta[3])/3) + 
  g2^4*MatMul[Ybc, Yb]*((-45*mh1)/2 - 63*mh1*Zeta[3]) + 
  g1^2*MatMul[Ybc, Yb, Ybc, Yb]*((28*mh1)/15 - (24*mh1*Zeta[3])/5) + 
  g1^4*MatMul[Ybc, Yb]*((-657*mh1)/50 - (12*mh2)/25 - (8*trace[mb])/25 - 
    (24*trace[me])/25 - (12*trace[ml])/25 - (4*trace[mq])/25 - 
    (32*trace[mt])/25 + (7*mh1*Zeta[3])/15) + g1^2*g2^2*MatMul[Ybc, Yb]*
   ((-41*mh1)/5 + (18*mh1*Zeta[3])/5) + g1^2*g3^2*MatMul[Ybc, Yb]*
   ((-152*mh1)/15 + (256*mh1*Zeta[3])/15) + g2^2*MatMul[Ybc, Yb, Ybc, Yb]*
   (-12*mh1 + 72*mh1*Zeta[3]) + g3^4*MatMul[Ytc, Yt]*
   ((16*mh2)/3 - (544*mh2*Zeta[3])/3) + g2^4*MatMul[Ytc, Yt]*
   ((-45*mh2)/2 - 63*mh2*Zeta[3]) + g1^2*MatMul[Ytc, Yt, Ytc, Yt]*
   ((44*mh2)/3 - 24*mh2*Zeta[3]) + g1^4*MatMul[Ytc, Yt]*
   ((-24*mh1)/25 - (3911*mh2)/150 - (16*trace[mb])/25 - (48*trace[me])/25 - 
    (24*trace[ml])/25 - (8*trace[mq])/25 - (64*trace[mt])/25 + 
    (143*mh2*Zeta[3])/75) + g1^2*g3^2*MatMul[Ytc, Yt]*
   ((-136*mh2)/5 + (256*mh2*Zeta[3])/15) + g1^2*g2^2*MatMul[Ytc, Yt]*
   ((-59*mh2)/5 + 18*mh2*Zeta[3]) + g2^2*MatMul[Ytc, Yt, Ytc, Yt]*
   (-12*mh2 + 72*mh2*Zeta[3]) + g3^2*MatMul[hbc, Yb]*
   (-16*trace[hb, Ybc] + 16*trace[he, Adj[Ye]] + 96*trace[hb, Ybc]*Zeta[3]) + 
  g3^2*M3*MatMul[Ybc, Yb]*(16*trace[hb, Ybc] + 16*trace[hbc, Yb] - 
    16*trace[he, Adj[Ye]] - 16*trace[hec, Ye] - 96*trace[hb, Ybc]*Zeta[3] - 
    96*trace[hbc, Yb]*Zeta[3]) + g3^2*MatMul[Ybc, hb]*
   (-16*trace[hbc, Yb] + 16*trace[hec, Ye] + 96*trace[hbc, Yb]*Zeta[3]) + 
  g1^2*MatMul[hbc, Yb]*((32*trace[hb, Ybc])/5 - (16*trace[he, Adj[Ye]])/5 - 
    (48*trace[hb, Ybc]*Zeta[3])/5 + (24*trace[he, Adj[Ye]]*Zeta[3])/5) + 
  g1^2*M1*MatMul[Ybc, Yb]*((-32*trace[hb, Ybc])/5 - (32*trace[hbc, Yb])/5 + 
    (16*trace[he, Adj[Ye]])/5 + (16*trace[hec, Ye])/5 + 
    (48*trace[hb, Ybc]*Zeta[3])/5 + (48*trace[hbc, Yb]*Zeta[3])/5 - 
    (24*trace[he, Adj[Ye]]*Zeta[3])/5 - (24*trace[hec, Ye]*Zeta[3])/5) + 
  g1^2*MatMul[Ybc, hb]*((32*trace[hbc, Yb])/5 - (16*trace[hec, Ye])/5 - 
    (48*trace[hbc, Yb]*Zeta[3])/5 + (24*trace[hec, Ye]*Zeta[3])/5) + 
  g1^2*MatMul[htc, Yt]*(4*trace[ht, Ytc] - (48*trace[ht, Ytc]*Zeta[3])/5) + 
  g3^2*MatMul[htc, Yt]*(-16*trace[ht, Ytc] + 96*trace[ht, Ytc]*Zeta[3]) + 
  g3^2*M3*MatMul[Ytc, Yt]*(16*trace[ht, Ytc] + 16*trace[htc, Yt] - 
    96*trace[ht, Ytc]*Zeta[3] - 96*trace[htc, Yt]*Zeta[3]) + 
  g1^2*MatMul[Ytc, ht]*(4*trace[htc, Yt] - (48*trace[htc, Yt]*Zeta[3])/5) + 
  g1^2*M1*MatMul[Ytc, Yt]*(-4*trace[ht, Ytc] - 4*trace[htc, Yt] + 
    (48*trace[ht, Ytc]*Zeta[3])/5 + (48*trace[htc, Yt]*Zeta[3])/5) + 
  g3^2*MatMul[Ytc, ht]*(-16*trace[htc, Yt] + 96*trace[htc, Yt]*Zeta[3]) + 
  g3^2*M3*MatMul[hbc, Yb]*(16*trace[Ybc, Yb] - 16*trace[Adj[Ye], Ye] - 
    96*trace[Ybc, Yb]*Zeta[3]) + g3^2*M3*MatMul[Ybc, hb]*
   (16*trace[Ybc, Yb] - 16*trace[Adj[Ye], Ye] - 96*trace[Ybc, Yb]*Zeta[3]) + 
  2*g3^2*mq*MatMul[Ybc, Yb]*(-8*trace[Ybc, Yb] + 8*trace[Adj[Ye], Ye] + 
    48*trace[Ybc, Yb]*Zeta[3]) + g3^2*mb*Yb*Ybc*(-16*trace[Ybc, Yb] + 
    16*trace[Adj[Ye], Ye] + 96*trace[Ybc, Yb]*Zeta[3]) + 
  g3^2*MatMul[hbc, hb]*(-16*trace[Ybc, Yb] + 16*trace[Adj[Ye], Ye] + 
    96*trace[Ybc, Yb]*Zeta[3]) + g3^2*M3^2*MatMul[Ybc, Yb]*
   (-32*trace[Ybc, Yb] + 32*trace[Adj[Ye], Ye] + 
    192*trace[Ybc, Yb]*Zeta[3]) + g3^2*M3*MatMul[htc, Yt]*
   (16*trace[Ytc, Yt] - 96*trace[Ytc, Yt]*Zeta[3]) + 
  g3^2*M3*MatMul[Ytc, ht]*(16*trace[Ytc, Yt] - 96*trace[Ytc, Yt]*Zeta[3]) + 
  g1^2*M1^2*MatMul[Ytc, Yt]*(8*trace[Ytc, Yt] - (96*trace[Ytc, Yt]*Zeta[3])/
     5) + g1^2*mt*Yt*Ytc*(4*trace[Ytc, Yt] - (48*trace[Ytc, Yt]*Zeta[3])/5) + 
  g1^2*MatMul[htc, ht]*(4*trace[Ytc, Yt] - (48*trace[Ytc, Yt]*Zeta[3])/5) + 
  2*g1^2*mq*MatMul[Ytc, Yt]*(2*trace[Ytc, Yt] - (24*trace[Ytc, Yt]*Zeta[3])/
     5) + g1^2*M1*MatMul[htc, Yt]*(-4*trace[Ytc, Yt] + 
    (48*trace[Ytc, Yt]*Zeta[3])/5) + g1^2*M1*MatMul[Ytc, ht]*
   (-4*trace[Ytc, Yt] + (48*trace[Ytc, Yt]*Zeta[3])/5) + 
  2*g3^2*mq*MatMul[Ytc, Yt]*(-8*trace[Ytc, Yt] + 48*trace[Ytc, Yt]*Zeta[3]) + 
  g3^2*mt*Yt*Ytc*(-16*trace[Ytc, Yt] + 96*trace[Ytc, Yt]*Zeta[3]) + 
  g3^2*MatMul[htc, ht]*(-16*trace[Ytc, Yt] + 96*trace[Ytc, Yt]*Zeta[3]) + 
  g3^2*M3^2*MatMul[Ytc, Yt]*(-32*trace[Ytc, Yt] + 
    192*trace[Ytc, Yt]*Zeta[3]) + g1^2*M1*MatMul[hbc, Yb]*
   ((-32*trace[Ybc, Yb])/5 + (16*trace[Adj[Ye], Ye])/5 + 
    (48*trace[Ybc, Yb]*Zeta[3])/5 - (24*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*M1*MatMul[Ybc, hb]*((-32*trace[Ybc, Yb])/5 + 
    (16*trace[Adj[Ye], Ye])/5 + (48*trace[Ybc, Yb]*Zeta[3])/5 - 
    (24*trace[Adj[Ye], Ye]*Zeta[3])/5) + 2*g1^2*mq*MatMul[Ybc, Yb]*
   ((16*trace[Ybc, Yb])/5 - (8*trace[Adj[Ye], Ye])/5 - 
    (24*trace[Ybc, Yb]*Zeta[3])/5 + (12*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*mb*Yb*Ybc*((32*trace[Ybc, Yb])/5 - (16*trace[Adj[Ye], Ye])/5 - 
    (48*trace[Ybc, Yb]*Zeta[3])/5 + (24*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*MatMul[hbc, hb]*((32*trace[Ybc, Yb])/5 - (16*trace[Adj[Ye], Ye])/5 - 
    (48*trace[Ybc, Yb]*Zeta[3])/5 + (24*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*M1^2*MatMul[Ybc, Yb]*((64*trace[Ybc, Yb])/5 - 
    (32*trace[Adj[Ye], Ye])/5 - (96*trace[Ybc, Yb]*Zeta[3])/5 + 
    (48*trace[Adj[Ye], Ye]*Zeta[3])/5) + g3^2*MatMul[Ybc, Yb]*
   (-16*trace[hb, hbc] + 16*trace[he, hec] - 32*mh1*trace[Ybc, Yb] + 
    32*mh1*trace[Adj[Ye], Ye] - 16*trace[Yb, Ybc, mb] - 
    16*trace[Ybc, Yb, mq] + 16*trace[Ye, Adj[Ye], me] + 
    16*trace[Adj[Ye], Ye, ml] + 96*trace[hb, hbc]*Zeta[3] + 
    192*mh1*trace[Ybc, Yb]*Zeta[3] + 96*trace[Yb, Ybc, mb]*Zeta[3] + 
    96*trace[Ybc, Yb, mq]*Zeta[3]) + g1^2*MatMul[Ytc, Yt]*
   (4*trace[ht, htc] + 8*mh2*trace[Ytc, Yt] + 4*trace[Yt, Ytc, mt] + 
    4*trace[Ytc, Yt, mq] - (48*trace[ht, htc]*Zeta[3])/5 - 
    (96*mh2*trace[Ytc, Yt]*Zeta[3])/5 - (48*trace[Yt, Ytc, mt]*Zeta[3])/5 - 
    (48*trace[Ytc, Yt, mq]*Zeta[3])/5) + g3^2*MatMul[Ytc, Yt]*
   (-16*trace[ht, htc] - 32*mh2*trace[Ytc, Yt] - 16*trace[Yt, Ytc, mt] - 
    16*trace[Ytc, Yt, mq] + 96*trace[ht, htc]*Zeta[3] + 
    192*mh2*trace[Ytc, Yt]*Zeta[3] + 96*trace[Yt, Ytc, mt]*Zeta[3] + 
    96*trace[Ytc, Yt, mq]*Zeta[3]) + g1^2*MatMul[Ybc, Yb]*
   ((32*trace[hb, hbc])/5 - (16*trace[he, hec])/5 + 
    (64*mh1*trace[Ybc, Yb])/5 - (32*mh1*trace[Adj[Ye], Ye])/5 + 
    (32*trace[Yb, Ybc, mb])/5 + (32*trace[Ybc, Yb, mq])/5 - 
    (16*trace[Ye, Adj[Ye], me])/5 - (16*trace[Adj[Ye], Ye, ml])/5 - 
    (48*trace[hb, hbc]*Zeta[3])/5 + (24*trace[he, hec]*Zeta[3])/5 - 
    (96*mh1*trace[Ybc, Yb]*Zeta[3])/5 + (48*mh1*trace[Adj[Ye], Ye]*Zeta[3])/
     5 - (48*trace[Yb, Ybc, mb]*Zeta[3])/5 - (48*trace[Ybc, Yb, mq]*Zeta[3])/
     5 + (24*trace[Ye, Adj[Ye], me]*Zeta[3])/5 + 
    (24*trace[Adj[Ye], Ye, ml]*Zeta[3])/5)}

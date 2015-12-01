{4*hb*hbc - (8*g1^2*M1^2)/15 - (32*g3^2*M3^2)/3 + 4*mq*Yb*Adj[Yb] + 
  4*mb*MatMul[Yb, Adj[Yb]] + 4*mh1*MatMul[Yb, Adj[Yb]], 
 (4*g1^2*hb*hbc)/5 + 12*g2^2*hb*hbc + (808*g1^4*M1^2)/75 + 
  (128*g1^2*g3^2*M1^2)/45 + (128*g1^2*g3^2*M1*M3)/45 + 
  (128*g1^2*g3^2*M3^2)/45 - (128*g3^4*M3^2)/3 - (4*g1^2*hbc*M1*Yb)/5 - 
  12*g2^2*hbc*M2*Yb + (4*g1^2*mq*Yb*Adj[Yb])/5 + 12*g2^2*mq*Yb*Adj[Yb] - 
  (4*g1^2*M1*MatMul[hb, Adj[Yb]])/5 - 12*g2^2*M2*MatMul[hb, Adj[Yb]] + 
  (8*g1^2*M1^2*MatMul[Yb, Adj[Yb]])/5 + 24*g2^2*M2^2*MatMul[Yb, Adj[Yb]] + 
  (4*g1^2*mb*MatMul[Yb, Adj[Yb]])/5 + 12*g2^2*mb*MatMul[Yb, Adj[Yb]] + 
  (4*g1^2*mh1*MatMul[Yb, Adj[Yb]])/5 + 12*g2^2*mh1*MatMul[Yb, Adj[Yb]] - 
  4*mb*MatMul[Yb, Adj[Yb]]^2 - 4*mt*MatMul[Yb, Adj[Yt]]*MatMul[Yt, Adj[Yb]] - 
  4*hbc*MatMul[hb, Adj[Yb], Yb] - 4*hbc*MatMul[hb, Adj[Yt], Yt] - 
  4*hbc*MatMul[Yb, Adj[Yb], hb] - 4*mq*Adj[Yb]*MatMul[Yb, Adj[Yb], Yb] - 
  4*hbc*MatMul[Yb, Adj[Yt], ht] - 4*mq*Adj[Yb]*MatMul[Yb, Adj[Yt], Yt] - 
  4*mq*Yb*MatMul[Adj[Yb], Yb, Adj[Yb]] - 
  4*mq*Yb*MatMul[Adj[Yt], Yt, Adj[Yb]] - 4*MatMul[hb, hbc, Yb, Adj[Yb]] - 
  4*MatMul[hb, htc, Yt, Adj[Yb]] - 4*MatMul[Yb, hbc, hb, Adj[Yb]] - 
  4*MatMul[Yb, htc, ht, Adj[Yb]] - 4*mb*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]] - 
  8*mh1*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]] - 
  4*mb*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]] + (-4*mh1 - 4*mh2)*
   MatMul[Yb, Adj[Yt], Yt, Adj[Yb]] + 
  g1^4*((4*mh1)/25 + (4*mh2)/25 + (8*trace[mb])/75 + (8*trace[me])/25 + 
    (4*trace[ml])/25 + (4*trace[mq])/75 + (32*trace[mt])/75) + 
  g3^4*((16*trace[mb])/3 + (32*trace[mq])/3 + (16*trace[mt])/3) + 
  Yb*(-12*hbc*trace[hb, Adj[Yb]] - 4*hbc*trace[he, Adj[Ye]]) + 
  MatMul[hb, Adj[Yb]]*(-12*trace[hbc, Yb] - 4*trace[hec, Ye]) + 
  mq*Yb*Adj[Yb]*(-12*trace[Adj[Yb], Yb] - 4*trace[Adj[Ye], Ye]) + 
  2*mb*MatMul[Yb, Adj[Yb]]*(-6*trace[Adj[Yb], Yb] - 2*trace[Adj[Ye], Ye]) + 
  hb*(-12*hbc*trace[Adj[Yb], Yb] - 4*hbc*trace[Adj[Ye], Ye]) + 
  MatMul[Yb, Adj[Yb]]*(-12*trace[hb, hbc] - 4*trace[he, hec] - 
    24*mh1*trace[Adj[Yb], Yb] - 8*mh1*trace[Adj[Ye], Ye] - 
    12*trace[Yb, Adj[Yb], mb] - 4*trace[Ye, Adj[Ye], me] - 
    12*trace[Adj[Yb], Yb, mq] - 4*trace[Adj[Ye], Ye, ml]), 
 (128*g3^2*mb*MatMul[Yb, Adj[Yb]]^2)/3 + 
  (128*g3^2*mt*MatMul[Yb, Adj[Yt]]*MatMul[Yt, Adj[Yb]])/3 + 
  (128*g3^2*hbc*MatMul[hb, Adj[Yb], Yb])/3 + 
  (128*g3^2*hbc*MatMul[hb, Adj[Yt], Yt])/3 + 
  (128*g3^2*hbc*MatMul[Yb, Adj[Yb], hb])/3 - 
  (128*g3^2*hbc*M3*MatMul[Yb, Adj[Yb], Yb])/3 + 
  (128*g3^2*mq*Adj[Yb]*MatMul[Yb, Adj[Yb], Yb])/3 + 
  (128*g3^2*hbc*MatMul[Yb, Adj[Yt], ht])/3 - 
  (128*g3^2*hbc*M3*MatMul[Yb, Adj[Yt], Yt])/3 + 
  (128*g3^2*mq*Adj[Yb]*MatMul[Yb, Adj[Yt], Yt])/3 + 
  (128*g3^2*mq*Yb*MatMul[Adj[Yb], Yb, Adj[Yb]])/3 - 
  4*mq*MatMul[Yb, Adj[Yt], Yt]*MatMul[Adj[Yb], Yb, Adj[Yb]] + 
  (128*g3^2*mq*Yb*MatMul[Adj[Yt], Yt, Adj[Yb]])/3 - 
  4*mq*MatMul[Yb, Adj[Yb], Yb]*MatMul[Adj[Yt], Yt, Adj[Yb]] + 
  12*mq*MatMul[Yb, Adj[Yt], Yt]*MatMul[Adj[Yt], Yt, Adj[Yb]] + 
  (128*g3^2*MatMul[hb, hbc, Yb, Adj[Yb]])/3 + 
  (128*g3^2*MatMul[hb, htc, Yt, Adj[Yb]])/3 - 
  (128*g3^2*M3*MatMul[hb, Adj[Yb], Yb, Adj[Yb]])/3 - 
  (128*g3^2*M3*MatMul[hb, Adj[Yt], Yt, Adj[Yb]])/3 + 
  (128*g3^2*MatMul[Yb, hbc, hb, Adj[Yb]])/3 - 
  (128*g3^2*M3*MatMul[Yb, hbc, Yb, Adj[Yb]])/3 + 
  (128*g3^2*MatMul[Yb, htc, ht, Adj[Yb]])/3 - 
  (128*g3^2*M3*MatMul[Yb, htc, Yt, Adj[Yb]])/3 - 
  (128*g3^2*M3*MatMul[Yb, Adj[Yb], hb, Adj[Yb]])/3 + 
  (256*g3^2*M3^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]])/3 + 
  (128*g3^2*mb*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]])/3 + 
  (256*g3^2*mh1*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]])/3 - 
  4*mt*MatMul[Yt, Adj[Yb]]*MatMul[Yb, Adj[Yb], Yb, Adj[Yt]] - 
  (128*g3^2*M3*MatMul[Yb, Adj[Yt], ht, Adj[Yb]])/3 + 
  (256*g3^2*M3^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]])/3 + 
  (128*g3^2*mb*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]])/3 + 
  g3^2*((128*mh1)/3 + (128*mh2)/3)*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]] - 
  8*mb*MatMul[Yb, Adj[Yb]]*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]] + 
  12*mt*MatMul[Yt, Adj[Yb]]*MatMul[Yb, Adj[Yt], Yt, Adj[Yt]] - 
  4*mt*MatMul[Yb, Adj[Yt]]*MatMul[Yt, Adj[Yb], Yb, Adj[Yb]] + 
  12*mt*MatMul[Yb, Adj[Yt]]*MatMul[Yt, Adj[Yt], Yt, Adj[Yb]] - 
  4*hbc*MatMul[hb, Adj[Yb], Yb, Adj[Yt], Yt] - 
  4*hbc*MatMul[hb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*hbc*MatMul[hb, Adj[Yt], Yt, Adj[Yt], Yt] - 
  4*hbc*MatMul[Yb, Adj[Yb], hb, Adj[Yt], Yt] - 
  4*hbc*MatMul[Yb, Adj[Yb], Yb, Adj[Yt], ht] - 
  4*mq*Adj[Yb]*MatMul[Yb, Adj[Yb], Yb, Adj[Yt], Yt] - 
  4*hbc*MatMul[Yb, Adj[Yt], ht, Adj[Yb], Yb] + 
  12*hbc*MatMul[Yb, Adj[Yt], ht, Adj[Yt], Yt] - 
  4*hbc*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], hb] - 
  4*mq*Adj[Yb]*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*hbc*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], ht] + 
  12*mq*Adj[Yb]*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt] - 
  4*mq*Yb*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb]] - 
  4*mq*Yb*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yb]] + 
  12*mq*Yb*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb]] - 
  4*MatMul[hb, hbc, Yb, Adj[Yt], Yt, Adj[Yb]] - 
  4*MatMul[hb, htc, Yt, Adj[Yb], Yb, Adj[Yb]] + 
  12*MatMul[hb, htc, Yt, Adj[Yt], Yt, Adj[Yb]] - 
  4*MatMul[hb, Adj[Yb], Yb, htc, Yt, Adj[Yb]] - 
  4*MatMul[hb, Adj[Yt], Yt, hbc, Yb, Adj[Yb]] + 
  12*MatMul[hb, Adj[Yt], Yt, htc, Yt, Adj[Yb]] - 
  4*MatMul[Yb, hbc, hb, Adj[Yt], Yt, Adj[Yb]] - 
  4*MatMul[Yb, hbc, Yb, Adj[Yt], ht, Adj[Yb]] - 
  4*MatMul[Yb, htc, ht, Adj[Yb], Yb, Adj[Yb]] + 
  12*MatMul[Yb, htc, ht, Adj[Yt], Yt, Adj[Yb]] - 
  4*MatMul[Yb, htc, Yt, Adj[Yb], hb, Adj[Yb]] + 
  12*MatMul[Yb, htc, Yt, Adj[Yt], ht, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yb], hb, htc, Yt, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yb], Yb, htc, ht, Adj[Yb]] - 
  4*mb*MatMul[Yb, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb]] + 
  (-8*mh1 - 4*mh2)*MatMul[Yb, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yt], ht, hbc, Yb, Adj[Yb]] + 
  12*MatMul[Yb, Adj[Yt], ht, htc, Yt, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yt], Yt, hbc, hb, Adj[Yb]] + 
  12*MatMul[Yb, Adj[Yt], Yt, htc, ht, Adj[Yb]] - 
  4*mb*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yb]] + 
  (-8*mh1 - 4*mh2)*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yb]] + 
  12*mb*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb]] + 
  (12*mh1 + 24*mh2)*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb]] + 
  g1^6*((-808*mh1)/375 - (808*mh2)/375 - (1616*trace[mb])/1125 - 
    (1616*trace[me])/375 - (808*trace[ml])/375 - (808*trace[mq])/1125 - 
    (6464*trace[mt])/1125) + g1^4*g3^2*((-64*mh1)/75 - (64*mh2)/75 - 
    (128*trace[mb])/225 - (128*trace[me])/75 - (64*trace[ml])/75 - 
    (64*trace[mq])/225 - (512*trace[mt])/225) + 
  g1^2*g3^4*((-64*trace[mb])/45 - (128*trace[mq])/45 - (64*trace[mt])/45) + 
  g3^6*((320*trace[mb])/9 + (640*trace[mq])/9 + (320*trace[mt])/9) + 
  MatMul[Yb, hbc, Yb, Adj[Yb]]*(12*trace[hb, Adj[Yb]] + 
    4*trace[he, Adj[Ye]]) + MatMul[Yb, Adj[Yb], Yb]*
   (12*hbc*trace[hb, Adj[Yb]] + 4*hbc*trace[he, Adj[Ye]]) + 
  MatMul[hb, Adj[Yb], Yb, Adj[Yb]]*(12*trace[hbc, Yb] + 4*trace[hec, Ye]) + 
  MatMul[Yb, Adj[Yb], hb, Adj[Yb]]*(12*trace[hbc, Yb] + 4*trace[hec, Ye]) + 
  MatMul[Yb, htc, Yt, Adj[Yb]]*(-12*trace[hb, Adj[Yb]] - 
    4*trace[he, Adj[Ye]] + 24*trace[ht, Adj[Yt]]) + 
  MatMul[Yb, Adj[Yt], Yt]*(-12*hbc*trace[hb, Adj[Yb]] - 
    4*hbc*trace[he, Adj[Ye]] + 24*hbc*trace[ht, Adj[Yt]]) + 
  g1^4*M1*((56*trace[hb, Adj[Yb]])/15 + (56*trace[hbc, Yb])/15 + 
    (24*trace[he, Adj[Ye]])/5 + (24*trace[hec, Ye])/5 + 
    (104*trace[ht, Adj[Yt]])/15 + (104*trace[htc, Yt])/15) + 
  MatMul[hb, Adj[Yt], Yt, Adj[Yb]]*(-12*trace[hbc, Yb] - 4*trace[hec, Ye] + 
    24*trace[htc, Yt]) + MatMul[Yb, Adj[Yt], ht, Adj[Yb]]*
   (-12*trace[hbc, Yb] - 4*trace[hec, Ye] + 24*trace[htc, Yt]) + 
  g3^4*M3*((320*trace[hb, Adj[Yb]])/3 + (320*trace[hbc, Yb])/3 + 
    (320*trace[ht, Adj[Yt]])/3 + (320*trace[htc, Yt])/3) + 
  2*mb*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*(6*trace[Adj[Yb], Yb] + 
    2*trace[Adj[Ye], Ye]) + mb*MatMul[Yb, Adj[Yb]]^2*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  mq*Adj[Yb]*MatMul[Yb, Adj[Yb], Yb]*(12*trace[Adj[Yb], Yb] + 
    4*trace[Adj[Ye], Ye]) + mq*Yb*MatMul[Adj[Yb], Yb, Adj[Yb]]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[hb, hbc, Yb, Adj[Yb]]*(12*trace[Adj[Yb], Yb] + 
    4*trace[Adj[Ye], Ye]) + MatMul[Yb, hbc, hb, Adj[Yb]]*
   (12*trace[Adj[Yb], Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[hb, Adj[Yb], Yb]*(12*hbc*trace[Adj[Yb], Yb] + 
    4*hbc*trace[Adj[Ye], Ye]) + MatMul[Yb, Adj[Yb], hb]*
   (12*hbc*trace[Adj[Yb], Yb] + 4*hbc*trace[Adj[Ye], Ye]) + 
  g3^4*M3^2*(-320*trace[Adj[Yb], Yb] - 320*trace[Adj[Yt], Yt]) + 
  g1^4*M1^2*((-56*trace[Adj[Yb], Yb])/5 - (72*trace[Adj[Ye], Ye])/5 - 
    (104*trace[Adj[Yt], Yt])/5) + 2*mb*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*
   (-6*trace[Adj[Yb], Yb] - 2*trace[Adj[Ye], Ye] + 12*trace[Adj[Yt], Yt]) + 
  mt*MatMul[Yb, Adj[Yt]]*MatMul[Yt, Adj[Yb]]*(-12*trace[Adj[Yb], Yb] - 
    4*trace[Adj[Ye], Ye] + 24*trace[Adj[Yt], Yt]) + 
  mq*Adj[Yb]*MatMul[Yb, Adj[Yt], Yt]*(-12*trace[Adj[Yb], Yb] - 
    4*trace[Adj[Ye], Ye] + 24*trace[Adj[Yt], Yt]) + 
  mq*Yb*MatMul[Adj[Yt], Yt, Adj[Yb]]*(-12*trace[Adj[Yb], Yb] - 
    4*trace[Adj[Ye], Ye] + 24*trace[Adj[Yt], Yt]) + 
  MatMul[hb, htc, Yt, Adj[Yb]]*(-12*trace[Adj[Yb], Yb] - 
    4*trace[Adj[Ye], Ye] + 24*trace[Adj[Yt], Yt]) + 
  MatMul[Yb, htc, ht, Adj[Yb]]*(-12*trace[Adj[Yb], Yb] - 
    4*trace[Adj[Ye], Ye] + 24*trace[Adj[Yt], Yt]) + 
  MatMul[hb, Adj[Yt], Yt]*(-12*hbc*trace[Adj[Yb], Yb] - 
    4*hbc*trace[Adj[Ye], Ye] + 24*hbc*trace[Adj[Yt], Yt]) + 
  MatMul[Yb, Adj[Yt], ht]*(-12*hbc*trace[Adj[Yb], Yb] - 
    4*hbc*trace[Adj[Ye], Ye] + 24*hbc*trace[Adj[Yt], Yt]) + 
  MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*(12*trace[hb, hbc] + 4*trace[he, hec] + 
    36*mh1*trace[Adj[Yb], Yb] + 12*mh1*trace[Adj[Ye], Ye] + 
    12*trace[Yb, Adj[Yb], mb] + 4*trace[Ye, Adj[Ye], me] + 
    12*trace[Adj[Yb], Yb, mq] + 4*trace[Adj[Ye], Ye, ml]) + 
  g3^4*(-64*trace[hb, hbc] - 64*trace[ht, htc] - 64*mh1*trace[Adj[Yb], Yb] - 
    64*mh2*trace[Adj[Yt], Yt] - 64*trace[Yb, Adj[Yb], mb] - 
    64*trace[Yt, Adj[Yt], mt] - 64*trace[Adj[Yb], Yb, mq] - 
    64*trace[Adj[Yt], Yt, mq]) + g1^4*((-56*trace[hb, hbc])/25 - 
    (72*trace[he, hec])/25 - (104*trace[ht, htc])/25 - 
    (56*mh1*trace[Adj[Yb], Yb])/25 - (72*mh1*trace[Adj[Ye], Ye])/25 - 
    (104*mh2*trace[Adj[Yt], Yt])/25 - (56*trace[Yb, Adj[Yb], mb])/25 - 
    (72*trace[Ye, Adj[Ye], me])/25 - (104*trace[Yt, Adj[Yt], mt])/25 - 
    (56*trace[Adj[Yb], Yb, mq])/25 - (72*trace[Adj[Ye], Ye, ml])/25 - 
    (104*trace[Adj[Yt], Yt, mq])/25) + MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*
   (-12*trace[hb, hbc] - 4*trace[he, hec] + 24*trace[ht, htc] - 
    24*mh1*trace[Adj[Yb], Yb] - 12*mh2*trace[Adj[Yb], Yb] - 
    8*mh1*trace[Adj[Ye], Ye] - 4*mh2*trace[Adj[Ye], Ye] + 
    24*mh1*trace[Adj[Yt], Yt] + 48*mh2*trace[Adj[Yt], Yt] - 
    12*trace[Yb, Adj[Yb], mb] - 4*trace[Ye, Adj[Ye], me] + 
    24*trace[Yt, Adj[Yt], mt] - 12*trace[Adj[Yb], Yb, mq] - 
    4*trace[Adj[Ye], Ye, ml] + 24*trace[Adj[Yt], Yt, mq]) + 
  Yb*(-72*hbc*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
    24*hbc*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb] - 
    24*hbc*trace[hb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
    8*hbc*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] + 
    24*hbc*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
    144*hbc*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
    48*hbc*trace[Adj[Ye], he, Adj[Ye], Ye] + 
    24*hbc*trace[Adj[Yt], ht, Adj[Yb], Yb]) + 
  MatMul[hb, Adj[Yb]]*(-72*trace[hbc, Yb]*trace[Adj[Yb], Yb] - 
    24*trace[hec, Ye]*trace[Adj[Yb], Yb] - 24*trace[hbc, Yb]*
     trace[Adj[Ye], Ye] - 8*trace[hec, Ye]*trace[Adj[Ye], Ye] + 
    24*trace[htc, Yt, Adj[Yb], Yb] + 144*trace[Adj[Yb], Yb, hbc, Yb] + 
    48*trace[Adj[Ye], Ye, hec, Ye] + 24*trace[Adj[Yt], Yt, hbc, Yb]) + 
  2*mb*MatMul[Yb, Adj[Yb]]*(-18*trace[Adj[Yb], Yb]^2 - 
    12*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 2*trace[Adj[Ye], Ye]^2 + 
    36*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 12*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    12*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + 
  mq*Yb*Adj[Yb]*(-36*trace[Adj[Yb], Yb]^2 - 24*trace[Adj[Yb], Yb]*
     trace[Adj[Ye], Ye] - 4*trace[Adj[Ye], Ye]^2 + 
    72*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 24*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    24*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + 
  hb*(-36*hbc*trace[Adj[Yb], Yb]^2 - 24*hbc*trace[Adj[Yb], Yb]*
     trace[Adj[Ye], Ye] - 4*hbc*trace[Adj[Ye], Ye]^2 + 
    72*hbc*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    24*hbc*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    24*hbc*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + 
  MatMul[Yb, Adj[Yb]]*(-72*trace[hb, Adj[Yb]]*trace[hbc, Yb] - 
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
    24*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq]) + g3^4*M3^2*MatMul[Yb, Adj[Yb]]*
   (64 - 2176*Zeta[3]) + g2^2*g3^4*M3^2*(720 - 864*Zeta[3]) + 
  g2^2*g3^4*M2*M3*(480 - 576*Zeta[3]) + g3^4*mq*Yb*Adj[Yb]*
   (32/3 - (1088*Zeta[3])/3) + g2^2*g3^4*M2^2*(288 - 288*Zeta[3]) + 
  g1^2*g3^4*M3^2*(1936/15 - (1056*Zeta[3])/5) + 
  g2^2*g3^2*M2*MatMul[hb, Adj[Yb]]*(176 - 192*Zeta[3]) + 
  g2^2*g3^2*M3*MatMul[hb, Adj[Yb]]*(176 - 192*Zeta[3]) + 
  2*g3^4*mb*MatMul[Yb, Adj[Yb]]*(16/3 - (544*Zeta[3])/3) + 
  g1^2*g3^4*M1*M3*(3616/45 - (704*Zeta[3])/5) + 
  g2^4*M2^2*MatMul[Yb, Adj[Yb]]*(-474 - 108*Zeta[3]) + 
  g1^4*g3^2*M1^2*(2912/75 - (2112*Zeta[3])/25) + 
  g1^6*M1^2*(227548/1125 - (9552*Zeta[3])/125) + 
  g2^2*M2^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*(36 - 72*Zeta[3]) + 
  g2^2*M2^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*(36 - 72*Zeta[3]) + 
  g1^2*g3^4*M1^2*(2336/45 - (352*Zeta[3])/5) + 
  g1^4*g3^2*M1*M3*(5824/225 - (1408*Zeta[3])/25) + 
  g2^2*mb*MatMul[Yb, Adj[Yb]]^2*(18 - 36*Zeta[3]) + 
  g2^2*mt*MatMul[Yb, Adj[Yt]]*MatMul[Yt, Adj[Yb]]*(18 - 36*Zeta[3]) + 
  g2^2*mq*Adj[Yb]*MatMul[Yb, Adj[Yb], Yb]*(18 - 36*Zeta[3]) + 
  g2^2*mq*Adj[Yb]*MatMul[Yb, Adj[Yt], Yt]*(18 - 36*Zeta[3]) + 
  g2^2*mq*Yb*MatMul[Adj[Yb], Yb, Adj[Yb]]*(18 - 36*Zeta[3]) + 
  g2^2*mq*Yb*MatMul[Adj[Yt], Yt, Adj[Yb]]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[hb, hbc, Yb, Adj[Yb]]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[hb, htc, Yt, Adj[Yb]]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[Yb, hbc, hb, Adj[Yb]]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[Yb, htc, ht, Adj[Yb]]*(18 - 36*Zeta[3]) + 
  g1^4*g3^2*M3^2*(3968/225 - (704*Zeta[3])/25) + 
  g1^4*g2^2*M1^2*(108/5 - (648*Zeta[3])/25) + 
  g1^2*g3^2*M1*MatMul[hb, Adj[Yb]]*(48/5 - (64*Zeta[3])/3) + 
  g1^2*g3^2*M3*MatMul[hb, Adj[Yb]]*(48/5 - (64*Zeta[3])/3) + 
  g2^4*mq*Yb*Adj[Yb]*(-87 - 18*Zeta[3]) + 
  2*g2^2*mb*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*(9 - 18*Zeta[3]) + 
  2*g2^2*mb*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*(9 - 18*Zeta[3]) + 
  g1^4*g2^2*M1*M2*(72/5 - (432*Zeta[3])/25) + 
  g1^2*g2^2*M1*MatMul[hb, Adj[Yb]]*(86/5 - (84*Zeta[3])/5) + 
  g1^2*g2^2*M2*MatMul[hb, Adj[Yb]]*(86/5 - (84*Zeta[3])/5) + 
  2*g2^4*mb*MatMul[Yb, Adj[Yb]]*(-87/2 - 9*Zeta[3]) + 
  g1^4*g2^2*M2^2*(216/25 - (216*Zeta[3])/25) + 
  g1^2*M1*MatMul[hb, Adj[Yt], Yt, Adj[Yb]]*(58/15 - (36*Zeta[3])/5) + 
  g1^2*M1*MatMul[Yb, htc, Yt, Adj[Yb]]*(58/15 - (36*Zeta[3])/5) + 
  g1^2*M1*MatMul[Yb, Adj[Yt], ht, Adj[Yb]]*(58/15 - (36*Zeta[3])/5) + 
  g1^2*M1*MatMul[hb, Adj[Yb], Yb, Adj[Yb]]*(2/3 - (12*Zeta[3])/5) + 
  g1^2*M1*MatMul[Yb, hbc, Yb, Adj[Yb]]*(2/3 - (12*Zeta[3])/5) + 
  g1^2*M1*MatMul[Yb, Adj[Yb], hb, Adj[Yb]]*(2/3 - (12*Zeta[3])/5) + 
  g1^4*M1^2*MatMul[Yb, Adj[Yb]]*(-674/5 - (28*Zeta[3])/25) + 
  g1^4*mq*Yb*Adj[Yb]*(-337/15 - (14*Zeta[3])/75) + 
  2*g1^4*mb*MatMul[Yb, Adj[Yb]]*(-337/30 - (7*Zeta[3])/75) + 
  g1^4*M1*MatMul[hb, Adj[Yb]]*(674/15 + (28*Zeta[3])/75) + 
  2*g1^2*mb*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*(-1/3 + (6*Zeta[3])/5) + 
  g1^2*mb*MatMul[Yb, Adj[Yb]]^2*(-2/3 + (12*Zeta[3])/5) + 
  g1^2*mq*Adj[Yb]*MatMul[Yb, Adj[Yb], Yb]*(-2/3 + (12*Zeta[3])/5) + 
  g1^2*mq*Yb*MatMul[Adj[Yb], Yb, Adj[Yb]]*(-2/3 + (12*Zeta[3])/5) + 
  g1^2*MatMul[hb, hbc, Yb, Adj[Yb]]*(-2/3 + (12*Zeta[3])/5) + 
  g1^2*MatMul[Yb, hbc, hb, Adj[Yb]]*(-2/3 + (12*Zeta[3])/5) + 
  2*g1^2*mb*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*(-29/15 + (18*Zeta[3])/5) + 
  g1^2*M1^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*(-4/3 + (24*Zeta[3])/5) + 
  g1^2*mt*MatMul[Yb, Adj[Yt]]*MatMul[Yt, Adj[Yb]]*(-58/15 + (36*Zeta[3])/5) + 
  g1^2*mq*Adj[Yb]*MatMul[Yb, Adj[Yt], Yt]*(-58/15 + (36*Zeta[3])/5) + 
  g1^2*mq*Yb*MatMul[Adj[Yt], Yt, Adj[Yb]]*(-58/15 + (36*Zeta[3])/5) + 
  g1^2*MatMul[hb, htc, Yt, Adj[Yb]]*(-58/15 + (36*Zeta[3])/5) + 
  g1^2*MatMul[Yb, htc, ht, Adj[Yb]]*(-58/15 + (36*Zeta[3])/5) + 
  2*g1^2*g2^2*mb*MatMul[Yb, Adj[Yb]]*(-43/5 + (42*Zeta[3])/5) + 
  2*g1^2*g3^2*mb*MatMul[Yb, Adj[Yb]]*(-24/5 + (32*Zeta[3])/3) + 
  2*mb*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb]]*(6 + 12*Zeta[3]) + 
  g1^2*M1^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*(-116/15 + (72*Zeta[3])/5) + 
  g1^2*g2^2*mq*Yb*Adj[Yb]*(-86/5 + (84*Zeta[3])/5) + 
  g1^2*g3^2*mq*Yb*Adj[Yb]*(-48/5 + (64*Zeta[3])/3) + 
  mq*MatMul[Yb, Adj[Yb], Yb]*MatMul[Adj[Yb], Yb, Adj[Yb]]*(12 + 24*Zeta[3]) + 
  2*mb*MatMul[Yb, Adj[Yb]]*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*
   (12 + 24*Zeta[3]) + mq*Adj[Yb]*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb]*
   (12 + 24*Zeta[3]) + mq*Yb*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb]]*
   (12 + 24*Zeta[3]) + MatMul[hb, hbc, Yb, Adj[Yb], Yb, Adj[Yb]]*
   (12 + 24*Zeta[3]) + MatMul[hb, Adj[Yb], Yb, hbc, Yb, Adj[Yb]]*
   (12 + 24*Zeta[3]) + MatMul[Yb, hbc, hb, Adj[Yb], Yb, Adj[Yb]]*
   (12 + 24*Zeta[3]) + MatMul[Yb, hbc, Yb, Adj[Yb], hb, Adj[Yb]]*
   (12 + 24*Zeta[3]) + MatMul[Yb, Adj[Yb], hb, hbc, Yb, Adj[Yb]]*
   (12 + 24*Zeta[3]) + MatMul[Yb, Adj[Yb], Yb, hbc, hb, Adj[Yb]]*
   (12 + 24*Zeta[3]) + g1^2*g2^2*M1^2*MatMul[Yb, Adj[Yb]]*
   (-172/5 + (168*Zeta[3])/5) + g1^2*g2^2*M1*M2*MatMul[Yb, Adj[Yb]]*
   (-172/5 + (168*Zeta[3])/5) + g1^2*g2^2*M2^2*MatMul[Yb, Adj[Yb]]*
   (-172/5 + (168*Zeta[3])/5) + g2^2*M2*MatMul[hb, Adj[Yb], Yb, Adj[Yb]]*
   (-18 + 36*Zeta[3]) + g2^2*M2*MatMul[hb, Adj[Yt], Yt, Adj[Yb]]*
   (-18 + 36*Zeta[3]) + g2^2*M2*MatMul[Yb, hbc, Yb, Adj[Yb]]*
   (-18 + 36*Zeta[3]) + g2^2*M2*MatMul[Yb, htc, Yt, Adj[Yb]]*
   (-18 + 36*Zeta[3]) + g2^2*M2*MatMul[Yb, Adj[Yb], hb, Adj[Yb]]*
   (-18 + 36*Zeta[3]) + g2^2*M2*MatMul[Yb, Adj[Yt], ht, Adj[Yb]]*
   (-18 + 36*Zeta[3]) + g2^4*M2*MatMul[hb, Adj[Yb]]*(174 + 36*Zeta[3]) + 
  g1^2*g3^2*M1^2*MatMul[Yb, Adj[Yb]]*(-96/5 + (128*Zeta[3])/3) + 
  g1^2*g3^2*M1*M3*MatMul[Yb, Adj[Yb]]*(-96/5 + (128*Zeta[3])/3) + 
  g1^2*g3^2*M3^2*MatMul[Yb, Adj[Yb]]*(-96/5 + (128*Zeta[3])/3) + 
  2*g2^2*g3^2*mb*MatMul[Yb, Adj[Yb]]*(-88 + 96*Zeta[3]) + 
  g2^2*g3^2*mq*Yb*Adj[Yb]*(-176 + 192*Zeta[3]) + 
  g2^2*g3^2*M2^2*MatMul[Yb, Adj[Yb]]*(-352 + 384*Zeta[3]) + 
  g2^2*g3^2*M2*M3*MatMul[Yb, Adj[Yb]]*(-352 + 384*Zeta[3]) + 
  g2^2*g3^2*M3^2*MatMul[Yb, Adj[Yb]]*(-352 + 384*Zeta[3]) + 
  g3^4*M3*MatMul[hb, Adj[Yb]]*(-64/3 + (2176*Zeta[3])/3) + 
  g3^6*M3^2*(20512/9 + 7680*Zeta[3]) + 
  g3^4*hb*((32*hbc)/3 - (1088*hbc*Zeta[3])/3) + 
  g2^2*g3^2*M2*Yb*(176*hbc - 192*hbc*Zeta[3]) + 
  g2^2*g3^2*M3*Yb*(176*hbc - 192*hbc*Zeta[3]) + 
  g2^2*MatMul[hb, Adj[Yb], Yb]*(18*hbc - 36*hbc*Zeta[3]) + 
  g2^2*MatMul[hb, Adj[Yt], Yt]*(18*hbc - 36*hbc*Zeta[3]) + 
  g2^2*MatMul[Yb, Adj[Yb], hb]*(18*hbc - 36*hbc*Zeta[3]) + 
  g2^2*MatMul[Yb, Adj[Yt], ht]*(18*hbc - 36*hbc*Zeta[3]) + 
  g1^2*g3^2*M1*Yb*((48*hbc)/5 - (64*hbc*Zeta[3])/3) + 
  g1^2*g3^2*M3*Yb*((48*hbc)/5 - (64*hbc*Zeta[3])/3) + 
  g2^4*hb*(-87*hbc - 18*hbc*Zeta[3]) + g1^2*g2^2*M1*Yb*
   ((86*hbc)/5 - (84*hbc*Zeta[3])/5) + g1^2*g2^2*M2*Yb*
   ((86*hbc)/5 - (84*hbc*Zeta[3])/5) + g1^2*M1*MatMul[Yb, Adj[Yt], Yt]*
   ((58*hbc)/15 - (36*hbc*Zeta[3])/5) + g1^2*M1*MatMul[Yb, Adj[Yb], Yb]*
   ((2*hbc)/3 - (12*hbc*Zeta[3])/5) + 
  g1^4*hb*((-337*hbc)/15 - (14*hbc*Zeta[3])/75) + 
  g1^4*M1*Yb*((674*hbc)/15 + (28*hbc*Zeta[3])/75) + 
  g1^2*MatMul[hb, Adj[Yb], Yb]*((-2*hbc)/3 + (12*hbc*Zeta[3])/5) + 
  g1^2*MatMul[Yb, Adj[Yb], hb]*((-2*hbc)/3 + (12*hbc*Zeta[3])/5) + 
  g1^2*MatMul[hb, Adj[Yt], Yt]*((-58*hbc)/15 + (36*hbc*Zeta[3])/5) + 
  g1^2*MatMul[Yb, Adj[Yt], ht]*((-58*hbc)/15 + (36*hbc*Zeta[3])/5) + 
  g1^2*g2^2*hb*((-86*hbc)/5 + (84*hbc*Zeta[3])/5) + 
  g1^2*g3^2*hb*((-48*hbc)/5 + (64*hbc*Zeta[3])/3) + 
  MatMul[hb, Adj[Yb], Yb, Adj[Yb], Yb]*(12*hbc + 24*hbc*Zeta[3]) + 
  MatMul[Yb, Adj[Yb], hb, Adj[Yb], Yb]*(12*hbc + 24*hbc*Zeta[3]) + 
  MatMul[Yb, Adj[Yb], Yb, Adj[Yb], hb]*(12*hbc + 24*hbc*Zeta[3]) + 
  g2^2*M2*MatMul[Yb, Adj[Yb], Yb]*(-18*hbc + 36*hbc*Zeta[3]) + 
  g2^2*M2*MatMul[Yb, Adj[Yt], Yt]*(-18*hbc + 36*hbc*Zeta[3]) + 
  g2^4*M2*Yb*(174*hbc + 36*hbc*Zeta[3]) + 
  g2^2*g3^2*hb*(-176*hbc + 192*hbc*Zeta[3]) + 
  g3^4*M3*Yb*((-64*hbc)/3 + (2176*hbc*Zeta[3])/3) + 
  g3^4*MatMul[Yb, Adj[Yb]]*((32*mh1)/3 - (1088*mh1*Zeta[3])/3) + 
  g2^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*(36*mh1 - 72*mh1*Zeta[3]) + 
  g2^4*MatMul[Yb, Adj[Yb]]*(-99*mh1 - 12*mh2 - 12*trace[ml] - 36*trace[mq] - 
    18*mh1*Zeta[3]) + g1^4*MatMul[Yb, Adj[Yb]]*
   ((-1721*mh1)/75 - (12*mh2)/25 - (8*trace[mb])/25 - (24*trace[me])/25 - 
    (12*trace[ml])/25 - (4*trace[mq])/25 - (32*trace[mt])/25 - 
    (14*mh1*Zeta[3])/75) + g1^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*
   ((-4*mh1)/3 + (24*mh1*Zeta[3])/5) + g1^2*g2^2*MatMul[Yb, Adj[Yb]]*
   ((-86*mh1)/5 + (84*mh1*Zeta[3])/5) + g1^2*g3^2*MatMul[Yb, Adj[Yb]]*
   ((-48*mh1)/5 + (64*mh1*Zeta[3])/3) + 
  MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb]]*(36*mh1 + 72*mh1*Zeta[3]) + 
  g2^2*g3^2*MatMul[Yb, Adj[Yb]]*(-176*mh1 + 192*mh1*Zeta[3]) + 
  g2^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*(18*mh1 + 18*mh2 - 36*mh1*Zeta[3] - 
    36*mh2*Zeta[3]) + g1^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*
   ((-58*mh1)/15 - (58*mh2)/15 + (36*mh1*Zeta[3])/5 + (36*mh2*Zeta[3])/5) + 
  g3^2*Yb*(-32*hbc*trace[hb, Adj[Yb]] + 32*hbc*trace[he, Adj[Ye]] + 
    192*hbc*trace[hb, Adj[Yb]]*Zeta[3]) + g3^2*M3*MatMul[Yb, Adj[Yb]]*
   (32*trace[hb, Adj[Yb]] + 32*trace[hbc, Yb] - 32*trace[he, Adj[Ye]] - 
    32*trace[hec, Ye] - 192*trace[hb, Adj[Yb]]*Zeta[3] - 
    192*trace[hbc, Yb]*Zeta[3]) + g3^2*MatMul[hb, Adj[Yb]]*
   (-32*trace[hbc, Yb] + 32*trace[hec, Ye] + 192*trace[hbc, Yb]*Zeta[3]) + 
  g2^2*Yb*(54*hbc*trace[hb, Adj[Yb]] + 18*hbc*trace[he, Adj[Ye]] - 
    108*hbc*trace[hb, Adj[Yb]]*Zeta[3] - 36*hbc*trace[he, Adj[Ye]]*Zeta[3]) + 
  g1^2*Yb*(14*hbc*trace[hb, Adj[Yb]] - 6*hbc*trace[he, Adj[Ye]] - 
    12*hbc*trace[hb, Adj[Yb]]*Zeta[3] + 12*hbc*trace[he, Adj[Ye]]*Zeta[3]) + 
  g2^2*MatMul[hb, Adj[Yb]]*(54*trace[hbc, Yb] + 18*trace[hec, Ye] - 
    108*trace[hbc, Yb]*Zeta[3] - 36*trace[hec, Ye]*Zeta[3]) + 
  g1^2*M1*MatMul[Yb, Adj[Yb]]*(-14*trace[hb, Adj[Yb]] - 14*trace[hbc, Yb] + 
    6*trace[he, Adj[Ye]] + 6*trace[hec, Ye] + 12*trace[hb, Adj[Yb]]*Zeta[3] + 
    12*trace[hbc, Yb]*Zeta[3] - 12*trace[he, Adj[Ye]]*Zeta[3] - 
    12*trace[hec, Ye]*Zeta[3]) + g1^2*MatMul[hb, Adj[Yb]]*
   (14*trace[hbc, Yb] - 6*trace[hec, Ye] - 12*trace[hbc, Yb]*Zeta[3] + 
    12*trace[hec, Ye]*Zeta[3]) + g2^2*M2*MatMul[Yb, Adj[Yb]]*
   (-54*trace[hb, Adj[Yb]] - 54*trace[hbc, Yb] - 18*trace[he, Adj[Ye]] - 
    18*trace[hec, Ye] + 108*trace[hb, Adj[Yb]]*Zeta[3] + 
    108*trace[hbc, Yb]*Zeta[3] + 36*trace[he, Adj[Ye]]*Zeta[3] + 
    36*trace[hec, Ye]*Zeta[3]) + g3^2*M3*MatMul[hb, Adj[Yb]]*
   (32*trace[Adj[Yb], Yb] - 32*trace[Adj[Ye], Ye] - 
    192*trace[Adj[Yb], Yb]*Zeta[3]) + 2*g3^2*mb*MatMul[Yb, Adj[Yb]]*
   (-16*trace[Adj[Yb], Yb] + 16*trace[Adj[Ye], Ye] + 
    96*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*mq*Yb*Adj[Yb]*
   (-32*trace[Adj[Yb], Yb] + 32*trace[Adj[Ye], Ye] + 
    192*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*M3^2*MatMul[Yb, Adj[Yb]]*
   (-64*trace[Adj[Yb], Yb] + 64*trace[Adj[Ye], Ye] + 
    384*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g3^2*M3*Yb*(32*hbc*trace[Adj[Yb], Yb] - 32*hbc*trace[Adj[Ye], Ye] - 
    192*hbc*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g3^2*hb*(-32*hbc*trace[Adj[Yb], Yb] + 32*hbc*trace[Adj[Ye], Ye] + 
    192*hbc*trace[Adj[Yb], Yb]*Zeta[3]) + g2^2*M2^2*MatMul[Yb, Adj[Yb]]*
   (108*trace[Adj[Yb], Yb] + 36*trace[Adj[Ye], Ye] - 
    216*trace[Adj[Yb], Yb]*Zeta[3] - 72*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g2^2*mq*Yb*Adj[Yb]*(54*trace[Adj[Yb], Yb] + 18*trace[Adj[Ye], Ye] - 
    108*trace[Adj[Yb], Yb]*Zeta[3] - 36*trace[Adj[Ye], Ye]*Zeta[3]) + 
  2*g2^2*mb*MatMul[Yb, Adj[Yb]]*(27*trace[Adj[Yb], Yb] + 
    9*trace[Adj[Ye], Ye] - 54*trace[Adj[Yb], Yb]*Zeta[3] - 
    18*trace[Adj[Ye], Ye]*Zeta[3]) + g1^2*M1*MatMul[hb, Adj[Yb]]*
   (-14*trace[Adj[Yb], Yb] + 6*trace[Adj[Ye], Ye] + 
    12*trace[Adj[Yb], Yb]*Zeta[3] - 12*trace[Adj[Ye], Ye]*Zeta[3]) + 
  2*g1^2*mb*MatMul[Yb, Adj[Yb]]*(7*trace[Adj[Yb], Yb] - 
    3*trace[Adj[Ye], Ye] - 6*trace[Adj[Yb], Yb]*Zeta[3] + 
    6*trace[Adj[Ye], Ye]*Zeta[3]) + g1^2*mq*Yb*Adj[Yb]*
   (14*trace[Adj[Yb], Yb] - 6*trace[Adj[Ye], Ye] - 
    12*trace[Adj[Yb], Yb]*Zeta[3] + 12*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g1^2*M1^2*MatMul[Yb, Adj[Yb]]*(28*trace[Adj[Yb], Yb] - 
    12*trace[Adj[Ye], Ye] - 24*trace[Adj[Yb], Yb]*Zeta[3] + 
    24*trace[Adj[Ye], Ye]*Zeta[3]) + g2^2*M2*MatMul[hb, Adj[Yb]]*
   (-54*trace[Adj[Yb], Yb] - 18*trace[Adj[Ye], Ye] + 
    108*trace[Adj[Yb], Yb]*Zeta[3] + 36*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g2^2*hb*(54*hbc*trace[Adj[Yb], Yb] + 18*hbc*trace[Adj[Ye], Ye] - 
    108*hbc*trace[Adj[Yb], Yb]*Zeta[3] - 36*hbc*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g1^2*M1*Yb*(-14*hbc*trace[Adj[Yb], Yb] + 6*hbc*trace[Adj[Ye], Ye] + 
    12*hbc*trace[Adj[Yb], Yb]*Zeta[3] - 12*hbc*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g1^2*hb*(14*hbc*trace[Adj[Yb], Yb] - 6*hbc*trace[Adj[Ye], Ye] - 
    12*hbc*trace[Adj[Yb], Yb]*Zeta[3] + 12*hbc*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g2^2*M2*Yb*(-54*hbc*trace[Adj[Yb], Yb] - 18*hbc*trace[Adj[Ye], Ye] + 
    108*hbc*trace[Adj[Yb], Yb]*Zeta[3] + 36*hbc*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g3^2*MatMul[Yb, Adj[Yb]]*(-32*trace[hb, hbc] + 32*trace[he, hec] - 
    64*mh1*trace[Adj[Yb], Yb] + 64*mh1*trace[Adj[Ye], Ye] - 
    32*trace[Yb, Adj[Yb], mb] + 32*trace[Ye, Adj[Ye], me] - 
    32*trace[Adj[Yb], Yb, mq] + 32*trace[Adj[Ye], Ye, ml] + 
    192*trace[hb, hbc]*Zeta[3] + 384*mh1*trace[Adj[Yb], Yb]*Zeta[3] + 
    192*trace[Yb, Adj[Yb], mb]*Zeta[3] + 192*trace[Adj[Yb], Yb, mq]*
     Zeta[3]) + g2^2*MatMul[Yb, Adj[Yb]]*(54*trace[hb, hbc] + 
    18*trace[he, hec] + 108*mh1*trace[Adj[Yb], Yb] + 
    36*mh1*trace[Adj[Ye], Ye] + 54*trace[Yb, Adj[Yb], mb] + 
    18*trace[Ye, Adj[Ye], me] + 54*trace[Adj[Yb], Yb, mq] + 
    18*trace[Adj[Ye], Ye, ml] - 108*trace[hb, hbc]*Zeta[3] - 
    36*trace[he, hec]*Zeta[3] - 216*mh1*trace[Adj[Yb], Yb]*Zeta[3] - 
    72*mh1*trace[Adj[Ye], Ye]*Zeta[3] - 108*trace[Yb, Adj[Yb], mb]*Zeta[3] - 
    36*trace[Ye, Adj[Ye], me]*Zeta[3] - 108*trace[Adj[Yb], Yb, mq]*Zeta[3] - 
    36*trace[Adj[Ye], Ye, ml]*Zeta[3]) + g1^2*MatMul[Yb, Adj[Yb]]*
   (14*trace[hb, hbc] - 6*trace[he, hec] + 28*mh1*trace[Adj[Yb], Yb] - 
    12*mh1*trace[Adj[Ye], Ye] + 14*trace[Yb, Adj[Yb], mb] - 
    6*trace[Ye, Adj[Ye], me] + 14*trace[Adj[Yb], Yb, mq] - 
    6*trace[Adj[Ye], Ye, ml] - 12*trace[hb, hbc]*Zeta[3] + 
    12*trace[he, hec]*Zeta[3] - 24*mh1*trace[Adj[Yb], Yb]*Zeta[3] + 
    24*mh1*trace[Adj[Ye], Ye]*Zeta[3] - 12*trace[Yb, Adj[Yb], mb]*Zeta[3] + 
    12*trace[Ye, Adj[Ye], me]*Zeta[3] - 12*trace[Adj[Yb], Yb, mq]*Zeta[3] + 
    12*trace[Adj[Ye], Ye, ml]*Zeta[3])}

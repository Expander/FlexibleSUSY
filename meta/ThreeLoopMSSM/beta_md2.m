{4*hb*hbc - (8*g1^2*M1^2)/15 - (32*g3^2*M3^2)/3 + 4*mq*Yb*Ybc + 
  4*mb*MatMul[Yb, Ybc] + 4*mh1*MatMul[Yb, Ybc], 
 (4*g1^2*hb*hbc)/5 + 12*g2^2*hb*hbc + (808*g1^4*M1^2)/75 + 
  (128*g1^2*g3^2*M1^2)/45 + (128*g1^2*g3^2*M1*M3)/45 + 
  (128*g1^2*g3^2*M3^2)/45 - (128*g3^4*M3^2)/3 - (4*g1^2*hbc*M1*Yb)/5 - 
  12*g2^2*hbc*M2*Yb + (4*g1^2*mq*Yb*Ybc)/5 + 12*g2^2*mq*Yb*Ybc - 
  (4*g1^2*M1*MatMul[hb, Ybc])/5 - 12*g2^2*M2*MatMul[hb, Ybc] + 
  (8*g1^2*M1^2*MatMul[Yb, Ybc])/5 + 24*g2^2*M2^2*MatMul[Yb, Ybc] + 
  (4*g1^2*mb*MatMul[Yb, Ybc])/5 + 12*g2^2*mb*MatMul[Yb, Ybc] + 
  (4*g1^2*mh1*MatMul[Yb, Ybc])/5 + 12*g2^2*mh1*MatMul[Yb, Ybc] - 
  4*mb*MatMul[Yb, Ybc]^2 - 4*mt*MatMul[Yb, Ytc]*MatMul[Yt, Ybc] - 
  4*hbc*MatMul[hb, Ybc, Yb] - 4*hbc*MatMul[hb, Ytc, Yt] - 
  4*hbc*MatMul[Yb, Ybc, hb] - 4*mq*Ybc*MatMul[Yb, Ybc, Yb] - 
  4*hbc*MatMul[Yb, Ytc, ht] - 4*mq*Ybc*MatMul[Yb, Ytc, Yt] - 
  4*mq*Yb*MatMul[Ybc, Yb, Ybc] - 4*mq*Yb*MatMul[Ytc, Yt, Ybc] - 
  4*MatMul[hb, hbc, Yb, Ybc] - 4*MatMul[hb, htc, Yt, Ybc] - 
  4*MatMul[Yb, hbc, hb, Ybc] - 4*MatMul[Yb, htc, ht, Ybc] - 
  4*mb*MatMul[Yb, Ybc, Yb, Ybc] - 8*mh1*MatMul[Yb, Ybc, Yb, Ybc] - 
  4*mb*MatMul[Yb, Ytc, Yt, Ybc] + (-4*mh1 - 4*mh2)*MatMul[Yb, Ytc, Yt, Ybc] + 
  g1^4*((4*mh1)/25 + (4*mh2)/25 + (8*trace[mb])/75 + (8*trace[me])/25 + 
    (4*trace[ml])/25 + (4*trace[mq])/75 + (32*trace[mt])/75) + 
  g3^4*((16*trace[mb])/3 + (32*trace[mq])/3 + (16*trace[mt])/3) + 
  Yb*(-12*hbc*trace[hb, Ybc] - 4*hbc*trace[he, Adj[Ye]]) + 
  MatMul[hb, Ybc]*(-12*trace[hbc, Yb] - 4*trace[hec, Ye]) + 
  mq*Yb*Ybc*(-12*trace[Ybc, Yb] - 4*trace[Adj[Ye], Ye]) + 
  2*mb*MatMul[Yb, Ybc]*(-6*trace[Ybc, Yb] - 2*trace[Adj[Ye], Ye]) + 
  hb*(-12*hbc*trace[Ybc, Yb] - 4*hbc*trace[Adj[Ye], Ye]) + 
  MatMul[Yb, Ybc]*(-12*trace[hb, hbc] - 4*trace[he, hec] - 
    24*mh1*trace[Ybc, Yb] - 8*mh1*trace[Adj[Ye], Ye] - 
    12*trace[Yb, Ybc, mb] - 12*trace[Ybc, Yb, mq] - 
    4*trace[Ye, Adj[Ye], me] - 4*trace[Adj[Ye], Ye, ml]), 
 (128*g3^2*mb*MatMul[Yb, Ybc]^2)/3 + 
  (128*g3^2*mt*MatMul[Yb, Ytc]*MatMul[Yt, Ybc])/3 + 
  (128*g3^2*hbc*MatMul[hb, Ybc, Yb])/3 + (128*g3^2*hbc*MatMul[hb, Ytc, Yt])/
   3 + (128*g3^2*hbc*MatMul[Yb, Ybc, hb])/3 - 
  (128*g3^2*hbc*M3*MatMul[Yb, Ybc, Yb])/3 + 
  (128*g3^2*mq*Ybc*MatMul[Yb, Ybc, Yb])/3 + 
  (128*g3^2*hbc*MatMul[Yb, Ytc, ht])/3 - 
  (128*g3^2*hbc*M3*MatMul[Yb, Ytc, Yt])/3 + 
  (128*g3^2*mq*Ybc*MatMul[Yb, Ytc, Yt])/3 + 
  (128*g3^2*mq*Yb*MatMul[Ybc, Yb, Ybc])/3 - 4*mq*MatMul[Yb, Ytc, Yt]*
   MatMul[Ybc, Yb, Ybc] + (128*g3^2*mq*Yb*MatMul[Ytc, Yt, Ybc])/3 - 
  4*mq*MatMul[Yb, Ybc, Yb]*MatMul[Ytc, Yt, Ybc] + 
  12*mq*MatMul[Yb, Ytc, Yt]*MatMul[Ytc, Yt, Ybc] + 
  (128*g3^2*MatMul[hb, hbc, Yb, Ybc])/3 + (128*g3^2*MatMul[hb, htc, Yt, Ybc])/
   3 - (128*g3^2*M3*MatMul[hb, Ybc, Yb, Ybc])/3 - 
  (128*g3^2*M3*MatMul[hb, Ytc, Yt, Ybc])/3 + 
  (128*g3^2*MatMul[Yb, hbc, hb, Ybc])/3 - 
  (128*g3^2*M3*MatMul[Yb, hbc, Yb, Ybc])/3 + 
  (128*g3^2*MatMul[Yb, htc, ht, Ybc])/3 - 
  (128*g3^2*M3*MatMul[Yb, htc, Yt, Ybc])/3 - 
  (128*g3^2*M3*MatMul[Yb, Ybc, hb, Ybc])/3 + 
  (256*g3^2*M3^2*MatMul[Yb, Ybc, Yb, Ybc])/3 + 
  (128*g3^2*mb*MatMul[Yb, Ybc, Yb, Ybc])/3 + 
  (256*g3^2*mh1*MatMul[Yb, Ybc, Yb, Ybc])/3 - 
  4*mt*MatMul[Yt, Ybc]*MatMul[Yb, Ybc, Yb, Ytc] - 
  (128*g3^2*M3*MatMul[Yb, Ytc, ht, Ybc])/3 + 
  (256*g3^2*M3^2*MatMul[Yb, Ytc, Yt, Ybc])/3 + 
  (128*g3^2*mb*MatMul[Yb, Ytc, Yt, Ybc])/3 + g3^2*((128*mh1)/3 + (128*mh2)/3)*
   MatMul[Yb, Ytc, Yt, Ybc] - 8*mb*MatMul[Yb, Ybc]*MatMul[Yb, Ytc, Yt, Ybc] + 
  12*mt*MatMul[Yt, Ybc]*MatMul[Yb, Ytc, Yt, Ytc] - 
  4*mt*MatMul[Yb, Ytc]*MatMul[Yt, Ybc, Yb, Ybc] + 
  12*mt*MatMul[Yb, Ytc]*MatMul[Yt, Ytc, Yt, Ybc] - 
  4*hbc*MatMul[hb, Ybc, Yb, Ytc, Yt] - 4*hbc*MatMul[hb, Ytc, Yt, Ybc, Yb] + 
  12*hbc*MatMul[hb, Ytc, Yt, Ytc, Yt] - 4*hbc*MatMul[Yb, Ybc, hb, Ytc, Yt] - 
  4*hbc*MatMul[Yb, Ybc, Yb, Ytc, ht] - 
  4*mq*Ybc*MatMul[Yb, Ybc, Yb, Ytc, Yt] - 
  4*hbc*MatMul[Yb, Ytc, ht, Ybc, Yb] + 12*hbc*MatMul[Yb, Ytc, ht, Ytc, Yt] - 
  4*hbc*MatMul[Yb, Ytc, Yt, Ybc, hb] - 
  4*mq*Ybc*MatMul[Yb, Ytc, Yt, Ybc, Yb] + 
  12*hbc*MatMul[Yb, Ytc, Yt, Ytc, ht] + 
  12*mq*Ybc*MatMul[Yb, Ytc, Yt, Ytc, Yt] - 
  4*mq*Yb*MatMul[Ybc, Yb, Ytc, Yt, Ybc] - 
  4*mq*Yb*MatMul[Ytc, Yt, Ybc, Yb, Ybc] + 
  12*mq*Yb*MatMul[Ytc, Yt, Ytc, Yt, Ybc] - 
  4*MatMul[hb, hbc, Yb, Ytc, Yt, Ybc] - 4*MatMul[hb, htc, Yt, Ybc, Yb, Ybc] + 
  12*MatMul[hb, htc, Yt, Ytc, Yt, Ybc] - 
  4*MatMul[hb, Ybc, Yb, htc, Yt, Ybc] - 4*MatMul[hb, Ytc, Yt, hbc, Yb, Ybc] + 
  12*MatMul[hb, Ytc, Yt, htc, Yt, Ybc] - 
  4*MatMul[Yb, hbc, hb, Ytc, Yt, Ybc] - 4*MatMul[Yb, hbc, Yb, Ytc, ht, Ybc] - 
  4*MatMul[Yb, htc, ht, Ybc, Yb, Ybc] + 
  12*MatMul[Yb, htc, ht, Ytc, Yt, Ybc] - 
  4*MatMul[Yb, htc, Yt, Ybc, hb, Ybc] + 
  12*MatMul[Yb, htc, Yt, Ytc, ht, Ybc] - 
  4*MatMul[Yb, Ybc, hb, htc, Yt, Ybc] - 4*MatMul[Yb, Ybc, Yb, htc, ht, Ybc] - 
  4*mb*MatMul[Yb, Ybc, Yb, Ytc, Yt, Ybc] + 
  (-8*mh1 - 4*mh2)*MatMul[Yb, Ybc, Yb, Ytc, Yt, Ybc] - 
  4*MatMul[Yb, Ytc, ht, hbc, Yb, Ybc] + 
  12*MatMul[Yb, Ytc, ht, htc, Yt, Ybc] - 
  4*MatMul[Yb, Ytc, Yt, hbc, hb, Ybc] + 
  12*MatMul[Yb, Ytc, Yt, htc, ht, Ybc] - 
  4*mb*MatMul[Yb, Ytc, Yt, Ybc, Yb, Ybc] + 
  (-8*mh1 - 4*mh2)*MatMul[Yb, Ytc, Yt, Ybc, Yb, Ybc] + 
  12*mb*MatMul[Yb, Ytc, Yt, Ytc, Yt, Ybc] + 
  (12*mh1 + 24*mh2)*MatMul[Yb, Ytc, Yt, Ytc, Yt, Ybc] + 
  g1^6*((-808*mh1)/375 - (808*mh2)/375 - (1616*trace[mb])/1125 - 
    (1616*trace[me])/375 - (808*trace[ml])/375 - (808*trace[mq])/1125 - 
    (6464*trace[mt])/1125) + g1^4*g3^2*((-64*mh1)/75 - (64*mh2)/75 - 
    (128*trace[mb])/225 - (128*trace[me])/75 - (64*trace[ml])/75 - 
    (64*trace[mq])/225 - (512*trace[mt])/225) + 
  g1^2*g3^4*((-64*trace[mb])/45 - (128*trace[mq])/45 - (64*trace[mt])/45) + 
  g3^6*((320*trace[mb])/9 + (640*trace[mq])/9 + (320*trace[mt])/9) + 
  MatMul[Yb, hbc, Yb, Ybc]*(12*trace[hb, Ybc] + 4*trace[he, Adj[Ye]]) + 
  MatMul[Yb, Ybc, Yb]*(12*hbc*trace[hb, Ybc] + 4*hbc*trace[he, Adj[Ye]]) + 
  MatMul[hb, Ybc, Yb, Ybc]*(12*trace[hbc, Yb] + 4*trace[hec, Ye]) + 
  MatMul[Yb, Ybc, hb, Ybc]*(12*trace[hbc, Yb] + 4*trace[hec, Ye]) + 
  MatMul[Yb, htc, Yt, Ybc]*(-12*trace[hb, Ybc] - 4*trace[he, Adj[Ye]] + 
    24*trace[ht, Ytc]) + MatMul[Yb, Ytc, Yt]*(-12*hbc*trace[hb, Ybc] - 
    4*hbc*trace[he, Adj[Ye]] + 24*hbc*trace[ht, Ytc]) + 
  g1^4*M1*((56*trace[hb, Ybc])/15 + (56*trace[hbc, Yb])/15 + 
    (24*trace[he, Adj[Ye]])/5 + (24*trace[hec, Ye])/5 + 
    (104*trace[ht, Ytc])/15 + (104*trace[htc, Yt])/15) + 
  MatMul[hb, Ytc, Yt, Ybc]*(-12*trace[hbc, Yb] - 4*trace[hec, Ye] + 
    24*trace[htc, Yt]) + MatMul[Yb, Ytc, ht, Ybc]*
   (-12*trace[hbc, Yb] - 4*trace[hec, Ye] + 24*trace[htc, Yt]) + 
  g3^4*M3*((320*trace[hb, Ybc])/3 + (320*trace[hbc, Yb])/3 + 
    (320*trace[ht, Ytc])/3 + (320*trace[htc, Yt])/3) + 
  g3^4*M3^2*(-320*trace[Ybc, Yb] - 320*trace[Ytc, Yt]) + 
  g1^4*M1^2*((-56*trace[Ybc, Yb])/5 - (104*trace[Ytc, Yt])/5 - 
    (72*trace[Adj[Ye], Ye])/5) + mt*MatMul[Yb, Ytc]*MatMul[Yt, Ybc]*
   (-12*trace[Ybc, Yb] + 24*trace[Ytc, Yt] - 4*trace[Adj[Ye], Ye]) + 
  mq*Ybc*MatMul[Yb, Ytc, Yt]*(-12*trace[Ybc, Yb] + 24*trace[Ytc, Yt] - 
    4*trace[Adj[Ye], Ye]) + mq*Yb*MatMul[Ytc, Yt, Ybc]*
   (-12*trace[Ybc, Yb] + 24*trace[Ytc, Yt] - 4*trace[Adj[Ye], Ye]) + 
  MatMul[hb, htc, Yt, Ybc]*(-12*trace[Ybc, Yb] + 24*trace[Ytc, Yt] - 
    4*trace[Adj[Ye], Ye]) + MatMul[Yb, htc, ht, Ybc]*
   (-12*trace[Ybc, Yb] + 24*trace[Ytc, Yt] - 4*trace[Adj[Ye], Ye]) + 
  2*mb*MatMul[Yb, Ytc, Yt, Ybc]*(-6*trace[Ybc, Yb] + 12*trace[Ytc, Yt] - 
    2*trace[Adj[Ye], Ye]) + 2*mb*MatMul[Yb, Ybc, Yb, Ybc]*
   (6*trace[Ybc, Yb] + 2*trace[Adj[Ye], Ye]) + 
  mb*MatMul[Yb, Ybc]^2*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  mq*Ybc*MatMul[Yb, Ybc, Yb]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  mq*Yb*MatMul[Ybc, Yb, Ybc]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[hb, hbc, Yb, Ybc]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[Yb, hbc, hb, Ybc]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[hb, Ytc, Yt]*(-12*hbc*trace[Ybc, Yb] + 24*hbc*trace[Ytc, Yt] - 
    4*hbc*trace[Adj[Ye], Ye]) + MatMul[Yb, Ytc, ht]*
   (-12*hbc*trace[Ybc, Yb] + 24*hbc*trace[Ytc, Yt] - 
    4*hbc*trace[Adj[Ye], Ye]) + MatMul[hb, Ybc, Yb]*
   (12*hbc*trace[Ybc, Yb] + 4*hbc*trace[Adj[Ye], Ye]) + 
  MatMul[Yb, Ybc, hb]*(12*hbc*trace[Ybc, Yb] + 4*hbc*trace[Adj[Ye], Ye]) + 
  g3^4*(-64*trace[hb, hbc] - 64*trace[ht, htc] - 64*mh1*trace[Ybc, Yb] - 
    64*mh2*trace[Ytc, Yt] - 64*trace[Yb, Ybc, mb] - 64*trace[Ybc, Yb, mq] - 
    64*trace[Yt, Ytc, mt] - 64*trace[Ytc, Yt, mq]) + 
  MatMul[Yb, Ytc, Yt, Ybc]*(-12*trace[hb, hbc] - 4*trace[he, hec] + 
    24*trace[ht, htc] - 24*mh1*trace[Ybc, Yb] - 12*mh2*trace[Ybc, Yb] + 
    24*mh1*trace[Ytc, Yt] + 48*mh2*trace[Ytc, Yt] - 
    8*mh1*trace[Adj[Ye], Ye] - 4*mh2*trace[Adj[Ye], Ye] - 
    12*trace[Yb, Ybc, mb] - 12*trace[Ybc, Yb, mq] - 
    4*trace[Ye, Adj[Ye], me] + 24*trace[Yt, Ytc, mt] + 
    24*trace[Ytc, Yt, mq] - 4*trace[Adj[Ye], Ye, ml]) + 
  g1^4*((-56*trace[hb, hbc])/25 - (72*trace[he, hec])/25 - 
    (104*trace[ht, htc])/25 - (56*mh1*trace[Ybc, Yb])/25 - 
    (104*mh2*trace[Ytc, Yt])/25 - (72*mh1*trace[Adj[Ye], Ye])/25 - 
    (56*trace[Yb, Ybc, mb])/25 - (56*trace[Ybc, Yb, mq])/25 - 
    (72*trace[Ye, Adj[Ye], me])/25 - (104*trace[Yt, Ytc, mt])/25 - 
    (104*trace[Ytc, Yt, mq])/25 - (72*trace[Adj[Ye], Ye, ml])/25) + 
  MatMul[Yb, Ybc, Yb, Ybc]*(12*trace[hb, hbc] + 4*trace[he, hec] + 
    36*mh1*trace[Ybc, Yb] + 12*mh1*trace[Adj[Ye], Ye] + 
    12*trace[Yb, Ybc, mb] + 12*trace[Ybc, Yb, mq] + 
    4*trace[Ye, Adj[Ye], me] + 4*trace[Adj[Ye], Ye, ml]) + 
  Yb*(-72*hbc*trace[hb, Ybc]*trace[Ybc, Yb] - 24*hbc*trace[he, Adj[Ye]]*
     trace[Ybc, Yb] - 24*hbc*trace[hb, Ybc]*trace[Adj[Ye], Ye] - 
    8*hbc*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] + 
    24*hbc*trace[hb, Ytc, Yt, Ybc] + 144*hbc*trace[Ybc, hb, Ybc, Yb] + 
    24*hbc*trace[Ytc, ht, Ybc, Yb] + 48*hbc*trace[Adj[Ye], he, Adj[Ye], 
      Ye]) + MatMul[hb, Ybc]*(-72*trace[hbc, Yb]*trace[Ybc, Yb] - 
    24*trace[hec, Ye]*trace[Ybc, Yb] - 24*trace[hbc, Yb]*trace[Adj[Ye], Ye] - 
    8*trace[hec, Ye]*trace[Adj[Ye], Ye] + 24*trace[htc, Yt, Ybc, Yb] + 
    144*trace[Ybc, Yb, hbc, Yb] + 24*trace[Ytc, Yt, hbc, Yb] + 
    48*trace[Adj[Ye], Ye, hec, Ye]) + 2*mb*MatMul[Yb, Ybc]*
   (-18*trace[Ybc, Yb]^2 - 12*trace[Ybc, Yb]*trace[Adj[Ye], Ye] - 
    2*trace[Adj[Ye], Ye]^2 + 36*trace[Ybc, Yb, Ybc, Yb] + 
    12*trace[Ytc, Yt, Ybc, Yb] + 12*trace[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  mq*Yb*Ybc*(-36*trace[Ybc, Yb]^2 - 24*trace[Ybc, Yb]*trace[Adj[Ye], Ye] - 
    4*trace[Adj[Ye], Ye]^2 + 72*trace[Ybc, Yb, Ybc, Yb] + 
    24*trace[Ytc, Yt, Ybc, Yb] + 24*trace[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  hb*(-36*hbc*trace[Ybc, Yb]^2 - 24*hbc*trace[Ybc, Yb]*trace[Adj[Ye], Ye] - 
    4*hbc*trace[Adj[Ye], Ye]^2 + 72*hbc*trace[Ybc, Yb, Ybc, Yb] + 
    24*hbc*trace[Ytc, Yt, Ybc, Yb] + 24*hbc*trace[Adj[Ye], Ye, Adj[Ye], 
      Ye]) + MatMul[Yb, Ybc]*(-72*trace[hb, Ybc]*trace[hbc, Yb] - 
    24*trace[hbc, Yb]*trace[he, Adj[Ye]] - 24*trace[hb, Ybc]*trace[hec, Ye] - 
    8*trace[he, Adj[Ye]]*trace[hec, Ye] - 72*trace[hb, hbc]*trace[Ybc, Yb] - 
    24*trace[he, hec]*trace[Ybc, Yb] - 108*mh1*trace[Ybc, Yb]^2 - 
    24*trace[hb, hbc]*trace[Adj[Ye], Ye] - 8*trace[he, hec]*
     trace[Adj[Ye], Ye] - 72*mh1*trace[Ybc, Yb]*trace[Adj[Ye], Ye] - 
    12*mh1*trace[Adj[Ye], Ye]^2 - 72*trace[Ybc, Yb]*trace[Yb, Ybc, mb] - 
    24*trace[Adj[Ye], Ye]*trace[Yb, Ybc, mb] - 72*trace[Ybc, Yb]*
     trace[Ybc, Yb, mq] - 24*trace[Adj[Ye], Ye]*trace[Ybc, Yb, mq] - 
    24*trace[Ybc, Yb]*trace[Ye, Adj[Ye], me] - 8*trace[Adj[Ye], Ye]*
     trace[Ye, Adj[Ye], me] - 24*trace[Ybc, Yb]*trace[Adj[Ye], Ye, ml] - 
    8*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, ml] + 
    24*trace[hb, htc, Yt, Ybc] + 144*trace[hbc, hb, Ybc, Yb] + 
    24*trace[hbc, hb, Ytc, Yt] + 48*trace[hec, he, Adj[Ye], Ye] + 
    24*trace[htc, ht, Ybc, Yb] + 144*trace[Ybc, hb, hbc, Yb] + 
    216*mh1*trace[Ybc, Yb, Ybc, Yb] + 24*trace[Ytc, ht, hbc, Yb] + 
    48*mh1*trace[Ytc, Yt, Ybc, Yb] + 24*mh2*trace[Ytc, Yt, Ybc, Yb] + 
    48*trace[Adj[Ye], he, hec, Ye] + 72*mh1*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    144*trace[Yb, Ybc, Yb, Ybc, mb] + 24*trace[Yb, Ytc, Yt, Ybc, mb] + 
    144*trace[Ybc, Yb, Ybc, Yb, mq] + 24*trace[Ybc, Yb, Ytc, Yt, mq] + 
    48*trace[Ye, Adj[Ye], Ye, Adj[Ye], me] + 24*trace[Yt, Ybc, Yb, Ytc, mt] + 
    24*trace[Ytc, Yt, Ybc, Yb, mq] + 48*trace[Adj[Ye], Ye, Adj[Ye], Ye, 
      ml]) + g3^4*M3^2*MatMul[Yb, Ybc]*(64 - 2176*Zeta[3]) + 
  g2^2*g3^4*M3^2*(720 - 864*Zeta[3]) + g2^2*g3^4*M2*M3*(480 - 576*Zeta[3]) + 
  g3^4*mq*Yb*Ybc*(32/3 - (1088*Zeta[3])/3) + 
  g2^2*g3^4*M2^2*(288 - 288*Zeta[3]) + g1^2*g3^4*M3^2*
   (1936/15 - (1056*Zeta[3])/5) + g2^2*g3^2*M2*MatMul[hb, Ybc]*
   (176 - 192*Zeta[3]) + g2^2*g3^2*M3*MatMul[hb, Ybc]*(176 - 192*Zeta[3]) + 
  2*g3^4*mb*MatMul[Yb, Ybc]*(16/3 - (544*Zeta[3])/3) + 
  g1^2*g3^4*M1*M3*(3616/45 - (704*Zeta[3])/5) + 
  g2^4*M2^2*MatMul[Yb, Ybc]*(-474 - 108*Zeta[3]) + 
  g1^4*g3^2*M1^2*(2912/75 - (2112*Zeta[3])/25) + 
  g1^6*M1^2*(227548/1125 - (9552*Zeta[3])/125) + 
  g2^2*M2^2*MatMul[Yb, Ybc, Yb, Ybc]*(36 - 72*Zeta[3]) + 
  g2^2*M2^2*MatMul[Yb, Ytc, Yt, Ybc]*(36 - 72*Zeta[3]) + 
  g1^2*g3^4*M1^2*(2336/45 - (352*Zeta[3])/5) + 
  g1^4*g3^2*M1*M3*(5824/225 - (1408*Zeta[3])/25) + 
  g2^2*mb*MatMul[Yb, Ybc]^2*(18 - 36*Zeta[3]) + 
  g2^2*mt*MatMul[Yb, Ytc]*MatMul[Yt, Ybc]*(18 - 36*Zeta[3]) + 
  g2^2*mq*Ybc*MatMul[Yb, Ybc, Yb]*(18 - 36*Zeta[3]) + 
  g2^2*mq*Ybc*MatMul[Yb, Ytc, Yt]*(18 - 36*Zeta[3]) + 
  g2^2*mq*Yb*MatMul[Ybc, Yb, Ybc]*(18 - 36*Zeta[3]) + 
  g2^2*mq*Yb*MatMul[Ytc, Yt, Ybc]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[hb, hbc, Yb, Ybc]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[hb, htc, Yt, Ybc]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[Yb, hbc, hb, Ybc]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[Yb, htc, ht, Ybc]*(18 - 36*Zeta[3]) + 
  g1^4*g3^2*M3^2*(3968/225 - (704*Zeta[3])/25) + 
  g1^4*g2^2*M1^2*(108/5 - (648*Zeta[3])/25) + g1^2*g3^2*M1*MatMul[hb, Ybc]*
   (48/5 - (64*Zeta[3])/3) + g1^2*g3^2*M3*MatMul[hb, Ybc]*
   (48/5 - (64*Zeta[3])/3) + g2^4*mq*Yb*Ybc*(-87 - 18*Zeta[3]) + 
  2*g2^2*mb*MatMul[Yb, Ybc, Yb, Ybc]*(9 - 18*Zeta[3]) + 
  2*g2^2*mb*MatMul[Yb, Ytc, Yt, Ybc]*(9 - 18*Zeta[3]) + 
  g1^4*g2^2*M1*M2*(72/5 - (432*Zeta[3])/25) + g1^2*g2^2*M1*MatMul[hb, Ybc]*
   (86/5 - (84*Zeta[3])/5) + g1^2*g2^2*M2*MatMul[hb, Ybc]*
   (86/5 - (84*Zeta[3])/5) + 2*g2^4*mb*MatMul[Yb, Ybc]*(-87/2 - 9*Zeta[3]) + 
  g1^4*g2^2*M2^2*(216/25 - (216*Zeta[3])/25) + 
  g1^2*M1*MatMul[hb, Ytc, Yt, Ybc]*(58/15 - (36*Zeta[3])/5) + 
  g1^2*M1*MatMul[Yb, htc, Yt, Ybc]*(58/15 - (36*Zeta[3])/5) + 
  g1^2*M1*MatMul[Yb, Ytc, ht, Ybc]*(58/15 - (36*Zeta[3])/5) + 
  g1^2*M1*MatMul[hb, Ybc, Yb, Ybc]*(2/3 - (12*Zeta[3])/5) + 
  g1^2*M1*MatMul[Yb, hbc, Yb, Ybc]*(2/3 - (12*Zeta[3])/5) + 
  g1^2*M1*MatMul[Yb, Ybc, hb, Ybc]*(2/3 - (12*Zeta[3])/5) + 
  g1^4*M1^2*MatMul[Yb, Ybc]*(-674/5 - (28*Zeta[3])/25) + 
  g1^4*mq*Yb*Ybc*(-337/15 - (14*Zeta[3])/75) + 2*g1^4*mb*MatMul[Yb, Ybc]*
   (-337/30 - (7*Zeta[3])/75) + g1^4*M1*MatMul[hb, Ybc]*
   (674/15 + (28*Zeta[3])/75) + 2*g1^2*mb*MatMul[Yb, Ybc, Yb, Ybc]*
   (-1/3 + (6*Zeta[3])/5) + g1^2*mb*MatMul[Yb, Ybc]^2*
   (-2/3 + (12*Zeta[3])/5) + g1^2*mq*Ybc*MatMul[Yb, Ybc, Yb]*
   (-2/3 + (12*Zeta[3])/5) + g1^2*mq*Yb*MatMul[Ybc, Yb, Ybc]*
   (-2/3 + (12*Zeta[3])/5) + g1^2*MatMul[hb, hbc, Yb, Ybc]*
   (-2/3 + (12*Zeta[3])/5) + g1^2*MatMul[Yb, hbc, hb, Ybc]*
   (-2/3 + (12*Zeta[3])/5) + 2*g1^2*mb*MatMul[Yb, Ytc, Yt, Ybc]*
   (-29/15 + (18*Zeta[3])/5) + g1^2*M1^2*MatMul[Yb, Ybc, Yb, Ybc]*
   (-4/3 + (24*Zeta[3])/5) + g1^2*mt*MatMul[Yb, Ytc]*MatMul[Yt, Ybc]*
   (-58/15 + (36*Zeta[3])/5) + g1^2*mq*Ybc*MatMul[Yb, Ytc, Yt]*
   (-58/15 + (36*Zeta[3])/5) + g1^2*mq*Yb*MatMul[Ytc, Yt, Ybc]*
   (-58/15 + (36*Zeta[3])/5) + g1^2*MatMul[hb, htc, Yt, Ybc]*
   (-58/15 + (36*Zeta[3])/5) + g1^2*MatMul[Yb, htc, ht, Ybc]*
   (-58/15 + (36*Zeta[3])/5) + 2*g1^2*g2^2*mb*MatMul[Yb, Ybc]*
   (-43/5 + (42*Zeta[3])/5) + 2*g1^2*g3^2*mb*MatMul[Yb, Ybc]*
   (-24/5 + (32*Zeta[3])/3) + 2*mb*MatMul[Yb, Ybc, Yb, Ybc, Yb, Ybc]*
   (6 + 12*Zeta[3]) + g1^2*M1^2*MatMul[Yb, Ytc, Yt, Ybc]*
   (-116/15 + (72*Zeta[3])/5) + g1^2*g2^2*mq*Yb*Ybc*
   (-86/5 + (84*Zeta[3])/5) + g1^2*g3^2*mq*Yb*Ybc*(-48/5 + (64*Zeta[3])/3) + 
  mq*MatMul[Yb, Ybc, Yb]*MatMul[Ybc, Yb, Ybc]*(12 + 24*Zeta[3]) + 
  2*mb*MatMul[Yb, Ybc]*MatMul[Yb, Ybc, Yb, Ybc]*(12 + 24*Zeta[3]) + 
  mq*Ybc*MatMul[Yb, Ybc, Yb, Ybc, Yb]*(12 + 24*Zeta[3]) + 
  mq*Yb*MatMul[Ybc, Yb, Ybc, Yb, Ybc]*(12 + 24*Zeta[3]) + 
  MatMul[hb, hbc, Yb, Ybc, Yb, Ybc]*(12 + 24*Zeta[3]) + 
  MatMul[hb, Ybc, Yb, hbc, Yb, Ybc]*(12 + 24*Zeta[3]) + 
  MatMul[Yb, hbc, hb, Ybc, Yb, Ybc]*(12 + 24*Zeta[3]) + 
  MatMul[Yb, hbc, Yb, Ybc, hb, Ybc]*(12 + 24*Zeta[3]) + 
  MatMul[Yb, Ybc, hb, hbc, Yb, Ybc]*(12 + 24*Zeta[3]) + 
  MatMul[Yb, Ybc, Yb, hbc, hb, Ybc]*(12 + 24*Zeta[3]) + 
  g1^2*g2^2*M1^2*MatMul[Yb, Ybc]*(-172/5 + (168*Zeta[3])/5) + 
  g1^2*g2^2*M1*M2*MatMul[Yb, Ybc]*(-172/5 + (168*Zeta[3])/5) + 
  g1^2*g2^2*M2^2*MatMul[Yb, Ybc]*(-172/5 + (168*Zeta[3])/5) + 
  g2^2*M2*MatMul[hb, Ybc, Yb, Ybc]*(-18 + 36*Zeta[3]) + 
  g2^2*M2*MatMul[hb, Ytc, Yt, Ybc]*(-18 + 36*Zeta[3]) + 
  g2^2*M2*MatMul[Yb, hbc, Yb, Ybc]*(-18 + 36*Zeta[3]) + 
  g2^2*M2*MatMul[Yb, htc, Yt, Ybc]*(-18 + 36*Zeta[3]) + 
  g2^2*M2*MatMul[Yb, Ybc, hb, Ybc]*(-18 + 36*Zeta[3]) + 
  g2^2*M2*MatMul[Yb, Ytc, ht, Ybc]*(-18 + 36*Zeta[3]) + 
  g2^4*M2*MatMul[hb, Ybc]*(174 + 36*Zeta[3]) + g1^2*g3^2*M1^2*MatMul[Yb, Ybc]*
   (-96/5 + (128*Zeta[3])/3) + g1^2*g3^2*M1*M3*MatMul[Yb, Ybc]*
   (-96/5 + (128*Zeta[3])/3) + g1^2*g3^2*M3^2*MatMul[Yb, Ybc]*
   (-96/5 + (128*Zeta[3])/3) + 2*g2^2*g3^2*mb*MatMul[Yb, Ybc]*
   (-88 + 96*Zeta[3]) + g2^2*g3^2*mq*Yb*Ybc*(-176 + 192*Zeta[3]) + 
  g2^2*g3^2*M2^2*MatMul[Yb, Ybc]*(-352 + 384*Zeta[3]) + 
  g2^2*g3^2*M2*M3*MatMul[Yb, Ybc]*(-352 + 384*Zeta[3]) + 
  g2^2*g3^2*M3^2*MatMul[Yb, Ybc]*(-352 + 384*Zeta[3]) + 
  g3^4*M3*MatMul[hb, Ybc]*(-64/3 + (2176*Zeta[3])/3) + 
  g3^6*M3^2*(20512/9 + 7680*Zeta[3]) + 
  g3^4*hb*((32*hbc)/3 - (1088*hbc*Zeta[3])/3) + 
  g2^2*g3^2*M2*Yb*(176*hbc - 192*hbc*Zeta[3]) + 
  g2^2*g3^2*M3*Yb*(176*hbc - 192*hbc*Zeta[3]) + 
  g2^2*MatMul[hb, Ybc, Yb]*(18*hbc - 36*hbc*Zeta[3]) + 
  g2^2*MatMul[hb, Ytc, Yt]*(18*hbc - 36*hbc*Zeta[3]) + 
  g2^2*MatMul[Yb, Ybc, hb]*(18*hbc - 36*hbc*Zeta[3]) + 
  g2^2*MatMul[Yb, Ytc, ht]*(18*hbc - 36*hbc*Zeta[3]) + 
  g1^2*g3^2*M1*Yb*((48*hbc)/5 - (64*hbc*Zeta[3])/3) + 
  g1^2*g3^2*M3*Yb*((48*hbc)/5 - (64*hbc*Zeta[3])/3) + 
  g2^4*hb*(-87*hbc - 18*hbc*Zeta[3]) + g1^2*g2^2*M1*Yb*
   ((86*hbc)/5 - (84*hbc*Zeta[3])/5) + g1^2*g2^2*M2*Yb*
   ((86*hbc)/5 - (84*hbc*Zeta[3])/5) + g1^2*M1*MatMul[Yb, Ytc, Yt]*
   ((58*hbc)/15 - (36*hbc*Zeta[3])/5) + g1^2*M1*MatMul[Yb, Ybc, Yb]*
   ((2*hbc)/3 - (12*hbc*Zeta[3])/5) + 
  g1^4*hb*((-337*hbc)/15 - (14*hbc*Zeta[3])/75) + 
  g1^4*M1*Yb*((674*hbc)/15 + (28*hbc*Zeta[3])/75) + 
  g1^2*MatMul[hb, Ybc, Yb]*((-2*hbc)/3 + (12*hbc*Zeta[3])/5) + 
  g1^2*MatMul[Yb, Ybc, hb]*((-2*hbc)/3 + (12*hbc*Zeta[3])/5) + 
  g1^2*MatMul[hb, Ytc, Yt]*((-58*hbc)/15 + (36*hbc*Zeta[3])/5) + 
  g1^2*MatMul[Yb, Ytc, ht]*((-58*hbc)/15 + (36*hbc*Zeta[3])/5) + 
  g1^2*g2^2*hb*((-86*hbc)/5 + (84*hbc*Zeta[3])/5) + 
  g1^2*g3^2*hb*((-48*hbc)/5 + (64*hbc*Zeta[3])/3) + 
  MatMul[hb, Ybc, Yb, Ybc, Yb]*(12*hbc + 24*hbc*Zeta[3]) + 
  MatMul[Yb, Ybc, hb, Ybc, Yb]*(12*hbc + 24*hbc*Zeta[3]) + 
  MatMul[Yb, Ybc, Yb, Ybc, hb]*(12*hbc + 24*hbc*Zeta[3]) + 
  g2^2*M2*MatMul[Yb, Ybc, Yb]*(-18*hbc + 36*hbc*Zeta[3]) + 
  g2^2*M2*MatMul[Yb, Ytc, Yt]*(-18*hbc + 36*hbc*Zeta[3]) + 
  g2^4*M2*Yb*(174*hbc + 36*hbc*Zeta[3]) + 
  g2^2*g3^2*hb*(-176*hbc + 192*hbc*Zeta[3]) + 
  g3^4*M3*Yb*((-64*hbc)/3 + (2176*hbc*Zeta[3])/3) + 
  g3^4*MatMul[Yb, Ybc]*((32*mh1)/3 - (1088*mh1*Zeta[3])/3) + 
  g2^2*MatMul[Yb, Ybc, Yb, Ybc]*(36*mh1 - 72*mh1*Zeta[3]) + 
  g2^4*MatMul[Yb, Ybc]*(-99*mh1 - 12*mh2 - 12*trace[ml] - 36*trace[mq] - 
    18*mh1*Zeta[3]) + g1^4*MatMul[Yb, Ybc]*((-1721*mh1)/75 - (12*mh2)/25 - 
    (8*trace[mb])/25 - (24*trace[me])/25 - (12*trace[ml])/25 - 
    (4*trace[mq])/25 - (32*trace[mt])/25 - (14*mh1*Zeta[3])/75) + 
  g1^2*MatMul[Yb, Ybc, Yb, Ybc]*((-4*mh1)/3 + (24*mh1*Zeta[3])/5) + 
  g1^2*g2^2*MatMul[Yb, Ybc]*((-86*mh1)/5 + (84*mh1*Zeta[3])/5) + 
  g1^2*g3^2*MatMul[Yb, Ybc]*((-48*mh1)/5 + (64*mh1*Zeta[3])/3) + 
  MatMul[Yb, Ybc, Yb, Ybc, Yb, Ybc]*(36*mh1 + 72*mh1*Zeta[3]) + 
  g2^2*g3^2*MatMul[Yb, Ybc]*(-176*mh1 + 192*mh1*Zeta[3]) + 
  g2^2*MatMul[Yb, Ytc, Yt, Ybc]*(18*mh1 + 18*mh2 - 36*mh1*Zeta[3] - 
    36*mh2*Zeta[3]) + g1^2*MatMul[Yb, Ytc, Yt, Ybc]*
   ((-58*mh1)/15 - (58*mh2)/15 + (36*mh1*Zeta[3])/5 + (36*mh2*Zeta[3])/5) + 
  g3^2*Yb*(-32*hbc*trace[hb, Ybc] + 32*hbc*trace[he, Adj[Ye]] + 
    192*hbc*trace[hb, Ybc]*Zeta[3]) + g3^2*M3*MatMul[Yb, Ybc]*
   (32*trace[hb, Ybc] + 32*trace[hbc, Yb] - 32*trace[he, Adj[Ye]] - 
    32*trace[hec, Ye] - 192*trace[hb, Ybc]*Zeta[3] - 
    192*trace[hbc, Yb]*Zeta[3]) + g3^2*MatMul[hb, Ybc]*
   (-32*trace[hbc, Yb] + 32*trace[hec, Ye] + 192*trace[hbc, Yb]*Zeta[3]) + 
  g2^2*Yb*(54*hbc*trace[hb, Ybc] + 18*hbc*trace[he, Adj[Ye]] - 
    108*hbc*trace[hb, Ybc]*Zeta[3] - 36*hbc*trace[he, Adj[Ye]]*Zeta[3]) + 
  g1^2*Yb*(14*hbc*trace[hb, Ybc] - 6*hbc*trace[he, Adj[Ye]] - 
    12*hbc*trace[hb, Ybc]*Zeta[3] + 12*hbc*trace[he, Adj[Ye]]*Zeta[3]) + 
  g2^2*MatMul[hb, Ybc]*(54*trace[hbc, Yb] + 18*trace[hec, Ye] - 
    108*trace[hbc, Yb]*Zeta[3] - 36*trace[hec, Ye]*Zeta[3]) + 
  g1^2*M1*MatMul[Yb, Ybc]*(-14*trace[hb, Ybc] - 14*trace[hbc, Yb] + 
    6*trace[he, Adj[Ye]] + 6*trace[hec, Ye] + 12*trace[hb, Ybc]*Zeta[3] + 
    12*trace[hbc, Yb]*Zeta[3] - 12*trace[he, Adj[Ye]]*Zeta[3] - 
    12*trace[hec, Ye]*Zeta[3]) + g1^2*MatMul[hb, Ybc]*
   (14*trace[hbc, Yb] - 6*trace[hec, Ye] - 12*trace[hbc, Yb]*Zeta[3] + 
    12*trace[hec, Ye]*Zeta[3]) + g2^2*M2*MatMul[Yb, Ybc]*
   (-54*trace[hb, Ybc] - 54*trace[hbc, Yb] - 18*trace[he, Adj[Ye]] - 
    18*trace[hec, Ye] + 108*trace[hb, Ybc]*Zeta[3] + 
    108*trace[hbc, Yb]*Zeta[3] + 36*trace[he, Adj[Ye]]*Zeta[3] + 
    36*trace[hec, Ye]*Zeta[3]) + g3^2*M3*MatMul[hb, Ybc]*
   (32*trace[Ybc, Yb] - 32*trace[Adj[Ye], Ye] - 192*trace[Ybc, Yb]*Zeta[3]) + 
  2*g3^2*mb*MatMul[Yb, Ybc]*(-16*trace[Ybc, Yb] + 16*trace[Adj[Ye], Ye] + 
    96*trace[Ybc, Yb]*Zeta[3]) + g3^2*mq*Yb*Ybc*(-32*trace[Ybc, Yb] + 
    32*trace[Adj[Ye], Ye] + 192*trace[Ybc, Yb]*Zeta[3]) + 
  g3^2*M3^2*MatMul[Yb, Ybc]*(-64*trace[Ybc, Yb] + 64*trace[Adj[Ye], Ye] + 
    384*trace[Ybc, Yb]*Zeta[3]) + g3^2*M3*Yb*(32*hbc*trace[Ybc, Yb] - 
    32*hbc*trace[Adj[Ye], Ye] - 192*hbc*trace[Ybc, Yb]*Zeta[3]) + 
  g3^2*hb*(-32*hbc*trace[Ybc, Yb] + 32*hbc*trace[Adj[Ye], Ye] + 
    192*hbc*trace[Ybc, Yb]*Zeta[3]) + g2^2*M2^2*MatMul[Yb, Ybc]*
   (108*trace[Ybc, Yb] + 36*trace[Adj[Ye], Ye] - 216*trace[Ybc, Yb]*Zeta[3] - 
    72*trace[Adj[Ye], Ye]*Zeta[3]) + g2^2*mq*Yb*Ybc*
   (54*trace[Ybc, Yb] + 18*trace[Adj[Ye], Ye] - 108*trace[Ybc, Yb]*Zeta[3] - 
    36*trace[Adj[Ye], Ye]*Zeta[3]) + 2*g2^2*mb*MatMul[Yb, Ybc]*
   (27*trace[Ybc, Yb] + 9*trace[Adj[Ye], Ye] - 54*trace[Ybc, Yb]*Zeta[3] - 
    18*trace[Adj[Ye], Ye]*Zeta[3]) + g1^2*M1*MatMul[hb, Ybc]*
   (-14*trace[Ybc, Yb] + 6*trace[Adj[Ye], Ye] + 12*trace[Ybc, Yb]*Zeta[3] - 
    12*trace[Adj[Ye], Ye]*Zeta[3]) + 2*g1^2*mb*MatMul[Yb, Ybc]*
   (7*trace[Ybc, Yb] - 3*trace[Adj[Ye], Ye] - 6*trace[Ybc, Yb]*Zeta[3] + 
    6*trace[Adj[Ye], Ye]*Zeta[3]) + g1^2*mq*Yb*Ybc*
   (14*trace[Ybc, Yb] - 6*trace[Adj[Ye], Ye] - 12*trace[Ybc, Yb]*Zeta[3] + 
    12*trace[Adj[Ye], Ye]*Zeta[3]) + g1^2*M1^2*MatMul[Yb, Ybc]*
   (28*trace[Ybc, Yb] - 12*trace[Adj[Ye], Ye] - 24*trace[Ybc, Yb]*Zeta[3] + 
    24*trace[Adj[Ye], Ye]*Zeta[3]) + g2^2*M2*MatMul[hb, Ybc]*
   (-54*trace[Ybc, Yb] - 18*trace[Adj[Ye], Ye] + 108*trace[Ybc, Yb]*Zeta[3] + 
    36*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g2^2*hb*(54*hbc*trace[Ybc, Yb] + 18*hbc*trace[Adj[Ye], Ye] - 
    108*hbc*trace[Ybc, Yb]*Zeta[3] - 36*hbc*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g1^2*M1*Yb*(-14*hbc*trace[Ybc, Yb] + 6*hbc*trace[Adj[Ye], Ye] + 
    12*hbc*trace[Ybc, Yb]*Zeta[3] - 12*hbc*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g1^2*hb*(14*hbc*trace[Ybc, Yb] - 6*hbc*trace[Adj[Ye], Ye] - 
    12*hbc*trace[Ybc, Yb]*Zeta[3] + 12*hbc*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g2^2*M2*Yb*(-54*hbc*trace[Ybc, Yb] - 18*hbc*trace[Adj[Ye], Ye] + 
    108*hbc*trace[Ybc, Yb]*Zeta[3] + 36*hbc*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g3^2*MatMul[Yb, Ybc]*(-32*trace[hb, hbc] + 32*trace[he, hec] - 
    64*mh1*trace[Ybc, Yb] + 64*mh1*trace[Adj[Ye], Ye] - 
    32*trace[Yb, Ybc, mb] - 32*trace[Ybc, Yb, mq] + 
    32*trace[Ye, Adj[Ye], me] + 32*trace[Adj[Ye], Ye, ml] + 
    192*trace[hb, hbc]*Zeta[3] + 384*mh1*trace[Ybc, Yb]*Zeta[3] + 
    192*trace[Yb, Ybc, mb]*Zeta[3] + 192*trace[Ybc, Yb, mq]*Zeta[3]) + 
  g2^2*MatMul[Yb, Ybc]*(54*trace[hb, hbc] + 18*trace[he, hec] + 
    108*mh1*trace[Ybc, Yb] + 36*mh1*trace[Adj[Ye], Ye] + 
    54*trace[Yb, Ybc, mb] + 54*trace[Ybc, Yb, mq] + 
    18*trace[Ye, Adj[Ye], me] + 18*trace[Adj[Ye], Ye, ml] - 
    108*trace[hb, hbc]*Zeta[3] - 36*trace[he, hec]*Zeta[3] - 
    216*mh1*trace[Ybc, Yb]*Zeta[3] - 72*mh1*trace[Adj[Ye], Ye]*Zeta[3] - 
    108*trace[Yb, Ybc, mb]*Zeta[3] - 108*trace[Ybc, Yb, mq]*Zeta[3] - 
    36*trace[Ye, Adj[Ye], me]*Zeta[3] - 36*trace[Adj[Ye], Ye, ml]*Zeta[3]) + 
  g1^2*MatMul[Yb, Ybc]*(14*trace[hb, hbc] - 6*trace[he, hec] + 
    28*mh1*trace[Ybc, Yb] - 12*mh1*trace[Adj[Ye], Ye] + 
    14*trace[Yb, Ybc, mb] + 14*trace[Ybc, Yb, mq] - 
    6*trace[Ye, Adj[Ye], me] - 6*trace[Adj[Ye], Ye, ml] - 
    12*trace[hb, hbc]*Zeta[3] + 12*trace[he, hec]*Zeta[3] - 
    24*mh1*trace[Ybc, Yb]*Zeta[3] + 24*mh1*trace[Adj[Ye], Ye]*Zeta[3] - 
    12*trace[Yb, Ybc, mb]*Zeta[3] - 12*trace[Ybc, Yb, mq]*Zeta[3] + 
    12*trace[Ye, Adj[Ye], me]*Zeta[3] + 12*trace[Adj[Ye], Ye, ml]*Zeta[3])}

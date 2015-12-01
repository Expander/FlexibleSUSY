{4*ht*htc - (32*g1^2*M1^2)/15 - (32*g3^2*M3^2)/3 + 4*mq*Yt*Adj[Yt] + 
  4*mh2*MatMul[Yt, Adj[Yt]] + 4*mt*MatMul[Yt, Adj[Yt]], 
 (-4*g1^2*ht*htc)/5 + 12*g2^2*ht*htc + (3424*g1^4*M1^2)/75 + 
  (512*g1^2*g3^2*M1^2)/45 + (512*g1^2*g3^2*M1*M3)/45 + 
  (512*g1^2*g3^2*M3^2)/45 - (128*g3^4*M3^2)/3 + (4*g1^2*htc*M1*Yt)/5 - 
  12*g2^2*htc*M2*Yt - (4*g1^2*mq*Yt*Adj[Yt])/5 + 12*g2^2*mq*Yt*Adj[Yt] + 
  (4*g1^2*M1*MatMul[ht, Adj[Yt]])/5 - 12*g2^2*M2*MatMul[ht, Adj[Yt]] - 
  4*mb*MatMul[Yb, Adj[Yt]]*MatMul[Yt, Adj[Yb]] - 
  (8*g1^2*M1^2*MatMul[Yt, Adj[Yt]])/5 + 24*g2^2*M2^2*MatMul[Yt, Adj[Yt]] - 
  (4*g1^2*mh2*MatMul[Yt, Adj[Yt]])/5 + 12*g2^2*mh2*MatMul[Yt, Adj[Yt]] - 
  (4*g1^2*mt*MatMul[Yt, Adj[Yt]])/5 + 12*g2^2*mt*MatMul[Yt, Adj[Yt]] - 
  4*mt*MatMul[Yt, Adj[Yt]]^2 - 4*htc*MatMul[ht, Adj[Yb], Yb] - 
  4*htc*MatMul[ht, Adj[Yt], Yt] - 4*htc*MatMul[Yt, Adj[Yb], hb] - 
  4*mq*Adj[Yt]*MatMul[Yt, Adj[Yb], Yb] - 4*htc*MatMul[Yt, Adj[Yt], ht] - 
  4*mq*Adj[Yt]*MatMul[Yt, Adj[Yt], Yt] - 
  4*mq*Yt*MatMul[Adj[Yb], Yb, Adj[Yt]] - 
  4*mq*Yt*MatMul[Adj[Yt], Yt, Adj[Yt]] - 4*MatMul[ht, hbc, Yb, Adj[Yt]] - 
  4*MatMul[ht, htc, Yt, Adj[Yt]] - 4*MatMul[Yt, hbc, hb, Adj[Yt]] - 
  4*MatMul[Yt, htc, ht, Adj[Yt]] + (-4*mh1 - 4*mh2)*
   MatMul[Yt, Adj[Yb], Yb, Adj[Yt]] - 4*mt*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]] - 
  8*mh2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]] - 
  4*mt*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]] + 
  g1^4*((16*mh1)/25 + (16*mh2)/25 + (32*trace[mb])/75 + (32*trace[me])/25 + 
    (16*trace[ml])/25 + (16*trace[mq])/75 + (128*trace[mt])/75) + 
  g3^4*((16*trace[mb])/3 + (32*trace[mq])/3 + (16*trace[mt])/3) - 
  12*htc*Yt*trace[ht, Adj[Yt]] - 12*MatMul[ht, Adj[Yt]]*trace[htc, Yt] - 
  12*ht*htc*trace[Adj[Yt], Yt] - 12*mq*Yt*Adj[Yt]*trace[Adj[Yt], Yt] - 
  12*mt*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  MatMul[Yt, Adj[Yt]]*(-12*trace[ht, htc] - 24*mh2*trace[Adj[Yt], Yt] - 
    12*trace[Yt, Adj[Yt], mt] - 12*trace[Adj[Yt], Yt, mq]), 
 (128*g3^2*mb*MatMul[Yb, Adj[Yt]]*MatMul[Yt, Adj[Yb]])/3 + 
  (128*g3^2*mt*MatMul[Yt, Adj[Yt]]^2)/3 + 
  (128*g3^2*htc*MatMul[ht, Adj[Yb], Yb])/3 + 
  (128*g3^2*htc*MatMul[ht, Adj[Yt], Yt])/3 + 
  (128*g3^2*htc*MatMul[Yt, Adj[Yb], hb])/3 - 
  (128*g3^2*htc*M3*MatMul[Yt, Adj[Yb], Yb])/3 + 
  (128*g3^2*mq*Adj[Yt]*MatMul[Yt, Adj[Yb], Yb])/3 + 
  (128*g3^2*htc*MatMul[Yt, Adj[Yt], ht])/3 - 
  (128*g3^2*htc*M3*MatMul[Yt, Adj[Yt], Yt])/3 + 
  (128*g3^2*mq*Adj[Yt]*MatMul[Yt, Adj[Yt], Yt])/3 + 
  (128*g3^2*mq*Yt*MatMul[Adj[Yb], Yb, Adj[Yt]])/3 + 
  12*mq*MatMul[Yt, Adj[Yb], Yb]*MatMul[Adj[Yb], Yb, Adj[Yt]] - 
  4*mq*MatMul[Yt, Adj[Yt], Yt]*MatMul[Adj[Yb], Yb, Adj[Yt]] + 
  (128*g3^2*mq*Yt*MatMul[Adj[Yt], Yt, Adj[Yt]])/3 - 
  4*mq*MatMul[Yt, Adj[Yb], Yb]*MatMul[Adj[Yt], Yt, Adj[Yt]] + 
  (128*g3^2*MatMul[ht, hbc, Yb, Adj[Yt]])/3 + 
  (128*g3^2*MatMul[ht, htc, Yt, Adj[Yt]])/3 - 
  (128*g3^2*M3*MatMul[ht, Adj[Yb], Yb, Adj[Yt]])/3 - 
  (128*g3^2*M3*MatMul[ht, Adj[Yt], Yt, Adj[Yt]])/3 + 
  12*mb*MatMul[Yt, Adj[Yb]]*MatMul[Yb, Adj[Yb], Yb, Adj[Yt]] - 
  4*mb*MatMul[Yt, Adj[Yb]]*MatMul[Yb, Adj[Yt], Yt, Adj[Yt]] + 
  (128*g3^2*MatMul[Yt, hbc, hb, Adj[Yt]])/3 - 
  (128*g3^2*M3*MatMul[Yt, hbc, Yb, Adj[Yt]])/3 + 
  (128*g3^2*MatMul[Yt, htc, ht, Adj[Yt]])/3 - 
  (128*g3^2*M3*MatMul[Yt, htc, Yt, Adj[Yt]])/3 - 
  (128*g3^2*M3*MatMul[Yt, Adj[Yb], hb, Adj[Yt]])/3 + 
  12*mb*MatMul[Yb, Adj[Yt]]*MatMul[Yt, Adj[Yb], Yb, Adj[Yb]] + 
  (256*g3^2*M3^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]])/3 + 
  g3^2*((128*mh1)/3 + (128*mh2)/3)*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]] + 
  (128*g3^2*mt*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]])/3 - 
  8*mt*MatMul[Yt, Adj[Yt]]*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]] - 
  (128*g3^2*M3*MatMul[Yt, Adj[Yt], ht, Adj[Yt]])/3 - 
  4*mb*MatMul[Yb, Adj[Yt]]*MatMul[Yt, Adj[Yt], Yt, Adj[Yb]] + 
  (256*g3^2*M3^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]])/3 + 
  (256*g3^2*mh2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]])/3 + 
  (128*g3^2*mt*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]])/3 + 
  12*htc*MatMul[ht, Adj[Yb], Yb, Adj[Yb], Yb] - 
  4*htc*MatMul[ht, Adj[Yb], Yb, Adj[Yt], Yt] - 
  4*htc*MatMul[ht, Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*htc*MatMul[Yt, Adj[Yb], hb, Adj[Yb], Yb] - 
  4*htc*MatMul[Yt, Adj[Yb], hb, Adj[Yt], Yt] + 
  12*htc*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], hb] + 
  12*mq*Adj[Yt]*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb] - 
  4*htc*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], ht] - 
  4*mq*Adj[Yt]*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt] - 
  4*htc*MatMul[Yt, Adj[Yt], ht, Adj[Yb], Yb] - 
  4*htc*MatMul[Yt, Adj[Yt], Yt, Adj[Yb], hb] - 
  4*mq*Adj[Yt]*MatMul[Yt, Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*mq*Yt*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yt]] - 
  4*mq*Yt*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yt]] - 
  4*mq*Yt*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt]] + 
  12*MatMul[ht, hbc, Yb, Adj[Yb], Yb, Adj[Yt]] - 
  4*MatMul[ht, hbc, Yb, Adj[Yt], Yt, Adj[Yt]] - 
  4*MatMul[ht, htc, Yt, Adj[Yb], Yb, Adj[Yt]] + 
  12*MatMul[ht, Adj[Yb], Yb, hbc, Yb, Adj[Yt]] - 
  4*MatMul[ht, Adj[Yb], Yb, htc, Yt, Adj[Yt]] - 
  4*MatMul[ht, Adj[Yt], Yt, hbc, Yb, Adj[Yt]] + 
  12*MatMul[Yt, hbc, hb, Adj[Yb], Yb, Adj[Yt]] - 
  4*MatMul[Yt, hbc, hb, Adj[Yt], Yt, Adj[Yt]] + 
  12*MatMul[Yt, hbc, Yb, Adj[Yb], hb, Adj[Yt]] - 
  4*MatMul[Yt, hbc, Yb, Adj[Yt], ht, Adj[Yt]] - 
  4*MatMul[Yt, htc, ht, Adj[Yb], Yb, Adj[Yt]] - 
  4*MatMul[Yt, htc, Yt, Adj[Yb], hb, Adj[Yt]] + 
  12*MatMul[Yt, Adj[Yb], hb, hbc, Yb, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yb], hb, htc, Yt, Adj[Yt]] + 
  12*MatMul[Yt, Adj[Yb], Yb, hbc, hb, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yb], Yb, htc, ht, Adj[Yt]] + 
  (24*mh1 + 12*mh2)*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yt]] + 
  12*mt*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yt]] + 
  (-4*mh1 - 8*mh2)*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yt]] - 
  4*mt*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yt], ht, hbc, Yb, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yt], Yt, hbc, hb, Adj[Yt]] + 
  (-4*mh1 - 8*mh2)*MatMul[Yt, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt]] - 
  4*mt*MatMul[Yt, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt]] + 
  g1^6*((-3424*mh1)/375 - (3424*mh2)/375 - (6848*trace[mb])/1125 - 
    (6848*trace[me])/375 - (3424*trace[ml])/375 - (3424*trace[mq])/1125 - 
    (27392*trace[mt])/1125) + g1^4*g3^2*((-256*mh1)/75 - (256*mh2)/75 - 
    (512*trace[mb])/225 - (512*trace[me])/75 - (256*trace[ml])/75 - 
    (256*trace[mq])/225 - (2048*trace[mt])/225) + 
  g1^2*g3^4*((-256*trace[mb])/45 - (512*trace[mq])/45 - (256*trace[mt])/45) + 
  g3^6*((320*trace[mb])/9 + (640*trace[mq])/9 + (320*trace[mt])/9) + 
  MatMul[Yt, hbc, Yb, Adj[Yt]]*(24*trace[hb, Adj[Yb]] + 
    8*trace[he, Adj[Ye]] - 12*trace[ht, Adj[Yt]]) + 
  12*htc*MatMul[Yt, Adj[Yt], Yt]*trace[ht, Adj[Yt]] + 
  12*MatMul[Yt, htc, Yt, Adj[Yt]]*trace[ht, Adj[Yt]] + 
  MatMul[Yt, Adj[Yb], Yb]*(24*htc*trace[hb, Adj[Yb]] + 
    8*htc*trace[he, Adj[Ye]] - 12*htc*trace[ht, Adj[Yt]]) + 
  MatMul[ht, Adj[Yb], Yb, Adj[Yt]]*(24*trace[hbc, Yb] + 8*trace[hec, Ye] - 
    12*trace[htc, Yt]) + MatMul[Yt, Adj[Yb], hb, Adj[Yt]]*
   (24*trace[hbc, Yb] + 8*trace[hec, Ye] - 12*trace[htc, Yt]) + 
  12*MatMul[ht, Adj[Yt], Yt, Adj[Yt]]*trace[htc, Yt] + 
  12*MatMul[Yt, Adj[Yt], ht, Adj[Yt]]*trace[htc, Yt] + 
  g1^4*M1*((224*trace[hb, Adj[Yb]])/15 + (224*trace[hbc, Yb])/15 + 
    (96*trace[he, Adj[Ye]])/5 + (96*trace[hec, Ye])/5 + 
    (416*trace[ht, Adj[Yt]])/15 + (416*trace[htc, Yt])/15) + 
  g3^4*M3*((320*trace[hb, Adj[Yb]])/3 + (320*trace[hbc, Yb])/3 + 
    (320*trace[ht, Adj[Yt]])/3 + (320*trace[htc, Yt])/3) + 
  g3^4*M3^2*(-320*trace[Adj[Yb], Yb] - 320*trace[Adj[Yt], Yt]) + 
  g1^4*M1^2*((-224*trace[Adj[Yb], Yb])/5 - (288*trace[Adj[Ye], Ye])/5 - 
    (416*trace[Adj[Yt], Yt])/5) + mb*MatMul[Yb, Adj[Yt]]*MatMul[Yt, Adj[Yb]]*
   (24*trace[Adj[Yb], Yb] + 8*trace[Adj[Ye], Ye] - 12*trace[Adj[Yt], Yt]) + 
  mq*Adj[Yt]*MatMul[Yt, Adj[Yb], Yb]*(24*trace[Adj[Yb], Yb] + 
    8*trace[Adj[Ye], Ye] - 12*trace[Adj[Yt], Yt]) + 
  mq*Yt*MatMul[Adj[Yb], Yb, Adj[Yt]]*(24*trace[Adj[Yb], Yb] + 
    8*trace[Adj[Ye], Ye] - 12*trace[Adj[Yt], Yt]) + 
  MatMul[ht, hbc, Yb, Adj[Yt]]*(24*trace[Adj[Yb], Yb] + 
    8*trace[Adj[Ye], Ye] - 12*trace[Adj[Yt], Yt]) + 
  MatMul[Yt, hbc, hb, Adj[Yt]]*(24*trace[Adj[Yb], Yb] + 
    8*trace[Adj[Ye], Ye] - 12*trace[Adj[Yt], Yt]) + 
  2*mt*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*(12*trace[Adj[Yb], Yb] + 
    4*trace[Adj[Ye], Ye] - 6*trace[Adj[Yt], Yt]) + 
  12*mt*MatMul[Yt, Adj[Yt]]^2*trace[Adj[Yt], Yt] + 
  12*htc*MatMul[ht, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*htc*MatMul[Yt, Adj[Yt], ht]*trace[Adj[Yt], Yt] + 
  12*mq*Adj[Yt]*MatMul[Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*mq*Yt*MatMul[Adj[Yt], Yt, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  12*MatMul[ht, htc, Yt, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  12*MatMul[Yt, htc, ht, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  12*mt*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  MatMul[ht, Adj[Yb], Yb]*(24*htc*trace[Adj[Yb], Yb] + 
    8*htc*trace[Adj[Ye], Ye] - 12*htc*trace[Adj[Yt], Yt]) + 
  MatMul[Yt, Adj[Yb], hb]*(24*htc*trace[Adj[Yb], Yb] + 
    8*htc*trace[Adj[Ye], Ye] - 12*htc*trace[Adj[Yt], Yt]) + 
  g3^4*(-64*trace[hb, hbc] - 64*trace[ht, htc] - 64*mh1*trace[Adj[Yb], Yb] - 
    64*mh2*trace[Adj[Yt], Yt] - 64*trace[Yb, Adj[Yb], mb] - 
    64*trace[Yt, Adj[Yt], mt] - 64*trace[Adj[Yb], Yb, mq] - 
    64*trace[Adj[Yt], Yt, mq]) + g1^4*((-224*trace[hb, hbc])/25 - 
    (288*trace[he, hec])/25 - (416*trace[ht, htc])/25 - 
    (224*mh1*trace[Adj[Yb], Yb])/25 - (288*mh1*trace[Adj[Ye], Ye])/25 - 
    (416*mh2*trace[Adj[Yt], Yt])/25 - (224*trace[Yb, Adj[Yb], mb])/25 - 
    (288*trace[Ye, Adj[Ye], me])/25 - (416*trace[Yt, Adj[Yt], mt])/25 - 
    (224*trace[Adj[Yb], Yb, mq])/25 - (288*trace[Adj[Ye], Ye, ml])/25 - 
    (416*trace[Adj[Yt], Yt, mq])/25) + MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*
   (24*trace[hb, hbc] + 8*trace[he, hec] - 12*trace[ht, htc] + 
    48*mh1*trace[Adj[Yb], Yb] + 24*mh2*trace[Adj[Yb], Yb] + 
    16*mh1*trace[Adj[Ye], Ye] + 8*mh2*trace[Adj[Ye], Ye] - 
    12*mh1*trace[Adj[Yt], Yt] - 24*mh2*trace[Adj[Yt], Yt] + 
    24*trace[Yb, Adj[Yb], mb] + 8*trace[Ye, Adj[Ye], me] - 
    12*trace[Yt, Adj[Yt], mt] + 24*trace[Adj[Yb], Yb, mq] + 
    8*trace[Adj[Ye], Ye, ml] - 12*trace[Adj[Yt], Yt, mq]) + 
  MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]*(12*trace[ht, htc] + 
    36*mh2*trace[Adj[Yt], Yt] + 12*trace[Yt, Adj[Yt], mt] + 
    12*trace[Adj[Yt], Yt, mq]) + 
  Yt*(-72*htc*trace[ht, Adj[Yt]]*trace[Adj[Yt], Yt] + 
    24*htc*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
    24*htc*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
    144*htc*trace[Adj[Yt], ht, Adj[Yt], Yt]) + 
  MatMul[ht, Adj[Yt]]*(-72*trace[htc, Yt]*trace[Adj[Yt], Yt] + 
    24*trace[htc, Yt, Adj[Yb], Yb] + 24*trace[Adj[Yt], Yt, hbc, Yb] + 
    144*trace[Adj[Yt], Yt, htc, Yt]) + 2*mt*MatMul[Yt, Adj[Yt]]*
   (-18*trace[Adj[Yt], Yt]^2 + 12*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    36*trace[Adj[Yt], Yt, Adj[Yt], Yt]) + 
  mq*Yt*Adj[Yt]*(-36*trace[Adj[Yt], Yt]^2 + 
    24*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    72*trace[Adj[Yt], Yt, Adj[Yt], Yt]) + 
  ht*(-36*htc*trace[Adj[Yt], Yt]^2 + 24*htc*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    72*htc*trace[Adj[Yt], Yt, Adj[Yt], Yt]) + 
  MatMul[Yt, Adj[Yt]]*(-72*trace[ht, Adj[Yt]]*trace[htc, Yt] - 
    72*trace[ht, htc]*trace[Adj[Yt], Yt] - 108*mh2*trace[Adj[Yt], Yt]^2 - 
    72*trace[Adj[Yt], Yt]*trace[Yt, Adj[Yt], mt] - 
    72*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, mq] + 
    24*trace[hb, htc, Yt, Adj[Yb]] + 24*trace[hbc, hb, Adj[Yt], Yt] + 
    24*trace[htc, ht, Adj[Yb], Yb] + 144*trace[htc, ht, Adj[Yt], Yt] + 
    24*trace[Adj[Yt], ht, hbc, Yb] + 144*trace[Adj[Yt], ht, htc, Yt] + 
    24*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    48*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    216*mh2*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
    24*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
    24*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
    144*trace[Yt, Adj[Yt], Yt, Adj[Yt], mt] + 
    24*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
    24*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
    144*trace[Adj[Yt], Yt, Adj[Yt], Yt, mq]) + g3^4*M3^2*MatMul[Yt, Adj[Yt]]*
   (64 - 2176*Zeta[3]) + g2^2*g3^4*M3^2*(720 - 864*Zeta[3]) + 
  g2^2*g3^4*M2*M3*(480 - 576*Zeta[3]) + g3^4*mq*Yt*Adj[Yt]*
   (32/3 - (1088*Zeta[3])/3) + g1^4*g3^2*M1^2*(8576/75 - (8448*Zeta[3])/25) + 
  g1^6*M1^2*(864496/1125 - (38208*Zeta[3])/125) + 
  g2^2*g3^4*M2^2*(288 - 288*Zeta[3]) + g1^4*g3^2*M1*M3*
   (17152/225 - (5632*Zeta[3])/25) + g1^2*g3^4*M3^2*
   (-176/15 - (1056*Zeta[3])/5) + g2^2*g3^2*M2*MatMul[ht, Adj[Yt]]*
   (176 - 192*Zeta[3]) + g2^2*g3^2*M3*MatMul[ht, Adj[Yt]]*
   (176 - 192*Zeta[3]) + 2*g3^4*mt*MatMul[Yt, Adj[Yt]]*
   (16/3 - (544*Zeta[3])/3) + g1^2*g3^4*M1*M3*(-1376/45 - (704*Zeta[3])/5) + 
  g1^4*g3^2*M3^2*(512/9 - (2816*Zeta[3])/25) + g2^4*M2^2*MatMul[Yt, Adj[Yt]]*
   (-474 - 108*Zeta[3]) + g1^4*g2^2*M1^2*(432/5 - (2592*Zeta[3])/25) + 
  g2^2*M2^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*(36 - 72*Zeta[3]) + 
  g2^2*M2^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]*(36 - 72*Zeta[3]) + 
  g1^2*g3^4*M1^2*(-32/9 - (352*Zeta[3])/5) + 
  g1^4*g2^2*M1*M2*(288/5 - (1728*Zeta[3])/25) + 
  g1^2*g3^2*M1^2*MatMul[Yt, Adj[Yt]]*(-32/15 - (896*Zeta[3])/15) + 
  g1^2*g3^2*M1*M3*MatMul[Yt, Adj[Yt]]*(-32/15 - (896*Zeta[3])/15) + 
  g1^2*g3^2*M3^2*MatMul[Yt, Adj[Yt]]*(-32/15 - (896*Zeta[3])/15) + 
  g1^4*M1^2*MatMul[Yt, Adj[Yt]]*(-4794/25 - (988*Zeta[3])/25) + 
  g2^2*mb*MatMul[Yb, Adj[Yt]]*MatMul[Yt, Adj[Yb]]*(18 - 36*Zeta[3]) + 
  g2^2*mt*MatMul[Yt, Adj[Yt]]^2*(18 - 36*Zeta[3]) + 
  g2^2*mq*Adj[Yt]*MatMul[Yt, Adj[Yb], Yb]*(18 - 36*Zeta[3]) + 
  g2^2*mq*Adj[Yt]*MatMul[Yt, Adj[Yt], Yt]*(18 - 36*Zeta[3]) + 
  g2^2*mq*Yt*MatMul[Adj[Yb], Yb, Adj[Yt]]*(18 - 36*Zeta[3]) + 
  g2^2*mq*Yt*MatMul[Adj[Yt], Yt, Adj[Yt]]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[ht, hbc, Yb, Adj[Yt]]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[ht, htc, Yt, Adj[Yt]]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[Yt, hbc, hb, Adj[Yt]]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[Yt, htc, ht, Adj[Yt]]*(18 - 36*Zeta[3]) + 
  g1^4*g2^2*M2^2*(864/25 - (864*Zeta[3])/25) + 
  g1^2*g2^2*M1*MatMul[ht, Adj[Yt]]*(134/5 - (156*Zeta[3])/5) + 
  g1^2*g2^2*M2*MatMul[ht, Adj[Yt]]*(134/5 - (156*Zeta[3])/5) + 
  g1^2*g3^2*mq*Yt*Adj[Yt]*(-16/15 - (448*Zeta[3])/15) + 
  g2^4*mq*Yt*Adj[Yt]*(-87 - 18*Zeta[3]) + 
  2*g2^2*mt*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*(9 - 18*Zeta[3]) + 
  2*g2^2*mt*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]*(9 - 18*Zeta[3]) + 
  2*g1^2*g3^2*mt*MatMul[Yt, Adj[Yt]]*(-8/15 - (224*Zeta[3])/15) + 
  g1^2*M1*MatMul[ht, Adj[Yt], Yt, Adj[Yt]]*(2/3 - 12*Zeta[3]) + 
  g1^2*M1*MatMul[Yt, htc, Yt, Adj[Yt]]*(2/3 - 12*Zeta[3]) + 
  g1^2*M1*MatMul[Yt, Adj[Yt], ht, Adj[Yt]]*(2/3 - 12*Zeta[3]) + 
  2*g2^4*mt*MatMul[Yt, Adj[Yt]]*(-87/2 - 9*Zeta[3]) + 
  g1^2*M1*MatMul[ht, Adj[Yb], Yb, Adj[Yt]]*(-38/15 - (36*Zeta[3])/5) + 
  g1^2*M1*MatMul[Yt, hbc, Yb, Adj[Yt]]*(-38/15 - (36*Zeta[3])/5) + 
  g1^2*M1*MatMul[Yt, Adj[Yb], hb, Adj[Yt]]*(-38/15 - (36*Zeta[3])/5) + 
  g1^4*mq*Yt*Adj[Yt]*(-799/25 - (494*Zeta[3])/75) + 
  2*g1^4*mt*MatMul[Yt, Adj[Yt]]*(-799/50 - (247*Zeta[3])/75) + 
  2*g1^2*mt*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*(19/15 + (18*Zeta[3])/5) + 
  2*g1^2*mt*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]*(-1/3 + 6*Zeta[3]) + 
  g1^2*mb*MatMul[Yb, Adj[Yt]]*MatMul[Yt, Adj[Yb]]*(38/15 + (36*Zeta[3])/5) + 
  g1^2*mq*Adj[Yt]*MatMul[Yt, Adj[Yb], Yb]*(38/15 + (36*Zeta[3])/5) + 
  g1^2*mq*Yt*MatMul[Adj[Yb], Yb, Adj[Yt]]*(38/15 + (36*Zeta[3])/5) + 
  g1^2*MatMul[ht, hbc, Yb, Adj[Yt]]*(38/15 + (36*Zeta[3])/5) + 
  g1^2*MatMul[Yt, hbc, hb, Adj[Yt]]*(38/15 + (36*Zeta[3])/5) + 
  g1^2*mt*MatMul[Yt, Adj[Yt]]^2*(-2/3 + 12*Zeta[3]) + 
  g1^2*mq*Adj[Yt]*MatMul[Yt, Adj[Yt], Yt]*(-2/3 + 12*Zeta[3]) + 
  g1^2*mq*Yt*MatMul[Adj[Yt], Yt, Adj[Yt]]*(-2/3 + 12*Zeta[3]) + 
  g1^2*MatMul[ht, htc, Yt, Adj[Yt]]*(-2/3 + 12*Zeta[3]) + 
  g1^2*MatMul[Yt, htc, ht, Adj[Yt]]*(-2/3 + 12*Zeta[3]) + 
  2*mt*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt]]*(6 + 12*Zeta[3]) + 
  g1^4*M1*MatMul[ht, Adj[Yt]]*(1598/25 + (988*Zeta[3])/75) + 
  g1^2*M1^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*(76/15 + (72*Zeta[3])/5) + 
  2*g1^2*g2^2*mt*MatMul[Yt, Adj[Yt]]*(-67/5 + (78*Zeta[3])/5) + 
  g1^2*M1^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]*(-4/3 + 24*Zeta[3]) + 
  mq*MatMul[Yt, Adj[Yt], Yt]*MatMul[Adj[Yt], Yt, Adj[Yt]]*(12 + 24*Zeta[3]) + 
  2*mt*MatMul[Yt, Adj[Yt]]*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]*
   (12 + 24*Zeta[3]) + mq*Adj[Yt]*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt]*
   (12 + 24*Zeta[3]) + mq*Yt*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt]]*
   (12 + 24*Zeta[3]) + MatMul[ht, htc, Yt, Adj[Yt], Yt, Adj[Yt]]*
   (12 + 24*Zeta[3]) + MatMul[ht, Adj[Yt], Yt, htc, Yt, Adj[Yt]]*
   (12 + 24*Zeta[3]) + MatMul[Yt, htc, ht, Adj[Yt], Yt, Adj[Yt]]*
   (12 + 24*Zeta[3]) + MatMul[Yt, htc, Yt, Adj[Yt], ht, Adj[Yt]]*
   (12 + 24*Zeta[3]) + MatMul[Yt, Adj[Yt], ht, htc, Yt, Adj[Yt]]*
   (12 + 24*Zeta[3]) + MatMul[Yt, Adj[Yt], Yt, htc, ht, Adj[Yt]]*
   (12 + 24*Zeta[3]) + g1^2*g3^2*M1*MatMul[ht, Adj[Yt]]*
   (16/15 + (448*Zeta[3])/15) + g1^2*g3^2*M3*MatMul[ht, Adj[Yt]]*
   (16/15 + (448*Zeta[3])/15) + g1^2*g2^2*mq*Yt*Adj[Yt]*
   (-134/5 + (156*Zeta[3])/5) + g2^2*M2*MatMul[ht, Adj[Yb], Yb, Adj[Yt]]*
   (-18 + 36*Zeta[3]) + g2^2*M2*MatMul[ht, Adj[Yt], Yt, Adj[Yt]]*
   (-18 + 36*Zeta[3]) + g2^2*M2*MatMul[Yt, hbc, Yb, Adj[Yt]]*
   (-18 + 36*Zeta[3]) + g2^2*M2*MatMul[Yt, htc, Yt, Adj[Yt]]*
   (-18 + 36*Zeta[3]) + g2^2*M2*MatMul[Yt, Adj[Yb], hb, Adj[Yt]]*
   (-18 + 36*Zeta[3]) + g2^2*M2*MatMul[Yt, Adj[Yt], ht, Adj[Yt]]*
   (-18 + 36*Zeta[3]) + g2^4*M2*MatMul[ht, Adj[Yt]]*(174 + 36*Zeta[3]) + 
  g1^2*g2^2*M1^2*MatMul[Yt, Adj[Yt]]*(-268/5 + (312*Zeta[3])/5) + 
  g1^2*g2^2*M1*M2*MatMul[Yt, Adj[Yt]]*(-268/5 + (312*Zeta[3])/5) + 
  g1^2*g2^2*M2^2*MatMul[Yt, Adj[Yt]]*(-268/5 + (312*Zeta[3])/5) + 
  2*g2^2*g3^2*mt*MatMul[Yt, Adj[Yt]]*(-88 + 96*Zeta[3]) + 
  g2^2*g3^2*mq*Yt*Adj[Yt]*(-176 + 192*Zeta[3]) + 
  g2^2*g3^2*M2^2*MatMul[Yt, Adj[Yt]]*(-352 + 384*Zeta[3]) + 
  g2^2*g3^2*M2*M3*MatMul[Yt, Adj[Yt]]*(-352 + 384*Zeta[3]) + 
  g2^2*g3^2*M3^2*MatMul[Yt, Adj[Yt]]*(-352 + 384*Zeta[3]) + 
  g3^4*M3*MatMul[ht, Adj[Yt]]*(-64/3 + (2176*Zeta[3])/3) + 
  g3^6*M3^2*(20512/9 + 7680*Zeta[3]) + 
  g3^4*ht*((32*htc)/3 - (1088*htc*Zeta[3])/3) + 
  g2^2*g3^2*M2*Yt*(176*htc - 192*htc*Zeta[3]) + 
  g2^2*g3^2*M3*Yt*(176*htc - 192*htc*Zeta[3]) + 
  g2^2*MatMul[ht, Adj[Yb], Yb]*(18*htc - 36*htc*Zeta[3]) + 
  g2^2*MatMul[ht, Adj[Yt], Yt]*(18*htc - 36*htc*Zeta[3]) + 
  g2^2*MatMul[Yt, Adj[Yb], hb]*(18*htc - 36*htc*Zeta[3]) + 
  g2^2*MatMul[Yt, Adj[Yt], ht]*(18*htc - 36*htc*Zeta[3]) + 
  g1^2*g2^2*M1*Yt*((134*htc)/5 - (156*htc*Zeta[3])/5) + 
  g1^2*g2^2*M2*Yt*((134*htc)/5 - (156*htc*Zeta[3])/5) + 
  g1^2*g3^2*ht*((-16*htc)/15 - (448*htc*Zeta[3])/15) + 
  g2^4*ht*(-87*htc - 18*htc*Zeta[3]) + g1^2*M1*MatMul[Yt, Adj[Yt], Yt]*
   ((2*htc)/3 - 12*htc*Zeta[3]) + g1^2*M1*MatMul[Yt, Adj[Yb], Yb]*
   ((-38*htc)/15 - (36*htc*Zeta[3])/5) + 
  g1^4*ht*((-799*htc)/25 - (494*htc*Zeta[3])/75) + 
  g1^2*MatMul[ht, Adj[Yb], Yb]*((38*htc)/15 + (36*htc*Zeta[3])/5) + 
  g1^2*MatMul[Yt, Adj[Yb], hb]*((38*htc)/15 + (36*htc*Zeta[3])/5) + 
  g1^2*MatMul[ht, Adj[Yt], Yt]*((-2*htc)/3 + 12*htc*Zeta[3]) + 
  g1^2*MatMul[Yt, Adj[Yt], ht]*((-2*htc)/3 + 12*htc*Zeta[3]) + 
  g1^4*M1*Yt*((1598*htc)/25 + (988*htc*Zeta[3])/75) + 
  MatMul[ht, Adj[Yt], Yt, Adj[Yt], Yt]*(12*htc + 24*htc*Zeta[3]) + 
  MatMul[Yt, Adj[Yt], ht, Adj[Yt], Yt]*(12*htc + 24*htc*Zeta[3]) + 
  MatMul[Yt, Adj[Yt], Yt, Adj[Yt], ht]*(12*htc + 24*htc*Zeta[3]) + 
  g1^2*g3^2*M1*Yt*((16*htc)/15 + (448*htc*Zeta[3])/15) + 
  g1^2*g3^2*M3*Yt*((16*htc)/15 + (448*htc*Zeta[3])/15) + 
  g1^2*g2^2*ht*((-134*htc)/5 + (156*htc*Zeta[3])/5) + 
  g2^2*M2*MatMul[Yt, Adj[Yb], Yb]*(-18*htc + 36*htc*Zeta[3]) + 
  g2^2*M2*MatMul[Yt, Adj[Yt], Yt]*(-18*htc + 36*htc*Zeta[3]) + 
  g2^4*M2*Yt*(174*htc + 36*htc*Zeta[3]) + 
  g2^2*g3^2*ht*(-176*htc + 192*htc*Zeta[3]) + 
  g3^4*M3*Yt*((-64*htc)/3 + (2176*htc*Zeta[3])/3) + 
  g3^4*MatMul[Yt, Adj[Yt]]*((32*mh2)/3 - (1088*mh2*Zeta[3])/3) + 
  g2^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]*(36*mh2 - 72*mh2*Zeta[3]) + 
  g2^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*(18*mh1 + 18*mh2 - 36*mh1*Zeta[3] - 
    36*mh2*Zeta[3]) + g1^2*g3^2*MatMul[Yt, Adj[Yt]]*
   ((-16*mh2)/15 - (448*mh2*Zeta[3])/15) + g2^4*MatMul[Yt, Adj[Yt]]*
   (-12*mh1 - 99*mh2 - 12*trace[ml] - 36*trace[mq] - 18*mh2*Zeta[3]) + 
  g1^4*MatMul[Yt, Adj[Yt]]*((12*mh1)/25 - (787*mh2)/25 + (8*trace[mb])/25 + 
    (24*trace[me])/25 + (12*trace[ml])/25 + (4*trace[mq])/25 + 
    (32*trace[mt])/25 - (494*mh2*Zeta[3])/75) + 
  g1^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*((38*mh1)/15 + (38*mh2)/15 + 
    (36*mh1*Zeta[3])/5 + (36*mh2*Zeta[3])/5) + 
  g1^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]*((-4*mh2)/3 + 24*mh2*Zeta[3]) + 
  g1^2*g2^2*MatMul[Yt, Adj[Yt]]*((-134*mh2)/5 + (156*mh2*Zeta[3])/5) + 
  MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt]]*(36*mh2 + 72*mh2*Zeta[3]) + 
  g2^2*g3^2*MatMul[Yt, Adj[Yt]]*(-176*mh2 + 192*mh2*Zeta[3]) + 
  g2^2*Yt*(54*htc*trace[ht, Adj[Yt]] - 108*htc*trace[ht, Adj[Yt]]*Zeta[3]) + 
  g1^2*Yt*(14*htc*trace[ht, Adj[Yt]] + (84*htc*trace[ht, Adj[Yt]]*Zeta[3])/
     5) + g3^2*Yt*(-32*htc*trace[ht, Adj[Yt]] + 192*htc*trace[ht, Adj[Yt]]*
     Zeta[3]) + g3^2*M3*MatMul[Yt, Adj[Yt]]*(32*trace[ht, Adj[Yt]] + 
    32*trace[htc, Yt] - 192*trace[ht, Adj[Yt]]*Zeta[3] - 
    192*trace[htc, Yt]*Zeta[3]) + g2^2*MatMul[ht, Adj[Yt]]*
   (54*trace[htc, Yt] - 108*trace[htc, Yt]*Zeta[3]) + 
  g1^2*M1*MatMul[Yt, Adj[Yt]]*(-14*trace[ht, Adj[Yt]] - 14*trace[htc, Yt] - 
    (84*trace[ht, Adj[Yt]]*Zeta[3])/5 - (84*trace[htc, Yt]*Zeta[3])/5) + 
  g1^2*MatMul[ht, Adj[Yt]]*(14*trace[htc, Yt] + (84*trace[htc, Yt]*Zeta[3])/
     5) + g2^2*M2*MatMul[Yt, Adj[Yt]]*(-54*trace[ht, Adj[Yt]] - 
    54*trace[htc, Yt] + 108*trace[ht, Adj[Yt]]*Zeta[3] + 
    108*trace[htc, Yt]*Zeta[3]) + g3^2*MatMul[ht, Adj[Yt]]*
   (-32*trace[htc, Yt] + 192*trace[htc, Yt]*Zeta[3]) + 
  g2^2*M2^2*MatMul[Yt, Adj[Yt]]*(108*trace[Adj[Yt], Yt] - 
    216*trace[Adj[Yt], Yt]*Zeta[3]) + g3^2*M3*MatMul[ht, Adj[Yt]]*
   (32*trace[Adj[Yt], Yt] - 192*trace[Adj[Yt], Yt]*Zeta[3]) + 
  g2^2*mq*Yt*Adj[Yt]*(54*trace[Adj[Yt], Yt] - 108*trace[Adj[Yt], Yt]*
     Zeta[3]) + 2*g2^2*mt*MatMul[Yt, Adj[Yt]]*(27*trace[Adj[Yt], Yt] - 
    54*trace[Adj[Yt], Yt]*Zeta[3]) + g1^2*M1*MatMul[ht, Adj[Yt]]*
   (-14*trace[Adj[Yt], Yt] - (84*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  2*g1^2*mt*MatMul[Yt, Adj[Yt]]*(7*trace[Adj[Yt], Yt] + 
    (42*trace[Adj[Yt], Yt]*Zeta[3])/5) + g1^2*mq*Yt*Adj[Yt]*
   (14*trace[Adj[Yt], Yt] + (84*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  g1^2*M1^2*MatMul[Yt, Adj[Yt]]*(28*trace[Adj[Yt], Yt] + 
    (168*trace[Adj[Yt], Yt]*Zeta[3])/5) + 2*g3^2*mt*MatMul[Yt, Adj[Yt]]*
   (-16*trace[Adj[Yt], Yt] + 96*trace[Adj[Yt], Yt]*Zeta[3]) + 
  g2^2*M2*MatMul[ht, Adj[Yt]]*(-54*trace[Adj[Yt], Yt] + 
    108*trace[Adj[Yt], Yt]*Zeta[3]) + g3^2*mq*Yt*Adj[Yt]*
   (-32*trace[Adj[Yt], Yt] + 192*trace[Adj[Yt], Yt]*Zeta[3]) + 
  g3^2*M3^2*MatMul[Yt, Adj[Yt]]*(-64*trace[Adj[Yt], Yt] + 
    384*trace[Adj[Yt], Yt]*Zeta[3]) + 
  g3^2*M3*Yt*(32*htc*trace[Adj[Yt], Yt] - 192*htc*trace[Adj[Yt], Yt]*
     Zeta[3]) + g2^2*ht*(54*htc*trace[Adj[Yt], Yt] - 
    108*htc*trace[Adj[Yt], Yt]*Zeta[3]) + 
  g1^2*M1*Yt*(-14*htc*trace[Adj[Yt], Yt] - 
    (84*htc*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  g1^2*ht*(14*htc*trace[Adj[Yt], Yt] + (84*htc*trace[Adj[Yt], Yt]*Zeta[3])/
     5) + g2^2*M2*Yt*(-54*htc*trace[Adj[Yt], Yt] + 
    108*htc*trace[Adj[Yt], Yt]*Zeta[3]) + 
  g3^2*ht*(-32*htc*trace[Adj[Yt], Yt] + 192*htc*trace[Adj[Yt], Yt]*Zeta[3]) + 
  g2^2*MatMul[Yt, Adj[Yt]]*(54*trace[ht, htc] + 108*mh2*trace[Adj[Yt], Yt] + 
    54*trace[Yt, Adj[Yt], mt] + 54*trace[Adj[Yt], Yt, mq] - 
    108*trace[ht, htc]*Zeta[3] - 216*mh2*trace[Adj[Yt], Yt]*Zeta[3] - 
    108*trace[Yt, Adj[Yt], mt]*Zeta[3] - 108*trace[Adj[Yt], Yt, mq]*
     Zeta[3]) + g1^2*MatMul[Yt, Adj[Yt]]*(14*trace[ht, htc] + 
    28*mh2*trace[Adj[Yt], Yt] + 14*trace[Yt, Adj[Yt], mt] + 
    14*trace[Adj[Yt], Yt, mq] + (84*trace[ht, htc]*Zeta[3])/5 + 
    (168*mh2*trace[Adj[Yt], Yt]*Zeta[3])/5 + 
    (84*trace[Yt, Adj[Yt], mt]*Zeta[3])/5 + 
    (84*trace[Adj[Yt], Yt, mq]*Zeta[3])/5) + g3^2*MatMul[Yt, Adj[Yt]]*
   (-32*trace[ht, htc] - 64*mh2*trace[Adj[Yt], Yt] - 
    32*trace[Yt, Adj[Yt], mt] - 32*trace[Adj[Yt], Yt, mq] + 
    192*trace[ht, htc]*Zeta[3] + 384*mh2*trace[Adj[Yt], Yt]*Zeta[3] + 
    192*trace[Yt, Adj[Yt], mt]*Zeta[3] + 192*trace[Adj[Yt], Yt, mq]*Zeta[3])}

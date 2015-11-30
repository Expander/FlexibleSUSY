{4*ht*htc - (32*g1^2*M1^2)/15 - (32*g3^2*M3^2)/3 + 4*mq*Yt*Ytc + 
  4*mh2*MatMul[Yt, Ytc] + 4*mt*MatMul[Yt, Ytc], 
 (-4*g1^2*ht*htc)/5 + 12*g2^2*ht*htc + (3424*g1^4*M1^2)/75 + 
  (512*g1^2*g3^2*M1^2)/45 + (512*g1^2*g3^2*M1*M3)/45 + 
  (512*g1^2*g3^2*M3^2)/45 - (128*g3^4*M3^2)/3 + (4*g1^2*htc*M1*Yt)/5 - 
  12*g2^2*htc*M2*Yt - (4*g1^2*mq*Yt*Ytc)/5 + 12*g2^2*mq*Yt*Ytc + 
  (4*g1^2*M1*MatMul[ht, Ytc])/5 - 12*g2^2*M2*MatMul[ht, Ytc] - 
  4*mb*MatMul[Yb, Ytc]*MatMul[Yt, Ybc] - (8*g1^2*M1^2*MatMul[Yt, Ytc])/5 + 
  24*g2^2*M2^2*MatMul[Yt, Ytc] - (4*g1^2*mh2*MatMul[Yt, Ytc])/5 + 
  12*g2^2*mh2*MatMul[Yt, Ytc] - (4*g1^2*mt*MatMul[Yt, Ytc])/5 + 
  12*g2^2*mt*MatMul[Yt, Ytc] - 4*mt*MatMul[Yt, Ytc]^2 - 
  4*htc*MatMul[ht, Ybc, Yb] - 4*htc*MatMul[ht, Ytc, Yt] - 
  4*mq*Yt*MatMul[Ybc, Yb, Ytc] - 4*htc*MatMul[Yt, Ybc, hb] - 
  4*mq*Ytc*MatMul[Yt, Ybc, Yb] - 4*htc*MatMul[Yt, Ytc, ht] - 
  4*mq*Ytc*MatMul[Yt, Ytc, Yt] - 4*mq*Yt*MatMul[Ytc, Yt, Ytc] - 
  4*MatMul[ht, hbc, Yb, Ytc] - 4*MatMul[ht, htc, Yt, Ytc] - 
  4*MatMul[Yt, hbc, hb, Ytc] - 4*MatMul[Yt, htc, ht, Ytc] + 
  (-4*mh1 - 4*mh2)*MatMul[Yt, Ybc, Yb, Ytc] - 4*mt*MatMul[Yt, Ybc, Yb, Ytc] - 
  8*mh2*MatMul[Yt, Ytc, Yt, Ytc] - 4*mt*MatMul[Yt, Ytc, Yt, Ytc] + 
  g1^4*((16*mh1)/25 + (16*mh2)/25 + (32*trace[mb])/75 + (32*trace[me])/25 + 
    (16*trace[ml])/25 + (16*trace[mq])/75 + (128*trace[mt])/75) + 
  g3^4*((16*trace[mb])/3 + (32*trace[mq])/3 + (16*trace[mt])/3) - 
  12*htc*Yt*trace[ht, Ytc] - 12*MatMul[ht, Ytc]*trace[htc, Yt] - 
  12*ht*htc*trace[Ytc, Yt] - 12*mq*Yt*Ytc*trace[Ytc, Yt] - 
  12*mt*MatMul[Yt, Ytc]*trace[Ytc, Yt] + MatMul[Yt, Ytc]*
   (-12*trace[ht, htc] - 24*mh2*trace[Ytc, Yt] - 12*trace[Yt, Ytc, mt] - 
    12*trace[Ytc, Yt, mq]), (128*g3^2*mb*MatMul[Yb, Ytc]*MatMul[Yt, Ybc])/3 + 
  (128*g3^2*mt*MatMul[Yt, Ytc]^2)/3 + (128*g3^2*htc*MatMul[ht, Ybc, Yb])/3 + 
  (128*g3^2*htc*MatMul[ht, Ytc, Yt])/3 + 
  (128*g3^2*mq*Yt*MatMul[Ybc, Yb, Ytc])/3 + 
  (128*g3^2*htc*MatMul[Yt, Ybc, hb])/3 - 
  (128*g3^2*htc*M3*MatMul[Yt, Ybc, Yb])/3 + 
  (128*g3^2*mq*Ytc*MatMul[Yt, Ybc, Yb])/3 + 12*mq*MatMul[Ybc, Yb, Ytc]*
   MatMul[Yt, Ybc, Yb] + (128*g3^2*htc*MatMul[Yt, Ytc, ht])/3 - 
  (128*g3^2*htc*M3*MatMul[Yt, Ytc, Yt])/3 + 
  (128*g3^2*mq*Ytc*MatMul[Yt, Ytc, Yt])/3 - 4*mq*MatMul[Ybc, Yb, Ytc]*
   MatMul[Yt, Ytc, Yt] + (128*g3^2*mq*Yt*MatMul[Ytc, Yt, Ytc])/3 - 
  4*mq*MatMul[Yt, Ybc, Yb]*MatMul[Ytc, Yt, Ytc] + 
  (128*g3^2*MatMul[ht, hbc, Yb, Ytc])/3 + (128*g3^2*MatMul[ht, htc, Yt, Ytc])/
   3 - (128*g3^2*M3*MatMul[ht, Ybc, Yb, Ytc])/3 - 
  (128*g3^2*M3*MatMul[ht, Ytc, Yt, Ytc])/3 + 12*mb*MatMul[Yt, Ybc]*
   MatMul[Yb, Ybc, Yb, Ytc] - 4*mb*MatMul[Yt, Ybc]*MatMul[Yb, Ytc, Yt, Ytc] + 
  (128*g3^2*MatMul[Yt, hbc, hb, Ytc])/3 - 
  (128*g3^2*M3*MatMul[Yt, hbc, Yb, Ytc])/3 + 
  (128*g3^2*MatMul[Yt, htc, ht, Ytc])/3 - 
  (128*g3^2*M3*MatMul[Yt, htc, Yt, Ytc])/3 - 
  (128*g3^2*M3*MatMul[Yt, Ybc, hb, Ytc])/3 + 12*mb*MatMul[Yb, Ytc]*
   MatMul[Yt, Ybc, Yb, Ybc] + (256*g3^2*M3^2*MatMul[Yt, Ybc, Yb, Ytc])/3 + 
  g3^2*((128*mh1)/3 + (128*mh2)/3)*MatMul[Yt, Ybc, Yb, Ytc] + 
  (128*g3^2*mt*MatMul[Yt, Ybc, Yb, Ytc])/3 - 8*mt*MatMul[Yt, Ytc]*
   MatMul[Yt, Ybc, Yb, Ytc] - (128*g3^2*M3*MatMul[Yt, Ytc, ht, Ytc])/3 - 
  4*mb*MatMul[Yb, Ytc]*MatMul[Yt, Ytc, Yt, Ybc] + 
  (256*g3^2*M3^2*MatMul[Yt, Ytc, Yt, Ytc])/3 + 
  (256*g3^2*mh2*MatMul[Yt, Ytc, Yt, Ytc])/3 + 
  (128*g3^2*mt*MatMul[Yt, Ytc, Yt, Ytc])/3 + 
  12*htc*MatMul[ht, Ybc, Yb, Ybc, Yb] - 4*htc*MatMul[ht, Ybc, Yb, Ytc, Yt] - 
  4*htc*MatMul[ht, Ytc, Yt, Ybc, Yb] + 
  12*mq*Yt*MatMul[Ybc, Yb, Ybc, Yb, Ytc] - 
  4*mq*Yt*MatMul[Ybc, Yb, Ytc, Yt, Ytc] + 
  12*htc*MatMul[Yt, Ybc, hb, Ybc, Yb] - 4*htc*MatMul[Yt, Ybc, hb, Ytc, Yt] + 
  12*htc*MatMul[Yt, Ybc, Yb, Ybc, hb] + 
  12*mq*Ytc*MatMul[Yt, Ybc, Yb, Ybc, Yb] - 
  4*htc*MatMul[Yt, Ybc, Yb, Ytc, ht] - 
  4*mq*Ytc*MatMul[Yt, Ybc, Yb, Ytc, Yt] - 
  4*htc*MatMul[Yt, Ytc, ht, Ybc, Yb] - 4*htc*MatMul[Yt, Ytc, Yt, Ybc, hb] - 
  4*mq*Ytc*MatMul[Yt, Ytc, Yt, Ybc, Yb] - 
  4*mq*Yt*MatMul[Ytc, Yt, Ybc, Yb, Ytc] + 
  12*MatMul[ht, hbc, Yb, Ybc, Yb, Ytc] - 
  4*MatMul[ht, hbc, Yb, Ytc, Yt, Ytc] - 4*MatMul[ht, htc, Yt, Ybc, Yb, Ytc] + 
  12*MatMul[ht, Ybc, Yb, hbc, Yb, Ytc] - 
  4*MatMul[ht, Ybc, Yb, htc, Yt, Ytc] - 4*MatMul[ht, Ytc, Yt, hbc, Yb, Ytc] + 
  12*MatMul[Yt, hbc, hb, Ybc, Yb, Ytc] - 
  4*MatMul[Yt, hbc, hb, Ytc, Yt, Ytc] + 
  12*MatMul[Yt, hbc, Yb, Ybc, hb, Ytc] - 
  4*MatMul[Yt, hbc, Yb, Ytc, ht, Ytc] - 4*MatMul[Yt, htc, ht, Ybc, Yb, Ytc] - 
  4*MatMul[Yt, htc, Yt, Ybc, hb, Ytc] + 
  12*MatMul[Yt, Ybc, hb, hbc, Yb, Ytc] - 
  4*MatMul[Yt, Ybc, hb, htc, Yt, Ytc] + 
  12*MatMul[Yt, Ybc, Yb, hbc, hb, Ytc] - 
  4*MatMul[Yt, Ybc, Yb, htc, ht, Ytc] + (24*mh1 + 12*mh2)*
   MatMul[Yt, Ybc, Yb, Ybc, Yb, Ytc] + 
  12*mt*MatMul[Yt, Ybc, Yb, Ybc, Yb, Ytc] + 
  (-4*mh1 - 8*mh2)*MatMul[Yt, Ybc, Yb, Ytc, Yt, Ytc] - 
  4*mt*MatMul[Yt, Ybc, Yb, Ytc, Yt, Ytc] - 
  4*MatMul[Yt, Ytc, ht, hbc, Yb, Ytc] - 4*MatMul[Yt, Ytc, Yt, hbc, hb, Ytc] + 
  (-4*mh1 - 8*mh2)*MatMul[Yt, Ytc, Yt, Ybc, Yb, Ytc] - 
  4*mt*MatMul[Yt, Ytc, Yt, Ybc, Yb, Ytc] + 
  g1^6*((-3424*mh1)/375 - (3424*mh2)/375 - (6848*trace[mb])/1125 - 
    (6848*trace[me])/375 - (3424*trace[ml])/375 - (3424*trace[mq])/1125 - 
    (27392*trace[mt])/1125) + g1^4*g3^2*((-256*mh1)/75 - (256*mh2)/75 - 
    (512*trace[mb])/225 - (512*trace[me])/75 - (256*trace[ml])/75 - 
    (256*trace[mq])/225 - (2048*trace[mt])/225) + 
  g1^2*g3^4*((-256*trace[mb])/45 - (512*trace[mq])/45 - (256*trace[mt])/45) + 
  g3^6*((320*trace[mb])/9 + (640*trace[mq])/9 + (320*trace[mt])/9) + 
  MatMul[Yt, hbc, Yb, Ytc]*(24*trace[hb, Ybc] + 8*trace[he, Adj[Ye]] - 
    12*trace[ht, Ytc]) + 12*htc*MatMul[Yt, Ytc, Yt]*trace[ht, Ytc] + 
  12*MatMul[Yt, htc, Yt, Ytc]*trace[ht, Ytc] + 
  MatMul[Yt, Ybc, Yb]*(24*htc*trace[hb, Ybc] + 8*htc*trace[he, Adj[Ye]] - 
    12*htc*trace[ht, Ytc]) + MatMul[ht, Ybc, Yb, Ytc]*
   (24*trace[hbc, Yb] + 8*trace[hec, Ye] - 12*trace[htc, Yt]) + 
  MatMul[Yt, Ybc, hb, Ytc]*(24*trace[hbc, Yb] + 8*trace[hec, Ye] - 
    12*trace[htc, Yt]) + 12*MatMul[ht, Ytc, Yt, Ytc]*trace[htc, Yt] + 
  12*MatMul[Yt, Ytc, ht, Ytc]*trace[htc, Yt] + 
  g1^4*M1*((224*trace[hb, Ybc])/15 + (224*trace[hbc, Yb])/15 + 
    (96*trace[he, Adj[Ye]])/5 + (96*trace[hec, Ye])/5 + 
    (416*trace[ht, Ytc])/15 + (416*trace[htc, Yt])/15) + 
  g3^4*M3*((320*trace[hb, Ybc])/3 + (320*trace[hbc, Yb])/3 + 
    (320*trace[ht, Ytc])/3 + (320*trace[htc, Yt])/3) + 
  g3^4*M3^2*(-320*trace[Ybc, Yb] - 320*trace[Ytc, Yt]) + 
  12*mt*MatMul[Yt, Ytc]^2*trace[Ytc, Yt] + 12*htc*MatMul[ht, Ytc, Yt]*
   trace[Ytc, Yt] + 12*htc*MatMul[Yt, Ytc, ht]*trace[Ytc, Yt] + 
  12*mq*Ytc*MatMul[Yt, Ytc, Yt]*trace[Ytc, Yt] + 
  12*mq*Yt*MatMul[Ytc, Yt, Ytc]*trace[Ytc, Yt] + 
  12*MatMul[ht, htc, Yt, Ytc]*trace[Ytc, Yt] + 12*MatMul[Yt, htc, ht, Ytc]*
   trace[Ytc, Yt] + 12*mt*MatMul[Yt, Ytc, Yt, Ytc]*trace[Ytc, Yt] + 
  g1^4*M1^2*((-224*trace[Ybc, Yb])/5 - (416*trace[Ytc, Yt])/5 - 
    (288*trace[Adj[Ye], Ye])/5) + 2*mt*MatMul[Yt, Ybc, Yb, Ytc]*
   (12*trace[Ybc, Yb] - 6*trace[Ytc, Yt] + 4*trace[Adj[Ye], Ye]) + 
  mb*MatMul[Yb, Ytc]*MatMul[Yt, Ybc]*(24*trace[Ybc, Yb] - 12*trace[Ytc, Yt] + 
    8*trace[Adj[Ye], Ye]) + mq*Yt*MatMul[Ybc, Yb, Ytc]*
   (24*trace[Ybc, Yb] - 12*trace[Ytc, Yt] + 8*trace[Adj[Ye], Ye]) + 
  mq*Ytc*MatMul[Yt, Ybc, Yb]*(24*trace[Ybc, Yb] - 12*trace[Ytc, Yt] + 
    8*trace[Adj[Ye], Ye]) + MatMul[ht, hbc, Yb, Ytc]*
   (24*trace[Ybc, Yb] - 12*trace[Ytc, Yt] + 8*trace[Adj[Ye], Ye]) + 
  MatMul[Yt, hbc, hb, Ytc]*(24*trace[Ybc, Yb] - 12*trace[Ytc, Yt] + 
    8*trace[Adj[Ye], Ye]) + MatMul[ht, Ybc, Yb]*(24*htc*trace[Ybc, Yb] - 
    12*htc*trace[Ytc, Yt] + 8*htc*trace[Adj[Ye], Ye]) + 
  MatMul[Yt, Ybc, hb]*(24*htc*trace[Ybc, Yb] - 12*htc*trace[Ytc, Yt] + 
    8*htc*trace[Adj[Ye], Ye]) + g3^4*(-64*trace[hb, hbc] - 
    64*trace[ht, htc] - 64*mh1*trace[Ybc, Yb] - 64*mh2*trace[Ytc, Yt] - 
    64*trace[Yb, Ybc, mb] - 64*trace[Ybc, Yb, mq] - 64*trace[Yt, Ytc, mt] - 
    64*trace[Ytc, Yt, mq]) + MatMul[Yt, Ytc, Yt, Ytc]*
   (12*trace[ht, htc] + 36*mh2*trace[Ytc, Yt] + 12*trace[Yt, Ytc, mt] + 
    12*trace[Ytc, Yt, mq]) + g1^4*((-224*trace[hb, hbc])/25 - 
    (288*trace[he, hec])/25 - (416*trace[ht, htc])/25 - 
    (224*mh1*trace[Ybc, Yb])/25 - (416*mh2*trace[Ytc, Yt])/25 - 
    (288*mh1*trace[Adj[Ye], Ye])/25 - (224*trace[Yb, Ybc, mb])/25 - 
    (224*trace[Ybc, Yb, mq])/25 - (288*trace[Ye, Adj[Ye], me])/25 - 
    (416*trace[Yt, Ytc, mt])/25 - (416*trace[Ytc, Yt, mq])/25 - 
    (288*trace[Adj[Ye], Ye, ml])/25) + MatMul[Yt, Ybc, Yb, Ytc]*
   (24*trace[hb, hbc] + 8*trace[he, hec] - 12*trace[ht, htc] + 
    48*mh1*trace[Ybc, Yb] + 24*mh2*trace[Ybc, Yb] - 12*mh1*trace[Ytc, Yt] - 
    24*mh2*trace[Ytc, Yt] + 16*mh1*trace[Adj[Ye], Ye] + 
    8*mh2*trace[Adj[Ye], Ye] + 24*trace[Yb, Ybc, mb] + 
    24*trace[Ybc, Yb, mq] + 8*trace[Ye, Adj[Ye], me] - 
    12*trace[Yt, Ytc, mt] - 12*trace[Ytc, Yt, mq] + 
    8*trace[Adj[Ye], Ye, ml]) + Yt*(-72*htc*trace[ht, Ytc]*trace[Ytc, Yt] + 
    24*htc*trace[hb, Ytc, Yt, Ybc] + 24*htc*trace[Ytc, ht, Ybc, Yb] + 
    144*htc*trace[Ytc, ht, Ytc, Yt]) + MatMul[ht, Ytc]*
   (-72*trace[htc, Yt]*trace[Ytc, Yt] + 24*trace[htc, Yt, Ybc, Yb] + 
    24*trace[Ytc, Yt, hbc, Yb] + 144*trace[Ytc, Yt, htc, Yt]) + 
  2*mt*MatMul[Yt, Ytc]*(-18*trace[Ytc, Yt]^2 + 12*trace[Ytc, Yt, Ybc, Yb] + 
    36*trace[Ytc, Yt, Ytc, Yt]) + mq*Yt*Ytc*(-36*trace[Ytc, Yt]^2 + 
    24*trace[Ytc, Yt, Ybc, Yb] + 72*trace[Ytc, Yt, Ytc, Yt]) + 
  ht*(-36*htc*trace[Ytc, Yt]^2 + 24*htc*trace[Ytc, Yt, Ybc, Yb] + 
    72*htc*trace[Ytc, Yt, Ytc, Yt]) + MatMul[Yt, Ytc]*
   (-72*trace[ht, Ytc]*trace[htc, Yt] - 72*trace[ht, htc]*trace[Ytc, Yt] - 
    108*mh2*trace[Ytc, Yt]^2 - 72*trace[Ytc, Yt]*trace[Yt, Ytc, mt] - 
    72*trace[Ytc, Yt]*trace[Ytc, Yt, mq] + 24*trace[hb, htc, Yt, Ybc] + 
    24*trace[hbc, hb, Ytc, Yt] + 24*trace[htc, ht, Ybc, Yb] + 
    144*trace[htc, ht, Ytc, Yt] + 24*trace[Ytc, ht, hbc, Yb] + 
    144*trace[Ytc, ht, htc, Yt] + 24*mh1*trace[Ytc, Yt, Ybc, Yb] + 
    48*mh2*trace[Ytc, Yt, Ybc, Yb] + 216*mh2*trace[Ytc, Yt, Ytc, Yt] + 
    24*trace[Yb, Ytc, Yt, Ybc, mb] + 24*trace[Ybc, Yb, Ytc, Yt, mq] + 
    24*trace[Yt, Ybc, Yb, Ytc, mt] + 144*trace[Yt, Ytc, Yt, Ytc, mt] + 
    24*trace[Ytc, Yt, Ybc, Yb, mq] + 144*trace[Ytc, Yt, Ytc, Yt, mq]) + 
  g3^4*M3^2*MatMul[Yt, Ytc]*(64 - 2176*Zeta[3]) + 
  g2^2*g3^4*M3^2*(720 - 864*Zeta[3]) + g2^2*g3^4*M2*M3*(480 - 576*Zeta[3]) + 
  g3^4*mq*Yt*Ytc*(32/3 - (1088*Zeta[3])/3) + 
  g1^4*g3^2*M1^2*(8576/75 - (8448*Zeta[3])/25) + 
  g1^6*M1^2*(864496/1125 - (38208*Zeta[3])/125) + 
  g2^2*g3^4*M2^2*(288 - 288*Zeta[3]) + g1^4*g3^2*M1*M3*
   (17152/225 - (5632*Zeta[3])/25) + g1^2*g3^4*M3^2*
   (-176/15 - (1056*Zeta[3])/5) + g2^2*g3^2*M2*MatMul[ht, Ytc]*
   (176 - 192*Zeta[3]) + g2^2*g3^2*M3*MatMul[ht, Ytc]*(176 - 192*Zeta[3]) + 
  2*g3^4*mt*MatMul[Yt, Ytc]*(16/3 - (544*Zeta[3])/3) + 
  g1^2*g3^4*M1*M3*(-1376/45 - (704*Zeta[3])/5) + 
  g1^4*g3^2*M3^2*(512/9 - (2816*Zeta[3])/25) + g2^4*M2^2*MatMul[Yt, Ytc]*
   (-474 - 108*Zeta[3]) + g1^4*g2^2*M1^2*(432/5 - (2592*Zeta[3])/25) + 
  g2^2*M2^2*MatMul[Yt, Ybc, Yb, Ytc]*(36 - 72*Zeta[3]) + 
  g2^2*M2^2*MatMul[Yt, Ytc, Yt, Ytc]*(36 - 72*Zeta[3]) + 
  g1^2*g3^4*M1^2*(-32/9 - (352*Zeta[3])/5) + 
  g1^4*g2^2*M1*M2*(288/5 - (1728*Zeta[3])/25) + 
  g1^2*g3^2*M1^2*MatMul[Yt, Ytc]*(-32/15 - (896*Zeta[3])/15) + 
  g1^2*g3^2*M1*M3*MatMul[Yt, Ytc]*(-32/15 - (896*Zeta[3])/15) + 
  g1^2*g3^2*M3^2*MatMul[Yt, Ytc]*(-32/15 - (896*Zeta[3])/15) + 
  g1^4*M1^2*MatMul[Yt, Ytc]*(-4794/25 - (988*Zeta[3])/25) + 
  g2^2*mb*MatMul[Yb, Ytc]*MatMul[Yt, Ybc]*(18 - 36*Zeta[3]) + 
  g2^2*mt*MatMul[Yt, Ytc]^2*(18 - 36*Zeta[3]) + 
  g2^2*mq*Yt*MatMul[Ybc, Yb, Ytc]*(18 - 36*Zeta[3]) + 
  g2^2*mq*Ytc*MatMul[Yt, Ybc, Yb]*(18 - 36*Zeta[3]) + 
  g2^2*mq*Ytc*MatMul[Yt, Ytc, Yt]*(18 - 36*Zeta[3]) + 
  g2^2*mq*Yt*MatMul[Ytc, Yt, Ytc]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[ht, hbc, Yb, Ytc]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[ht, htc, Yt, Ytc]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[Yt, hbc, hb, Ytc]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[Yt, htc, ht, Ytc]*(18 - 36*Zeta[3]) + 
  g1^4*g2^2*M2^2*(864/25 - (864*Zeta[3])/25) + g1^2*g2^2*M1*MatMul[ht, Ytc]*
   (134/5 - (156*Zeta[3])/5) + g1^2*g2^2*M2*MatMul[ht, Ytc]*
   (134/5 - (156*Zeta[3])/5) + g1^2*g3^2*mq*Yt*Ytc*
   (-16/15 - (448*Zeta[3])/15) + g2^4*mq*Yt*Ytc*(-87 - 18*Zeta[3]) + 
  2*g2^2*mt*MatMul[Yt, Ybc, Yb, Ytc]*(9 - 18*Zeta[3]) + 
  2*g2^2*mt*MatMul[Yt, Ytc, Yt, Ytc]*(9 - 18*Zeta[3]) + 
  2*g1^2*g3^2*mt*MatMul[Yt, Ytc]*(-8/15 - (224*Zeta[3])/15) + 
  g1^2*M1*MatMul[ht, Ytc, Yt, Ytc]*(2/3 - 12*Zeta[3]) + 
  g1^2*M1*MatMul[Yt, htc, Yt, Ytc]*(2/3 - 12*Zeta[3]) + 
  g1^2*M1*MatMul[Yt, Ytc, ht, Ytc]*(2/3 - 12*Zeta[3]) + 
  2*g2^4*mt*MatMul[Yt, Ytc]*(-87/2 - 9*Zeta[3]) + 
  g1^2*M1*MatMul[ht, Ybc, Yb, Ytc]*(-38/15 - (36*Zeta[3])/5) + 
  g1^2*M1*MatMul[Yt, hbc, Yb, Ytc]*(-38/15 - (36*Zeta[3])/5) + 
  g1^2*M1*MatMul[Yt, Ybc, hb, Ytc]*(-38/15 - (36*Zeta[3])/5) + 
  g1^4*mq*Yt*Ytc*(-799/25 - (494*Zeta[3])/75) + 
  2*g1^4*mt*MatMul[Yt, Ytc]*(-799/50 - (247*Zeta[3])/75) + 
  2*g1^2*mt*MatMul[Yt, Ybc, Yb, Ytc]*(19/15 + (18*Zeta[3])/5) + 
  2*g1^2*mt*MatMul[Yt, Ytc, Yt, Ytc]*(-1/3 + 6*Zeta[3]) + 
  g1^2*mb*MatMul[Yb, Ytc]*MatMul[Yt, Ybc]*(38/15 + (36*Zeta[3])/5) + 
  g1^2*mq*Yt*MatMul[Ybc, Yb, Ytc]*(38/15 + (36*Zeta[3])/5) + 
  g1^2*mq*Ytc*MatMul[Yt, Ybc, Yb]*(38/15 + (36*Zeta[3])/5) + 
  g1^2*MatMul[ht, hbc, Yb, Ytc]*(38/15 + (36*Zeta[3])/5) + 
  g1^2*MatMul[Yt, hbc, hb, Ytc]*(38/15 + (36*Zeta[3])/5) + 
  g1^2*mt*MatMul[Yt, Ytc]^2*(-2/3 + 12*Zeta[3]) + 
  g1^2*mq*Ytc*MatMul[Yt, Ytc, Yt]*(-2/3 + 12*Zeta[3]) + 
  g1^2*mq*Yt*MatMul[Ytc, Yt, Ytc]*(-2/3 + 12*Zeta[3]) + 
  g1^2*MatMul[ht, htc, Yt, Ytc]*(-2/3 + 12*Zeta[3]) + 
  g1^2*MatMul[Yt, htc, ht, Ytc]*(-2/3 + 12*Zeta[3]) + 
  2*mt*MatMul[Yt, Ytc, Yt, Ytc, Yt, Ytc]*(6 + 12*Zeta[3]) + 
  g1^4*M1*MatMul[ht, Ytc]*(1598/25 + (988*Zeta[3])/75) + 
  g1^2*M1^2*MatMul[Yt, Ybc, Yb, Ytc]*(76/15 + (72*Zeta[3])/5) + 
  2*g1^2*g2^2*mt*MatMul[Yt, Ytc]*(-67/5 + (78*Zeta[3])/5) + 
  g1^2*M1^2*MatMul[Yt, Ytc, Yt, Ytc]*(-4/3 + 24*Zeta[3]) + 
  mq*MatMul[Yt, Ytc, Yt]*MatMul[Ytc, Yt, Ytc]*(12 + 24*Zeta[3]) + 
  2*mt*MatMul[Yt, Ytc]*MatMul[Yt, Ytc, Yt, Ytc]*(12 + 24*Zeta[3]) + 
  mq*Ytc*MatMul[Yt, Ytc, Yt, Ytc, Yt]*(12 + 24*Zeta[3]) + 
  mq*Yt*MatMul[Ytc, Yt, Ytc, Yt, Ytc]*(12 + 24*Zeta[3]) + 
  MatMul[ht, htc, Yt, Ytc, Yt, Ytc]*(12 + 24*Zeta[3]) + 
  MatMul[ht, Ytc, Yt, htc, Yt, Ytc]*(12 + 24*Zeta[3]) + 
  MatMul[Yt, htc, ht, Ytc, Yt, Ytc]*(12 + 24*Zeta[3]) + 
  MatMul[Yt, htc, Yt, Ytc, ht, Ytc]*(12 + 24*Zeta[3]) + 
  MatMul[Yt, Ytc, ht, htc, Yt, Ytc]*(12 + 24*Zeta[3]) + 
  MatMul[Yt, Ytc, Yt, htc, ht, Ytc]*(12 + 24*Zeta[3]) + 
  g1^2*g3^2*M1*MatMul[ht, Ytc]*(16/15 + (448*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[ht, Ytc]*(16/15 + (448*Zeta[3])/15) + 
  g1^2*g2^2*mq*Yt*Ytc*(-134/5 + (156*Zeta[3])/5) + 
  g2^2*M2*MatMul[ht, Ybc, Yb, Ytc]*(-18 + 36*Zeta[3]) + 
  g2^2*M2*MatMul[ht, Ytc, Yt, Ytc]*(-18 + 36*Zeta[3]) + 
  g2^2*M2*MatMul[Yt, hbc, Yb, Ytc]*(-18 + 36*Zeta[3]) + 
  g2^2*M2*MatMul[Yt, htc, Yt, Ytc]*(-18 + 36*Zeta[3]) + 
  g2^2*M2*MatMul[Yt, Ybc, hb, Ytc]*(-18 + 36*Zeta[3]) + 
  g2^2*M2*MatMul[Yt, Ytc, ht, Ytc]*(-18 + 36*Zeta[3]) + 
  g2^4*M2*MatMul[ht, Ytc]*(174 + 36*Zeta[3]) + g1^2*g2^2*M1^2*MatMul[Yt, Ytc]*
   (-268/5 + (312*Zeta[3])/5) + g1^2*g2^2*M1*M2*MatMul[Yt, Ytc]*
   (-268/5 + (312*Zeta[3])/5) + g1^2*g2^2*M2^2*MatMul[Yt, Ytc]*
   (-268/5 + (312*Zeta[3])/5) + 2*g2^2*g3^2*mt*MatMul[Yt, Ytc]*
   (-88 + 96*Zeta[3]) + g2^2*g3^2*mq*Yt*Ytc*(-176 + 192*Zeta[3]) + 
  g2^2*g3^2*M2^2*MatMul[Yt, Ytc]*(-352 + 384*Zeta[3]) + 
  g2^2*g3^2*M2*M3*MatMul[Yt, Ytc]*(-352 + 384*Zeta[3]) + 
  g2^2*g3^2*M3^2*MatMul[Yt, Ytc]*(-352 + 384*Zeta[3]) + 
  g3^4*M3*MatMul[ht, Ytc]*(-64/3 + (2176*Zeta[3])/3) + 
  g3^6*M3^2*(20512/9 + 7680*Zeta[3]) + 
  g3^4*ht*((32*htc)/3 - (1088*htc*Zeta[3])/3) + 
  g2^2*g3^2*M2*Yt*(176*htc - 192*htc*Zeta[3]) + 
  g2^2*g3^2*M3*Yt*(176*htc - 192*htc*Zeta[3]) + 
  g2^2*MatMul[ht, Ybc, Yb]*(18*htc - 36*htc*Zeta[3]) + 
  g2^2*MatMul[ht, Ytc, Yt]*(18*htc - 36*htc*Zeta[3]) + 
  g2^2*MatMul[Yt, Ybc, hb]*(18*htc - 36*htc*Zeta[3]) + 
  g2^2*MatMul[Yt, Ytc, ht]*(18*htc - 36*htc*Zeta[3]) + 
  g1^2*g2^2*M1*Yt*((134*htc)/5 - (156*htc*Zeta[3])/5) + 
  g1^2*g2^2*M2*Yt*((134*htc)/5 - (156*htc*Zeta[3])/5) + 
  g1^2*g3^2*ht*((-16*htc)/15 - (448*htc*Zeta[3])/15) + 
  g2^4*ht*(-87*htc - 18*htc*Zeta[3]) + g1^2*M1*MatMul[Yt, Ytc, Yt]*
   ((2*htc)/3 - 12*htc*Zeta[3]) + g1^2*M1*MatMul[Yt, Ybc, Yb]*
   ((-38*htc)/15 - (36*htc*Zeta[3])/5) + 
  g1^4*ht*((-799*htc)/25 - (494*htc*Zeta[3])/75) + 
  g1^2*MatMul[ht, Ybc, Yb]*((38*htc)/15 + (36*htc*Zeta[3])/5) + 
  g1^2*MatMul[Yt, Ybc, hb]*((38*htc)/15 + (36*htc*Zeta[3])/5) + 
  g1^2*MatMul[ht, Ytc, Yt]*((-2*htc)/3 + 12*htc*Zeta[3]) + 
  g1^2*MatMul[Yt, Ytc, ht]*((-2*htc)/3 + 12*htc*Zeta[3]) + 
  g1^4*M1*Yt*((1598*htc)/25 + (988*htc*Zeta[3])/75) + 
  MatMul[ht, Ytc, Yt, Ytc, Yt]*(12*htc + 24*htc*Zeta[3]) + 
  MatMul[Yt, Ytc, ht, Ytc, Yt]*(12*htc + 24*htc*Zeta[3]) + 
  MatMul[Yt, Ytc, Yt, Ytc, ht]*(12*htc + 24*htc*Zeta[3]) + 
  g1^2*g3^2*M1*Yt*((16*htc)/15 + (448*htc*Zeta[3])/15) + 
  g1^2*g3^2*M3*Yt*((16*htc)/15 + (448*htc*Zeta[3])/15) + 
  g1^2*g2^2*ht*((-134*htc)/5 + (156*htc*Zeta[3])/5) + 
  g2^2*M2*MatMul[Yt, Ybc, Yb]*(-18*htc + 36*htc*Zeta[3]) + 
  g2^2*M2*MatMul[Yt, Ytc, Yt]*(-18*htc + 36*htc*Zeta[3]) + 
  g2^4*M2*Yt*(174*htc + 36*htc*Zeta[3]) + 
  g2^2*g3^2*ht*(-176*htc + 192*htc*Zeta[3]) + 
  g3^4*M3*Yt*((-64*htc)/3 + (2176*htc*Zeta[3])/3) + 
  g3^4*MatMul[Yt, Ytc]*((32*mh2)/3 - (1088*mh2*Zeta[3])/3) + 
  g2^2*MatMul[Yt, Ytc, Yt, Ytc]*(36*mh2 - 72*mh2*Zeta[3]) + 
  g2^2*MatMul[Yt, Ybc, Yb, Ytc]*(18*mh1 + 18*mh2 - 36*mh1*Zeta[3] - 
    36*mh2*Zeta[3]) + g1^2*g3^2*MatMul[Yt, Ytc]*
   ((-16*mh2)/15 - (448*mh2*Zeta[3])/15) + g2^4*MatMul[Yt, Ytc]*
   (-12*mh1 - 99*mh2 - 12*trace[ml] - 36*trace[mq] - 18*mh2*Zeta[3]) + 
  g1^4*MatMul[Yt, Ytc]*((12*mh1)/25 - (787*mh2)/25 + (8*trace[mb])/25 + 
    (24*trace[me])/25 + (12*trace[ml])/25 + (4*trace[mq])/25 + 
    (32*trace[mt])/25 - (494*mh2*Zeta[3])/75) + 
  g1^2*MatMul[Yt, Ybc, Yb, Ytc]*((38*mh1)/15 + (38*mh2)/15 + 
    (36*mh1*Zeta[3])/5 + (36*mh2*Zeta[3])/5) + g1^2*MatMul[Yt, Ytc, Yt, Ytc]*
   ((-4*mh2)/3 + 24*mh2*Zeta[3]) + g1^2*g2^2*MatMul[Yt, Ytc]*
   ((-134*mh2)/5 + (156*mh2*Zeta[3])/5) + MatMul[Yt, Ytc, Yt, Ytc, Yt, Ytc]*
   (36*mh2 + 72*mh2*Zeta[3]) + g2^2*g3^2*MatMul[Yt, Ytc]*
   (-176*mh2 + 192*mh2*Zeta[3]) + g2^2*Yt*(54*htc*trace[ht, Ytc] - 
    108*htc*trace[ht, Ytc]*Zeta[3]) + 
  g1^2*Yt*(14*htc*trace[ht, Ytc] + (84*htc*trace[ht, Ytc]*Zeta[3])/5) + 
  g3^2*Yt*(-32*htc*trace[ht, Ytc] + 192*htc*trace[ht, Ytc]*Zeta[3]) + 
  g3^2*M3*MatMul[Yt, Ytc]*(32*trace[ht, Ytc] + 32*trace[htc, Yt] - 
    192*trace[ht, Ytc]*Zeta[3] - 192*trace[htc, Yt]*Zeta[3]) + 
  g2^2*MatMul[ht, Ytc]*(54*trace[htc, Yt] - 108*trace[htc, Yt]*Zeta[3]) + 
  g1^2*M1*MatMul[Yt, Ytc]*(-14*trace[ht, Ytc] - 14*trace[htc, Yt] - 
    (84*trace[ht, Ytc]*Zeta[3])/5 - (84*trace[htc, Yt]*Zeta[3])/5) + 
  g1^2*MatMul[ht, Ytc]*(14*trace[htc, Yt] + (84*trace[htc, Yt]*Zeta[3])/5) + 
  g2^2*M2*MatMul[Yt, Ytc]*(-54*trace[ht, Ytc] - 54*trace[htc, Yt] + 
    108*trace[ht, Ytc]*Zeta[3] + 108*trace[htc, Yt]*Zeta[3]) + 
  g3^2*MatMul[ht, Ytc]*(-32*trace[htc, Yt] + 192*trace[htc, Yt]*Zeta[3]) + 
  g2^2*M2^2*MatMul[Yt, Ytc]*(108*trace[Ytc, Yt] - 
    216*trace[Ytc, Yt]*Zeta[3]) + g3^2*M3*MatMul[ht, Ytc]*
   (32*trace[Ytc, Yt] - 192*trace[Ytc, Yt]*Zeta[3]) + 
  g2^2*mq*Yt*Ytc*(54*trace[Ytc, Yt] - 108*trace[Ytc, Yt]*Zeta[3]) + 
  2*g2^2*mt*MatMul[Yt, Ytc]*(27*trace[Ytc, Yt] - 54*trace[Ytc, Yt]*Zeta[3]) + 
  g1^2*M1*MatMul[ht, Ytc]*(-14*trace[Ytc, Yt] - (84*trace[Ytc, Yt]*Zeta[3])/
     5) + 2*g1^2*mt*MatMul[Yt, Ytc]*(7*trace[Ytc, Yt] + 
    (42*trace[Ytc, Yt]*Zeta[3])/5) + g1^2*mq*Yt*Ytc*
   (14*trace[Ytc, Yt] + (84*trace[Ytc, Yt]*Zeta[3])/5) + 
  g1^2*M1^2*MatMul[Yt, Ytc]*(28*trace[Ytc, Yt] + (168*trace[Ytc, Yt]*Zeta[3])/
     5) + 2*g3^2*mt*MatMul[Yt, Ytc]*(-16*trace[Ytc, Yt] + 
    96*trace[Ytc, Yt]*Zeta[3]) + g2^2*M2*MatMul[ht, Ytc]*
   (-54*trace[Ytc, Yt] + 108*trace[Ytc, Yt]*Zeta[3]) + 
  g3^2*mq*Yt*Ytc*(-32*trace[Ytc, Yt] + 192*trace[Ytc, Yt]*Zeta[3]) + 
  g3^2*M3^2*MatMul[Yt, Ytc]*(-64*trace[Ytc, Yt] + 
    384*trace[Ytc, Yt]*Zeta[3]) + g3^2*M3*Yt*(32*htc*trace[Ytc, Yt] - 
    192*htc*trace[Ytc, Yt]*Zeta[3]) + 
  g2^2*ht*(54*htc*trace[Ytc, Yt] - 108*htc*trace[Ytc, Yt]*Zeta[3]) + 
  g1^2*M1*Yt*(-14*htc*trace[Ytc, Yt] - (84*htc*trace[Ytc, Yt]*Zeta[3])/5) + 
  g1^2*ht*(14*htc*trace[Ytc, Yt] + (84*htc*trace[Ytc, Yt]*Zeta[3])/5) + 
  g2^2*M2*Yt*(-54*htc*trace[Ytc, Yt] + 108*htc*trace[Ytc, Yt]*Zeta[3]) + 
  g3^2*ht*(-32*htc*trace[Ytc, Yt] + 192*htc*trace[Ytc, Yt]*Zeta[3]) + 
  g2^2*MatMul[Yt, Ytc]*(54*trace[ht, htc] + 108*mh2*trace[Ytc, Yt] + 
    54*trace[Yt, Ytc, mt] + 54*trace[Ytc, Yt, mq] - 
    108*trace[ht, htc]*Zeta[3] - 216*mh2*trace[Ytc, Yt]*Zeta[3] - 
    108*trace[Yt, Ytc, mt]*Zeta[3] - 108*trace[Ytc, Yt, mq]*Zeta[3]) + 
  g1^2*MatMul[Yt, Ytc]*(14*trace[ht, htc] + 28*mh2*trace[Ytc, Yt] + 
    14*trace[Yt, Ytc, mt] + 14*trace[Ytc, Yt, mq] + 
    (84*trace[ht, htc]*Zeta[3])/5 + (168*mh2*trace[Ytc, Yt]*Zeta[3])/5 + 
    (84*trace[Yt, Ytc, mt]*Zeta[3])/5 + (84*trace[Ytc, Yt, mq]*Zeta[3])/5) + 
  g3^2*MatMul[Yt, Ytc]*(-32*trace[ht, htc] - 64*mh2*trace[Ytc, Yt] - 
    32*trace[Yt, Ytc, mt] - 32*trace[Ytc, Yt, mq] + 
    192*trace[ht, htc]*Zeta[3] + 384*mh2*trace[Ytc, Yt]*Zeta[3] + 
    192*trace[Yt, Ytc, mt]*Zeta[3] + 192*trace[Ytc, Yt, mq]*Zeta[3])}

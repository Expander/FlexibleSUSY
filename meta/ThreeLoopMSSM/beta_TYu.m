{(-13*g1^2*ht)/15 - 3*g2^2*ht - (16*g3^2*ht)/3 + (26*g1^2*M1*Yt)/15 + 
  6*g2^2*M2*Yt + (32*g3^2*M3*Yt)/3 + MatMul[ht, Ybc, Yb] + 
  5*MatMul[ht, Ytc, Yt] + 2*MatMul[Yt, Ybc, hb] + 4*MatMul[Yt, Ytc, ht] + 
  6*Yt*trace[ht, Ytc] + 3*ht*trace[Ytc, Yt], 
 (2743*g1^4*ht)/450 + g1^2*g2^2*ht + (15*g2^4*ht)/2 + (136*g1^2*g3^2*ht)/45 + 
  8*g2^2*g3^2*ht - (16*g3^4*ht)/9 - (5486*g1^4*M1*Yt)/225 - 
  2*g1^2*g2^2*M1*Yt - (272*g1^2*g3^2*M1*Yt)/45 - 2*g1^2*g2^2*M2*Yt - 
  30*g2^4*M2*Yt - 16*g2^2*g3^2*M2*Yt - (272*g1^2*g3^2*M3*Yt)/45 - 
  16*g2^2*g3^2*M3*Yt + (64*g3^4*M3*Yt)/9 + (2*g1^2*MatMul[ht, Ybc, Yb])/5 + 
  12*g2^2*MatMul[ht, Ytc, Yt] + (4*g1^2*MatMul[Yt, Ybc, hb])/5 - 
  (4*g1^2*M1*MatMul[Yt, Ybc, Yb])/5 + (6*g1^2*MatMul[Yt, Ytc, ht])/5 + 
  6*g2^2*MatMul[Yt, Ytc, ht] - (4*g1^2*M1*MatMul[Yt, Ytc, Yt])/5 - 
  12*g2^2*M2*MatMul[Yt, Ytc, Yt] - 2*MatMul[ht, Ybc, Yb, Ybc, Yb] - 
  4*MatMul[ht, Ybc, Yb, Ytc, Yt] - 6*MatMul[ht, Ytc, Yt, Ytc, Yt] - 
  4*MatMul[Yt, Ybc, hb, Ybc, Yb] - 4*MatMul[Yt, Ybc, hb, Ytc, Yt] - 
  4*MatMul[Yt, Ybc, Yb, Ybc, hb] - 2*MatMul[Yt, Ybc, Yb, Ytc, ht] - 
  8*MatMul[Yt, Ytc, ht, Ytc, Yt] - 6*MatMul[Yt, Ytc, Yt, Ytc, ht] + 
  MatMul[Yt, Ybc, Yb]*(-6*trace[hb, Ybc] - 2*trace[he, Adj[Ye]]) + 
  (8*g1^2*Yt*trace[ht, Ytc])/5 + 32*g3^2*Yt*trace[ht, Ytc] - 
  18*MatMul[Yt, Ytc, Yt]*trace[ht, Ytc] + (4*g1^2*ht*trace[Ytc, Yt])/5 + 
  16*g3^2*ht*trace[Ytc, Yt] - (8*g1^2*M1*Yt*trace[Ytc, Yt])/5 - 
  32*g3^2*M3*Yt*trace[Ytc, Yt] - 15*MatMul[ht, Ytc, Yt]*trace[Ytc, Yt] - 
  12*MatMul[Yt, Ytc, ht]*trace[Ytc, Yt] + MatMul[Yt, Ybc, hb]*
   (-6*trace[Ybc, Yb] - 2*trace[Adj[Ye], Ye]) + 
  MatMul[ht, Ybc, Yb]*(-3*trace[Ybc, Yb] - trace[Adj[Ye], Ye]) + 
  Yt*(-6*trace[hb, Ytc, Yt, Ybc] - 6*trace[Ytc, ht, Ybc, Yb] - 
    36*trace[Ytc, ht, Ytc, Yt]) + ht*(-3*trace[Ytc, Yt, Ybc, Yb] - 
    9*trace[Ytc, Yt, Ytc, Yt]), (-8*g1^2*g2^2*g3^2*ht)/5 + 
  (16*g1^2*g2^2*g3^2*M1*Yt)/5 + (16*g1^2*g2^2*g3^2*M2*Yt)/5 + 
  (16*g1^2*g2^2*g3^2*M3*Yt)/5 - 4*g2^2*g3^2*MatMul[ht, Ybc, Yb] - 
  8*g2^2*g3^2*MatMul[Yt, Ybc, hb] + 8*g2^2*g3^2*M2*MatMul[Yt, Ybc, Yb] + 
  8*g2^2*g3^2*M3*MatMul[Yt, Ybc, Yb] + (64*g3^2*MatMul[ht, Ybc, Yb, Ybc, Yb])/
   3 + (128*g3^2*MatMul[ht, Ybc, Yb, Ytc, Yt])/3 + 
  64*g3^2*MatMul[ht, Ytc, Yt, Ytc, Yt] + 
  (128*g3^2*MatMul[Yt, Ybc, hb, Ybc, Yb])/3 + 
  (128*g3^2*MatMul[Yt, Ybc, hb, Ytc, Yt])/3 + 
  (128*g3^2*MatMul[Yt, Ybc, Yb, Ybc, hb])/3 - 
  (128*g3^2*M3*MatMul[Yt, Ybc, Yb, Ybc, Yb])/3 + 
  (64*g3^2*MatMul[Yt, Ybc, Yb, Ytc, ht])/3 - 
  (128*g3^2*M3*MatMul[Yt, Ybc, Yb, Ytc, Yt])/3 + 
  (20*g1^2*MatMul[Yt, Ytc, ht, Ytc, Yt])/3 + 
  12*g2^2*MatMul[Yt, Ytc, ht, Ytc, Yt] + 
  (256*g3^2*MatMul[Yt, Ytc, ht, Ytc, Yt])/3 + 
  64*g3^2*MatMul[Yt, Ytc, Yt, Ytc, ht] - 
  (20*g1^2*M1*MatMul[Yt, Ytc, Yt, Ytc, Yt])/3 - 
  12*g2^2*M2*MatMul[Yt, Ytc, Yt, Ytc, Yt] - 
  (256*g3^2*M3*MatMul[Yt, Ytc, Yt, Ytc, Yt])/3 + 
  12*MatMul[ht, Ybc, Yb, Ybc, Yb, Ytc, Yt] + 
  4*MatMul[ht, Ybc, Yb, Ytc, Yt, Ybc, Yb] - 
  4*MatMul[ht, Ybc, Yb, Ytc, Yt, Ytc, Yt] + 
  12*MatMul[Yt, Ybc, hb, Ybc, Yb, Ytc, Yt] + 
  8*MatMul[Yt, Ybc, hb, Ytc, Yt, Ybc, Yb] - 
  4*MatMul[Yt, Ybc, hb, Ytc, Yt, Ytc, Yt] + 
  12*MatMul[Yt, Ybc, Yb, Ybc, hb, Ytc, Yt] + 
  6*MatMul[Yt, Ybc, Yb, Ybc, Yb, Ytc, ht] + 
  8*MatMul[Yt, Ybc, Yb, Ytc, ht, Ybc, Yb] - 
  4*MatMul[Yt, Ybc, Yb, Ytc, ht, Ytc, Yt] + 
  8*MatMul[Yt, Ybc, Yb, Ytc, Yt, Ybc, hb] - 
  2*MatMul[Yt, Ybc, Yb, Ytc, Yt, Ytc, ht] + 
  4*MatMul[Yt, Ytc, ht, Ybc, Yb, Ytc, Yt] + 
  4*MatMul[Yt, Ytc, Yt, Ybc, hb, Ytc, Yt] + 
  6*MatMul[Yt, Ytc, Yt, Ybc, Yb, Ytc, ht] + MatMul[Yt, Ybc, Yb, Ybc, Yb]*
   (12*trace[hb, Ybc] + 4*trace[he, Adj[Ye]]) + 
  g2^2*MatMul[Yt, Ybc, Yb]*(36*trace[hb, Ybc] + 12*trace[he, Adj[Ye]]) + 
  MatMul[Yt, Ybc, Yb, Ytc, Yt]*(24*trace[hb, Ybc] + 8*trace[he, Adj[Ye]] - 
    12*trace[ht, Ytc]) + 24*MatMul[Yt, Ytc, Yt, Ytc, Yt]*trace[ht, Ytc] + 
  18*MatMul[ht, Ytc, Yt, Ytc, Yt]*trace[Ytc, Yt] + 
  24*MatMul[Yt, Ytc, ht, Ytc, Yt]*trace[Ytc, Yt] + 
  18*MatMul[Yt, Ytc, Yt, Ytc, ht]*trace[Ytc, Yt] + 
  g2^2*M2*MatMul[Yt, Ybc, Yb]*(-36*trace[Ybc, Yb] - 12*trace[Adj[Ye], Ye]) + 
  MatMul[ht, Ybc, Yb, Ybc, Yb]*(6*trace[Ybc, Yb] + 2*trace[Adj[Ye], Ye]) + 
  MatMul[Yt, Ybc, hb, Ybc, Yb]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[Yt, Ybc, Yb, Ybc, hb]*(12*trace[Ybc, Yb] + 4*trace[Adj[Ye], Ye]) + 
  MatMul[Yt, Ybc, Yb, Ytc, ht]*(12*trace[Ybc, Yb] - 6*trace[Ytc, Yt] + 
    4*trace[Adj[Ye], Ye]) + g2^2*MatMul[ht, Ybc, Yb]*
   (18*trace[Ybc, Yb] + 6*trace[Adj[Ye], Ye]) + 
  MatMul[ht, Ybc, Yb, Ytc, Yt]*(24*trace[Ybc, Yb] - 12*trace[Ytc, Yt] + 
    8*trace[Adj[Ye], Ye]) + MatMul[Yt, Ybc, hb, Ytc, Yt]*
   (24*trace[Ybc, Yb] - 12*trace[Ytc, Yt] + 8*trace[Adj[Ye], Ye]) + 
  g2^2*MatMul[Yt, Ybc, hb]*(36*trace[Ybc, Yb] + 12*trace[Adj[Ye], Ye]) + 
  MatMul[Yt, Ytc, Yt]*(-108*trace[ht, Ytc]*trace[Ytc, Yt] + 
    36*trace[hb, Ytc, Yt, Ybc] + 36*trace[Ytc, ht, Ybc, Yb] + 
    216*trace[Ytc, ht, Ytc, Yt]) + MatMul[Yt, Ytc, ht]*
   (-36*trace[Ytc, Yt]^2 + 24*trace[Ytc, Yt, Ybc, Yb] + 
    72*trace[Ytc, Yt, Ytc, Yt]) + MatMul[ht, Ytc, Yt]*
   (-45*trace[Ytc, Yt]^2 + 30*trace[Ytc, Yt, Ybc, Yb] + 
    90*trace[Ytc, Yt, Ytc, Yt]) + MatMul[Yt, Ybc, Yb]*
   (-36*trace[hb, Ybc]*trace[Ybc, Yb] - 12*trace[he, Adj[Ye]]*
     trace[Ybc, Yb] - 12*trace[hb, Ybc]*trace[Adj[Ye], Ye] - 
    4*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] + 12*trace[hb, Ytc, Yt, Ybc] + 
    72*trace[Ybc, hb, Ybc, Yb] + 12*trace[Ytc, ht, Ybc, Yb] + 
    24*trace[Adj[Ye], he, Adj[Ye], Ye]) + MatMul[ht, Ybc, Yb]*
   (-9*trace[Ybc, Yb]^2 - 6*trace[Ybc, Yb]*trace[Adj[Ye], Ye] - 
    trace[Adj[Ye], Ye]^2 + 18*trace[Ybc, Yb, Ybc, Yb] + 
    6*trace[Ytc, Yt, Ybc, Yb] + 6*trace[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  MatMul[Yt, Ybc, hb]*(-18*trace[Ybc, Yb]^2 - 12*trace[Ybc, Yb]*
     trace[Adj[Ye], Ye] - 2*trace[Adj[Ye], Ye]^2 + 
    36*trace[Ybc, Yb, Ybc, Yb] + 12*trace[Ytc, Yt, Ybc, Yb] + 
    12*trace[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  g3^6*M3*Yt*(-10880/9 - 3840*Zeta[3]) + g2^6*M2*Yt*(-1035 - 1890*Zeta[3]) + 
  g3^4*MatMul[ht, Ytc, Yt]*(40/3 - (1360*Zeta[3])/3) + 
  g3^4*MatMul[Yt, Ytc, ht]*(32/3 - (1088*Zeta[3])/3) + 
  g2^4*g3^2*ht*(140 - 216*Zeta[3]) + g2^2*g3^2*M2*MatMul[Yt, Ytc, Yt]*
   (184 - 192*Zeta[3]) + g2^2*g3^2*M3*MatMul[Yt, Ytc, Yt]*
   (184 - 192*Zeta[3]) + g3^4*MatMul[Yt, Ybc, hb]*(16/3 - (544*Zeta[3])/3) + 
  g2^2*g3^4*ht*(68 - 144*Zeta[3]) + g3^4*MatMul[ht, Ybc, Yb]*
   (8/3 - (272*Zeta[3])/3) + g2^4*MatMul[Yt, Ytc, ht]*(-66 - 72*Zeta[3]) + 
  g2^4*MatMul[Yt, Ybc, hb]*(-45/2 - 63*Zeta[3]) + 
  g2^4*MatMul[ht, Ytc, Yt]*(-393/4 - (99*Zeta[3])/2) + 
  g1^2*g2^2*M1*MatMul[Yt, Ytc, Yt]*(193/5 - (246*Zeta[3])/5) + 
  g1^2*g2^2*M2*MatMul[Yt, Ytc, Yt]*(193/5 - (246*Zeta[3])/5) + 
  g1^4*g3^2*ht*(5308/225 - (1144*Zeta[3])/25) + 
  g2^2*M2*MatMul[Yt, Ybc, Yb, Ybc, Yb]*(6 - 36*Zeta[3]) + 
  g2^2*MatMul[ht, Ybc, Yb, Ytc, Yt]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[Yt, Ybc, hb, Ytc, Yt]*(18 - 36*Zeta[3]) + 
  g1^2*g3^4*ht*(436/45 - (176*Zeta[3])/5) + g2^4*MatMul[ht, Ybc, Yb]*
   (-45/4 - (63*Zeta[3])/2) + g1^2*g3^2*MatMul[ht, Ytc, Yt]*
   (-44/3 - (64*Zeta[3])/3) + g1^6*ht*(352097/6750 - (2587*Zeta[3])/125) + 
  g2^2*MatMul[Yt, Ybc, Yb, Ytc, ht]*(9 - 18*Zeta[3]) + 
  g2^2*MatMul[ht, Ytc, Yt, Ytc, Yt]*(15 - 18*Zeta[3]) + 
  g1^2*g3^2*M1*MatMul[Yt, Ybc, Yb]*(152/15 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[Yt, Ybc, Yb]*(152/15 - (256*Zeta[3])/15) + 
  g1^2*g2^4*ht*(17/2 - (81*Zeta[3])/5) + 
  g1^4*g2^2*ht*(379/50 - (351*Zeta[3])/25) + 
  g1^2*M1*MatMul[Yt, Ybc, Yb, Ytc, Yt]*(-38/15 - (36*Zeta[3])/5) + 
  g1^2*MatMul[Yt, Ytc, Yt, Ytc, ht]*(7 - 6*Zeta[3]) + 
  g1^4*MatMul[ht, Ytc, Yt]*(-2671/60 - (169*Zeta[3])/30) + 
  g1^2*g2^2*M1*MatMul[Yt, Ybc, Yb]*(41/5 - (18*Zeta[3])/5) + 
  g1^2*g2^2*M2*MatMul[Yt, Ybc, Yb]*(41/5 - (18*Zeta[3])/5) + 
  g1^2*MatMul[Yt, Ybc, hb, Ybc, Yb]*(14/15 - (12*Zeta[3])/5) + 
  g1^2*MatMul[Yt, Ybc, Yb, Ybc, hb]*(14/15 - (12*Zeta[3])/5) + 
  g1^4*MatMul[Yt, Ytc, ht]*(-3082/75 - (104*Zeta[3])/75) + 
  g1^2*MatMul[ht, Ybc, Yb, Ybc, Yb]*(7/15 - (6*Zeta[3])/5) + 
  g1^4*M1*MatMul[Yt, Ybc, Yb]*(633/25 - (14*Zeta[3])/15) + 
  g1^4*MatMul[ht, Ybc, Yb]*(-633/100 + (7*Zeta[3])/30) + 
  g1^4*MatMul[Yt, Ybc, hb]*(-633/50 + (7*Zeta[3])/15) + 
  6*MatMul[ht, Ybc, Yb, Ybc, Yb, Ybc, Yb]*Zeta[3] + 
  12*MatMul[Yt, Ybc, hb, Ybc, Yb, Ybc, Yb]*Zeta[3] + 
  12*MatMul[Yt, Ybc, Yb, Ybc, hb, Ybc, Yb]*Zeta[3] + 
  12*MatMul[Yt, Ybc, Yb, Ybc, Yb, Ybc, hb]*Zeta[3] + 
  g1^2*g2^2*MatMul[ht, Ybc, Yb]*(-41/10 + (9*Zeta[3])/5) + 
  g1^2*g3^2*MatMul[Yt, Ytc, ht]*(-416/15 + (32*Zeta[3])/15) + 
  g1^2*M1*MatMul[Yt, Ybc, Yb, Ybc, Yb]*(-14/15 + (12*Zeta[3])/5) + 
  g1^2*g2^2*MatMul[Yt, Ybc, hb]*(-41/5 + (18*Zeta[3])/5) + 
  g1^2*MatMul[Yt, Ybc, Yb, Ytc, ht]*(19/15 + (18*Zeta[3])/5) + 
  g1^2*MatMul[ht, Ytc, Yt, Ytc, Yt]*(3 + 6*Zeta[3]) + 
  g1^2*MatMul[ht, Ybc, Yb, Ytc, Yt]*(38/15 + (36*Zeta[3])/5) + 
  g1^2*MatMul[Yt, Ybc, hb, Ytc, Yt]*(38/15 + (36*Zeta[3])/5) + 
  g1^2*g3^2*MatMul[ht, Ybc, Yb]*(-76/15 + (128*Zeta[3])/15) + 
  g1^4*M1*MatMul[Yt, Ytc, Yt]*(8561/75 + (234*Zeta[3])/25) + 
  g1^2*g3^2*M1*MatMul[Yt, Ytc, Yt]*(424/15 + (64*Zeta[3])/5) + 
  g1^2*g3^2*M3*MatMul[Yt, Ytc, Yt]*(424/15 + (64*Zeta[3])/5) + 
  g1^2*g3^2*MatMul[Yt, Ybc, hb]*(-152/15 + (256*Zeta[3])/15) + 
  g2^2*MatMul[ht, Ybc, Yb, Ybc, Yb]*(-3 + 18*Zeta[3]) + 
  g2^2*MatMul[Yt, Ytc, Yt, Ytc, ht]*(3 + 18*Zeta[3]) + 
  MatMul[Yt, Ytc, Yt, Ytc, Yt, Ytc, ht]*(6 + 24*Zeta[3]) + 
  g1^4*g2^2*M2*Yt*(-379/25 + (702*Zeta[3])/25) + 
  MatMul[ht, Ytc, Yt, Ytc, Yt, Ytc, Yt]*(12 + 30*Zeta[3]) + 
  g1^2*g2^4*M1*Yt*(-17 + (162*Zeta[3])/5) + g1^2*g2^2*MatMul[Yt, Ytc, ht]*
   (-126/5 + (168*Zeta[3])/5) + g2^2*M2*MatMul[Yt, Ybc, Yb, Ytc, Yt]*
   (-18 + 36*Zeta[3]) + g2^2*MatMul[Yt, Ybc, hb, Ybc, Yb]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Yt, Ybc, Yb, Ybc, hb]*(-6 + 36*Zeta[3]) + 
  MatMul[Yt, Ytc, ht, Ytc, Yt, Ytc, Yt]*(12 + 36*Zeta[3]) + 
  MatMul[Yt, Ytc, Yt, Ytc, ht, Ytc, Yt]*(12 + 36*Zeta[3]) + 
  g1^2*g2^2*MatMul[ht, Ytc, Yt]*(-327/10 + (201*Zeta[3])/5) + 
  g1^4*g2^2*M1*Yt*(-758/25 + (1404*Zeta[3])/25) + 
  g1^2*g2^4*M2*Yt*(-34 + (324*Zeta[3])/5) + 
  g1^2*g3^4*M1*Yt*(-872/45 + (352*Zeta[3])/5) + 
  g1^4*g3^2*M3*Yt*(-10616/225 + (2288*Zeta[3])/25) + 
  g2^2*g3^2*MatMul[Yt, Ytc, ht]*(-96 + 96*Zeta[3]) + 
  g1^6*M1*Yt*(-352097/1125 + (15522*Zeta[3])/125) + 
  g2^4*M2*MatMul[Yt, Ybc, Yb]*(45 + 126*Zeta[3]) + 
  g1^2*g3^4*M3*Yt*(-1744/45 + (704*Zeta[3])/5) + 
  g2^4*M2*MatMul[Yt, Ytc, Yt]*(219 + 162*Zeta[3]) + 
  g1^4*g3^2*M1*Yt*(-21232/225 + (4576*Zeta[3])/25) + 
  g2^2*g3^2*MatMul[ht, Ytc, Yt]*(-180 + 192*Zeta[3]) + 
  g2^2*g3^4*M2*Yt*(-136 + 288*Zeta[3]) + g2^6*ht*(345/2 + 315*Zeta[3]) + 
  g3^4*M3*MatMul[Yt, Ybc, Yb]*(-32/3 + (1088*Zeta[3])/3) + 
  g2^4*g3^2*M3*Yt*(-280 + 432*Zeta[3]) + g2^2*g3^4*M3*Yt*
   (-272 + 576*Zeta[3]) + g3^6*ht*(5440/27 + 640*Zeta[3]) + 
  g2^4*g3^2*M2*Yt*(-560 + 864*Zeta[3]) + g3^4*M3*MatMul[Yt, Ytc, Yt]*
   (-32 + 1088*Zeta[3]) + g3^2*MatMul[Yt, Ybc, Yb]*
   (-16*trace[hb, Ybc] + 16*trace[he, Adj[Ye]] + 96*trace[hb, Ybc]*Zeta[3]) + 
  g1^2*MatMul[Yt, Ybc, Yb]*((32*trace[hb, Ybc])/5 - 
    (16*trace[he, Adj[Ye]])/5 - (48*trace[hb, Ybc]*Zeta[3])/5 + 
    (24*trace[he, Adj[Ye]]*Zeta[3])/5) + 
  g2^4*Yt*(-90*trace[hb, Ybc] - 30*trace[he, Adj[Ye]] - 
    (315*trace[ht, Ytc])/2 - 189*trace[ht, Ytc]*Zeta[3]) + 
  g2^2*MatMul[Yt, Ytc, Yt]*(90*trace[ht, Ytc] - 108*trace[ht, Ytc]*Zeta[3]) + 
  g3^4*Yt*((-320*trace[hb, Ybc])/3 - (640*trace[ht, Ytc])/3 - 
    32*trace[ht, Ytc]*Zeta[3]) + g1^4*Yt*((-182*trace[hb, Ybc])/15 - 
    (78*trace[he, Adj[Ye]])/5 - (171*trace[ht, Ytc])/2 - 
    (13*trace[ht, Ytc]*Zeta[3])/5) + g1^2*MatMul[Yt, Ytc, Yt]*
   (18*trace[ht, Ytc] + (36*trace[ht, Ytc]*Zeta[3])/5) + 
  g1^2*g2^2*Yt*((-57*trace[ht, Ytc])/5 + (126*trace[ht, Ytc]*Zeta[3])/5) + 
  g1^2*g3^2*Yt*((-248*trace[ht, Ytc])/3 + (416*trace[ht, Ytc]*Zeta[3])/5) + 
  g2^2*g3^2*Yt*(-264*trace[ht, Ytc] + 288*trace[ht, Ytc]*Zeta[3]) + 
  g3^2*MatMul[Yt, Ytc, Yt]*(-48*trace[ht, Ytc] + 
    288*trace[ht, Ytc]*Zeta[3]) + g3^2*M3*MatMul[Yt, Ybc, Yb]*
   (16*trace[Ybc, Yb] - 16*trace[Adj[Ye], Ye] - 96*trace[Ybc, Yb]*Zeta[3]) + 
  g3^2*MatMul[ht, Ybc, Yb]*(-8*trace[Ybc, Yb] + 8*trace[Adj[Ye], Ye] + 
    48*trace[Ybc, Yb]*Zeta[3]) + g3^2*MatMul[Yt, Ybc, hb]*
   (-16*trace[Ybc, Yb] + 16*trace[Adj[Ye], Ye] + 96*trace[Ybc, Yb]*Zeta[3]) + 
  g3^2*M3*MatMul[Yt, Ytc, Yt]*(48*trace[Ytc, Yt] - 
    288*trace[Ytc, Yt]*Zeta[3]) + g2^2*g3^2*M2*Yt*
   (264*trace[Ytc, Yt] - 288*trace[Ytc, Yt]*Zeta[3]) + 
  g2^2*g3^2*M3*Yt*(264*trace[Ytc, Yt] - 288*trace[Ytc, Yt]*Zeta[3]) + 
  g2^2*MatMul[ht, Ytc, Yt]*(72*trace[Ytc, Yt] - 108*trace[Ytc, Yt]*Zeta[3]) + 
  g2^4*ht*(-45*trace[Ybc, Yb] - (315*trace[Ytc, Yt])/4 - 
    15*trace[Adj[Ye], Ye] - (189*trace[Ytc, Yt]*Zeta[3])/2) + 
  g1^2*g3^2*M1*Yt*((248*trace[Ytc, Yt])/3 - (416*trace[Ytc, Yt]*Zeta[3])/5) + 
  g1^2*g3^2*M3*Yt*((248*trace[Ytc, Yt])/3 - (416*trace[Ytc, Yt]*Zeta[3])/5) + 
  g2^2*MatMul[Yt, Ytc, ht]*(63*trace[Ytc, Yt] - 54*trace[Ytc, Yt]*Zeta[3]) + 
  g1^2*g2^2*M1*Yt*((57*trace[Ytc, Yt])/5 - (126*trace[Ytc, Yt]*Zeta[3])/5) + 
  g1^2*g2^2*M2*Yt*((57*trace[Ytc, Yt])/5 - (126*trace[Ytc, Yt]*Zeta[3])/5) + 
  g3^4*ht*((-160*trace[Ybc, Yb])/3 - (320*trace[Ytc, Yt])/3 - 
    16*trace[Ytc, Yt]*Zeta[3]) + g1^2*M1*MatMul[Yt, Ytc, Yt]*
   (-18*trace[Ytc, Yt] - (36*trace[Ytc, Yt]*Zeta[3])/5) + 
  g1^4*ht*((-91*trace[Ybc, Yb])/15 - (171*trace[Ytc, Yt])/4 - 
    (39*trace[Adj[Ye], Ye])/5 - (13*trace[Ytc, Yt]*Zeta[3])/10) + 
  g1^2*MatMul[Yt, Ytc, ht]*(11*trace[Ytc, Yt] - (6*trace[Ytc, Yt]*Zeta[3])/
     5) + g1^4*M1*Yt*((364*trace[Ybc, Yb])/15 + 171*trace[Ytc, Yt] + 
    (156*trace[Adj[Ye], Ye])/5 + (26*trace[Ytc, Yt]*Zeta[3])/5) + 
  g1^2*MatMul[ht, Ytc, Yt]*(16*trace[Ytc, Yt] + 12*trace[Ytc, Yt]*Zeta[3]) + 
  g1^2*g2^2*ht*((-57*trace[Ytc, Yt])/10 + (63*trace[Ytc, Yt]*Zeta[3])/5) + 
  g1^2*g3^2*ht*((-124*trace[Ytc, Yt])/3 + (208*trace[Ytc, Yt]*Zeta[3])/5) + 
  g3^4*M3*Yt*((640*trace[Ybc, Yb])/3 + (1280*trace[Ytc, Yt])/3 + 
    64*trace[Ytc, Yt]*Zeta[3]) + g2^2*M2*MatMul[Yt, Ytc, Yt]*
   (-90*trace[Ytc, Yt] + 108*trace[Ytc, Yt]*Zeta[3]) + 
  g2^2*g3^2*ht*(-132*trace[Ytc, Yt] + 144*trace[Ytc, Yt]*Zeta[3]) + 
  g3^2*MatMul[Yt, Ytc, ht]*(-32*trace[Ytc, Yt] + 
    192*trace[Ytc, Yt]*Zeta[3]) + g3^2*MatMul[ht, Ytc, Yt]*
   (-40*trace[Ytc, Yt] + 240*trace[Ytc, Yt]*Zeta[3]) + 
  g2^4*M2*Yt*(180*trace[Ybc, Yb] + 315*trace[Ytc, Yt] + 
    60*trace[Adj[Ye], Ye] + 378*trace[Ytc, Yt]*Zeta[3]) + 
  g1^2*M1*MatMul[Yt, Ybc, Yb]*((-32*trace[Ybc, Yb])/5 + 
    (16*trace[Adj[Ye], Ye])/5 + (48*trace[Ybc, Yb]*Zeta[3])/5 - 
    (24*trace[Adj[Ye], Ye]*Zeta[3])/5) + g1^2*MatMul[ht, Ybc, Yb]*
   ((16*trace[Ybc, Yb])/5 - (8*trace[Adj[Ye], Ye])/5 - 
    (24*trace[Ybc, Yb]*Zeta[3])/5 + (12*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*MatMul[Yt, Ybc, hb]*((32*trace[Ybc, Yb])/5 - 
    (16*trace[Adj[Ye], Ye])/5 - (48*trace[Ybc, Yb]*Zeta[3])/5 + 
    (24*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g3^2*Yt*(48*trace[hb, Ytc, Yt, Ybc] + 48*trace[Ytc, ht, Ybc, Yb] + 
    288*trace[Ytc, ht, Ytc, Yt] - 96*trace[hb, Ytc, Yt, Ybc]*Zeta[3] - 
    96*trace[Ytc, ht, Ybc, Yb]*Zeta[3] - 576*trace[Ytc, ht, Ytc, Yt]*
     Zeta[3]) + g1^2*Yt*((12*trace[hb, Ytc, Yt, Ybc])/5 + 
    (12*trace[Ytc, ht, Ybc, Yb])/5 + (228*trace[Ytc, ht, Ytc, Yt])/5 + 
    (12*trace[hb, Ytc, Yt, Ybc]*Zeta[3])/5 + 
    (12*trace[Ytc, ht, Ybc, Yb]*Zeta[3])/5 - 
    (72*trace[Ytc, ht, Ytc, Yt]*Zeta[3])/5) + 
  g2^2*Yt*(36*trace[hb, Ytc, Yt, Ybc] + 36*trace[Ytc, ht, Ybc, Yb] + 
    36*trace[Ytc, ht, Ytc, Yt] + 216*trace[Ytc, ht, Ytc, Yt]*Zeta[3]) + 
  g3^2*ht*(24*trace[Ytc, Yt, Ybc, Yb] + 72*trace[Ytc, Yt, Ytc, Yt] - 
    48*trace[Ytc, Yt, Ybc, Yb]*Zeta[3] - 144*trace[Ytc, Yt, Ytc, Yt]*
     Zeta[3]) + g2^2*M2*Yt*(-36*trace[Ytc, Yt, Ybc, Yb] - 
    18*trace[Ytc, Yt, Ytc, Yt] - 108*trace[Ytc, Yt, Ytc, Yt]*Zeta[3]) + 
  g1^2*ht*((6*trace[Ytc, Yt, Ybc, Yb])/5 + (57*trace[Ytc, Yt, Ytc, Yt])/5 + 
    (6*trace[Ytc, Yt, Ybc, Yb]*Zeta[3])/5 - 
    (18*trace[Ytc, Yt, Ytc, Yt]*Zeta[3])/5) + 
  g1^2*M1*Yt*((-12*trace[Ytc, Yt, Ybc, Yb])/5 - (114*trace[Ytc, Yt, Ytc, Yt])/
     5 - (12*trace[Ytc, Yt, Ybc, Yb]*Zeta[3])/5 + 
    (36*trace[Ytc, Yt, Ytc, Yt]*Zeta[3])/5) + 
  g2^2*ht*(18*trace[Ytc, Yt, Ybc, Yb] + 9*trace[Ytc, Yt, Ytc, Yt] + 
    54*trace[Ytc, Yt, Ytc, Yt]*Zeta[3]) + 
  g3^2*M3*Yt*(-48*trace[Ytc, Yt, Ybc, Yb] - 144*trace[Ytc, Yt, Ytc, Yt] + 
    96*trace[Ytc, Yt, Ybc, Yb]*Zeta[3] + 288*trace[Ytc, Yt, Ytc, Yt]*
     Zeta[3]) + Yt*(36*trace[Ybc, Yb]*trace[hb, Ytc, Yt, Ybc] + 
    12*trace[Adj[Ye], Ye]*trace[hb, Ytc, Yt, Ybc] + 
    36*trace[Ybc, Yb]*trace[Ytc, ht, Ybc, Yb] + 12*trace[Adj[Ye], Ye]*
     trace[Ytc, ht, Ybc, Yb] + 216*trace[Ytc, Yt]*trace[Ytc, ht, Ytc, Yt] + 
    36*trace[hb, Ybc]*trace[Ytc, Yt, Ybc, Yb] + 12*trace[he, Adj[Ye]]*
     trace[Ytc, Yt, Ybc, Yb] + 108*trace[ht, Ytc]*trace[Ytc, Yt, Ytc, Yt] + 
    18*trace[Ybc, hb, Ytc, Yt, Ybc, Yb] + 
    18*trace[Ybc, Yb, Ytc, ht, Ybc, Yb] + 
    18*trace[Ytc, Yt, Ybc, hb, Ybc, Yb] + 
    18*trace[Ytc, Yt, Ytc, ht, Ytc, Yt] + 
    108*trace[Ytc, Yt, Ytc, ht, Ytc, Yt]*Zeta[3]) + 
  ht*(18*trace[Ybc, Yb]*trace[Ytc, Yt, Ybc, Yb] + 
    6*trace[Adj[Ye], Ye]*trace[Ytc, Yt, Ybc, Yb] + 
    54*trace[Ytc, Yt]*trace[Ytc, Yt, Ytc, Yt] + 
    9*trace[Ybc, Yb, Ytc, Yt, Ybc, Yb] + 3*trace[Ytc, Yt, Ytc, Yt, Ytc, Yt] + 
    18*trace[Ytc, Yt, Ytc, Yt, Ytc, Yt]*Zeta[3])}

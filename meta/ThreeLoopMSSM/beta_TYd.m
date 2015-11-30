{(-7*g1^2*hb)/15 - 3*g2^2*hb - (16*g3^2*hb)/3 + (14*g1^2*M1*Yb)/15 + 
  6*g2^2*M2*Yb + (32*g3^2*M3*Yb)/3 + 5*MatMul[hb, Ybc, Yb] + 
  MatMul[hb, Ytc, Yt] + 4*MatMul[Yb, Ybc, hb] + 2*MatMul[Yb, Ytc, ht] + 
  Yb*(6*trace[hb, Ybc] + 2*trace[he, Adj[Ye]]) + 
  hb*(3*trace[Ybc, Yb] + trace[Adj[Ye], Ye]), 
 (287*g1^4*hb)/90 + g1^2*g2^2*hb + (15*g2^4*hb)/2 + (8*g1^2*g3^2*hb)/9 + 
  8*g2^2*g3^2*hb - (16*g3^4*hb)/9 - (574*g1^4*M1*Yb)/45 - 2*g1^2*g2^2*M1*Yb - 
  (16*g1^2*g3^2*M1*Yb)/9 - 2*g1^2*g2^2*M2*Yb - 30*g2^4*M2*Yb - 
  16*g2^2*g3^2*M2*Yb - (16*g1^2*g3^2*M3*Yb)/9 - 16*g2^2*g3^2*M3*Yb + 
  (64*g3^4*M3*Yb)/9 + (6*g1^2*MatMul[hb, Ybc, Yb])/5 + 
  12*g2^2*MatMul[hb, Ybc, Yb] + (4*g1^2*MatMul[hb, Ytc, Yt])/5 + 
  (6*g1^2*MatMul[Yb, Ybc, hb])/5 + 6*g2^2*MatMul[Yb, Ybc, hb] - 
  (8*g1^2*M1*MatMul[Yb, Ybc, Yb])/5 - 12*g2^2*M2*MatMul[Yb, Ybc, Yb] + 
  (8*g1^2*MatMul[Yb, Ytc, ht])/5 - (8*g1^2*M1*MatMul[Yb, Ytc, Yt])/5 - 
  6*MatMul[hb, Ybc, Yb, Ybc, Yb] - 4*MatMul[hb, Ytc, Yt, Ybc, Yb] - 
  2*MatMul[hb, Ytc, Yt, Ytc, Yt] - 8*MatMul[Yb, Ybc, hb, Ybc, Yb] - 
  6*MatMul[Yb, Ybc, Yb, Ybc, hb] - 4*MatMul[Yb, Ytc, ht, Ybc, Yb] - 
  4*MatMul[Yb, Ytc, ht, Ytc, Yt] - 2*MatMul[Yb, Ytc, Yt, Ybc, hb] - 
  4*MatMul[Yb, Ytc, Yt, Ytc, ht] + 32*g3^2*Yb*trace[hb, Ybc] + 
  MatMul[Yb, Ybc, Yb]*(-18*trace[hb, Ybc] - 6*trace[he, Adj[Ye]]) + 
  g1^2*Yb*((-4*trace[hb, Ybc])/5 + (12*trace[he, Adj[Ye]])/5) - 
  6*MatMul[Yb, Ytc, Yt]*trace[ht, Ytc] + 16*g3^2*hb*trace[Ybc, Yb] - 
  32*g3^2*M3*Yb*trace[Ybc, Yb] - 3*MatMul[hb, Ytc, Yt]*trace[Ytc, Yt] - 
  6*MatMul[Yb, Ytc, ht]*trace[Ytc, Yt] + MatMul[hb, Ybc, Yb]*
   (-15*trace[Ybc, Yb] - 5*trace[Adj[Ye], Ye]) + 
  MatMul[Yb, Ybc, hb]*(-12*trace[Ybc, Yb] - 4*trace[Adj[Ye], Ye]) + 
  g1^2*M1*Yb*((4*trace[Ybc, Yb])/5 - (12*trace[Adj[Ye], Ye])/5) + 
  g1^2*hb*((-2*trace[Ybc, Yb])/5 + (6*trace[Adj[Ye], Ye])/5) + 
  Yb*(-6*trace[hb, Ytc, Yt, Ybc] - 36*trace[Ybc, hb, Ybc, Yb] - 
    6*trace[Ytc, ht, Ybc, Yb] - 12*trace[Adj[Ye], he, Adj[Ye], Ye]) + 
  hb*(-9*trace[Ybc, Yb, Ybc, Yb] - 3*trace[Ytc, Yt, Ybc, Yb] - 
    3*trace[Adj[Ye], Ye, Adj[Ye], Ye]), (-8*g1^2*g2^2*g3^2*hb)/5 + 
  (16*g1^2*g2^2*g3^2*M1*Yb)/5 + (16*g1^2*g2^2*g3^2*M2*Yb)/5 + 
  (16*g1^2*g2^2*g3^2*M3*Yb)/5 - 4*g2^2*g3^2*MatMul[hb, Ytc, Yt] - 
  8*g2^2*g3^2*MatMul[Yb, Ytc, ht] + 8*g2^2*g3^2*M2*MatMul[Yb, Ytc, Yt] + 
  8*g2^2*g3^2*M3*MatMul[Yb, Ytc, Yt] + 64*g3^2*MatMul[hb, Ybc, Yb, Ybc, Yb] + 
  (128*g3^2*MatMul[hb, Ytc, Yt, Ybc, Yb])/3 + 
  (64*g3^2*MatMul[hb, Ytc, Yt, Ytc, Yt])/3 + 
  (4*g1^2*MatMul[Yb, Ybc, hb, Ybc, Yb])/15 + 
  12*g2^2*MatMul[Yb, Ybc, hb, Ybc, Yb] + 
  (256*g3^2*MatMul[Yb, Ybc, hb, Ybc, Yb])/3 + 
  64*g3^2*MatMul[Yb, Ybc, Yb, Ybc, hb] - 
  (4*g1^2*M1*MatMul[Yb, Ybc, Yb, Ybc, Yb])/15 - 
  12*g2^2*M2*MatMul[Yb, Ybc, Yb, Ybc, Yb] - 
  (256*g3^2*M3*MatMul[Yb, Ybc, Yb, Ybc, Yb])/3 + 
  (128*g3^2*MatMul[Yb, Ytc, ht, Ybc, Yb])/3 + 
  (128*g3^2*MatMul[Yb, Ytc, ht, Ytc, Yt])/3 + 
  (64*g3^2*MatMul[Yb, Ytc, Yt, Ybc, hb])/3 - 
  (128*g3^2*M3*MatMul[Yb, Ytc, Yt, Ybc, Yb])/3 + 
  (128*g3^2*MatMul[Yb, Ytc, Yt, Ytc, ht])/3 - 
  (128*g3^2*M3*MatMul[Yb, Ytc, Yt, Ytc, Yt])/3 - 
  4*MatMul[hb, Ytc, Yt, Ybc, Yb, Ybc, Yb] + 
  4*MatMul[hb, Ytc, Yt, Ybc, Yb, Ytc, Yt] + 
  12*MatMul[hb, Ytc, Yt, Ytc, Yt, Ybc, Yb] + 
  4*MatMul[Yb, Ybc, hb, Ytc, Yt, Ybc, Yb] + 
  4*MatMul[Yb, Ybc, Yb, Ytc, ht, Ybc, Yb] + 
  6*MatMul[Yb, Ybc, Yb, Ytc, Yt, Ybc, hb] - 
  4*MatMul[Yb, Ytc, ht, Ybc, Yb, Ybc, Yb] + 
  8*MatMul[Yb, Ytc, ht, Ybc, Yb, Ytc, Yt] + 
  12*MatMul[Yb, Ytc, ht, Ytc, Yt, Ybc, Yb] - 
  4*MatMul[Yb, Ytc, Yt, Ybc, hb, Ybc, Yb] + 
  8*MatMul[Yb, Ytc, Yt, Ybc, hb, Ytc, Yt] - 
  2*MatMul[Yb, Ytc, Yt, Ybc, Yb, Ybc, hb] + 
  8*MatMul[Yb, Ytc, Yt, Ybc, Yb, Ytc, ht] + 
  12*MatMul[Yb, Ytc, Yt, Ytc, ht, Ybc, Yb] + 
  6*MatMul[Yb, Ytc, Yt, Ytc, Yt, Ybc, hb] + MatMul[Yb, Ybc, Yb, Ybc, Yb]*
   (24*trace[hb, Ybc] + 8*trace[he, Adj[Ye]]) + 
  36*g2^2*MatMul[Yb, Ytc, Yt]*trace[ht, Ytc] + 
  12*MatMul[Yb, Ytc, Yt, Ytc, Yt]*trace[ht, Ytc] + 
  MatMul[Yb, Ytc, Yt, Ybc, Yb]*(-12*trace[hb, Ybc] - 4*trace[he, Adj[Ye]] + 
    24*trace[ht, Ytc]) + 18*g2^2*MatMul[hb, Ytc, Yt]*trace[Ytc, Yt] + 
  36*g2^2*MatMul[Yb, Ytc, ht]*trace[Ytc, Yt] - 36*g2^2*M2*MatMul[Yb, Ytc, Yt]*
   trace[Ytc, Yt] + 6*MatMul[hb, Ytc, Yt, Ytc, Yt]*trace[Ytc, Yt] + 
  12*MatMul[Yb, Ytc, ht, Ytc, Yt]*trace[Ytc, Yt] + 
  12*MatMul[Yb, Ytc, Yt, Ytc, ht]*trace[Ytc, Yt] + 
  MatMul[hb, Ytc, Yt, Ybc, Yb]*(-12*trace[Ybc, Yb] + 24*trace[Ytc, Yt] - 
    4*trace[Adj[Ye], Ye]) + MatMul[Yb, Ytc, ht, Ybc, Yb]*
   (-12*trace[Ybc, Yb] + 24*trace[Ytc, Yt] - 4*trace[Adj[Ye], Ye]) + 
  MatMul[Yb, Ytc, Yt, Ybc, hb]*(-6*trace[Ybc, Yb] + 12*trace[Ytc, Yt] - 
    2*trace[Adj[Ye], Ye]) + MatMul[hb, Ybc, Yb, Ybc, Yb]*
   (18*trace[Ybc, Yb] + 6*trace[Adj[Ye], Ye]) + 
  MatMul[Yb, Ybc, Yb, Ybc, hb]*(18*trace[Ybc, Yb] + 6*trace[Adj[Ye], Ye]) + 
  MatMul[Yb, Ybc, hb, Ybc, Yb]*(24*trace[Ybc, Yb] + 8*trace[Adj[Ye], Ye]) + 
  MatMul[Yb, Ytc, Yt]*(-36*trace[ht, Ytc]*trace[Ytc, Yt] + 
    12*trace[hb, Ytc, Yt, Ybc] + 12*trace[Ytc, ht, Ybc, Yb] + 
    72*trace[Ytc, ht, Ytc, Yt]) + MatMul[hb, Ytc, Yt]*
   (-9*trace[Ytc, Yt]^2 + 6*trace[Ytc, Yt, Ybc, Yb] + 
    18*trace[Ytc, Yt, Ytc, Yt]) + MatMul[Yb, Ytc, ht]*
   (-18*trace[Ytc, Yt]^2 + 12*trace[Ytc, Yt, Ybc, Yb] + 
    36*trace[Ytc, Yt, Ytc, Yt]) + MatMul[Yb, Ybc, Yb]*
   (-108*trace[hb, Ybc]*trace[Ybc, Yb] - 36*trace[he, Adj[Ye]]*
     trace[Ybc, Yb] - 36*trace[hb, Ybc]*trace[Adj[Ye], Ye] - 
    12*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] + 36*trace[hb, Ytc, Yt, Ybc] + 
    216*trace[Ybc, hb, Ybc, Yb] + 36*trace[Ytc, ht, Ybc, Yb] + 
    72*trace[Adj[Ye], he, Adj[Ye], Ye]) + MatMul[Yb, Ybc, hb]*
   (-36*trace[Ybc, Yb]^2 - 24*trace[Ybc, Yb]*trace[Adj[Ye], Ye] - 
    4*trace[Adj[Ye], Ye]^2 + 72*trace[Ybc, Yb, Ybc, Yb] + 
    24*trace[Ytc, Yt, Ybc, Yb] + 24*trace[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  MatMul[hb, Ybc, Yb]*(-45*trace[Ybc, Yb]^2 - 30*trace[Ybc, Yb]*
     trace[Adj[Ye], Ye] - 5*trace[Adj[Ye], Ye]^2 + 
    90*trace[Ybc, Yb, Ybc, Yb] + 30*trace[Ytc, Yt, Ybc, Yb] + 
    30*trace[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  g3^6*M3*Yb*(-10880/9 - 3840*Zeta[3]) + g2^6*M2*Yb*(-1035 - 1890*Zeta[3]) + 
  g3^4*MatMul[hb, Ybc, Yb]*(40/3 - (1360*Zeta[3])/3) + 
  g3^4*MatMul[Yb, Ybc, hb]*(32/3 - (1088*Zeta[3])/3) + 
  g2^4*g3^2*hb*(140 - 216*Zeta[3]) + g2^2*g3^2*M2*MatMul[Yb, Ybc, Yb]*
   (184 - 192*Zeta[3]) + g2^2*g3^2*M3*MatMul[Yb, Ybc, Yb]*
   (184 - 192*Zeta[3]) + g3^4*MatMul[Yb, Ytc, ht]*(16/3 - (544*Zeta[3])/3) + 
  g2^2*g3^4*hb*(68 - 144*Zeta[3]) + g3^4*MatMul[hb, Ytc, Yt]*
   (8/3 - (272*Zeta[3])/3) + g2^4*MatMul[Yb, Ybc, hb]*(-66 - 72*Zeta[3]) + 
  g2^4*MatMul[Yb, Ytc, ht]*(-45/2 - 63*Zeta[3]) + 
  g2^4*MatMul[hb, Ybc, Yb]*(-393/4 - (99*Zeta[3])/2) + 
  g1^2*g3^2*M1*MatMul[Yb, Ybc, Yb]*(296/15 - (192*Zeta[3])/5) + 
  g1^2*g3^2*M3*MatMul[Yb, Ybc, Yb]*(296/15 - (192*Zeta[3])/5) + 
  g2^2*M2*MatMul[Yb, Ytc, Yt, Ytc, Yt]*(6 - 36*Zeta[3]) + 
  g2^2*MatMul[hb, Ytc, Yt, Ybc, Yb]*(18 - 36*Zeta[3]) + 
  g2^2*MatMul[Yb, Ytc, ht, Ybc, Yb]*(18 - 36*Zeta[3]) + 
  g1^2*g3^4*hb*(212/9 - (176*Zeta[3])/5) + g2^4*MatMul[hb, Ytc, Yt]*
   (-45/4 - (63*Zeta[3])/2) + g1^4*g3^2*hb*(3892/225 - (616*Zeta[3])/25) + 
  g1^2*g2^2*M1*MatMul[Yb, Ybc, Yb]*(127/5 - (102*Zeta[3])/5) + 
  g1^2*g2^2*M2*MatMul[Yb, Ybc, Yb]*(127/5 - (102*Zeta[3])/5) + 
  g2^2*MatMul[Yb, Ytc, Yt, Ybc, hb]*(9 - 18*Zeta[3]) + 
  g1^2*g2^2*M1*MatMul[Yb, Ytc, Yt]*(59/5 - 18*Zeta[3]) + 
  g1^2*g2^2*M2*MatMul[Yb, Ytc, Yt]*(59/5 - 18*Zeta[3]) + 
  g2^2*MatMul[hb, Ybc, Yb, Ybc, Yb]*(15 - 18*Zeta[3]) + 
  g1^2*g3^2*M1*MatMul[Yb, Ytc, Yt]*(136/5 - (256*Zeta[3])/15) + 
  g1^2*g3^2*M3*MatMul[Yb, Ytc, Yt]*(136/5 - (256*Zeta[3])/15) + 
  g1^2*g2^4*hb*(17/2 - (81*Zeta[3])/5) + g1^2*MatMul[Yb, Ytc, ht, Ytc, Yt]*
   (22/3 - 12*Zeta[3]) + g1^2*MatMul[Yb, Ytc, Yt, Ytc, ht]*
   (22/3 - 12*Zeta[3]) + g1^6*hb*(194651/6750 - (1393*Zeta[3])/125) + 
  g1^4*g2^2*hb*(109/50 - (189*Zeta[3])/25) + 
  g1^2*M1*MatMul[Yb, Ytc, Yt, Ybc, Yb]*(58/15 - (36*Zeta[3])/5) + 
  g1^2*MatMul[hb, Ytc, Yt, Ytc, Yt]*(11/3 - 6*Zeta[3]) + 
  g1^4*M1*MatMul[Yb, Ytc, Yt]*(3767/75 - (286*Zeta[3])/75) + 
  g1^2*MatMul[Yb, Ybc, Yb, Ybc, hb]*(3/5 - (6*Zeta[3])/5) + 
  g1^4*M1*MatMul[Yb, Ybc, Yb]*(5269/75 - (14*Zeta[3])/25) + 
  g1^4*MatMul[hb, Ybc, Yb]*(-8639/300 + (7*Zeta[3])/150) + 
  g1^4*MatMul[Yb, Ybc, hb]*(-1792/75 + (28*Zeta[3])/75) + 
  g1^4*MatMul[hb, Ytc, Yt]*(-3767/300 + (143*Zeta[3])/150) + 
  6*MatMul[hb, Ytc, Yt, Ytc, Yt, Ytc, Yt]*Zeta[3] + 
  12*MatMul[Yb, Ytc, ht, Ytc, Yt, Ytc, Yt]*Zeta[3] + 
  12*MatMul[Yb, Ytc, Yt, Ytc, ht, Ytc, Yt]*Zeta[3] + 
  12*MatMul[Yb, Ytc, Yt, Ytc, Yt, Ytc, ht]*Zeta[3] + 
  g1^2*MatMul[hb, Ybc, Yb, Ybc, Yb]*(-1/5 + (6*Zeta[3])/5) + 
  g1^4*MatMul[Yb, Ytc, ht]*(-3767/150 + (143*Zeta[3])/75) + 
  g1^2*MatMul[Yb, Ytc, Yt, Ybc, hb]*(-29/15 + (18*Zeta[3])/5) + 
  g1^2*MatMul[hb, Ytc, Yt, Ybc, Yb]*(-58/15 + (36*Zeta[3])/5) + 
  g1^2*MatMul[Yb, Ytc, ht, Ybc, Yb]*(-58/15 + (36*Zeta[3])/5) + 
  g1^2*g3^2*MatMul[hb, Ytc, Yt]*(-68/5 + (128*Zeta[3])/15) + 
  g1^2*g2^2*MatMul[hb, Ytc, Yt]*(-59/10 + 9*Zeta[3]) + 
  g1^2*g2^2*MatMul[Yb, Ybc, hb]*(-84/5 + 12*Zeta[3]) + 
  g1^2*M1*MatMul[Yb, Ytc, Yt, Ytc, Yt]*(-22/3 + 12*Zeta[3]) + 
  g1^4*g2^2*M2*Yb*(-109/25 + (378*Zeta[3])/25) + 
  g1^2*g3^2*MatMul[Yb, Ytc, ht]*(-136/5 + (256*Zeta[3])/15) + 
  g1^2*g2^2*MatMul[Yb, Ytc, ht]*(-59/5 + 18*Zeta[3]) + 
  g2^2*MatMul[hb, Ytc, Yt, Ytc, Yt]*(-3 + 18*Zeta[3]) + 
  g2^2*MatMul[Yb, Ybc, Yb, Ybc, hb]*(3 + 18*Zeta[3]) + 
  g1^2*g2^2*MatMul[hb, Ybc, Yb]*(-213/10 + (93*Zeta[3])/5) + 
  MatMul[Yb, Ybc, Yb, Ybc, Yb, Ybc, hb]*(6 + 24*Zeta[3]) + 
  g1^2*g3^2*MatMul[Yb, Ybc, hb]*(-224/15 + (416*Zeta[3])/15) + 
  g1^2*g3^2*MatMul[hb, Ybc, Yb]*(-44/3 + (448*Zeta[3])/15) + 
  MatMul[hb, Ybc, Yb, Ybc, Yb, Ybc, Yb]*(12 + 30*Zeta[3]) + 
  g1^4*g2^2*M1*Yb*(-218/25 + (756*Zeta[3])/25) + 
  g1^2*g2^4*M1*Yb*(-17 + (162*Zeta[3])/5) + 
  g2^2*M2*MatMul[Yb, Ytc, Yt, Ybc, Yb]*(-18 + 36*Zeta[3]) + 
  g2^2*MatMul[Yb, Ytc, ht, Ytc, Yt]*(-6 + 36*Zeta[3]) + 
  g2^2*MatMul[Yb, Ytc, Yt, Ytc, ht]*(-6 + 36*Zeta[3]) + 
  MatMul[Yb, Ybc, hb, Ybc, Yb, Ybc, Yb]*(12 + 36*Zeta[3]) + 
  MatMul[Yb, Ybc, Yb, Ybc, hb, Ybc, Yb]*(12 + 36*Zeta[3]) + 
  g1^4*g3^2*M3*Yb*(-7784/225 + (1232*Zeta[3])/25) + 
  g1^2*g2^4*M2*Yb*(-34 + (324*Zeta[3])/5) + 
  g1^6*M1*Yb*(-194651/1125 + (8358*Zeta[3])/125) + 
  g1^2*g3^4*M1*Yb*(-424/9 + (352*Zeta[3])/5) + g2^2*g3^2*MatMul[Yb, Ybc, hb]*
   (-96 + 96*Zeta[3]) + g1^4*g3^2*M1*Yb*(-15568/225 + (2464*Zeta[3])/25) + 
  g2^4*M2*MatMul[Yb, Ytc, Yt]*(45 + 126*Zeta[3]) + 
  g1^2*g3^4*M3*Yb*(-848/9 + (704*Zeta[3])/5) + g2^4*M2*MatMul[Yb, Ybc, Yb]*
   (219 + 162*Zeta[3]) + g2^2*g3^2*MatMul[hb, Ybc, Yb]*(-180 + 192*Zeta[3]) + 
  g2^2*g3^4*M2*Yb*(-136 + 288*Zeta[3]) + g2^6*hb*(345/2 + 315*Zeta[3]) + 
  g3^4*M3*MatMul[Yb, Ytc, Yt]*(-32/3 + (1088*Zeta[3])/3) + 
  g2^4*g3^2*M3*Yb*(-280 + 432*Zeta[3]) + g2^2*g3^4*M3*Yb*
   (-272 + 576*Zeta[3]) + g3^6*hb*(5440/27 + 640*Zeta[3]) + 
  g2^4*g3^2*M2*Yb*(-560 + 864*Zeta[3]) + g3^4*M3*MatMul[Yb, Ybc, Yb]*
   (-32 + 1088*Zeta[3]) + g3^4*Yb*((-640*trace[hb, Ybc])/3 - 
    (320*trace[ht, Ytc])/3 - 32*trace[hb, Ybc]*Zeta[3]) + 
  g1^2*g3^2*Yb*((-568*trace[hb, Ybc])/15 + (224*trace[hb, Ybc]*Zeta[3])/5) + 
  g2^2*g3^2*Yb*(-264*trace[hb, Ybc] + 288*trace[hb, Ybc]*Zeta[3]) + 
  g3^2*MatMul[Yb, Ybc, Yb]*(-48*trace[hb, Ybc] + 48*trace[he, Adj[Ye]] + 
    288*trace[hb, Ybc]*Zeta[3]) + g2^4*Yb*((-315*trace[hb, Ybc])/2 - 
    (105*trace[he, Adj[Ye]])/2 - 90*trace[ht, Ytc] - 
    189*trace[hb, Ybc]*Zeta[3] - 63*trace[he, Adj[Ye]]*Zeta[3]) + 
  g2^2*MatMul[Yb, Ybc, Yb]*(90*trace[hb, Ybc] + 30*trace[he, Adj[Ye]] - 
    108*trace[hb, Ybc]*Zeta[3] - 36*trace[he, Adj[Ye]]*Zeta[3]) + 
  g1^4*Yb*((-63*trace[hb, Ybc])/2 - (633*trace[he, Adj[Ye]])/10 - 
    (182*trace[ht, Ytc])/15 - (77*trace[hb, Ybc]*Zeta[3])/25 + 
    (81*trace[he, Adj[Ye]]*Zeta[3])/25) + g1^2*MatMul[Yb, Ybc, Yb]*
   ((102*trace[hb, Ybc])/5 - (46*trace[he, Adj[Ye]])/5 - 
    (108*trace[hb, Ybc]*Zeta[3])/5 + (84*trace[he, Adj[Ye]]*Zeta[3])/5) + 
  g1^2*g2^2*Yb*((-3*trace[hb, Ybc])/5 - (81*trace[he, Adj[Ye]])/5 - 
    18*trace[hb, Ybc]*Zeta[3] + (162*trace[he, Adj[Ye]]*Zeta[3])/5) + 
  g1^2*MatMul[Yb, Ytc, Yt]*(4*trace[ht, Ytc] - (48*trace[ht, Ytc]*Zeta[3])/
     5) + g3^2*MatMul[Yb, Ytc, Yt]*(-16*trace[ht, Ytc] + 
    96*trace[ht, Ytc]*Zeta[3]) + g2^2*g3^2*M2*Yb*
   (264*trace[Ybc, Yb] - 288*trace[Ybc, Yb]*Zeta[3]) + 
  g2^2*g3^2*M3*Yb*(264*trace[Ybc, Yb] - 288*trace[Ybc, Yb]*Zeta[3]) + 
  g3^2*M3*MatMul[Yb, Ybc, Yb]*(48*trace[Ybc, Yb] - 48*trace[Adj[Ye], Ye] - 
    288*trace[Ybc, Yb]*Zeta[3]) + g1^2*g3^2*M1*Yb*
   ((568*trace[Ybc, Yb])/15 - (224*trace[Ybc, Yb]*Zeta[3])/5) + 
  g1^2*g3^2*M3*Yb*((568*trace[Ybc, Yb])/15 - (224*trace[Ybc, Yb]*Zeta[3])/
     5) + g3^4*hb*((-320*trace[Ybc, Yb])/3 - (160*trace[Ytc, Yt])/3 - 
    16*trace[Ybc, Yb]*Zeta[3]) + g1^2*g3^2*hb*((-284*trace[Ybc, Yb])/15 + 
    (112*trace[Ybc, Yb]*Zeta[3])/5) + 
  g3^4*M3*Yb*((1280*trace[Ybc, Yb])/3 + (640*trace[Ytc, Yt])/3 + 
    64*trace[Ybc, Yb]*Zeta[3]) + g2^2*g3^2*hb*(-132*trace[Ybc, Yb] + 
    144*trace[Ybc, Yb]*Zeta[3]) + g3^2*MatMul[Yb, Ybc, hb]*
   (-32*trace[Ybc, Yb] + 32*trace[Adj[Ye], Ye] + 
    192*trace[Ybc, Yb]*Zeta[3]) + g3^2*MatMul[hb, Ybc, Yb]*
   (-40*trace[Ybc, Yb] + 40*trace[Adj[Ye], Ye] + 
    240*trace[Ybc, Yb]*Zeta[3]) + g3^2*M3*MatMul[Yb, Ytc, Yt]*
   (16*trace[Ytc, Yt] - 96*trace[Ytc, Yt]*Zeta[3]) + 
  g1^2*MatMul[Yb, Ytc, ht]*(4*trace[Ytc, Yt] - (48*trace[Ytc, Yt]*Zeta[3])/
     5) + g1^2*MatMul[hb, Ytc, Yt]*(2*trace[Ytc, Yt] - 
    (24*trace[Ytc, Yt]*Zeta[3])/5) + g1^2*M1*MatMul[Yb, Ytc, Yt]*
   (-4*trace[Ytc, Yt] + (48*trace[Ytc, Yt]*Zeta[3])/5) + 
  g3^2*MatMul[hb, Ytc, Yt]*(-8*trace[Ytc, Yt] + 48*trace[Ytc, Yt]*Zeta[3]) + 
  g3^2*MatMul[Yb, Ytc, ht]*(-16*trace[Ytc, Yt] + 96*trace[Ytc, Yt]*Zeta[3]) + 
  g2^2*MatMul[hb, Ybc, Yb]*(72*trace[Ybc, Yb] + 24*trace[Adj[Ye], Ye] - 
    108*trace[Ybc, Yb]*Zeta[3] - 36*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g1^2*g2^2*M1*Yb*((3*trace[Ybc, Yb])/5 + (81*trace[Adj[Ye], Ye])/5 + 
    18*trace[Ybc, Yb]*Zeta[3] - (162*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*g2^2*M2*Yb*((3*trace[Ybc, Yb])/5 + (81*trace[Adj[Ye], Ye])/5 + 
    18*trace[Ybc, Yb]*Zeta[3] - (162*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g2^4*hb*((-315*trace[Ybc, Yb])/4 - 45*trace[Ytc, Yt] - 
    (105*trace[Adj[Ye], Ye])/4 - (189*trace[Ybc, Yb]*Zeta[3])/2 - 
    (63*trace[Adj[Ye], Ye]*Zeta[3])/2) + g2^2*MatMul[Yb, Ybc, hb]*
   (63*trace[Ybc, Yb] + 21*trace[Adj[Ye], Ye] - 54*trace[Ybc, Yb]*Zeta[3] - 
    18*trace[Adj[Ye], Ye]*Zeta[3]) + g1^2*M1*MatMul[Yb, Ybc, Yb]*
   ((-102*trace[Ybc, Yb])/5 + (46*trace[Adj[Ye], Ye])/5 + 
    (108*trace[Ybc, Yb]*Zeta[3])/5 - (84*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^4*M1*Yb*(63*trace[Ybc, Yb] + (364*trace[Ytc, Yt])/15 + 
    (633*trace[Adj[Ye], Ye])/5 + (154*trace[Ybc, Yb]*Zeta[3])/25 - 
    (162*trace[Adj[Ye], Ye]*Zeta[3])/25) + 
  g1^4*hb*((-63*trace[Ybc, Yb])/4 - (91*trace[Ytc, Yt])/15 - 
    (633*trace[Adj[Ye], Ye])/20 - (77*trace[Ybc, Yb]*Zeta[3])/50 + 
    (81*trace[Adj[Ye], Ye]*Zeta[3])/50) + g1^2*MatMul[Yb, Ybc, hb]*
   ((67*trace[Ybc, Yb])/5 - (31*trace[Adj[Ye], Ye])/5 - 
    (78*trace[Ybc, Yb]*Zeta[3])/5 + (54*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*MatMul[hb, Ybc, Yb]*((86*trace[Ybc, Yb])/5 - 
    (38*trace[Adj[Ye], Ye])/5 - (84*trace[Ybc, Yb]*Zeta[3])/5 + 
    (72*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*g2^2*hb*((-3*trace[Ybc, Yb])/10 - (81*trace[Adj[Ye], Ye])/10 - 
    9*trace[Ybc, Yb]*Zeta[3] + (81*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g2^2*M2*MatMul[Yb, Ybc, Yb]*(-90*trace[Ybc, Yb] - 30*trace[Adj[Ye], Ye] + 
    108*trace[Ybc, Yb]*Zeta[3] + 36*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g2^4*M2*Yb*(315*trace[Ybc, Yb] + 180*trace[Ytc, Yt] + 
    105*trace[Adj[Ye], Ye] + 378*trace[Ybc, Yb]*Zeta[3] + 
    126*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g3^2*Yb*(48*trace[hb, Ytc, Yt, Ybc] + 288*trace[Ybc, hb, Ybc, Yb] + 
    48*trace[Ytc, ht, Ybc, Yb] - 96*trace[hb, Ytc, Yt, Ybc]*Zeta[3] - 
    576*trace[Ybc, hb, Ybc, Yb]*Zeta[3] - 96*trace[Ytc, ht, Ybc, Yb]*
     Zeta[3]) + g3^2*hb*(72*trace[Ybc, Yb, Ybc, Yb] + 
    24*trace[Ytc, Yt, Ybc, Yb] - 144*trace[Ybc, Yb, Ybc, Yb]*Zeta[3] - 
    48*trace[Ytc, Yt, Ybc, Yb]*Zeta[3]) + 
  g3^2*M3*Yb*(-144*trace[Ybc, Yb, Ybc, Yb] - 48*trace[Ytc, Yt, Ybc, Yb] + 
    288*trace[Ybc, Yb, Ybc, Yb]*Zeta[3] + 96*trace[Ytc, Yt, Ybc, Yb]*
     Zeta[3]) + g1^2*Yb*((-24*trace[hb, Ytc, Yt, Ybc])/5 + 
    12*trace[Ybc, hb, Ybc, Yb] - (24*trace[Ytc, ht, Ybc, Yb])/5 + 
    36*trace[Adj[Ye], he, Adj[Ye], Ye] + (84*trace[hb, Ytc, Yt, Ybc]*Zeta[3])/
     5 + (216*trace[Ybc, hb, Ybc, Yb]*Zeta[3])/5 + 
    (84*trace[Ytc, ht, Ybc, Yb]*Zeta[3])/5 - 
    (216*trace[Adj[Ye], he, Adj[Ye], Ye]*Zeta[3])/5) + 
  g2^2*Yb*(36*trace[hb, Ytc, Yt, Ybc] + 36*trace[Ybc, hb, Ybc, Yb] + 
    36*trace[Ytc, ht, Ybc, Yb] + 12*trace[Adj[Ye], he, Adj[Ye], Ye] + 
    216*trace[Ybc, hb, Ybc, Yb]*Zeta[3] + 72*trace[Adj[Ye], he, Adj[Ye], Ye]*
     Zeta[3]) + g2^2*M2*Yb*(-18*trace[Ybc, Yb, Ybc, Yb] - 
    36*trace[Ytc, Yt, Ybc, Yb] - 6*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
    108*trace[Ybc, Yb, Ybc, Yb]*Zeta[3] - 36*trace[Adj[Ye], Ye, Adj[Ye], Ye]*
     Zeta[3]) + g1^2*hb*(3*trace[Ybc, Yb, Ybc, Yb] - 
    (12*trace[Ytc, Yt, Ybc, Yb])/5 + 9*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    (54*trace[Ybc, Yb, Ybc, Yb]*Zeta[3])/5 + 
    (42*trace[Ytc, Yt, Ybc, Yb]*Zeta[3])/5 - 
    (54*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3])/5) + 
  g2^2*hb*(9*trace[Ybc, Yb, Ybc, Yb] + 18*trace[Ytc, Yt, Ybc, Yb] + 
    3*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 54*trace[Ybc, Yb, Ybc, Yb]*Zeta[3] + 
    18*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3]) + 
  g1^2*M1*Yb*(-6*trace[Ybc, Yb, Ybc, Yb] + (24*trace[Ytc, Yt, Ybc, Yb])/5 - 
    18*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
    (108*trace[Ybc, Yb, Ybc, Yb]*Zeta[3])/5 - 
    (84*trace[Ytc, Yt, Ybc, Yb]*Zeta[3])/5 + 
    (108*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3])/5) + 
  Yb*(36*trace[Ytc, Yt]*trace[hb, Ytc, Yt, Ybc] + 
    216*trace[Ybc, Yb]*trace[Ybc, hb, Ybc, Yb] + 72*trace[Adj[Ye], Ye]*
     trace[Ybc, hb, Ybc, Yb] + 108*trace[hb, Ybc]*trace[Ybc, Yb, Ybc, Yb] + 
    36*trace[he, Adj[Ye]]*trace[Ybc, Yb, Ybc, Yb] + 
    36*trace[Ytc, Yt]*trace[Ytc, ht, Ybc, Yb] + 36*trace[ht, Ytc]*
     trace[Ytc, Yt, Ybc, Yb] + 72*trace[Ybc, Yb]*trace[Adj[Ye], he, Adj[Ye], 
      Ye] + 24*trace[Adj[Ye], Ye]*trace[Adj[Ye], he, Adj[Ye], Ye] + 
    36*trace[hb, Ybc]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    12*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    18*trace[hb, Ytc, Yt, Ytc, Yt, Ybc] + 
    18*trace[Ybc, Yb, Ybc, hb, Ybc, Yb] + 
    18*trace[Ytc, ht, Ytc, Yt, Ybc, Yb] + 
    18*trace[Ytc, Yt, Ytc, ht, Ybc, Yb] + 
    6*trace[Adj[Ye], Ye, Adj[Ye], he, Adj[Ye], Ye] + 
    108*trace[Ybc, Yb, Ybc, hb, Ybc, Yb]*Zeta[3] + 
    36*trace[Adj[Ye], Ye, Adj[Ye], he, Adj[Ye], Ye]*Zeta[3]) + 
  hb*(54*trace[Ybc, Yb]*trace[Ybc, Yb, Ybc, Yb] + 
    18*trace[Adj[Ye], Ye]*trace[Ybc, Yb, Ybc, Yb] + 
    18*trace[Ytc, Yt]*trace[Ytc, Yt, Ybc, Yb] + 18*trace[Ybc, Yb]*
     trace[Adj[Ye], Ye, Adj[Ye], Ye] + 6*trace[Adj[Ye], Ye]*
     trace[Adj[Ye], Ye, Adj[Ye], Ye] + 3*trace[Ybc, Yb, Ybc, Yb, Ybc, Yb] + 
    9*trace[Ytc, Yt, Ytc, Yt, Ybc, Yb] + trace[Adj[Ye], Ye, Adj[Ye], Ye, 
     Adj[Ye], Ye] + 18*trace[Ybc, Yb, Ybc, Yb, Ybc, Yb]*Zeta[3] + 
    6*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3])}

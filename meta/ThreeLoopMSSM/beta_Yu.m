{(-13*g1^2*Yt)/15 - 3*g2^2*Yt - (16*g3^2*Yt)/3 + MatMul[Yt, Ybc, Yb] + 
  3*MatMul[Yt, Ytc, Yt] + 3*Yt*trace[Ytc, Yt], 
 (2743*g1^4*Yt)/450 + g1^2*g2^2*Yt + (15*g2^4*Yt)/2 + (136*g1^2*g3^2*Yt)/45 + 
  8*g2^2*g3^2*Yt - (16*g3^4*Yt)/9 + (2*g1^2*MatMul[Yt, Ybc, Yb])/5 + 
  (2*g1^2*MatMul[Yt, Ytc, Yt])/5 + 6*g2^2*MatMul[Yt, Ytc, Yt] - 
  2*MatMul[Yt, Ybc, Yb, Ybc, Yb] - 2*MatMul[Yt, Ybc, Yb, Ytc, Yt] - 
  4*MatMul[Yt, Ytc, Yt, Ytc, Yt] + (4*g1^2*Yt*trace[Ytc, Yt])/5 + 
  16*g3^2*Yt*trace[Ytc, Yt] - 9*MatMul[Yt, Ytc, Yt]*trace[Ytc, Yt] + 
  MatMul[Yt, Ybc, Yb]*(-3*trace[Ybc, Yb] - trace[Adj[Ye], Ye]) + 
  Yt*(-3*trace[Ytc, Yt, Ybc, Yb] - 9*trace[Ytc, Yt, Ytc, Yt]), 
 (-8*g1^2*g2^2*g3^2*Yt)/5 - 4*g2^2*g3^2*MatMul[Yt, Ybc, Yb] + 
  (64*g3^2*MatMul[Yt, Ybc, Yb, Ybc, Yb])/3 + 
  (64*g3^2*MatMul[Yt, Ybc, Yb, Ytc, Yt])/3 + 
  (10*g1^2*MatMul[Yt, Ytc, Yt, Ytc, Yt])/3 + 
  6*g2^2*MatMul[Yt, Ytc, Yt, Ytc, Yt] + 
  (128*g3^2*MatMul[Yt, Ytc, Yt, Ytc, Yt])/3 + 
  6*MatMul[Yt, Ybc, Yb, Ybc, Yb, Ytc, Yt] + 
  4*MatMul[Yt, Ybc, Yb, Ytc, Yt, Ybc, Yb] - 
  2*MatMul[Yt, Ybc, Yb, Ytc, Yt, Ytc, Yt] + 
  2*MatMul[Yt, Ytc, Yt, Ybc, Yb, Ytc, Yt] + 12*MatMul[Yt, Ytc, Yt, Ytc, Yt]*
   trace[Ytc, Yt] + MatMul[Yt, Ybc, Yb, Ybc, Yb]*
   (6*trace[Ybc, Yb] + 2*trace[Adj[Ye], Ye]) + MatMul[Yt, Ybc, Yb, Ytc, Yt]*
   (12*trace[Ybc, Yb] - 6*trace[Ytc, Yt] + 4*trace[Adj[Ye], Ye]) + 
  g2^2*MatMul[Yt, Ybc, Yb]*(18*trace[Ybc, Yb] + 6*trace[Adj[Ye], Ye]) + 
  MatMul[Yt, Ytc, Yt]*(-27*trace[Ytc, Yt]^2 + 18*trace[Ytc, Yt, Ybc, Yb] + 
    54*trace[Ytc, Yt, Ytc, Yt]) + MatMul[Yt, Ybc, Yb]*
   (-9*trace[Ybc, Yb]^2 - 6*trace[Ybc, Yb]*trace[Adj[Ye], Ye] - 
    trace[Adj[Ye], Ye]^2 + 18*trace[Ybc, Yb, Ybc, Yb] + 
    6*trace[Ytc, Yt, Ybc, Yb] + 6*trace[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  g3^4*MatMul[Yt, Ytc, Yt]*(8 - 272*Zeta[3]) + 
  g2^4*g3^2*Yt*(140 - 216*Zeta[3]) + g2^2*g3^4*Yt*(68 - 144*Zeta[3]) + 
  g3^4*MatMul[Yt, Ybc, Yb]*(8/3 - (272*Zeta[3])/3) + 
  g1^4*g3^2*Yt*(5308/225 - (1144*Zeta[3])/25) + 
  g2^4*MatMul[Yt, Ytc, Yt]*(-219/4 - (81*Zeta[3])/2) + 
  g1^2*g3^4*Yt*(436/45 - (176*Zeta[3])/5) + g2^4*MatMul[Yt, Ybc, Yb]*
   (-45/4 - (63*Zeta[3])/2) + g1^6*Yt*(352097/6750 - (2587*Zeta[3])/125) + 
  g2^2*MatMul[Yt, Ybc, Yb, Ytc, Yt]*(9 - 18*Zeta[3]) + 
  g1^2*g2^4*Yt*(17/2 - (81*Zeta[3])/5) + 
  g1^4*g2^2*Yt*(379/50 - (351*Zeta[3])/25) + g1^2*g3^2*MatMul[Yt, Ytc, Yt]*
   (-212/15 - (32*Zeta[3])/5) + g1^4*MatMul[Yt, Ytc, Yt]*
   (-8561/300 - (117*Zeta[3])/50) + g1^2*MatMul[Yt, Ybc, Yb, Ybc, Yb]*
   (7/15 - (6*Zeta[3])/5) + g1^4*MatMul[Yt, Ybc, Yb]*
   (-633/100 + (7*Zeta[3])/30) + 6*MatMul[Yt, Ybc, Yb, Ybc, Yb, Ybc, Yb]*
   Zeta[3] + g1^2*g2^2*MatMul[Yt, Ybc, Yb]*(-41/10 + (9*Zeta[3])/5) + 
  g1^2*MatMul[Yt, Ybc, Yb, Ytc, Yt]*(19/15 + (18*Zeta[3])/5) + 
  g1^2*g3^2*MatMul[Yt, Ybc, Yb]*(-76/15 + (128*Zeta[3])/15) + 
  g2^2*MatMul[Yt, Ybc, Yb, Ybc, Yb]*(-3 + 18*Zeta[3]) + 
  MatMul[Yt, Ytc, Yt, Ytc, Yt, Ytc, Yt]*(6 + 18*Zeta[3]) + 
  g1^2*g2^2*MatMul[Yt, Ytc, Yt]*(-193/10 + (123*Zeta[3])/5) + 
  g2^2*g3^2*MatMul[Yt, Ytc, Yt]*(-92 + 96*Zeta[3]) + 
  g2^6*Yt*(345/2 + 315*Zeta[3]) + g3^6*Yt*(5440/27 + 640*Zeta[3]) + 
  g3^2*MatMul[Yt, Ybc, Yb]*(-8*trace[Ybc, Yb] + 8*trace[Adj[Ye], Ye] + 
    48*trace[Ybc, Yb]*Zeta[3]) + g2^4*Yt*(-45*trace[Ybc, Yb] - 
    (315*trace[Ytc, Yt])/4 - 15*trace[Adj[Ye], Ye] - 
    (189*trace[Ytc, Yt]*Zeta[3])/2) + g2^2*MatMul[Yt, Ytc, Yt]*
   (45*trace[Ytc, Yt] - 54*trace[Ytc, Yt]*Zeta[3]) + 
  g3^4*Yt*((-160*trace[Ybc, Yb])/3 - (320*trace[Ytc, Yt])/3 - 
    16*trace[Ytc, Yt]*Zeta[3]) + g1^4*Yt*((-91*trace[Ybc, Yb])/15 - 
    (171*trace[Ytc, Yt])/4 - (39*trace[Adj[Ye], Ye])/5 - 
    (13*trace[Ytc, Yt]*Zeta[3])/10) + g1^2*MatMul[Yt, Ytc, Yt]*
   (9*trace[Ytc, Yt] + (18*trace[Ytc, Yt]*Zeta[3])/5) + 
  g1^2*g2^2*Yt*((-57*trace[Ytc, Yt])/10 + (63*trace[Ytc, Yt]*Zeta[3])/5) + 
  g1^2*g3^2*Yt*((-124*trace[Ytc, Yt])/3 + (208*trace[Ytc, Yt]*Zeta[3])/5) + 
  g2^2*g3^2*Yt*(-132*trace[Ytc, Yt] + 144*trace[Ytc, Yt]*Zeta[3]) + 
  g3^2*MatMul[Yt, Ytc, Yt]*(-24*trace[Ytc, Yt] + 
    144*trace[Ytc, Yt]*Zeta[3]) + g1^2*MatMul[Yt, Ybc, Yb]*
   ((16*trace[Ybc, Yb])/5 - (8*trace[Adj[Ye], Ye])/5 - 
    (24*trace[Ybc, Yb]*Zeta[3])/5 + (12*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g3^2*Yt*(24*trace[Ytc, Yt, Ybc, Yb] + 72*trace[Ytc, Yt, Ytc, Yt] - 
    48*trace[Ytc, Yt, Ybc, Yb]*Zeta[3] - 144*trace[Ytc, Yt, Ytc, Yt]*
     Zeta[3]) + g1^2*Yt*((6*trace[Ytc, Yt, Ybc, Yb])/5 + 
    (57*trace[Ytc, Yt, Ytc, Yt])/5 + (6*trace[Ytc, Yt, Ybc, Yb]*Zeta[3])/5 - 
    (18*trace[Ytc, Yt, Ytc, Yt]*Zeta[3])/5) + 
  g2^2*Yt*(18*trace[Ytc, Yt, Ybc, Yb] + 9*trace[Ytc, Yt, Ytc, Yt] + 
    54*trace[Ytc, Yt, Ytc, Yt]*Zeta[3]) + 
  Yt*(18*trace[Ybc, Yb]*trace[Ytc, Yt, Ybc, Yb] + 
    6*trace[Adj[Ye], Ye]*trace[Ytc, Yt, Ybc, Yb] + 
    54*trace[Ytc, Yt]*trace[Ytc, Yt, Ytc, Yt] + 
    9*trace[Ybc, Yb, Ytc, Yt, Ybc, Yb] + 3*trace[Ytc, Yt, Ytc, Yt, Ytc, Yt] + 
    18*trace[Ytc, Yt, Ytc, Yt, Ytc, Yt]*Zeta[3])}

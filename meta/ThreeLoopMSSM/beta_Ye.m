{(-9*g1^2*Ye)/5 - 3*g2^2*Ye + 3*MatMul[Ye, Adj[Ye], Ye] + 
  Ye*(3*trace[Ybc, Yb] + trace[Adj[Ye], Ye]), 
 (27*g1^4*Ye)/2 + (9*g1^2*g2^2*Ye)/5 + (15*g2^4*Ye)/2 + 
  6*g2^2*MatMul[Ye, Adj[Ye], Ye] - 4*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  16*g3^2*Ye*trace[Ybc, Yb] + MatMul[Ye, Adj[Ye], Ye]*
   (-9*trace[Ybc, Yb] - 3*trace[Adj[Ye], Ye]) + 
  g1^2*Ye*((-2*trace[Ybc, Yb])/5 + (6*trace[Adj[Ye], Ye])/5) + 
  Ye*(-9*trace[Ybc, Yb, Ybc, Yb] - 3*trace[Ytc, Yt, Ybc, Yb] - 
    3*trace[Adj[Ye], Ye, Adj[Ye], Ye]), 
 (54*g1^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye])/5 + 
  6*g2^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye]*(12*trace[Ybc, Yb] + 
    4*trace[Adj[Ye], Ye]) + MatMul[Ye, Adj[Ye], Ye]*
   (-27*trace[Ybc, Yb]^2 - 18*trace[Ybc, Yb]*trace[Adj[Ye], Ye] - 
    3*trace[Adj[Ye], Ye]^2 + 54*trace[Ybc, Yb, Ybc, Yb] + 
    18*trace[Ytc, Yt, Ybc, Yb] + 18*trace[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  g2^4*g3^2*Ye*(180 - 216*Zeta[3]) + 
  g1^4*g3^2*Ye*(396/5 - (2376*Zeta[3])/25) + 
  g1^6*Ye*(24993/250 - (5373*Zeta[3])/125) + g2^4*MatMul[Ye, Adj[Ye], Ye]*
   (-219/4 - (81*Zeta[3])/2) + g1^4*g2^2*Ye*(837/50 - (729*Zeta[3])/25) + 
  g1^2*g2^4*Ye*(9/2 - (81*Zeta[3])/5) + g1^4*MatMul[Ye, Adj[Ye], Ye]*
   (-5751/100 - (729*Zeta[3])/50) + 
  MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye]*(6 + 18*Zeta[3]) + 
  g1^2*g2^2*MatMul[Ye, Adj[Ye], Ye]*(-351/10 + (243*Zeta[3])/5) + 
  g2^6*Ye*(345/2 + 315*Zeta[3]) + g3^4*Ye*((-160*trace[Ybc, Yb])/3 - 
    16*trace[Ybc, Yb]*Zeta[3]) + g1^2*g3^2*Ye*((-284*trace[Ybc, Yb])/15 + 
    (112*trace[Ybc, Yb]*Zeta[3])/5) + 
  g2^2*g3^2*Ye*(-132*trace[Ybc, Yb] + 144*trace[Ybc, Yb]*Zeta[3]) + 
  g3^2*MatMul[Ye, Adj[Ye], Ye]*(-96*trace[Ybc, Yb] + 
    144*trace[Ybc, Yb]*Zeta[3]) + g2^4*Ye*((-315*trace[Ybc, Yb])/4 - 
    45*trace[Ytc, Yt] - (105*trace[Adj[Ye], Ye])/4 - 
    (189*trace[Ybc, Yb]*Zeta[3])/2 - (63*trace[Adj[Ye], Ye]*Zeta[3])/2) + 
  g2^2*MatMul[Ye, Adj[Ye], Ye]*(45*trace[Ybc, Yb] + 15*trace[Adj[Ye], Ye] - 
    54*trace[Ybc, Yb]*Zeta[3] - 18*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g1^4*Ye*((-301*trace[Ybc, Yb])/12 - (117*trace[Ytc, Yt])/5 - 
    (873*trace[Adj[Ye], Ye])/20 - (77*trace[Ybc, Yb]*Zeta[3])/50 + 
    (81*trace[Adj[Ye], Ye]*Zeta[3])/50) + g1^2*MatMul[Ye, Adj[Ye], Ye]*
   ((147*trace[Ybc, Yb])/5 + (9*trace[Adj[Ye], Ye])/5 - 
    (18*trace[Ybc, Yb]*Zeta[3])/5 + (54*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*g2^2*Ye*((-3*trace[Ybc, Yb])/10 - (81*trace[Adj[Ye], Ye])/10 - 
    9*trace[Ybc, Yb]*Zeta[3] + (81*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g3^2*Ye*(72*trace[Ybc, Yb, Ybc, Yb] + 24*trace[Ytc, Yt, Ybc, Yb] - 
    144*trace[Ybc, Yb, Ybc, Yb]*Zeta[3] - 48*trace[Ytc, Yt, Ybc, Yb]*
     Zeta[3]) + g1^2*Ye*(3*trace[Ybc, Yb, Ybc, Yb] - 
    (12*trace[Ytc, Yt, Ybc, Yb])/5 + 9*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    (54*trace[Ybc, Yb, Ybc, Yb]*Zeta[3])/5 + 
    (42*trace[Ytc, Yt, Ybc, Yb]*Zeta[3])/5 - 
    (54*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3])/5) + 
  g2^2*Ye*(9*trace[Ybc, Yb, Ybc, Yb] + 18*trace[Ytc, Yt, Ybc, Yb] + 
    3*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 54*trace[Ybc, Yb, Ybc, Yb]*Zeta[3] + 
    18*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3]) + 
  Ye*(54*trace[Ybc, Yb]*trace[Ybc, Yb, Ybc, Yb] + 
    18*trace[Adj[Ye], Ye]*trace[Ybc, Yb, Ybc, Yb] + 
    18*trace[Ytc, Yt]*trace[Ytc, Yt, Ybc, Yb] + 18*trace[Ybc, Yb]*
     trace[Adj[Ye], Ye, Adj[Ye], Ye] + 6*trace[Adj[Ye], Ye]*
     trace[Adj[Ye], Ye, Adj[Ye], Ye] + 3*trace[Ybc, Yb, Ybc, Yb, Ybc, Yb] + 
    9*trace[Ytc, Yt, Ytc, Yt, Ybc, Yb] + trace[Adj[Ye], Ye, Adj[Ye], Ye, 
     Adj[Ye], Ye] + 18*trace[Ybc, Yb, Ybc, Yb, Ybc, Yb]*Zeta[3] + 
    6*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3])}

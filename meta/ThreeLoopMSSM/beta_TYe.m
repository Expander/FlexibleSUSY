{(-9*g1^2*he)/5 - 3*g2^2*he + (18*g1^2*M1*Ye)/5 + 6*g2^2*M2*Ye + 
  5*MatMul[he, Adj[Ye], Ye] + 4*MatMul[Ye, Adj[Ye], he] + 
  Ye*(6*trace[hb, Adj[Yb]] + 2*trace[he, Adj[Ye]]) + 
  he*(3*trace[Adj[Yb], Yb] + trace[Adj[Ye], Ye]), 
 (27*g1^4*he)/2 + (9*g1^2*g2^2*he)/5 + (15*g2^4*he)/2 - 54*g1^4*M1*Ye - 
  (18*g1^2*g2^2*M1*Ye)/5 - (18*g1^2*g2^2*M2*Ye)/5 - 30*g2^4*M2*Ye - 
  (6*g1^2*MatMul[he, Adj[Ye], Ye])/5 + 12*g2^2*MatMul[he, Adj[Ye], Ye] + 
  (6*g1^2*MatMul[Ye, Adj[Ye], he])/5 + 6*g2^2*MatMul[Ye, Adj[Ye], he] - 
  12*g2^2*M2*MatMul[Ye, Adj[Ye], Ye] - 
  6*MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye] - 
  8*MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye] - 
  6*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he] + 32*g3^2*Ye*trace[hb, Adj[Yb]] + 
  MatMul[Ye, Adj[Ye], Ye]*(-18*trace[hb, Adj[Yb]] - 6*trace[he, Adj[Ye]]) + 
  g1^2*Ye*((-4*trace[hb, Adj[Yb]])/5 + (12*trace[he, Adj[Ye]])/5) + 
  16*g3^2*he*trace[Adj[Yb], Yb] - 32*g3^2*M3*Ye*trace[Adj[Yb], Yb] + 
  MatMul[he, Adj[Ye], Ye]*(-15*trace[Adj[Yb], Yb] - 5*trace[Adj[Ye], Ye]) + 
  MatMul[Ye, Adj[Ye], he]*(-12*trace[Adj[Yb], Yb] - 4*trace[Adj[Ye], Ye]) + 
  g1^2*M1*Ye*((4*trace[Adj[Yb], Yb])/5 - (12*trace[Adj[Ye], Ye])/5) + 
  g1^2*he*((-2*trace[Adj[Yb], Yb])/5 + (6*trace[Adj[Ye], Ye])/5) + 
  Ye*(-6*trace[hb, Adj[Yt], Yt, Adj[Yb]] - 
    36*trace[Adj[Yb], hb, Adj[Yb], Yb] - 12*trace[Adj[Ye], he, Adj[Ye], Ye] - 
    6*trace[Adj[Yt], ht, Adj[Yb], Yb]) + 
  he*(-9*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
    3*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 3*trace[Adj[Yt], Yt, Adj[Yb], Yb]), 
 (108*g1^2*MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye])/5 + 
  12*g2^2*MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye] - 
  (108*g1^2*M1*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye])/5 - 
  12*g2^2*M2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye]*(24*trace[hb, Adj[Yb]] + 
    8*trace[he, Adj[Ye]]) + MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye]*
   (18*trace[Adj[Yb], Yb] + 6*trace[Adj[Ye], Ye]) + 
  MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he]*(18*trace[Adj[Yb], Yb] + 
    6*trace[Adj[Ye], Ye]) + MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye]*
   (24*trace[Adj[Yb], Yb] + 8*trace[Adj[Ye], Ye]) + 
  MatMul[Ye, Adj[Ye], Ye]*(-108*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
    36*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb] - 36*trace[hb, Adj[Yb]]*
     trace[Adj[Ye], Ye] - 12*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] + 
    36*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
    216*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
    72*trace[Adj[Ye], he, Adj[Ye], Ye] + 
    36*trace[Adj[Yt], ht, Adj[Yb], Yb]) + MatMul[Ye, Adj[Ye], he]*
   (-36*trace[Adj[Yb], Yb]^2 - 24*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
    4*trace[Adj[Ye], Ye]^2 + 72*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    24*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    24*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + MatMul[he, Adj[Ye], Ye]*
   (-45*trace[Adj[Yb], Yb]^2 - 30*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
    5*trace[Adj[Ye], Ye]^2 + 90*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    30*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    30*trace[Adj[Yt], Yt, Adj[Yb], Yb]) + g2^6*M2*Ye*(-1035 - 1890*Zeta[3]) + 
  g2^4*g3^2*he*(180 - 216*Zeta[3]) + g1^2*g2^2*M1*MatMul[Ye, Adj[Ye], Ye]*
   (351/5 - (486*Zeta[3])/5) + g1^2*g2^2*M2*MatMul[Ye, Adj[Ye], Ye]*
   (351/5 - (486*Zeta[3])/5) + g1^4*g3^2*he*(396/5 - (2376*Zeta[3])/25) + 
  g2^4*MatMul[Ye, Adj[Ye], he]*(-66 - 72*Zeta[3]) + 
  g2^4*MatMul[he, Adj[Ye], Ye]*(-393/4 - (99*Zeta[3])/2) + 
  g1^6*he*(24993/250 - (5373*Zeta[3])/125) + g1^4*MatMul[he, Adj[Ye], Ye]*
   (-8757/100 - (1539*Zeta[3])/50) + 
  g1^4*g2^2*he*(837/50 - (729*Zeta[3])/25) + 
  g2^2*MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye]*(15 - 18*Zeta[3]) + 
  g1^2*g2^4*he*(9/2 - (81*Zeta[3])/5) + g1^4*MatMul[Ye, Adj[Ye], he]*
   (-2124/25 - (324*Zeta[3])/25) + g1^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he]*
   (99/5 - (54*Zeta[3])/5) + g1^2*MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye]*
   (63/5 + (54*Zeta[3])/5) + g2^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he]*
   (3 + 18*Zeta[3]) + MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], he]*
   (6 + 24*Zeta[3]) + MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye]*
   (12 + 30*Zeta[3]) + g1^2*g2^4*M1*Ye*(-9 + (162*Zeta[3])/5) + 
  MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye, Adj[Ye], Ye]*(12 + 36*Zeta[3]) + 
  MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he, Adj[Ye], Ye]*(12 + 36*Zeta[3]) + 
  g1^4*g2^2*M2*Ye*(-837/25 + (1458*Zeta[3])/25) + 
  g1^4*M1*MatMul[Ye, Adj[Ye], Ye]*(5751/25 + (1458*Zeta[3])/25) + 
  g1^2*g2^2*MatMul[Ye, Adj[Ye], he]*(-216/5 + (324*Zeta[3])/5) + 
  g1^2*g2^4*M2*Ye*(-18 + (324*Zeta[3])/5) + g1^2*g2^2*MatMul[he, Adj[Ye], Ye]*
   (-621/10 + 81*Zeta[3]) + g1^4*g2^2*M1*Ye*(-1674/25 + (2916*Zeta[3])/25) + 
  g2^4*M2*MatMul[Ye, Adj[Ye], Ye]*(219 + 162*Zeta[3]) + 
  g1^4*g3^2*M3*Ye*(-792/5 + (4752*Zeta[3])/25) + 
  g1^6*M1*Ye*(-74979/125 + (32238*Zeta[3])/125) + 
  g2^6*he*(345/2 + 315*Zeta[3]) + g1^4*g3^2*M1*Ye*
   (-1584/5 + (9504*Zeta[3])/25) + g2^4*g3^2*M3*Ye*(-360 + 432*Zeta[3]) + 
  g2^4*g3^2*M2*Ye*(-720 + 864*Zeta[3]) + 
  g3^4*Ye*((-320*trace[hb, Adj[Yb]])/3 - 32*trace[hb, Adj[Yb]]*Zeta[3]) + 
  g1^2*g3^2*Ye*((-568*trace[hb, Adj[Yb]])/15 + 
    (224*trace[hb, Adj[Yb]]*Zeta[3])/5) + 
  g2^2*g3^2*Ye*(-264*trace[hb, Adj[Yb]] + 288*trace[hb, Adj[Yb]]*Zeta[3]) + 
  g3^2*MatMul[Ye, Adj[Ye], Ye]*(-192*trace[hb, Adj[Yb]] + 
    288*trace[hb, Adj[Yb]]*Zeta[3]) + 
  g2^4*Ye*((-315*trace[hb, Adj[Yb]])/2 - (105*trace[he, Adj[Ye]])/2 - 
    90*trace[ht, Adj[Yt]] - 189*trace[hb, Adj[Yb]]*Zeta[3] - 
    63*trace[he, Adj[Ye]]*Zeta[3]) + g2^2*MatMul[Ye, Adj[Ye], Ye]*
   (90*trace[hb, Adj[Yb]] + 30*trace[he, Adj[Ye]] - 
    108*trace[hb, Adj[Yb]]*Zeta[3] - 36*trace[he, Adj[Ye]]*Zeta[3]) + 
  g1^4*Ye*((-301*trace[hb, Adj[Yb]])/6 - (873*trace[he, Adj[Ye]])/10 - 
    (234*trace[ht, Adj[Yt]])/5 - (77*trace[hb, Adj[Yb]]*Zeta[3])/25 + 
    (81*trace[he, Adj[Ye]]*Zeta[3])/25) + g1^2*MatMul[Ye, Adj[Ye], Ye]*
   ((294*trace[hb, Adj[Yb]])/5 + (18*trace[he, Adj[Ye]])/5 - 
    (36*trace[hb, Adj[Yb]]*Zeta[3])/5 + (108*trace[he, Adj[Ye]]*Zeta[3])/5) + 
  g1^2*g2^2*Ye*((-3*trace[hb, Adj[Yb]])/5 - (81*trace[he, Adj[Ye]])/5 - 
    18*trace[hb, Adj[Yb]]*Zeta[3] + (162*trace[he, Adj[Ye]]*Zeta[3])/5) + 
  g3^2*M3*MatMul[Ye, Adj[Ye], Ye]*(192*trace[Adj[Yb], Yb] - 
    288*trace[Adj[Yb], Yb]*Zeta[3]) + g2^2*g3^2*M2*Ye*
   (264*trace[Adj[Yb], Yb] - 288*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g2^2*g3^2*M3*Ye*(264*trace[Adj[Yb], Yb] - 288*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g1^2*g3^2*M1*Ye*((568*trace[Adj[Yb], Yb])/15 - 
    (224*trace[Adj[Yb], Yb]*Zeta[3])/5) + 
  g1^2*g3^2*M3*Ye*((568*trace[Adj[Yb], Yb])/15 - 
    (224*trace[Adj[Yb], Yb]*Zeta[3])/5) + 
  g3^4*he*((-160*trace[Adj[Yb], Yb])/3 - 16*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g1^2*g3^2*he*((-284*trace[Adj[Yb], Yb])/15 + 
    (112*trace[Adj[Yb], Yb]*Zeta[3])/5) + 
  g3^4*M3*Ye*((640*trace[Adj[Yb], Yb])/3 + 64*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g2^2*g3^2*he*(-132*trace[Adj[Yb], Yb] + 144*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g3^2*MatMul[Ye, Adj[Ye], he]*(-128*trace[Adj[Yb], Yb] + 
    192*trace[Adj[Yb], Yb]*Zeta[3]) + g3^2*MatMul[he, Adj[Ye], Ye]*
   (-160*trace[Adj[Yb], Yb] + 240*trace[Adj[Yb], Yb]*Zeta[3]) + 
  g2^2*MatMul[he, Adj[Ye], Ye]*(72*trace[Adj[Yb], Yb] + 
    24*trace[Adj[Ye], Ye] - 108*trace[Adj[Yb], Yb]*Zeta[3] - 
    36*trace[Adj[Ye], Ye]*Zeta[3]) + g1^2*g2^2*M1*Ye*
   ((3*trace[Adj[Yb], Yb])/5 + (81*trace[Adj[Ye], Ye])/5 + 
    18*trace[Adj[Yb], Yb]*Zeta[3] - (162*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*g2^2*M2*Ye*((3*trace[Adj[Yb], Yb])/5 + (81*trace[Adj[Ye], Ye])/5 + 
    18*trace[Adj[Yb], Yb]*Zeta[3] - (162*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g2^4*he*((-315*trace[Adj[Yb], Yb])/4 - (105*trace[Adj[Ye], Ye])/4 - 
    45*trace[Adj[Yt], Yt] - (189*trace[Adj[Yb], Yb]*Zeta[3])/2 - 
    (63*trace[Adj[Ye], Ye]*Zeta[3])/2) + g1^2*M1*MatMul[Ye, Adj[Ye], Ye]*
   ((-294*trace[Adj[Yb], Yb])/5 - (18*trace[Adj[Ye], Ye])/5 + 
    (36*trace[Adj[Yb], Yb]*Zeta[3])/5 - (108*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g2^2*MatMul[Ye, Adj[Ye], he]*(63*trace[Adj[Yb], Yb] + 
    21*trace[Adj[Ye], Ye] - 54*trace[Adj[Yb], Yb]*Zeta[3] - 
    18*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g1^4*M1*Ye*((301*trace[Adj[Yb], Yb])/3 + (873*trace[Adj[Ye], Ye])/5 + 
    (468*trace[Adj[Yt], Yt])/5 + (154*trace[Adj[Yb], Yb]*Zeta[3])/25 - 
    (162*trace[Adj[Ye], Ye]*Zeta[3])/25) + 
  g1^4*he*((-301*trace[Adj[Yb], Yb])/12 - (873*trace[Adj[Ye], Ye])/20 - 
    (117*trace[Adj[Yt], Yt])/5 - (77*trace[Adj[Yb], Yb]*Zeta[3])/50 + 
    (81*trace[Adj[Ye], Ye]*Zeta[3])/50) + g1^2*MatMul[Ye, Adj[Ye], he]*
   ((187*trace[Adj[Yb], Yb])/5 + (9*trace[Adj[Ye], Ye])/5 - 
    (78*trace[Adj[Yb], Yb]*Zeta[3])/5 + (54*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*g2^2*he*((-3*trace[Adj[Yb], Yb])/10 - (81*trace[Adj[Ye], Ye])/10 - 
    9*trace[Adj[Yb], Yb]*Zeta[3] + (81*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  g1^2*MatMul[he, Adj[Ye], Ye]*((254*trace[Adj[Yb], Yb])/5 + 
    (18*trace[Adj[Ye], Ye])/5 + (24*trace[Adj[Yb], Yb]*Zeta[3])/5 + 
    (108*trace[Adj[Ye], Ye]*Zeta[3])/5) + g2^2*M2*MatMul[Ye, Adj[Ye], Ye]*
   (-90*trace[Adj[Yb], Yb] - 30*trace[Adj[Ye], Ye] + 
    108*trace[Adj[Yb], Yb]*Zeta[3] + 36*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g2^4*M2*Ye*(315*trace[Adj[Yb], Yb] + 105*trace[Adj[Ye], Ye] + 
    180*trace[Adj[Yt], Yt] + 378*trace[Adj[Yb], Yb]*Zeta[3] + 
    126*trace[Adj[Ye], Ye]*Zeta[3]) + 
  g2^2*Ye*(36*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
    36*trace[Adj[Yb], hb, Adj[Yb], Yb] + 12*trace[Adj[Ye], he, Adj[Ye], Ye] + 
    36*trace[Adj[Yt], ht, Adj[Yb], Yb] + 216*trace[Adj[Yb], hb, Adj[Yb], Yb]*
     Zeta[3] + 72*trace[Adj[Ye], he, Adj[Ye], Ye]*Zeta[3]) + 
  g2^2*M2*Ye*(-18*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
    6*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 36*trace[Adj[Yt], Yt, Adj[Yb], Yb] - 
    108*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] - 
    36*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3]) + 
  g2^2*he*(9*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    3*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 18*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    54*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
    18*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3]) + 
  g3^2*Ye*(48*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
    288*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
    48*trace[Adj[Yt], ht, Adj[Yb], Yb] - 96*trace[hb, Adj[Yt], Yt, Adj[Yb]]*
     Zeta[3] - 576*trace[Adj[Yb], hb, Adj[Yb], Yb]*Zeta[3] - 
    96*trace[Adj[Yt], ht, Adj[Yb], Yb]*Zeta[3]) + 
  g1^2*Ye*((-24*trace[hb, Adj[Yt], Yt, Adj[Yb]])/5 + 
    12*trace[Adj[Yb], hb, Adj[Yb], Yb] + 36*trace[Adj[Ye], he, Adj[Ye], Ye] - 
    (24*trace[Adj[Yt], ht, Adj[Yb], Yb])/5 + 
    (84*trace[hb, Adj[Yt], Yt, Adj[Yb]]*Zeta[3])/5 + 
    (216*trace[Adj[Yb], hb, Adj[Yb], Yb]*Zeta[3])/5 - 
    (216*trace[Adj[Ye], he, Adj[Ye], Ye]*Zeta[3])/5 + 
    (84*trace[Adj[Yt], ht, Adj[Yb], Yb]*Zeta[3])/5) + 
  g3^2*he*(72*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    24*trace[Adj[Yt], Yt, Adj[Yb], Yb] - 144*trace[Adj[Yb], Yb, Adj[Yb], Yb]*
     Zeta[3] - 48*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3]) + 
  g1^2*M1*Ye*(-6*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
    18*trace[Adj[Ye], Ye, Adj[Ye], Ye] + (24*trace[Adj[Yt], Yt, Adj[Yb], Yb])/
     5 - (108*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3])/5 + 
    (108*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3])/5 - 
    (84*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3])/5) + 
  g1^2*he*(3*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    9*trace[Adj[Ye], Ye, Adj[Ye], Ye] - (12*trace[Adj[Yt], Yt, Adj[Yb], Yb])/
     5 + (54*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3])/5 - 
    (54*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3])/5 + 
    (42*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3])/5) + 
  g3^2*M3*Ye*(-144*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
    48*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 288*trace[Adj[Yb], Yb, Adj[Yb], Yb]*
     Zeta[3] + 96*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3]) + 
  Ye*(36*trace[Adj[Yt], Yt]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
    216*trace[Adj[Yb], Yb]*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
    72*trace[Adj[Ye], Ye]*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
    108*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    36*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    72*trace[Adj[Yb], Yb]*trace[Adj[Ye], he, Adj[Ye], Ye] + 
    24*trace[Adj[Ye], Ye]*trace[Adj[Ye], he, Adj[Ye], Ye] + 
    36*trace[hb, Adj[Yb]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    12*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    36*trace[Adj[Yt], Yt]*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
    36*trace[ht, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    18*trace[hb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb]] + 
    18*trace[Adj[Yb], Yb, Adj[Yb], hb, Adj[Yb], Yb] + 
    6*trace[Adj[Ye], Ye, Adj[Ye], he, Adj[Ye], Ye] + 
    18*trace[Adj[Yt], ht, Adj[Yt], Yt, Adj[Yb], Yb] + 
    18*trace[Adj[Yt], Yt, Adj[Yt], ht, Adj[Yb], Yb] + 
    108*trace[Adj[Yb], Yb, Adj[Yb], hb, Adj[Yb], Yb]*Zeta[3] + 
    36*trace[Adj[Ye], Ye, Adj[Ye], he, Adj[Ye], Ye]*Zeta[3]) + 
  he*(54*trace[Adj[Yb], Yb]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    18*trace[Adj[Ye], Ye]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    18*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    6*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    18*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    3*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb] + 
    trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
    9*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], Yb] + 
    18*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
    6*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3])}

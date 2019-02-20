(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

(* loop functions from Appendix A of arxiv:1407.4081 *)

F1[x_] := x Log[x^2] / (x^2 - 1);

F2[x_] := 6 x^2 (2 - 2 x^2 + (1 + x^2) Log[x^2]) / (x^2 - 1)^3;

F3[x_] := 2 x (5 (1 - x^2) + (1 + 4 x^2) Log[x^2]) / (3 (x^2 - 1)^2);

F4[x_] := 2 x (x^2 - 1 - Log[x^2]) / (x^2 - 1)^2;

F5[x_] := 3 x (1 - x^4 + 2 x^2 Log[x^2]) / (1 - x^2)^3;

F6[x_] := (x^2 - 3)/(4 (1 - x^2)) + x^2 (x^2 - 2) Log[x^2] / (2 (1 - x^2)^2);

F7[x_] := -3 (x^4 - 6 x^2 + 1)/(2 (x^2 - 1)^2) + 3 x^4 (x^2 - 3) Log[x^2] / (x^2 - 1)^3;

F8[x1_, x2_] := -2 + 2/(x1^2 - x2^2) (x1^4 Log[x1^2]/(x1^2 - 1) - x2^4 Log[x2^2]/(x2^2 - 1));

F9[x1_, x2_] := 2/(x1^2 - x2^2) (x1^2 Log[x1^2]/(x1^2 - 1) - x2^2 Log[x2^2]/(x2^2 - 1));

f[r_] := F5[r];

g[r_] := F7[r];

f1[r_] := 6 (r^2 + 3) r^2 / (7 (r^2 - 1)^2) + 6 (r^2 - 5) r^4 Log[r^2] / (7 (r^2 - 1)^3);

f2[r_] := 2 (r^2 + 11) r^2 / (9 (r^2 - 1)^2) + 2 (5 r^2 - 17) r^4 Log[r^2] / (9 (r^2 - 1)^3);

f3[r_] := 2 (r^4 + 9 r^2 + 2) / (3 (r^2 - 1)^2) + 2 (r^4 - 7 r^2 - 6) r^2 Log[r^2] / (3 (r^2 - 1)^3);

f4[r_] := 2 (5 r^4 + 25 r^2 + 6) / (7 (r^2 - 1)^2) + 2 (r^4 - 19 r^2 - 18) r^2 Log[r^2] / (7 (r^2 - 1)^3);

f5[r1_, r2_] := (3/4) (
    (1 + (r1 + r2)^2 - r1^2 r2^2)/((r1^2 - 1) (r2^2 - 1))
    + (r1^3 (r1^2 + 1) Log[r1^2])/((r1^2 - 1)^2 (r1 - r2))
    - (r2^3 (r2^2 + 1) Log[r2^2])/((r1 - r2) (r2^2 - 1)^2)
);

f6[r1_, r2_] := (6/7) (
    (r1^2 + r2^2 + r1 r2 - r1^2 r2^2)/((r1^2 - 1) (r2^2 - 1))
    + ((r1)^5 Log[r1^2])/((r1^2 - 1)^2 (r1 - r2))
    - ((r2)^5 Log[r2^2])/((r1 - r2) (r2^2 - 1)^2)
);

f7[r1_, r2_] := 6 (
    (1 + r1 r2)/((r1^2 - 1) (r2^2 - 1))
    + ((r1)^3 Log[r1^2])/((r1^2 - 1)^2 (r1 - r2))
    - ((r2)^3 Log[r2^2])/((r1 - r2) (r2^2 - 1)^2)
);

f8[r1_, r2_] := (3/2) (
    (r1 + r2)/((r1^2 - 1) (r2^2 - 1))
    + ((r1)^4 Log[r1^2])/((r1^2 - 1)^2 (r1 - r2))
    - ((r2)^4 Log[r2^2])/((r1 - r2) (r2^2 - 1)^2)
);

(* Delta and Phi function for 2-loop threshold corrections *)

Cl2[x_] := Im[Li2[Exp[I x]]];

Li2[x_] := PolyLog[2, x];

Delta[x_, y_, z_] := x^2 + y^2 + z^2 - 2 (x y + x z + y z);

Lambda2[u_, v_] := (1 - u - v)^2 - 4 u v;

PhiPos[u_, v_] :=
    Module[{lam = Sqrt[Lambda2[u, v]]},
           1/lam (
               - Log[u] Log[v]
               + 2 Log[(1 - lam + u - v)/2] Log[(1 - lam - u + v)/2]
               - 2 Li2[(1 - lam + u - v)/2]
               - 2 Li2[(1 - lam - u + v)/2]
               + Pi^2/3
           )
          ];

PhiNeg[u_, v_] :=
    Module[{lam = Sqrt[-Lambda2[u, v]]},
           2/lam (
               + Cl2[2 ArcCos[(1 + u - v)/(2 Sqrt[u])]]
               + Cl2[2 ArcCos[(1 - u + v)/(2 Sqrt[v])]]
               + Cl2[2 ArcCos[(-1 + u + v)/(2 Sqrt[u v])]]
           )
          ];

PhiUV[u_, v_] :=
    Module[{lam = Lambda2[u, v]},
           If[lam > 0,
              Which[
                  u <= 1 && v <= 1, PhiPos[u, v],
                  u >= 1 && v/u <= 1, PhiPos[1/u, v/u]/u,
                  True, PhiPos[1/v, u/v]/v
                   ],
              PhiNeg[u, v]
             ]
          ];

(* Phi[x,y,z] function from
 * Davydychev and Tausk, Nucl. Phys. B397 (1993) 23.
 * Arguments are interpreted as squared mass parameters.
 *)
Phi[x_, y_, z_] := PhiUV[x/z, y/z];

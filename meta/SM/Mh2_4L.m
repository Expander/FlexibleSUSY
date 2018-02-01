(* effective potential from (5.3) of [1508.00912] *)
(* k = 1/(4 Pi)^2 *)
V4L[h_] := h^4 gs^6 yt^4 k^4 (a4 l^4 + a3 l^3 + a2 l^2 + a1 l + a0) /.
    { l -> Log[h^2*yt^2/(2*Q^2)] };

(* tadpole / v *)
t4L = -1/v D[V4L[h], h] /. { h -> v } /. { v -> Sqrt[2] mt/yt } // Simplify;

(* self energy *)
s4L = -V4L''[h] /. { h -> Sqrt[2] mt/yt } // Simplify;

(* Delta Mh^2 at 4-loop O(at*as^3) *)
DMh24L = t4L - s4L /. {
    a4 -> aL4, (* 1380 *)
    a3 -> aL3, (* -9144 *)
    a2 -> aL2, (* 27699.06 *)
    a1 -> aL1, (* -54056.36 *)
    a0 -> aL0  (* 59366.97 *)
} // Simplify;

(* coefficients from Eq.(4.39) *)
aL4 = 1380;

aL3 = -9144;

aL2 = 30584 - 2400 Zeta[3];

aL1 = 27680/3 Zeta[3] - 63200/9 Zeta[5] - 1547146/27 - 208/9 Pi^4 +
    640/3 Log[2]^2 (Log[2]^2 - Pi^2) + 5120 PolyLog[4,1/2];

aL0 = 13820381/270 + 1747112 Zeta[3]/45 + 1984 Zeta[5]/9 -
    40288/9 Zeta[3]^2 - 298894/1215 Pi^4 - 1780/243 Pi^6 +
    5888/135 Log[2]^5 - 5888/81 Pi^2 Log[2]^3 - 36064/405 Pi^4 Log[2] +
    78464/81 Log[2]^2 (Log[2]^2 - Pi^2) + 627712/27 PolyLog[4,1/2] -
    47104/9 PolyLog[5,1/2];

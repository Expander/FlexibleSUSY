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
    a4 -> 1380,
    a3 -> -9144,
    a2 -> 27699.06,
    a1 -> -54056.36,
    a0 -> 59366.97
} // Simplify;

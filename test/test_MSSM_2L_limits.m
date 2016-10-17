(* Script for converting MSSM 2L Higgs corrections from Fortran to
   Mathematica form and derivation of the limits s2t -> 0 and mst1 ->
   mst2.
 *)

repl = {t -> mt^2, g -> mg^2, q -> Q^2, T1 -> mst1^2, T2 -> mst2^2, T -> T2, phi -> phiExpr};
A = Xt - mu/tanb; (* notice sign convention for mu! *)
Xt = (T1 - T2) s2t/(2 mt);
del = g^2 + t^2 + T1^2 - 2 (g t + g T1 + t T1);

(* phi[x,y,z] function from Pietro *)
phiExpr[x_, y_, z_] :=
    Module[{u = x/z, v = y/z, lambda, xp, xm},
           lambda = Sqrt[(1 - u - v)^2 - 4 u v];
           xp = 1/2 (1 + (u - v) - lambda);
           xm = 1/2 (1 - (u - v) - lambda);
           1/lambda (2 Log[xp] Log[xm] - Log[u] Log[v] - 
                     2 (PolyLog[2, xp] + PolyLog[2, xm]) + Pi^2/3)
          ];

(*****************************)
(* CP-even tadpoles O(at*as) *)
(*****************************)

(* Fortran expressions from Pietro Slavich *)

F2lStr = "4*mg*mt*(1+4*c2t**2)/s2t
       $     -(2*(T1-T2)+4*mg*mt/s2t)*Log(g/q)*Log(t/q)
       $     -2*(4-s2t**2)*(T1-T2)
       $     +(4*T1*T2-s2t**2*(T1+T2)**2)/(T1-T2)*Log(T1/q)*Log(T2/q)
       $     + strF2lc(t,mg,T1,T2,s2t,c2t,q)
       $     - strF2lc(t,mg,T2,T1,-s2t,c2t,q)";

G2lStr = "5*mg/mt*s2t*(T1-T2)-10*(T1+T2-2*t)-4*g
       $     + 12*t*(Log(t/q)**2-2*Log(t/q))
       $     +(4*g-s2t*mg/mt*(T1-T2))*Log(g/q)*Log(t/q)
       $     +s2t**2*(T1+T2)*Log(T1/q)*Log(T2/q)
       $     + strG2lc(t,mg,T1,T2,s2t,q)
       $     + strG2lc(t,mg,T2,T1,-s2t,q)";

strF2lcStr = "(4*(g+t+2*T1)-s2t**2*(3*T1+T2)-4*s2t*mg*mt
       $     -16*mg*mt*T1*c2t**2/s2t/(T1-T2))*Log(T1/q)       
       $     +T1/(T1-T2)*(s2t**2*(T1+T2)-2*(2*T1-T2))*Log(T1/q)**2
       $     +2*(T1-g-t+mg*mt*s2t
       $     +2*c2t**2*mg*mt*T1/s2t/(T1-T2))*Log(g*t/q**2)*Log(T1/q)
       $     +4*mg*mt*c2t**2*(t-g)/s2t/(T1-T2)*Log(t/g)*Log(T1/q)
       $     +((2*del+4*g*t)/T1-2*mg*mt*s2t/T1*(g+t-T1)
       $     +4*c2t**2*mg*mt/T1/(T1-T2)/s2t*del)*phi(g,t,T1)";

strG2lcStr = "(4*(g+t+2*T1)+s2t**2*(T1-T2)
       $     -4*mg/mt*s2t*(t+T1))*Log(T1/q) 
       $     +(mg/mt*s2t*(5*t-g+T1)-2*(g+2*t))*Log(t/q)*Log(T1/q)
       $     +(mg/mt*s2t*(g-t+T1)-2*g)*Log(g/q)*Log(T1/q)
       $     -(2+s2t**2)*T1*Log(T1/q)**2
       $     +(2*g/T1*(g+t-T1-2*mg*mt*s2t)
       $     +mg/mt*s2t*del/T1)*phi(g,t,T1)";

exprStr = ("(" <> # <> ")")& /@ {F2lStr, G2lStr, strF2lcStr, strG2lcStr};

NoBracket[c_String] := !MemberQ[{"(", ")"}, c];

MakeFunc[head_String] :=
    head <> "(" ~~ q__?NoBracket ~~ ")" :> head <> "[" <> q <> "]";

exprInputForm = \
    StringReplace[
        StringReplace[exprStr, {"$" -> "", "**" -> "^"}],
        {MakeFunc["Log"],
         MakeFunc["strF2lc"],
         MakeFunc["strG2lc"],
         MakeFunc["phi"]}
    ];

expr = ToExpression /@ exprInputForm;

F2l = expr[[1]];
G2l = expr[[2]];
strF2lc[t_, mg_, T1_, T2_, s2t_, c2t_, q_] := Evaluate[expr[[3]]];
strG2lc[t_, mg_, T1_, T2_, s2t_, q_] := Evaluate[expr[[4]]];

k = 4 gs^2/(16 Pi^2)^2;

S1 = mt*mu/tanb*s2t*F2l;
S1 = S1/v1^2;

S2 = mt*A*s2t*F2l + 2*mt^2*G2l;
S2 = S2/v2^2;

S1 = k*S1;
S2 = k*S2;

(* CP-even tadpoles O(at*as) *)
tadpoles = {S1, S2};

(* Limits *)
(* s2t -> 0 *)

tadpolesS2t0 =
    Simplify[Limit[tadpoles, s2t -> 0] /. c2t -> 1 /. Xt -> 0];

(* s2t -> 0 and mst1 -> mst2 *)
tadpolesS2t0T1eqT2 =
    Simplify[Normal /@
             Series[tadpolesS2t0 /. phi -> phiExpr /. T1 -> T2 + eps,
                    {eps, 0, 0}] /. T2 -> T];

Print[CForm /@ FullSimplify[tadpolesS2t0]];
Print[CForm /@ Simplify[tadpolesS2t0T1eqT2]];

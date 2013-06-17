<< models/fmssm/fmssm_lattice_defs.m
<< meta/writeNRGE.m

gY := Sqrt[3/5] g1;

m2stop := {
   {m2Q[[3,3]] + vu^2 norm[Yu[[3,3]]] + 1/12 (3 g2^2 - gY^2) (vd^2-vu^2),
    vu Conjugate[TAu[[3,3]]] - vd Conjugate[Yu[[3,3]]] mu},
   {vu TAu[[3,3]] - vd Yu[[3,3]] Conjugate[mu],
    m2U[[3,3]] + vu^2 norm[Yu[[3,3]]] + 1/3 gY^2 (vd^2-vu^2)}};

hTable[mat_, nr_] := Table[Which[
    i===j, Re[mat[i,j]],
    i < j, mat[i,j],
    i > j, Conjugate[mat[j,i]]], {i,nr}, {j,nr}];

hTableT[mat_, nr_] := Transpose[hTable[mat, nr]];

headername =
    StringReplace[filename, "_dependence.cpp" ~~ EndOfString -> ".hpp"];

WriteString[filename,
  "#include <vector>\n",
  "#include ",InputForm[headername],"\n",
  "\n",
  "\n"
];

writeNBCDeps[filename, {
    M1i, M2i, M3i, A0i,
    _m2QiT, _m2UiT, _m2DiT, _m2LiT, _m2EiT, _m2NiT,
    _AuiT, _AdiT, _AniT, _AeiT,
    _YuiT, _YdiT, _YniT, _YeiT
},{
{"Fmssm_constraint_on_ms_n::dependence", Det[m2stop] - (scale0 Exp[t])^4},
{"Fmssm_constraint_on_gauge_couplings_n_::dependence", {g1 - g1i, g2 - g2i, g3 - g3i}},
{"Fmssm_constraint_on_yukawas_n_::dependence", {
    Yu - Table[YuiT[j,i], {i,3}, {j,3}],
    Yd - Table[YdiT[j,i], {i,3}, {j,3}],
    Ye - Table[YeiT[j,i], {i,3}, {j,3}]}},
{"Fmssm_constraint_on_ewsb_n_::dependence", {
    b vd/vu + (g2^2+gY^2) (vd^4-vu^4)/(4(vd^2+vu^2)) - (norm[mu]+m2Hu),
    b vu/vd + (g2^2+gY^2) (vu^4-vd^4)/(4(vd^2+vu^2)) - (norm[mu]+m2Hd),
    Im[mu], Im[b]}}
},{TAu, TAd, TAe}];

Close[filename];

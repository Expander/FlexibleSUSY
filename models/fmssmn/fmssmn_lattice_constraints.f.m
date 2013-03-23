<< models/fmssmn/fmssmn_lattice_defs.m

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

writeBCs[filename,
"g1i,g2i,g3i,\n" <>
"     &YuiT,YdiT,YniT,YeiT,\n" <>
"     &m2Hui,m2Hdi,m2QiT,m2UiT,m2DiT,m2LiT,m2NiT,m2EiT,\n" <>
"     &AuiT,AdiT,AniT,AeiT,\n" <>
"     &M1i,M2i,M3i,\n" <>
"     &vu,vd,\n",
"      double precision g1i,g2i,g3i\n" <>
"      double complex   YuiT(3,3),YdiT(3,3),YniT(3,3),YeiT(3,3)\n" <>
"      double precision m2Hui,m2Hdi\n" <>
"      double complex   m2QiT(3,3),m2UiT(3,3),m2DiT(3,3),\n" <>
"     &                 m2LiT(3,3),m2NiT(3,3),m2EiT(3,3)\n" <>
"      double complex   AuiT(3,3),AdiT(3,3),AniT(3,3),AeiT(3,3)\n" <>
"      double complex   M1i,M2i,M3i\n" <>
"      double precision vu,vd\n",
{
    YuiT, YdiT, YniT, YeiT,
    m2QiT, m2UiT, m2DiT, m2LiT, m2NiT, m2EiT,
    AuiT, AdiT, AniT, AeiT
},{
    M1i, M2i, M3i,
    _m2QiT, _m2UiT, _m2DiT, _m2LiT, _m2EiT, _m2NiT,
    _AuiT, _AdiT, _AniT, _AeiT,
    _YuiT, _YdiT, _YniT, _YeiT
},{
{"Fmssmn_mx", g1 - g2},
{"Fmssmn_higgs_masses", {m2Hu - m2Hui, m2Hd - m2Hdi}},
{"Fmssmn_gaugino_masses", {M1 - M1i, M2 - M2i, M3 - M3i}},
{"Fmssmn_sfermion_masses", {
    m2Q - hTableT[m2QiT, 3],
    m2U - hTableT[m2UiT, 3],
    m2D - hTableT[m2DiT, 3],
    m2L - hTableT[m2LiT, 3],
    m2E - hTableT[m2EiT, 3],
    m2N - hTableT[m2NiT, 3]}},
{"Fmssmn_trilinears", {
    TAu - Yu Table[AuiT[j,i], {i,3}, {j,3}],
    TAd - Yd Table[AdiT[j,i], {i,3}, {j,3}],
    TAn - Yn Table[AniT[j,i], {i,3}, {j,3}],
    TAe - Ye Table[AeiT[j,i], {i,3}, {j,3}]}},
{"Fmssmn_ms", Det[m2stop] - (scale0 Exp[t])^4},
{"Fmssmn_gauge_couplings", {g1 - g1i, g2 - g2i, g3 - g3i}},
{"Fmssmn_yude", {
    Yu - Table[YuiT[j,i], {i,3}, {j,3}],
    Yd - Table[YdiT[j,i], {i,3}, {j,3}],
    Ye - Table[YeiT[j,i], {i,3}, {j,3}]}},
{"Fmssmn_yn", Yn - Table[YniT[j,i], {i,3}, {j,3}]},
{"Fmssmn_ewsb", {
    b vd/vu + (g2^2+gY^2) (vd^4-vu^4)/(4(vd^2+vu^2)) - (norm[mu]+m2Hu),
    b vu/vd + (g2^2+gY^2) (vu^4-vd^4)/(4(vd^2+vu^2)) - (norm[mu]+m2Hd),
    Im[mu], Im[b]}}
}];

Close[filename];

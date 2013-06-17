<< models/fmssmn/fmssmn_lattice_defs.m
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
    M1i, M2i, M3i,
    _m2QiT, _m2UiT, _m2DiT, _m2LiT, _m2EiT, _m2NiT,
    _AuiT, _AdiT, _AniT, _AeiT,
    _YuiT, _YdiT, _YniT, _YeiT
},{
{"Fmssmn_constraint_on_yn_n_::dependence", Yn - Table[YniT[j,i], {i,3}, {j,3}]}
},
{Yn,m2N,TAn}];

Close[filename];

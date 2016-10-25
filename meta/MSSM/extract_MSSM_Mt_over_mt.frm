* Converts 2-loop relation between Mt and mt from
* [arXiv:hep-ph/0210258] to Mathematica form.
*
* Usage : form extract_MSSM_Mt_over_mt.frm
* Input : mt.res
* Output: mt.m

Symbol cs2t, cs2b, sn2t, sn4t, sn2b, sn4b, zt2;
Symbol mt, mmt, mb, mmb;
Symbol mgl, mmgl, mmst1, mmst2, mmsb1, mmsb2, mmsusy, mmu;
Function fin, den, Log;

#include mt.res
format mathematica;
#write <mt.m> "%e",resmt

.end

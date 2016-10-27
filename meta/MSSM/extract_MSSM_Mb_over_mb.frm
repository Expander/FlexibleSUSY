* Converts 2-loop relation between Mb and mb from
* [arXiv:hep-ph/0210258] to Mathematica form.
*
* Usage : form extract_MSSM_Mb_over_mb.frm
* Input : mb.res
* Output: mb.m

Symbol cs2t, cs2b, sn2t, sn4t, sn2b, sn4b, zt2;
Symbol mt, mmt, mb, mmb;
Symbol mgl, mmgl, mmst1, mmst2, mmsb1, mmsb2, mmsusy, mmu;
Function fin, den, Log;

#include mb.res
format mathematica;
#write <mb.m> "%e",resmb

.end

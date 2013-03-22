<< ../fmssm/fmssm_lattice_defs.m

lowert    = t  	 /. x -> w;
lowerg1   = g1 	 /. x -> w;
lowerg2   = g2 	 /. x -> w;
lowerg3   = g3 	 /. x -> w;
lowerYu   = Yu 	 /. x -> w;
lowerYd   = Yd 	 /. x -> w;
lowerYe   = Ye 	 /. x -> w;
lowermu   = mu 	 /. x -> w;
lowerb    = b  	 /. x -> w;
lowerM1   = M1 	 /. x -> w;
lowerM2   = M2 	 /. x -> w;
lowerM3   = M3 	 /. x -> w;
lowerm2Hu = m2Hu /. x -> w;
lowerm2Hd = m2Hd /. x -> w;
lowerm2Q  = m2Q  /. x -> w;
lowerm2U  = m2U  /. x -> w;
lowerm2D  = m2D  /. x -> w;
lowerm2L  = m2L  /. x -> w;
lowerm2E  = m2E  /. x -> w;
lowerTAu  = TAu  /. x -> w;
lowerTAd  = TAd  /. x -> w;
lowerTAe  = TAe  /. x -> w;

Clear[Yn,m2N,TAn];

<< fmssmn_lattice_defs.m

writeMCs[filename, "\n", "", {}, {}, {
{"Fmssm_fmssmn_gauge_couplings", {
    lowerg1   - g1,
    lowerg2   - g2,
    lowerg3   - g3 }},
{"Fmssm_fmssmn_Yukawas", {
    lowerYu   - Yu,
    lowerYd   - Yd,
    lowerYe   - Ye }},
{"Fmssm_fmssmn_mu_b", {
    lowermu   - mu,
    lowerb    - b }},
{"Fmssm_fmssmn_gaugino_masses", {
    lowerM1   - M1,
    lowerM2   - M2,
    lowerM3   - M3 }},
{"Fmssm_fmssmn_Higgs_masses", {
    lowerm2Hu - m2Hu,
    lowerm2Hd - m2Hd }},
{"Fmssm_fmssmn_sfermion_masses", {
    lowerm2Q  - m2Q,
    lowerm2U  - m2U,
    lowerm2D  - m2D,
    lowerm2L  - m2L,
    lowerm2E  - m2E }},
{"Fmssm_fmssmn_trilinears", {
    lowerTAu  - TAu,
    lowerTAd  - TAd,
    lowerTAe  - TAe }}
},
{Yn,m2N,TAn}];

Close[filename];

Needs["SUSYHD`"];
SetSMparameters[MtPole, AlphaS];
scale=Q;
tanBeta=TB;
Xtt=Xt;
gauginoM1=scale;
gauginoM2=scale;
gauginoM3=scale;
higgsMu=scale;
higgsAt=(higgsMu/tanBeta + Xtt scale);
mq3=scale;
mu3=scale;
md3=scale;
mq2=scale;
mu2=scale;
md2=scale;
mq1=scale;
mu1=scale;
md1=scale;
ml3=scale;
me3=scale;
ml2=scale;
me2=scale;
ml1=scale;
me1=scale;
mA=scale;
mass = MHiggs[{tanBeta, gauginoM1, gauginoM2, gauginoM3, higgsMu, higgsAt, mq3, mu3, md3, mq2, mu2, md2, mq1, mu1, md1, ml3, me3, ml2, me2, ml1, me1, mA}, Rscale->scale, scheme->"DRbar", hiOrd->{1,1,1,1}, numerical->True, split->False];
dmass= \[CapitalDelta]MHiggs[{tanBeta, gauginoM1, gauginoM2, gauginoM3, higgsMu, higgsAt, mq3, mu3, md3, mq2, mu2, md2, mq1, mu1, md1, ml3, me3, ml2, me2, ml1, me1, mA}, Rscale->scale, scheme->"DRbar", sources->{1,1,1}, numerical->False];

Write[Streams["stderr"], ToString[mass] <> " " <> ToString[dmass]];

Quit[0]

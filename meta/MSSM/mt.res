* Two-loop O(\alpha_s^2) MSSM corrections to the pole masses of heavy quarks
* A.Bednyakov, A.Onishchenko, V.Velizhanin, O.Veretin
* *
* This file contains the two-loop MSSM corrections to the relation
* between pole and running masses of t-quark (defined in Eq.(57))
* *
* zt2, zt3 - Zeta-functions
* mmsb1 - first  sbottom mass square
* mmsb2 - second sbottom mass square
* mmst1 - first  stop mass square
* mmst2 - second stop mass square
* mt - b-quark mass
* mmt - b-quark mass square
* mgl - gluino mass
* mmgl - gluino mass square
* csb, cs2b, cs4b - cosine mixing angle for sbottom sector
* snb, sn2b, sn4b - sine mixing angle for sbottom sector
* cst, cs2t, cs4t - cosine mixing angle for stop sector
* snt, sn2t, sn4t - sine mixing angle for stop sectort
* mmsusy - two first genration squarks mass square
* mmu - scale square
* *
* functions fin(mm1,mm2) defined as
* if mm1>mm2
*  fin(mm1,mm2)=(-7/2 - (7*mm2)/(2*mm1) +
*   PolyLog(mm2/mm1) - (mm2*PolyLog(mm2/mm1))/mm1 -
*   (3*mm2*Log(mm1))/mm1 - Log(mm1)*Log(mm1 - mm2) +
*   (mm2*Log(mm1)*Log(mm1 - mm2))/mm1 + (3*mm2*Log(mm2))/mm1 -
*   Log(mm1)*Log(mm2) + (2*mm2*Log(mm1)*Log(mm2))/mm1 +
*   Log(mm1 - mm2)*Log(mm2) - (mm2*Log(mm1 - mm2)*Log(mm2))/mm1 - PI^2/4 +
*   (mm2*PI^2)/(12*mm1) + Log(mm1)^2 - (3*mm2*Log(mm1)^2)/(2*mm1) -
*   (mm2*Log(mm2)^2)/(2*mm1));
* if mm1<mm2
*  fin(mm1,mm2)=(-7/2 - (7*mm2)/(2*mm1) -
*   PolyLog(mm1/mm2) + (mm2*PolyLog(mm1/mm2))/mm1 -
*   (3*mm2*Log(mm1))/mm1 + (3*mm2*Log(mm2))/mm1 +
*   (mm2*Log(mm1)*Log(mm2))/mm1 - Log(mm1)*Log(-mm1 + mm2) +
*   (mm2*Log(mm1)*Log(-mm1 + mm2))/mm1 + Log(mm2)*Log(-mm1 + mm2) -
*   (mm2*Log(mm2)*Log(-mm1 + mm2))/mm1 + PI^2/12 -
*   (mm2*PI^2)/(4*mm1) + Log(mm1)^2/2 - (mm2*Log(mm1)^2)/mm1 -
*   Log(mm2)^2/2);
* if mm1=mm2
*    fin(mm1,mm1)=-7-PI^2/6;
*
* where PolyLog(x)=PolyLog(x,2)=Spence(1-x) dilogarifm function

Local resmt =

       + cs2t^2 * (
          - 640/9
          - 128/9*zt2
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,1)*sn2t * (
          + 32/3*mmsusy*mt^-1*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,1) * (
          + 16/3*mmst1
          - 8*mmsusy
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,2)*sn2t * (
          + 32/3*mmst1*mmsusy*mt^-1*mgl
          - 32/3*mmst1^2*mt^-1*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,2) * (
          - 56/3*mmst1*mmsusy
          + 56/3*mmst1^2
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,3) * (
          - 32/3*mmst1^2*mmsusy
          + 32/3*mmst1^3
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,1)*sn2t * (
          - 32/3*mmsusy*mt^-1*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,1) * (
          + 16/3*mmst2
          - 8*mmsusy
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,2)*sn2t * (
          - 32/3*mmst2*mmsusy*mt^-1*mgl
          + 32/3*mmst2^2*mt^-1*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,2) * (
          - 56/3*mmst2*mmsusy
          + 56/3*mmst2^2
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,3) * (
          - 32/3*mmst2^2*mmsusy
          + 32/3*mmst2^3
          )

       + fin(mmgl,mmsusy) * (
          - 16/3
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,1)*sn2t * (
          + 4/3*mmsb1*mt^-1*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,1) * (
          - 1/3*mmsb1
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb1*mmst1*mt^-1*mgl
          + 4/3*mmsb1^2*mt^-1*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,2) * (
          + mmsb1*mmst1
          - mmsb1^2
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,3) * (
          + 4/3*mmsb1*mmst1^2
          - 4/3*mmsb1^2*mmst1
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,1)*sn2t * (
          - 4/3*mmsb1*mt^-1*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,1) * (
          - 1/3*mmsb1
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb1*mmst2*mt^-1*mgl
          - 4/3*mmsb1^2*mt^-1*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,2) * (
          + mmsb1*mmst2
          - mmsb1^2
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,3) * (
          + 4/3*mmsb1*mmst2^2
          - 4/3*mmsb1^2*mmst2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,1)*sn2t * (
          + 4/3*mmsb2*mt^-1*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,1) * (
          - 1/3*mmsb2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb2*mmst1*mt^-1*mgl
          + 4/3*mmsb2^2*mt^-1*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,2) * (
          + mmsb2*mmst1
          - mmsb2^2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,3) * (
          + 4/3*mmsb2*mmst1^2
          - 4/3*mmsb2^2*mmst1
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,1)*sn2t * (
          - 4/3*mmsb2*mt^-1*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,1) * (
          - 1/3*mmsb2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb2*mmst2*mt^-1*mgl
          - 4/3*mmsb2^2*mt^-1*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,2) * (
          + mmsb2*mmst2
          - mmsb2^2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,3) * (
          + 4/3*mmsb2*mmst2^2
          - 4/3*mmsb2^2*mmst2
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1)*sn2t * (
          + 88/9*mmst1*mt^-1*mgl
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1)*cs2t^2 * (
          + 5/3*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          - 154/9*mmst1^2
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 34/9*mmst1^2
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1) * (
          + 22/9*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,2) * (
          + 12*mmst1^2
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,3) * (
          + 16/3*mmst1^3
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1)*sn2t * (
          - 16/9*mmst1*mt^-1*mgl
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1)*cs2t^2 * (
          - 5/3*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          + 26/9*mmst1^2
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 34/9*mmst1^2
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1) * (
          + 16/9*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmst1*mmst2*mt^-1*mgl
          - 4/3*mmst1^2*mt^-1*mgl
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,2)*cs2t^2 * (
          + 26/9*mmst1*mmst2
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,2) * (
          - 17/9*mmst1*mmst2
          - mmst1^2
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,3) * (
          + 4/3*mmst1*mmst2^2
          - 4/3*mmst1^2*mmst2
          )

       + fin(mmst1,mmgl)*den(mmst1 - mmst2,1)*cs2t^2 * (
          - 128/9*mmst1
          )

       + fin(mmst1,mmsb1)*den(mmgl - mmst1,1) * (
          - 2/3*mmst1
          )

       + fin(mmst1,mmsb1)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb1*mmst1*mt^-1*mgl
          + 4/3*mmst1^2*mt^-1*mgl
          )

       + fin(mmst1,mmsb1)*den(mmgl - mmst1,2) * (
          + mmsb1*mmst1
          - 7/3*mmst1^2
          )

       + fin(mmst1,mmsb1)*den(mmgl - mmst1,3) * (
          + 4/3*mmsb1*mmst1^2
          - 4/3*mmst1^3
          )

       + fin(mmst1,mmsb2)*den(mmgl - mmst1,1) * (
          - 2/3*mmst1
          )

       + fin(mmst1,mmsb2)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb2*mmst1*mt^-1*mgl
          + 4/3*mmst1^2*mt^-1*mgl
          )

       + fin(mmst1,mmsb2)*den(mmgl - mmst1,2) * (
          + mmsb2*mmst1
          - 7/3*mmst1^2
          )

       + fin(mmst1,mmsb2)*den(mmgl - mmst1,3) * (
          + 4/3*mmsb2*mmst1^2
          - 4/3*mmst1^3
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,1)*sn2t * (
          - 4/9*mmst1*mt^-1*mgl
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,1)*cs2t^2 * (
          - 11/9*mmst1
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          - 8/9*mmst1^2
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,1) * (
          + mmst1
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmst1*mmst2*mt^-1*mgl
          + 4/3*mmst1^2*mt^-1*mgl
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,2)*cs2t^2 * (
          - 26/9*mmst1^2
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,2) * (
          + mmst1*mmst2
          + 5/9*mmst1^2
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,3) * (
          + 4/3*mmst1^2*mmst2
          - 4/3*mmst1^3
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,1)*sn2t * (
          + 4/9*mmst1*mt^-1*mgl
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,1)*cs2t^2 * (
          - 11/9*mmst1
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          + 8/9*mmst1^2
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,1) * (
          + 1/9*mmst1
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,2)*sn2t * (
          - 4/3*mmst1*mmst2*mt^-1*mgl
          + 4/3*mmst1^2*mt^-1*mgl
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,2)*cs2t^2 * (
          - 26/9*mmst1*mmst2
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,2) * (
          + 5/9*mmst1*mmst2
          + mmst1^2
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,3) * (
          - 4/3*mmst1*mmst2^2
          + 4/3*mmst1^2*mmst2
          )

       + fin(mmst1,mmsusy)*den(mmgl - mmst1,1) * (
          - 16/3*mmst1
          )

       + fin(mmst1,mmsusy)*den(mmgl - mmst1,2)*sn2t * (
          - 32/3*mmst1*mmsusy*mt^-1*mgl
          + 32/3*mmst1^2*mt^-1*mgl
          )

       + fin(mmst1,mmsusy)*den(mmgl - mmst1,2) * (
          + 8*mmst1*mmsusy
          - 56/3*mmst1^2
          )

       + fin(mmst1,mmsusy)*den(mmgl - mmst1,3) * (
          + 32/3*mmst1^2*mmsusy
          - 32/3*mmst1^3
          )

       + fin(mmst2,mmgl)*cs2t^2 * (
          - 128/9
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1)*sn2t * (
          + 16/9*mmst2*mt^-1*mgl
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1)*cs2t^2 * (
          + 26/9*mmst1
          + 11/9*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          - 26/9*mmst1^2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 34/9*mmst1^2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1) * (
          - 34/9*mmst1
          - 2*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmst1*mmst2*mt^-1*mgl
          + 4/3*mmst2^2*mt^-1*mgl
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,2)*cs2t^2 * (
          + 26/9*mmst1*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,2) * (
          - 17/9*mmst1*mmst2
          - mmst2^2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,3) * (
          - 4/3*mmst1*mmst2^2
          + 4/3*mmst1^2*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1)*sn2t * (
          - 88/9*mmst2*mt^-1*mgl
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1)*cs2t^2 * (
          - 154/9*mmst1
          - 139/9*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          + 154/9*mmst1^2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 34/9*mmst1^2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1) * (
          + 34/9*mmst1
          + 56/9*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,2) * (
          + 12*mmst2^2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,3) * (
          + 16/3*mmst2^3
          )

       + fin(mmst2,mmgl)*den(mmst1 - mmst2,1)*cs2t^2 * (
          + 128/9*mmst1
          )

       + fin(mmst2,mmsb1)*den(mmgl - mmst2,1) * (
          - 2/3*mmst2
          )

       + fin(mmst2,mmsb1)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb1*mmst2*mt^-1*mgl
          - 4/3*mmst2^2*mt^-1*mgl
          )

       + fin(mmst2,mmsb1)*den(mmgl - mmst2,2) * (
          + mmsb1*mmst2
          - 7/3*mmst2^2
          )

       + fin(mmst2,mmsb1)*den(mmgl - mmst2,3) * (
          + 4/3*mmsb1*mmst2^2
          - 4/3*mmst2^3
          )

       + fin(mmst2,mmsb2)*den(mmgl - mmst2,1) * (
          - 2/3*mmst2
          )

       + fin(mmst2,mmsb2)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb2*mmst2*mt^-1*mgl
          - 4/3*mmst2^2*mt^-1*mgl
          )

       + fin(mmst2,mmsb2)*den(mmgl - mmst2,2) * (
          + mmsb2*mmst2
          - 7/3*mmst2^2
          )

       + fin(mmst2,mmsb2)*den(mmgl - mmst2,3) * (
          + 4/3*mmsb2*mmst2^2
          - 4/3*mmst2^3
          )

       + fin(mmst2,mmsusy)*den(mmgl - mmst2,1) * (
          - 16/3*mmst2
          )

       + fin(mmst2,mmsusy)*den(mmgl - mmst2,2)*sn2t * (
          + 32/3*mmst2*mmsusy*mt^-1*mgl
          - 32/3*mmst2^2*mt^-1*mgl
          )

       + fin(mmst2,mmsusy)*den(mmgl - mmst2,2) * (
          + 8*mmst2*mmsusy
          - 56/3*mmst2^2
          )

       + fin(mmst2,mmsusy)*den(mmgl - mmst2,3) * (
          + 32/3*mmst2^2*mmsusy
          - 32/3*mmst2^3
          )

       + Log(mmgl)*cs2t^2 * (
          + 128/3
          )

       + Log(mmgl)^2*cs2t^2 * (
          - 64/9
          )

       + Log(mmgl)^2*den(mmgl - mmst1,1)*sn2t * (
          + 2/3*mmsb1*mt^-1*mgl
          + 2/3*mmsb2*mt^-1*mgl
          + 10*mmst1*mt^-1*mgl
          + 2/3*mmst2*mt^-1*mgl
          + 16/3*mmsusy*mt^-1*mgl
          )

       + Log(mmgl)^2*den(mmgl - mmst1,1)*cs2t^2 * (
          + 53/3*mmst1
          )

       + Log(mmgl)^2*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/3*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)^2*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          - 223/9*mmst1^2
          )

       + Log(mmgl)^2*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 43/9*mmst1^2
          )

       + Log(mmgl)^2*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2)*sn2t * (
          - 16/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmgl)^2*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          + 32/9*mmst1^3
          )

       + Log(mmgl)^2*den(mmgl - mmst1,1)*den(mmst1 - mmst2,3) * (
          - 16/9*mmst1^4
          )

       + Log(mmgl)^2*den(mmgl - mmst1,1) * (
          - 1/2*mmsb1
          - 1/2*mmsb2
          - 73/6*mmst1
          - 1/2*mmst2
          + 4/3*mmsusy
          )

       + Log(mmgl)^2*den(mmgl - mmst1,2)*sn2t * (
          + 2/3*mmsb1*mmst1*mt^-1*mgl
          + 2/3*mmsb2*mmst1*mt^-1*mgl
          + 2/3*mmst1*mmst2*mt^-1*mgl
          - 16*mmst1*mmsusy*mt^-1*mgl
          - 122/9*mmst1^2*mt^-1*mgl
          + 32/3*mmsusy^2*mt^-1*mgl
          )

       + Log(mmgl)^2*den(mmgl - mmst1,2)*cs2t^2 * (
          + 53/6*mmst1^2
          )

       + Log(mmgl)^2*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmgl)^2*den(mmgl - mmst1,2)*den(mmst1 - mmst2,2) * (
          + 8/9*mmst1^4
          )

       + Log(mmgl)^2*den(mmgl - mmst1,2) * (
          - 7/6*mmsb1*mmst1
          - 7/6*mmsb2*mmst1
          - 7/6*mmst1*mmst2
          + 52/3*mmst1*mmsusy
          + 151/9*mmst1^2
          - 8*mmsusy^2
          )

       + Log(mmgl)^2*den(mmgl - mmst1,3)*sn2t * (
          - 8/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmgl)^2*den(mmgl - mmst1,3) * (
          - 2/3*mmsb1*mmst1^2
          - 2/3*mmsb2*mmst1^2
          - 32/3*mmst1*mmsusy^2
          - 2/3*mmst1^2*mmst2
          + 16*mmst1^2*mmsusy
          + 46/3*mmst1^3
          )

       + Log(mmgl)^2*den(mmgl - mmst1,4) * (
          + 4/9*mmst1^4
          )

       + Log(mmgl)^2*den(mmgl - mmst2,1)*sn2t * (
          - 2/3*mmsb1*mt^-1*mgl
          - 2/3*mmsb2*mt^-1*mgl
          + 2/9*mmst1*mt^-1*mgl
          - 98/9*mmst2*mt^-1*mgl
          - 16/3*mmsusy*mt^-1*mgl
          )

       + Log(mmgl)^2*den(mmgl - mmst2,1)*cs2t^2 * (
          - 223/9*mmst1
          - 64/9*mmst2
          )

       + Log(mmgl)^2*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t * (
          - 8/3*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)^2*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          + 223/9*mmst1^2
          )

       + Log(mmgl)^2*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 43/9*mmst1^2
          )

       + Log(mmgl)^2*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2)*sn2t * (
          + 16/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmgl)^2*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          - 32/9*mmst1^3
          )

       + Log(mmgl)^2*den(mmgl - mmst2,1)*den(mmst1 - mmst2,3) * (
          + 16/9*mmst1^4
          )

       + Log(mmgl)^2*den(mmgl - mmst2,1) * (
          - 1/2*mmsb1
          - 1/2*mmsb2
          + 109/18*mmst1
          - 101/18*mmst2
          + 4/3*mmsusy
          )

       + Log(mmgl)^2*den(mmgl - mmst2,2)*sn2t * (
          - 2/3*mmsb1*mmst2*mt^-1*mgl
          - 2/3*mmsb2*mmst2*mt^-1*mgl
          - 14/9*mmst1*mmst2*mt^-1*mgl
          - 8/9*mmst1^2*mt^-1*mgl
          + 16*mmst2*mmsusy*mt^-1*mgl
          + 38/3*mmst2^2*mt^-1*mgl
          - 32/3*mmsusy^2*mt^-1*mgl
          )

       + Log(mmgl)^2*den(mmgl - mmst2,2)*cs2t^2 * (
          + 53/6*mmst2^2
          )

       + Log(mmgl)^2*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmgl)^2*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          - 32/9*mmst1^3
          )

       + Log(mmgl)^2*den(mmgl - mmst2,2)*den(mmst1 - mmst2,2) * (
          + 8/9*mmst1^4
          )

       + Log(mmgl)^2*den(mmgl - mmst2,2) * (
          - 7/6*mmsb1*mmst2
          - 7/6*mmsb2*mmst2
          + 11/18*mmst1*mmst2
          + 8/3*mmst1^2
          + 52/3*mmst2*mmsusy
          + 53/3*mmst2^2
          - 8*mmsusy^2
          )

       + Log(mmgl)^2*den(mmgl - mmst2,3)*sn2t * (
          + 8/9*mmst2^3*mt^-1*mgl
          )

       + Log(mmgl)^2*den(mmgl - mmst2,3) * (
          - 2/3*mmsb1*mmst2^2
          - 2/3*mmsb2*mmst2^2
          - 2/3*mmst1*mmst2^2
          - 32/3*mmst2*mmsusy^2
          + 16*mmst2^2*mmsusy
          + 46/3*mmst2^3
          )

       + Log(mmgl)^2*den(mmgl - mmst2,4) * (
          + 4/9*mmst2^4
          )

       + Log(mmgl)^2 * (
          - 166/9
          )

       + Log(mmgl)*Log(mmsb1)*den(mmgl - mmst1,1)*sn2t * (
          - 4/3*mmsb1*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmsb1)*den(mmgl - mmst1,1) * (
          + mmsb1
          - 4/3*mmst1
          )

       + Log(mmgl)*Log(mmsb1)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb1*mmst1*mt^-1*mgl
          + 8/3*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmsb1)*den(mmgl - mmst1,2) * (
          + 7/3*mmsb1*mmst1
          - 14/3*mmst1^2
          )

       + Log(mmgl)*Log(mmsb1)*den(mmgl - mmst1,3) * (
          + 4/3*mmsb1*mmst1^2
          - 8/3*mmst1^3
          )

       + Log(mmgl)*Log(mmsb1)*den(mmgl - mmst2,1)*sn2t * (
          + 4/3*mmsb1*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmsb1)*den(mmgl - mmst2,1) * (
          + mmsb1
          - 4/3*mmst2
          )

       + Log(mmgl)*Log(mmsb1)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb1*mmst2*mt^-1*mgl
          - 8/3*mmst2^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmsb1)*den(mmgl - mmst2,2) * (
          + 7/3*mmsb1*mmst2
          - 14/3*mmst2^2
          )

       + Log(mmgl)*Log(mmsb1)*den(mmgl - mmst2,3) * (
          + 4/3*mmsb1*mmst2^2
          - 8/3*mmst2^3
          )

       + Log(mmgl)*Log(mmsb1) * (
          + 4/3
          )

       + Log(mmgl)*Log(mmsb2)*den(mmgl - mmst1,1)*sn2t * (
          - 4/3*mmsb2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmsb2)*den(mmgl - mmst1,1) * (
          + mmsb2
          - 4/3*mmst1
          )

       + Log(mmgl)*Log(mmsb2)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb2*mmst1*mt^-1*mgl
          + 8/3*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmsb2)*den(mmgl - mmst1,2) * (
          + 7/3*mmsb2*mmst1
          - 14/3*mmst1^2
          )

       + Log(mmgl)*Log(mmsb2)*den(mmgl - mmst1,3) * (
          + 4/3*mmsb2*mmst1^2
          - 8/3*mmst1^3
          )

       + Log(mmgl)*Log(mmsb2)*den(mmgl - mmst2,1)*sn2t * (
          + 4/3*mmsb2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmsb2)*den(mmgl - mmst2,1) * (
          + mmsb2
          - 4/3*mmst2
          )

       + Log(mmgl)*Log(mmsb2)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb2*mmst2*mt^-1*mgl
          - 8/3*mmst2^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmsb2)*den(mmgl - mmst2,2) * (
          + 7/3*mmsb2*mmst2
          - 14/3*mmst2^2
          )

       + Log(mmgl)*Log(mmsb2)*den(mmgl - mmst2,3) * (
          + 4/3*mmsb2*mmst2^2
          - 8/3*mmst2^3
          )

       + Log(mmgl)*Log(mmsb2) * (
          + 4/3
          )

       + Log(mmgl)*Log(mmst1)*sn2t * (
          - 208/9*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*cs2t^2 * (
          - 64/9
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,1)*sn2t * (
          - 28/3*mmst1*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,1)*cs2t^2 * (
          - 121/9*mmst1
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/3*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          + 16/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*cs2t^2
       * (
          + 287/9*mmst1^2
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          - 43/9*mmst1^2
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          + 16/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          - 32/9*mmst1^3
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,3) * (
          + 16/9*mmst1^4
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,1) * (
          - 34/3*mmst1
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,2)*sn2t * (
          + 172/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,2)*sn4t*cs2t * (
          - 8/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,2)*cs2t^2 * (
          - 11/9*mmst1^2
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,2)*den(mmst1 - mmst2,2) * (
          - 8/9*mmst1^4
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,2) * (
          - 130/3*mmst1^2
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,3)*sn2t * (
          + 16/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,3)*cs2t^2 * (
          + 16/9*mmst1^3
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,3) * (
          - 212/9*mmst1^3
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst1,4) * (
          - 8/9*mmst1^4
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,1)*sn2t * (
          + 20/9*mmst1*mt^-1*mgl
          + 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 8/9*mmst1*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,1)*cs2t^2 * (
          + 5/3*mmst1
          - 22/9*mmst2
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/3*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          - 16/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*cs2t^2
       * (
          - 31/9*mmst1^2
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          + 43/9*mmst1^2
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          - 16/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          + 32/9*mmst1^3
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,3) * (
          - 16/9*mmst1^4
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,1) * (
          - 34/9*mmst1
          + 2/9*mmst2
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,2)*sn2t * (
          + 4*mmst1*mmst2*mt^-1*mgl
          + 8/9*mmst1^2*mt^-1*mgl
          - 8/3*mmst2^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,2)*cs2t^2 * (
          - 32/9*mmst1*mmst2
          - 52/9*mmst2^2
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          + 32/9*mmst1^3
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,2) * (
          - 8/9*mmst1^4
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,2) * (
          + 37/9*mmst1*mmst2
          - 8/3*mmst1^2
          + 10/9*mmst2^2
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,3)*cs2t^2 * (
          - 16/9*mmst1*mmst2^2
          )

       + Log(mmgl)*Log(mmst1)*den(mmgl - mmst2,3) * (
          + 28/9*mmst1*mmst2^2
          - 8/3*mmst2^3
          )

       + Log(mmgl)*Log(mmst1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          + 256/9*mmgl
          + 256/9*mmst1
          )

       + Log(mmgl)*Log(mmst1) * (
          + 52/9
          )

       + Log(mmgl)*Log(mmst2)*sn2t * (
          + 208/9*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*cs2t^2 * (
          + 64/3
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,1)*sn2t * (
          - 28/9*mmst2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,1)*sn4t*cs2t * (
          + 16/9*mmst1*mt^-1*mgl
          + 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,1)*cs2t^2 * (
          - 53/9*mmst1
          - 16/9*mmst2
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/3*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          - 16/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*cs2t^2
       * (
          + 31/9*mmst1^2
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          - 43/9*mmst1^2
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          + 16/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          - 32/9*mmst1^3
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,3) * (
          + 16/9*mmst1^4
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,1) * (
          + 61/9*mmst1
          + 25/9*mmst2
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,2)*sn2t * (
          - 28/9*mmst1*mmst2*mt^-1*mgl
          + 32/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,2)*cs2t^2 * (
          - 32/9*mmst1*mmst2
          - 52/9*mmst1^2
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,2)*den(mmst1 - mmst2,2) * (
          - 8/9*mmst1^4
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,2) * (
          + 53/9*mmst1*mmst2
          + 2*mmst1^2
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,3)*cs2t^2 * (
          - 16/9*mmst1^2*mmst2
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst1,3) * (
          + 28/9*mmst1^2*mmst2
          - 8/3*mmst1^3
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,1)*sn2t * (
          - 8/9*mmst1*mt^-1*mgl
          + 92/9*mmst2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 16/9*mmst1*mt^-1*mgl
          - 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,1)*cs2t^2 * (
          + 287/9*mmst1
          + 166/9*mmst2
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/3*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          + 16/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*cs2t^2
       * (
          - 287/9*mmst1^2
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          + 43/9*mmst1^2
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          - 16/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          + 32/9*mmst1^3
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,3) * (
          - 16/9*mmst1^4
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,1) * (
          - 59/9*mmst1
          - 161/9*mmst2
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,2)*sn2t * (
          + 8/9*mmst1*mmst2*mt^-1*mgl
          + 8/9*mmst1^2*mt^-1*mgl
          - 164/9*mmst2^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,2)*sn4t*cs2t * (
          + 8/9*mmst2^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,2)*cs2t^2 * (
          - 11/9*mmst2^2
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          + 32/9*mmst1^3
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,2) * (
          - 8/9*mmst1^4
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,2) * (
          - 16/9*mmst1*mmst2
          - 8/3*mmst1^2
          - 398/9*mmst2^2
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,3)*sn2t * (
          - 16/9*mmst2^3*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,3)*cs2t^2 * (
          + 16/9*mmst2^3
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,3) * (
          - 212/9*mmst2^3
          )

       + Log(mmgl)*Log(mmst2)*den(mmgl - mmst2,4) * (
          - 8/9*mmst2^4
          )

       + Log(mmgl)*Log(mmst2)*den(mmst1 - mmst2,1)*cs2t^2 * (
          - 256/9*mmgl
          - 256/9*mmst1
          )

       + Log(mmgl)*Log(mmst2) * (
          + 52/9
          )

       + Log(mmgl)*Log(mmsusy)*den(mmgl - mmst1,1)*sn2t * (
          - 32/3*mmsusy*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmsusy)*den(mmgl - mmst1,1) * (
          - 8/3*mmsusy
          )

       + Log(mmgl)*Log(mmsusy)*den(mmgl - mmst1,2)*sn2t * (
          + 32*mmst1*mmsusy*mt^-1*mgl
          - 64/3*mmsusy^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmsusy)*den(mmgl - mmst1,2) * (
          - 104/3*mmst1*mmsusy
          + 16*mmsusy^2
          )

       + Log(mmgl)*Log(mmsusy)*den(mmgl - mmst1,3) * (
          + 64/3*mmst1*mmsusy^2
          - 32*mmst1^2*mmsusy
          )

       + Log(mmgl)*Log(mmsusy)*den(mmgl - mmst2,1)*sn2t * (
          + 32/3*mmsusy*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmsusy)*den(mmgl - mmst2,1) * (
          - 8/3*mmsusy
          )

       + Log(mmgl)*Log(mmsusy)*den(mmgl - mmst2,2)*sn2t * (
          - 32*mmst2*mmsusy*mt^-1*mgl
          + 64/3*mmsusy^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmsusy)*den(mmgl - mmst2,2) * (
          - 104/3*mmst2*mmsusy
          + 16*mmsusy^2
          )

       + Log(mmgl)*Log(mmsusy)*den(mmgl - mmst2,3) * (
          + 64/3*mmst2*mmsusy^2
          - 32*mmst2^2*mmsusy
          )

       + Log(mmgl)*Log(mmt)*den(mmgl - mmst1,1)*sn2t * (
          + 16/3*mmst1*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmt)*den(mmgl - mmst1,1) * (
          - 16/3*mmst1
          )

       + Log(mmgl)*Log(mmt)*den(mmgl - mmst1,2) * (
          - 8/3*mmst1^2
          )

       + Log(mmgl)*Log(mmt)*den(mmgl - mmst2,1)*sn2t * (
          - 16/3*mmst2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmt)*den(mmgl - mmst2,1) * (
          - 16/3*mmst2
          )

       + Log(mmgl)*Log(mmt)*den(mmgl - mmst2,2) * (
          - 8/3*mmst2^2
          )

       + Log(mmgl)*Log(mmt) * (
          - 40/3
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst1,1)*sn2t * (
          - 16*mmst1*mt^-1*mgl
          + 16/9*mmst2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1*mt^-1*mgl
          - 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst1,1)*cs2t^2 * (
          - 16*mmst1
          + 16/9*mmst2
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*cs2t^2
       * (
          + 128/9*mmst1^2
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst1,1) * (
          + 332/9*mmst1
          - 16/9*mmst2
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst1,2)*sn2t * (
          + 16/9*mmst1*mmst2*mt^-1*mgl
          - 8/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst1,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2*mt^-1*mgl
          + 8/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst1,2)*cs2t^2 * (
          + 32/9*mmst1*mmst2
          - 32/3*mmst1^2
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst1,2) * (
          - 32/9*mmst1*mmst2
          + 178/9*mmst1^2
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst1,3)*cs2t^2 * (
          + 16/9*mmst1^2*mmst2
          - 16/9*mmst1^3
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst1,3) * (
          - 16/9*mmst1^2*mmst2
          + 8/9*mmst1^3
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst2,1)*sn2t * (
          - 16/9*mmst1*mt^-1*mgl
          + 16*mmst2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 8/9*mmst1*mt^-1*mgl
          + 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst2,1)*cs2t^2 * (
          + 16*mmst1
          - 16/9*mmst2
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*cs2t^2
       * (
          - 128/9*mmst1^2
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst2,1) * (
          - 16/9*mmst1
          + 332/9*mmst2
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst2,2)*sn2t * (
          - 16/9*mmst1*mmst2*mt^-1*mgl
          + 8/9*mmst2^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst2,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2*mt^-1*mgl
          - 8/9*mmst2^2*mt^-1*mgl
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst2,2)*cs2t^2 * (
          + 32/9*mmst1*mmst2
          - 32/3*mmst2^2
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst2,2) * (
          - 32/9*mmst1*mmst2
          + 178/9*mmst2^2
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst2,3)*cs2t^2 * (
          + 16/9*mmst1*mmst2^2
          - 16/9*mmst2^3
          )

       + Log(mmgl)*Log(mmu)*den(mmgl - mmst2,3) * (
          - 16/9*mmst1*mmst2^2
          + 8/9*mmst2^3
          )

       + Log(mmgl)*Log(mmu) * (
          + 36
          )

       + Log(mmgl)*den(mmgl - mmst1,1)*sn2t * (
          - 8/3*mmsb1*mt^-1*mgl
          - 8/3*mmsb2*mt^-1*mgl
          - 176/3*mmst1*mt^-1*mgl
          - 8/9*mmst2*mt^-1*mgl
          + 128/3*mmsusy*mt^-1*mgl
          )

       + Log(mmgl)*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1*mt^-1*mgl
          - 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmgl)*den(mmgl - mmst1,1)*cs2t^2 * (
          - 428/9*mmst1
          + 16/9*mmst2
          )

       + Log(mmgl)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          + 796/9*mmst1^2
          )

       + Log(mmgl)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          - 76/3*mmst1^2
          )

       + Log(mmgl)*den(mmgl - mmst1,1) * (
          + 2*mmsb1
          + 2*mmsb2
          + 290/9*mmst1
          + 2/9*mmst2
          - 16*mmsusy
          )

       + Log(mmgl)*den(mmgl - mmst1,2)*sn2t * (
          - 8/3*mmsb1*mmst1*mt^-1*mgl
          - 8/3*mmsb2*mmst1*mt^-1*mgl
          - 8/9*mmst1*mmst2*mt^-1*mgl
          - 64/3*mmst1*mmsusy*mt^-1*mgl
          + 56/9*mmst1^2*mt^-1*mgl
          + 32*mmsusy^2*mt^-1*mgl
          )

       + Log(mmgl)*den(mmgl - mmst1,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2*mt^-1*mgl
          + 8/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmgl)*den(mmgl - mmst1,2)*cs2t^2 * (
          + 32/9*mmst1*mmst2
          - 238/9*mmst1^2
          )

       + Log(mmgl)*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1) * (
          - 8/9*mmst1^3
          )

       + Log(mmgl)*den(mmgl - mmst1,2) * (
          + 14/3*mmsb1*mmst1
          + 14/3*mmsb2*mmst1
          + 10/9*mmst1*mmst2
          + 16/3*mmst1*mmsusy
          - 290/9*mmst1^2
          - 24*mmsusy^2
          )

       + Log(mmgl)*den(mmgl - mmst1,3)*cs2t^2 * (
          + 16/9*mmst1^2*mmst2
          - 16/9*mmst1^3
          )

       + Log(mmgl)*den(mmgl - mmst1,3) * (
          + 8/3*mmsb1*mmst1^2
          + 8/3*mmsb2*mmst1^2
          - 32*mmst1*mmsusy^2
          + 8/9*mmst1^2*mmst2
          + 64/3*mmst1^2*mmsusy
          - 200/9*mmst1^3
          )

       + Log(mmgl)*den(mmgl - mmst2,1)*sn2t * (
          + 8/3*mmsb1*mt^-1*mgl
          + 8/3*mmsb2*mt^-1*mgl
          + 8/9*mmst1*mt^-1*mgl
          + 176/3*mmst2*mt^-1*mgl
          - 128/3*mmsusy*mt^-1*mgl
          )

       + Log(mmgl)*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 8/9*mmst1*mt^-1*mgl
          + 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmgl)*den(mmgl - mmst2,1)*cs2t^2 * (
          + 812/9*mmst1
          + 368/9*mmst2
          )

       + Log(mmgl)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          - 796/9*mmst1^2
          )

       + Log(mmgl)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          + 76/3*mmst1^2
          )

       + Log(mmgl)*den(mmgl - mmst2,1) * (
          + 2*mmsb1
          + 2*mmsb2
          - 226/9*mmst1
          + 62/9*mmst2
          - 16*mmsusy
          )

       + Log(mmgl)*den(mmgl - mmst2,2)*sn2t * (
          + 8/3*mmsb1*mmst2*mt^-1*mgl
          + 8/3*mmsb2*mmst2*mt^-1*mgl
          + 8/9*mmst1*mmst2*mt^-1*mgl
          + 64/3*mmst2*mmsusy*mt^-1*mgl
          - 56/9*mmst2^2*mt^-1*mgl
          - 32*mmsusy^2*mt^-1*mgl
          )

       + Log(mmgl)*den(mmgl - mmst2,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2*mt^-1*mgl
          - 8/9*mmst2^2*mt^-1*mgl
          )

       + Log(mmgl)*den(mmgl - mmst2,2)*cs2t^2 * (
          + 32/9*mmst1*mmst2
          - 238/9*mmst2^2
          )

       + Log(mmgl)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          + 8/9*mmst1^3
          )

       + Log(mmgl)*den(mmgl - mmst2,2) * (
          + 14/3*mmsb1*mmst2
          + 14/3*mmsb2*mmst2
          + 2/9*mmst1*mmst2
          - 8/9*mmst1^2
          + 16/3*mmst2*mmsusy
          - 298/9*mmst2^2
          - 24*mmsusy^2
          )

       + Log(mmgl)*den(mmgl - mmst2,3)*cs2t^2 * (
          + 16/9*mmst1*mmst2^2
          - 16/9*mmst2^3
          )

       + Log(mmgl)*den(mmgl - mmst2,3) * (
          + 8/3*mmsb1*mmst2^2
          + 8/3*mmsb2*mmst2^2
          + 8/9*mmst1*mmst2^2
          - 32*mmst2*mmsusy^2
          + 64/3*mmst2^2*mmsusy
          - 200/9*mmst2^3
          )

       + Log(mmgl) * (
          + 232/3
          )

       + Log(mmsb1)^2*den(mmgl - mmst1,1)*sn2t * (
          + 2/3*mmsb1*mt^-1*mgl
          )

       + Log(mmsb1)^2*den(mmgl - mmst1,1) * (
          - 7/6*mmsb1
          + 2/3*mmst1
          )

       + Log(mmsb1)^2*den(mmgl - mmst1,2)*sn2t * (
          + 8/3*mmsb1*mmst1*mt^-1*mgl
          - 4/3*mmsb1^2*mt^-1*mgl
          - 4/3*mmst1^2*mt^-1*mgl
          )

       + Log(mmsb1)^2*den(mmgl - mmst1,2) * (
          - 4*mmsb1*mmst1
          + mmsb1^2
          + 7/3*mmst1^2
          )

       + Log(mmsb1)^2*den(mmgl - mmst1,3) * (
          - 8/3*mmsb1*mmst1^2
          + 4/3*mmsb1^2*mmst1
          + 4/3*mmst1^3
          )

       + Log(mmsb1)^2*den(mmgl - mmst2,1)*sn2t * (
          - 2/3*mmsb1*mt^-1*mgl
          )

       + Log(mmsb1)^2*den(mmgl - mmst2,1) * (
          - 7/6*mmsb1
          + 2/3*mmst2
          )

       + Log(mmsb1)^2*den(mmgl - mmst2,2)*sn2t * (
          - 8/3*mmsb1*mmst2*mt^-1*mgl
          + 4/3*mmsb1^2*mt^-1*mgl
          + 4/3*mmst2^2*mt^-1*mgl
          )

       + Log(mmsb1)^2*den(mmgl - mmst2,2) * (
          - 4*mmsb1*mmst2
          + mmsb1^2
          + 7/3*mmst2^2
          )

       + Log(mmsb1)^2*den(mmgl - mmst2,3) * (
          - 8/3*mmsb1*mmst2^2
          + 4/3*mmsb1^2*mmst2
          + 4/3*mmst2^3
          )

       + Log(mmsb1)^2 * (
          - 1/3
          )

       + Log(mmsb1)*Log(mmst1)*den(mmgl - mmst1,1) * (
          + 4/3*mmsb1
          )

       + Log(mmsb1)*Log(mmst1)*den(mmgl - mmst1,2)*sn2t * (
          - 4*mmsb1*mmst1*mt^-1*mgl
          + 8/3*mmsb1^2*mt^-1*mgl
          )

       + Log(mmsb1)*Log(mmst1)*den(mmgl - mmst1,2) * (
          + 17/3*mmsb1*mmst1
          - 2*mmsb1^2
          )

       + Log(mmsb1)*Log(mmst1)*den(mmgl - mmst1,3) * (
          + 4*mmsb1*mmst1^2
          - 8/3*mmsb1^2*mmst1
          )

       + Log(mmsb1)*Log(mmst2)*den(mmgl - mmst2,1) * (
          + 4/3*mmsb1
          )

       + Log(mmsb1)*Log(mmst2)*den(mmgl - mmst2,2)*sn2t * (
          + 4*mmsb1*mmst2*mt^-1*mgl
          - 8/3*mmsb1^2*mt^-1*mgl
          )

       + Log(mmsb1)*Log(mmst2)*den(mmgl - mmst2,2) * (
          + 17/3*mmsb1*mmst2
          - 2*mmsb1^2
          )

       + Log(mmsb1)*Log(mmst2)*den(mmgl - mmst2,3) * (
          + 4*mmsb1*mmst2^2
          - 8/3*mmsb1^2*mmst2
          )

       + Log(mmsb1)*Log(mmt) * (
          - 2/3
          )

       + Log(mmsb1)*den(mmgl - mmst1,1)*sn2t * (
          + 8/3*mmsb1*mt^-1*mgl
          )

       + Log(mmsb1)*den(mmgl - mmst1,1) * (
          + 2*mmst1
          )

       + Log(mmsb1)*den(mmgl - mmst1,2)*sn2t * (
          + 4*mmsb1^2*mt^-1*mgl
          - 4*mmst1^2*mt^-1*mgl
          )

       + Log(mmsb1)*den(mmgl - mmst1,2) * (
          + 4/3*mmsb1*mmst1
          - 3*mmsb1^2
          + 7*mmst1^2
          )

       + Log(mmsb1)*den(mmgl - mmst1,3) * (
          - 4*mmsb1^2*mmst1
          + 4*mmst1^3
          )

       + Log(mmsb1)*den(mmgl - mmst2,1)*sn2t * (
          - 8/3*mmsb1*mt^-1*mgl
          )

       + Log(mmsb1)*den(mmgl - mmst2,1) * (
          + 2*mmst2
          )

       + Log(mmsb1)*den(mmgl - mmst2,2)*sn2t * (
          - 4*mmsb1^2*mt^-1*mgl
          + 4*mmst2^2*mt^-1*mgl
          )

       + Log(mmsb1)*den(mmgl - mmst2,2) * (
          + 4/3*mmsb1*mmst2
          - 3*mmsb1^2
          + 7*mmst2^2
          )

       + Log(mmsb1)*den(mmgl - mmst2,3) * (
          - 4*mmsb1^2*mmst2
          + 4*mmst2^3
          )

       + Log(mmsb1) * (
          + 1/9
          )

       + Log(mmsb2)^2*den(mmgl - mmst1,1)*sn2t * (
          + 2/3*mmsb2*mt^-1*mgl
          )

       + Log(mmsb2)^2*den(mmgl - mmst1,1) * (
          - 7/6*mmsb2
          + 2/3*mmst1
          )

       + Log(mmsb2)^2*den(mmgl - mmst1,2)*sn2t * (
          + 8/3*mmsb2*mmst1*mt^-1*mgl
          - 4/3*mmsb2^2*mt^-1*mgl
          - 4/3*mmst1^2*mt^-1*mgl
          )

       + Log(mmsb2)^2*den(mmgl - mmst1,2) * (
          - 4*mmsb2*mmst1
          + mmsb2^2
          + 7/3*mmst1^2
          )

       + Log(mmsb2)^2*den(mmgl - mmst1,3) * (
          - 8/3*mmsb2*mmst1^2
          + 4/3*mmsb2^2*mmst1
          + 4/3*mmst1^3
          )

       + Log(mmsb2)^2*den(mmgl - mmst2,1)*sn2t * (
          - 2/3*mmsb2*mt^-1*mgl
          )

       + Log(mmsb2)^2*den(mmgl - mmst2,1) * (
          - 7/6*mmsb2
          + 2/3*mmst2
          )

       + Log(mmsb2)^2*den(mmgl - mmst2,2)*sn2t * (
          - 8/3*mmsb2*mmst2*mt^-1*mgl
          + 4/3*mmsb2^2*mt^-1*mgl
          + 4/3*mmst2^2*mt^-1*mgl
          )

       + Log(mmsb2)^2*den(mmgl - mmst2,2) * (
          - 4*mmsb2*mmst2
          + mmsb2^2
          + 7/3*mmst2^2
          )

       + Log(mmsb2)^2*den(mmgl - mmst2,3) * (
          - 8/3*mmsb2*mmst2^2
          + 4/3*mmsb2^2*mmst2
          + 4/3*mmst2^3
          )

       + Log(mmsb2)^2 * (
          - 1/3
          )

       + Log(mmsb2)*Log(mmst1)*den(mmgl - mmst1,1) * (
          + 4/3*mmsb2
          )

       + Log(mmsb2)*Log(mmst1)*den(mmgl - mmst1,2)*sn2t * (
          - 4*mmsb2*mmst1*mt^-1*mgl
          + 8/3*mmsb2^2*mt^-1*mgl
          )

       + Log(mmsb2)*Log(mmst1)*den(mmgl - mmst1,2) * (
          + 17/3*mmsb2*mmst1
          - 2*mmsb2^2
          )

       + Log(mmsb2)*Log(mmst1)*den(mmgl - mmst1,3) * (
          + 4*mmsb2*mmst1^2
          - 8/3*mmsb2^2*mmst1
          )

       + Log(mmsb2)*Log(mmst2)*den(mmgl - mmst2,1) * (
          + 4/3*mmsb2
          )

       + Log(mmsb2)*Log(mmst2)*den(mmgl - mmst2,2)*sn2t * (
          + 4*mmsb2*mmst2*mt^-1*mgl
          - 8/3*mmsb2^2*mt^-1*mgl
          )

       + Log(mmsb2)*Log(mmst2)*den(mmgl - mmst2,2) * (
          + 17/3*mmsb2*mmst2
          - 2*mmsb2^2
          )

       + Log(mmsb2)*Log(mmst2)*den(mmgl - mmst2,3) * (
          + 4*mmsb2*mmst2^2
          - 8/3*mmsb2^2*mmst2
          )

       + Log(mmsb2)*Log(mmt) * (
          - 2/3
          )

       + Log(mmsb2)*den(mmgl - mmst1,1)*sn2t * (
          + 8/3*mmsb2*mt^-1*mgl
          )

       + Log(mmsb2)*den(mmgl - mmst1,1) * (
          + 2*mmst1
          )

       + Log(mmsb2)*den(mmgl - mmst1,2)*sn2t * (
          + 4*mmsb2^2*mt^-1*mgl
          - 4*mmst1^2*mt^-1*mgl
          )

       + Log(mmsb2)*den(mmgl - mmst1,2) * (
          + 4/3*mmsb2*mmst1
          - 3*mmsb2^2
          + 7*mmst1^2
          )

       + Log(mmsb2)*den(mmgl - mmst1,3) * (
          - 4*mmsb2^2*mmst1
          + 4*mmst1^3
          )

       + Log(mmsb2)*den(mmgl - mmst2,1)*sn2t * (
          - 8/3*mmsb2*mt^-1*mgl
          )

       + Log(mmsb2)*den(mmgl - mmst2,1) * (
          + 2*mmst2
          )

       + Log(mmsb2)*den(mmgl - mmst2,2)*sn2t * (
          - 4*mmsb2^2*mt^-1*mgl
          + 4*mmst2^2*mt^-1*mgl
          )

       + Log(mmsb2)*den(mmgl - mmst2,2) * (
          + 4/3*mmsb2*mmst2
          - 3*mmsb2^2
          + 7*mmst2^2
          )

       + Log(mmsb2)*den(mmgl - mmst2,3) * (
          - 4*mmsb2^2*mmst2
          + 4*mmst2^3
          )

       + Log(mmsb2) * (
          + 1/9
          )

       + Log(mmst1)*sn2t * (
          + 280/9*mt^-1*mgl
          )

       + Log(mmst1)*cs2t^2 * (
          + 64/9
          )

       + Log(mmst1)^2*sn2t * (
          + 8*mt^-1*mgl
          )

       + Log(mmst1)^2*den(mmgl - mmst1,1)*sn2t * (
          - 2/9*mmst1*mt^-1*mgl
          - 4/9*mmst2*mt^-1*mgl
          )

       + Log(mmst1)^2*den(mmgl - mmst1,1)*sn4t*cs2t * (
          + 8/9*mmst1*mt^-1*mgl
          )

       + Log(mmst1)^2*den(mmgl - mmst1,1)*cs2t^2 * (
          - 14/9*mmst1
          - 11/9*mmst2
          )

       + Log(mmst1)^2*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn4t*cs2t * (
          - 16/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmst1)^2*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          - 77/9*mmst1^2
          )

       + Log(mmst1)^2*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 13/9*mmst1^2
          )

       + Log(mmst1)^2*den(mmgl - mmst1,1) * (
          - 2/3*mmsb1
          - 2/3*mmsb2
          + 419/18*mmst1
          + mmst2
          - 16/3*mmsusy
          )

       + Log(mmst1)^2*den(mmgl - mmst1,2)*sn2t * (
          + 2*mmsb1*mmst1*mt^-1*mgl
          - 4/3*mmsb1^2*mt^-1*mgl
          + 2*mmsb2*mmst1*mt^-1*mgl
          - 4/3*mmsb2^2*mt^-1*mgl
          + 2*mmst1*mmst2*mt^-1*mgl
          + 16*mmst1*mmsusy*mt^-1*mgl
          - 86/9*mmst1^2*mt^-1*mgl
          - 4/3*mmst2^2*mt^-1*mgl
          - 32/3*mmsusy^2*mt^-1*mgl
          )

       + Log(mmst1)^2*den(mmgl - mmst1,2)*sn4t*cs2t * (
          + 8/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmst1)^2*den(mmgl - mmst1,2)*cs2t^2 * (
          - 26/9*mmst1*mmst2
          - 85/18*mmst1^2
          )

       + Log(mmst1)^2*den(mmgl - mmst1,2) * (
          - 17/6*mmsb1*mmst1
          + mmsb1^2
          - 17/6*mmsb2*mmst1
          + mmsb2^2
          + 1/18*mmst1*mmst2
          - 68/3*mmst1*mmsusy
          + 92/3*mmst1^2
          + mmst2^2
          + 8*mmsusy^2
          )

       + Log(mmst1)^2*den(mmgl - mmst1,3)*sn2t * (
          - 8/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmst1)^2*den(mmgl - mmst1,3)*cs2t^2 * (
          - 16/9*mmst1^3
          )

       + Log(mmst1)^2*den(mmgl - mmst1,3) * (
          - 2*mmsb1*mmst1^2
          + 4/3*mmsb1^2*mmst1
          - 2*mmsb2*mmst1^2
          + 4/3*mmsb2^2*mmst1
          + 4/3*mmst1*mmst2^2
          + 32/3*mmst1*mmsusy^2
          - 2*mmst1^2*mmst2
          - 16*mmst1^2*mmsusy
          + 110/9*mmst1^3
          )

       + Log(mmst1)^2*den(mmgl - mmst1,4) * (
          + 4/9*mmst1^4
          )

       + Log(mmst1)^2*den(mmgl - mmst2,1)*sn2t * (
          - 2/3*mmst1*mt^-1*mgl
          )

       + Log(mmst1)^2*den(mmgl - mmst2,1)*cs2t^2 * (
          - 13/9*mmst1
          )

       + Log(mmst1)^2*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          + 13/9*mmst1^2
          )

       + Log(mmst1)^2*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 13/9*mmst1^2
          )

       + Log(mmst1)^2*den(mmgl - mmst2,1) * (
          + 17/18*mmst1
          )

       + Log(mmst1)^2*den(mmgl - mmst2,2) * (
          - 2/3*mmst1*mmst2
          )

       + Log(mmst1)^2*den(mmst1 - mmst2,1)*cs2t^2 * (
          - 128/9*mmgl
          - 64/9*mmst1
          )

       + Log(mmst1)^2 * (
          + 41/9
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,1)*sn2t * (
          - 8/9*mmst1*mt^-1*mgl
          + 8/3*mmst2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 16/9*mmst1*mt^-1*mgl
          - 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,1)*cs2t^2 * (
          + 5/9*mmst1
          + 38/9*mmst2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/3*mmst1^2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          + 16/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*cs2t^2
       * (
          - 5/9*mmst1^2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 17/9*mmst1^2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          - 16/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          + 32/9*mmst1^3
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,3) * (
          - 16/9*mmst1^4
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,1) * (
          - 11/3*mmst1
          - 34/9*mmst2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,2)*sn2t * (
          - 20/9*mmst1*mmst2*mt^-1*mgl
          - 8/9*mmst1^2*mt^-1*mgl
          + 8/3*mmst2^2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,2)*cs2t^2 * (
          + 28/3*mmst1*mmst2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,2)*den(mmst1 - mmst2,2) * (
          + 8/9*mmst1^4
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,2) * (
          - 11/3*mmst1*mmst2
          - 8/9*mmst1^2
          - 2*mmst2^2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,3)*cs2t^2 * (
          + 16/9*mmst1^2*mmst2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst1,3) * (
          - 8/3*mmst1*mmst2^2
          + 20/9*mmst1^2*mmst2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,1)*sn2t * (
          - 8/9*mmst1*mt^-1*mgl
          - 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,1)*cs2t^2 * (
          + 11/9*mmst1
          + 22/9*mmst2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/3*mmst1^2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          + 16/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*cs2t^2
       * (
          + 5/9*mmst1^2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 17/9*mmst1^2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          + 16/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          - 32/9*mmst1^3
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,3) * (
          + 16/9*mmst1^4
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,1) * (
          + 17/9*mmst1
          - 2/9*mmst2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,2)*sn2t * (
          - 4*mmst1*mmst2*mt^-1*mgl
          - 8/9*mmst1^2*mt^-1*mgl
          + 8/3*mmst2^2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,2)*cs2t^2 * (
          + 32/9*mmst1*mmst2
          + 52/9*mmst2^2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/9*mmst1^3*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          - 32/9*mmst1^3
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,2) * (
          + 8/9*mmst1^4
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,2) * (
          - 25/9*mmst1*mmst2
          + 8/3*mmst1^2
          - 10/9*mmst2^2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,3)*cs2t^2 * (
          + 16/9*mmst1*mmst2^2
          )

       + Log(mmst1)*Log(mmst2)*den(mmgl - mmst2,3) * (
          - 28/9*mmst1*mmst2^2
          + 8/3*mmst2^3
          )

       + Log(mmst1)*Log(mmsusy)*den(mmgl - mmst1,1) * (
          + 32/3*mmsusy
          )

       + Log(mmst1)*Log(mmsusy)*den(mmgl - mmst1,2)*sn2t * (
          - 32*mmst1*mmsusy*mt^-1*mgl
          + 64/3*mmsusy^2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmsusy)*den(mmgl - mmst1,2) * (
          + 136/3*mmst1*mmsusy
          - 16*mmsusy^2
          )

       + Log(mmst1)*Log(mmsusy)*den(mmgl - mmst1,3) * (
          - 64/3*mmst1*mmsusy^2
          + 32*mmst1^2*mmsusy
          )

       + Log(mmst1)*Log(mmt)*den(mmgl - mmst1,1)*sn2t * (
          - 16/3*mmst1*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmt)*den(mmgl - mmst1,1) * (
          + 16/3*mmst1
          )

       + Log(mmst1)*Log(mmt)*den(mmgl - mmst1,2) * (
          + 8/3*mmst1^2
          )

       + Log(mmst1)*Log(mmt) * (
          - 2/3
          )

       + Log(mmst1)*Log(mmu)*sn2t * (
          + 64/9*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmu)*cs2t^2 * (
          + 64/9
          )

       + Log(mmst1)*Log(mmu)*den(mmgl - mmst1,1)*sn2t * (
          + 16*mmst1*mt^-1*mgl
          - 16/9*mmst2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmu)*den(mmgl - mmst1,1)*sn4t*cs2t * (
          + 8/9*mmst1*mt^-1*mgl
          + 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmu)*den(mmgl - mmst1,1)*cs2t^2 * (
          + 16*mmst1
          - 16/9*mmst2
          )

       + Log(mmst1)*Log(mmu)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*cs2t^2
       * (
          - 128/9*mmst1^2
          )

       + Log(mmst1)*Log(mmu)*den(mmgl - mmst1,1) * (
          - 332/9*mmst1
          + 16/9*mmst2
          )

       + Log(mmst1)*Log(mmu)*den(mmgl - mmst1,2)*sn2t * (
          - 16/9*mmst1*mmst2*mt^-1*mgl
          + 8/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmu)*den(mmgl - mmst1,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2*mt^-1*mgl
          - 8/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmst1)*Log(mmu)*den(mmgl - mmst1,2)*cs2t^2 * (
          - 32/9*mmst1*mmst2
          + 32/3*mmst1^2
          )

       + Log(mmst1)*Log(mmu)*den(mmgl - mmst1,2) * (
          + 32/9*mmst1*mmst2
          - 178/9*mmst1^2
          )

       + Log(mmst1)*Log(mmu)*den(mmgl - mmst1,3)*cs2t^2 * (
          - 16/9*mmst1^2*mmst2
          + 16/9*mmst1^3
          )

       + Log(mmst1)*Log(mmu)*den(mmgl - mmst1,3) * (
          + 16/9*mmst1^2*mmst2
          - 8/9*mmst1^3
          )

       + Log(mmst1)*Log(mmu)*den(mmst1 - mmst2,1)*cs2t^2 * (
          - 128/9*mmst1
          )

       + Log(mmst1)*Log(mmu) * (
          - 128/9
          )

       + Log(mmst1)*den(mmgl - mmst1,1)*sn2t * (
          + 172/3*mmst1*mt^-1*mgl
          - 28/9*mmst2*mt^-1*mgl
          )

       + Log(mmst1)*den(mmgl - mmst1,1)*sn4t*cs2t * (
          + 16/9*mmst1*mt^-1*mgl
          + 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmst1)*den(mmgl - mmst1,1)*cs2t^2 * (
          + 229/9*mmst1
          - 49/9*mmst2
          )

       + Log(mmst1)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t * (
          - 8/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmst1)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          - 718/9*mmst1^2
          )

       + Log(mmst1)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 34/3*mmst1^2
          )

       + Log(mmst1)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          - 8/9*mmst1^3
          )

       + Log(mmst1)*den(mmgl - mmst1,1) * (
          - 2*mmsb1
          - 2*mmsb2
          - 13/3*mmst1
          + 43/9*mmst2
          - 16*mmsusy
          )

       + Log(mmst1)*den(mmgl - mmst1,2)*sn2t * (
          + 8/3*mmsb1*mmst1*mt^-1*mgl
          - 4*mmsb1^2*mt^-1*mgl
          + 8/3*mmsb2*mmst1*mt^-1*mgl
          - 4*mmsb2^2*mt^-1*mgl
          + 8/9*mmst1*mmst2*mt^-1*mgl
          + 64/3*mmst1*mmsusy*mt^-1*mgl
          + 52/9*mmst1^2*mt^-1*mgl
          - 4*mmst2^2*mt^-1*mgl
          - 32*mmsusy^2*mt^-1*mgl
          )

       + Log(mmst1)*den(mmgl - mmst1,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2*mt^-1*mgl
          - 8/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmst1)*den(mmgl - mmst1,2)*cs2t^2 * (
          - 110/9*mmst1*mmst2
          + 16*mmst1^2
          )

       + Log(mmst1)*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1) * (
          + 8/9*mmst1^3
          )

       + Log(mmst1)*den(mmgl - mmst1,2) * (
          - 6*mmsb1*mmst1
          + 3*mmsb1^2
          - 6*mmsb2*mmst1
          + 3*mmsb2^2
          + 56/9*mmst1*mmst2
          - 48*mmst1*mmsusy
          + 187/9*mmst1^2
          + 3*mmst2^2
          + 24*mmsusy^2
          )

       + Log(mmst1)*den(mmgl - mmst1,3)*cs2t^2 * (
          - 16/9*mmst1^2*mmst2
          + 16/9*mmst1^3
          )

       + Log(mmst1)*den(mmgl - mmst1,3) * (
          - 8/3*mmsb1*mmst1^2
          + 4*mmsb1^2*mmst1
          - 8/3*mmsb2*mmst1^2
          + 4*mmsb2^2*mmst1
          + 4*mmst1*mmst2^2
          + 32*mmst1*mmsusy^2
          - 8/9*mmst1^2*mmst2
          - 64/3*mmst1^2*mmsusy
          + 92/9*mmst1^3
          )

       + Log(mmst1)*den(mmgl - mmst2,1)*sn2t * (
          - 16/3*mmst1*mt^-1*mgl
          )

       + Log(mmst1)*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 8/9*mmst1*mt^-1*mgl
          )

       + Log(mmst1)*den(mmgl - mmst2,1)*cs2t^2 * (
          - 6*mmst1
          )

       + Log(mmst1)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmst1)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          + 26/3*mmst1^2
          )

       + Log(mmst1)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 34/3*mmst1^2
          )

       + Log(mmst1)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          + 8/9*mmst1^3
          )

       + Log(mmst1)*den(mmgl - mmst2,1) * (
          + 52/9*mmst1
          )

       + Log(mmst1)*den(mmgl - mmst2,2)*cs2t^2 * (
          + 16/9*mmst1*mmst2
          )

       + Log(mmst1)*den(mmgl - mmst2,2) * (
          - 40/9*mmst1*mmst2
          )

       + Log(mmst1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          - 128/3*mmgl
          - 640/9*mmst1
          )

       + Log(mmst1) * (
          + 5/9
          )

       + Log(mmst2)*sn2t * (
          - 280/9*mt^-1*mgl
          )

       + Log(mmst2)*cs2t^2 * (
          - 64
          )

       + Log(mmst2)^2*sn2t * (
          - 8*mt^-1*mgl
          )

       + Log(mmst2)^2*cs2t^2 * (
          - 64/9
          )

       + Log(mmst2)^2*den(mmgl - mmst1,1)*sn2t * (
          + 4/9*mmst1*mt^-1*mgl
          + 2/9*mmst2*mt^-1*mgl
          )

       + Log(mmst2)^2*den(mmgl - mmst1,1)*cs2t^2 * (
          + 8/3*mmst1
          - 11/9*mmst2
          )

       + Log(mmst2)^2*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          - 13/9*mmst1^2
          )

       + Log(mmst2)^2*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 13/9*mmst1^2
          )

       + Log(mmst2)^2*den(mmgl - mmst1,1) * (
          - 14/9*mmst1
          + 1/2*mmst2
          )

       + Log(mmst2)^2*den(mmgl - mmst1,2)*sn2t * (
          + 8/3*mmst1*mmst2*mt^-1*mgl
          - 4/3*mmst1^2*mt^-1*mgl
          - 4/3*mmst2^2*mt^-1*mgl
          )

       + Log(mmst2)^2*den(mmgl - mmst1,2)*cs2t^2 * (
          - 26/9*mmst1*mmst2
          + 26/9*mmst1^2
          )

       + Log(mmst2)^2*den(mmgl - mmst1,2) * (
          - 10/9*mmst1*mmst2
          - 5/9*mmst1^2
          + mmst2^2
          )

       + Log(mmst2)^2*den(mmgl - mmst1,3) * (
          + 4/3*mmst1*mmst2^2
          - 8/3*mmst1^2*mmst2
          + 4/3*mmst1^3
          )

       + Log(mmst2)^2*den(mmgl - mmst2,1)*sn2t * (
          + 2/3*mmst2*mt^-1*mgl
          )

       + Log(mmst2)^2*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 16/9*mmst1*mt^-1*mgl
          + 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmst2)^2*den(mmgl - mmst2,1)*cs2t^2 * (
          - 77/9*mmst1
          - 34/3*mmst2
          )

       + Log(mmst2)^2*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn4t*cs2t * (
          - 16/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmst2)^2*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          + 77/9*mmst1^2
          )

       + Log(mmst2)^2*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 13/9*mmst1^2
          )

       + Log(mmst2)^2*den(mmgl - mmst2,1) * (
          - 2/3*mmsb1
          - 2/3*mmsb2
          + 13/9*mmst1
          + 149/6*mmst2
          - 16/3*mmsusy
          )

       + Log(mmst2)^2*den(mmgl - mmst2,2)*sn2t * (
          - 2*mmsb1*mmst2*mt^-1*mgl
          + 4/3*mmsb1^2*mt^-1*mgl
          - 2*mmsb2*mmst2*mt^-1*mgl
          + 4/3*mmsb2^2*mt^-1*mgl
          + 2/3*mmst1*mmst2*mt^-1*mgl
          - 16*mmst2*mmsusy*mt^-1*mgl
          + 74/9*mmst2^2*mt^-1*mgl
          + 32/3*mmsusy^2*mt^-1*mgl
          )

       + Log(mmst2)^2*den(mmgl - mmst2,2)*sn4t*cs2t * (
          - 8/9*mmst2^2*mt^-1*mgl
          )

       + Log(mmst2)^2*den(mmgl - mmst2,2)*cs2t^2 * (
          - 137/18*mmst2^2
          )

       + Log(mmst2)^2*den(mmgl - mmst2,2) * (
          - 17/6*mmsb1*mmst2
          + mmsb1^2
          - 17/6*mmsb2*mmst2
          + mmsb2^2
          + 1/2*mmst1*mmst2
          - 68/3*mmst2*mmsusy
          + 281/9*mmst2^2
          + 8*mmsusy^2
          )

       + Log(mmst2)^2*den(mmgl - mmst2,3)*sn2t * (
          + 8/9*mmst2^3*mt^-1*mgl
          )

       + Log(mmst2)^2*den(mmgl - mmst2,3)*cs2t^2 * (
          - 16/9*mmst2^3
          )

       + Log(mmst2)^2*den(mmgl - mmst2,3) * (
          - 2*mmsb1*mmst2^2
          + 4/3*mmsb1^2*mmst2
          - 2*mmsb2*mmst2^2
          + 4/3*mmsb2^2*mmst2
          + 2/3*mmst1*mmst2^2
          + 32/3*mmst2*mmsusy^2
          - 16*mmst2^2*mmsusy
          + 98/9*mmst2^3
          )

       + Log(mmst2)^2*den(mmgl - mmst2,4) * (
          + 4/9*mmst2^4
          )

       + Log(mmst2)^2*den(mmst1 - mmst2,1)*cs2t^2 * (
          + 128/9*mmgl
          + 64/9*mmst1
          )

       + Log(mmst2)^2 * (
          + 41/9
          )

       + Log(mmst2)*Log(mmsusy)*den(mmgl - mmst2,1) * (
          + 32/3*mmsusy
          )

       + Log(mmst2)*Log(mmsusy)*den(mmgl - mmst2,2)*sn2t * (
          + 32*mmst2*mmsusy*mt^-1*mgl
          - 64/3*mmsusy^2*mt^-1*mgl
          )

       + Log(mmst2)*Log(mmsusy)*den(mmgl - mmst2,2) * (
          + 136/3*mmst2*mmsusy
          - 16*mmsusy^2
          )

       + Log(mmst2)*Log(mmsusy)*den(mmgl - mmst2,3) * (
          - 64/3*mmst2*mmsusy^2
          + 32*mmst2^2*mmsusy
          )

       + Log(mmst2)*Log(mmt)*den(mmgl - mmst2,1)*sn2t * (
          + 16/3*mmst2*mt^-1*mgl
          )

       + Log(mmst2)*Log(mmt)*den(mmgl - mmst2,1) * (
          + 16/3*mmst2
          )

       + Log(mmst2)*Log(mmt)*den(mmgl - mmst2,2) * (
          + 8/3*mmst2^2
          )

       + Log(mmst2)*Log(mmt) * (
          - 2/3
          )

       + Log(mmst2)*Log(mmu)*sn2t * (
          - 64/9*mt^-1*mgl
          )

       + Log(mmst2)*Log(mmu)*cs2t^2 * (
          - 64/9
          )

       + Log(mmst2)*Log(mmu)*den(mmgl - mmst2,1)*sn2t * (
          + 16/9*mmst1*mt^-1*mgl
          - 16*mmst2*mt^-1*mgl
          )

       + Log(mmst2)*Log(mmu)*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1*mt^-1*mgl
          - 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmst2)*Log(mmu)*den(mmgl - mmst2,1)*cs2t^2 * (
          - 16*mmst1
          + 16/9*mmst2
          )

       + Log(mmst2)*Log(mmu)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*cs2t^2
       * (
          + 128/9*mmst1^2
          )

       + Log(mmst2)*Log(mmu)*den(mmgl - mmst2,1) * (
          + 16/9*mmst1
          - 332/9*mmst2
          )

       + Log(mmst2)*Log(mmu)*den(mmgl - mmst2,2)*sn2t * (
          + 16/9*mmst1*mmst2*mt^-1*mgl
          - 8/9*mmst2^2*mt^-1*mgl
          )

       + Log(mmst2)*Log(mmu)*den(mmgl - mmst2,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2*mt^-1*mgl
          + 8/9*mmst2^2*mt^-1*mgl
          )

       + Log(mmst2)*Log(mmu)*den(mmgl - mmst2,2)*cs2t^2 * (
          - 32/9*mmst1*mmst2
          + 32/3*mmst2^2
          )

       + Log(mmst2)*Log(mmu)*den(mmgl - mmst2,2) * (
          + 32/9*mmst1*mmst2
          - 178/9*mmst2^2
          )

       + Log(mmst2)*Log(mmu)*den(mmgl - mmst2,3)*cs2t^2 * (
          - 16/9*mmst1*mmst2^2
          + 16/9*mmst2^3
          )

       + Log(mmst2)*Log(mmu)*den(mmgl - mmst2,3) * (
          + 16/9*mmst1*mmst2^2
          - 8/9*mmst2^3
          )

       + Log(mmst2)*Log(mmu)*den(mmst1 - mmst2,1)*cs2t^2 * (
          + 128/9*mmst1
          )

       + Log(mmst2)*Log(mmu) * (
          - 128/9
          )

       + Log(mmst2)*den(mmgl - mmst1,1)*sn2t * (
          + 4/9*mmst1*mt^-1*mgl
          + 52/9*mmst2*mt^-1*mgl
          )

       + Log(mmst2)*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmst2)*den(mmgl - mmst1,1)*cs2t^2 * (
          + 37/3*mmst1
          + 19/3*mmst2
          )

       + Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          - 26/3*mmst1^2
          )

       + Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 14*mmst1^2
          )

       + Log(mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          + 8/9*mmst1^3
          )

       + Log(mmst2)*den(mmgl - mmst1,1) * (
          - 137/9*mmst1
          - 23/3*mmst2
          )

       + Log(mmst2)*den(mmgl - mmst1,2)*sn2t * (
          - 4*mmst1^2*mt^-1*mgl
          + 4*mmst2^2*mt^-1*mgl
          )

       + Log(mmst2)*den(mmgl - mmst1,2)*cs2t^2 * (
          + 94/9*mmst1*mmst2
          + 26/3*mmst1^2
          )

       + Log(mmst2)*den(mmgl - mmst1,2) * (
          - 82/9*mmst1*mmst2
          - 5/3*mmst1^2
          - 3*mmst2^2
          )

       + Log(mmst2)*den(mmgl - mmst1,3) * (
          - 4*mmst1*mmst2^2
          + 4*mmst1^3
          )

       + Log(mmst2)*den(mmgl - mmst2,1)*sn2t * (
          + 8/3*mmst1*mt^-1*mgl
          - 520/9*mmst2*mt^-1*mgl
          )

       + Log(mmst2)*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1*mt^-1*mgl
          - 16/9*mmst2*mt^-1*mgl
          )

       + Log(mmst2)*den(mmgl - mmst2,1)*cs2t^2 * (
          - 734/9*mmst1
          - 152/3*mmst2
          )

       + Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t * (
          - 8/9*mmst1^2*mt^-1*mgl
          )

       + Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          + 718/9*mmst1^2
          )

       + Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 14*mmst1^2
          )

       + Log(mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          - 8/9*mmst1^3
          )

       + Log(mmst2)*den(mmgl - mmst2,1) * (
          - 2*mmsb1
          - 2*mmsb2
          + 50/3*mmst1
          + 52/9*mmst2
          - 16*mmsusy
          )

       + Log(mmst2)*den(mmgl - mmst2,2)*sn2t * (
          - 8/3*mmsb1*mmst2*mt^-1*mgl
          + 4*mmsb1^2*mt^-1*mgl
          - 8/3*mmsb2*mmst2*mt^-1*mgl
          + 4*mmsb2^2*mt^-1*mgl
          - 8/9*mmst1*mmst2*mt^-1*mgl
          - 64/3*mmst2*mmsusy*mt^-1*mgl
          - 16/9*mmst2^2*mt^-1*mgl
          + 32*mmsusy^2*mt^-1*mgl
          )

       + Log(mmst2)*den(mmgl - mmst2,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2*mt^-1*mgl
          + 8/9*mmst2^2*mt^-1*mgl
          )

       + Log(mmst2)*den(mmgl - mmst2,2)*cs2t^2 * (
          - 32/9*mmst1*mmst2
          + 74/3*mmst2^2
          )

       + Log(mmst2)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          - 8/9*mmst1^3
          )

       + Log(mmst2)*den(mmgl - mmst2,2) * (
          - 6*mmsb1*mmst2
          + 3*mmsb1^2
          - 6*mmsb2*mmst2
          + 3*mmsb2^2
          + 22/9*mmst1*mmst2
          + 8/9*mmst1^2
          - 48*mmst2*mmsusy
          + 20*mmst2^2
          + 24*mmsusy^2
          )

       + Log(mmst2)*den(mmgl - mmst2,3)*cs2t^2 * (
          - 16/9*mmst1*mmst2^2
          + 16/9*mmst2^3
          )

       + Log(mmst2)*den(mmgl - mmst2,3) * (
          - 8/3*mmsb1*mmst2^2
          + 4*mmsb1^2*mmst2
          - 8/3*mmsb2*mmst2^2
          + 4*mmsb2^2*mmst2
          - 8/9*mmst1*mmst2^2
          + 32*mmst2*mmsusy^2
          - 64/3*mmst2^2*mmsusy
          + 128/9*mmst2^3
          )

       + Log(mmst2)*den(mmst1 - mmst2,1)*cs2t^2 * (
          + 128/3*mmgl
          + 640/9*mmst1
          )

       + Log(mmst2) * (
          + 5/9
          )

       + Log(mmsusy)^2*den(mmgl - mmst1,1)*sn2t * (
          + 16/3*mmsusy*mt^-1*mgl
          )

       + Log(mmsusy)^2*den(mmgl - mmst1,1) * (
          - 4*mmsusy
          )

       + Log(mmsusy)^2*den(mmgl - mmst1,2) * (
          - 16/3*mmst1*mmsusy
          )

       + Log(mmsusy)^2*den(mmgl - mmst2,1)*sn2t * (
          - 16/3*mmsusy*mt^-1*mgl
          )

       + Log(mmsusy)^2*den(mmgl - mmst2,1) * (
          - 4*mmsusy
          )

       + Log(mmsusy)^2*den(mmgl - mmst2,2) * (
          - 16/3*mmst2*mmsusy
          )

       + Log(mmsusy)^2 * (
          + 8/3
          )

       + Log(mmsusy)*Log(mmt) * (
          - 16/3
          )

       + Log(mmsusy)*den(mmgl - mmst1,1)*sn2t * (
          - 128/3*mmsusy*mt^-1*mgl
          )

       + Log(mmsusy)*den(mmgl - mmst1,1) * (
          + 32*mmsusy
          )

       + Log(mmsusy)*den(mmgl - mmst1,2) * (
          + 128/3*mmst1*mmsusy
          )

       + Log(mmsusy)*den(mmgl - mmst2,1)*sn2t * (
          + 128/3*mmsusy*mt^-1*mgl
          )

       + Log(mmsusy)*den(mmgl - mmst2,1) * (
          + 32*mmsusy
          )

       + Log(mmsusy)*den(mmgl - mmst2,2) * (
          + 128/3*mmst2*mmsusy
          )

       + Log(mmsusy) * (
          + 152/9
          )

       + Log(mmt)*Log(mmu) * (
          + 64/3
          )

       + Log(mmt)*den(mmgl - mmst1,1) * (
          + 8/3*mmst1
          )

       + Log(mmt)*den(mmgl - mmst2,1) * (
          + 8/3*mmst2
          )

       + Log(mmt) * (
          + 8
          )

       + Log(mmu)*cs2t^2 * (
          + 128/9
          )

       + Log(mmu)^2 * (
          - 130/9
          )

       + Log(mmu)*den(mmgl - mmst1,1)*sn2t * (
          + 8/9*mmst1*mt^-1*mgl
          - 16/9*mmst2*mt^-1*mgl
          )

       + Log(mmu)*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1*mt^-1*mgl
          + 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmu)*den(mmgl - mmst1,1)*cs2t^2 * (
          + 88/9*mmst1
          - 8/3*mmst2
          )

       + Log(mmu)*den(mmgl - mmst1,1) * (
          - 58/3*mmst1
          + 8/3*mmst2
          )

       + Log(mmu)*den(mmgl - mmst1,2)*cs2t^2 * (
          - 16/9*mmst1*mmst2
          + 16/9*mmst1^2
          )

       + Log(mmu)*den(mmgl - mmst1,2) * (
          + 16/9*mmst1*mmst2
          - 8/9*mmst1^2
          )

       + Log(mmu)*den(mmgl - mmst2,1)*sn2t * (
          + 16/9*mmst1*mt^-1*mgl
          - 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmu)*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1*mt^-1*mgl
          + 8/9*mmst2*mt^-1*mgl
          )

       + Log(mmu)*den(mmgl - mmst2,1)*cs2t^2 * (
          - 8/3*mmst1
          + 88/9*mmst2
          )

       + Log(mmu)*den(mmgl - mmst2,1) * (
          + 8/3*mmst1
          - 58/3*mmst2
          )

       + Log(mmu)*den(mmgl - mmst2,2)*cs2t^2 * (
          - 16/9*mmst1*mmst2
          + 16/9*mmst2^2
          )

       + Log(mmu)*den(mmgl - mmst2,2) * (
          + 16/9*mmst1*mmst2
          - 8/9*mmst2^2
          )

       + Log(mmu) * (
          - 932/9
          )

       + den(mmgl - mmst1,1)*sn2t * (
          + 4/3*mmsb1*mt^-1*mgl*zt2
          + 8*mmsb1*mt^-1*mgl
          + 4/3*mmsb2*mt^-1*mgl*zt2
          + 8*mmsb2*mt^-1*mgl
          + 88/9*mmst1*mt^-1*mgl*zt2
          + 676/9*mmst1*mt^-1*mgl
          + 4/3*mmst2*mt^-1*mgl*zt2
          + 56/9*mmst2*mt^-1*mgl
          + 32/3*mmsusy*mt^-1*mgl*zt2
          + 64*mmsusy*mt^-1*mgl
          )

       + den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1*mt^-1*mgl
          + 8/9*mmst2*mt^-1*mgl
          )

       + den(mmgl - mmst1,1)*cs2t^2 * (
          + 41/9*mmst1*zt2
          + 439/9*mmst1
          - 8/3*mmst2
          )

       + den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          - 20*mmst1^2*zt2
          - 140*mmst1^2
          )

       + den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 20/3*mmst1^2*zt2
          + 428/9*mmst1^2
          )

       + den(mmgl - mmst1,1) * (
          - mmsb1*zt2
          - 6*mmsb1
          - mmsb2*zt2
          - 6*mmsb2
          + 50/9*mmst1*zt2
          + 61/9*mmst1
          - mmst2*zt2
          - 10/3*mmst2
          - 8*mmsusy*zt2
          - 48*mmsusy
          )

       + den(mmgl - mmst1,2)*cs2t^2 * (
          - 16/9*mmst1*mmst2
          + 16/9*mmst1^2
          )

       + den(mmgl - mmst1,2) * (
          - 4/3*mmsb1*mmst1*zt2
          - 8*mmsb1*mmst1
          - 4/3*mmsb2*mmst1*zt2
          - 8*mmsb2*mmst1
          - 4/3*mmst1*mmst2*zt2
          - 56/9*mmst1*mmst2
          - 32/3*mmst1*mmsusy*zt2
          - 64*mmst1*mmsusy
          + 44/3*mmst1^2*zt2
          + 868/9*mmst1^2
          )

       + den(mmgl - mmst1,3) * (
          + 16/3*mmst1^3*zt2
          + 112/3*mmst1^3
          )

       + den(mmgl - mmst2,1)*sn2t * (
          - 4/3*mmsb1*mt^-1*mgl*zt2
          - 8*mmsb1*mt^-1*mgl
          - 4/3*mmsb2*mt^-1*mgl*zt2
          - 8*mmsb2*mt^-1*mgl
          - 4/3*mmst1*mt^-1*mgl*zt2
          - 56/9*mmst1*mt^-1*mgl
          - 88/9*mmst2*mt^-1*mgl*zt2
          - 676/9*mmst2*mt^-1*mgl
          - 32/3*mmsusy*mt^-1*mgl*zt2
          - 64*mmsusy*mt^-1*mgl
          )

       + den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1*mt^-1*mgl
          + 8/9*mmst2*mt^-1*mgl
          )

       + den(mmgl - mmst2,1)*cs2t^2 * (
          - 20*mmst1*zt2
          - 428/3*mmst1
          - 139/9*mmst2*zt2
          - 821/9*mmst2
          )

       + den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*cs2t^2 * (
          + 20*mmst1^2*zt2
          + 140*mmst1^2
          )

       + den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 20/3*mmst1^2*zt2
          - 428/9*mmst1^2
          )

       + den(mmgl - mmst2,1) * (
          - mmsb1*zt2
          - 6*mmsb1
          - mmsb2*zt2
          - 6*mmsb2
          + 17/3*mmst1*zt2
          + 398/9*mmst1
          + 110/9*mmst2*zt2
          + 163/3*mmst2
          - 8*mmsusy*zt2
          - 48*mmsusy
          )

       + den(mmgl - mmst2,2)*cs2t^2 * (
          - 16/9*mmst1*mmst2
          + 16/9*mmst2^2
          )

       + den(mmgl - mmst2,2) * (
          - 4/3*mmsb1*mmst2*zt2
          - 8*mmsb1*mmst2
          - 4/3*mmsb2*mmst2*zt2
          - 8*mmsb2*mmst2
          - 4/3*mmst1*mmst2*zt2
          - 56/9*mmst1*mmst2
          - 32/3*mmst2*mmsusy*zt2
          - 64*mmst2*mmsusy
          + 44/3*mmst2^2*zt2
          + 868/9*mmst2^2
          )

       + den(mmgl - mmst2,3) * (
          + 16/3*mmst2^3*zt2
          + 112/3*mmst2^3
          )

       - 77/2
          + 8/9*zt2
         ;



(* FORMULAE FOR THE TWO-LOOP O(at^2) CORRECTIONS TO THE HIGGS MASSES *)
(*                                                                   *)
(* WRITTEN BY P. SLAVICH                                             *)
(*                                                                   *)  
(* BASED ON A. BRIGNOLE, G. DEGRASSI, P. SLAVICH AND F. ZWIRNER,     *)
(*                                                                   *) 
(* HEP-PH/0112177.                                                   *)


ct2 = (1+c2t)/2;st2=(1-c2t)/2;

(* F1 *)

f1ab = (
(2*BL*mu2/t*(mu2-BL+t)/delta[BL,mu2,t]
+(BL+mu2-t)/t)*phi[BL,mu2,t]
+2*A0*cb^2*(A0-6*t)/(A0-4*t)/t*phi[A0,t,t]
-2*cb^2*Li2[1-A0/t]
-2-(2+sb^2)/3*Pi^2
+ Log[t/q]*(
(4*(BL-mu2-10*t)*t+A0*(mu2-BL+(6+4 sb^2)*t))/(A0-4*t)/t
+1/delta[BL,mu2,t]*((BL-mu2)^3/t
+(2*mu2^2+2*BL*mu2+5*BL*t+mu2*t-4*BL^2-2*t^2)))
+Log[A0/q]*(4*A0*cb^2/(A0-4*t))
+Log[BL/q]*(-BL/t+BL*(-BL+mu2+t)^2/t/delta[BL,mu2,t])
+Log[mu2/q]*(mu2/t+mu2*(t^2-(BL-mu2)^2)/t/delta[BL,mu2,t])
+Nc*(Log[t/q]^2-Log[T1/q]^2/2-Log[T2/q]^2/2)
+Log[BL/q]*Log[mu2/q]-Log[BL/q]*Log[t/q]-Log[mu2/q]*Log[t/q]
-3*Log[t/q]^2-2*Log[T1/q]^2-2*Log[T2/q]^2);

f1c = (
2*mu2^2*(mu2+t-T1)/T1/delta[T1,mu2,t]*phi[mu2,t,T1]
-A0^2*(1+c2t^2)*cb^2*Yt^2
/2/T2/delta[A0,T1,T2]*phi[A0,T1,T2]
-cb^2*A0/T1*((2*Sqrt[t]+s2t*Yt)/Sqrt[t]
+(2*Sqrt[t]+s2t*Yt)^2/2/(A0-4*T1))*phi[A0,T1,T1]
-cb^2/T1*phi[A0,BL,T1]*(
2*A0*BL*(ct2*t+s2t*Yt*Sqrt[t]+st2*Yt^2)/delta[A0,T1,BL]
+(A0+BL-T1)*((1+c2t)*Sqrt[t]+s2t*Yt)/2/Sqrt[t])
+sb^2*(1+c2t+s2t*Xt/Sqrt[t])*Li2[1-BL/T1]
+(1-c2t)*Li2[1-mu2/T1]
+s2t^2*(T1-T2)^2/4/T1/T2*Nc+Pi^2/3
-3*s2t*(sb^2*Xt+cb^2*Yt)*t/Sqrt[t]/T1
-(sb^2*Xt^2+cb^2*Yt^2)*(6-2*c2t)/4/T1
+(3-c2t)*mu2/T1-(3-c2t)*cb^2*A0/2/T1
-(1+c2t)*t/2/T1-(1-c2t)*BL/2/T1
-(1+c2t^2)*(T1^2+T2^2)/4/T1/T2+5/2-s2t^2/2
+Log[t/q]*(1+mu2/t-t/T1-T1/t
+(-(mu2-T1)^3/t + 4*mu2^2 + 5*mu2*t + 2*t^2
+mu2^2*t/T1-t^3/T1-2*mu2*T1-2*T1^2)/delta[T1,mu2,t])
+Log[mu2/q]*(mu2*((-2+c2t)*t+T1)/t/T1
-mu2*((T1-t)^3+2*mu2*(t-T1)*T1+mu2^2*(t+T1))
/t/T1/delta[mu2,t,T1])
+Log[BL/q]*((1-c2t)/2*BL/T1
-cb^2*BL*(A0-BL+T1)/T1/delta[A0,T1,BL]*
(ct2*t+s2t*Yt*Sqrt[t]+st2*Yt^2)
+sb^2*BL*(ct2*t+s2t*Xt*Sqrt[t]+st2*Xt^2)/(BL-T1)/T1)
+Log[A0/q]*((3-c2t)*A0*cb^2/2/T1
+A0*cb^2*(2*Sqrt[t]+s2t*Yt)^2/2/(A0-4*T1)/T1
+A0*(1+c2t^2)*cb^2*(A0*(T1+T2)-(T1-T2)^2)*Yt^2
/4/T1/T2/delta[A0,T1,T2]
+A0*cb^2*(A0-BL-T1)/T1/delta[A0,BL,T1]*
(ct2*t+s2t*Yt*Sqrt[t]+st2*Yt^2))
+Log[T2/q]*((1+c2t^2)*T2/4/T1-Nc*s2t^2*T2/4/T1
-cb^2*(1+c2t^2)*Yt^2/4/T2
+sb^2*(1+c2t^2)*Xt^2/4/T1
+cb^2*(1+c2t^2)*Yt^2/4/T1/T2/delta[A0,T1,T2]*
(A0^2*T1-A0*(2*T1^2+5*T1*T2+T2^2)+(T1-T2)^2*(T1+T2)))
+Log[T1/q]*(cb^2*(ct2*t+s2t*Yt*Sqrt[t]+st2*Yt^2)*
(1/T1-(A0+BL-T1)/delta[A0,T1,BL])
+cb^2*(1+c2t^2)*Yt^2/4/T1/T2/delta[A0,T1,T2]*
(A0^2*T2-A0*(T1^2+5*T1*T2+2*T2^2)+(T1-T2)^2*(T1+T2))
+1/delta[T1,mu2,t]*((mu2-T1)^2*T1/t-(mu2-t)^2*(mu2+t)/T1
+6*mu2^2+4*mu2*t+2*t^2-3*mu2*T1-2*T1^2)
+cb^2*(A0-8*T1)/2/T1/(A0-4*T1)*(2*Sqrt[t]+s2t*Yt)^2
+sb^2*(BL-2*T1)/T1/(BL-T1)*
(ct2*t+s2t*Xt*Sqrt[t]+st2*Xt^2)
-(1-c2t)*(mu2-2*T1)/2/T1-s2t^2*Nc*(T1-2*T2)/4/T2
+sb^2*c2t*(1-c2t)*Xt^2/2/T1
+sb^2*(1+c2t^2)*Xt^2/4/T2
+sb^2*s2t*(t-3*T1)/Sqrt[t]/T1*Xt
-3*cb^2*s2t*(t+T1)/Sqrt[t]/T1*Yt
-cb^2*(5-2*c2t-c2t^2)*Yt^2/4/T1
-(3+c2t-8*sb^2)*t/2/T1
+(3-c2t)/2/T1*mu2+(1+c2t^2)*T1/4/T2
-T1/t+(-14+c2t-c2t^2)/2)
+(Nc+1)*s2t^2/4*
(3*Log[T1/q]^2-2*Log[T1/q]*Log[T2/q]-Log[T2/q]^2)
+s2t*(6*sb^2*Xt+5*cb^2*Yt)/2/Sqrt[t]*Log[T1/q]^2
+(9+sb^2-cb^2*c2t)/2*Log[T1/q]^2
+Log[T1/q]*Log[T2/q]-2*Log[t/q]*Log[T1/q]
+cb^2*(1+c2t+s2t*Yt/Sqrt[t])/2*(Log[A0/q]*Log[T1/q]
+Log[BL/q]*Log[T1/q]-Log[A0/q]*Log[BL/q])
     );

f1c1 = f1c /. {T1->T2,T2->T1,s2t->-s2t,c2t->-c2t};
f1c1 = f1c1 /.  phi[A0,T2,T1] -> T1/T2 phi[A0,T1,T2];

F1 = f1ab + f1c + f1c1;


(* F2 *)

f2ab = (-(3+Nc)/2 (Log[T1/q]^2-Log[T2/q]^2));

f2c = (
4*mu2^2*t/T1/delta[mu2,t,T1]*phi[mu2,t,T1]
+(A0*c2t^2*Yt^2/(T1-T2)/T2 
+(1+c2t^2)/2*(T1/T2-1)*Yt^2*A0/delta[A0,T1,T2])*cb^2*phi[A0,T1,T2]
-(2*A0*c2t^2*Sqrt[t]*Yt/s2t/T1/(T1-T2)
+s2t*(A0*Yt/2/Sqrt[t]/T1+Yt*A0*(A0-4*T1)/2/Sqrt[t]/T1/(T1-T2))
+ A0*(2*Sqrt[t]+s2t*Yt)^2/2/T1/(A0-4*T1)
+ A0/T1+A0*c2t^2*Yt^2/T1/(T1-T2))*cb^2*phi[A0,T1,T1]
-(2*A0*BL/T1*(ct2*t+s2t*Sqrt[t]*Yt+st2*Yt^2)/delta[A0,BL,T1]
+c2t^2*(A0+BL-T1)*Yt*Sqrt[t]/s2t/T1/(T1-T2)
+s2t*((A0+BL-T1)*Yt/4/Sqrt[t]/T1+delta[A0,BL,T1]*Yt/2/Sqrt[t]/T1/(T1-T2))
+ct2*(A0+BL-T1)/2/T1+c2t*Yt^2*(A0+BL-T1)/2/T1/(T1-T2)
+c2t*delta[A0,BL,T1]/2/T1/(T1-T2)-c2t*t*(A0+BL-T1)/2/T1/(T1-T2))
*cb^2*phi[A0,BL,T1]
+(s2t*(BL-T1)*Xt/Sqrt[t]/(T1-T2)+s2t*Xt/2/Sqrt[t]+ct2
-c2t*(t+T1-BL)/(T1-T2)+c2t*Xt*(2*c2t*Sqrt[t]+s2t*Xt)/s2t/(T1-T2))
*sb^2*Li2[1-BL/T1]
+(1-c2t-2*c2t*(mu2-T1)/(T1-T2))*Li2[1-mu2/T1]
-21*s2t*T1*At/2/(T1-T2)/Sqrt[t]+c2t*T1/2/(T1-T2)
+3*s2t*At/4/Sqrt[t]-3*s2t*Sqrt[t]*At/T1
+c2t*(2*BL+2*A0*cb^2-4*mu2-2*t+T1+2*(sb^2*Xt^2+cb^2*Yt^2))/4/T1
+(-5+c2t^2)/4/T1*(sb^2*Xt^2+cb^2*Yt^2)
-(2*BL+6*A0*cb^2-12*mu2+2*t+(1+c2t^2-Nc*s2t^2)*T2)/4/T1
+((1+c2t^2-Nc*s2t^2)*T1+(1+c2t^2)*(sb^2*Xt^2+cb^2*Yt^2))/4/T2
-t/T1*Log[t/q]+(2*mu2*(t-T1)*T1+mu2^2*(t+T1)-(t-T1)^3)
/T1/delta[T1,mu2,t]*Log[t*T1/q^2]-Log[mu2/q]*(2-c2t)*mu2/T1
+mu2*((t-T1)^2-mu2^2)/T1/delta[T1,mu2,t]*Log[mu2*T1/q^2]
+Log[BL/q]*(st2*BL/T1+sb^2*BL/T1/(BL-T1)*(ct2*t+s2t*Sqrt[t]*Xt+st2*Xt^2))
-cb^2*BL*(A0-BL+T1)*(ct2*t+s2t*Sqrt[t]*Yt+st2*Yt^2)
/T1/delta[A0,T1,BL]*Log[BL*T1/q^2]
+Log[A0/q]*((3-c2t)*cb^2*A0/2/T1
+A0*cb^2*(2*Sqrt[t]+s2t*Yt)^2/2/(A0-4*T1)/T1)
-cb^2*A0*(1+c2t^2)*(A0-T1-T2)*(T1-T2)*Yt^2
/4/T1/T2/delta[A0,T1,T2]*Log[A0*T1/q^2]
+cb^2*A0*(A0-BL-T1)*(ct2*t+s2t*Sqrt[t]*Yt+st2*Yt^2)
/T1/delta[A0,T1,BL]*Log[A0*T1/q^2]
+Log[T2/q]*(c2t^2*(1+Nc)*(T1+T2)/(T1-T2)
+c2t^2*(sb^2*Xt^2+cb^2*Yt^2)/(T1-T2)
-(1+c2t^2)*sb^2*(T1-T2)*Xt^2/4/T1/T2+(1+c2t^2-s2t^2*Nc)*T2/4/T1
+(1+c2t^2)*(sb^2*Xt^2+cb^2*Yt^2)/4/T2)
-(1+c2t^2)*cb^2*Yt^2*(A0^2*T1+(T1-T2)^3-A0*(2*T1^2+3*T1*T2-T2^2))
/4/T1/T2/delta[A0,T1,T2]*Log[T2*T1/q^2]
+Log[T1/q]*(4*mu2*(mu2+t-T1)/delta[T1,mu2,t]
+2*(1+c2t^2)*cb^2*Yt^2*A0*(A0-T1-3T2)/4/T2/delta[A0,T1,T2]
-2*cb^2*(A0+BL-T1)*(ct2*t+s2t*Sqrt[t]*Yt+st2*Yt^2)/delta[A0,T1,BL]
+cb^2*(2*Sqrt[t]+s2t*Yt)^2/2/T1-2*cb^2*(2*Sqrt[t]+s2t*Yt)^2/(A0-4*T1)
+sb^2*(ct2*t+s2t*Sqrt[t]*Xt+st2*Xt^2)/T1
-sb^2*(ct2*t+s2t*Sqrt[t]*Xt+st2*Xt^2)/(BL-T1)
-st2*mu2/T1+1-c2t-6*c2t^2*Sqrt[t]*At/s2t/(T1-T2)
+9*s2t*T1*At/Sqrt[t]/(T1-T2)-c2t*(1+c2t)*(sb^2*Xt^2+cb^2*Yt^2)/(T1-T2)
-(1+c2t^2)*sb^2*(T1-T2)*Xt^2/4/T1/T2
-c2t*(BL+A0*cb^2-2*mu2-t+(1+3*c2t*(1+Nc))*T1-c2t*(1+Nc)*T2)/(T1-T2)
-3*s2t*At/2/Sqrt[t]+s2t*Sqrt[t]*(sb^2*Xt-3*cb^2*Yt)/T1
-(3-2*c2t+c2t^2)/4/T1*(sb^2*Xt^2+cb^2*Yt^2)+s2t^2/2/T1*(sb^2*Xt^2-cb^2*Yt^2)
-6+c2t+(Nc+1)*s2t^2/2+((3-c2t)*mu2-(c2t+8*cb^2-5)*t)/2/T1
-(2*c2t^2+(1-Nc)*s2t^2)*T1/4/T2)
+Log[T2/q]^2*(3*c2t^2*Sqrt[t]*At/2/s2t/(T1-T2)
+c2t*(1-2*c2t)*(sb^2Xt^2+cb^2Yt^2)/4/(T1-T2)+3/8
+c2t*(2*BL+2*A0*cb^2-4*mu2-2*t+(1-2*c2t*(1+Nc))*T1
+(1-6*c2t*(1+Nc))*T2)/8/(T1-T2))
+Log[T1/q]^2*(3*s2t/2/Sqrt[t]*At+s2t*cb^2*A0*Yt/(T1-T2)/Sqrt[t]
+s2t*(BL-6*T1)*At/2/(T1-T2)/Sqrt[t]
+9*c2t^2*Sqrt[t]/2/s2t*At/(T1-T2)
+3*c2t*(1+2*c2t)*(sb^2Xt^2+cb^2Yt^2)/4/(T1-T2)
+c2t*(6*BL+6*A0*cb^2-12*mu2-6*t+(7+26*c2t*(1+Nc))*T1
-(1+2*c2t*(1+Nc))*T2)/8/(T1-T2)+25/8-c2t/2+s2t^2*(1+Nc))
-(s2t*(2*A0+2*BL-T1-T2)*Yt/4/Sqrt[t]/(T1-T2)
+(1+c2t)/4+c2t^2*Sqrt[t]*Yt/s2t/(T1-T2)
+c2t*(A0+BL-t-T1+Yt^2)/2/(T1-T2))*cb^2*Log[A0/T1]*Log[BL/T1]
-cb^2*s2t*2*A0*Yt/(T1-T2)/Sqrt[t]*Log[A0/q]*Log[T1/q]
-s2t*BL*At/(T1-T2)/Sqrt[t]*Log[BL/q]*Log[T1/q]
-(c2t^2*(1+Nc)*(T1+T2)/(T1-T2)
+c2t^2*(sb^2Xt^2+cb^2Yt^2)/(T1-T2))*Log[T1/q]*Log[T2/q]
);


f2c1 = f2c /. {T1->T2,T2->T1,s2t->-s2t,c2t->-c2t};
f2c1 = f2c1 /.  phi[A0,T2,T1] -> T1/T2 phi[A0,T1,T2];

F2 = f2ab + f2c - f2c1;


(* F3 *)

f3ab = (2+Nc)/2(2-Log[T1/q]-Log[T2/q])(2-(T1+T2)/(T1-T2)Log[T1/T2]);

f3c = (
(2*mu2*t*(mu2+t-T1)/T1/delta[T1,mu2,t]
-(4*mu2*t+2*delta[T1,mu2,t])/T1/(T1-T2)
-(mu2+t-T1)/T1)*phi[mu2,t,T1]
+(A0*(1+c2t^2)*(A0-2*(T1+T2))*Yt^2/2/T2/delta[A0,T1,T2]
+4*A0*c2t^2*(A0-2*(T1+T2))*Yt^2/T2/(T1-T2)^2
-(1-3 c2t^2)*Yt^2/2/T2)*cb^2*phi[A0,T1,T2]
+(-A0*(2*Sqrt[t]+s2t*Yt)^2/2/T1/(A0-4*T1)
+A0*(2*Sqrt[t]+s2t*Yt)^2/2/T1/(T1-T2)
-2*A0*c2t^2*Yt^2*(2*A0-7*T1-T2)/T1/(T1-T2)^2
+2*A0*(1-3*c2t^2)*Sqrt[t]*(A0-4*T1)*Yt/s2t/T1/(T1-T2)^2
-4*A0*c2t^2*Sqrt[t]*Yt/s2t/T1/(T1-T2))*cb^2*phi[A0,T1,T1]
+(-2*A0*BL*(ct2*t+s2t*Sqrt[t]*Yt+st2*Yt^2)/T1/delta[A0,T1,BL]
+2*(1-3*c2t^2)*Sqrt[t]*Yt*delta[A0,T1,BL]/s2t/T1/(T1-T2)^2
-2*c2t^2*Sqrt[t]*Yt*(A0+BL-T1)/s2t/T1/(T1-T2)
+3*c2t*delta[A0,T1,BL]*(t-Yt^2)/T1/(T1-T2)^2
+(A0+BL-T1)*(c2t*(t-Yt^2)+(ct2*t+s2t*Sqrt[t]*Yt+st2*Yt^2))
/T1/(T1-T2))*cb^2*phi[A0,BL,T1]
-(1-3*c2t^2)*sb^2*Xt^2/(T1-T2)*Li2[1-T2/T1]
+(1-c2t+2*(1-3c2t)*(mu2-T1)/(T1-T2)
-6*c2t*(mu2-T1)^2/(T1-T2)^2)*Li2[1-mu2/T1]
+(-4*(1-3*c2t^2)*sb^2*Sqrt[t]*(BL-T1)*Xt/s2t/(T1-T2)^2
+4*c2t^2*sb^2*Sqrt[t]*Xt/s2t/(T1-T2)
-2*sb^2*(ct2*t+s2t*Sqrt[t]*Xt+st2*Xt^2)/(T1-T2)
-2*sb^2*c2t*(3*BL-2*T1-T2)(t-Xt^2)/(T1-T2)^2)*Li2[1-BL/T1]
+2*Sqrt[t]*At*(3*(7-18*c2t^2)*T1-9*c2t^2*T2)/s2t/(T1-T2)^2
-3*s2t*Sqrt[t]*At*(4*T1-T2)/T1/(T1-T2)+2*c2t*At^2/(T1-T2)
-21*c2t*T1*(sb^2*Xt^2+cb^2*Yt^2)/(T1-T2)^2
+5*c2t/2/(T1-T2)*(sb^2*Xt^2+cb^2*Yt^2)
-(6-2*c2t)/4/T1*(sb^2*Xt^2+cb^2*Yt^2)
+3*mu2/T1+c2t*mu2*(T1+T2)/T1/(T1-T2)+18*c2t*mu2*T1/(T1-T2)^2
+3*cb^2*A0*c2t/(T1-T2)-cb^2*A0*(3-c2t)/2/T1
-12*(cb^2*A0+BL)*T1*c2t/(T1-T2)^2+3*c2t*BL/(T1-T2)-(1-c2t)*BL/2/T1
+15*c2t*T1*t/(T1-T2)^2-3*c2t*t/2/(T1-T2)-(1+c2t)*t/2/T1
-s2t^2/2*(1+Nc)-c2t^2*(1+Nc)*(T1-T2)^2/4/T1/T2
-9*(1+Nc)*T1*c2t^2/(T1-T2)-3*T1*(2*T1+3*T2)/(T1-T2)^2*c2t
-(1-Nc)*(T1-T2)^2/4/T1/T2-(28+6*Nc)*T1/2/(T1-T2)
+Log[t/q]*(-t/T1)-Log[t*T1/q^2]*t*((t-T1)^2-mu2^2)/T1/delta[T1,mu2,t]
+Log[mu2/q]*(-6*mu2*(3*mu2-T2)*c2t/(T1-T2)^2-(2-c2t)*mu2/T1)
-Log[mu2*T1/q^2]*mu2*(T1^2+mu2^2-t^2-2*mu2*T1)/T1/delta[T1,mu2,t]
+Log[BL/q]*(12*BL*c2t*T1/(T1-T2)^2
+sb^2*BL*(ct2*t+s2t*Sqrt[t]*Xt+st2*Xt^2)/T1/(BL-T1)
+(1-c2t)*BL/2/T1-3*c2t*BL/(T1-T2))
-Log[BL*T1/q^2]*cb^2*BL*(A0-BL+T1)/T1/delta[A0,T1,BL]
*(ct2*t+s2t*Sqrt[t]*Yt+st2*Yt^2)
+Log[A0/q]*A0*cb^2*((2*Sqrt[t]+s2t*Yt)^2/2/(A0-4*T1)/T1
+12*c2t*T1/(T1-T2)^2+(3-c2t)/2/T1-3*c2t/(T1-T2))
-Log[A0*T1/q^2]*A0*cb^2*(Yt^2*(1+c2t^2)/4/T1/T2/delta[A0,T1,T2]
*((T1+T2)^2-A0*(T1+T2)+4*T1*T2)-(A0-BL-T1)/T1/delta[A0,T1,BL]
*(ct2*t+s2t*Sqrt[t]*Yt+st2*Yt^2))
+Log[T2/q]*(12*c2t^2*T2*Sqrt[t]*At/s2t/(T1-T2)^2
-6*s2t*Sqrt[t]*T2*At/(T1-T2)^2
-3*(1-c2t)*cb^2*A0*T2/(T1-T2)^2-(1-3*c2t)*T2*BL/(T1-T2)^2
+6*(1-c2t)*T2*mu2/(T1-T2)^2-(1+3*c2t)*T2*t/(T1-T2)^2
+(sb^2*Xt^2+cb^2*Yt^2)*(-(2-3*c2t+23*c2t^2)*T2/(T1-T2)^2
-(1+5*c2t^2)/2/(T1-T2)-(1+c2t^2)/4/T2)
+(1+c2t^2)*sb^2*Xt^2*(T1^2-4*T1*T2-T2^2)/4/T1/T2/(T1-T2)
-(1+Nc)*c2t^2/4/(T1-T2)^2*(9*T1^2+32*T1*T2+19*T2^2)
-(1+Nc)*(T1-T2)/4/T1*c2t^2+T2*(2*T1+T2)/(T1-T2)^2*c2t
-T2*((7+3*Nc)*T1+(1-Nc)*T2)/2/(T1-T2)^2
-((1-Nc)*T1-(3-Nc)*T2)/2/(T1-T2)+(1-Nc)*T2/4/T1)
+Log[T2*T1/q^2]*(1+c2t^2)*cb^2*Yt^2/4/T1/T2/delta[A0,T1,T2]
*(A0^2*T1+(T1-T2)^3-2*T2*(T1^2-T2^2)-A0*(2*T1^2+T1*T2+T2^2))
+Log[T1/q]*(2*(mu2+t-T1)^2/delta[T1,mu2,t]
-cb^2*(1+c2t^2)*Yt^2/2/T2/delta[A0,T1,T2]
*(A0^2-4*T2*(T1-T2)-A0*(T1+3*T2))
-2*cb^2*(A0+BL-T1)*(ct2*t+s2t*Sqrt[t]*Yt+st2*Yt^2)/delta[A0,T1,BL]
+cb^2*(2*Sqrt[t]+s2t*Yt)^2/2/T1-2*cb^2*(2*Sqrt[t]+s2t*Yt)^2/(A0-4*T1)
+sb^2*(ct2*t+s2t*Sqrt[t]*Xt+st2*Xt^2)/T1
-sb^2*(ct2*t+s2t*Sqrt[t]*Xt+st2*Xt^2)/(BL-T1)-(1-c2t)/2*mu2/T1+1-c2t
+12*Sqrt[t]*At*((-3+7*c2t^2)*T1+c2t^2*T2)/s2t/(T1-T2)^2
+6*Sqrt[t]*(2*T1-T2)*s2t*At/(T1-T2)^2+sb^2*s2t^2*Xt^2/T1
+s2t*Sqrt[t]/T1*(sb^2*Xt-3*cb^2*Yt)
+(1+c2t^2)*sb^2*Xt^2*(T1^2+4*T1*T2-T2^2)/4/T1/T2/(T1-T2)
+(sb^2*Xt^2+cb^2*Yt^2)*((2+15*c2t+23*c2t^2)*T1/(T1-T2)^2
+(5-6*c2t-5*c2t^2)/2/(T1-T2)-(5-2*c2t-c2t^2)/4/T1)
+3*A0*cb^2*(2*(1+c2t)*T1-(1-c2t)*T2)/(T1-T2)^2
+BL*(2*(1+3*c2t)*T1-(1-3*c2t)*T2)/(T1-T2)^2
+t*(2*(1-6*c2t)*T1-(1+3*c2t)*T2)/(T1-T2)^2
+(5-c2t-8*cb^2)*t/2/T1+(3-c2t)*mu2/2/T1
-2*(3+12*c2t)*mu2*T1/(T1-T2)^2-6*(1-c2t)*mu2/(T1-T2)
-9/2+s2t^2/2*(1+Nc)+3*c2t/2+c2t*T1*(4*T1+11*T2)/(T1-T2)^2
+(1+Nc)*(65*T1^2+4*T1*T2-9*T2^2)*c2t^2/4/(T1-T2)^2
+(1+Nc)*(T1-T2)*c2t^2/4/T2+(1-Nc)*T1/4/T2
+T1*((7+3*Nc)*T1+(1-Nc)*T2)/2/(T1-T2)^2
+(5*(5+Nc)*T1+(1-Nc)*T2)/2/(T1-T2))
+(2*Log[mu2/q]^2-Log[mu2/T1]^2)*3*c2t*mu2^2/(T1-T2)^2
+(2*Log[BL/q]^2-Log[BL/T1]^2)*
(2*BL*(1-3*c2t^2)*Sqrt[t]*At/s2t/(T1-T2)^2-BL/2/(T1-T2)
-3*BL*c2t*(T1+T2-2*t+2*sb^2*Xt^2+2*cb^2*Yt^2)/2/(T1-T2)^2)
+(2*Log[A0/q]^2-Log[A0/T1]^2)*cb^2*
(4*A0*(1-3*c2t^2)*Sqrt[t]*Yt/s2t/(T1-T2)^2-3*A0/2/(T1-T2)
-3*A0*c2t*(T1+T2-2*t+2*Yt^2)/2/(T1-T2)^2)
+Log[t/T1]*Log[mu2/T1]*(T1+T2-2*t-2*mu2)/(T1-T2)
+Log[T2/T1]*Log[A0/T1]*Yt^2*cb^2*
(8*A0*c2t^2+(1-3*c2t^2)*(T1-T2))/2/(T1-T2)^2
+Log[BL/T1]*Log[A0/T1]*cb^2*(2*Sqrt[t]*Yt/(T1-T2)^2/s2t
*((A0+BL)*(1-3*c2t^2)-(1-2*c2t^2)*T1+c2t^2*T2)
+s2t*Sqrt[t]*Yt/(T1-T2)+(t+Yt^2)/2/(T1-T2)
+3*c2t*(2*A0+2*BL-T1-T2)*(t-Yt^2)/2/(T1-T2)^2)
+Log[T1/q]*Log[T2/q]*(2*c2t^2*(1+Nc)*(T1+T2)^2/(T1-T2)^2
+((1+5*c2t^2)*T1-(1-11*c2t^2)*T2)/2/(T1-T2)^2*(sb^2*Xt^2+cb^2*Yt^2))
+Log[T1/q]^2*(9*c2t^2*Sqrt[t]*At/s2t/(T1-T2)
+6*Sqrt[t]*T1*(2-5*c2t^2)*At/s2t/(T1-T2)^2
-3*Sqrt[t]*s2t*(2*T1-T2)*At/(T1-T2)^2
-((5+5*c2t+13*c2t^2)*T1-(3-4*c2t-10*c2t^2)*T2)
/2/(T1-T2)^2*(sb^2*Xt^2+cb^2*Yt^2)
+mu2*((6+5*c2t)*T1-(3-4*c2t)*T2)/(T1-T2)^2
+9/4-c2t+5/4*s2t^2*(1+Nc)
+c2t*(t-BL-A0*cb^2)*(5*T1+4*T2)/2/(T1-T2)^2
-(t+BL+3*A0*cb^2)*(2*T1-T2)/2/(T1-T2)^2
+c2t/4-9*T1*T2*c2t/2/(T1-T2)^2
-c2t^2*(1+Nc)*(8*T1^2+16*T1*T2-T2^2)/2/(T1-T2)^2
-T1*((14+5*Nc)*T1-(10+4*Nc)*T2)/2/(T1-T2)^2)
+Log[T2/q]^2*(-3*c2t^2*Sqrt[t]*(T1+T2)*At/s2t/(T1-T2)^2
+3*s2t*Sqrt[t]*T2*At/(T1-T2)^2
+(2*T2-c2t(T1+2*T2)+c2t^2*(2*T1+5*T2))
/2/(T1-T2)^2*(sb^2*Xt^2+cb^2*Yt^2)
+T2*(3*cb^2*A0+BL-6*mu2+t)/2/(T1-T2)^2
-c2t*(cb^2*A0+BL-mu2-t)*(T1+2*T2)/2/(T1-T2)^2
+c2t*mu2*(T1+2*T2)/2/(T1-T2)^2-3/4+s2t^2/4*(1+Nc)
+c2t^2*((T1+T2)^2+3*T2^2)/2/(T1-T2)^2*(1+Nc)
-c2t*((T1+T2)^2+2*T1*T2)/4/(T1-T2)^2
+T2*(2*(1+Nc)*T1+(2-Nc)*T2)/2/(T1-T2)^2)
);

f3c1 = f3c /. {T1->T2,T2->T1,s2t->-s2t,c2t->-c2t};
f3c1 = f3c1 /.  phi[A0,T2,T1] -> T1/T2 phi[A0,T1,T2];

F3 = f3ab + f3c + f3c1;

(* FA *)

faab = (5+2*Nc)*(T1*(1-Log[T1/q]+Log[T1/q]^2/2)-T2(1-Log[T2/q]+Log[T2/q]^2/2));

fac1 = (
-(delta[mu2,t,T1]+2*t*mu2)/T1*phi[mu2,t,T1]
+cb^2*Yt^2*(A0*c2t^2*(A0-2*(T1+T2))/2/(T1-T2)/T2 
-s2t^2*(T1-T2)/4/T2)*phi[A0,T1,T2]
+cb^2*(-A0*Sqrt[t]*Yt*(A0-4*T1)/s2t/At/T1/(T1-T2)*(Xt-s2t^2*At)
-A0*c2t^2*Yt^2*(A0-4*T1)/2/T1/(T1-T2)
+A0*(2*Sqrt[t]+s2t*Yt)^2/4/T1)*phi[A0,T1,T1]
+cb^2(-Sqrt[t]*Yt*delta[A0,BL,T1]*(Xt-s2t^2*At)/s2t/At/(T1-T2)/T1
+c2t*delta[A0,BL,T1]*(t-Yt^2)/2/T1/(T1-T2)
+(A0+BL-T1)*(ct2*t+s2t*Sqrt[t]*Yt+st2*Yt^2)/2/T1)*phi[A0,BL,T1]
-s2t^2*sb^2*Xt^2/2*Li2[1-T2/T1]
+(mu2-T1)*(1-c2t*(1+(mu2-T1)/(T1-T2)))*Li2[1-mu2/T1]
+(2*sb^2*Sqrt[t]*Xt*(Xt-s2t^2*At)*(BL-T1)/s2t/At/(T1-T2)
-c2t*sb^2*(t-Xt^2)*(BL-T1)/(T1-T2)
-sb^2*(ct2*t+s2t*Sqrt[t]*Xt+st2*Xt^2))*Li2[1-BL/T1]
-21*Sqrt[t]*T1*(Xt-s2t^2*At)/s2t/(T1-T2)
-Nc*c2t^2*T1-9*s2t*Sqrt[t]*At/2+c2t/2*BL
-c2t*T1/(T1-T2)*(A0*cb^2+2*BL-2*t-3*mu2+2*sb^2*Xt^2+2*cb^2*Yt^2)
-T1/2-T1*(2*T1+5*T2)/2/(T1-T2)*c2t-T1*c2t^2-(3+Nc)*(T1-T2)
+mu2*c2t*T2/(T1-T2)*Log[mu2/q]
+BL*c2t*(3*T1+T2)/2/(T1-T2)*Log[BL/q]
+A0*cb^2*c2t*(3*T1+T2)/2/(T1-T2)*Log[A0/q]
+Log[T2/q]*(-3*s2t*Sqrt[t]*T2*At/(T1-T2)
+c2t*T2/2/(T1-T2)*(BL+cb^2*A0-2*mu2-t+T1+sb^2*Xt^2+cb^2*Yt^2)
+((s2t^2-2)*T1+3*(3*s2t^2-4)*T2)/4/(T1-T2)*(sb^2*Xt^2+cb^2*Yt^2)
+Nc*(T1+T2)/4/(T1-T2)*(s2t^2*(T1+2*T2)-4*T2)
+s2t^2/4/(T1-T2)*(T1^2+3*T1*T2+2*T2^2)
-(T1^2-3*T1*T2+9*T2^2)/2/(T1-T2)-(3+Nc)*T2
-T2*(BL+3*A0*cb^2-6*mu2+t)/2/(T1-T2))
+Log[T1/q]*(18*Sqrt[t]*T1*(Xt-At*s2t^2)/(T1-T2)/s2t
+3*s2t*Sqrt[t]*(2*T1-T2)/(T1-T2)*At
+((9+8*c2t+9*c2t^2)*T1+(-5+2*c2t+c2t^2)*T2)/
4/(T1-T2)*(sb^2*Xt^2+cb^2*Yt^2)
+Nc*(4*(3-2*s2t^2)*T1^2+(-4+s2t^2)*T1*T2+s2t^2*T2^2)/4/(T1-T2)
+c2t^2*(8*T1^2-T1*T2-T2^2)/4/(T1-T2)
+cb^2*A0*(2*(3+c2t)*T1-(3-c2t)*T2)/2/(T1-T2)
+BL*(2*(1+c2t)*T1-(1-c2t)*T2)/2/(T1-T2)
-mu2*(3*(2+c2t)*T1-(3-c2t)*T2)/(T1-T2)
+t*(2*(1-2*c2t)*T1-(1+c2t)*T2)/2/(T1-T2)
+c2t*T1*(T1+4*T2)/2/(T1-T2)+(3+Nc)*T1
+(20*T1^2-11*T1*T2-T2^2)/4/(T1-T2))
+T1*Log[mu2/q]*Log[t/q]+(mu2+t-T1)*Log[T1/q]*Log[t/q]
+(mu2+t-T1+c2t*mu2^2/(T1-T2))*Log[T1/q]*Log[mu2/q]
-c2t*BL*(T1+T2)/4/(T1-T2)*Log[BL/q]^2
+Log[BL/q]*Log[T1/q]*(cb^2*Sqrt[t]*(A0-BL-T1)*Xt*Yt/At/s2t/(T1-T2)
-2*BL*sb^2*Sqrt[t]*Xt^2/s2t/At/(T1-T2)
+2*BL*Sqrt[t]*s2t^2*At/s2t/(T1-T2)
-s2t*cb^2*Sqrt[t]*Yt*(A0+BL-T1)/(T1-T2)-cb^2*Sqrt[t]*Yt*s2t/2
-c2t*cb^2*(A0+BL-T1)*(t-Yt^2)/2/(T1-T2)
+c2t*BL*(t-T1-sb^2*Xt^2-cb^2*Yt^2)/(T1-T2)
-(1-c2t)/2*BL-cb^2*((1+c2t)*t+(1-c2t)*Yt^2)/4)
-c2t*A0*cb^2*(T1+T2)/4/(T1-T2)*Log[A0/q]^2
+cb^2*Log[A0/q]*Log[T1/q]*(
-Sqrt[t]*(3*A0-BL+T1)*Yt*(Xt-s2t^2*At)/s2t/At/(T1-T2)
-Sqrt[t]*Yt*s2t/2+c2t^2*Yt^2*(T1-T2-2*A0)/4/(T1-T2)
-(6*A0+t+2*Yt^2)/4-c2t*(2*BL-2*A0-T1-T2)*(t-Yt^2)/4/(T1-T2)
-c2t*A0*(T1+T2)/2/(T1-T2))
+cb^2*Log[A0/q]*Log[BL/q]*(Sqrt[t]*T1*Yt*(Xt-s2t^2*At)/s2t/At/(T1-T2)
+s2t*Sqrt[t]*Yt/2-c2t*(t-Yt^2)*(T1+T2)/2/(T1-T2))
+cb^2*(A0*c2t^2*Yt^2/2/(T1-T2)+s2t^2*Yt^2/4)*Log[A0/q]*Log[T2/q]
+Log[T2/q]^2*(3*s2t*Sqrt[t]*T2*At/2/(T1-T2)+
(2-c2t+c2t^2)*T2/4/(T1-T2)*(sb^2*Xt^2+cb^2*Yt^2)
+Nc*(1+c2t^2)*T2^2/4/(T1-T2)-s2t^2*T2^2/4/(T1-T2)
+BL*T2*(1-c2t)/4/(T1-T2)-(3-c2t)*mu2*T2/2/(T1-T2)
+(3-c2t)*A0*T2*cb^2/4/(T1-T2)+(1+c2t)*t*T2/4/(T1-T2)
+T2*(9*T2-(4+c2t)*T1)/4/(T1-T2)+(3+Nc)*T2/2)
+Log[T1/q]*Log[T2/q]*(
Nc*(4*T1*T2-s2t^2*(T1+T2)^2)/4/(T1-T2)+(1+c2t^2)/4*(T1-T2)
+cb^2*c2t^2*Yt^2*(T1+T2-A0)/2/(T1-T2)
+c2t^2*T2*(T1+sb^2*Xt^2)/(T1-T2))
+Log[T1/q]^2*(cb^2*Sqrt[t]*Yt*(A0-BL+T1)*(Xt-s2t^2*At)/s2t/At/(T1-T2)
+Sqrt[t]*(BL-6*T1)*(Xt-s2t^2*At)/s2t/(T1-T2)
+cb^2*s2t*Sqrt[t]*Yt/2-3*At*Sqrt[t]*(2*T1-T2)*s2t/2/(T1-T2)
+Nc*T1*(3*(-2+s2t^2)*T1+2*s2t^2*T2)/4/(T1-T2)
-sb^2*Xt^2/4/(T1-T2)*((5+2*c2t+3*c2t^2)*T1+(-3+c2t+2*c2t^2)*T2-2*BL*c2t)
-cb^2*Yt^2/4/(T1-T2)*((3+c2t+4*c2t^2)*T1-s2t^2*T2-2*A0*c2t^2)
-(3+c2t)*cb^2*A0*T1/4/(T1-T2)-T1*BL*(1+c2t)/4/(T1-T2)
-t*c2t*sb^2*BL/2/(T1-T2)
+mu2/2/(T1-T2)*(2*(2+c2t)*T1-(1-c2t)*T2)-c2t*mu2^2/2/(T1-T2)
+t/4/(T1-T2)*((-5+c2t-(1-c2t)*sb^2)*T1+(4+(1+c2t)*sb^2)*T2)
-T1*((4+3*c2t^2)*T1+c2t*(3+2*c2t)*T2)/4/(T1-T2)-(3+Nc)*T1/2)
);

FA = faab + fac1 - (fac1 /.{T1->T2,T2->T1,s2t->-s2t,c2t->-c2t});

(* extra stuff *)

DmuF2 = -Nc/4*(Log[T1/q]^2-Log[T2/q]^2);
DmuF3 = -Nc/4*(Log[T1/q]+Log[T2/q]-2)*(2-(T1+T2)/(T1-T2)*Log[T1/T2]);

DAtF2 = -(3+Nc)/2*(Log[T1/q]^2-Log[T2/q]^2);
DAtF3 = -(3+Nc)/2*(Log[T1/q]+Log[T2/q]-2)*(2-(T1+T2)/(T1-T2)*Log[T1/T2]);

DM12 = ht^2*Nc*(s2t/4*mu*Sqrt[t]*(Log[T1/q]^2-Log[T2/q]^2)
+s2t^2/8*mu*At*(Log[T1/q]+Log[T2/q]-2)*(2-(T1+T2)/(T1-T2)*Log[T1/T2]));

DM22 = ht^2*Nc*(t*(Log[T1/q]^2+Log[T2/q]^2-2*Log[t/q]^2)
+s2t*At*Sqrt[t]*(Log[T1/q]^2-Log[T2/q]^2)+
s2t^2/4*At^2*(Log[T1/q]+Log[T2/q]-2)*(2-(T1+T2)/(T1-T2)*Log[T1/T2]));


(* DRBAR matrix elements *)

k = ht^2 Nc/(16 Pi^2)^2

M11 = k * ht^2 * mu2 * s2t^2/2 * (F3 + 2*DmuF3);

M12 = k * (ht^2 * At * mu * s2t^2/2 * (F3 + DmuF3 + DAtF3) + 
       ht^2 * mu * Sqrt[t] * s2t * (F2 + DmuF2) + DM12);

M22 = k * (ht^2 * At^2 * s2t^2/2 * (F3 + 2*DAtF3) +
       2 * ht^2 * At * Sqrt[t] * s2t * (F2 + DAtF2) +
       2 * ht^2 * t * F1 + DM22); 

DMA = k * ht^2/sb/cb * mu * At /(T1-T2) * FA;

(* tools *)

(*

delta[x_,y_,z_] := x^2 + y^2 + z^2 - 2 (x y + x z + y z);

Li2[x_]:= PolyLog[2,x];

lam[x_,y_,z_] := Sqrt[(1- x/z - y/z)^2 - 4 x y/z^2];
xp[x_,y_,z_]  := 1/2 ( 1 + x/z -y/z - lam[x,y,z]);            
xm[x_,y_,z_]  := 1/2 ( 1 - x/z +y/z - lam[x,y,z]);            
phi2[x_,y_,z_] :=  (2*Log[xp[x,y,z]]*Log[xm[x,y,z]] -  Log[x/z]*Log[y/z] -
                   2*(Li2[xp[x,y,z]] +Li2[xm[x,y,z]]) +
                   Pi^2/3)/lam[x,y,z] ;
phi[x_,y_,z_] := Which[ x <= z && y <= z,     phi2[x,y,z], 
                        z <= x && y <= x, z/x phi2[z,y,x], 
                        z <= y && x <= y, z/y phi2[z,x,y]]

*)

(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

Get["utils/load-FlexibleSUSY.m"];
Start["SM"];
Needs["TestSuite`", "TestSuite.m"];

Print["testing FSAntifield[] ..."];

TestEquality[FSAntiField[hh], hh];
TestEquality[FSAntiField[Ah], Ah];
TestEquality[FSAntiField[Hp], conj[Hp]];
TestEquality[FSAntiField[VP], VP];
TestEquality[FSAntiField[VZ], VZ];
TestEquality[FSAntiField[VWp], conj[VWp]];
TestEquality[FSAntiField[gP], conj[gP]];
TestEquality[FSAntiField[gZ], conj[gZ]];
TestEquality[FSAntiField[gWp], conj[gWp]];
TestEquality[FSAntiField[gWpC], conj[gWpC]];
TestEquality[FSAntiField[Fv], bar[Fv]];
TestEquality[FSAntiField[Fe], bar[Fe]];
TestEquality[FSAntiField[Fu], bar[Fu]];
TestEquality[FSAntiField[Fd], bar[Fd]];

Print["testing GetSMParticles[] ..."];

smParticles = {VG, gG, Hp, Fv, Ah, hh, VP, VZ, gP, gZ, gWp, gWpC, Fd, Fu, Fe, VWp};

TestEquality[Sort[TreeMasses`GetSMParticles[]], Sort[smParticles]];
TestEquality[Sort[TreeMasses`GetParticles[]]  , Sort[smParticles]];

Print["testing IsParticle[] ..."];

TestEquality[TreeMasses`IsParticle /@ smParticles, Array[True&, Length[smParticles]]];

conjSMParticles = {VG, conj[gG], conj[Hp], bar[Fv], Ah, hh, VP, VZ,
                   conj[gP], conj[gZ], conj[gWp], conj[gWpC], bar[Fd],
                   bar[Fu], bar[Fe], conj[VWp]};

TestEquality[Sort[FSAntiField /@ smParticles], Sort[conjSMParticles]];

TestEquality[TreeMasses`IsParticle /@ FSAntiField /@ smParticles,
             Array[True&, Length[conjSMParticles]]];

Print["testing IsSMParticle[] ..."];

TestEquality[TreeMasses`IsSMParticle /@ TreeMasses`GetSMParticles[],
             Array[True&, Length[TreeMasses`GetSMParticles[]]]];

TestEquality[TreeMasses`IsSMParticle /@ FSAntiField /@ TreeMasses`GetSMParticles[],
             Array[True&, Length[TreeMasses`GetSMParticles[]]]];

Print["testing IsSMParticleElementwise[] ..."];

TestEquality[IsSMParticleElementwise[Fu], {True, True, True}];
TestEquality[IsSMParticleElementwise[Fd], {True, True, True}];
TestEquality[IsSMParticleElementwise[Fe], {True, True, True}];
TestEquality[IsSMParticleElementwise[Fv], {True, True, True}];
TestEquality[IsSMParticleElementwise[VP], {True}];
TestEquality[IsSMParticleElementwise[VZ], {True}];
TestEquality[IsSMParticleElementwise[VWp], {True}];
TestEquality[IsSMParticleElementwise[hh, IncludeHiggs -> True], {True}];
TestEquality[IsSMParticleElementwise[hh, IncludeHiggs -> False], {True}];
TestEquality[IsSMParticleElementwise[Ah], {True}];
TestEquality[IsSMParticleElementwise[Hp], {True}];

PrintTestSummary[];

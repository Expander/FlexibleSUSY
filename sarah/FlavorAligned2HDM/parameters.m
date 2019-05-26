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

ParameterDefinitions = { 

{g1,        { Description -> "Hypercharge-Coupling"}},
{g2,        { Description -> "Left-Coupling"}},
{g3,        { Description -> "Strong-Coupling"}},    
{AlphaS,    {Description -> "Alpha Strong"}},	
{e,         { Description -> "electric charge"}}, 

{Gf,        { Description -> "Fermi's constant"}},
{aEWinv,    { Description -> "inverse weak coupling constant at mZ"}},

{Yu,        { Description -> "Up-Yukawa-Coupling",
			 DependenceNum ->  None}}, 
             									
{Yd,        { Description -> "Down-Yukawa-Coupling",
			  DependenceNum ->  None}},
             									
{Ye,        { Description -> "Lepton-Yukawa-Coupling",
			  DependenceNum ->  None}}, 

(* {Yu,        { Description -> "Up-Yukawa-Coupling",
			 DependenceNum ->  Sqrt[2]/v2* {{Mass[Fu,1],0,0},
             									{0, Mass[Fu,2],0},
             									{0, 0, Mass[Fu,3]}}}}, 
             									
{Yd,        { Description -> "Down-Yukawa-Coupling",
			  DependenceNum ->  Sqrt[2]/v1* {{Mass[Fd,1],0,0},
             									{0, Mass[Fd,2],0},
             									{0, 0, Mass[Fd,3]}}}},
             									
{Ye,        { Description -> "Lepton-Yukawa-Coupling",
			  DependenceNum ->  Sqrt[2]/v1* {{Mass[Fe,1],0,0},
             									{0, Mass[Fe,2],0},
             									{0, 0, Mass[Fe,3]}}}}, *)


{epsU,       { LaTeX -> "\\epsilon_u",
                        OutputName -> epsU,
                        LesHouches -> epsU }},

{epsD,       { LaTeX -> "\\epsilon_d",
                        OutputName -> epsD,
                        LesHouches -> epsD }},

{epsE,       { LaTeX -> "\\epsilon_e",
                        OutputName -> epsE,
                        LesHouches -> epsE }},

(* {xiE,        { LaTeX -> "\\xi_e",
               OutputName -> xiE,
               LesHouches -> {HMIX,41}}},
{xiD,        { LaTeX -> "\\xi_d",
               OutputName -> xiD,
               LesHouches -> {HMIX,42}}},
{xiU,        { LaTeX -> "\\xi_u",
               OutputName -> xiU,
               LesHouches -> {HMIX,43}}}, *)

                                                                            
{Lambda1,    { LaTeX -> "\\lambda_1",
               OutputName -> Lam1,
               LesHouches -> {HMIX,31}}},
{Lambda2,    { LaTeX -> "\\lambda_2",
               OutputName -> Lam2,
               LesHouches -> {HMIX,32}}},
{Lambda3,    { LaTeX -> "\\lambda_3",
               OutputName -> Lam3,
               LesHouches -> {HMIX,33}}},
{Lambda4,    { LaTeX -> "\\lambda_4",
               OutputName -> Lam4,
               LesHouches -> {HMIX,34}}},
{Lambda5,    { LaTeX -> "\\lambda_5",
               OutputName -> Lam5,
               LesHouches -> {HMIX,35}}},

{Lambda6,    { LaTeX -> "\\lambda_6",
               OutputName -> Lam6,
               LesHouches -> {HMIX,36}}},

{Lambda7,    { LaTeX -> "\\lambda_7",
               OutputName -> Lam7,
               LesHouches -> {HMIX,37}}},


{M112,    {    LaTeX -> "m^2_1",
               OutputName -> M112,
               LesHouches -> {HMIX,20}}},


{M222,    {    LaTeX -> "m^2_2",
               OutputName -> M222,
               LesHouches -> {HMIX,21}}},

{M122,    {    LaTeX -> "m^2_{12}",
               OutputName -> M122,
               LesHouches -> {HMIX,22}}},


{v1,        { Description -> "Down-VEV", LaTeX -> "v_1"}}, 
{v2,        { Description -> "Up-VEV", LaTeX -> "v_2"}},       
{v,         { Description -> "EW-VEV", DependenceSPheno -> Sqrt[v1^2 + v2^2] }},
             
{\[Beta],   { Description -> "Pseudo Scalar mixing angle"  }},             
{TanBeta,   { Description -> "Tan Beta" }},              
{\[Alpha],  { Description -> "Scalar mixing angle" }},  

{ZH,        { Description->"Scalar-Mixing-Matrix"}},

{ZP,        { Description->"Charged-Mixing-Matrix"}},  

{eta,       { Real -> True,
	      LaTeX -> "\\eta",
              OutputName -> "eta",
              LesHouches -> {HMIX,500} }},
              

{ThetaW,    { Description -> "Weinberg-Angle"}}, 

{ZZ, {Description ->   "Photon-Z Mixing Matrix"}},
{ZW, {Description -> "W Mixing Matrix" }},


{Vu,        {Description ->"Left-Up-Mixing-Matrix"}},
{Vd,        {Description ->"Left-Down-Mixing-Matrix"}},
{Uu,        {Description ->"Right-Up-Mixing-Matrix"}},
{Ud,        {Description ->"Right-Down-Mixing-Matrix"}}, 
{Ve,        {Description ->"Left-Lepton-Mixing-Matrix"}},
{Ue,        {Description ->"Right-Lepton-Mixing-Matrix"}}

 }; 
 


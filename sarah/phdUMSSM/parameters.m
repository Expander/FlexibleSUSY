ParameterDefinitions = { 

{g1,        { Description -> "Hypercharge-Coupling"}},
{g2,        { Description -> "Left-Coupling"}},
{g3,        { Description -> "Strong-Coupling"}},    
{AlphaS,    {Description -> "Alpha Strong"}},	
{Gf,        { Description -> "Fermi's constant"}},
{aEWinv,    { Description -> "inverse weak coupling constant at mZ"}},

{e,         { Description -> "electric charge"}}, 

{Yu,        { Description -> "Up-Yukawa-Coupling"   }}, 
{Yd,        { Description -> "Down-Yukawa-Coupling"}},
{Ye,        { Description -> "Lepton-Yukawa-Coupling"}}, 

{T[Ye],     { Description -> "Trilinear-Lepton-Coupling"}},
{T[Yd],     { Description -> "Trilinear-Down-Coupling"}}, 
{T[Yu],     { Description -> "Trilinear-Up-Coupling"}}, 

{\[Mu],     { Description -> "Mu-parameter"}}, 
{B[\[Mu]],  { Description -> "Bmu-parameter"}},        

{mq2,       { Description -> "Softbreaking left Squark Mass"}},
{me2,       { Description -> "Softbreaking right Slepton Mass"}},
{ml2,       { Description -> "Softbreaking left Slepton Mass"}},
{mu2,       { Description -> "Softbreaking right Up-Squark Mass"}},
{md2,       { Description -> "Softbreaking right Down-Squark Mass"}},
{mHd2,      { Description -> "Softbreaking Down-Higgs Mass"}}, 
{mHu2,      { Description -> "Softbreaking Up-Higgs Mass"}}, 

{MassB,     { Description -> "Bino Mass parameter"}},
{MassWB,    { Description -> "Wino Mass parameter"}},
{MassG,     { Description -> "Gluino Mass parameter"}},
                           
{vd,        { Description -> "Down-VEV"}}, 
{vu,        { Description -> "Up-VEV"}},       
{v,         { Description -> "EW-VEV"}},
             
{\[Beta],   { Description -> "Pseudo Scalar mixing angle"  }},             
{TanBeta,   { Description -> "Tan Beta" }},              
{\[Alpha],  { Description -> "Scalar mixing angle" }},                          
             
    
{ZD,        { Description -> "Down-Squark-Mixing-Matrix" }},             
{ZU,        { Description -> "Up-Squark-Mixing-Matrix"}},             
{ZE,        { Description -> "Slepton-Mixing-Matrix"}}, 
{ZV,        { Description->"Sneutrino Mixing-Matrix"}},   
{ZH,        { Description->"Scalar-Mixing-Matrix", 
              Real -> False, 
              Dependence -> None,
              DependenceNum -> None,
              DependenceOptional -> None,
              LesHouches -> NMHMIX      }},
{ZA,        { Description->"Pseudo-Scalar-Mixing-Matrix", 
               Real ->False,
               Dependence -> None,
              DependenceNum -> None,
              DependenceOptional -> None,
               LesHouches -> NMAMIX     }},
{ZP,        { Description->"Charged-Mixing-Matrix"}},                      
                                          
 
{ZN,        { Description->"Neutralino Mixing-Matrix",
              LesHouches -> NMNMIX   }}, 

{UP,        { Description->"Chargino-plus Mixing-Matrix"}}, 
{UM,        { Description->"Chargino-minus Mixing-Matrix"}}, 

{ZEL,       { Description ->"Left-Lepton-Mixing-Matrix"}},
{ZER,       { Description ->"Right-Lepton-Mixing-Matrix" }},                          
{ZDL,       { Description ->"Left-Down-Mixing-Matrix"}},                       
{ZDR,       { Description ->"Right-Down-Mixing-Matrix"}},              
{ZUL,       { Description ->"Left-Up-Mixing-Matrix"}},                        
{ZUR,       { Description ->"Right-Up-Mixing-Matrix"}},           
              
{ThetaW,    { Description -> "Weinberg-Angle"}},              

                         
{PhaseGlu,  { Description -> "Gluino-Phase" }},   
                                                                      
{\[Kappa],   {Description -> "Singlet Self-Interaction"}},                               
{T[\[Kappa]],  { Description -> "Softbreaking Singlet Self-Interaction" }}, 
{\[Lambda],   { Description -> "Singlet-Higgs-Interaction"   }},                               
{T[\[Lambda]],  {Description -> "Softbreaking Singlet-Higgs-Interaction"}},        
             
{ms2,       { Description -> "Softbreaking Singlet Mass" }},
{vS,        { Description -> "Singlet-VEV"}},


(* ------------------------------------------- *)
(* Expression for ThetaWp has to be added!     *)
(* ------------------------------------------- *)

{ZZ, {Description ->   "Photon-Z-Z' Mixing Matrix"}},
{ZW, {Description -> "W Mixing Matrix" }},
{ZfW, {Description ->    "Wino Mixing Matrix"}},

{ThetaWp,  { Description -> "Theta'",
              DependenceNum -> 0  }},
              
              
{gp,       {Description -> "U(1)' Gauge Coupling"}},


{Qq,        { LaTeX -> "Q_q",
              Real -> True,
		OutputName -> Qq,
			  LesHouches -> {XCharge,1}}},             
{Ql,        { LaTeX -> "Q_l",
              Real -> True,
		OutputName -> Ql,
			  LesHouches -> {XCharge,2}}},
{QHu,        { LaTeX -> "Q_{H_u}",
              Real -> True,
		OutputName -> QHu,
			  LesHouches -> {XCharge,3}}},
{QHd,        { LaTeX -> "Q_{H_d}",
              Real -> True,
		OutputName -> QHd,
			  LesHouches -> {XCharge,4}}},
{Qd,        { LaTeX -> "Q_d",
              Real -> True,
		OutputName -> Qd,
			  LesHouches -> {XCharge,5}}},
{Qu,        { LaTeX -> "Q_u",
              Real -> True,
		OutputName -> Qu,
			  LesHouches -> {XCharge,6}}},             
{Qe,        { LaTeX -> "Q_e",
              Real -> True,
		OutputName -> Qe,
			  LesHouches -> {XCharge,7}}},
{Qs,        { LaTeX -> "Q_s",
              Real -> True,
		OutputName -> Qs,
			  LesHouches -> {XCharge,8}}},

{ThetaZ,    { LaTeX -> "\\Theta_Z",
			  OutputName -> ThZ}},

{MassU,    { LaTeX -> "M_U",
             Form -> Scalar,
             OutputName -> M4,
             LesHouches -> {MSOFT,4} }}

 }; 
 



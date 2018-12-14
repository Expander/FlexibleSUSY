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
    {AlphaS,    { Description -> "Alpha Strong"}},
    {e,         { Description -> "electric charge"}},
    {Gf,        { Description -> "Fermi's constant"}},
    {aEWinv,    { Description -> "inverse weak coupling constant at mZ"}},
    {Yu,        { Description -> "Up-Yukawa-Coupling",
                  DependenceNum ->  Sqrt[2]/v* {{Mass[Fu,1],0,0},
                                                {0, Mass[Fu,2],0},
                                                {0, 0, Mass[Fu,3]}}}},
    {Yd,        { Description -> "Down-Yukawa-Coupling",
                  DependenceNum ->  Sqrt[2]/v* {{Mass[Fd,1],0,0},
                                                {0, Mass[Fd,2],0},
                                                {0, 0, Mass[Fd,3]}}}},
    {Ye,        { Description -> "Lepton-Yukawa-Coupling",
                  DependenceNum ->  Sqrt[2]/v* {{Mass[Fe,1],0,0},
                                                {0, Mass[Fe,2],0},
                                                {0, 0, Mass[Fe,3]}}}},
    {mu2,       { Description -> "SM Mu Parameter",
                  LaTeX -> "m^2",
                  OutputName -> mu2 }},
    {\[Lambda], { Description -> "SM Higgs Selfcouplings",
                  DependenceNum -> Mass[hh]^2/(2 v^2)}},
    {v,         { Description -> "EW-VEV",
                  DependenceNum -> Sqrt[4*Mass[VWp]^2/(g2^2)],
                  DependenceSPheno -> None  }},
    {ThetaW,    { Description -> "Weinberg-Angle",
                  DependenceNum -> ArcSin[Sqrt[1 - Mass[VWp]^2/Mass[VZ]^2]]  }},
    {ZZ,        { Description -> "Photon-Z Mixing Matrix"}},
    {ZW,        { Description -> "W Mixing Matrix",
                  Dependence ->   1/Sqrt[2] {{1, 1},
                                             {\[ImaginaryI],-\[ImaginaryI]}} }},
    {Vu,        { Description ->"Left-Up-Mixing-Matrix"}},
    {Vd,        { Description ->"Left-Down-Mixing-Matrix"}},
    {Uu,        { Description ->"Right-Up-Mixing-Matrix"}},
    {Ud,        { Description ->"Right-Down-Mixing-Matrix"}},
    {Ve,        { Description ->"Left-Lepton-Mixing-Matrix"}},
    {Ue,        { Description ->"Right-Lepton-Mixing-Matrix"}},

(* split MSSM parameters *)

    {ZN,        { Description -> "Neutralino Mixing-Matrix" }},
    {UP,        { Description -> "Chargino-plus Mixing-Matrix" }},
    {UM,        { Description -> "Chargino-minus Mixing-Matrix" }},
    {ZE,        { Description -> "Slepton-Mixing-Matrix" }},
    {ZV,        { Description -> "Sneutrino Mixing-Matrix" }},
    {MassB,     { Description -> "Bino Mass parameter" }},
    {MassWB,    { Description -> "Wino Mass parameter" }},
    {\[Mu],     { Description -> "Mu-parameter",
                  Real -> True,
                  LaTeX -> "\\mu_{\\text{MSSM}}",
                  OutputName -> Mu }},
    {gYu,       { Description -> "Higgs-Bino-Up-Higgsino-Coupling",
                  Real -> True,
                  OutputName -> gYu,
                  LaTeX -> "\\tilde{g}_{1\\text{u}}",
                  LesHouches -> {SplitMSSM,1} }},
    {g2u,       { Description -> "Higgs-Wino-Up-Higgsino-Coupling",
                  Real -> True,
                  OutputName -> g2u,
                  LaTeX -> "\\tilde{g}_{2\\text{u}}",
                  LesHouches -> {SplitMSSM,2} }},
    {gYd,       { Description -> "Higgs-Bino-Down-Higgsino-Coupling",
                  Real -> True,
                  OutputName -> gYd,
                  LaTeX -> "\\tilde{g}_{1\\text{d}}",
                  LesHouches -> {SplitMSSM,3} }},
    {g2d,       { Description -> "Higgs-Wino-Down-Higgsino-Coupling",
                  Real -> True,
                  OutputName -> g2d,
                  LaTeX -> "\\tilde{g}_{2\\text{d}}",
                  LesHouches -> {SplitMSSM,4} }},
    {me2,       { Description -> "Softbreaking right Slepton Mass" }},
    {ml2,       { Description -> "Softbreaking left Slepton Mass" }},
    {Te,        { Description -> "Trilinear-Lepton-Coupling" }},

    {g2llll,    { Description -> "Quartic-Slepton-Coupling",
                  Real -> True,
                  OutputName -> g2llll,
                  LaTeX -> "\\tilde{g}_{2,\\tilde{l}\\tilde{l}\\tilde{l}\\tilde{l}}",
                  LesHouches -> {SplitMSSM,5} }},

    {g2llhh,    { Description -> "Quartic-Slepton-Coupling",
                  Real -> True,
                  OutputName -> g2llhh,
                  LaTeX -> "\\tilde{g}_{2,\\tilde{l}\\tilde{l}\\tilde{l}\\tilde{l}}",
                  LesHouches -> {SplitMSSM,6} }},

    {g2hhhh,    { Description -> "Quartic-Slepton-Coupling",
                  Real -> True,
                  OutputName -> g2hhhh,
                  LaTeX -> "\\tilde{g}_{2,\\tilde{l}\\tilde{l}\\tilde{l}\\tilde{l}}",
                  LesHouches -> {SplitMSSM,7} }},

    {g1llll,    { Description -> "Quartic-Slepton-Coupling",
                  Real -> True,
                  OutputName -> g1llll,
                  LaTeX -> "\\tilde{g}_{2,\\tilde{l}\\tilde{l}\\tilde{l}\\tilde{l}}",
                  LesHouches -> {SplitMSSM,8} }},

    {g1eeee,    { Description -> "Quartic-Slepton-Coupling",
                  Real -> True,
                  OutputName -> g1eeee,
                  LaTeX -> "\\tilde{g}_{2,\\tilde{l}\\tilde{l}\\tilde{l}\\tilde{l}}",
                  LesHouches -> {SplitMSSM,9} }},

    {g1hhhh,    { Description -> "Quartic-Slepton-Coupling",
                  Real -> True,
                  OutputName -> g1hhhh,
                  LaTeX -> "\\tilde{g}_{2,\\tilde{l}\\tilde{l}\\tilde{l}\\tilde{l}}",
                  LesHouches -> {SplitMSSM,10} }},

    {g1llhh,    { Description -> "Quartic-Slepton-Coupling",
                  Real -> True,
                  OutputName -> g1llhh,
                  LaTeX -> "\\tilde{g}_{2,\\tilde{l}\\tilde{l}\\tilde{l}\\tilde{l}}",
                  LesHouches -> {SplitMSSM,11} }},

    {g1eehh,    { Description -> "Quartic-Slepton-Coupling",
                  Real -> True,
                  OutputName -> g1eehh,
                  LaTeX -> "\\tilde{g}_{2,\\tilde{l}\\tilde{l}\\tilde{l}\\tilde{l}}",
                  LesHouches -> {SplitMSSM,12} }},

    {g1llee,    { Description -> "Quartic-Slepton-Coupling",
                  Real -> True,
                  OutputName -> g1llee,
                  LaTeX -> "\\tilde{g}_{2,\\tilde{l}\\tilde{l}\\tilde{l}\\tilde{l}}",
                  LesHouches -> {SplitMSSM,13} }},

    {g2slwl,    { Description -> "Quartic-Slepton-Coupling",
                  Real -> True,
                  OutputName -> g2slwl,
                  LaTeX -> "\\tilde{g}_{2,\\tilde{l}\\tilde{l}\\tilde{l}\\tilde{l}}",
                  LesHouches -> {SplitMSSM,14} }},

    {g1slbl,    { Description -> "Quartic-Slepton-Coupling",
                  Real -> True,
                  OutputName -> g1slbl,
                  LaTeX -> "\\tilde{g}_{2,\\tilde{l}\\tilde{l}\\tilde{l}\\tilde{l}}",
                  LesHouches -> {SplitMSSM,15} }},

    {g1sebe,    { Description -> "Quartic-Slepton-Coupling",
                  Real -> True,
                  OutputName -> g1sebe,
                  LaTeX -> "\\tilde{g}_{2,\\tilde{l}\\tilde{l}\\tilde{l}\\tilde{l}}",
                  LesHouches -> {SplitMSSM,16} }},

    {gylele,    { Description -> "Quartic-Slepton-Coupling",
                  Real -> True,
                  OutputName -> gylele,
                  LaTeX -> "\\tilde{g}_{2,\\tilde{l}\\tilde{l}\\tilde{l}\\tilde{l}}",
                  LesHouches -> {SplitMSSM,17} }},

    {gylehh,    { Description -> "Quartic-Slepton-Coupling",
                  Real -> True,
                  OutputName -> gylehh,
                  LaTeX -> "\\tilde{g}_{2,\\tilde{l}\\tilde{l}\\tilde{l}\\tilde{l}}",
                  LesHouches -> {SplitMSSM,18} }},

    {gydsle,    { Description -> "Triple-down-Higgsino-Slepton-Electron-Coupling",
                  Real -> True,
                  OutputName -> gYedsle,
                  LaTeX -> "\\tilde{g}_{y_e,H_d\\tilde{l}e}",
                  LesHouches -> {SplitMSSM,19} }},

    {gydlse,    { Description -> "Triple-down-Higgsino-Electron-Right-Slepton-Coupling",
                  Real -> True,
                  OutputName -> gYedlse,
                  LaTeX -> "\\tilde{g}_{y_e,H_dl\\tilde{e}}",
                  LesHouches -> {SplitMSSM,20} }}
};

BeginPackage["EffectiveCouplings`", {"SARAH`", "CConversion`", "SelfEnergies`", "TreeMasses`", "Vertices`", "Observables`"}];

CalculateEffectiveCouplings::usage="";

Begin["`Private`"];

GetPhoton[] := SARAH`VectorP;
GetGluon[] := SARAH`VectorG;

CreateEffectiveCouplingName[pIn_, pOut_] :=
    "eff_Cp" <> CConversion`ToValidCSymbolString[pIn] <> CConversion`ToValidCSymbolString[pOut] <> CConversion`ToValidCSymbolString[pOut];

GetAllowedCouplings[] :=
    Module[{dim,
            valid = {FlexibleSUSYObservable`CpHiggsPhotonPhoton,
                     FlexibleSUSYObservable`CpHiggsGluonGluon,
                     FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton,
                     FlexibleSUSYObservable`CpPseudoScalarGluonGluon}
           },
           If[FreeQ[TreeMasses`GetParticles[], SARAH`HiggsBoson],
              valid = DeleteCases[valid, a_ /; (a === FlexibleSUSYObservable`CpHiggsPhotonPhoton ||
                                                a === FlexibleSUSYObservable`CpHiggsGluonGluon)];
             ];
           If[FreeQ[TreeMasses`GetParticles[], SARAH`PseudoScalar],
              valid = DeleteCases[valid, a_ /; (a === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton ||
                                                a === FlexibleSUSYObservable`CpPseudoScalarGluonGluon)];
             ];
           valid
          ];

GetExternalStates[coupling_] :=
    Module[{particle, vectorBoson},
           Which[coupling === FlexibleSUSYObservable`CpHiggsPhotonPhoton,
                 particle = SARAH`HiggsBoson;
                 vectorBoson = SARAH`VectorP;,
                 coupling === FlexibleSUSYObservable`CpHiggsGluonGluon,
                 particle = SARAH`HiggsBoson;
                 vectorBoson = SARAH`VectorG;,
                 coupling === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton,
                 particle = SARAH`PseudoScalar;
                 vectorBoson = SARAH`VectorP;,
                 coupling === FlexibleSUSYObservable`CpPseudoScalarGluonGluon,
                 particle = SARAH`PseudoScalar;
                 vectorBoson = SARAH`VectorG;,
                 True, particle = Null; vectorBoson = Null
                ];
           {particle, vectorBoson}
          ];

CreateEffectiveCouplingPrototype[coupling_] :=
    Module[{particle, vectorBoson, dim, type, name, result = ""},
           {particle, vectorBoson} = GetExternalStates[coupling];
           If[particle =!= Null && vectorBoson =!= Null,
              dim = TreeMasses`GetDimension[particle];
              type = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
              name = CreateEffectiveCouplingName[particle, vectorBoson];
              result = type <> " " <> name <> If[dim == 1, "();\n", "(unsigned gO1);\n"];
             ];
           result
          ];

RunToDecayingParticleScale[particle_] :=
    Module[{savedMass, body, result},
           savedMass = CConversion`RValueToCFormString[FlexibleSUSY`M[particle]];
           If[TreeMasses`GetDimension[particle] != 1,
              savedMass = savedMass <> "(gO1)";
             ];
           body = "model.run_to(" <> savedMass <> ");\nmodel.calculate_DRbar_masses();\n";
           "if (rg_improve && scale != " <> savedMass <> ") {\n"
           <> TextFormatting`IndentText[body] <> "}\n"
          ];

ResetSavedParameters[] :=
    Module[{body},
           body = "model.set_scale(scale);\nmodel.set(saved_parameters);\n";
           "if (model.get_scale() != scale) {\n" <> TextFormatting`IndentText[body] <> "}\n"
          ];

CreateEffectiveCouplingFunction[coupling_] :=
    Module[{particle, vectorBoson, dim, type, name, savedMass, mass, body = "", result = ""},
           {particle, vectorBoson} = GetExternalStates[coupling];
           If[particle =!= Null && vectorBoson =!= Null,
              dim = TreeMasses`GetDimension[particle];
              type = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
              name = CreateEffectiveCouplingName[particle, vectorBoson];
              result = type <> " " <> FlexibleSUSY`FSModelName <> "_effective_couplings::"
                       <> name;
              mass = ToValidCSymbolString[FlexibleSUSY`M[particle]];
              savedMass = "const auto " <> mass <> " = PHYSICAL(" <> mass <> ");\n";
              If[dim == 1,
                 result = result <> "()\n";,
                 result = result <> "(unsigned gO1)\n{\n";
                 mass = mass <> "(gO1)";
                ];
              body = "const double scale = model.get_scale();\n"
                     <> "const Eigen::ArrayXd saved_parameters(model.get());\n"
                     <> savedMass <> "\n"
                     <> RunToDecayingParticleScale[particle] <> "\n"
                     <> type <> " result;\n\n";
              body = body <> ResetSavedParameters[] <> "\n";
              body = body <> "return result;\n";
              result = result <> TextFormatting`IndentText[body] <> "}\n";
             ];
           result
          ];

CreateEffectiveCouplingsPrototypes[couplings_List] :=
    Module[{result = ""},
           (result = result <> CreateEffectiveCouplingPrototype[#])& /@ couplings;
           result
          ];

CreateEffectiveCouplingsFunctions[couplings_List] :=
    Module[{result = ""},
           (result = result <> CreateEffectiveCouplingFunction[#] <> "\n")& /@ couplings;
           result
          ];

CalculateEffectiveCouplings[] :=
    Module[{couplings, prototypes = "", functions = ""},
           couplings = GetAllowedCouplings[];
           (* @todo find needed couplings *)
           prototypes = prototypes <> CreateEffectiveCouplingsPrototypes[couplings];
           functions = functions <> CreateEffectiveCouplingsFunctions[couplings];
           {prototypes, functions}
          ];

(* @todo the following is to be rewritten using a more general approach to getting the vertices *)
(*
GetElectricCharge[state_] := SARAH`getElectricCharge[state];

(* @todo this way of getting the particles coupling to photons and gluons is very slow *)
IsCoupledToPhoton[state_] :=
    Module[{charge, result},
           charge = {GetElectricCharge[state]};
           If[TreeMasses`IsScalar[state] || TreeMasses`IsVector[state],
              If[!NumericQ[charge],
                 charge = Cases[SARAH`Vertex[{Susyno`LieGroups`conj[state],
                                              state, GetPhoton[]}, UseDependences -> True][[2,1]], _?NumberQ];
                ];
              result = charge =!= {} && Susyno`LieGroups`conj[state] =!= state && charge =!= {0};,
              If[!NumericQ[charge],
                 charge = Cases[SARAH`Vertex[{SARAH`bar[state], state, GetPhoton[]},
                                             UseDependences -> True][[2,1]], _?NumberQ];
                ];
              result = charge =!= {} && SARAH`bar[state] =!= state && charge =!= {0};
             ];
           result
          ];

IsCoupledToGluon[state_] :=
    Module[{charge, result},
           If[TreeMasses`IsFermion[state],
              charge = SARAH`Vertex[{SARAH`bar[state], state, GetGluon[]}][[2,1]];,
              charge = SARAH`Vertex[{Susyno`LieGroups`conj[state], state, GetGluon[]}][[2,1]];
             ];
           charge =!= 0
          ];

GetParticlesCoupledToPhoton[particles_List] := Select[particles, (IsCoupledToPhoton[#] && !TreeMasses`IsGhost[#])&];

GetParticlesCoupledToGluon[particles_List] := Select[particles, (IsCoupledToGluon[#] && !TreeMasses`IsGhost[#])&];

(* hacked (temporary) solution: get the two-body decays this way *)
GetTwoBodyDecays[state_, nPointFunctions_, eigenstates_:FlexibleSUSY`FSEigenstates] :=
    Module[{dim, selfEnergy},
           dim = TreeMasses`GetDimension[state, eigenstates];
           If[dim === 1,
              selfEnergy = First[Cases[nPointFunctions, SelfEnergies`FSSelfEnergy[state, _]]];,
              selfEnergy = First[Cases[nPointFunctions, SelfEnergies`FSSelfEnergy[state[__], _]]];
             ];
           twoBodyDecays = DeleteDuplicates[Vertices`GetParticleList /@
                                            Join[Cases[selfEnergy, SARAH`Cp[p1_, p2_, p3_]
                                                       :> (SARAH`Cp @@ Vertices`ToRotatedField /@ {p1, p2, p3}), {0, Infinity}],
                                                 Cases[selfEnergy, SARAH`Cp[p1_, p2_, p3_][PL]
                                                       :> (SARAH`Cp @@ Vertices`ToRotatedField /@ {p1, p2, p3})[PL], {0, Infinity}],
                                                 Cases[selfEnergy, SARAH`Cp[p1_, p2_, p3_][PR]
                                                       :> (SARAH`Cp @@ Vertices`ToRotatedField /@ {p1, p2, p3})[PR], {0,Infinity}]]
                                            /. p_[{__}] -> p /. SARAH`Cp[p__][PL | PR] -> SARAH`Cp[p]];
           twoBodyDecays = DeleteCases[twoBodyDecays, SARAH`Cp[p1_, p2_, p3_] /; Or @@ (TreeMasses`IsGhost /@ {p1, p2, p3})];
           twoBodyDecays /. SARAH`Cp[p1_, p2_, p3_] -> {p2, p3, SARAH`Cp[p1, p2, p3]}
          ];

GetEffectiveCouplingToPhoton[state_, nPointFunctions_, eigenstates_:FlexibleSUSY`FSEigenstates] :=
    Module[{particles = TreeMasses`GetParticles[eigenstates],
            chargedParticles, decayParticles, neededCouplings},
           chargedParticles = GetParticlesCoupledToPhoton[particles];
           decayParticles = GetTwoBodyDecays[state, nPointFunctions, eigenstates];

          ];

GetEffectiveCouplingToGluon[state_, eigenstates_:FlexibleSUSY`FSEigenstates] :=
    Module[{particles = TreeMasses`GetParticles[eigenstates],
            coloredParticles, decayParticles, neededCouplings = {}},
           coloredParticles = Cases[GetParticlesCoupledToGluon[particles], Except[GetGluon[]]];
          ];

*)

End[];

EndPackage[];

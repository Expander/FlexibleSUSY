BeginPackage["EffectiveCouplings`", {"SARAH`", "CConversion`", "SelfEnergies`", "TreeMasses`", "Vertices`", "Observables`"}];

CalculateEffectiveCouplings::usage="";

Begin["`Private`"];

GetPhoton[] := SARAH`VectorP;
GetGluon[] := SARAH`VectorG;

CreateEffectiveCouplingName[pIn_, pOut_] :=
    "eff_Cp" <> CConversion`ToValidCSymbolString[pIn] <> CConversion`ToValidCSymbolString[pOut] <> CConversion`ToValidCSymbolString[pOut];

IsImplementedCoupling[obs_ /; obs === FlexibleSUSYObservable`CpHiggsPhotonPhoton] := True;
IsImplementedCoupling[obs_ /; obs === FlexibleSUSYObservable`CpHiggsGluonGluon] := True;
IsImplementedCoupling[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] := True;
IsImplementedCoupling[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarGluonGluon] := True;
IsImplementedCoupling[_] := False;

FilterOutInvalidCouplings[couplings_List] :=
    Module[{dim, valid},
           valid = couplings;
           If[MemberQ[couplings, FlexibleSUSYObservable`CpHiggsPhotonPhoton] ||
              MemberQ[couplings, FlexibleSUSYObservable`CpHiggsGluonGluon],
              dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson]
              If[FreeQ[TreeMasses`GetParticles[], SARAH`HiggsBoson] ||
                 TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson] == 0,
                 Print["Warning: no physical Higgs boson found."];
                 Print["         Effective couplings for Higgs boson will not"];
                 Print["         be calculated."];
                 valid = DeleteCases[valid, a_ /; (a === FlexibleSUSYObservable`CpHiggsPhotonPhoton ||
                                                   a === FlexibleSUSYObservable`CpHiggsGluonGluon)];
                ];
             ];
           If[MemberQ[couplings, FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] ||
              MemberQ[couplings, FlexibleSUSYObservable`CpPseudoScalarGluonGluon],
              If[FreeQ[TreeMasses`GetParticles[], SARAH`PseudoScalar] ||
                 TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar] == 0,
                 Print["Warning: no physical pseudoscalar boson found."];
                 Print["         Effective couplings for pseudoscalar boson will not"];
                 Print["         be calculated."];
                 valid = DeleteCases[valid, a_ /; (a === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton ||
                                                   a === FlexibleSUSYObservable`CpPseudoScalarGluonGluon)]; 
                ];
             ];
           valid
          ];

CreateEffectiveCouplingPrototype[coupling_] :=
    Module[{particle, vectorBoson, dim, type, name},
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
                 True, Return[""]
                ];
           dim = TreeMasses`GetDimensionWithoutGoldstones[particle];
           type = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
           name = CreateEffectiveCouplingName[particle, vectorBoson];
           type <> " " <> name <> If[dim == 1, "();\n", "(unsigned gO1);\n"]
          ];

CreateEffectiveCouplingsPrototypes[couplings_List] :=
    Module[{result = ""},
           (result = result <> CreateEffectiveCouplingPrototype[#])& /@ couplings;
           result
          ];

CreateEffectiveCouplingsFunctions[couplings_List] :=
    Module[{result = ""},
           result
          ];

CalculateEffectiveCouplings[observables_List] :=
    Module[{couplings, prototypes = "", functions = ""},
           couplings = DeleteDuplicates[Cases[observables, p_ /; IsImplementedCoupling[p], {0, Infinity}]];
           couplings = FilterOutInvalidCouplings[observables];
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

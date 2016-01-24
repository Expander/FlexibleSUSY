BeginPackage["EffectiveCouplings`", {"SARAH`", "CConversion`", "SelfEnergies`", "TreeMasses`", "TextFormatting`", "Vertices`", "Observables`"}];

InitializeEffectiveCouplings::usage="";
GetNeededVerticesList::usage="";
CreateEffectiveCouplings::usage="";

Begin["`Private`"];

GetAllowedCouplingsForModel[] :=
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

(* @todo much of this is either identical or very similar to required functions
   for the GMuonMinus2 branch, and will be required again for decays, it might
   be a good idea to make them general (e.g. defined in Vertices?).             *)
AntiParticle[p_] := If[TreeMasses`IsScalar[p] || TreeMasses`IsVector[p],
                       Susyno`LieGroups`conj[p],
                       SARAH`bar[p]];

NonZeroVertexQ[vertex_] := MemberQ[vertex[[2 ;;]][[All, 1]], Except[0]];

(* @todo extend to multiple non-Abelian groups *)
IsColorOrLorentzIndex[index_] := StringMatchQ[ToString @ index, "ct" ~~ __] ||
                                 StringMatchQ[ToString @ index, "lt" ~~ __];
StripColorAndLorentzIndices[p_Symbol] := p;
StripColorAndLorentzIndices[SARAH`bar[p_]] := SARAH`bar[StripColorAndLorentzIndices[p]];
StripColorAndLorentzIndices[Susyno`LieGroups`conj[p_]] := Susyno`LieGroups`conj[StripColorAndLorentzIndices[p]];
StripColorAndLorentzIndices[p_] :=
    Module[{remainingIndices},
           remainingIndices = Select[p[[1]], (!IsColorOrLorentzIndex[#])&];
           If[Length[remainingIndices] === 0,
              Head[p],
              Head[p][remainingIndices]
             ]
          ];
SetAttributes[StripColorAndLorentzIndices, {Listable}];

(* @todo this is very slow because each possible vertex must be calculated,
   this can be improved by either pre-calculating all vertices or saving
   previous results (e.g. define along the lines of f[p] := f[p] = ...)    *)
GetTwoBodyDecays[particle_] :=
    Module[{i, j, allParticles, combinations, vertex, fields, couplings,
            candidate, found = {}, decays = {}},
           allParticles = Select[TreeMasses`GetParticles[], !TreeMasses`IsGhost[#]&];
           combinations = Table[Sort[{AntiParticle[particle],
                                      AntiParticle[allParticles[[i]]],
                                      allParticles[[j]]}],
                                {i, 1, Length[allParticles]}, {j, 1, Length[allParticles]}];
           combinations = DeleteDuplicates[Flatten[combinations, 1]];
           For[i = 1, i <= Length[combinations], i++,
               vertex = SARAH`Vertex[combinations[[i]], UseDependences -> True];
               If[NonZeroVertexQ[vertex],
                  fields = First[vertex];
                  coupling = Rest[vertex];
                  If[Length[coupling] > 1,
                     coupling = (SARAH`Cp @@ (StripColorAndLorentzIndices @ fields))[SARAH`PL];,
                     coupling = SARAH`Cp @@ (StripColorAndLorentzIndices @ fields);
                    ];
                  candidate = Append[DeleteCases[fields /. head_[{__}] :> head, p_ /; p === AntiParticle[particle], {0, Infinity}, 1], coupling];
                  If[FreeQ[found, C[candidate[[1]], candidate[[2]]]] &&
                     FreeQ[found, C[AntiParticle[candidate[[1]]], AntiParticle[candidate[[2]]]]] ||
                     AntiParticle[particle] =!= particle,
                     If[((!TreeMasses`IsMassless[candidate[[1]]] || !TreeMasses`IsVector[candidate[[1]]]) &&
                        (!TreeMasses`IsMassless[candidate[[2]]] || !TreeMasses`IsVector[candidate[[2]]])) ||
                        !FreeQ[SARAH`AllowDecaysMasslessVectors, SARAH`RE[particle]],
                        decays = Append[decays, candidate];
                        found = Append[found, C[candidate[[1]], candidate[[2]]]];
                       ];
                    ];
                 ];
              ];
           decays
          ];

GetElectricCharge[p_] := SARAH`getElectricCharge[p];

GetParticlesCouplingToVectorBoson[vector_] :=
    Module[{i, charge, allParticles, particles = {}},
           allParticles = Select[TreeMasses`GetParticles[], !TreeMasses`IsGhost[#]&];
           For[i = 1, i <= Length[allParticles], i++,
               If[vector === SARAH`VectorG,
                  (* @note could use defined functions in e.g. TreeMasses, plus check
                     for undefined group factors, but will do it this way for now
                     to ensure consistency with SARAH                                  *)
                  charge = SARAH`Vertex[{AntiParticle[allParticles[[i]]],
                                         allParticles[[i]], SARAH`VectorG},
                                        UseDependences -> True][[2,1]];
                  If[charge =!= 0,
                     particles = Append[particles, allParticles[[i]]];
                    ];,
                  charge = GetElectricCharge[allParticles[[i]]];
                  If[NumericQ[charge],
                     charge = {charge},
                     charge = Cases[SARAH`Vertex[{AntiParticle[allParticles[[i]]],
                                                  allParticles[[i]], SARAH`VectorP},
                                                 UseDependences -> True][[2,1]], _?NumberQ];
                    ];
                  If[charge =!= {} && AntiParticle[allParticles[[i]]] =!= allParticles[[i]] &&
                     charge =!= {0},
                     particles = Append[particles, allParticles[[i]]];
                    ];
                 ];
              ];
           particles
          ];

InitializeEffectiveCouplings[] :=
    Module[{i, couplings, particle, vectorBoson,
            allParticles = {}, allVectorBosons = {},
            twoBodyDecays, vectorBosonInteractions,
            neededTwoBodyDecays, neededVectorBosonInteractions,
            neededCoups, result = {}},
           couplings = GetAllowedCouplingsForModel[];
           {allParticles, allVectorBosons} = DeleteDuplicates /@ {(#[[1]])& /@ (GetExternalStates[#]& /@ couplings),
                                                                  (#[[2]])& /@ (GetExternalStates[#]& /@ couplings)};
           twoBodyDecays = {#, GetTwoBodyDecays[#]}& /@ allParticles;
           vectorBosonInteractions = {#, GetParticlesCouplingToVectorBoson[#]}& /@ allVectorBosons;
           For[i = 1, i <= Length[couplings], i++,
                {particle, vectorBoson} = GetExternalStates[couplings[[i]]];
                neededTwoBodyDecays = First[Select[twoBodyDecays, (#[[1]] === particle)&]];
                neededVectorBosonInteractions = First[Select[vectorBosonInteractions, (#[[1]] === vectorBoson)&]];
                neededCoups = Select[neededTwoBodyDecays[[2]],
                                    (MemberQ[neededVectorBosonInteractions[[2]], #[[1]]] ||
                                     MemberQ[neededVectorBosonInteractions[[2]], #[[2]]])&];
                result = Append[result, {couplings[[i]], #[[3]]& /@ neededCoups}];
              ];
           result
          ];

GetNeededVerticesList[couplings_List] :=
    {Null[Null, Join[(#[[2]])& /@ couplings]]};

CreateEffectiveCouplingName[pIn_, pOut_] :=
    "eff_Cp" <> CConversion`ToValidCSymbolString[pIn] <> CConversion`ToValidCSymbolString[pOut] <> CConversion`ToValidCSymbolString[pOut];

(* @todo these are basically identical to those in SelfEnergies,
   it would be better to reuse the definitions there if possible  *)
GetParticleIndicesInCoupling[SARAH`Cp[a__]] := Flatten[Cases[{a}, List[__], Infinity]];

GetParticleIndicesInCoupling[SARAH`Cp[a__][_]] := GetParticleIndices[SARAH`Cp[a]];

CreateEffectiveCouplingPrototype[coupling_] :=
    Module[{couplingName = coupling[[1]], particle, vectorBoson,
            dim, type, name, result = ""},
           {particle, vectorBoson} = GetExternalStates[couplingName];
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
    Module[{couplingName = coupling[[1]], neededCouplings = coupling[[2]],
            particle, vectorBoson, dim, type, name, savedMass, mass, body = "", result = ""},
           {particle, vectorBoson} = GetExternalStates[couplingName];
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
                     <> savedMass
                     <> "const double scale_factor = 0.25 * Sqr(" <> mass <> ");\n\n"
                     <> RunToDecayingParticleScale[particle] <> "\n"
                     <> type <> " result;\n\n";

              (* @todo convert expression to code *)

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

CreateEffectiveCouplings[couplings_List] :=
    Module[{prototypes = "", functions = ""},
           prototypes = prototypes
                        <> CreateEffectiveCouplingsPrototypes[couplings];
           functions = functions
                       <> CreateEffectiveCouplingsFunctions[couplings];
           {prototypes, functions}
          ];

End[];

EndPackage[];

BeginPackage["EffectiveCouplings`", {"SARAH`", "CConversion`", "SelfEnergies`", "TreeMasses`", "Vertices`", "Observables`"}];

CalculateEffectiveCouplings::usage="";

Begin["`Private`"];

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

(* @todo much of this is either identical or very similar to required functions
   for the GMuonMinus2 branch, and will be required again for decays, it might
   be a good idea to make them general (e.g. defined in Vertices?).             *)
AntiParticle[p_] := If[TreeMasses`IsScalar[p] || TreeMasses`IsVector[p],
                       Susyno`LieGroups`conj[p],
                       SARAH`bar[p]];

NonZeroVertexQ[vertex_] := MemberQ[vertex[[2 ;;]][[All, 1]], Except[0]];

(* @todo this is very slow because each possible vertex must be calculated,
   this can be improved by either pre-calculating all vertices or saving
   previous results (e.g. define along the lines of f[p] := f[p] = ...)    *)
GetTwoBodyDecays[particle_] :=
    Module[{i, j, allParticles, combinations, vertex,
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
                  candidate = Append[DeleteCases[vertex[[1]] /. head_[{__}] :> head, p_ /; p === AntiParticle[particle], {0, Infinity}, 1], SARAH`Cp @@ (vertex[[1]] /. head_[{__}] :> head)];
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

GetEffectiveCouplingExpression[coupling_] :=
    Module[{particle, vectorBoson, twoBodyDecays,
            vectorBosonInteractions, expr, neededCouplings},
           Which[couplings[[i]] === FlexibleSUSYObservable`CpHiggsPhotonPhoton,
                 particle = SARAH`HiggsBoson;
                 vectorBoson = SARAH`VectorP;,
                 couplings[[i]] === FlexibleSUSYObservable`CpHiggsGluonGluon,
                 particle = SARAH`HiggsBoson;
                 vectorBoson = SARAH`VectorG;,
                 couplings[[i]] === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton,
                 particle = SARAH`PseudoScalar;
                 vectorBoson = SARAH`VectorP;,
                 couplings[[i]] === FlexibleSUSYObservable`CpPseudoScalarGluonGluon,
                 particle = SARAH`PseudoScalar;
                 vectorBoson = SARAH`VectorG;,
                 True,
                 Print["Error: unsupported coupling ", couplings[[i]]];
                 Quit[1];
                ];

           twoBodyDecays = GetTwoBodyDecays[particle];
           vectorBosonInteractions = GetParticlesCouplingToVectorBoson[vectorBoson];
           neededCouplings = Select[twoBodyDecays,
                                    (MemberQ[vectorBosonInteractions, #[[1]]] ||
                                     MemberQ[vectorBosonInteractions, #[[2]]])&];

           {{coupling, expr}, neededCouplings}
          ];

GetEffectiveCouplingsExpressions[couplings_] :=
    Module[{i, particle, vectorBoson, expr, neededCoups,
            expressions = {}, neededCouplings = {}},
           For[i = 1, i <= Length[couplings], i++,
                {expr, neededCoups} = GetEffectiveCouplingExpression[couplings[[i]]];
                expressions = Append[expressions, expr];
                neededCouplings = DeleteDuplicates[Join[neededCouplings, neededCoups]];
              ];
           {expressions, neededCouplings}
          ];

CreateEffectiveCouplingName[pIn_, pOut_] :=
    "eff_Cp" <> CConversion`ToValidCSymbolString[pIn] <> CConversion`ToValidCSymbolString[pOut] <> CConversion`ToValidCSymbolString[pOut];

(* @todo these are basically identical to those in SelfEnergies,
   it would be better to reuse the definitions there if possible  *)
GetParticleIndicesInCoupling[SARAH`Cp[a__]] := Flatten[Cases[{a}, List[__], Infinity]];

GetParticleIndicesInCoupling[SARAH`Cp[a__][_]] := GetParticleIndices[SARAH`Cp[a]];

CreateNeededCouplingSymbol[coupling_] :=
    Module[{symbol, indices},
           indices = GetParticleIndicesInCoupling[coupling];
           symbol = CConversion`ToValidCSymbol[coupling /. a_[List[__]] :> a];
           symbol[Sequence @@ indices]
          ];

CreateNeededCouplingName[coupling_] :=
    CConversion`ToValidCSymbolString[GetHead[CreateNeededCouplingSymbol[coupling]]];

CreateNeededCouplingFunction[coupling_, expr_] :=
    Module[{i, indices, functionName, initialValue, type, typeStr,
            body, prototype, definition},
           indices = GetParticleIndicesInCoupling[coupling];
           functionName = CreateNeededCouplingName[coupling];
           functionName = functionName <> "(";
           For[i = 1, i <= Length[indices], i++,
               If[i > 1, functionName = functionName <> ", ";];
               functionName = functionName <> "unsigned ";
               (* variable names must not be integers *)
               If[!IntegerQ[indices[[i]]] && !FreeQ[expr, indices[[i]]],
                  functionName = functionName <> CConversion`ToValidCSymbolString[indices[[i]]];
                 ];
              ];
           functionName = functionName <> ")";
           If[Parameters`IsRealExpression[expr],
              type = CConversion`ScalarType[CConversion`realScalarCType];    initialValue = " = 0.0";,
              type = CConversion`ScalarType[CConversion`complexScalarCType]; initialValue = "";];
           typeStr = CConversion`CreateCType[type];
           prototype = typeStr <> " " <> functionName <> " const;\n";
           definition = typeStr <> " " <> FlexibleSUSY`FSModelName <> "_effective_couplings::"
                        <> functionName <> " const\n{\n";
           body = Parameters`CreateLocalConstRefs[expr] <> "\n" <>
                  typeStr <> " result" <> initialValue <> ";\n\n";
           body = body <> TreeMasses`ExpressionToString[expr, "result"];
           body = body <> "\nreturn result;\n";
           body = IndentText[WrapLines[body]];
           definition = definition <> body <> "}\n";
           {prototype, definition}
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

CreateEffectiveCouplingFunction[coupling_, expr_] :=
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
                     <> savedMass
                     <> "const double scale_factor = 0.25 * Sqr(" <> mass <> ");\n\n";
                     <> RunToDecayingParticleScale[particle] <> "\n"
                     <> type <> " result;\n\n";

              (* @todo convert expression to code *)

              body = body <> ResetSavedParameters[] <> "\n";
              body = body <> "return result;\n";
              result = result <> TextFormatting`IndentText[body] <> "}\n";
             ];
           result
          ];

CreateNeededCouplingsFunctions[couplings_List] :=
    Module[{i, proto, func, prototypes = "", functions = ""},
           For[i = 1, i <= Length[couplings], i++,
               {proto, func} = CreateNeededCouplingFunction[couplings[[i,1]], couplings[[i,2]]];
               prototypes = prototypes <> proto;
               functions = functions <> func <> "\n";
              ];
           {prototypes, functions}
          ];

CreateEffectiveCouplingsPrototypes[couplings_List] :=
    Module[{result = ""},
           (result = result <> CreateEffectiveCouplingPrototype[#])& /@ couplings;
           result
          ];

CreateEffectiveCouplingsFunctions[couplings_List] :=
    Module[{result = ""},
           (result = result <> CreateEffectiveCouplingFunction[#[[1]], #[[2]]] <> "\n")& /@ couplings;
           result
          ];

CalculateEffectiveCouplings[] :=
    Module[{couplings, couplingsAndExprs, neededCouplings,
            neededPrototypes, neededFunctions,
            prototypes = "", functions = ""},
           couplings = GetAllowedCouplings[];
           {couplingsAndExprs, neededCouplings} = GetEffectiveCouplingsExpressions[couplings];
           {neededPrototypes, neededFunctions} = CreateNeededCouplingsFunctions[neededCouplings];
           prototypes = prototypes <> neededPrototypes
                        <> CreateEffectiveCouplingsPrototypes[couplings];
           functions = functions <> neededFunctions
                       <> CreateEffectiveCouplingsFunctions[couplingAndExprs];
           {prototypes, functions}
          ];

End[];

EndPackage[];

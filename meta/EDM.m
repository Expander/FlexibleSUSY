BeginPackage["EDM`", {"SARAH`", "TextFormatting`", "TreeMasses`", "LoopMasses`", "Vertices`"}];

Test::usage="";

SetEDMParticles::usage="Set the particles for which the EDMs shall be calculated.";

CreateParticles::usage="Returns the c++ code that contains all particle classes";
CreateEDMParticleFunctions::usage="Returns the c++ code that contains all EDM particle functions";
CreateDiagrams::usage="Returns the c++ code that contains all relevant diagram classes";
CreateVertexFunctionData::usage="Returns the c++ code that contains all relevant vertex function data";

CreateCalculation::usage="Returns the c++ code that performs the actual calculation the magnetic moment";

CreateDefinitions::usage="Returns the c++ that contains all function definitions"

NPointFunctions::usage="Returns a list of all n point functions that are needed. Actually it is a list of fake functions to extract vertex functions...";

(******** IMPORTANT NOTES:
 If you add new kinds of vertices (e.g for new diagram types):
 - Add the new types to vertexTypes
 - Expand CouplingsForParticles[] and VertexTypeForParticles[] accordingly
 - Write the c++ class for the new vertex type

 When adding support for new diagram types, do the following:
 - Add the new types to contributingDiagramTypes
 - Write new overloads for CreateDiagramEvaluatorClass[], ContributingDiagramsOfType[] and VerticesForDiagram[]
 - Write the necessary c++ code: loop functions, DiagramEvaluator<> specialisations
 **********)

(* TODO: privatize interface again Begin["`Private`"]; *)

(************* Begin public interface *******************)

edmParticles = Null;
SetEDMParticles[particles_List] := (edmParticles = particles;)

IsSMParticle[particle_] :=
    SARAH`SMQ[particle] || TreeMasses`IsSMGoldstone[particle];

CreateSMParticleFlags[particle_] :=
    Module[{result = "", i,
            numberOfGenerations = TreeMasses`GetDimension[particle]},
           For[i = 1, i <= numberOfGenerations, i++,
               If[i > 1, result = result <> ", ";];
               If[IsSMParticle[particle[i]] === True ||
                  IsSMParticle[particle] === True,
                  result = result <> "true";,
                  result = result <> "false";
                 ];
              ];
           "{ " <> result <> " }"
          ];

(* Create c++ classes for all particles *)
CreateParticles[] :=
    Module[{particles, code},
           (* Get a list of all particles *)
           particles = TreeMasses`GetParticles[];

           code = ("// Particles (SARAH-style)\n" <>
                   "struct Particle {};\n\n" <>

                   StringJoin @ Riffle[("struct " <> ParticleToCXXName[#] <>
                                        ": public Particle {\n" <>
                                        TextFormatting`IndentText[
                                            "static const unsigned numberOfGenerations = " <>
                                            ToString @ TreeMasses`GetDimension[#] <> ";\n" <>
                                            "static const bool is_sm_particle[numberOfGenerations];"
                                                                 ] <> "\n};\n" <>
                                        "const bool " <> ParticleToCXXName[#] <> "::is_sm_particle[" <>
                                        ParticleToCXXName[#] <> "::numberOfGenerations] = " <>
                                        CreateSMParticleFlags[#] <> ";\n"
                                        &) /@ particles, "\n"] <> "\n\n" <>
                   "// Special particle families\n" <>
                   "typedef " <> ParticleToCXXName @ GetPhoton[] <> " Photon;\n" <>
                   "typedef " <> ParticleToCXXName @ GetMuonFamily[] <> " MuonFamily;\n\n" <>

                   "// AntiFields\n" <>
                   "template<class P> struct anti : public Particle\n" <>
                   "{\n" <>
                   IndentText @
                   ("static const unsigned numberOfGenerations = P::numberOfGenerations;\n" <>
                    "typedef anti<P> type;\n") <>
                   "};\n" <>
                   "template<class P> struct anti<anti<P>> { typedef P type; };\n\n" <>

                   "// Particles that are their own AntiFields\n" <>
                   StringJoin @ Riffle[("template<> struct " <>
                                        "anti<" <> ParticleToCXXName[#] <> ">" <>
                                        " { typedef " <> ParticleToCXXName[#] <> " type; };"
                                        &) /@ Select[particles, (# == AntiField[#] &)],
                                       "\n"]
                  );

           Return[code];
          ];

CreateChargeGetters[] :=
    "template<class Particle>\n" <>
    "double charge( EvaluationContext & );\n\n" <>
    "template<class Particle>\n" <>
    "double charge( unsigned, EvaluationContext & );\n\n" <>
    StringJoin @ Riffle[
        (Module[{photonVertexParticles, particleDim = TreeMasses`GetDimension[#]},
           photonVertexParticles = {GetPhoton[], #, AntiField[#]};
           "template<>\n" <>
           "double charge<" <> ParticleToCXXName[#] <> ">" <>
           If[particleDim === 1, "( ",
              "( unsigned index, "] <> "EvaluationContext& context )\n{\n" <>
           IndentText[
               "typedef VertexFunction<" <>
               StringJoin @ Riffle[ParticleToCXXName /@ photonVertexParticles, ", "] <>
               "> VF;\n" <>
               If[particleDim === 1,
                  "const std::array<unsigned, 0> indices{};\n",
                  "const std::array<unsigned, 2> indices{ index, index };\n"
                 ] <>
               "return VF::vertex(indices, context).left().real();"
           ] <> "\n}"
          ] &) /@ edmParticles,
                        "\n\n"]

CreateDiagrams[] :=
    Module[{diagramTypes, diagramTypeHeads, code},
           diagrams = contributingDiagramTypes;
           diagramHeads = DeleteDuplicates @ (Head /@ diagrams);

           code = "// The different diagram types that contribute to the muon magnetic moment\n";
           code = (code <>
                   StringJoin @ Riffle[("template<unsigned> class " <> SymbolName[#] <> ";" &)
                                       /@ diagramHeads, "\n"] <>
                   "\n\n");

           code = (code <> "// Indexed diagram types\n" <>
                   StringJoin @ Riffle[("template<> class " <> SymbolName[Head[#]] <>
                                        "<" <> ToString @ #[[1]] <> "> {};" &)
                                       /@ diagrams, "\n"]);

           code = (code <> "\n\n" <>
                   StringJoin @ Riffle[CreateDiagramEvaluatorClass /@ contributingDiagramTypes, "\n\n"]);

           Return[code];
          ];

CreateDiagramEvaluatorClass[type_OneLoopDiagram] :=
    ("template<class EDMParticle, class PhotonEmitter, class ExchangeParticle>\n" <>
     "struct DiagramEvaluator<OneLoopDiagram<" <>
     ToString @ type[[1]] <>
     ">, EDMParticle, PhotonEmitter, ExchangeParticle>\n" <>
     "{ static double value(EvaluationContext& context); };");

CreateVertexFunctionData[vertexRules_List] := CreateVertices[vertexRules][[1]];

calculationCode = Null;
CreateCalculation[] :=
    Module[{code, evaluators},
           (* If we have been here before return the old result *)
           If[calculationCode =!= Null, Return[calculationCode]];
           
           evaluators = ConcreteDiagramEvaluators[];
           
           code = "/********** EDM.m generated calculation code **********/\n\n" <>
                  "template<class Particle> double edm( void );\n\n";
           
           code = code <> StringJoin @ Riffle[
                  Module[{pEvaluators = Cases[evaluators, {#, ev_List} -> ev]},
                   "template<> double edm<" <> ParticleToCXXName[#] <> ">( void )\n" <>
                   "{\n" <>
                   IndentText["EvaluationContext context{ model };\n" <>
                              "double val = 0.0\n\n;" <>
                              StringJoin @ Riffle[("val += " <> # <> "::value(context);" &) /@ pEvaluators,
                                                  "\n"] <> "\n\n" <>
                              "return val;"
                   ] <>
                   "}"] & /@ evaluators, "\n\n"];

           calculationCode = code;
           Return[code];
          ];

CreateDefinitions[vertexRules_List] :=
    (CreateEvaluationContextSpecializations[] <> "\n\n" <>
     CreateVertices[vertexRules][[2]]);

nPointFunctions = Null;
NPointFunctions[] :=
    Module[{contributingDiagrams, vertices},
           If[nPointFunctions =!= Null, Return[nPointFunctions]];

           contributingDiagrams = ContributingDiagrams[];

           vertices = Flatten[VerticesForDiagram /@ contributingDiagrams, 1];
           AppendTo[vertices, StripLorentzIndices @ MemoizingVertex[{GetPhoton[], GetMuonFamily[], SARAH`bar[GetMuonFamily[]]}][[1]]];
           vertices = DeleteDuplicates[vertices];

           vertices = (OrderParticles[#, Ordering[(Vertices`StripFieldIndices /@ #)]] &) /@ vertices;
           vertices = DeleteDuplicates[vertices,
                                       (Vertices`StripFieldIndices[#1] === Vertices`StripFieldIndices[#2] &)];

           nPointFunctions = Flatten[(Null[Null, #] &) /@ ((CouplingsForParticles[#] &) /@ vertices)];
           Return[nPointFunctions];
          ];

(**************** End public interface *****************)

(* Effectively generate all mass calculation functions *)
CreateEvaluationContextSpecializations[] :=
Module[{particles, code},
       particles = TreeMasses`GetParticles[];
       particles = Select[particles, (! TreeMasses`IsGhost[#] &)];

       code = (StringJoin @
               Riffle[("template<> double EvaluationContext::mass<" <> ToString[#] <> ">(" <>
                       If[TreeMasses`GetDimension[#] === 1, "", "unsigned index"] <> ") const\n" <>
                       "{ return model.get_M" <> ParticleToCXXName[#] <>
                       If[TreeMasses`GetDimension[#] === 1, "()", "(index)"] <> "; }"
                       &) /@ particles, "\n\n"]);

       Return[code];
       ];

(************************ Begin helper routines *******************************)

GetPhoton[] := SARAH`Photon;

IsLorentzIndex[index_] := StringMatchQ[ToString @ index, "lt" ~~ __];

StripLorentzIndices[p_Symbol] := p;
StripLorentzIndices[SARAH`bar[p_]] := SARAH`bar[StripLorentzIndices[p]];
StripLorentzIndices[Susyno`LieGroups`conj[p_]] := Susyno`LieGroups`conj[StripLorentzIndices[p]];
StripLorentzIndices[p_] := Module[{remainingIndices},
                                  remainingIndices = Select[p[[1]], (!IsLorentzIndex[#] &)];
                                  If[Length[remainingIndices] === 0, Head[p],
                                     Head[p][remainingIndices]]
                                  ];
SetAttributes[StripLorentzIndices, {Listable}];

(* Return a string corresponding to the c++ class name of the particle.
 Note that "bar" and "conj" get turned into anti<...>::type! *)
ParticleToCXXName[p_] := SymbolName[p];
ParticleToCXXName[SARAH`bar[p_]] := "anti<" <> SymbolName[p] <> ">::type";
ParticleToCXXName[Susyno`LieGroups`conj[p_]] := "anti<" <> SymbolName[p] <> ">::type";

(* Return a string corresponding to the name of the particle.
 Note that "bar" and "conj" are left as they are! *)
ParticleToSARAHString[p_] := SymbolName[p];
ParticleToSARAHString[SARAH`bar[p_]] := "bar" <> SymbolName[p];
ParticleToSARAHString[Susyno`LieGroups`conj[p_]] := "conj" <> SymbolName[p];

subIndexPattern = (ReplacePart[SARAH`subIndizes[[All, 1]], 0 -> Alternatives] -> ___);
AddIndexPattern[SARAH`bar[p_]] := SARAH`bar[AddIndexPattern[p]];
AddIndexPattern[Susyno`LieGroups`conj[p_]] := Susyno`LieGroups`conj[AddIndexPattern[p]];
AddIndexPattern[particle_] := SARAH`getFull[SARAH`getBlank[particle]] /. subIndexPattern;

CachedVertex[particles_List] :=
    Module[{
        vertexPattern = ReplacePart[({#, ___} &) /@
                                        Permutations[AddIndexPattern /@ particles],
                                        0 -> Alternatives],
            vertexList = Symbol["SARAH`VertexList" <> ToString @ Length[particles]]},
           FirstCase[vertexList, vertexPattern]
           ];

(* Returns the name of the coupling function that FlexibleSUSY generates for
 a specific vertex in a canonical order! *)
NameOfCouplingFunction[particles_List] :=
((* FIXME: Not upwards compatible if naming conventions change *)
 "Cp" <> StringJoin @ (ParticleToSARAHString /@ Sort[particles]));

(********************** End helper routines **************************)

(* The different vertex types that are supported.
 They have the same names as their c++ counterparts. *)
vertexTypes = {
    SingleComponentedVertex,
    LeftAndRightComponentedVertex
};

(* The different diagram types that should be taken into consideration *)
(* They need to be called DIAGRAMTYPENAME[_Integer]! See CreateDiagramClasses[] below. *)
(* There is no bounds check done on the integers, so they have to fit
 into a standard c++ unsigned (!) int *)
contributingDiagramTypes = {
    OneLoopDiagram[0]
};

(* Find all diagrams of the type type_, testing all corresponding combinations of particles *)
(* IMPORTANT: Return value should have the format
 {{edmParticle1, {Diagram[DIAGRAMTYPENAME[_Integer], Particles___], Diagram[...], ...}},
  {edmParticle2, {...}},
  ...} *)

ContributingDiagramsOfType[OneLoopDiagram[0]] :=
    Module[{edmParticle = #, diagrams = SARAH`InsFields[
                           {{C[#, SARAH`AntiField[SARAH`FieldToInsert[1]],
                               SARAH`AntiField[SARAH`FieldToInsert[2]]],
                             C[SARAH`FieldToInsert[1], GetPhoton[],
                               SARAH`AntiField[SARAH`FieldToInsert[1]]],
                             C[SARAH`FieldToInsert[1], SARAH`FieldToInsert[2],
                               SARAH`AntiField[#]]},
                            {SARAH`FieldToInsert[1], SARAH`FieldToInsert[2]}}]},
           {edmParticle, DeleteDuplicates[(Module[{photonEmitter = #[[2,1]],
                                                   exchangeParticle = #[[2,2]]},
                  Diagram[OneLoopDiagram[0], edmParticle, photonEmitter, exchangeParticle]]
           &) /@ diagrams]}] & /@ edmParticles;

(* Returns the necessary c++ code corresponding to the vertices that need to be calculated.
 The returned value is a list {prototypes, definitions}. *)
createdVertices = Null;
CreateVertices[vertexRules_List] :=
    Module[{contributingDiagrams, vertices,
            vertexClassesPrototypes, vertexClassesDefinitions},
           If[createdVertices =!= Null, Return[createdVertices]];
           contributingDiagrams = ContributingDiagrams[];

           vertices = DeleteDuplicates @ Flatten[VerticesForDiagram /@
                                                 Flatten @ contributingDiagrams[[All, 2]], 1];

           {vertexClassesPrototypes, vertexClassesDefinitions} = Transpose @
               ((CreateVertexFunction[#, vertexRules] &) /@ vertices);
           vertexClassesPrototypes = Cases[vertexClassesPrototypes, Except[""]];
           vertexClassesDefinitions = Cases[vertexClassesDefinitions, Except[""]];

           createdVertices = {vertexClassesPrototypes, vertexClassesDefinitions};
           createdVertices = (StringJoin @ Riffle[#, "\n\n"] &) /@ createdVertices;

           Return[createdVertices];
          ];

(* Returns the vertices that are present in the specified diagram.
 This function should be overloaded for future diagram types.
 IMPORTANT: Lorentz indices have to be stripped away (They are unnecessary anyway) *)
VerticesForDiagram[Diagram[loopDiagram_OneLoopDiagram, edmParticle_, photonEmitter_, exchangeParticle_]] :=
    Module[{edmVertex1, photonVertex, edmVertex2},
           edmVertex1 = CachedVertex[{edmParticle, AntiField[photonEmitter], AntiField[exchangeParticle]}];
           photonVertex = CachedVertex[{photonEmitter, GetPhoton[], AntiField[photonEmitter]}];
           edmVertex2 = CachedVertex[{photonEmitter, exchangeParticle, AntiField[edmParticle]}];
           
           edmVertex1 = StripLorentzIndices @ edmVertex1[[1]];
           photonVertex = StripLorentzIndices @ photonVertex[[1]];
           edmVertex2 = StripLorentzIndices @ edmVertex2[[1]];

           Return[{edmVertex1, photonVertex, edmVertex2}];
           ];

(* Returns the vertex type for a vertex with a given list of particles *)
VertexTypeForParticles[particles_List] :=
    Module[{strippedParticles, scalars, vectors, fermions, scalarCount, vectorCount, fermionCount},
           strippedParticles = Vertices`StripFieldIndices /@ particles;

           scalars = Select[strippedParticles, (TreeMasses`IsScalar[#] || TreeMasses`IsScalar[AntiField[#]] &)];
           vectors = Select[strippedParticles, (TreeMasses`IsVector[#] || TreeMasses`IsVector[AntiField[#]] &)];
           fermions = Select[strippedParticles, (TreeMasses`IsFermion[#] || TreeMasses`IsFermion[AntiField[#]] &)];

           scalarCount = Length[scalars];
           vectorCount = Length[vectors];
           fermionCount = Length[fermions];

           If[fermionCount === 2 && scalarCount === 1 && vectorCount === 0,
              Return[LeftAndRightComponentedVertex]];
           If[fermionCount === 2 && scalarCount === 0 && vectorCount === 1,
              If[fermions[[1]] === AntiField[fermions[[2]]],
                 Return[LeftAndRightComponentedVertex]]];
           If[fermionCount === 0 && scalarCount === 2 && vectorCount === 1,
              Return[SingleComponentedVertex]];

           Return[Null];
           ];

(* Returns the different SARAH`Cp coupling parts for a vertex with a given list of particles *)
CouplingsForParticles[particles_List] :=
    Module[{vertexType, couplings},
           vertexType = VertexTypeForParticles[particles];
           couplings = {ReplacePart[particles, 0 -> SARAH`Cp]};

           couplings = Switch[vertexType,
                              SingleComponentedVertex, couplings,
                              LeftAndRightComponentedVertex, {couplings[[1]][SARAH`PL], couplings[[1]][SARAH`PR]}];
           Return[couplings];
           ];

(* Creates the actual c++ code for a vertex with given particles.
 This involves creating the VertexFunctionData<> code as well as
 the VertexFunction<> code. You should never need to change this code! *)
vertexFunctions = {};
CreateVertexFunction[indexedParticles_List, vertexRules_List] :=
    Module[{prototypes, definitions = "", ordering, particles, orderedParticles,
             orderedIndexedParticles, addSpacing = True},
            particles = Vertices`StripFieldIndices /@ indexedParticles;
            If[MemberQ[vertexFunctions, particles], Return[{"",""}]];

            ordering = Ordering[particles];
            orderedParticles = particles[[ordering]];
            orderedIndexedParticles = OrderParticles[indexedParticles, ordering];

            If[MemberQ[vertexFunctions, orderedParticles] === True,
               (* There is already an entry *)
               prototypes = "";
               addSpacing = False,
               (* There is no entry yet, create it *)
               {prototypes, definitions} = CreateOrderedVertexFunction[orderedIndexedParticles, vertexRules];
               AppendTo[vertexFunctions, orderedParticles];
               ];

            If[ordering === Table[i, {i, 1, Length[ordering]}],
               Return[{prototypes, definitions}]];

            orderedVertexFunction = ("VertexFunction<" <>
                                     StringJoin @ Riffle[ParticleToCXXName /@ orderedParticles, ", "] <>
                                     ">");

            prototypes = (prototypes <> If[addSpacing, "\n\n", ""] <>
                          "template<> struct VertexFunctionData<" <>
                          StringJoin @ Riffle[ParticleToCXXName /@ particles, ", "] <>
                          ">\n" <>
                          "{\n" <>
                          IndentText @
                          ("static const bool is_permutation = true;\n" <>
                           "typedef " <> orderedVertexFunction <> " orig_type;\n" <>
                           "typedef boost::mpl::vector_c<unsigned, " <>
                           StringJoin @ Riffle[ToString /@ (Ordering[ordering] - 1), ", "] <>
                           "> particlePermutation;\n"
                           ) <>
                          "};");

            AppendTo[vertexFunctions, particles];
            Return[{prototypes, definitions}];
          ];

(* Creates local declarations of field indices, whose values are taken
   from the elements of `arrayName'.
 *)
DeclareIndices[indexedParticles_List, arrayName_String] :=
    Module[{p, total = 0, fieldIndexList, decl = ""},
           DeclareIndex[idx_, num_Integer, an_String] := (
               "const unsigned " <> CConversion`ToValidCSymbolString[idx] <>
               " = " <> an <> "[" <> ToString[num] <> "];\n");
           For[p = 1, p <= Length[indexedParticles], p++,
               fieldIndexList = FieldIndexList[indexedParticles[[p]]];
               decl = decl <> StringJoin[DeclareIndex[#, total++, arrayName]& /@ fieldIndexList];
              ];
           Assert[total == Total[Length[FieldIndexList[#]]& /@ indexedParticles]];
           decl
          ];

(* ParsedVertex structure:
 ParsedVertex[
              {numP1Indices, numP2Indices, ...},
              {{minIndex1, minIndex2, ...}, {maxIndex1+1, maxIndex2+1, ...}},
              VertexClassName,
              VertexFunctionBody
              ]

 Getters are available! Given below ParseVertex[]
 *)

(* The heart of the algorithm! From the particle content, determine all
 necessary information. *)
ParseVertex[indexedParticles_List, vertexRules_List] :=
    Module[{particles, numberOfIndices, declareIndices,
        parsedVertex, vertexClassName, vertexFunctionBody,
        sarahParticles, particleInfo, indexBounds, expr, exprL, exprR},
           numberOfIndices = ((Length @ FieldIndexList[#] &) /@ indexedParticles);
           particles = Vertices`StripFieldIndices /@ indexedParticles;
           declareIndices = DeclareIndices[indexedParticles, "indices"];

           vertexClassName = SymbolName[VertexTypeForParticles[particles]];
           vertexFunctionBody = Switch[vertexClassName,
                                       "SingleComponentedVertex",
                                       expr = (SARAH`Cp @@ indexedParticles) /. vertexRules;
                                       expr = TreeMasses`ReplaceDependenciesReverse[expr];
                                       "std::complex<double> result;\n\n" <>
                                       declareIndices <>
                                       Parameters`CreateLocalConstRefs[expr] <> "\n" <>
                                       TreeMasses`ExpressionToString[expr, "result"] <> "\n" <>
                                       "return vertex_type(result);",

                                       "LeftAndRightComponentedVertex",
                                       exprL = SARAH`Cp[Sequence @@ indexedParticles][SARAH`PL] /. vertexRules;
                                       exprR = SARAH`Cp[Sequence @@ indexedParticles][SARAH`PR] /. vertexRules;
                                       exprL = TreeMasses`ReplaceDependenciesReverse[exprL];
                                       exprR = TreeMasses`ReplaceDependenciesReverse[exprR];
                                       "std::complex<double> left, right;\n\n" <>
                                       declareIndices <>
                                       Parameters`CreateLocalConstRefs[exprL + exprR] <> "\n" <>
                                       TreeMasses`ExpressionToString[exprL, "left"] <> "\n" <>
                                       TreeMasses`ExpressionToString[exprR, "right"] <> "\n" <>
                                       "return vertex_type(left, right);"];

           sarahParticles = SARAH`getParticleName /@ particles;
           particleInfo = Flatten[(Cases[SARAH`Particles[FlexibleSUSY`FSEigenstates], {#, ___}] &) /@
                                  sarahParticles, 1];

           (* INFO: I do not think this ever occurs... *)
           particleInfo = DeleteCases[particleInfo, {SARAH`generation, 1}, {3}];
           particleInfo = DeleteCases[particleInfo, {SARAH`lorentz, _}, {3}];

           indexBounds = (With[{particleIndex = #},
                               (If[#[[1]] === SARAH`generation,
                                   {particleInfo[[particleIndex, 2]]-1, particleInfo[[particleIndex, 3]]},
                                   {1, #[[2]]}]
                                &) /@ particleInfo[[particleIndex, 5]]]
                          &) /@ Table[i, {i, Length[particles]}];
           indexBounds = Cases[Flatten[indexBounds, 1], Except[{}]];

           If[indexBounds === {},
              indexBounds = {{},{}},
              indexBounds = Transpose @ indexBounds];

           parsedVertex = ParsedVertex[numberOfIndices,
                                       indexBounds,
                                       vertexClassName,
                                       vertexFunctionBody];

           Return[parsedVertex];
           ];

(** Getters to the ParsedVertex structure **)
NumberOfIndices[parsedVertex_ParsedVertex] := Total[parsedVertex[[1]]];
NumberOfIndices[parsedVertex_ParsedVertex, pIndex_Integer] := parsedVertex[[1, pIndex]];

IndexBounds[parsedVertex_ParsedVertex] := parsedVertex[[2]];

VertexClassName[parsedVertex_ParsedVertex] := parsedVertex[[3]];
VertexFunctionBody[parsedVertex_ParsedVertex] := parsedVertex[[4]];
(** End getters **)

(* Create the c++ code for a canonically ordered vertex *)
CreateOrderedVertexFunction[orderedIndexedParticles_List, vertexRules_List] :=
    Module[{prototype, definition, orderedParticles, dataClassName, functionClassName,
            parsedVertex, particleIndexStartF, particleIndexStart, indexBounds},
            orderedParticles = Vertices`StripFieldIndices /@ orderedIndexedParticles;
            parsedVertex = ParseVertex[orderedIndexedParticles, vertexRules];
            dataClassName = "VertexFunctionData<" <> StringJoin @ Riffle[ParticleToCXXName /@ orderedParticles, ", "] <> ">";
            functionClassName = "VertexFunction<" <> StringJoin @ Riffle[ParticleToCXXName /@ orderedParticles, ", "] <> ">";

            particleIndexStartF[1] = 0;
            particleIndexStartF[pIndex_] := particleIndexStartF[pIndex-1] + NumberOfIndices[parsedVertex, pIndex-1];
            particleIndexStartF[Length[orderedParticles]+1] = NumberOfIndices[parsedVertex];

            particleIndexStart = Table[particleIndexStartF[i], {i, 1, Length[orderedParticles] + 1}];

            prototype = ("template<> struct " <> dataClassName <> "\n" <>
                         "{\n" <>
                         IndentText @
                         ("static const bool is_permutation = false;\n" <>
                          "typedef IndexBounds<" <> ToString @ NumberOfIndices[parsedVertex] <> "> index_bounds;\n" <>
                          "typedef " <> VertexClassName[parsedVertex] <> " vertex_type;\n" <>
                          "typedef boost::mpl::vector_c<unsigned, " <>
                             StringJoin @ Riffle[ToString /@ particleIndexStart, ", "] <>
                          "> particleIndexStart;\n" <>
                          "static const index_bounds indexB;\n"
                          ) <>
                         "};");

            indexBounds = IndexBounds[parsedVertex];

            If[NumberOfIndices[parsedVertex] =!= 0,
               prototype = (prototype <> "\n" <>
                            "const " <> dataClassName <> "::index_bounds " <> dataClassName <> "::indexB = { " <>
                            "{ " <> StringJoin @ Riffle[ToString /@ indexBounds[[1]], ", "] <> " }, " <>
                            "{ " <> StringJoin @ Riffle[ToString /@ indexBounds[[2]], ", "] <> " } };"
                            );];
            definition = ("template<> template<> " <> functionClassName <> "::vertex_type\n" <>
                          functionClassName <> "::vertex(const indices_type &indices, EvaluationContext &context)\n" <>
                          "{\n" <>
                          IndentText @ VertexFunctionBody[parsedVertex] <> "\n" <>
                          "}");

            Return[{prototype, definition}];
            ];

(* Find all contributing diagrams *)
cachedContributingDiagrams = Null;
ContributingDiagrams[] :=
       Module[{diagrams},
           If[cachedContributingDiagrams =!= Null, Return[cachedContributingDiagrams]];

           cachedContributingDiagrams = Flatten[(ContributingDiagramsOfType[#] &)
                                                /@ contributingDiagramTypes
                                                , 1];
           cachedContributingDiagrams = ({#, Union @
                   ReplacePart[Cases[cachedContributingDiagrams,
                         {#, diagrams_List} -> diagrams], 0 -> Sequence]} &) /@ edmParticles;
           
           Return[cachedContributingDiagrams];
          ];

(* Returns a list of all concrete diagram evaluators
 format: {{edmParticle1, {"DiagramEvaluator<OneLoopDiagram<1>, Fe, VP>", "...", ... }},
          {edmParticle2, {"...", ... }},
          ...}
 that need to be invoked in our calculation *)
ConcreteDiagramEvaluators[] :=
     ({#[[1]],
         (("DiagramEvaluator<" <> SymbolName @ Head @ #[[1]] <> "<" <>
           ToString @ #[[1,1]] <> ">, " <>
           StringJoin @ (Riffle[ParticleToCXXName /@ ReplacePart[#[[2;;]], 0 -> List], ", "]) <>
           ">" &)
          /@ #[[2]]) } &) /@ ContributingDiagrams[];

(* TODO: End[]; *)

EndPackage[];

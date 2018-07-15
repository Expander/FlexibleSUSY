(* ::Package:: *)

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

BeginPackage["FFMassiveVFormFactors`",
    {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "CXXDiagrams`", "ColorMathInterface`"}
];

FFMassiveVFormFactorsCreateInterface::usage = "";
f::usage = "";
MassiveVIndices::usage = "";

Begin["Private`"];

MassiveVIndices[V_] :=
   Module[{name = CXXNameOfField[V], dim = TreeMasses`GetDimension[V]},
      (* we expect that the gauge group from which massive vector boson come from is broke *)
      Assert[dim === 1];

      "template<>\n" <>
      "typename field_indices<" <> name <> ">::type\n" <>
      "default_indices_for_spectator<" <> name <> "> (void)\n" <>
      "{\n" <>
         "return {};\n" <>
      "}\n"
   ];

FFMassiveVFormFactorsCreateInterface[inFermion_, outFermion_, spectator_, loopParticles_List] :=
    Module[{prototype, definition,
            numberOfIndices1 = CXXDiagrams`NumberOfFieldIndices[inFermion],
            numberOfIndices2 = CXXDiagrams`NumberOfFieldIndices[outFermion],
            numberOfIndices3 = CXXDiagrams`NumberOfFieldIndices[spectator]},

      prototype =
         "std::valarray<std::complex<double>> calculate_" <> CXXNameOfField[inFermion] <>
            "_" <> CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[spectator] <> "_form_factors" <>
            " (\n" <>
            If[TreeMasses`GetDimension[inFermion] =!= 1,
               "   int generationIndex1, ",
               " "
            ] <>
            If[TreeMasses`GetDimension[outFermion] =!= 1,
               " int generationIndex2, ",
               " "
            ] <>
            "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model )";

      definition =
          prototype <> "{\n" <>
            IndentText[
               FlexibleSUSY`FSModelName <> "_mass_eigenstates model_ = model;\n" <>
               "EvaluationContext context{ model_ };\n" <>
               "std::array<int, " <> ToString @ numberOfIndices1 <> "> indices1 = {" <>
                     (* TODO: Specify indices correctly *)
                       If[TreeMasses`GetDimension[inFermion] =!= 1,
                          " generationIndex1" <>
                          If[numberOfIndices1 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices1-1}],
                             ""] <> " ",
                          If[numberOfIndices1 =!= 0,
                             StringJoin @ Riffle[Table[" 0", {numberOfIndices1}], ","] <> " ",
                             ""]
                         ] <> "};\n" <>
                   "std::array<int, " <> ToString @ numberOfIndices2 <>
                     "> indices2 = {" <>
                       If[TreeMasses`GetDimension[outFermion] =!= 1,
                          " generationIndex2" <>
                          If[numberOfIndices2 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices2-1}],
                             ""] <> " ",
                          If[numberOfIndices2 =!= 0,
                             StringJoin @ Riffle[Table[" 0", {numberOfIndices2}], ","] <> " ",
                             ""]
                         ] <> "};\n\n" <>

               "std::valarray<std::complex<double>> val {0.0, 0.0};\n\n" <>

               StringJoin[
                  ("val += std::complex<double> {" <> (ToString @ N[#[[2,1]], 16]) <> "} * FFMassiveVVertexCorrectionFS<" <>
                   StringRiffle[CXXDiagrams`CXXNameOfField /@ {inFermion, outFermion, spectator, #[[1,1]], #[[1,2]]}, ","]  <>
                   ">::value(indices1, indices2, context);\n") & /@ loopParticles
               ] <> "\n" <>

               "return val;"

            ] <> "\n}\n\n";

    {prototype <> ";", definition}
  ];

(* if t is TreeMasses`IsScalar then returns list of scalars and anti-scalar, etc. *)
getParticlesOfType[t_] :=
    DeleteDuplicates@Join[#, CXXDiagrams`LorentzConjugate /@ #]& @
      Select[TreeMasses`GetParticles[], t];

vertexNonZero[vertex_] :=
    Transpose[Drop[vertex, 1]][[1]] =!= {0,0};

(* if a diagram exists, return a color factor and a list of particles in vertices, otherwise return an empty list *)
singleDiagram[inFermion_, outFermion_, spectator_, F_?TreeMasses`IsFermion, S_?TreeMasses`IsScalar] :=
   Module[{FBarFjSBar, FiBarFS, SBarSVBar, FBarFVBar, v1, v2, v3, v4,colorIndexAssociation},

      (* if the electric charge of an incomind particle doesn't equal to the sum of charges of outgoing ones,
         return an {} *)
      If[TreeMasses`GetElectricCharge[inFermion] =!= Plus @@ (TreeMasses`GetElectricCharge /@ {S,F}),
         Return[{}]
      ];

      colorIndexAssociation = Association[{}];
      If[TreeMasses`ColorChargedQ[#],
            If[TreeMasses`GetDimension[#] === 1,
               AssociateTo[colorIndexAssociation, # -> {Unique["ct"]}],
               AssociateTo[colorIndexAssociation, # -> {Unique["gt"], Unique["ct"]}]
         ]
      ]& /@
          {inFermion, CXXDiagrams`LorentzConjugate[outFermion], CXXDiagrams`LorentzConjugate[spectator], F, S};

      (*Assert[ Keys[colorIndexAssociation]]*)

      Print["========================="];
      (*Print[colorIndexAssociation];*)
      (*Print[AddIndices[F, colorIndexAssociation], " ", AddIndices[CXXDiagrams`LorentzConjugate[F], colorIndexAssociation]];*)

      v1 = {CXXDiagrams`LorentzConjugate[F], inFermion, CXXDiagrams`LorentzConjugate[S]};
      FBarFjSBar = SARAH`Vertex[AddIndices[#, colorIndexAssociation]& /@ v1];
      v2 = {CXXDiagrams`LorentzConjugate[outFermion], F, S};
      FiBarFS = SARAH`Vertex[AddIndices[#, colorIndexAssociation]& /@ v2];
      v3 = {CXXDiagrams`LorentzConjugate[S], S, CXXDiagrams`LorentzConjugate[spectator]};
      SBarSVBar = SARAH`Vertex[AddIndices[#, colorIndexAssociation]& /@ v3];
      v4 = {CXXDiagrams`LorentzConjugate[F], F, CXXDiagrams`LorentzConjugate[spectator]};
      FBarFVBar = SARAH`Vertex[AddIndices[#, colorIndexAssociation]& /@ v4];

      If[
        vertexNonZero[FBarFjSBar]
            && vertexNonZero[FiBarFS]
            && (vertexNonZero[SBarSVBar] || vertexNonZero[FBarFVBar]),
         Print[ColorMathInterface`DropColorles[{FBarFjSBar, FiBarFS, SBarSVBar, FBarFVBar}]];
         Print[{FBarFjSBar, FiBarFS, SBarSVBar, FBarFVBar} //
             ColorMathInterface`DropColorles (*//
             ColorMathInterface`TakeOnlyColor*)
         ];
        {1, {v1, v2, v3, v4}},
        {}
      ]
   ];

AddIndices[field_, ass_] :=
    Module[{temp, kupa},
    temp = Lookup[ass, field, {}];
    (*Print[field, " ", temp];*)
    kupa = If[temp === {},
       field,
       field /. {h_Symbol[f_Symbol] :> h[f[temp]],
          f_Symbol /; (f =!= bar && f =!= conj && f =!= List) :> f[temp]}
    ];
       (*Print["kupa", kupa];*)
       kupa
    ];

f[inFermion_, outFermion_, spectator_] :=
   Module[{scalars, fermions, internalParticles = {}, temp},

      scalars = getParticlesOfType[TreeMasses`IsScalar];
      fermions = getParticlesOfType[TreeMasses`IsFermion];

      Map[
         (temp = singleDiagram[inFermion, outFermion, spectator, #[[1]], #[[2]]];
         If[temp =!= {},
            AppendTo[internalParticles, {#, temp}];
            ])&,
         Tuples[{fermions, scalars}]
      ];

      internalParticles
   ];

(* evaluate single diagram *)
(*
CXXEvaluatorsForLeptonPairAndDiagramFromGraph[inFermion_, outFermion_, spectator_, diagram_, vertexCorrectionGraph] := 
    Module[{Emitter, exchangeParticle, colorFactor, colorFactorStr},

        Emitter = CXXDiagrams`LorentzConjugate[diagram[[4,3]]]; (* Edge between vertices 4 and 6 (3rd edge of vertex 4) *)
        exchangeParticle = diagram[[4,2]]; (* Edge between vertices 4 and 5 (2nd edge of vertex 4) *)
    
        colorFactor = getChargeFactor[
 {
  {
   Cp[inFermion, exchangeParticle, AntiField[Emitter]],
   Cp[spectator, Emitter, AntiField[Emitter]],
   Cp[AntiField[outFermion], AntiField[exchangeParticle], 
    Emitter]
   },
  {
   External[1] -> inFermion, External[2] -> AntiField[outFermion], 
   External[3] -> spectator,
   Internal[1] -> Emitter, Internal[2] -> exchangeParticle, 
   Internal[3] -> AntiField[Emitter]
   }
  },
 {
  {{inFermion, ex1}, {exchangeParticle, 
    in2}, {AntiField[Emitter], in1}},
  {{spectator, ex3}, {Emitter, in3}, {AntiField[Emitter], in1}},
  {{AntiField[outFermion], ex2}, {AntiField[exchangeParticle], 
    in2}, {Emitter, in3}}
  }
 ];

        colorFactorStr = "std::complex<double> " <>
            ToString @ (N[#, 16]& /@ (ReIm[colorFactor]/EvaluateColorStruct[Emitter, exchangeParticle]));

        colorFactorStr <> " * " <> CXXEvaluator[{inFermion, outFermion, spectator}, {Emitter, exchangeParticle}]
    ];
*)

(* TODO: add other topologies? *)

End[];
EndPackage[];

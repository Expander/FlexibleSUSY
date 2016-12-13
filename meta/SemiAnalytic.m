
BeginPackage["SemiAnalytic`", {"SARAH`", "CConversion`", "Constraint`", "Parameters`",
                               "TextFormatting`", "WriteOut`"}];

Begin["`Private`"];

IsDimensionOne[par_] :=
    Module[{dimOnePars},
           dimOnePars = { SARAH`BetaTijk };
           If[SARAH`SupersymmetricModel,
              dimOnePars = Append[dimOnePars, SARAH`BetaMi];,
              dimOnePars = Append[dimOnePars, SARAH`BetaMuij];
             ];
           dimOnePars = (Parameters`StripIndices[#[[1]]])& /@ (Join @@ dimOnePars);
           MemberQ[dimOnePars, Parameters`StripIndices[par]]
          ];

IsDiracGauginoMass[par_] :=
    Module[{diracMasses = {}},
           If[SARAH`SupersymmetricModel,
              diracMasses = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaDGi;
             ];
           MemberQ[diracMasses, Parameters`StripIndices[par]]
          ];

IsScalarMass[par_] :=
    Module[{scalarMasses},
           If[SARAH`SupersymmetricModel,
              scalarMasses = Parameters`StripIndices[#[[1]]]& /@ SARAH`Betam2ij;,
              scalarMasses = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaBij;
             ];
           MemberQ[scalarMasses, Parameters`StripIndices[par]]
          ];

IsSoftBilinear[par_] :=
    Module[{softBilinears = {}},
           If[SARAH`SupersymmetricModel,
              softBilinears = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaBij;
             ];
           MemberQ[softBilinears, Parameters`StripIndices[par]]
          ];

IsSoftLinear[par_] :=
    Module[{softLinears = {}},
           If[SARAH`SupersymmetricModel,
              softLinears = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaLSi;
             ];
           MemberQ[softLinears, Parameters`StripIndices[par]]
          ];

IsAllowedSemiAnalyticParameter[par_] :=
    Or[IsDimensionOne[par],
       IsDiracGauginoMass[par],
       IsScalarMass[par],
       IsSoftBilinear[par],
       IsSoftLinear[par]];

End[];

EndPackage[];

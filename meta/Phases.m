
BeginPackage["Phases`", {"TextFormatting`", "CConversion`"}];

ConvertSarahPhases::usage="convert SARAH phases to our format";

CreatePhasesDefinition::usage="creates definitions of the phases";

CreatePhasesGetterSetters::usage="creates function definitions
for phase getters and setters";

CreatePhasesInitialization::usage="creates initialization list of
phases"

CreatePhaseName::usage = "creates the name of a SARAH phase";

GetArg::usage = "returs the real arg of a given phase";

Begin["`Private`"];

GetArg[Exp[phase_]] := GetArg[phase];
GetArg[I phase_] := phase;
GetArg[phase_] := phase;

CreatePhaseName[phase_] :=
    CConversion`ToValidCSymbolString[GetArg[phase]];

ConvertSarahPhases[phases_List] :=
    DeleteDuplicates[(#[[2]])& /@ phases];

GetPhaseType[Exp[_]] := CConversion`ScalarType[CConversion`realScalarCType];
GetPhaseType[_]      := CConversion`ScalarType[CConversion`complexScalarCType];

CreatePhaseInit[type_ /; type === CConversion`ScalarType[CConversion`realScalarCType]] :=
    "(0)";

CreatePhaseInit[type_ /; type === CConversion`ScalarType[CConversion`complexScalarCType]] :=
    "(1,0)";

CreatePhasesDefinition[phases_List] :=
    Module[{result = "", k},
           For[k = 1, k <= Length[phases], k++,
               result = result <> CConversion`CreateCType[
                            GetPhaseType[phases[[k]]]] <> " " <>
                        CreatePhaseName[phases[[k]]] <> ";\n";
              ];
           Return[result];
          ];

CreatePhasesGetterSetters[phases_List] :=
    Module[{result = "", k, phaseName, type},
           For[k = 1, k <= Length[phases], k++,
               phaseName = CreatePhaseName[phases[[k]]];
               type = GetPhaseType[phases[[k]]];
               result = result <>
                        CConversion`CreateInlineSetter[phaseName, type] <>
                        CConversion`CreateInlineGetter[phaseName, type];
              ];
           Return[result];
          ];

CreatePhasesInitialization[phases_List] :=
    Module[{result = "", k, phaseName, type},
           For[k = 1, k <= Length[phases], k++,
               phaseName = CreatePhaseName[phases[[k]]];
               type = GetPhaseType[phases[[k]]];
               result = result <> ", " <> phaseName <> CreatePhaseInit[type];
              ];
           Return[result];
          ];

End[];

EndPackage[];

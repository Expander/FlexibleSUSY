BeginPackage["FlexibleTower`", {"CConversion`"}];

GetModelTypes::usage = "Returns comma-separated list of model types";

Begin["`Private`"];

GetModelType[name_String, templatePar_String] :=
    name <> "<" <> templatePar <> ">";

GetSLHAModelType[name_String, templatePar_String] :=
    GetModelType[name <> "_slha", GetModelType[name, templatePar]];

GetModelTypes[] :=
    Utils`StringJoinWithSeparator[
        Which[FlexibleSUSY`FlexibleEFTHiggs === True,
              { GetSLHAModelType[FlexibleSUSY`FSModelName, "Solver_type"],
                GetModelType["standard_model::StandardModel", "Solver_type"] },
              True,
              { GetSLHAModelType[FlexibleSUSY`FSModelName, "Solver_type"] }
             ],
        ", "
    ];

End[];

EndPackage[];

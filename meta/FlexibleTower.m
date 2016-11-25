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
              { GetSLHAModelType[FlexibleSUSY`FSModelName, "algorithm_type"],
                GetModelType["standard_model::StandardModel", "algorithm_type"] },
              True,
              { GetSLHAModelType[FlexibleSUSY`FSModelName, "algorithm_type"] }
             ],
        ", "
    ];

End[];

EndPackage[];

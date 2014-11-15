$flexiblesusyOutputDir   = "models/MSSMNoFVatMGUT";
$flexiblesusyConfigDir   = FileNameJoin[{Directory[], "config"}];
$flexiblesusyMetaDir     = FileNameJoin[{Directory[], "meta"}];
$flexiblesusyTemplateDir = FileNameJoin[{Directory[], "templates"}];
AppendTo[$Path, $flexiblesusyMetaDir];

Needs["SARAH`"];
Needs["FlexibleSUSY`", FileNameJoin[{$flexiblesusyMetaDir, "FlexibleSUSY.m"}]];

workingDirectory = Directory[];
SARAH`SARAH[OutputDirectory] = FileNameJoin[{workingDirectory, "Output"}];
SARAH`SARAH[InputDirectories] = {
    ToFileName[{$sarahDir, "Models"}],
    FileNameJoin[{workingDirectory, "sarah"}]
};

Print["Current working directory: ", workingDirectory];
Print["SARAH output directory: ", SARAH`SARAH[OutputDirectory]];
Print["Starting model: ", ToString[HoldForm[Start["MSSMNoFV"]]]];

Start["MSSMNoFV"];

MakeFlexibleSUSY[InputFile -> FileNameJoin[{"models/MSSMNoFVatMGUT","FlexibleSUSY.m"}]];

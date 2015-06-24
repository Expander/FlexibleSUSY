(* ::Package:: *)

$flexiblesusyDir         = "/home/bachi/Programme/FlexibleSUSY";
$flexiblesusyOutputDir   = FileNameJoin[{$flexiblesusyDir, "markus_tests"}];
$flexiblesusyConfigDir   = FileNameJoin[{$flexiblesusyDir, "config"}];
$flexiblesusyMetaDir     = FileNameJoin[{$flexiblesusyDir, "meta"}];
$flexiblesusyTemplateDir = FileNameJoin[{$flexiblesusyDir, "templates"}];

AppendTo[$Path, $flexiblesusyMetaDir];

(* Needs["SARAH`"]; *)
Needs["FlexibleSUSY`", FileNameJoin[{$flexiblesusyMetaDir, "FlexibleSUSY.m"}]];

(* SARAH`SARAH[OutputDirectory] = FileNameJoin[{$flexiblesusyDir, "Output"}];
SARAH`SARAH[InputDirectories] = {
    FileNameJoin[{$flexiblesusyDir, "sarah"}],
    ToFileName[{$sarahDir, "Models"}]
}; *)

Print["Let's start!"];

FlexibleSUSY`FSModelName = "ModelName";

FlexibleSUSY`Private`WriteWeinbergAngleClass[
               {{FileNameJoin[{Global`$flexiblesusyTemplateDir, "weinberg_angle.hpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_weinberg_angle.hpp"}]},
                {FileNameJoin[{Global`$flexiblesusyTemplateDir, "weinberg_angle.cpp.in"}],
                 FileNameJoin[{Global`$flexiblesusyOutputDir, FlexibleSUSY`FSModelName <> "_weinberg_angle.cpp"}]}
               }
                                      ];

Print["Work done!"];

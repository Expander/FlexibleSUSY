(*
   This script loads SARAH and FlexibleSUSY.
   Usage:

     Get["utils/load-FlexibleSUSY.m"];
*)

AppendTo[$Path, FileNameJoin[{Directory[], "meta"}]];

Needs["SARAH`"];
Needs["FlexibleSUSY`", FileNameJoin[{Directory[], "meta", "FlexibleSUSY.m"}]];

SARAH`SARAH[OutputDirectory] = FileNameJoin[{Directory[], "Output"}];

SARAH`SARAH[InputDirectories] = {
    FileNameJoin[{Directory[], "sarah"}],
    ToFileName[{$sarahDir, "Models"}]
};

Print["Current working directory: ", Directory[]];
Print["SARAH output directory   : ", SARAH`SARAH[OutputDirectory]];
Print["SARAH models directories : ", SARAH`SARAH[InputDirectories]];

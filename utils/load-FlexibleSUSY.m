(*
   This script loads SARAH and FlexibleSUSY.
   Usage:

     Get["load-FlexibleSUSY.m"];
*)

AppendTo[$Path, FileNameJoin[{Directory[], "meta"}]];

Needs["SARAH`"];
Needs["FlexibleSUSY`", FileNameJoin[{Directory[], "meta", "FlexibleSUSY.m"}]];


BeginPackage["AnomalousDimension`", {"SARAH`", "TextFormatting`", "CConversion`", "TreeMasses`"}];

AnomalousDimension[];

CreateAnomDimFunctions::usage="";
CreateAnomDimPrototypes::usage="";
ConvertSarahAnomDim::usage="";

Begin["Private`"];

GetName[AnomalousDimension[name_, type_, anom_List]] := name;

GetType[AnomalousDimension[name_, type_, anom_List]] := type;

GetAnomDim1Loop[AnomalousDimension[name_, type_, anom_List]] := anom[[1]];

GetAnomDim2Loop[AnomalousDimension[name_, type_, anom_List]] := anom[[2]];

GetAllAnomDims[AnomalousDimension[name_, type_, anom_List]] := anom;

CreateValidAnomDimName[names_List /; Length[names] == 2] :=
    Module[{name1, name2},
           name1 = GetHead[names[[1]]];
           name2 = GetHead[names[[2]]];
           Return[Symbol[ToString[name1] <> ToString[name2]]];
          ];

CreateValidAnomDimName[names_] :=
    Module[{},
           Print["Error: don't know how to create name for anomalous ",
                 "dimension from symbol: ", names];
           Quit[];
          ];

(* Converts SARAH anomalous dimensions
 *
 * SARAH format:
 *   { {name1, name2}, one-loop, two-loop }
 *
 * Our format:
 *   AnomalousDimension[cName, type, {one-loop, two-loop}]
 *
 * @param gij list of SARAH-like formated anomalous dimensions
 *)
ConvertSarahAnomDim[gij_List] :=
    Module[{lst = {}, adim, i, name, type, dim},
           For[i = 1, i <= Length[gij], i++,
               adim = gij[[i]];
               (* adim[[1]] == {name1,name2}, adim[[2]] == 1-loop anom. dim *)
               name = CreateValidAnomDimName[adim[[1]]];
               If[FreeQ[adim[[2]], a_[i1,i2]],
                  type = CConversion`ScalarType["double"];,
                  dim = TreeMasses`GetDimension[adim[[1,1]]];
                  type = CConversion`MatrixType["Eigen::Matrix<double," <>
                                                ToString[dim] <> "," <>
                                                ToString[dim] <> ">", dim, dim];
                 ];
               AppendTo[lst, AnomalousDimension[name, type, Drop[adim, 1]]];
              ];
           Return[lst];
          ];

(* create anomalous dimension prototypes *)
CreateAnomDimPrototypes[anomDim_AnomalousDimension] :=
    Module[{prototypes, name, type},
           name = ToValidCSymbolString[GetName[anomDim]];
           type = GetType[anomDim];
           prototypes = CreateGetterPrototype[name, GetCParameterType[type]];
           Return[prototypes];
          ];

CreateAnomDimPrototypes[anomDim_List] :=
    Module[{prototypes = ""},
           (prototypes = prototypes <> CreateAnomDimPrototypes[#])& /@ anomDim;
           Return[prototypes];
          ];

CreateAnomDimFunction[anomDim_AnomalousDimension] :=
    Module[{def = "", body = "", type, name, unitMatrix},
           type = GetType[anomDim];
           name = ToValidCSymbolString[GetName[anomDim]];
           unitMatrix = CreateUnitMatrix[type];
           body = "const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);\n" <>
                  "const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;\n";
           body = body <> CreateDefaultDefinition["anomDim", type] <> ";\n";
           (* one-loop *)
           body = body <> "\nanomDim = " <>
                  RValueToCFormString[(CConversion`oneOver16PiSqr * GetAnomDim1Loop[anomDim])
                                /. { Kronecker[i1,i2] -> unitMatrix }
                                /. { a_[i1,i2] :> a }];
           body = body <> ";\n";
           (* two-loop *)
           If[Length[GetAllAnomDims[anomDim]] > 1,
              body = body <> "\nif (displayLoops() > 1) {\n";
              body = body <> "   anomDim += " <>
                     RValueToCFormString[(CConversion`twoLoop * GetAnomDim2Loop[anomDim])
                                   /. { Kronecker[i1,i2] -> unitMatrix }
                                   /. { a_[i1,i2] :> a }];
              body = body <> ";\n}\n";
             ];
           body = body <> "\nreturn anomDim;\n";
           body = IndentText[body];
           def  = GetCParameterType[type] <> " CLASSNAME::get_" <>
                  name <> "() const\n{\n" <> body <> "}\n\n";
           Return[def];
          ];

CreateAnomDimFunctions[anomDim_List] :=
    Module[{functions = ""},
           (functions = functions <> CreateAnomDimFunction[#])& /@ anomDim;
           Return[functions];
          ];

End[];

EndPackage[];

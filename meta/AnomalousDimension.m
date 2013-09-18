
BeginPackage["AnomalousDimension`", {"SARAH`", "TextFormatting`", "CConversion`", "TreeMasses`", "Parameters`"}];

AnomalousDimension[];

CreateAnomDimFunctions::usage="";
CreateAnomDimPrototypes::usage="";
ConvertSarahAnomDim::usage="";

Begin["`Private`"];

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
           Quit[1];
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
               If[FreeQ[adim[[2]], a_[Susyno`LieGroups`i1,SARAH`i2]],
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
    Module[{def, body, type, name, unitMatrix,
            exprOneLoop, exprTwoLoop, inputParsDecl},
           type = GetType[anomDim];
           name = ToValidCSymbolString[GetName[anomDim]];
           unitMatrix = CreateUnitMatrix[type];
           (* one-loop *)
           exprOneLoop = (CConversion`oneOver16PiSqr * GetAnomDim1Loop[anomDim]) /.
                         { Kronecker[Susyno`LieGroups`i1,SARAH`i2] -> unitMatrix,
                           a_[Susyno`LieGroups`i1,SARAH`i2] :> a };
           body = "\nanomDim = " <> RValueToCFormString[exprOneLoop] <> ";\n";
           (* two-loop *)
           If[Length[GetAllAnomDims[anomDim]] > 1,
              exprTwoLoop = (CConversion`twoLoop * GetAnomDim2Loop[anomDim]) /.
                            { Kronecker[Susyno`LieGroups`i1,SARAH`i2] -> unitMatrix,
                              a_[Susyno`LieGroups`i1,SARAH`i2] :> a };
              body = body <> "\nif (get_loops() > 1) {\n" <>
                     IndentText["anomDim += " <> RValueToCFormString[exprTwoLoop]] <>
                     ";\n}\n";
             ];
           inputParsDecl = Parameters`CreateLocalConstRefsForInputParameters[exprOneLoop + exprTwoLoop];
           body = "const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;\n" <>
                  CreateDefaultDefinition["anomDim", type] <> ";\n" <>
                  inputParsDecl <>
                  body <>
                  "\nreturn anomDim;\n";
           def  = GetCParameterType[type] <> " CLASSNAME::get_" <>
                  name <> "() const\n{\n" <> IndentText[body] <> "}\n\n";
           Return[def];
          ];

CreateAnomDimFunctions[anomDim_List] :=
    Module[{functions = ""},
           (functions = functions <> CreateAnomDimFunction[#])& /@ anomDim;
           Return[functions];
          ];

End[];

EndPackage[];

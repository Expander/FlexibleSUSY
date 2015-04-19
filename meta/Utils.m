
BeginPackage["Utils`"];

StringJoinWithSeparator::usage = "Joins a list of strings with a given separator string";

Zip::usage = "Combines two lists to a list of touples.
Example:

In[]:= Zip[{a,b,c},{d,e,f}]
Out[]= {{a, d}, {b, e}, {c, f}}
";

StringZip::usage = "Combines lists of strings to a list of concated
strings.  Example:

In[]:= StringZip[{\"a\",\"b\"},{\"d\",\"e\"}]
Out[]= {\"ad\", \"be\"}}
";

StringZipWithSeparator::usage = "Combines lists of strings to a list
of concated strings using a separator.  Example:

In[]:= StringZipWithSeparator[{\"a\",\"b\"},{\"d\",\"e\"}, \"_\"]
Out[]= {\"a_d\", \"b_e\"}}
";

FSGetOption::usage = "Returns the value of an option from a list of
options.

Example:

In[]:= opts = {o1 -> x, o2 -> y}

In[]:= FSGetOption[opts, o1]
Out[]= x

In[]:= FSGetOption[opts, o2]
Out[]= y

In[]:= FSGetOption[opts, o3]
Error: option o3 not found
";

FSSetOption::usage = "Returns given list of options, where the
occurrence of the given rule is replaced (if it exists) or added (if
it does not exist).

Example:

In[]:= opts = {o1 -> x, o2 -> y}

In[]:= FSSetOption[opts, o1 -> z]
Out[]= {o2 -> y, o1 -> z}

In[]:= FSSetOption[opts, o3 -> z]
Out[]= {o1 -> x, o2 -> y, o3 -> z}
";

FSImportString::usage = "Returns the content of a file in form of a string.  If the file does not exist, \"unknown\" is returned.";

Begin["`Private`"];

StringJoinWithSeparator[list_List, separator_String, transformer_:ToString] :=
    Module[{result = "", i},
           For[i = 1, i <= Length[list], i++,
               If[i > 1, result = result <> separator;];
               result = result <> transformer[list[[i]]];
              ];
           Return[result];
          ];

Zip[list1_List, list2_List] :=
    MapThread[List, {list1, list2}];

StringZip[lists___List] :=
    MapThread[StringJoin, {lists}];

StringZipWithSeparator[lists___List, separator_String] :=
    MapThread[StringJoinWithSeparator[{##},separator]&, {lists}];

FSGetOption[opts_List, opt_] :=
    Module[{values},
           values = Cases[opts, (Rule[opt, value_] | RuleDelayed[opt, value_]) :> value];
           Switch[Length[values],
                  0, Print["Error: option ", opt, " not found"];
                     Null,
                  1, values[[1]],
                  _, Print["Warning: option ", opt, " not unique"];
                     values[[1]]
                 ]
          ];

FSSetOption[opts_List, rule:Rule[opt_, value_]] :=
    Join[
        Cases[opts, Except[Rule[opt,_] | RuleDelayed[opt,_]]],
        {rule}
        ];

FSSetOption[opts_List, rule:RuleDelayed[opt_, value_]] :=
    Join[
        Cases[opts, Except[Rule[opt,_] | RuleDelayed[opt,_]]],
        {rule}
        ];

FSImportString[fileName_String] :=
    Module[{str = Import[fileName, "String"]},
           If[str =!= $Failed,
              str,
              "unknown"
             ]
          ];

End[];

EndPackage[];


BeginPackage["Utils`"];

Zip::usage = "Combines two lists to a list of touples.
Example:

In[]:= Zip[{a,b,c},{d,e,f}]
Out[]= {{a, d}, {b, e}, {c, f}}
";

StringZip::usage = "Combines two lists of strings to a list of concated strings.
Example:

In[]:= StringZip[{\"a\",\"b\"},{\"d\",\"e\"}]
Out[]= {\"ad\", \"be\"}}
";

Begin["`Private`"];

Zip[list1_List, list2_List] :=
    MapThread[List, {list1, list2}];

StringZip[list1_List, list2_List] :=
    MapThread[StringJoin, {list1, list2}];

End[];

EndPackage[];

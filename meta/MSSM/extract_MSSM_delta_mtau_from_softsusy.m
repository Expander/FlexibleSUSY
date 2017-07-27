(* Converts the 2-loop SUSY-QCD corrections O(alpha_s^2) to the DR-bar
   tau Yukawa coupling in the MSSM from GiNaC form to Mathematica
   form.

   The GiNaC expression in "dmtauas2.expr" has been extracted from
   SOFTSUSY 4.0.1 by adding the following C++ code snippet into the
   file
   src/two_loop_thresholds/two_loop_archives/tau_corrections.cpp
   right after the cache has been filled:

   if (cache.size() == 1) {
      ofstream out("dmtauas2.expr");
      out << cache[0] << endl;
   }

   Note: The expression does not include the 2-loop factor 1/(4 Pi)^4 .
 *)

str = Import["dmtauas2.expr", "String"];
ex  = ToExpression[StringReplace[str, "--" -> "+"], TraditionalForm];
ex >> "dmtauas2.m"

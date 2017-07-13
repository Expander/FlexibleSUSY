(* Converts the 2-loop SUSY corrections to the DR-bar alpha_s coupling
   in the MSSM from GiNaC form to Mathematica form.

   The GiNaC expression in "das2.expr" has been extracted from
   SOFTSUSY 4.0.2 by adding the following C++ code snippet into the
   file
   src/two_loop_thresholds/gs_corrections.cpp
   right after the cache has been filled:


   if (cache.size() == 1) {
      ofstream out("das2.expr");
      out << cache[0] << endl;
   }

   Note: The expression does not include the 2-loop factor 1/(4 Pi)^4 .
 *)

str = Import["das2.expr", "String"];
ex  = ToExpression[str, TraditionalForm];
ex >> "das2.m"

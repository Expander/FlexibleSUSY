(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

BeginPackage["LatticeUtils`", {
    "SARAH`",
    "TreeMasses`"}]

FixDiagonalization::usage;
SingleCase::usage;
Done::usage;
DoneLn::usage;
MajoranaQ::usage;
MajoranaMassMatrixQ::usage;

Begin["`Private`"]

FixDiagonalization[fsMassMatrices_List] := FixDiagonalization/@fsMassMatrices;

(*
   Diagonalization conventions of SARAH:

      SVD: m = u^T diag v,
	where u and v are the 1st and the 2nd mixing matrices from
	DEFINITION[_][MatterSector]

      hermitian: m = z^dagger diag z

   According to the SARAH documentation, the specification of the
   neutraliino mass matrix is indistinguishable from that of a
   hermitian matrix even though it must be diagonalized as

      symmetric: m = u^T diag u

   This leads to the following amendment:
 *)
FixDiagonalization[TreeMasses`FSMassMatrix[m_, f_, z_Symbol]?MajoranaMassMatrixQ] :=
    TreeMasses`FSMassMatrix[m, f, {z, z}];

FixDiagonalization[m_TreeMasses`FSMassMatrix] := m;

MajoranaMassMatrixQ[TreeMasses`FSMassMatrix[_?MatrixQ, _, _]?MajoranaMassQ] := True;

MajoranaMassMatrixQ[_TreeMasses`FSMassMatrix] := False;

MajoranaMassQ[TreeMasses`FSMassMatrix[_, _?MajoranaQ, _]] := True;

MajoranaMassQ[_TreeMasses`FSMassMatrix] := False;

MajoranaQ[field_] := MemberQ[SARAH`MajoranaPart, field];

SingleCase[args__] := Module[{
	cases = Cases[args]
    },
    Assert[Length[cases] === 1];
    First[cases]
];

SetAttributes[Done, HoldFirst];

Done[exp_, msg__] := Module[{
	result,
	time
    },
    WriteString["stdout", msg];
    result = Timing[exp];
    If[(time = Round[First[result] 1*^3]) === 0,
       time = ToString[Round[First[result] 1*^6]] <> " us",
       time = ToString[time] <> " ms"];
    WriteString["stdout", time];
    Last[result]
];

SetAttributes[DoneLn, HoldFirst];

DoneLn[exp_, msg__] := Module[{
	result = Done[exp, msg]
    },
    WriteString["stdout", "\n"];
    result
];

End[] (* `Private` *)

EndPackage[]

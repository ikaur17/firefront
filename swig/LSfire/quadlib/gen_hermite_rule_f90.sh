#!/bin/bash
#
#
#  http://people.sc.fsu.edu/~jburkardt/cpp_src/gen_hermite_rule/gen_hermite_rule.html
#
#
# GEN_HERMITE_RULE is a C++ program which generates a specific generalized Gauss-Hermite quadrature rule, based on user input.
#
# The rule is written to three files for easy use as input to other programs.
#
# The generalized Gauss Hermite quadrature rule is used as follows:
#
#         Integral ( -oo < x < +oo ) |x-a|^alpha * exp( - b * ( x - a)^2 ) f(x) dx
#      
# is to be approximated by
#         Sum ( 1 <= i <= order ) w(i) * f(x(i))
#      
# Usage:
#
# gen_hermite_rule order alpha a b filename
#
# where
#
# - order is the number of points in the quadrature rule.
# - alpha is the parameter for the generalized Gauss-Hermite quadrature rule. The value of alpha may be any real value greater than -1.0. Specifying alpha=0.0 results in the basic (non-generalized) rule.
# - a is the center point (default 0);
# - b is the scale factor (default 1);
# - filename specifies the names of the output files: filename_w.txt, filename_x.txt, and filename_r.txt, containing the weights, abscissas, and interval limits.

gfortran -c -g gen_hermite_rule.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gen_hermite_rule.f90"
  exit
fi
rm compiler.txt
#
gfortran gen_hermite_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gen_hermite_rule.o"
  exit
fi
rm gen_hermite_rule.o
#
chmod ugo+x a.out
mv a.out gen_hermite_rule_f90.exe 
echo "Executable created: gen_hermite_rule_f90.exe "
#mv a.out ~/bin/$ARCH/gen_hermite_rule
#
#echo "Executable installed as ~/bin/$ARCH/gen_hermite_rule"

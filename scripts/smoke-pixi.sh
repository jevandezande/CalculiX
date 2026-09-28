#!/bin/sh
set -eu

project_dir=$(pwd)
smoke_dir=$(mktemp -d)
trap 'rm -rf "$smoke_dir"' EXIT HUP INT TERM
cp test/spring1.inp "$smoke_dir/spring1.inp"
sed 's/^\*STATIC$/\*STATIC,SOLVER=SPOOLES/' test/oneel.inp > "$smoke_dir/oneel.inp"
cp test/spring4.inp "$smoke_dir/spring4.inp"
(
  cd "$smoke_dir"
  "$project_dir/build/pixi/CalculiX" -i spring1 > spring1.log 2>&1
  "$project_dir/build/pixi/CalculiX" -i oneel > oneel.log 2>&1
  "$project_dir/build/pixi/CalculiX" -i spring4 > spring4.log 2>&1
)
grep -q 'Job finished' "$smoke_dir/spring1.log"
grep -q '2  1.000000E-01' "$smoke_dir/spring1.dat"
grep -q 'symmetric spooles solver' "$smoke_dir/oneel.log"
grep -q 'Job finished' "$smoke_dir/oneel.log"
grep -q 'Calculating the eigenvalues and the eigenmodes' "$smoke_dir/spring4.log"
grep -q 'E I G E N V A L U E   O U T P U T' "$smoke_dir/spring4.dat"
printf 'Static, SPOOLES, and ARPACK smoke tests passed\n'

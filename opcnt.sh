#!/bin/bash

for op in vsubsd vaddsd vmulsd vfmadd[1-3]{3}sd vfmsub[1-3]{3}sd vfnmadd[1-3]{3}sd vfnmsub[1-3]{3}sd vdivsd vmovsd; do
  echo ${op}
  grep -E ${op} $1 | wc -l
done

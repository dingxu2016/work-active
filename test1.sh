#!/bin/bash

dt=0.03
v0=0.5
iseed=96
natom=1024
phi=0.6
step=200000
shear_ratio=0.05

run=/data1/dingxu/active_shear/one-fix_one-shear/build/test
cd $run

/data1/dingxu/active_shear/one-fix_one-shear/build/main<<EOF

$dt
$v0
$iseed
$natom
$phi
$step
$shear_ratio

EOF

exit 0;

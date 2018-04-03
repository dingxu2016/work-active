#!/bin/bash

dt=0.01
v0=0.1
iseed=81
natom=4096
phi=0.6
step=200000

run=/data1/dingxu/active_shear/one-fix_one-shear/build/test
cd $run

/data1/dingxu/active_shear/one-fix_one-shear/build/main<<EOF

$dt
$v0
$iseed
$natom
$phi
$step

EOF

exit 0;

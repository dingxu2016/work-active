#!/bin/bash

dt=0.01
v0=0.1
iseed=88
natom=4096
phi=0.75
step=200000

run=/data1/dingxu/active_shear/fix_all-boundary/build/test
cd $run

/data1/dingxu/active_shear/fix_all-boundary/build/main<<EOF

$dt
$v0
$iseed
$natom
$phi
$step

EOF

exit 0;

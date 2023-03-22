#!/bin/bash
echo "hello, Bshall!"
cd /albedo/work/user/chshao001/WRFPD
./compile em_real >& log.compilewrf
./compile em_tropical_cyclone >& log.compilecyc
cd /albedo/work/user/chshao001/WRFPD/test/em_tropical_cyclone
sbatch run_wrf.cmd

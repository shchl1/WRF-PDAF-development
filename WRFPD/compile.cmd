#!/bin/bash

cd /albedo/work/user/chshao001/WRFPD
./compile em_real >& log.compilewrf
./compile em_tropical_cyclone >& log.compilecyc


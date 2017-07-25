#!/usr/bin/python
import sys
import os
import string

for i in range(0,10):
         dt=0.001+i*0.001
         command='./md2d -nsteps 500000 -natoms 10 -dt %f -e 3 -box 20 20 -periodic 0 > n10_nstep500000_dt'%dt+str(dt)+'.energy'
         os.system(command)

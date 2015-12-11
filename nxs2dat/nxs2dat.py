#! /usr/bin/env python

import os
import sys

sys.path.append("/opt/Mantid/bin")
from mantid.simpleapi import *

nxsdir=sys.argv[1]
if len(sys.argv)>=3:
    asciidir=sys.argv[2]
else:
    asciidir=os.getcwd()
config.setDataSearchDirs(nxsdir)
config['defaultsave.directory']=asciidir

if not os.path.exists(asciidir):
    os.makedirs(asciidir)
    print "Info: "+asciidir+" does not exist and will be created."


for file in os.listdir(nxsdir):
    if file.endswith(".nxs"):
        LoadNexusProcessed(Filename=file, OutputWorkspace='nxs2dat', LoadHistory=False)
        #Rebin(InputWorkspace='nxs2dat', OutputWorkspace='nxs2dat', Params='-2,0.015,5,-0.005,1000', PreserveEvents=False)
        asciifile=file.strip(".nxs")+".dat"
        SaveAscii(InputWorkspace='nxs2dat',Filename=asciifile,Separator='Space')


from os.path import exists
from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('shaper') #

import sys

inputPassName="hits"
nEv=400000

if len(sys.argv) > 2 :
    timeSample=int(sys.argv[2])
else :
    timeSample=21
    
from LDMX.TrigScint.trigScint import *


shaper =dataShaper("tShaped")
shaper.input_pass_name=inputPassName
shaper.input_pass_name2="unpack"
#tbClustersUp.input_collection="TestBeamHitsUp"
shaper.pad_time=0.
shaper.time_tolerance=50.
shaper.verbosity=0
p.sequence = [
    shaper
    ]
p.histogramFile = sys.argv[1].replace(".root","_dataShaped.root")

#generate on the fly
p.inputFiles = [sys.argv[1]]
#p.outputFiles = [ sys.argv[1].replace(".root", "_shaped.root") ]
p.maxEvents = nEv

p.termLogLevel = 2
p.logFileName = sys.argv[1].replace(".root", "_dataShaped.log") 
p.logFileLevel=0


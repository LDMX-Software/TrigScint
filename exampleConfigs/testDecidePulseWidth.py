import json 
from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('plot') #

import sys

#inputPassName="hits"
nEv=3000

if len(sys.argv) > 2 :
    startSample=int(sys.argv[2])
#else :
    #startSample=12


#from LDMX.TrigScint.trigScint import TestBeamHitAnalyzer
from LDMX.TrigScint.trigScint import TestBeamDecideWidth


# ------------------- all set; setup in detail, and run with these settings ---------------

#tsEv=TestBeamHitAnalyzer("plotMaker")
tsEv=TestBeamDecideWidth("Charge")
tsEv.inputPassName=''
tsEv.nInstrumentedChannels=12
# now in default config, too, but with test beam values :
#these are derived as the mean of gaussian fits to the "event pedestal" (average over middle two quartiles) for each channel
#tsEv.startSample=startSample
tsEv.inputCollection = "QIEsamplesUp"
#tsEv.gain=2e6
tsEv.gain=[1.89894e+06,1.97174e+06,1.90704e+06,1.9782e+06,1.95579e+06,1.95556e+06,1.88863e+06,1.84093e+06,1.84093e+06,1.88551e+06,1.90385e+06,1.90495e+06]                 

p.sequence = [
    tsEv
    ]


#generate on the fly
p.inputFiles = [sys.argv[1]]
outname=sys.argv[1]
outname=outname.replace(".root", "_width.root")
#p.outputFiles = [ outname ]

p.histogramFile = outname #.replace(".root"

p.maxEvents = nEv

p.logFileName=outname.replace(".root",".log")
p.termLogLevel = 2
p.logFileLevel=1#0

json.dumps(p.parameterDump(), indent=2)

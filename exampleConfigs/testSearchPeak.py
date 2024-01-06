import json 
from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('plot') #

import sys

#inputPassName="hits"
nEv=80000

if len(sys.argv) > 3 :
    startSample=int(sys.argv[3])
#else :
    #startSample=12


#from LDMX.TrigScint.trigScint import TestBeamHitAnalyzer
from LDMX.TrigScint.trigScint import TestBeamSearchPeak


# ------------------- all set; setup in detail, and run with these settings ---------------

#tsEv=TestBeamHitAnalyzer("plotMaker")
tsEv=TestBeamSearchPeak("Charge")
tsEv.inputPassName=''
tsEv.nInstrumentedChannels=12
# now in default config, too, but with test beam values :
#these are derived as the mean of gaussian fits to the "event pedestal" (average over middle two quartiles) for each channel
#tsEv.startSample=startSample
tsEv.inputCollection = "QIEsamplesUp"
#tsEv.gain=2e6
#tsEv.gain=[1.89894e+06,1.97174e+06,1.90704e+06,1.9782e+06,1.95579e+06,1.95556e+06,1.88863e+06,1.84093e+06,1.84093e+06,1.88551e+06,1.90385e+06,1.90495e+06]                 
gainList=[1871449.488636,1938213.246753,1889102.768362,1947630.000000,1923191.472393,1919736.424242,1862303.588235,1830186.931818,1844377.142857,1875542.807018,1887319.775281,1889952.458101]
pedList = [2.425455,-0.654356,-1.064245,3.278855,2.382266,-1.818987,2.441944,-1.604241,0.633108,-0.812120,1.946567,1.634450]
tsEv.gain = gainList
tsEv.pedestals = pedList
p.sequence = [
    tsEv
    ]


#generate on the fly
p.inputFiles = [sys.argv[1]]
outname=sys.argv[2]
#outname=outname.replace(".root", "_width.root")
#p.outputFiles = [ outname ]
#p.outputFiles = [sys.argv[2]]
p.histogramFile = outname #.replace(".root"

p.maxEvents = nEv

p.logFileName=outname.replace(".root",".log")
p.termLogLevel = 2
p.logFileLevel=1#0

json.dumps(p.parameterDump(), indent=2)

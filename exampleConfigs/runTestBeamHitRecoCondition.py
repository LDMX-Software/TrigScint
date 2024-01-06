from os.path import exists
from os import path
from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('hits') #

import sys

inputPassName="conv"
nEv=40000

#if len(sys.argv) > 2 :
    #timeSample=int(sys.argv[2])
#else :
    #timeSample=15
    
from LDMX.TrigScint.trigScint import TestBeamHitProducerCondition
#import LDMX.TrigScint.ts_hardcoded_conditions
import LDMX.TrigScint.ts_testbeam_conditions

nChannels=12
#gainList=[2e6]*nChannels
#3gainList=[1.89894e+06,1.97174e+06,1.90704e+06,1.9782e+06,1.95579e+06,1.95556e+06,1.88863e+06,1.84093e+06,1.84093e+06,1.88551e+06,1.90385e+06,1.90495e+06]
#now if there is a gain file, use that instead to read in the gain for each channel
gainFileName=sys.argv[1].replace("_linearize.root", "_gains.txt")
gainFileName=gainFileName.replace("_adcTrig", "")  #not derived for adcTrig events 
startSampleList=[12]*nChannels
startSampleFileName=gainFileName.replace("gains", "width")
#defaultPedFileName=dataPath+"/"+defaultRun+"_peds.txt"

if not exists(startSampleFileName) :  
    startSampleList=[12]*nChannels
    #pedFileName=defaultPedFileName

if exists(startSampleFileName) :
    with open(startSampleFileName) as f:
        for line in f.readlines() :
            line=line.split(',')  #values are comma separated, one channel per line: channelNB, ped
            startSampleList[ int(line[0].strip()) ] = int(line[1].strip())

print("Using this list of start samples:")
print(startSampleList)

tbHitsUp  =TestBeamHitProducerCondition("tbHits")
tbHitsUp.input_pass_name=inputPassName
tbHitsUp.input_collection="QIEsamplesUp"
#tbHitsUp.pedestals=pedList
#tbHitsUp.gain=gainList 
tbHitsUp.startSample=startSampleList
tbHitsUp.pulseWidth=5 #5 
tbHitsUp.pulseWidthLYSO = 5#7 
tbHitsUp.doCleanHits=True
tbHitsUp.nInstrumentedChannels=12
p.sequence = [
    tbHitsUp
    ]


#generate on the fly
p.inputFiles = [sys.argv[1]]
p.outputFiles = [ sys.argv[1].replace(".root", "_18_5_5_hits_full_general.root") ]
p.maxEvents = nEv

p.termLogLevel = 2

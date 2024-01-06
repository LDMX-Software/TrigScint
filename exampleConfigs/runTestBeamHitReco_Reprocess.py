from os.path import exists
from os import path
import numpy as np
from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('hits') #
#p = ldmxcfg.Process('run')
import sys

inputPassName="conv"
nEv=700000
#nEv=4
#if len(sys.argv) > 2 :
    #timeSample=int(sys.argv[2])
#else :
    #timeSample=15
    
from LDMX.TrigScint.trigScint import TestBeamHitProducer


nChannels=12
_,RunID_Gain,GEG0,GEG1,GEG2,GEG3,GEG4,GEG5,GEG6,GEG7,GEG8,GEG9,GEG10,GEG11 = np.loadtxt("/home/dhruvanshu/LDMX_Analysis_Files/ParametersExtracted_ReprocessStartSamples/GainData_Poly.txt",unpack=True,delimiter=',',dtype=str)
RunIDint=RunID_Gain.astype(np.int32)
ParGain = [GEG0.astype(np.float64),GEG1.astype(np.float64),GEG2.astype(np.float64),GEG3.astype(np.float64),GEG4.astype(np.float64),GEG5.astype(np.float64),GEG6.astype(np.float64),GEG7.astype(np.float64),GEG8.astype(np.float64),GEG9.astype(np.float64),GEG10.astype(np.float64),GEG11.astype(np.float64)]
gainList=[1871449.488636,1938213.246753,1889102.768362,1947630.000000,1923191.472393,1919736.424242,1862303.588235,1830186.931818,1844377.142857,1875542.807018,1887319.775281,1889952.458101]
print(type(sys.argv[2]))
for ir in range(len(RunIDint)):
	if RunIDint[ir] == int(sys.argv[2]):
		gainList=[ParGain[0][ir],ParGain[1][ir],ParGain[2][ir],ParGain[3][ir],ParGain[4][ir],ParGain[5][ir],ParGain[6][ir],ParGain[7][ir],ParGain[8][ir],ParGain[9][ir],ParGain[10][ir],ParGain[11][ir]]
		
#gainList=[1902940.000000,1975320.000000,1912200.000000,1989050.000000,1961480.000000,1954340.000000,1888030.000000,1850080.000000,113605.000000,1889230.000000,1911360.000000,1910540.000000]
print("Using this list of gains:")
print(gainList)

pedList = [2.425455,-0.654356,-1.064245,3.278855,2.382266,-1.818987,2.441944,-1.604241,0.633108,-0.812120,1.946567,1.634450]
print("Using this list of pedestals:")
print(pedList)

startSampleList=[14]*nChannels
print("Here:")
#startSampleFileName=gainFileName.replace("gains", "width")
#defaultPedFileName=dataPath+"/"+defaultRun+"_peds.txt"
startSampleFileName=sys.argv[1].replace(".root", "_linearize_width.txt")
#startSampleFileName=sys.argv[3]
print(startSampleFileName)

if not exists(startSampleFileName) :  
    startSampleList=[14]*nChannels
    #pedFileName=defaultPedFileName

if exists(startSampleFileName) :
    with open(startSampleFileName) as f:
        for line in f.readlines() :
            line=line.split(',')  #values are comma separated, one channel per line: channelNB, ped
            startSampleList[ int(line[0].strip()) ] = int(line[1].strip())

#startSampleList=[17,17,17,17,17,17,15,15,16,15,17,17]
#startSampleList=[6,6,6,6,6,6,6,6,6,6,6,6]
startSampleList=[19,19,19,19,19,19,16,16,16,16,19,19]

print("Using this list of start samples:")
print(startSampleList)

#MIPresponseList = [0.999006,1.0028,0.987226,0.934034,1.05343,0.941591,1.03258,1.02464,0,1.0454,1.05341,1.04235]
MIPresponseList = [1.0]*nChannels
print("Using this list of MIP responses:")
print(MIPresponseList)

tbHitsUp=TestBeamHitProducer("tbHits")
tbHitsUp.input_pass_name=inputPassName
tbHitsUp.input_collection="QIEsamplesUp"
tbHitsUp.pedestals=pedList
tbHitsUp.gain=gainList 
tbHitsUp.startSample=startSampleList
tbHitsUp.MIPresponse=MIPresponseList
tbHitsUp.pulseWidth=1 #5 
tbHitsUp.pulseWidthLYSO = 1#7 
tbHitsUp.doCleanHits=True
tbHitsUp.nInstrumentedChannels=12
p.sequence = [
    tbHitsUp
    ]


#generate on the fly
p.inputFiles = [sys.argv[1]]
#p.outputFiles = [sys.argv[2]]
p.outputFiles = [ sys.argv[1].replace(".root", "_only_peak_timesample.root") ]
#p.outputFiles = [ sys.argv[1].replace(".root", "_v2.root") ]
p.maxEvents = nEv
#p.run = 179

p.termLogLevel = 2

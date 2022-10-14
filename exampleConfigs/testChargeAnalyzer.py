from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('plot')
import sys
from LDMX.TrigScint.trigScint import ChargeAnalyzer
from LDMX.TrigScint.trigScint import TestBeamHitAnalyzer
from LDMX.TrigScint.trigScint import ADCAnalyzer
#from ROOT import TTree, TFile
#p.maxEvents = 1
#p.maxEvents = 15071
#p = ldmxcfg.Process('plot')
tbh_analyzer = TestBeamHitAnalyzer('Hits')
tbh_analyzer.inputCollection = "testBeamHitsUp"
tbh_analyzer.inputPassName = ""

Q_analyzer = ChargeAnalyzer('Charge')
Q_analyzer.inputCollection = "QIEsamplesUp"
Q_analyzer.inputPassName = ""

ADC_analyzer = ADCAnalyzer('ADC')
ADC_analyzer.inputCollection = "decodedQIEUp"
ADC_analyzer.inputPassName = ""

p.sequence = [ADC_analyzer,Q_analyzer,tbh_analyzer]
p.inputFiles = [sys.argv[1]]
p.histogramFile = sys.argv[1].replace(".root","_multiple_analysed.root")
print("Analyzed")

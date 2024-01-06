#Python config file used to run ChargeAnalyzerCrazy.cxx which plots the charge vs timesamples based on hit flags for a particular event of a run. Event number can be passed as an input argument to
#this script
import json
import sys
from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('plot') 

from LDMX.TrigScint.trigScint import ChargeAnalyzerCrazy
from LDMX.TrigScint.trigScint import QualityFlagAnalyzer
#inputPassName="conv"
#p.maxEvents = 1
Q_analyzer = ChargeAnalyzerCrazy('Charge')
Q_analyzer.inputCollection = "QIEsamplesUp"
Q_analyzer.inputPassName = ""
Q_analyzer.EventDisplayNumber = int(sys.argv[2])

F_analyzer = QualityFlagAnalyzer('Flag')
F_analyzer.inputCollection = "QIEsamplesUp"
F_analyzer.inputPassName = ""
#F_analyzer.EventDisplayNumber = int(sys.argv[2])

p.sequence = [Q_analyzer, F_analyzer]
p.inputFiles = [sys.argv[1]]
p.histogramFile = sys.argv[1].replace(".root","_multiple_analysed_crazy.root")
json.dumps(p.parameterDump(), indent=2)



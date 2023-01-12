import json
import sys
from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('plot') 

from LDMX.TrigScint.trigScint import ChargeAnalyzerCrazy

#inputPassName="conv"
#p.maxEvents = 1

Q_analyzer = ChargeAnalyzerCrazy('Charge')
Q_analyzer.inputCollection = "QIEsamplesUp"
Q_analyzer.inputPassName = ""

p.sequence = [Q_analyzer]
p.inputFiles = [sys.argv[1]]
p.histogramFile = sys.argv[1].replace(".root","_multiple_analysed_crazy.root")
json.dumps(p.parameterDump(), indent=2)



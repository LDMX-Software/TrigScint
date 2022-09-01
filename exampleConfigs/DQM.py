from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('DQM')
import sys
from LDMX.TrigScint.trigScint import *

p.maxEvents = 1000
tbh_analyzer = TestBeamHitAnalyzer('test')
tbh_analyzer.inputCollection = "testBeamHitsUp"
tbh_analyzer.inputPassName = ""
p.sequence = [tbh_analyzer]
p.inputFiles = [sys.argv[1]]
p.histogramFile = sys.argv[1].replace(".root","_DQM.root")

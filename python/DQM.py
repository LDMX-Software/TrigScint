from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('bbb')
import sys
from LDMX.TrigScint.trigScint import *
dqm = TestBeamDQM('aaa')
dqm.inputCollection='decodedQIEUp'
#dqm.inputPassName='unpack'
p.maxEvents = 10
p.sequence = [dqm]

p.inputFile = [sys.argv[1]]
p.outputFiles = [sys.argv[1].replace(".root","_DQM.root")]
p.histogramFile = "myHistFile.root"

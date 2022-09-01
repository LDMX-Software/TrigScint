from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('DQM')
import sys
from LDMX.TrigScint.trigScint import *

p.maxEvents = 1000
p.sequence = [TestBeamProducer()]
p.inputFile = [sys.argv[1]]
p.outputFiles = [sys.argv[1].replace(".root","_DQM.root")]

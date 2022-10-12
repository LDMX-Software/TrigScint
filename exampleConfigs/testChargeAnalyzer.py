from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('plot')
import sys
from LDMX.TrigScint.trigScint import ChargeAnalyzer
#from ROOT import TTree, TFile
#p.maxEvents = 1
#p.maxEvents = 15071
tbh_analyzer = ChargeAnalyzer('test')
tbh_analyzer.inputCollection = "QIEsamplesUp"
tbh_analyzer.inputPassName = ""
p.sequence = [tbh_analyzer]
p.inputFiles = [sys.argv[1]]
#p.inputFiles = ['/home/dhruvanshu/ldmx-sw/unpacked_ldmx_captan_out_13-04-2022_00-00-38__187_reformat_30timeSamplesFrom0_linearize_hits_clusters.root']
p.histogramFile = sys.argv[1].replace(".root","_DQM_Qanalyzed.root")

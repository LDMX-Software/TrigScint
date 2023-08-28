import sys
from LDMX.Framework import ldmxcfg
from LDMX.SimCore import generators
from LDMX.SimCore import simulator


thisPassName="sim"
p = ldmxcfg.Process(thisPassName)

gunZpos=1100 #3000  #mm -- define as positive here, for file naming; set sign below
detV=2        #detector geometry version number 
beamXsmear=7.5 #mm
beamYsmear=20 #mm
noisePerEvent=1.  #average number of PEs from SiPM noise per event (gets scaled by nTimeSamples to be constant)
startSample=17.
INDEX=0
if len(sys.argv) > 1 :
    INDEX=int(sys.argv[1])
NUM=0
if len(sys.argv) > 2 :
    NUM=int(sys.argv[2]) 
#DETECTORS=['ldmx-hcal-prototype-v'+str(detV)+'.0Gap0','ldmx-hcal-prototype-v'+str(detV)+'.0Gap1','ldmx-hcal-prototype-v'+str(detV)+'.0Gap2','ldmx-hcal-prototype-v'+str(detV)+'.0Gap3']
DETECTORS=["ldmx-det-v14"]

#['ldmx-hcal-prototype-v'+str(detV)+'.0Shift0','ldmx-hcal-prototype-v'+str(detV)+'.0Shift1','ldmx-hcal-prototype-v'+str(detV)+'.0Shift2','ldmx-hcal-prototype-v'+str(detV)+'.0Shift3']

#if len(sys.argv) > 1 :
#    nTimeSamples=int(sys.argv[1])
#else :
nTimeSamples=30 #config default is 5 
#if len(sys.argv) > 2 :
#    elecNoise=float(sys.argv[2])
#else :
elecNoise=1.5 #config default is 1.5
#if len(sys.argv) > 3 :
#    kExpo=float(sys.argv[3])
#else :
kExpo=0.1   #config default is 0.1
    
p.run = 10
p.maxEvents = 20000
#'testbeamSim_zNeg'+str(gunZpos)+'mm_beamSpot'+str(beamXsmear)+'x'+str(beamYsmear)+'mm_'+str(nTimeSamples)+'tSamp_eNoise'+str(elecNoise)+'_tauInv'+str(kExpo)+'_detV'+str(detV)+'_'+str(p.maxEvents)+'ev.root']
p.outputFiles = ["Shnap"+str(INDEX)+str(NUM)+".root"]
print("Producing output file: "+p.outputFiles[0])

gunZpos=-float(gunZpos)  #get sign right, and make floats, to use as parameters
beamXsmear=float(beamXsmear)
beamYsmear=float(beamYsmear)


#SINGLE PARTICLE GUN

#Theta=[0.0,0.5235,1.047197,1.570796]
gun = generators.gun('particle_gun')
gun.particle = 'e-'
import math
gun.direction=[0.,0.,1.0]#[-1.0*math.sin(Theta[INDEX]),0.,1.0*math.cos(Theta[INDEX])]
#0., 0., 1.]
#gunZpos*math.sin(Theta[INDEX])
gun.position=[0.,0.,gunZpos]#[10.0*math.sin(Theta[INDEX]),0.,-500+10.0*math.cos(Theta[INDEX])]#[ -597.5*math.sin(Theta[INDEX]) , 0.,gunZpos+597.5-597.5*math.cos(Theta[INDEX]) ]#*math.cos(Theta[INDEX])]
gun.energy = 4.   #gev

#print(gun.direction[0])
#print(gun.position[0])

#MULTIPLE PARTICLE GUN

#INDEX=0
mpgGen = generators.multi( "mgpGen" ) # this is the line that actually creates the generator
mpgGen.vertex = [0.,0.,gunZpos]#[ -27.926, 0., -700 ] # mm
mpgGen.nParticles = 2
mpgGen.pdgID = 11
#doPoisson=False
#mpgGen.enablePoisson = doPoisson
import math
theta = 0.0 #math.radians(Theta[INDEX])
mpgGen.momentum = [ 4000.*math.sin(theta) , 0, 4000.*math.cos(theta) ]


simulation = simulator.simulator('test_TS')
simulation.rootPrimaryGenUseSeed = True
simulation.generators=[gun]
simulation.setDetector("ldmx-hcal-measured2-v2.0")#"ldmx-hcal-prototype-v2.0")#"ldmx-hcal-measured-v2.0")
simulation.beamSpotSmear = [beamXsmear, beamYsmear, 0] #mm, at start position


from LDMX.Hcal import HcalGeometry
import LDMX.Ecal.EcalGeometry
from LDMX.TrigScint.trigScint import TrigScintQIEDigiProducer
from LDMX.TrigScint.trigScint import TrigScintRecHitProducer
from LDMX.TrigScint.trigScint import TrigScintClusterProducer

tsDigis = TrigScintQIEDigiProducer.pad1()
#tsDigis.randomSeed(NUM)
tsDigis.number_of_strips = 12
tsDigis.input_collection = "TriggerPadUpSimHits"
tsDigis.mean_noise=float(noisePerEvent/nTimeSamples)
tsDigis.maxts = nTimeSamples
tsDigis.elec_noise = elecNoise
tsDigis.expo_k = kExpo
#tsDigis.sipm_gain = 2.e6
tsDigis.zeroSupp_in_pe=0.5
tsDigis.toff_overall = 25.*startSample
tsDigis.pe_per_mip = 110. #from fit to 2-hit clusters in data run 183, April4 2310, all plastic

tsRecHits = TrigScintRecHitProducer.pad1()
tsRecHits.pe_per_mip = tsDigis.pe_per_mip #pick it up here too
#tsRecHits.gain = tsDigis.gain # same


tsCl = TrigScintClusterProducer.pad1()
tsCl.input_collection = tsRecHits.output_collection
tsCl.pad_time = 100.
tsCl.time_tolerance = 999.
#tsCl.verbosity = 3
tsCl.clustering_threshold = 30.  #to add in neighboring
tsCl.seed_threshold = 40.



from LDMX.TrigScint.trigScint import EventReadoutProducer

tsEv=EventReadoutProducer("eventLinearizer")
tsEv.input_pass_name=thisPassName
tsEv.input_collection=tsDigis.output_collection
tsEv.time_shift=0 #timeOffset

nChannels=12
gainList=[tsDigis.sipm_gain]*nChannels
pedList=[tsDigis.pedestal]*nChannels


from LDMX.TrigScint.trigScint import QIEAnalyzer

tsAna=QIEAnalyzer("plotMaker")
tsAna.inputPassName=thisPassName
tsAna.startSample=0
tsAna.pedestals=pedList
tsAna.gain=gainList

outname=p.outputFiles[0].replace(".root", "_plots.root")

from LDMX.TrigScint.trigScint import dataShaper

shaper=dataShaper("tShaped")
shaper.input_pass_name="sim"#inputPassName
shaper.input_pass_name2="sim"#"unpack"
shaper.input_pass_name3="sim"#"hits"
#tbClustersUp.input_collection="TestBeamHitsUp"
shaper.pad_time=0.
shaper.time_tolerance=50.
shaper.verbosity=0

p.histogramFile = outname 

p.sequence=[simulation,
                  tsDigis,
                  tsRecHits,
                  tsCl,
                  tsEv,
                  tsAna,
                  shaper
                  ]


p.termLogLevel = 3 #0

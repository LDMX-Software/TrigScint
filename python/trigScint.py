"""Configuration for Trigger Scintillator digitization, cluster, and track producers

Sets all parameters to reasonable defaults.

Examples
--------
    from LDMX.TrigScint.trigScint import TrigScintDigiProducer
    p.sequence.extend([ TrigScintDigiProducer.up() , TrigScintDigiProducer.down() , TrigScintDigiProducer.tagger() ])
    from LDMX.TrigScint.trigScint import TrigScintClusterProducer
    p.sequence.extend([ TrigScintClusterProducer.up() , TrigScintClusterProducer.down() , TrigScintClusterProducer.tagger() ])

"""

from LDMX.Framework import ldmxcfg

class TrigScintDigiProducer(ldmxcfg.Producer) :
    """Configuration for digitizer for Trigger Scintillators"""

    def __init__(self,name) :
        super().__init__(name,'trigscint::TrigScintDigiProducer','TrigScint')

        self.mean_noise = 0.02
        self.number_of_strips = 50
        self.number_of_arrays = 1
        self.mev_per_mip = 0.4
        self.pe_per_mip = 100.
        self.input_collection="TriggerPadUpSimHits"
        self.input_pass_name="" #take any pass
        self.output_collection="trigScintDigisUp"
        import time
        self.randomSeed = int(time.time())
        self.verbose = False

    def up() :
        """Get the digitizer for the trigger pad upstream of target"""
        digi = TrigScintDigiProducer( 'trigScintDigisUp' )
        digi.input_collection = 'TriggerPadUpSimHits'
        digi.output_collection= 'trigScintDigisUp'
        return digi

    def down() :
        """Get the digitizer for the trigger pad downstream of target"""
        digi = TrigScintDigiProducer( 'trigScintDigisDn' )
        digi.input_collection = 'TriggerPadDownSimHits'
        digi.output_collection= 'trigScintDigisDn'
        return digi

    def tagger() :
        """Get the digitizer for the trigger pad upstream of tagger"""
        digi = TrigScintDigiProducer( 'trigScintDigisTag' )
        digi.input_collection = 'TriggerPadTaggerSimHits'
        digi.output_collection= 'trigScintDigisTag'
        return digi


class TrigScintQIEDigiProducer(ldmxcfg.Producer) :
    """Configuration for digitizer for Trigger Scintillators's QIE chip"""

    def __init__(self,name) :
        super().__init__(name,'trigscint::TrigScintQIEDigiProducer','TrigScint')

        self.mean_noise = 0.02
        self.number_of_strips = 50
        self.number_of_arrays = 1
        self.mev_per_mip = 0.4
        self.pe_per_mip = 100.
        self.input_collection="TriggerPadUpSimHits"
        self.input_pass_name="" #take any pass
        self.output_collection="trigScintQIEDigisUp"
        self.input_pulse_shape="Expo" # Name of the input pulse class
        self.expo_k=0.1          # Inverse of decay time of piece-wise exponential 
        self.expo_tmax=5.0       # Time at which piece-wise exponential peaks
        self.maxts=5             # No. of time samples to analyze
        self.toff_overall = 55.0 # Global time offset
        self.tdc_thr = 3.4       # Threshold current in uA for TDC latch
        self.pedestal= 6.0       # QIE pedestal value (in fC)
        self.elec_noise = 1.5    # Electronic noise (in fC)
        self.sipm_gain = 1.e6    # SiPM Gain
        self.qie_sf = 40.        # QIE sampling frequency in MHz

        import time
        self.verbose = False

    def up() :
        """Get the digitizer for the trigger pad upstream of target"""
        digi = TrigScintQIEDigiProducer( 'trigScintQIEDigisUp' )
        digi.input_collection = 'TriggerPadUpSimHits'
        digi.output_collection= 'trigScintQIEDigisUp'
        return digi

    def down() :
        """Get the digitizer for the trigger pad downstream of target"""
        digi = TrigScintQIEDigiProducer( 'trigScintQIEDigisDn' )
        digi.input_collection = 'TriggerPadDownSimHits'
        digi.output_collection= 'trigScintQIEDigisDn'
        return digi

    def tagger() :
        """Get the digitizer for the trigger pad upstream of tagger"""
        digi = TrigScintQIEDigiProducer( 'trigScintQIEDigisTag' )
        digi.input_collection = 'TriggerPadTaggerSimHits'
        digi.output_collection= 'trigScintQIEDigisTag'
        return digi



class TrigScintRecHitProducer(ldmxcfg.Producer) :
    """Configuration for rechit producer for Trigger Scintillators"""

    def __init__(self,name) :
        super().__init__(name,'trigscint::TrigScintRecHitProducer','TrigScint')

        self .mev_per_mip = 0.4   #\
                                  # >>>both are for converting edep to PEs 
        self.pe_per_mip = 100.    #/
        self.pedestal= 6.0        # QIE pedestal value (in fC)
        self.gain = 1.e6      # SiPM Gain
        self.input_collection="trigScintQIEDigisUp"
        self.input_pass_name=""   #take any pass
        self.output_collection="trigScintRecHitsUp"
        self.verbose = False
        self.sample_of_interest=2 # Sample of interest. Range 0 to 3
        self.En_Reco_Option = 0   # Toggle Energy reconstruction algorithm

        self.input_pulse_shape="Expo" # Name of the input pulse class
        self.expo_k=0.1          # Inverse of decay time of piece-wise exponential
        self.expo_tmax=5.0       # Time at which piece-wise exponential peaks

    def up() : 
        """Get the rechit producer for upstream pad"""
        rechit = TrigScintRecHitProducer( 'trigScintRecHitsUp' )
        rechit.input_collection  = 'trigScintQIEDigisUp'
        rechit.output_collection = 'trigScintRecHitsUp'
        return rechit

    def down() : 
        """Get the rechit producer for downstream pad"""
        rechit = TrigScintRecHitProducer( 'trigScintRecHitsDown' )
        rechit.input_collection  = 'trigScintQIEDigisDn'
        rechit.output_collection = 'trigScintRecHitsDn'
        return rechit

    def tagger() : 
        """Get the rechit producer for tagger pad"""
        rechit = TrigScintRecHitProducer( 'trigScintRecHitsTag' )
        rechit.input_collection  = 'trigScintQIEDigisTag'
        rechit.output_collection = 'trigScintRecHitsTag'
        return rechit

class TrigScintClusterProducer(ldmxcfg.Producer) :
    """Configuration for cluster producer for Trigger Scintillators"""

    def __init__(self,name) :
        super().__init__(name,'trigscint::TrigScintClusterProducer','TrigScint')

        self.max_cluster_width = 2
        self.clustering_threshold = 0.  #to add in neighboring channels
        self.seed_threshold = 30.
        self.pad_time = 0.
        self.time_tolerance = 0.5
        self.input_collection="trigScintDigisTag"
        self.input_pass_name="" #take any pass
        self.output_collection="TriggerPadTaggerClusters"
        self.verbosity = 0

    def up() :
        """Get the cluster producer for the trigger pad upstream of target"""
        cluster = TrigScintClusterProducer( 'trigScintClustersUp' )
        cluster.input_collection = 'trigScintDigisUp'
        cluster.output_collection= 'TriggerPadUpClusters'
        cluster.pad_time= 0.
        return cluster

    def down() :
        """Get the cluster producer for the trigger pad downstream of target"""
        cluster = TrigScintClusterProducer( 'trigScintClustersDown' )
        cluster.input_collection = 'trigScintDigisDn'
        cluster.output_collection= 'TriggerPadDownClusters'
        cluster.pad_time= 0.
        return cluster

    def tagger() :
        """Get the cluster producer for the trigger pad upstream of tagger"""
        cluster = TrigScintClusterProducer( 'trigScintClustersTag' )
        cluster.input_collection = 'trigScintDigisTag'
        cluster.output_collection= 'TriggerPadTaggerClusters'
        cluster.pad_time= -2.
        return cluster


from LDMX.Framework import ldmxcfg

class TrigScintTrackProducer(ldmxcfg.Producer) :
    """Configuration for track producer for Trigger Scintillators"""

    def __init__(self,name) :
        super().__init__(name,'trigscint::TrigScintTrackProducer','TrigScint')

        self.delta_max = 0.75
        self.tracking_threshold = 0.  #to add in neighboring channels
        self.seeding_collection = "TriggerPadTaggerClusters"
        self.further_input_collections = ["TriggerPadUpClusters","TriggerPadDownClusters"]
        self.input_pass_name="" #take any pass
        self.output_collection="TriggerPadTracks"
        self.verbosity = 0

trigScintTrack = TrigScintTrackProducer( "trigScintTrack" )



#include "TrigScint/dataShaper.h"

#include <iterator>  // std::next
#include <map>

namespace trigscint {

void dataShaper::configure(framework::config::Parameters &ps) {
  seed_ = ps.getParameter<double>("seed_threshold");
  minThr_ = ps.getParameter<double>("clustering_threshold");
  maxWidth_ = ps.getParameter<int>("max_cluster_width");
  input_collection_ = ps.getParameter<std::string>("input_collection");
  input_collection2_ = ps.getParameter<std::string>("input_collection2");
  input_collection3_ = ps.getParameter<std::string>("input_collection3");
  passName_ = ps.getParameter<std::string>("input_pass_name");
  passName2_ = ps.getParameter<std::string>("input_pass_name2");
  passName3_ = ps.getParameter<std::string>("input_pass_name3");
  output_collection_ = ps.getParameter<std::string>("output_collection");
  verbose_ = ps.getParameter<int>("verbosity");

  timeTolerance_ = ps.getParameter<double>("time_tolerance");
  padTime_ = ps.getParameter<double>("pad_time");
  if (verbose_) {
    ldmx_log(info) << "In TestBeamClusterProducer: configure done!";
    ldmx_log(info) << "Got parameters: \nSeed threshold:   " << seed_
                   << "\nClustering threshold: " << minThr_
                   << "\nMax cluster width: " << maxWidth_
                   << "\nExpected pad hit time: " << padTime_
                   << "\nMax hit time delay: " << timeTolerance_
				   << "\n\t doCleanHits = " << doCleanHits_
                   << "\nInput collection:     " << input_collection_
                   << "\nInput pass name:     " << passName_
                   << "\nOutput collection:    " << output_collection_
                   << "\nVerbosity: " << verbose_;
  }

  return;
}

void dataShaper::produce(framework::Event &event) {
 if (verbose_) {
    ldmx_log(debug)
        << "dataShaper: produce() starts! Event number: "
        << event.getEventHeader().getEventNumber();
  }

  // looper over digi hits and aggregate energy depositions for each detID
  SimQIE qie;
  //float maxEnergy=0;
  for(int I=0; I<17; I++){
    integratedCharge_[I]=0.0;
    PENumber_[I]=0.0;
    TimeSampleWMaxCharge_[I]=-400.0;
    if(I<5){
        populated_[I]=0;
        chargeInCluster_[I]=0;
        hitNumInCluster_[I]=0;
        centroidOfCluster_[I]=0;
    }
  }
  overflow_=0;
  const auto digisO{event.getCollection<trigscint::TrigScintQIEDigis>(input_collection2_, passName2_)};
  for (const auto &digiO : digisO) {
    //digiO.getChanID();
    auto adc{digiO.getADC()};
    auto tdc{digiO.getADC()};
    //std::cout<<adc.size()<<std::endl;
    float charge=0;
    float maxEnergy=-1000;
    for(int I=0;I<30;I++){
      adcsList_[I]=adc[I];
      tdcsList_[I]=tdc[I];
      QList_[I]=qie.ADC2Q(adc[I]);
      if(QList_[I]>maxEnergy){
        TimeSampleWMaxCharge_[(int)digiO.getChanID()]=I;
        maxEnergy=QList_[I];
      }
    }
  }
  const auto digis{
      event.getCollection<trigscint::TestBeamHit>(input_collection3_, passName3_)};
  //TrigScintCluster
  const auto clusters{
      event.getCollection<ldmx::TrigScintCluster>(input_collection_, passName_)};
  HitNum_=0.0;
  for(const auto &hit : digis){
      if (hit.getPE()>minThr_) {        
        HitNum_+=1;
      }
      int ID = hit.getBarID();
      integratedCharge_[ID]+=hit.getQ();
      PENumber_[ID]+=hit.getPE();
  }
  CluNum_=0.0;
  for(const auto &cluster : clusters){
    CluNum_+=1;
    if(CluNum_>=5.0){overflow_=1.0;} 
    populated_[(int)CluNum_]=1;
    chargeInCluster_[(int)CluNum_]=cluster.getEnergy();
    hitNumInCluster_[(int)CluNum_]=cluster.getNHits();
    centroidOfCluster_[(int)CluNum_]=cluster.getCentroid();
  }      
  
  if(CluNum_>=5.0){overflow_=1.0;}else{
      populated_[(int)CluNum_]=1;
      chargeInCluster_[(int)CluNum_]=valE_;
      hitNumInCluster_[(int)CluNum_]=v_addedIndices_.size();
      centroidOfCluster_[(int)CluNum_]=centroid_;
      CluNum_+=1.0;
  }
  // over channels

  //if (trigScintClusters.size() > 0)
  //  event.add(output_collection_, trigScintClusters);
  tree_->Fill();
  // book keep which channels have already been added to a cluster

  return;
}


void dataShaper::onFileOpen() {
  ldmx_log(debug) << "Opening file!";

  return;
}

void dataShaper::onFileClose() {
  ldmx_log(debug) << "Closing file!";

  return;
}

void dataShaper::onProcessStart() {
  ldmx_log(debug) << "Process starts!";
  getHistoDirectory();
  tree_ = new TTree("Helper","Flattened and decoded raw WR data");
  tree_->Branch("IntCharge",integratedCharge_);
  tree_->Branch("PENumber",PENumber_);
  tree_->Branch("MaxTime",TimeSampleWMaxCharge_);
  tree_->Branch("CluNumb",&CluNum_);
  tree_->Branch("HitsAboveEThresh",&HitNum_);
  tree_->Branch("ClusterPop",populated_);
  tree_->Branch("chargeInCluster",chargeInCluster_);
  tree_->Branch("hitNumInCluster",hitNumInCluster_);
  tree_->Branch("centroidOfCluster",centroidOfCluster_);
  
  tree_->Branch("adcList",adcsList_);
  tree_->Branch("tdcList",tdcsList_);
  tree_->Branch("QList",QList_);
  tree_->Branch("Overflow",&overflow_);
  std::cout<<"Made all branch"<<std::endl;
  return;
}

void dataShaper::onProcessEnd() {
  ldmx_log(debug) << "Process ends!";

  return;
}

}  // namespace trigscint

DECLARE_PRODUCER_NS(trigscint, dataShaper);

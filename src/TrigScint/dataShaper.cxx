
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
  std::cout<<"Did I get here 1"<<std::endl;
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
    for(int J=0; J<30;J++){
        switch(I){
         case 0: 
            ChargeAtTimeForChan0_[J]=0;
         case 1: 
            ChargeAtTimeForChan1_[J]=0;
         case 2: 
            ChargeAtTimeForChan2_[J]=0;
         case 3: 
            ChargeAtTimeForChan3_[J]=0;
         case 4: 
            ChargeAtTimeForChan4_[J]=0;
         case 5: 
            ChargeAtTimeForChan5_[J]=0;
         case 6: 
            ChargeAtTimeForChan6_[J]=0;
         case 7: 
            ChargeAtTimeForChan7_[J]=0;
         case 8: 
            ChargeAtTimeForChan8_[J]=0;
         case 9: 
            ChargeAtTimeForChan9_[J]=0;
         case 10: 
            ChargeAtTimeForChan10_[J]=0;
         case 11: 
            ChargeAtTimeForChan11_[J]=0;
      }
    }
  }
  std::cout<<"Did I get here 2"<<std::endl;
  overflow_=0;
  //I COMMENTED THIS OUT 5/23/23
  /*const auto digisO{event.getCollection<trigscint::TrigScintQIEDigis>(input_collection2_, passName2_)};
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
      switch((int)digiO.getChanID()){
         case 0: 
            ChargeAtTimeForChan0_[I]+=QList_[I];
         case 1: 
            ChargeAtTimeForChan1_[I]+=QList_[I];
         case 2: 
            ChargeAtTimeForChan2_[I]+=QList_[I];
         case 3: 
            ChargeAtTimeForChan3_[I]+=QList_[I];
         case 4: 
            ChargeAtTimeForChan4_[I]+=QList_[I];
         case 5: 
            ChargeAtTimeForChan5_[I]+=QList_[I];
         case 6: 
            ChargeAtTimeForChan6_[I]+=QList_[I];
         case 7: 
            ChargeAtTimeForChan7_[I]+=QList_[I];
         case 8: 
            ChargeAtTimeForChan8_[I]+=QList_[I];
         case 9: 
            ChargeAtTimeForChan9_[I]+=QList_[I];
         case 10: 
            ChargeAtTimeForChan10_[I]+=QList_[I];
         case 11: 
            ChargeAtTimeForChan11_[I]+=QList_[I];
      }
      
    }
  }*/
  std::cout<<"Did I get here 3"<<std::endl;
  
  //I COMMENTED THIS OUT 5/23/23
  //const auto clusters{
  //    event.getCollection<ldmx::TrigScintCluster>(input_collection_, passName_)};
  

  std::cout<<"Did I get here 4"<<std::endl;
  std::cout<<input_collection3_<<std::endl;
  const auto digis{
      //4/16/2023 INCLUDE IF YOU ARE USING TESTBEAM HITS
      //was trigscint::TestBeamHit
  event.getCollection<ldmx::TrigScintHit>(input_collection3_, "sim")};//passName3_)};
  //TrigScintCluster
  std::cout<<"Did I get here 5"<<std::endl;
  HitNum_=0.0;
  for(const auto &hit : digis){
      //if(hit.getQualityFlag()<1){continue;}
      if (hit.getPE()>minThr_) {        
        HitNum_+=1;
      }
      int ID = hit.getBarID();
      
      //4/16/2023 INCLUDE IF YOU ARE USING TESTBEAM HITS
      //integratedCharge_[ID]+=hit.getQ();
      
      
      PENumber_[ID]+=hit.getPE();
      //std::cout<<"Channel ID: "<<ID<<", Integrated Charge: "<<integratedCharge_[ID]<<" "<<hit.getQ()<<std::endl;
  }
   std::cout<<"Did I get here 6"<<std::endl;
  /*for(int II=0;II<12;II++){
      HitEfficiency_[II]=0;
      if((II==0)or(II==11)or(II=7)or(II=9)){
        if((II==0)or(II==9)){
            if(PENumber_[II+1]>30){
                if(PENumber_[II]>30){
                    HitEfficiency_[II]=1;
                }else{
                    HitEfficiency_[II]=-1;
                
            }
        }else{
            if(PENumber_[II-1]>30){
                if(PENumber_[II]>30){
                    HitEfficiency_[II]=1;
                }else{
                    HitEfficiency_[II]=-1;
                }
            } 
        }
      }else{
          if((PENumber_[II+1]>30)and(PENumber_[II-1]>30)){
              if(PENumber_[II]>30){
                HitEfficiency_[II]=1;
              }else{
                HitEfficiency_[II]=-1;
              }
          }
      }
  }*/
  //for(int i=0;i<=12;i++){
  //  std::cout<<"Final Charge at "<<i<<" : "<<integratedCharge_[i]<<std::endl;
  //}
  
   
   //I COMMENTED THIS OUT 5/23/23
  /* std::cout<<"Did I get here 7"<<std::endl;
  CluNum_=0.0;
  for(const auto &cluster : clusters){
    CluNum_+=1;
    if(CluNum_>=5.0){overflow_=1.0;} 
    populated_[(int)CluNum_]=1;
    chargeInCluster_[(int)CluNum_]=cluster.getEnergy();
    hitNumInCluster_[(int)CluNum_]=cluster.getNHits();
    centroidOfCluster_[(int)CluNum_]=cluster.getCentroid();
  }      
  std::cout<<"Did I get here 8"<<std::endl;
  if(CluNum_>=5.0){overflow_=1.0;}else{
      populated_[(int)CluNum_]=1;
      chargeInCluster_[(int)CluNum_]=valE_;
      hitNumInCluster_[(int)CluNum_]=v_addedIndices_.size();
      centroidOfCluster_[(int)CluNum_]=centroid_;
      CluNum_+=1.0;
  }
  */
  
  
  // over channels
  std::cout<<"Did I get here 9"<<std::endl;
  //if (trigScintClusters.size() > 0)
  //  event.add(output_collection_, trigScintClusters);
  tree_->Fill();
  // book keep which channels have already been added to a cluster
  std::cout<<"Did I get here 10"<<std::endl; 
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
  tree_->Branch("IntCharge",integratedCharge_,"IntCharge[17]/F");
  tree_->Branch("PENumber",PENumber_,"PENumber[17]/F");
  tree_->Branch("MaxTime",TimeSampleWMaxCharge_,"MaxTime[17]/F");
  tree_->Branch("ChargeAtTimeForChan0",ChargeAtTimeForChan0_,"ChargeAtTimeForChan0[30]/F");
  tree_->Branch("ChargeAtTimeForChan1",ChargeAtTimeForChan1_,"ChargeAtTimeForChan1[30]/F");
  tree_->Branch("ChargeAtTimeForChan2",ChargeAtTimeForChan2_,"ChargeAtTimeForChan2[30]/F");
  tree_->Branch("ChargeAtTimeForChan3",ChargeAtTimeForChan3_,"ChargeAtTimeForChan3[30]/F");
  tree_->Branch("ChargeAtTimeForChan4",ChargeAtTimeForChan4_,"ChargeAtTimeForChan4[30]/F");
  tree_->Branch("ChargeAtTimeForChan5",ChargeAtTimeForChan5_,"ChargeAtTimeForChan5[30]/F");
  tree_->Branch("ChargeAtTimeForChan6",ChargeAtTimeForChan6_,"ChargeAtTimeForChan6[30]/F");
  tree_->Branch("ChargeAtTimeForChan7",ChargeAtTimeForChan7_,"ChargeAtTimeForChan7[30]/F");
  tree_->Branch("ChargeAtTimeForChan8",ChargeAtTimeForChan8_,"ChargeAtTimeForChan8[30]/F");
  tree_->Branch("ChargeAtTimeForChan9",ChargeAtTimeForChan9_,"ChargeAtTimeForChan9[30]/F");
  tree_->Branch("ChargeAtTimeForChan10",ChargeAtTimeForChan10_,"ChargeAtTimeForChan10[30]/F");
  tree_->Branch("ChargeAtTimeForChan11",ChargeAtTimeForChan11_,"ChargeAtTimeForChan11[30]/F"); 
  //tree_->Branch("HitEfficiency",HitEfficiency_,"HitEfficiency[12]/F");

  tree_->Branch("CluNumb",&CluNum_);
  tree_->Branch("HitsAboveEThresh",&HitNum_);
  tree_->Branch("ClusterPop",populated_,"ClusterPop[5]/F");
  tree_->Branch("chargeInCluster",chargeInCluster_,"chargeInCluster[5]/F");
  tree_->Branch("hitNumInCluster",hitNumInCluster_,"hitNumInCluster[5]/F");
  tree_->Branch("centroidOfCluster",centroidOfCluster_,"centroidOfCluster[5]/F");
  tree_->Branch("adcList",adcsList_,"adcList[30]/F");
  tree_->Branch("tdcList",tdcsList_,"tdcList[30]/F");
  tree_->Branch("QList",QList_,"QList[30]/F");
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

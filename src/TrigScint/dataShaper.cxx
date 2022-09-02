
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
  passName_ = ps.getParameter<std::string>("input_pass_name");
  passName2_ = ps.getParameter<std::string>("input_pass_name2");
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
  const auto digisO{event.getCollection<trigscint::TrigScintQIEDigis>(
      input_collection2_, passName2_)};
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
      event.getCollection<trigscint::TestBeamHit>(input_collection_, passName_)};

  if (verbose_) {
    ldmx_log(debug) << "Got digi collection " << input_collection_ << "_"
                    << passName_ << " with " << digis.size() << " entries ";
  }

  // TODO remove this once verified that the noise overlap bug is gone
  bool doDuplicate = true;

  // 1. store all the channel digi content in channel order
  auto iDigi{0}; 
  CluNum_=0.0;
  HitNum_=0.0;
  for (const auto &digi : digis) {
    // these are unordered hits, and this collection is zero-suppressed
    // map the index of the digi to the channel index
	ldmx_log(debug) << "Digi has PE count " << digi.getPE() 
					<< " and energy " << digi.getEnergy();

	if (doCleanHits_ && digi.getQualityFlag() ) {
	  ldmx_log(debug) << "Skipping hit with non-zero quality flag  " << digi.getQualityFlag();
	  continue;
	}
	
    if (digi.getPE() >
        minThr_) {  // cut on a min threshold (for a non-seeding hit to be added
                    // to seeded clusters) already here
      HitNum_+=1;
      int ID = digi.getBarID();
      integratedCharge_[ID]+=digi.getQ();
      PENumber_[ID]+=digi.getPE();
      //if(digi.getEnergy()>maxEnergy){
//	TimeSampleWMaxCharge_[ID]=digi.getTime();
//	maxEnergy=digi.getEnergy();
      //}
	  if ( ID > 100) { //test beam has some uninstrumented channels (could also consider setting these to 0 in hit producer)
		//ldmx_log(debug) << "Skipping channel with bar ID = " << ID <<" > " << maxChannelID_ << " (max instrumented nb)";
		continue;
	  }
      // first check if there is a (pure) noise hit at this channel,  and
      // replace it in that case. this is a protection against a problem that
      // shouldn't be there in the first place.
      if (doDuplicate && hitChannelMap_.find((ID)) != hitChannelMap_.end()) {
        int idx = ID;
        std::map<int, int>::iterator itr = hitChannelMap_.find(idx);
        double oldVal = digis.at(itr->second).getPE();
        if (verbose_) {
          ldmx_log(debug) << "Got duplicate digis for channel " << idx
                          << ", with already inserted value " << oldVal
                          << " and new " << digi.getPE();
        }
        if (digi.getPE() > oldVal) {
          hitChannelMap_.erase(itr->first);
          if (verbose_) {
            ldmx_log(debug)
                << "Skipped duplicate digi with smaller value for channel "
                << idx;
          }
        }
      }

      // don't add in late hits 
      if (digi.getTime() > padTime_ + timeTolerance_ )
		continue;

      hitChannelMap_.insert(std::pair<int, int>(ID, iDigi));
      // the channel number is the key, the digi list index is the value

      if (verbose_) {
        ldmx_log(debug) << "Mapping digi hit nb " << iDigi
                        << " with energy = " << digi.getEnergy()
                        << " MeV, nPE = " << digi.getPE() << " > " << minThr_
                        << " to key/channel " << ID;
      }
    }
    iDigi++;
  }

  // 2. now step through all the channels in the map and cluster the hits

  std::map<int, int>::iterator itr;

  // Create the container to hold the digitized trigger scintillator hits.
  //std::vector<ldmx::TrigScintCluster> trigScintClusters;

  // loop over channels
  for (itr = hitChannelMap_.begin(); itr != hitChannelMap_.end(); ++itr) {
    // this hit may have disappeared
    if (hitChannelMap_.find(itr->first) == hitChannelMap_.end()) {
      if (verbose_ > 1) {
        ldmx_log(debug) << "Attempting to use removed hit at channel "
                        << itr->first << "; skipping.";
      }
      continue;
    }

    // i don't like this but for now, erasing elements in the map leads, as it
    // turns out, to edge cases where i miss out on hits or run into
    // non-existing indices. so while what i do below means that i don't need to
    // erase hits, i'd rather find a way to do that and skip this book keeping:
    bool hasUsed = false;
    for (const auto &index : v_usedIndices_) {
      if (index == itr->first) {
        if (verbose_ > 1) {
          ldmx_log(warn) << "Attempting to re-use hit at channel " << itr->first
                         << "; skipping.";
        }
        hasUsed = true;
      }
    }
    if (hasUsed) continue;
    if (verbose_ > 1) {
      ldmx_log(debug) << "\t At hit with channel nb " << itr->first << ".";
    }

    if (hitChannelMap_.size() ==
        0)  // we removed them all..? shouldn't ever happen
    {
      if (verbose_)
        ldmx_log(warn) << "Time flies, and all clusters have already been "
                          "removed! Unclear how we even got here; interfering "
                          "here to get out of the loop. ";
      break;
    }

    trigscint::TestBeamHit digi = (trigscint::TestBeamHit)digis.at(itr->second);

    // skip all until hit a seed
    if (digi.getPE() >= seed_) {
      if (verbose_ > 1) {
        ldmx_log(debug) << "Seeding cluster with channel " << itr->first
                        << "; content " << digi.getPE();
      }

      // 1.  add seeding hit to cluster

      addHit(itr->first, (trigscint::TestBeamHit)digi);

      if (verbose_ > 1) {
        ldmx_log(debug) << "\t itr is pointing at hit with channel nb "
                        << itr->first << ".";
      }

      // ----- first look back one step

      // we have added the hit from the neighbouring channel to the list only if
      // it's above clustering threshold so no check needed now
      std::map<int, int>::iterator itrBack =
          hitChannelMap_.find(itr->first - 1);

      bool hasBacked = false;

      if (itrBack !=
          hitChannelMap_.end()) {  // there is an entry for the previous
                                   // channel, so it had content above threshold
        // but it wasn't enough to seed a cluster. so, unambiguous that it
        // should be added here because it's its only chance to get in.

        // need to check again for backwards hits
        hasUsed = false;
        for (const auto &index : v_usedIndices_) {
          if (index == itrBack->first) {
            if (verbose_ > 1) {
              ldmx_log(warn) << "Attempting to re-use hit at channel "
                             << itrBack->first << "; skipping.";
            }
            hasUsed = true;
          }
        }
        if (!hasUsed) {
          digi = (trigscint::TestBeamHit)digis.at(itrBack->second);

          // 2. add seed-1 to cluster
          addHit(itrBack->first, digi);
          hasBacked = true;

          if (verbose_ > 1) {
            ldmx_log(debug) << "Added -1 channel " << itrBack->first
                            << " to cluster; content " << digi.getPE();
            ldmx_log(debug) << "\t itr is pointing at hit with channel nb "
                            << itr->first << ".";
          }

        }  // if seed-1 wasn't used already
      }    // there exists a lower, unused neighbour

      // 3. check next: if seed+1 exists && seed +2 exists,
      //    3a. if seed-1 is in already, this is a case for a split, at seed. go
      //    directly to check on seed-2, don't add more here. 3b. else. addHit
      //    (seed+1) 3c. if seed+3 exists, this is a split, at seed+1. don't add
      //    more here. 3d. else addHit(seed+2)
      // 4. if seed+1 and !seed+2 --> go to addHit(seed+1)

      // --- now, step 3, 4: look ahead 1 step from seed

      if (v_addedIndices_.size() < maxWidth_) {
        // (in principle these don't need to be different iterators, but it
        // makes the logic easier to follow)
        std::map<int, int>::iterator itrNeighb =
            hitChannelMap_.find(itr->first + 1);
        if (itrNeighb !=
            hitChannelMap_.end()) {  // there is an entry for the next channel,
                                     // so it had content above threshold
          // seed+1 exists
          // check if there is sth in position seed+2
          if (hitChannelMap_.find(itrNeighb->first + 1) !=
              hitChannelMap_.end()) {  // a hit with that key exists, so seed+1
                                       // and seed+2 exist
            if (!hasBacked) {  // there is no seed-1 in the cluster. room for at
                               // least seed+1, and for seed+2 only if there is
                               // no seed+3
              // 3b
              digi = (trigscint::TestBeamHit)digis.at(itrNeighb->second);
              addHit(itrNeighb->first, digi);

              if (verbose_ > 1) {
                ldmx_log(debug)
                    << "No -1 hit. Added +1 channel " << itrNeighb->first
                    << " to cluster; content " << digi.getPE();
                ldmx_log(debug) << "\t itr is pointing at hit with channel nb "
                                << itr->first << ".";
              }

              if (v_addedIndices_.size() < maxWidth_) {
                if (hitChannelMap_.find(itrNeighb->first + 2) ==
                    hitChannelMap_
                        .end()) {  // no seed+3. also no seed-1. so add seed+2
                  // 3d.  add seed+2 to the cluster
                  itrNeighb = hitChannelMap_.find(itr->first + 2);
                  digi = (trigscint::TestBeamHit)digis.at(itrNeighb->second);
                  addHit(itrNeighb->first, digi);
                  if (verbose_ > 1) {
                    ldmx_log(debug)
                        << "No +3 hit. Added +2 channel " << itrNeighb->first
                        << " to cluster; content " << digi.getPE();
                    ldmx_log(debug)
                        << "\t itr is pointing at hit with channel nb "
                        << itr->first << ".";
                  }
                }

              }   // if no seed+3 --> added seed+2
            }     // if seed-1 wasn't added
          }       // if seed+2 exists. then already added seed+1.
          else {  // so: if not, then we need to add seed+1 here. (step 4)
            digi = (trigscint::TestBeamHit)digis.at(
                itrNeighb->second);  // itrNeighb hasn't moved since there was
                                     // no seed+2
            addHit(itrNeighb->first, digi);

            if (verbose_ > 1) {
              ldmx_log(debug)
                  << "Added +1 channel " << itrNeighb->first
                  << " as last channel to cluster; content " << digi.getPE();
              ldmx_log(debug) << "\t itr is pointing at hit with channel nb "
                              << itr->first << ".";
            }
          }
        }  // if seed+1 exists
        // 5. at this point, if clusterSize is 2 hits and seed+1 didn't exist,
        // we can afford to walk back one more step and add whatever junk was
        // there (we know it's not a seed)
        else if (hasBacked &&
                 hitChannelMap_.find(itrBack->first - 1) !=
                     hitChannelMap_
                         .end()) {  // seed-1 has been added, but not seed+1,
                                    // and there is a hit in seed-2
          itrBack = hitChannelMap_.find(itr->first - 2);
          digi = (trigscint::TestBeamHit)digis.at(itrBack->second);
          addHit(itrBack->first, digi);

          if (verbose_ > 1) {
            ldmx_log(debug) << "Added -2 channel " << itrBack->first
                            << " to cluster; content " << digi.getPE();
          }
          if (verbose_ > 1) {
            ldmx_log(debug) << "\t itr is pointing at hit with channel nb "
                            << itr->first << ".";
          }

        }  // check if add in seed -2

      }  // if adding another hit, going forward, was allowed

      // done adding hits to cluster. calculate centroid
      centroid_ /= val_;  // final weighting step: divide by total
      centroid_ -= 1;     // shift back to actual channel center

      //ldmx::TrigScintCluster cluster;

      if (verbose_ > 1) {
        ldmx_log(debug) << "Now have " << v_addedIndices_.size()
                        << " hits in the cluster ";
      }
      //cluster.setSeed(v_addedIndices_.at(0));
      //cluster.setIDs(v_addedIndices_);
      //cluster.setNHits(v_addedIndices_.size());
      //cluster.setCentroid(centroid_);
      //cluster.setEnergy(valE_);
      //cluster.setPE(val_);
      //cluster.setTime(time_ / val_);
      //cluster.setBeamEfrac(beamE_ / valE_);
      if(CluNum_>=5.0){overflow_=1.0;}else{
      populated_[(int)CluNum_]=1;
      chargeInCluster_[(int)CluNum_]=valE_;
      hitNumInCluster_[(int)CluNum_]=v_addedIndices_.size();
      centroidOfCluster_[(int)CluNum_]=centroid_;
      CluNum_+=1.0;
      }
      //trigScintClusters.push_back(cluster);

      //if (verbose_) cluster.Print();
      //std::cout<<"HELLLLLLLLLLLLLO"<<std::endl;
      centroid_ = 0;
      val_ = 0;
      valE_ = 0;
      beamE_ = 0;
      time_ = 0;
      v_addedIndices_.resize(
          0);  // book keep which channels have already been added to a cluster

      if (verbose_ > 1) {
        ldmx_log(debug)
            << "\t Finished processing of seeding hit with channel nb "
            << itr->first << ".";
      }

    }  // if content enough to seed a cluster

    if (hitChannelMap_.begin() == hitChannelMap_.end()) {
      if (verbose_)
        ldmx_log(warn) << "Time flies, and all clusters have already been "
                          "removed! Interfering here to get out of the loop. ";
      break;
    }
  }  // over channels

  //if (trigScintClusters.size() > 0)
  //  event.add(output_collection_, trigScintClusters);
  tree_->Fill();
  hitChannelMap_.clear();
  v_usedIndices_.resize(
      0);  // book keep which channels have already been added to a cluster

  return;
}

void dataShaper::addHit(uint idx, trigscint::TestBeamHit hit) {
  float ampl = hit.getPE();
  val_ += ampl;
  float energy = hit.getEnergy();
  valE_ += energy;

  centroid_ += (idx + 1) * ampl;  // need non-zero weight of channel 0. shifting
                                  // centroid back by 1 in the end
  // this number gets divided by val at the end
  v_addedIndices_.push_back(idx);

  beamE_ += hit.getBeamEfrac() * energy;
  if (hit.getTime() > -990.) {
    time_ += hit.getTime() * ampl;
  }

  v_usedIndices_.push_back(idx);
  /*    // not working properly, but i'd prefer this type of solution
        hitChannelMap_.erase( idx ) ;
        if (verbose_ > 1 ) {
        ldmx_log(debug) << "Removed used hit " << idx << " from list";
        }
        if ( hitChannelMap_.find( idx) != hitChannelMap_.end() )
        std::cerr << "----- WARNING! Hit still present in map after removal!! ";
  */
  if (verbose_ > 1) {
    ldmx_log(debug) << "   In addHit, adding hit at " << idx
                    << " with amplitude " << ampl
                    << ", updating cluster to current centroid "
                    << centroid_ / val_ - 1 << " and energy " << val_
                    << ". index vector now ends with "
                    << v_addedIndices_.back();
  }

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

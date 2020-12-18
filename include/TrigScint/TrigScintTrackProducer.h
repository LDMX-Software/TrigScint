
#ifndef EVENTPROC_TRIGSCINTTRACKPRODUCER_H
#define EVENTPROC_TRIGSCINTTRACKPRODUCER_H

// LDMX Framework
#include "TrigScint/Event/TrigScintCluster.h"
#include "TrigScint/Event/TrigScintTrack.h"

namespace trigscint {

/**
 * @class TrigScintTrackProducer
 * @brief making tracks from trigger scintillator clusters
 */
class TrigScintTrackProducer : public framework::Producer {
 public:
  TrigScintTrackProducer(const std::string &name, framework::Process&process)
      : Producer(name, process) {}

  virtual void configure(framework::config::Parameters &ps);

  virtual void produce(framework::Event &event);

  virtual void onFileOpen();

  virtual void onFileClose();

  virtual void onProcessStart();

  virtual void onProcessEnd();

 private:
  // collection of produced tracks
  std::vector<trigscint::event::TrigScintTrack> tracks_;

  // add a cluster to a track
  trigscint::event::TrigScintTrack makeTrack(std::vector<trigscint::event::TrigScintCluster> clusters);

  // maximum difference (in channel number space) between track seed and cluster
  // in the next pad tolerated to form a track
  double maxDelta_{0.};

  // producer specific verbosity
  int verbose_{0};

  // collection used to seed the tracks
  std::string seeding_collection_;

  // other cluster collections used in track making
  std::vector<std::string> input_collections_;

  // output collection (tracks)
  std::string output_collection_;

  // specific pass name to use for track making
  std::string passName_{""};

  // track centroid in units of channel nb (will not be content weighted)
  float centroid_{0.};

  // track residual in units of channel nb (will not be content weighted)
  float residual_{0.};
};
}  // namespace trigscint

#endif  // EVENTPROC_TRIGSCINTTRACKPRODUCER_H

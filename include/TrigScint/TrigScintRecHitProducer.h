/**
 * @file TrigScintRecHitProducer.h
 * @brief Class that builds recHits
 * @author Andrew Whitbeck, TTU
 */

#ifndef TRIGSCINT_TRIGSCINTDIGIPRODUCER_H
#define TRIGSCINT_TRIGSCINTDIGIPRODUCER_H

/*~~~~~~~~~~*/
/*   ROOT   */
/*~~~~~~~~~~*/
#include "TRandom3.h"
#include "TLinearFitter.h"
#include "TVectorD.h"

// LDMX
#include "DetDescr/TrigScintID.h"
#include "Recon/Event/EventConstants.h"
#include "Tools/NoiseGenerator.h"
#include "TrigScint/Event/TrigScintHit.h"
#include "TrigScint/Event/TrigScintQIEDigis.h"

/*~~~~~~~~~~~~~~~*/
/*   Framework   */
/*~~~~~~~~~~~~~~~*/
#include "Framework/Configure/Parameters.h"
#include "Framework/EventProcessor.h"

/*~~~~~~~~~~~*/
/* TrigScint */
/*~~~~~~~~~~~*/
#include "TrigScint/SimQIE.h"

namespace trigscint {

/**
 * @class TrigScintRecHitProducer
 * @brief Organizes digis into TrigScintHits, linearizes TDC
 * and ADC info, and converts amplitudes to PEs

 */
class TrigScintRecHitProducer : public framework::Producer {
 public:
  TrigScintRecHitProducer(const std::string& name, framework::Process& process);

  ~TrigScintRecHitProducer();

  /**
   * Callback for the processor to configure itself from the given set
   * of parameters.
   *
   * @param parameters ParameterSet for configuration.
   */
  void configure(framework::config::Parameters& parameters) final override;

  void produce(framework::Event& event);

 private:
  /**
   * Reconstruct true charge deposited in each time sample
   * @param adc array of adcs for give event, cell
   * @param tdc array of tdcs for give event, cell
   * @param sample_of_interest sample of interest
   */
  Double_t ChargeReconstruction(std::vector<int>adc
                                ,std::vector<int>tdc
                                ,int sample_of_interest=2);

  /// Class to set the verbosity level.
  // TODO: Make use of the global verbose parameter.
  bool verbose_{false};

  /// Name of the input collection containing the sim hits
  std::string inputCollection_;

  /// Name of the pass that the input collection is on (empty string means take
  /// any pass)
  std::string inputPassName_;

  /// Name of the output collection that will be used to stored the
  /// digitized trigger scintillator hits
  std::string outputCollection_;

  /// SiPM gain
  double gain_{1e6};

  /// QIE pedestal
  double pedestal_{6.0};

  /// Total MeV per MIP
  double mevPerMip_{1.40};

  /// Total number of photoelectrons per MIP
  double pePerMip_{13.5};

  /// Sample of interest
  int sample_of_interest_{2};

  /**
   * Energy Reconstruction option
   * 0 -> Default. (summation of charges)
   * 1 -> Pulse fitting over 5 ts
   */
  int En_Reco_Option_{0};

  /// Input pulse shape for fitting
  std::string input_pulse_shape_;

  /// Input pulse parameters for fitting
  std::vector<float> pulse_params_;
};

}  // namespace trigscint

#endif

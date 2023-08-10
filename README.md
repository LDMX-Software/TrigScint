# TrigScint

(An Update with the Files Relevant to running Hit Effiicencies)
In this branch I have put all the files needed to run Hit Efficiency metrics on MC and real data. Here they are as follows:

## TrigScint/{src/TrigScint/dataShaper.cxx,include/TrigScint/dataShaper.h}:
This was the old ntupplizer in the NTuppler branch but configured to get the properties from MC and Real Data we need for the hit efficiency metrics.

## TrigScint/util/ChangedPlotterModOverlay.C
Once one obtains the result of running datashaper over clustered hits, you can run this file through ChangedPlotterModOverlay.C. This ROOT Macro has two
methods (ChangedPlotterGeom and ChanegedPlotterMod) which plot the hit density per bar and the hit efficiency metric repsectively. They take on file name arguments, which means one can get either from any MC or data file one has in mind.

## TrigScint/util/cpyBash.h
This bash script was used to take test stand data at earlier steps in reconstruction and run it up to the datashaper. The file structre of events taken after April 1rst subtly changed, meaning this would need to be modified to run on later runs. I am working on this right now.

## OtherChanges
The main other changes I have implemented while running the Hit Efficiency script was to geometries. This occurs in Detectors/data/ldmx-hcal-prototype-v2.0 and other similarly named files. I will include this in a branch of ldmx-sw similarly named to this branch of TrigScint (i.e. NTupplerHitEff or something).

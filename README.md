# TrigScint


## Instructions for performing reconstruction on test beam data

### Step 1: post process raw file:
To make the event sizes more uniform in the raw data files, a post-processing step is necessary.  If you skip this step
the unpacking step will likely complain about events not having the correct size.  This step can be performed by running
the following command:
`python3 TrigScint/util/decode_2fibers_toRAW_fromBin.py -i ../ts/data/ldmx_captan_out_17-04-2022_10-51-06__169_reformat_30timeSamplesFrom0.dat -o test -n 30 -n 1000 -f 0 -p > test.log`
the `-i` argument is the input raw (binary) file.  `-n` is the number of time samples to include.  `-f` is the index of the first time sample to include.
`-p` mean don't preselect events.  This will produce a very large log file :( and a `.txt` file that will be the input to the
next step.  
### Step 2:  Unpack the raw data into LDMX-SW events
To convert the raw data format into LDMX-SW events to be processed within the LDMX-SW framework, you need to run the
unpacker on the post-processed raw data.  You can do this by running the following command:
`ldmx fire TrigScint/exampleConfigs/runRawUnpacker.py ../ts/data/ldmx_captan_out_17-04-2022_10-51-06__169_reformat_30timeSamplesFrom0.txt test_unpacked.root 24`
The first argument is the input file raw data file from step 1.  The second argument is the output root file that will be
created.  The third argument is the number of time samples to unpack. 
### Step 3: Digi producer
In this step we will decode the raw data into formatted in to QIE digis.  To perform this step, you should
run the following command: 
`ldmx fire TrigScint/exampleConfigs/runQIEDecode.py test_unpacked.root test_digi.root "raw" 24 TrigScint/data/channelMap_identity_16channels.txt`
The first argument is the input root file from step 2; the second is the output file name; the third is ????; the fourth
is the channel mapping file. 
### Step 4: Linearize hits
`ldmx fire TrigScint/exampleConfigs/runTSLinearizer.py test_digi.root 2`
### Step 5: Hit reconstruction
To perform this step, run the following command: 
`ldmx fire TrigScint/exampleConfigs/runTestBeamHitReco.py test_digi_linearize.root 14`
The first argument is the input root file and the second command is the location of the pulse within dataframe.  


#!/usr/bin/scl enable devtoolset-8 -- /bin/bash
HOMEDERP=/sdf/group/ldmx/users/rodwyer1/ldmxStuff/nfsDirectory/rodwyer1/ldmx-swTry2
source ${HOMEDERP}/ldmx-sw/scripts/ldmx-env.sh
COUNTER=1
for file in /sdf/group/ldmx/data/TS-data/test_stand_data/rootfiles/*_hits.root
do
    #cp ${file} "./file10420No"${COUNTER}".evio"
    echo ${file}
    if [ $COUNTER -lt 100 ]
    then
        cutfile=${file:55} 
        echo ${file}
        cp ${file} .
        ldmx fire ${HOMEDERP}/ldmx-sw/TrigScint/exampleConfigs/runTestBeamClustering.py ${cutfile} #&> log1${COUNTER}.txt
        wait $!
        rm ${cutfile}
        echo ${cutfile::-5}
        ldmx fire ${HOMEDERP}/ldmx-sw/TrigScint/exampleConfigs/dataShaper.py ${cutfile::-5}_clusters.root #&> log${COUNTER}.txt
        wait $!
        rm ${cutfile::-5}_clusters.root
        rm ${cutfile::-5}_clusters.log 
        rm ${cutfile::-5}_clusters_dataShaped.log  
        wait $!
        #wait $!
        #rm ${cutfile}
    fi
    COUNTER=$((COUNTER+1))
    #echo "./file10420No"${COUNTER}".evio"
done


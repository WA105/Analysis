#!/bin/bash

##########################################################################################################################
## Create a text file with the arguments for the condor_subit file for swiss_highway_data_2018_jun_24_reco_functions.py ##
##########################################################################################################################

#path to my files
PathToRecoFiles="/eos/experiment/wa105/offline/LArSoft/Data/Reco/2018_June_24/ROOT/recofull"
#path to the output directory where Condor will copy the files
OutputDirectory=""
#name of the .txt files holding the argument list for the Condor
OutputList="arguments.txt"

if [[ -f $OutputList ]];
  then
    rm $OutputList
  fi

touch $OutputList

for run in 840; #enter here the range of run numbers you want to save to the OutputList, 633..1199
do
  for file in $PathToRecoFiles/$run/*.root;
  do

      RecoFile=$file

		    if [[ -f $RecoFile ]];
		      then
		          echo -i $RecoFile -o $OutputDirectory >> $OutputList
	        fi

  done

done

echo "All done"

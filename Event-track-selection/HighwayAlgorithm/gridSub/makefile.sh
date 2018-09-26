#!/bin/bash

##########################################################################################################################
## Create a text file with the arguments of swiss_highway.py for the submit.sub file ##
##########################################################################################################################

#path to the input files
PathToRecoFiles="/eos/experiment/wa105/offline/LArSoft/Data/Reco/2018_June_24/ROOT/recofull"
#path to the output directory where HTCondor will copy the files
OutputDirectory=""
#name of the .txt files holding the argument list for submit.sub
OutputList="arguments.txt"

if [[ -f $OutputList ]];
  then
    rm $OutputList
  fi

touch $OutputList

for run in 840; #enter here the range of run numbers you want to process (every subrun of the specified runs will be processed)
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

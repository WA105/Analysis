#!/bin/bash

PathToRecoFiles="/eos/experiment/wa105/offline/LArSoft/MC/MC6/ROOT/recofull"
PathToG4Files="/eos/experiment/wa105/offline/LArSoft/MC/MC6/ROOT/g4detsim"

OutputList="arguments.txt"

rm $OutputList
touch $OutputList

for i in {0..500} #enter here the range of run numbers you want to save to the OutputList, 633..1199
do
    G4File=$PathToG4Files"/$i-G4Detsim-Parser.root"
    RecoFile=$PathToRecoFiles"/$i-RecoFull-Parser.root"

		if [[ -f $RecoFile && -f $G4File ]];
		then
		    echo -s $G4File -r $RecoFile -n 840 >> $OutputList
	  fi

done

echo "All done"

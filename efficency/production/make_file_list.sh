#!/bin/bash

export version=$1
export PathToRecoFiles="/eos/user/a/ascarpel/stitch/MCB${version}/"
#export PathToRecoFiles="/eos/experiment/wa105/offline/LArSoft/MC/MCB/MCB1/MuonsCN/ROOT/recofast/"
export PathToG4Files="/eos/experiment/wa105/offline/LArSoft/MC/MCB/MCB${version}/MuonsCN/ROOT/g4detsim/"

OutputList="arguments_MCB${version}.txt"

rm $OutputList
touch $OutputList

for i in {0..1000} #enter here the range of run numbers you want to save to the OutputList, 633..1199
do
    G4File=$PathToG4Files"/$i-G4Detsim-Parser.root"
    RecoFile=$PathToRecoFiles"/$i-RecoFast-Parser.root"
    #RecoFile=$PathToRecoFiles"/3x1x1dp_reco_${i}_parser.root"

		if [[ -f $RecoFile && -f $G4File ]];
		then
        echo 'processing file ' ${i}
		    echo -s $G4File -r $RecoFile -o recoEfficiency.root >> $OutputList
	  fi

done

echo "All done"

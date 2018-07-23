# Highway algorithm #

## Prerequisites ##
In order to perform efficently swiss_highway_data_2018_jun_24_reco_functions.py expects some packages:
  *numpy
  *matplotlib
  *scipy
  *root
  *root_numpy

## Run the algorithm ##
To run "Highway algorithm" just do:
```
python swiss_highway_data_2018_jun_24_reco_functions.py -i /path/to/my/file/file.root -o /my/output/directory/
```
The script expects an input file name in the format: run-subrun-RecoFull-Parser.root. For the moment the only reconstruction version supported is '2018_June_24'

## Run the algorithm on a SWAN session ##
In order to run the highway yourself, copy the files "swiss_highway.ipynb" and "swiss_highway_data_2018_jun_24_reco_functions.py" to your CERNBox. Then open the ipynb script in [SWAN](https://swan.cern.ch) and specify the runs you want the highway to run over in the "run_list".

The default settings (which you can all change in the .ipynb script) are:
minimum track length: 20 cm
small rectangle width: 3.5 cm
large rectangle width: 10 cm
minimum box length (in case of multiple boxes): 50 cm
maximum number of tracks in event: 10000000 (effectively infinity)

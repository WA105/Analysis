# Highway algorithm #

## Prerequisites ##
In order to perform efficently swiss_highway.py expects some packages:
  * [numpy]()
  * [scipy]()
  * [root]()
  * [root_numpy]()

## Run the algorithm ##
The following tutorial will consider only the case when swiss_highway.py is run on the lxplus machines at CERN.

### Prepare the environment ###

### Run the algorithm on lxplus ###
To run "Highway algorithm" just do:
```
python swiss_highway.py -i /path/to/my/file/file.root -o /my/output/directory/
```
The script expects an input file name in the format: run-subrun-RecoFull-Parser.root. For the moment the only reconstruction version supported is '2018_June_24'

### Run the algorithm on Tier0 ###


## Highway algorithm on a SWAN session ##
In order to run the highway yourself, copy the files "swiss_highway.ipynb" and "swiss_highway.py" to your CERNBox. Then open the ipynb script in [SWAN](https://swan.cern.ch) and specify the runs you want the highway to run over in the "run_list".

The default settings (which you can all change in the .ipynb script) are:
* minimum track length: 20 cm
* small rectangle width: 3.5 cm
* large rectangle width: 10 cm
* minimum box length (in case of multiple boxes): 50 cm
* maximum number of tracks in event: 10000000 (effectively infinity)

## Disclamer ##

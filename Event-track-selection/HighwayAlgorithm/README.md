# Highway algorithm #

## Prerequisites ##
In order to perform efficently, swiss_highway.py expects some packages:
  * [numpy](http://www.numpy.org/)
  * [scipy](https://www.scipy.org/)
  * [ROOT](https://root.cern.ch/pyroot) (should be automatically configured sourcing the latest root version on `/afs` - see later)
  * [root_numpy](http://scikit-hep.org/root_numpy/)
  * [python-dateutil](https://dateutil.readthedocs.io/en/stable/)

## Run the algorithm ##
The following tutorial is valid for the lxplus machines at CERN.

### Prepare the environment ###
Some of the packages necessary are unfortunately not available natively on lxplus, nor support the default version of python. However, it is possible to create a virtual enviroment where to configure all the necessary packages. If it is the first time you are accessing your virtual enviromenent you may probably have to install or reinstall some packages to solve version conflict errors. The procedure might be cumbersome. In order to create a virtual environment with the correct version of python type after logging into an lxplus machine at CERN: 
```
scl enable python27 bash
```
You can verify to be in the correct environment by just checking the current running version of python. If it is 2.7 or higher, then the environment has been properly configured.  
Once logged in, you need to source the binaries of a recent ( > 6.0 ) version of ROOT available on `/afs/`: <a href="https://root.cern.ch/content/release-60606" target="_blank"> https://root.cern.ch/content/release-60606 </a>. It is particular 
importat to have the correct version of ROOT configured in order to install root_numpy properly and be able to use ROOT within python. 
To install the packages you are missing you can use the `pip` command:
```
pip install --user packagename
```
similarly to disinstall a package
```
pip uninstall packagename
```
or to upgrade it
```
pip install --upgrade --user packagename
```
To exit this environment simply type `exit`. The correct ROOT version must be sourced at every relogin on the environment, while it is not necessary to reinstall the packages. 

### Run the algorithm on lxplus ###
To run "swiss_highway.py" just do:
```
python swiss_highway.py -i /path/to/my/file/file.root -o /my/output/directory/
```
The script expects an input file name in the format: run-subrun-RecoFull-Parser.root. For the moment the only reconstruction version supported is '2018_June_24'

### Run the algorithm as batch job ###
In the folder gridSub are available two scripts to run "swiss_highway.py" as batch job on the NP02 computing grid. The scripts are intended as examples, every uses must adapt them following the comments withing the scripts themselves. 
* **makefile.sh** will create a text file with the list of argument required by "swiss_highway.py". Each line of this file corresponds to a different jobs. 
* **submit.sub** is the submission file for the batch service.
In order to submit your jobs do:
```
condor_submit submit.sub
```
To fetch information about the status of the jobs:
```
condor_q username
```
For more information about the commands available see the [HTCondor manual for CERN users](http://batchdocs.web.cern.ch/batchdocs/local/quick.html)  

## Highway algorithm on a SWAN session ##
One can run the algorithm directly on `/eos/` using the [SWAN](https://swan.cern.ch) service. In order to do that copy the files "swiss_highway.ipynb" and "swiss_highway.py" to your CERNBox. Then open the ipynb script in SWAN and specify the runs you want the highway to run over in the "run_list".

The default settings (which you can all change in the .ipynb script) are:
* minimum track length: 20 cm
* small rectangle width: 3.5 cm
* large rectangle width: 10 cm
* minimum box length (in case of multiple boxes): 50 cm
* maximum number of tracks in event: 10000000 (effectively infinity)

## Contacts ##
For every question about the algorithm contact <a href="mailto:caspar.maria.schloesser@cern.ch" target="_blank"> caspar.maria.schloesser@cern.ch </a>. For questions about the job submission and the environment creation contact For any question contact: <a href="mailto:andrea.scarpelli@cern.ch" target="_blank"> andrea.scarpelli@cern.ch </a> 

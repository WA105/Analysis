# Highway algorithm .npy files

## Prerequisites ##
  * numpy

## Accessing the highway output data in python ##
If you are using python, you can use the .npy files instead of the .txt files to access the highway output data. To get the charge ratio for a particlular run and subrun open the file in python, do the following:

highway_output = numpy.load('box_charge_ratios_run840_subrun0.npy')

This will contain a list of dictionaries with the following keys: ['run', 'subrun', 'event', 'track', 'charge_ratios_view0', 'charge_ratios_view1']

The entries for 'run', 'subrun', 'event' and 'track' are integers, the entries for 'charge_ratios_viewX' are lists where the number of entries depends on the length of the track and the number of boxes used.

**Event and muon track selection**:
1. select events and tracks according to the output txt file from the [Highway track selection algorithm](../Event-track-selection/HighwayAlgorithm/)

2. additional cut on track-length: take only tracks with more than 50 cm

3. depending on the analysis remmber to remove track hits belonging to corner LEMs (operated at a lower gain). 

4. for computation of electron lifetime keep only tracks which travers anode and cathode

**ADC to charge converstion:** use 55 (ADC x ticks)/fC for view 0 and 67 (ADC x ticks) /fC for view 1.
This is based on charge injection in the range of 14-150 fC during cold operation of the detector. See presentation from July 6th-2018 on indico. 

**Electron lifetime**: use a common value of 7 ms. This is based on best fit to most of the data (see below) and presentation from June 29th-2018 on indico.
![alt text](Lifetime_all_runs.png)

**Monte Carlo**: 
See details of productions on the [Twiki](https://twiki.cern.ch/twiki/bin/view/Sandbox/MonteCarloSamples3x1x1)


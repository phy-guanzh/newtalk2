# Talk06 -- get your own plots in a control region

mainly talk about how to plot based on your samples.

## Preparation 

```export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_29
cd CMSSW_10_6_29/src
cmsenv
https://github.com/phy-guanzh/newtalk2.git
cd newtalk2/talk06plot
```
need to make sure that the samples path is correct

```
vim loopPlot.C 

#change samples files path from
vector<TString> pathToTrees = {"../Common/rootfiles/"};
#to 
vector<TString> pathToTrees = {"/eos/user/z/zguan/www/talk/"};
#or your own path
```
Details:

in EDBRHistoMaker.h, you can check the selections we used or we can call it Control Region phase space defination.
Can you find the common selections for data and MC?

which lepton selections are used in data/MC?


## Run
```commandline
source make.sh 

./execute_name
```

to check the plots, just enter fig-output_18 or any year you want to check.


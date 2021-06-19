import FWCore.ParameterSet.Config as cms

from SVJ.Production.svjHelper import svjHelper

from collections import OrderedDict
from copy import deepcopy
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import numpy as np

# implementation of recursive loop over any number of dimensions
# creates grid of all possible combinations of parameter values
def varyAll(pos,paramlist,sig,sigs):
    param = paramlist[pos][0]
    vals = paramlist[pos][1]
    for v in vals:
        stmp = sig[:]+[v]
        # check if last param
        if pos+1==len(paramlist):
            sigs.add(tuple(stmp))
        else:
            varyAll(pos+1,paramlist,stmp,sigs)

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-y","--year", dest="year", type=int, default=2018, help="which year to simulate (specifies generator tune)")
parser.add_argument("-Z","--mZprime", dest="mZprime", type=float, default=[3000], nargs='+', help="mZprime value(s)")
parser.add_argument("-D","--mDark", dest="mDark", type=float, default=[20], nargs='+', help="mDark value(s)")
parser.add_argument("-R","--rinv", dest="rinv", type=float, default=[0.3], nargs='+', help="rinv value(s)")
parser.add_argument("-A","--alpha", dest="alpha", type=str, default=["peak"], nargs='+', help="alpha value(s)")
parser.add_argument("-n","--name", dest="name", type=str, default="Scan", help="custom name for scan")
args = parser.parse_args()

# specification of tunes for each year
if args.year==2016:
    tune_loc = "Configuration.Generator.Pythia8CUEP8M1Settings_cfi"
    tune_block = "pythia8CUEP8M1SettingsBlock"
    tune_suff = "TuneCUETP8M1_13TeV_pythia8"
    gen_tag = "cms.InputTag('generator')"
elif args.year==2017 or args.year==2018:
    tune_loc = "Configuration.Generator.MCTunes2017.PythiaCP2Settings_cfi"
    tune_block = "pythia8CP2SettingsBlock"
    tune_suff = "TuneCP2_13TeV_pythia8"
    gen_tag = "cms.InputTag('generator','unsmeared')"
else:
    parser.error("Unknown year: "+str(args.year))

# complete set of parameter values
params = OrderedDict([
    ("mZprime", args.mZprime),
    ("mDark", args.mDark),
    ("rinv", args.rinv),
    ("alpha", args.alpha),
])

# set to accumulate all scan points
sigs = set()

# full 4D scan
varyAll(0,list(params.iteritems()),[],sigs)

# format first part of output config
first_part = """
import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from {0} import * 
from Configuration.Generator.PSweightsPythia.PythiaPSweightsSettings_cfi import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.),
    RandomizedParameters = cms.VPSet(),
)
""".format(tune_loc)

# append process parameters for each model point
helper = svjHelper()
points = []
for point in sorted(sigs):
    mZprime = point[0]
    mDark = point[1]
    rinv = point[2]
    alpha = point[3]

    weight = 1.0
    # down-weight rinv=0 b/c all events pass filter
    if rinv==0.0:
        weight = 0.5

    helper.setModel("s",mZprime,mDark,rinv,alpha)
    
    pdict = {
        'weight': weight,
        'processParameters': helper.getPythiaSettings(),
        'name': helper.getOutName(outpre="SVJ"),
    }
    points.append(pdict)

# format last part of config (loop over all points)
last_part = """
for point in points:
    basePythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock, 
        {0},
        pythia8PSweightsSettingsBlock,
        processParameters = cms.vstring(point['processParameters']),
        parameterSets = cms.vstring(
            'pythia8CommonSettings',
            '{1}',
            'pythia8PSweightsSettings',
            'processParameters',
        )
    )

    generator.RandomizedParameters.append(
        cms.PSet(
            ConfigWeight = cms.double(point['weight']),
            ConfigDescription = cms.string(point['name']),
            PythiaParameters = basePythiaParameters,
        ),
    )

darkhadronZ2filter = cms.EDFilter("MCParticleModuloFilter",
    moduleLabel = {2},
    particleIDs = cms.vint32(51,53),
    multipleOf = cms.uint32(4),
    absID = cms.bool(True),
)

darkquarkFilter = cms.EDFilter("MCParticleModuloFilter",
    moduleLabel = {2},
    particleIDs = cms.vint32(4900101),
    multipleOf = cms.uint32(2),
    absID = cms.bool(True),
    min = cms.uint32(2),
    status = cms.int32(23),
)

ProductionFilterSequence = cms.Sequence(generator+darkhadronZ2filter+darkquarkFilter)
""".format(tune_block,tune_block.replace("Block",""),gen_tag)

with open("SVJ_{}_{}_{}_cff.py".format(args.name,args.year,tune_suff),'w') as ofile:
    ofile.write(first_part)
    ofile.write("\npoints = "+str(points)+"\n")
    ofile.write(last_part)

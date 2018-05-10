#!/usr/bin/python

print("Importing libraries")

import matplotlib
matplotlib.use('Agg') # To save files remotely.  Must be before importing matplotlib.pyplot or pylab!
from optparse import OptionParser
from collections import defaultdict
import os
import sys
import gzip
from tqdm import tqdm
import shutil
import time
# imports to run the model
import cobra
import PyDSTool as dst
import matplotlib.pyplot as plt
import pandas as pd
#import seaborn as sns
#sns.set_style('whitegrid')
import copy
import sys
sys.path.append('src')
import simulator
import utils


DFBA_DIR = 'dfba/data'
RESULTSDIR = 'outputs/'
# I tried simulating their models, but the reactions don't follow the standard nomenclature
#SBML_DIR = "%s/human-metabolic-atlas" % (DFBA_DIR)
#SBML_DIR = "%s/average-european-diet" % (DFBA_DIR)
# TODO make the Data Directory into an options
DATA_DIR = '/data/amogh-jalihal/Problem-Solving/'
#ALL_SBML_DIR = "/data/amogh-jalihal/Problem-Solving/average-european-diets/sbml"
DIET_DIR = "%s/diet-definitions" % (DFBA_DIR)

SBML_Dict = {
        'unconstrained':         '%s/unconstrained' % (DFBA_DIR),
        'western':               '%s/western-diet' % (DFBA_DIR),
        'european':              '%s/average-european-diet' % (DFBA_DIR),
        'high-fiber':            '%s/high-fiber-diet' % (DFBA_DIR),
#        'human-metabolic-atlas': '%s/human-metabolic-atlas' % (DFBA_DIR),
}

ALL_SBML_DIR = {
        'unconstrained':  '%s/unconstrained/sbml' % (DATA_DIR),
        'western':        '%s/western-diet/sbml' % (DATA_DIR),
        'european':       '%s/average-european-diets/sbml' % (DATA_DIR),
        'high-fiber':     '%s/high-fiber-diet/sbml' % (DATA_DIR),
#        'human-metabolic-atlas': '%s/human-metabolic-atlas' % (DFBA_DIR),
}

DietDict = {
    'HighFiber':     '%s/VMH_HighFiber.tsv' % (DIET_DIR),
    'HighProtein':   '%s/VMH_HighProtein.tsv' % (DIET_DIR),
    'Vegetarian':    '%s/VMH_Vegetarian.tsv' % (DIET_DIR),
    'Vegan':         '%s/VMH_Vegan.tsv' % (DIET_DIR),
    'Type2Diabetes': '%s/VMH_Type2Diabetes.tsv' % (DIET_DIR),
    'Unhealthy':     '%s/VMH_Unhealthy.tsv' % (DIET_DIR),
}


def main(species, abundances={}, out_pref=None, constraints='european', diet="HighFiber", max_iters=10, cobraonly=False, with_essential=False):

    #model_file_template = "%s/%%s.xml" % (SBML_DIR)
    # start with the default 0.01 for all species
    for s in species:
        if s not in abundances:
            abundances[s] = 0.01

    #count = 0.0
    SpeciesDict = {}
    for s in species:
        #model_file = model_file_template % (s)
        model_file = "%s/%s.xml" % (SBML_Dict[constraints], s)
        if not os.path.isfile(model_file):
            #print("%s doesn't exist. copying it from %s" % (model_file, ALL_SBML_DIR))
            print("%s doesn't exist. using it from %s" % (model_file, ALL_SBML_DIR[constraints]))
            model_file = "%s/%s.xml" % (ALL_SBML_DIR[constraints], s)
        SpeciesDict[s] = {'File': model_file, 'initAbundance': abundances[s], 'orig_name': s}
        #count += 1

    out_file_pref = "%s%dmodels-%s-%s-%diters-" % (out_pref, len(species), constraints, diet, max_iters)
    out_dir = os.path.dirname(out_file_pref)
    utils.checkDir(out_dir)

    simulate_models(species, SpeciesDict, diet=diet, out_file_pref=out_file_pref, max_iters=max_iters, cobraonly=cobraonly, with_essential=with_essential)


def simulate_models(species, SpeciesDict, constraints='european', diet="HighFiber", out_file_pref=None, max_iters=10,cobraonly=False, with_essential=False):
    """
    *species*: list of species to simulate together
    *abundances*: abundances for those species
    """

    Diet = pd.read_csv(DietDict[diet], sep='\t')
    Output, SpecDict, modDef = simulator.simulateCommunity(SpeciesDict, Diet, MaxIter=max_iters, cobraonly=cobraonly,with_essential=with_essential)

    species_abundances = defaultdict(list)
    for sp in SpecDict.keys():
        orig_name = SpecDict[sp]['orig_name']
        Name = SpeciesDict[sp]['Name']
        for P in Output:
            points = list(P[Name])
            # keep track of the abundance at the end of each iteration and write it to a file
            species_abundances[orig_name].append(points[-1])
    df = pd.DataFrame(species_abundances)
    log_file = out_file_pref + 'logs.tsv'
    
    print("Writing the abundances to %s" % (log_file))
    
    df.to_csv(log_file, sep='\t')

    simulator.plotBiomass(SpecDict, Output)
    plt.title("Constraints: %s, Diet: %s" % (constraints, diet))

    biomass_file = out_file_pref + 'biomass.pdf'
    print("writing %s" % (biomass_file))
    plt.savefig(biomass_file, bbox_inches='tight')
    plt.savefig(biomass_file.replace('.pdf', '.png'), bbox_inches='tight')
    
    if not cobraonly:
        fig = plt.figure(figsize=(20,15))
        simulator.plotImmuneResponse(SpecDict, Output)
        fig.suptitle("Constraints: %s, Diet: %s" % (constraints, diet))
        immune_response_file = out_file_pref + 'immune_response.pdf'
        print("writing %s" % (immune_response_file))
        plt.savefig(immune_response_file, bbox_inches='tight')
        plt.savefig(immune_response_file.replace('.pdf', '.png'), bbox_inches='tight')

        plt.figure()
        simulator.plotBiomassInMucosa(SpecDict, Output)
        plt.title("Diet: %s" % (diet))
        immune_response_file = out_file_pref + 'biomass-mucosa.pdf'
        print("writing %s" % (immune_response_file))
        plt.savefig(immune_response_file, bbox_inches='tight')
        plt.savefig(immune_response_file.replace('.pdf', '.png'), bbox_inches='tight')


def parse_args(args):
    usage = '%s [options]\n' % (args[0])
    parser = OptionParser(usage=usage)
    parser.add_option('','--species-file',type='string', 
                      help="List of species to simulate together. Not needed if abundances_file is given")
    parser.add_option('','--abundances-file',type='string', 
                      help="File containing the species name in the first column, and the abundance to use in the second")
    parser.add_option('','--constraints',type='string', default='european',
                      help="SBML constraints to use. Default='european'. Options are: '%s'" % ("', '".join(sorted(SBML_Dict.keys()))))
    parser.add_option('','--diet',type='string', default='HighFiber',
                      help="Diet to use Default='HighFiber'. Options are: '%s'" % ("', '".join(sorted(DietDict.keys()))))
    parser.add_option('','--out-pref',type='string', default="outputs/dfba/test",
                      help="output prefix for plot")
    parser.add_option('','--max-iters',type='int', default=10,
                      help="number of iterations to run the simulation")
    parser.add_option('','--with-essential',action="store_true", default=False,
                      help="Find exchange reactions essential to the growth of the microbes, and add those to the diet")
    parser.add_option('','--cobraonly',action="store_true", default=False,
                    help="Simulate the FBA models only, not the immune system.")
    parser.add_option('','--shigella', type='float',
                    help="Add Shigella flexneri (pathogen) to the list of species with a specified abundance (gdw). Usually 0.01 is the default.")
    parser.add_option('','--indv-sp',action="store_true", default=False,
                      help="Simulate each of the species individually")

    (opts, args) = parser.parse_args()

    return opts, args


if __name__ == "__main__":
    opts, args = parse_args(sys.argv)

    if opts.diet not in DietDict:
        print("Error: %s not an available diet. Options are: '%s'" % ("', '".join(sorted(DietDict.keys()))))
        sys.exit(1)
    if opts.constraints not in SBML_Dict:
        print("Error: %s not an available diet. Options are: '%s'" % ("', '".join(sorted(SBML_Dict.keys()))))
        sys.exit(1)

    #abundances = {s: 0.01 for s in species}
    species = None
    abundances = {}
    
    if opts.species_file is not None:
        df = pd.read_csv(opts.species_file, sep='\t')
        species = set(df['species'].tolist())
        print(species)
            
    elif opts.abundances_file is not None:
        df = pd.read_csv(opts.abundances_file, sep='\t')
        print(df)
        abundances = dict(zip(df['species'], df['abundance']))
        if species is None:
            species = set(df['species'].tolist())
    else:
        ListOfMostAbundantSpecies = [
            'Bacteroides_thetaiotaomicron_VPI_5482',
            'Bacteroides_sp_2_2_4',
            #'Parabacteroides_johnsonii_DSM_18315',
            #'Prevotella_oralis_ATCC_33269'
        ]
        species = set(ListOfMostAbundantSpecies)
        print("No species specified. Using defaults:")
        print(species)

    if opts.shigella is not None:
        species.add('Shigella_flexneri_2002017')
        abundances['Shigella_flexneri_2002017'] = opts.shigella

    if opts.indv_sp:
        for s in species:
            out_pref = opts.out_pref + "/individual/%s" % (s)
            main(set([s]), constraints=opts.constraints, diet=opts.diet,
                 out_pref=out_pref, max_iters=opts.max_iters, cobraonly=opts.cobraonly,
                 with_essential=opts.with_essential)
    else:
        #main(ListOfMostAbundantSpecies, opts.out_pref, abundances, max_iters=opts.max_iters)
        main(species, abundances=abundances, constraints=opts.constraints, diet=opts.diet,
             out_pref=opts.out_pref, max_iters=opts.max_iters, cobraonly=opts.cobraonly,
             with_essential=opts.with_essential)

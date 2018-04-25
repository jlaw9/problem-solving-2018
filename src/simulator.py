import cobra
import PyDSTool as dst
import pandas as pd
import matplotlib.pyplot as plt
import copy
import sys
import pdb
import time


def cleanupname(name):
    """
     The reaction names in the model files 
     don't have brackets or parentheses. I replaced
     those found in the mediaFluxes file.
     """
    # name = name.replace('[', '_LPAREN_')
    # name = name.replace(']', '_RPAREN_')
    # name = name.replace('(', '_LPAREN_')
    # name = name.replace(')', '_RPAREN_')

    name = name.replace('[', '__40__')
    name = name.replace(']', '__41__')
    name = name.replace('(', '__40__')
    name = name.replace(')', '__41__')


    # name = name.replace('[', '__91__')
    # name = name.replace(']', '__93__')
    # name = name.replace('(', '__91__')
    # name = name.replace(')', '__93__')
    return name

def defineDFBAModel(SpeciesDict , MediaDF):
    print("Defining Dynamical model... \n")
    ParDef = dict()
    VarDef = dict()
    ICS = dict()
    exchange_list = []
    mediaDerivedComponents = {}
    for i, row in MediaDF.iterrows():
        N = cleanupname(row.Reaction)
        mediaDerivedComponents[N] = row['Flux Value'] / (24.0*60.0) # Per minute

    for species in SpeciesDict.keys():
        print("\nReading species " + str(species))
        SpeciesDict[species]['SpeciesModel'] = cobra.io.read_sbml_model(SpeciesDict[species]['File'])
        SpeciesDict[species]['OriginalLB'] = {r.id:r.lower_bound/10.0 for r in SpeciesDict[species]['SpeciesModel'].exchanges}
        # SpeciesDict[species]['OriginalLB'] = {r.id:r.lower_bound for r in SpeciesDict[species]['SpeciesModel'].exchanges}
        SpeciesDict[species]['solution'] = SpeciesDict[species]['SpeciesModel'].optimize()
        SpeciesDict[species]['Name'] = SpeciesDict[species]['SpeciesModel'].name.split(' ')[0] + '_' \
        + SpeciesDict[species]['SpeciesModel'].name.split(' ')[1].replace('.','')
        exchange_list += SpeciesDict[species]['SpeciesModel'].exchanges
        Name=SpeciesDict[species]['Name']
        ParDef['mu' + '_' + Name] = SpeciesDict[species]['solution'].objective_value/60
        VarDef[Name] =  'mu_' + Name + ' * ' + Name + ' - ' + 'Dilution * ' + Name ### Biomass
        ICS[Name] = SpeciesDict[species]['initAbundance']

    ParDef['Dilution'] = 0.002

    all_exchanges = set()

    for ex in exchange_list:
        all_exchanges.add(ex.id)

    for rid in all_exchanges:
        VarDef[rid] = '- Dilution * ' + rid
        ICS[rid] = 0.1 #10.0

        if rid in mediaDerivedComponents.keys():
            ParDef[rid + '_influx'] = mediaDerivedComponents[rid]
            VarDef[rid] += ' + ' +  rid + '_influx'

        for species in SpeciesDict.keys():
            if 'h2o' in rid: # Check to see if a unique metabolite is represented only once
                print(species, rid)
            if rid in [species_r.id for species_r in SpeciesDict[species]['SpeciesModel'].exchanges]:
                Name = SpeciesDict[species]['Name']
                ParDef[rid + '_' + Name] = SpeciesDict[species]['solution'].fluxes[rid]/60.0
                VarDef[rid] += ' + ' +  rid + '_' + Name + ' * ' + Name

    ModelDef = dst.args(name='Comunity',
                        pars=ParDef,
                        varspecs=VarDef,
                        ics=ICS)
    ModelDS = dst.Vode_ODEsystem(ModelDef)
    print("Done!")
    return (SpeciesDict, ModelDef, ModelDS)

# Functions for model updates

def recomputeLowerBounds(SpeciesDict, PrevSteadyState, Kmax):
    for species in SpeciesDict.keys():
        for rid in [rxn.id for rxn in SpeciesDict[species]['SpeciesModel'].exchanges]:
            SpeciesDict[species]['SpeciesModel'].reactions.get_by_id(rid) \
                                                   .lower_bound = \
                                                                  SpeciesDict[species]['OriginalLB'][rid] \
                                                                  ,* PrevSteadyState[rid]/(Kmax+PrevSteadyState[rid])
    return SpeciesDict

def updateFluxParameters(SpeciesDict, ModelDS, PrevSteadyState):
    ParDef = {}
    ICS = {}
    for species in SpeciesDict:
        solution = SpeciesDict[species]['SpeciesModel'].optimize()
        Name = SpeciesDict[species]['Name']
        ParDef['mu_' + Name] = solution.objective_value/60.0
        ICS[Name] = PrevSteadyState[Name]
        exchanges = [r.id for r in SpeciesDict[species]['SpeciesModel'].exchanges]
        for rid in exchanges:
            # Control for cobra fl
            # Because very small non-zero solutions may come up despite 0 LB
            if abs(solution.fluxes[rid]/60.0) < 1e-12: 
                solution.fluxes[rid] = 0
            ParDef[rid + '_' + Name] = solution.fluxes[rid]/60.0
            ICS[rid] = PrevSteadyState[rid]
    ModelDS.set(pars=ParDef, ics=ICS)
    return ModelDS

def update(SpeciesDict, ModelDS, PrevSteadyState, Kmax):
    UpdatedSpeciesDict = recomputeLowerBounds(SpeciesDict,
                                              PrevSteadyState, Kmax)

    UpdatedDynamicModel = updateFluxParameters(UpdatedSpeciesDict,
                                               ModelDS,
                                               PrevSteadyState)
    return(UpdatedSpeciesDict, UpdatedDynamicModel)

def get_ss(PointSet):
    SSPoints={}
    for k in PointSet.keys():
        SSPoints[k]=PointSet[k][-1]
    return(SSPoints)

def checkNegativeMetabolites(PointSet, StoreNegatives):
    IndexStop = len(PointSet['t'])

    for variable in PointSet.keys():
        if any(PointSet[variable] < 0.0): # checking only final Tpoint, b/c monotonic
            varIndex = next((index for index,value in enumerate(PointSet[variable]) if value < 0), None)
            if varIndex < IndexStop:
                # Update the index for the first negative crossing
                IndexStop = varIndex

    if IndexStop < len(PointSet['t']) and IndexStop > 0:
        P_tilFirstNeg={}
        if len(PointSet[variable] > IndexStop+5):
            Extension = 5
        elif len(PointSet[variable] > IndexStop+2):
            Extension = 2
        else:
            Extension = 0
        for variable in PointSet.keys():
            P_tilFirstNeg[variable]=PointSet[variable][:IndexStop]
            if PointSet[variable][IndexStop + Extension] < 0.0:
                P_tilFirstNeg[variable][IndexStop - 1] = 0
                StoreNegatives.add(variable)
                print('\t' + str(variable) +  'is 0 at ' + str( PointSet['t'][IndexStop]))

        P_tilFirstNeg['t'] = PointSet['t'][:IndexStop]
        PointSet = P_tilFirstNeg
    return(PointSet,StoreNegatives)

def plotBiomass(SpeciesDict, AllPoints):
    TimePoints={}
    TimePoints['t'] =[]

    for P in AllPoints:
        TimePoints['t'] += list(P['t'])

    for sp in SpeciesDict.keys():
        Name = SpeciesDict[sp]['Name']
        TimePoints[Name] = []
        for P in AllPoints:
            TimePoints[Name]+=list(P[Name])

    for k in TimePoints.keys():
        if k != 't':
            plt.plot(TimePoints['t'], TimePoints[k], label = k)

        plt.xlabel('Time (minutes)')
    plt.ylabel('gdw')
    plt.legend(bbox_to_anchor=(1.2,1.2))


def plotMetabolites(AllPoints):
    TimePoints={}
    TimePoints['t'] =[]

    for P in AllPoints:
        TimePoints['t'] += list(P['t'])

    for v in P.keys():
        TimePoints[v] = []
        for P in AllPoints:
            TimePoints[v]+=list(P[v])

    for k in TimePoints.keys():
        if k != 't':
            plt.plot(TimePoints['t'], TimePoints[v])

        plt.xlabel('Time (minutes)')
    plt.ylabel('mmol')
    plt.legend()

def simulateCommunity(SpeciesDict, Diet, TEND=2000, MaxIter=200, Kmax=0.01, InitialValues = {}):
    """
    Simulates the microbial community.
    Arguments:
        1. Dictionary with paths to cobra with the following structure:
           SpeciesDict = {'sp_id':{'File':'PATH/TO/COBRA.MODEL','initAbundance':0.0}}
        2. Diet with medium input definition
    Returns:
        List of PointSet objects and updated SpeciesDict
    """
    if not InitialValues:
        SpeciesDict, Definition, ModelDS = defineDFBAModel(SpeciesDict, Diet)
        InitialValues = {k:[v] for (k,v) in Definition.ics.iteritems()}
    AllPoints = []
    StoreNegatives = set()
    P = InitialValues
    T0= 0
    TSPAN = 60
    IndexStop = 1 
    i=0

    clockstart = time.clock()
    while T0 < TEND and i < MaxIter:
        i+=1
        print(str(i) + ' ' + str(T0))
        SpeciesDict, ModelDS = update(SpeciesDict, ModelDS, get_ss(P), Kmax)

        if T0+TSPAN > TEND:
            TSPAN = TEND - T0

        ModelDS.set(tdata=[T0, T0 + TSPAN])
        P = ModelDS.compute('test').sample() 
        OldT = P['t'][-1]
        # Initialize
        P, StoreNegatives = checkNegativeMetabolites(P, StoreNegatives) 
        T0 = P['t'][-1]
        if OldT != T0:
            TSPAN = 1.0
        else:
            TSPAN = 60
        AllPoints.append(P)

    print("This took " + str(time.clock() - clockstart) + "s")
    return(AllPoints, SpeciesDict)



############################################################
# Sample Inputs

# SpeciesDict = {
#     'Sp1': {'File': '../dfba/data/average-european-diet/Bacteroides_sp_1_1_14.xml',
#             'initAbundance': 0.01},
#     'Sp2': {'File': '../dfba/data/average-european-diet/Ruminococcus_flavefaciens_FD_1.xml',
#             'initAbundance': 0.01},
#     'Sp3': {'File': '../dfba/data/average-european-diet/Lactobacillus_brevis_ATCC_367.xml',
#             'initAbundance': 0.01},
#     'Sp4': {'File': '../dfba/data/average-european-diet/Mycobacterium_avium_subsp_avium_ATCC_25291.xml',
#             'initAbundance': 0.01},
#     'Sp5': {'File': '../dfba/data/average-european-diet/Actinomyces_viscosus_C505.xml',
#             'initAbundance': 0.01},
#     'Sp6': {'File': '../dfba/data/average-european-diet/Exiguobacterium_aurantiacum_DSM_6208.xml',
#             'initAbundance': 0.01},
#     'Sp7': {'File': '../dfba/data/average-european-diet/Arcanobacterium_haemolyticum_DSM_20595.xml',
#             'initAbundance': 0.01},
#     'Sp8': {'File': '../dfba/data/average-european-diet/Streptococcus_intermedius_JTH08.xml',
#             'initAbundance': 0.01},
#     'Sp9': {'File': '../dfba/data/average-european-diet/Bifidobacterium_longum_infantis_ATCC_15697.xml',
#             'initAbundance': 0.01},
#     'Sp10': {'File': '../dfba/data/average-european-diet/Desulfovibrio_piger_ATCC_29098.xml',
#              'initAbundance': 0.01},
#     # 'Sp11': {'File': '../dfba/data/average-european-diet/Escherichia_coli_O157_H7_str_Sakai.xml',
#     #          'initAbundance': 0.01},

# } 

# Diet = pd.read_csv('../dfba/data/VMH_HighFiber.tsv', sep='\t')

#####################################################################

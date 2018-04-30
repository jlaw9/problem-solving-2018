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

def defineDFBAModel(SpeciesDict , MediaDF, cobraonly):
    print("Defining Dynamical model... \n")
    ParDef = dict()
    VarDef = dict()
    ICS = dict()
    exchange_list = []
    mediaDerivedComponents = {}
    for i, row in MediaDF.iterrows():
        N = cleanupname(row.Reaction)
        mediaDerivedComponents[N] = row['Flux Value'] / (24.0*60.0) # Per minute
    if not cobraonly:
        variable_dict = {
            # 'B_M':'B_MSource - (k_AD*B_M)/(k_3+B_M) - (k_AT*R_E*B_M)/(alpha_EM+R_E)-epsilon*B_M',
            # 'B':'max(0,epsilon*B_M-T) -k_5*P*B',
            'epsilon':'(epsilon_0-epsilon)/tau_p + f*P*(epsilon_max-epsilon)',
            'R_E':'(k_1*P)/(1+alpha_RE*I_E)-mu_RE*R_E +(1/(1+alpha_RE*I_E))',#*(a_1*B_M*T_I)/(gamma_1+B_M)',
            'I_E':'(k_IE*R_E)/(gamma_IE+R_E)-mu_IE*I_E', # +alpha_11*B_M',
            'P':'(k_PE*max(0,R_E-T_RE))/(1+gamma_PE*I_E)-mu_4*P ' # +(k_PM*B)/(gamma_12+B)
            
        }
        
        parameter_dict = {
            # 'B_MSource':2.2e6,#1.5e6,#
            'k_AD':1.5e6/60,#
            'k_3':6e6,
            'k_AT':0.03/60,
            'alpha_EM':0.18,
            'epsilon_0':0.1,
            'epsilon_max':0.21,
            'tau_p':24*60,
            'f':0.5,
            'a_1':0.1,#
            'gamma_1':5e6,
            'k_1':0.5/60,
            'alpha_RE':2,
            'mu_RE':0.1/60,
            'k_IE':50,#
            'gamma_IE':10,
            'alpha_11':0.1,
            'mu_IE':1,
            'T':1.1e6,
            'k_5':25,
            'k_PM':0.8,
            'gamma_12':1.2e6,
            'k_PE':0.002/60,
            'T_RE':0.65,
            'gamma_PE':1,
            'mu_4':0.05/60,
            'T_I':1e3 #????
        }
        #notebook
        initial_conditions = {
            #'B_M':1,
            'epsilon':0.1,
            'R_E':0,
            'I_E':0,
            #'B':1,
            'P':1
        }
        ParDef.update(parameter_dict)
        ICS.update(initial_conditions)
        VarDef.update(variable_dict)
        
        ParDef['K_diff'] = 1e-2
        
    for species in SpeciesDict.keys():
        print("\nReading species " + str(species))
        SpeciesDict[species]['SpeciesModel'] = cobra.io.read_sbml_model(SpeciesDict[species]['File'])
        SpeciesDict[species]['OriginalLB'] = {r.id:r.lower_bound/10.0 for r in SpeciesDict[species]['SpeciesModel'].exchanges}
        # SpeciesDict[species]['OriginalLB'] = {r.id:r.lower_bound for r in SpeciesDict[species]['SpeciesModel'].exchanges}
        SpeciesDict[species]['solution'] = SpeciesDict[species]['SpeciesModel'].optimize()
        SpeciesDict[species]['Name'] = SpeciesDict[species]['SpeciesModel'].name.split(' ')[0] + '_' \
                                       + SpeciesDict[species]['SpeciesModel'].name.split(' ')[1].replace('.','')
        SpeciesDict[species]['exchanges'] = [r.id for r in SpeciesDict[species]['SpeciesModel'].exchanges]
        exchange_list += SpeciesDict[species]['exchanges']
        
        Name=SpeciesDict[species]['Name']
        ICS[Name] = SpeciesDict[species]['initAbundance']
        ParDef['mu' + '_' + Name] = SpeciesDict[species]['solution'].objective_value/60
        VarDef[Name] =  'mu_' + Name + ' * ' + Name + ' - ' + 'Dilution * ' + Name
        if not cobraonly:
            VarDef[Name] += '- K_diff *' + Name ### Biomass, 20% diffuse to mucosa
            # 10^12 is a placeholder for mass/cell
            VarDef[Name + '_M'] = 'K_diff* ' + Name +'* 10^7 - (k_AD * ' + Name + '_M) / (k_3+' + Name + '_M)'\
                                  ' - (k_AT*R_E*' + Name + '_M)/(alpha_EM+R_E)-epsilon*' + Name + '_M'
            ICS[Name + '_M'] = 0.1*SpeciesDict[species]['initAbundance']*1e5# 0.0 # This should be non zero
            
            VarDef[Name + '_BT'] = 'max(0,epsilon*' + Name + '_M - T_'+ Name + ') - k_5 * P * ' + Name + '_BT'
            ICS[Name + '_BT'] = 0.0
            ParDef['T_' + Name] = 1.1e6
            
            VarDef['P'] += '+ (k_PM *' + Name + '_BT / (gamma_12 + ' + Name + '_BT ))'
            VarDef['R_E'] += '* (a_1 * '+ Name + '_M * T_I/(gamma_1 + '+Name +'_M))'
            VarDef['I_E'] += ' + alpha_11 * ' + Name + '_M'


    ParDef['Dilution'] = 0.002
  
    all_exchanges = set()
    
    for ex in exchange_list:
        all_exchanges.add(ex)
        
    for rid in all_exchanges:
        VarDef[rid] = '- Dilution * ' + rid
        ICS[rid] = 0.1 #10.0

        if rid in mediaDerivedComponents.keys():
            ParDef[rid + '_influx'] = mediaDerivedComponents[rid]
            VarDef[rid] += ' + ' +  rid + '_influx'
            
        for species in SpeciesDict.keys():
            if 'h2o' in rid: # Check to see if a unique metabolite is represented only once
                print(species, rid)
            if rid in SpeciesDict[species]['exchanges']:
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
        for rid in SpeciesDict[species]['exchanges']:
            SpeciesDict[species]['SpeciesModel'].reactions.get_by_id(rid) \
                                                          .lower_bound = \
                                                                         SpeciesDict[species]['OriginalLB'][rid] \
                                                                         * PrevSteadyState[rid]/(Kmax+PrevSteadyState[rid])
    return SpeciesDict
        
def updateFluxParameters(SpeciesDict, ModelDS, PrevSteadyState, cobraonly):
    ParDef = {}
    ICS = {}
    if not cobraonly:
        ICS['P'] = PrevSteadyState['P']
        ICS['R_E'] = PrevSteadyState['R_E']
        ICS['I_E'] = PrevSteadyState['I_E']
        ICS['epsilon'] = PrevSteadyState['epsilon']

    for species in SpeciesDict:
        solution = SpeciesDict[species]['SpeciesModel'].optimize()
        Name = SpeciesDict[species]['Name']
        ParDef['mu_' + Name] = solution.objective_value/60.0
        ICS[Name] = PrevSteadyState[Name]
        if not cobraonly:
            ICS[Name + '_M'] = PrevSteadyState[Name+'_M']
            ICS[Name + '_BT'] = PrevSteadyState[Name+'_BT']
            
        for rid in SpeciesDict[species]['exchanges']:
            # Control for cobra fl
            # Because very small non-zero solutions may come up despite 0 LB
            if abs(solution.fluxes[rid]/60.0) < 1e-12: 
                solution.fluxes[rid] = 0
            ParDef[rid + '_' + Name] = solution.fluxes[rid]/60.0
            ICS[rid] = PrevSteadyState[rid]
            ModelDS.set(pars=ParDef, ics=ICS)
    return ModelDS

def update(SpeciesDict, ModelDS, PrevSteadyState, Kmax, cobraonly):
    UpdatedSpeciesDict = recomputeLowerBounds(SpeciesDict,
                                              PrevSteadyState, Kmax)

    UpdatedDynamicModel = updateFluxParameters(UpdatedSpeciesDict,
                                               ModelDS,
                                               PrevSteadyState,
                                               cobraonly)
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
#        if len(PointSet[variable] > IndexStop+5):
#            Extension = 5
#        elif len(PointSet[variable] > IndexStop+2):
#            Extension = 2
#        else:
#            Extension = 0
        for variable in PointSet.keys():
            P_tilFirstNeg[variable]=PointSet[variable][:IndexStop]
            if PointSet[variable][IndexStop] < 0.0: #+ Extension]
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

def plotBiomassInMucosa(SpeciesDict, AllPoints):
    TimePoints={}
    TimePoints['t'] =[]
    
    for P in AllPoints:
        TimePoints['t'] += list(P['t'])
        
    for sp in SpeciesDict.keys():
        Name = SpeciesDict[sp]['Name']
        TimePoints[Name+'_M'] = []
        for P in AllPoints:
            TimePoints[Name+'_M']+=list(P[Name+'_M'])

    for k in TimePoints.keys():
        if k != 't':
            plt.plot(TimePoints['t'], TimePoints[k], label = k)

        plt.xlabel('Time (minutes)')
    plt.ylabel('gdw')
    plt.legend(bbox_to_anchor=(1.2,1.2))

def plotBiomassInBT(SpeciesDict, AllPoints):
    TimePoints={}
    TimePoints['t'] =[]

    for P in AllPoints:
        TimePoints['t'] += list(P['t'])

    for sp in SpeciesDict.keys():
        Name = SpeciesDict[sp]['Name']
        TimePoints[Name+'_BT'] = []
        for P in AllPoints:
            TimePoints[Name+'_BT']+=list(P[Name+'_BT'])

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

def simulateCommunity(SpeciesDict, Diet, TEND=2000, MaxIter=200, Kmax=0.01, InitialValues = {}, cobraonly=False):
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
        SpeciesDict, Definition, ModelDS = defineDFBAModel(SpeciesDict, Diet,cobraonly)
        InitialValues = {k:[v] for (k,v) in Definition.ics.iteritems()}
    AllPoints = []
    StoreNegatives = set()
    P = InitialValues
    T0 = 0
    TSPAN = 60
    IndexStop = 1 
    i = 0

    clockstart = time.clock()
    while T0 < TEND and i < MaxIter:
        i+=1
        print(str(i) + ' ' + str(T0))
        SpeciesDict, ModelDS = update(SpeciesDict, ModelDS, get_ss(P), Kmax, cobraonly)

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
    print("\t\tSTATISTICS")
    print("Variables\t\t"+str(len(Definition.varspecs.keys())))
    print("Parameters\t\t"+str(len(Definition.pars.keys())))

    return(AllPoints, SpeciesDict, Definition)



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

import os
import sys
import pandas as pd
RELPATH = os.getcwd()
SRCPATH = RELPATH + "/../../src/"
sys.path.append(SRCPATH)
import simulator
import matplotlib.pyplot as plt

ModelDefinitions = ['high-fiber-diet',
                    'western-diet',
                    #'average-european-diet',
                    'unconstrained']

ListOfMostAbundantSpecies = [
    'Bacteroides_thetaiotaomicron_VPI_5482',
    'Bacteroides_sp_2_2_4',
    'Parabacteroides_johnsonii_DSM_18315',
    'Prevotella_oralis_ATCC_33269',
    'Eubacterium_eligens_ATCC_27750',
    'Slackia_exigua_ATCC_700122',
    'Dorea_longicatena_DSM_13814',
    'Clostridium_bartlettii_DSM_16795',
    'Streptococcus_sp_I_P16',
    'Blautia_hydrogenotrophica_DSM_10507',
    'Brevundimonas_bacteroides_DSM_4726',
    'Clostridium_hylemonae_DSM_15053',
    'Sutterella_wadsworthensis_3_1_45B',
    'Enterobacteriaceae_bacterium_9_2_54FAA',
    'Bacillus_megaterium_DSM319',
    'Peptostreptococcus_stomatis_DSM_17678', # average-european-diet causes problems
    'Brachybacterium_paraconglomeratum_LC44',
    'Neisseria_elongata_subsp_glycolytica_ATCC_29315',
    'Rothia_aeria_F0474',
    'Staphylococcus_hominis_subsp_hominis_C80'
]
i = 0


DietDict = {
    'HighFiber': '../../dfba/data/diet-definitions/VMH_HighFiber.tsv',
    'HighProtein': '../../dfba/data/diet-definitions/VMH_HighProtein.tsv',
    'Vegetarian': '../../dfba/data/diet-definitions/VMH_Vegetarian.tsv',
    'Vegan': '../../dfba/data/diet-definitions/VMH_Vegan.tsv',
    'Type2Diabetes': '../../dfba/data/diet-definitions/VMH_Type2Diabetes.tsv',
    'Unhealthy': '../../dfba/data/diet-definitions/VMH_Unhealthy.tsv'
}

def characterizeSingleSpecies(ModelType, SingleSpeciesDict, DietDict):
    plt.figure()
    sp_id = SingleSpeciesDict.keys()[0]
    for diet in DietDict.keys():
        print('\t\t' + diet)
        DietDF = pd.read_csv(DietDict[diet], sep='\t')
        Output, SpecDict, Definition = simulator.simulateCommunity(SingleSpeciesDict, DietDF, MaxIter=20, cobraonly=True)
        T = []
        Biomass = []
        for P in Output:
            T += list(P['t'])
            Biomass += list(P[SpecDict[sp_id]['Name']])
        plt.plot(T, Biomass, label=diet)
    plt.legend()
    plt.title(SingleSpeciesDict[sp_id]['Name'])
    plt.xlabel('min')
    plt.ylabel('gdw')
    plt.savefig('./figs/'+ ModelType + '_'  + SingleSpeciesDict[sp_id]['Name'] + '.pdf', dpi=300)

for sp in ListOfMostAbundantSpecies:
    print(sp)
    for m in ModelDefinitions:
        print('\t' + m)
        characterizeSingleSpecies(m, {'1':{'File':RELPATH + '/../../dfba/data/' + m + '/' + sp + '.xml', 'initAbundance':0.01}}, DietDict)

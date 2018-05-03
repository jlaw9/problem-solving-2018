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
                    'unconstrained',
                    'average-european-diet',
]

# ListOfRandomAbundantSpecies = [
#     'Bacteroides_thetaiotaomicron_VPI_5482',
#     'Bacteroides_sp_2_2_4',
#     'Parabacteroides_johnsonii_DSM_18315',
#     'Prevotella_oralis_ATCC_33269',
#     'Eubacterium_eligens_ATCC_27750',
#     'Slackia_exigua_ATCC_700122',
#     'Dorea_longicatena_DSM_13814',
#     'Clostridium_bartlettii_DSM_16795',
#     'Streptococcus_sp_I_P16',
#     'Blautia_hydrogenotrophica_DSM_10507',
#     'Brevundimonas_bacteroides_DSM_4726',
#     'Clostridium_hylemonae_DSM_15053',
#     'Sutterella_wadsworthensis_3_1_45B',
#     'Enterobacteriaceae_bacterium_9_2_54FAA',
#     'Bacillus_megaterium_DSM319',
#     'Peptostreptococcus_stomatis_DSM_17678', # average-european-diet causes problems
#     'Brachybacterium_paraconglomeratum_LC44',
#     'Neisseria_elongata_subsp_glycolytica_ATCC_29315',
#     'Rothia_aeria_F0474',
#     'Staphylococcus_hominis_subsp_hominis_C80'
# ]

ListOfMostAbundantSpecies = [
'Bacteroides_thetaiotaomicron_VPI_5482',
'Bacteroides_uniformis_ATCC_8492',
'Bacteroides_vulgatus_ATCC_8482',
'Bacteroides_massiliensis_B846dnLKV334',
'Eubacterium_eligens_ATCC_27750',
'Campylobacter_showae_CSUNSWCD',
'Ruminococcus_bromii_L2_63',
'Ruminococcus_sp_5_1_39BFAA',
'Bifidobacterium_longum_longum_JCM_1217',
'Roseburia_intestinalis_L1_82',
'Akkermansia_muciniphila_ATCC_BAA_835',
'Eubacterium_rectale_M104_1',
'Sutterella_wadsworthensis_3_1_45B',
'Escherichia_coli_O157_H7_str_Sakai',
'Escherichia_albertii_KF1',
'Streptococcus_parasanguinis_ATCC_903',
'Dialister_succinatiphilus_YIT_11850',
'Haemophilus_parainfluenzae_T3T1',
'Dorea_formicigenerans_ATCC_27755',
'Faecalibacterium_prausnitzii_A2_165',
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
    plt.savefig('./figs_abundance/'+ ModelType + '_'  + SingleSpeciesDict[sp_id]['Name'] + '.pdf', dpi=300)

for sp in ListOfMostAbundantSpecies:
    print(sp)
    for m in ModelDefinitions:
        print('\t' + m)
        characterizeSingleSpecies(m, {'1':{'File':RELPATH + '/../../dfba/data/' + m + '/' + sp + '.xml', 'initAbundance':0.01}}, DietDict)

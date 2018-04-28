import pandas as pd
#RELPATH = os.getcwd()
#SRCPATH = RELPATH + "/../../src/"
#sys.path.append(SRCPATH)
import cobra
#import simulator
import matplotlib.pyplot as plt
from tqdm import tqdm

ListOfRandomSpecies = [
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
'Peptostreptococcus_stomatis_DSM_17678',
'Brachybacterium_paraconglomeratum_LC44',
'Neisseria_elongata_subsp_glycolytica_ATCC_29315',
'Rothia_aeria_F0474',
'Staphylococcus_hominis_subsp_hominis_C80',
]

for sp in ListOfRandomSpecies:
    print("Reading file " + sp)
    model = cobra.io.read_sbml_model("../../dfba/data/average-european-diet/" + sp + ".xml")
    AllExchanges = [r.id for r in model.exchanges]
    print("\tThe number of exchange reactions are " + str(len(AllExchanges)))
    EssentialGenes = [r.id for r in cobra.flux_analysis.variability.find_essential_genes(model)]
    print("\tThe number of essential genes are " + str(len(EssentialGenes)))
    Intersection = [rxn for rxn in EssentialGenes if rxn in AllExchanges]
    print("\tEssential exchanges for " + sp + " are:")
    for i in Intersection:
        print(i)

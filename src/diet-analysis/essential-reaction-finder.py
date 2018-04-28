import pandas as pd
import cobra
import copy
import matplotlib.pyplot as plt
from tqdm import tqdm


ListOfRandomSpecies = [
'Bacteroides_thetaiotaomicron_VPI_5482',
'Bacteroides_sp_2_2_4',
# 'Parabacteroides_johnsonii_DSM_18315',
# 'Prevotella_oralis_ATCC_33269',
# 'Eubacterium_eligens_ATCC_27750',
# 'Slackia_exigua_ATCC_700122',
# 'Dorea_longicatena_DSM_13814',
# 'Clostridium_bartlettii_DSM_16795',
# 'Streptococcus_sp_I_P16',
# 'Blautia_hydrogenotrophica_DSM_10507',
# 'Brevundimonas_bacteroides_DSM_4726',
# 'Clostridium_hylemonae_DSM_15053',
# 'Sutterella_wadsworthensis_3_1_45B',
# # # 'Enterobacteriaceae_bacterium_9_2_54FAA',
# # 'Bacillus_megaterium_DSM319',
# # 'Peptostreptococcus_stomatis_DSM_17678',
# # 'Brachybacterium_paraconglomeratum_LC44',
# 'Neisseria_elongata_subsp_glycolytica_ATCC_29315',
# 'Rothia_aeria_F0474',
# 'Staphylococcus_hominis_subsp_hominis_C80', 
]

betternames = [n.split('_')[0][0] + '.' + n.split('_')[1] for n in ListOfRandomSpecies]
stats = {}
types = ['western-diet', 'average-european-diet', 'high-fiber-diet', 'unconstrained']
for t in types:
    for sp in ListOfRandomSpecies:
        if sp not in stats.keys():
            stats[sp] = {}
        print(sp)

        model = cobra.io.read_sbml_model("../../dfba/data/"+ t +"/" + sp + ".xml")
        AllExchanges = [r.id for r in model.exchanges]
        solution = model.optimize()
        igr = solution.objective_value
        print("\tInitial growth rate:" + str(igr))
        DietSensitivity = {}
        
        count = 0
        for r in AllExchanges:
            if r not in DietSensitivity.keys():
                originallb = model.reactions.get_by_id(r).lower_bound
                model.reactions.get_by_id(r).lower_bound = 0.0
                s = model.optimize()
                gr = s.objective_value
                if abs(igr-gr)/igr > 0.01:
                    DietSensitivity[r] = abs(igr-gr)/igr
                    model.reactions.get_by_id(r).lower_bound = originallb
            count += 1

        solution = model.optimize()
        
        print("\tWith minimal influxes: " + str(solution.objective_value))
        print("\tThe minimal set contains " + str(len(DietSensitivity.keys())) + " metabolites")

        stats[sp][t]=len(DietSensitivity.keys())
        
for s in stats:
    print(s + '\t' + str(stats[s][types[0]]) + '\t' +  str(stats[s][types[1]]) + '\t' +  str(stats[s][types[2]]))

pos = -0.2
for t in types:
    plotthis = {sp:stats[sp][t] for sp in stats.keys()}
    plt.bar([r + pos for r in range(len(plotthis))], plotthis.values(),width=0.15, align='center', label=t)
    pos +=0.2
plt.xticks(range(len(plotthis)), betternames,rotation=30)
plt.legend()
plt.show()
    

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

ListOfColors=[
    'aqua','olive', 'black', 'chartreuse',
    'coral', 'crimson', 'darkgreen', 'green',
    'indigo', 'grey',
    'aqua','olive', 'black', 'chartreuse',
    'coral', 'crimson', 'darkgreen', 'green',
    'indigo', 'grey']

Styles = ['--','--','--','--','--','--','--','--','--','--',
'-','-','-','-','-','-','-','-','-','-']
DietDict = {
    'HighFiber': '../../dfba/data/diet-definitions/VMH_HighFiber.tsv',
    'HighProtein': '../../dfba/data/diet-definitions/VMH_HighProtein.tsv',
    'Vegetarian': '../../dfba/data/diet-definitions/VMH_Vegetarian.tsv',
    'Vegan': '../../dfba/data/diet-definitions/VMH_Vegan.tsv',
    'Type2Diabetes': '../../dfba/data/diet-definitions/VMH_Type2Diabetes.tsv',
    'Unhealthy': '../../dfba/data/diet-definitions/VMH_Unhealthy.tsv'
}

def plotBiomass(ModelType, SpeciesDict, AllPoints, diet):
    TimePoints={}
    TimePoints['t'] =[]
    plotspecs = {}
    for P in AllPoints:
        TimePoints['t'] += list(P['t'])

    for sp in SpeciesDict.keys():
        Name = SpeciesDict[sp]['Name']
        TimePoints[Name] = []
        plotspecs[Name] = {}
        plotspecs[Name]['color'] = 'xkcd:' + SpeciesDict[sp]['plotspec']['color']
        plotspecs[Name]['style'] = SpeciesDict[sp]['plotspec']['style']
        for P in AllPoints:
            TimePoints[Name]+=list(P[Name])

    for k in TimePoints.keys():
        if k != 't':
            plt.plot(TimePoints['t'], TimePoints[k], label = k,color=plotspecs[k]['color'], linestyle=plotspecs[k]['style'])

    plt.xlabel('Time (minutes)')
    plt.ylabel('gdw')
    plt.title(diet)
    plt.legend(bbox_to_anchor=(1.2,1.2))
    plt.savefig('./figs_community/'+ ModelType + '_'  + diet +'.pdf', bbox_inches='tight', dpi=300)

for m in ModelDefinitions:
    for diet in DietDict.keys():
        SpeciesDict = {}
        count = 0
        for s in ListOfMostAbundantSpecies:
            SpeciesDict['sp_' + str(count)] = {'File': RELPATH +  '/../../dfba/data/' + m + '/' + s + '.xml', 'initAbundance':0.01,'plotspec':{'color':ListOfColors[count],'style':Styles[count]}}
            count += 1
        Diet = pd.read_csv(RELPATH +'/' + DietDict[diet], sep='\t')
        Output, SpecDict, modDef = simulator.simulateCommunity(SpeciesDict, Diet, MaxIter=20, cobraonly=True)
        plotBiomass(m, SpecDict, Output,diet)
        del Diet
        del Output
        del SpecDict
        del modDef

# count = 0
i=0
x=[1,2,3,4]
for s in ListOfMostAbundantSpecies:
    plt.plot(x,x, label=s.split('_')[0]+ '_' + s.split('_')[1], color=ListOfColors[i], linestyle=Styles[i])
    i+=1
    
plt.legend()
plt.savefig('./figs_community/legend.svg',dpi=300)

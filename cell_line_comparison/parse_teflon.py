import os
import pysam
import matplotlib.pyplot as plt
import numpy as np
from statistics import mean,median,stdev
import pickle
import vcf as pyvcf
from copy import deepcopy
import sys

from functions import *

chroms = ['2L','2R','3L','3R','4','X']
autosomes = ['2L','2R','3L','3R']
X = ['X']


teflon_TEid_map_file = 'teflonRef.prep_TF/teflonRef.hier'

teflon_file = 'teflon_output/genotypes/our.genotypes.txt'

teflon_results_dir = 'teflon_output/genotypes'

output_dir = 'outputs'

global_vars_modules.global_idx_iterator = 0 # global_vars_modules is a dummy-class from my functions.py

############ SCRIPT START
teflon_TEid_map = {}
with open(teflon_TEid_map_file,'r') as f:
    for ln,line in enumerate(f):
        if ln == 0: continue
        line = line.strip('\n')
        line = line.split()
        if line[0].find('teid') == -1: continue
    
        TEid_raw,family,typee = line
        TEid = TEid_raw.split('#')[0]
        teflon_TEid_map[TEid] = [TEid_raw,family,typee]

SAMPLE_TEid_stats = {}
IDE_LTR_counts = {'pan':{},'reqRO':{}}
for file_ in sorted(os.listdir(teflon_results_dir)):
#if 1:
    #file_ = 'asd'
    if file_[0] == '.': continue
    teflon_file = teflon_results_dir+'/'+file_    
    
    rname_GTs = {}
    rname_AOs = {}
    rname_TEtypes = {}
    ref_TE_IDs = {}
    typees = {}
    fams = {}
    with open(teflon_file,'r') as f:
        for line in f:
            line = line.strip('\n')
            line = line.split()
            # presence reads = supp BP, absence = "where the alignment spans the breakpoints" (from paper)
            rname,BP_left,BP_right,family,typee,strand,ref_TE_ID,BP_left_hasSoftClip,BP_right_hasSoftClip,count_supp,count_absence,count_ambig,genotype,TE_id = line
            genotype = float(genotype)
            count_supp,count_absence,count_ambig = int(count_supp),int(count_absence),int(count_ambig)
            if BP_left.isdigit():        BP_left = int(BP_left)
            if BP_right.isdigit():       BP_right = int(BP_right)
            
            if not rname in chroms: continue
        
            init(SAMPLE_TEid_stats,TE_id,{})
            SAMPLE_TEid_stats[TE_id][file_] = {'count_supp':count_supp,'genotype':genotype,'rname':rname,'rcoords':[BP_left,BP_right],'ref_TE_ID':ref_TE_ID,'type':typee,'count_absence':count_absence}
            
            #if ref_TE_ID != '-': continue
            #if count_supp < 2: continue
            #if count_absence < 2: continue
            #if sum([count_supp,count_absence,count_ambig]) >= 15: continue
            if not typee == 'LTR': continue
        
            init(ref_TE_IDs,ref_TE_ID,0)
            ref_TE_IDs[ref_TE_ID] += 1
            
            init(typees,typee,0)
            typees[typee] += 1
            
            init(fams,family,0)
            fams[family] += 1
            
            #call_stats = teflon_TEid_map[TE_id]
            #call_type = call_stats[2]
            #init(rname_TEtypes,call_type,0)
            #rname_TEtypes[call_type] += 1
        
            init(rname_GTs,rname,[])
            rname_GTs[rname].append(genotype)
            
            init(rname_AOs,rname,[])
            rname_AOs[rname].append(count_supp)
            
    if 0:
        for rname,vals in rname_GTs.items():
            filt0 = vals
            filt1 = filterArr( filt0, 0, 'threshold' )
            filt2 = filterArr( filt1, 200, 'cutoff' )
            fig = plt.figure(figsize=(11,11))
            plt.hist( filt2 ,bins=200)
            plt.title(file_.split('.')[0] + ' || '+ rname+' || numVals='+str(len(filt2)))
            break
            
    typees_fracs = {}
    for typee,abu in typees.items():
        typees_fracs[typee] = round(abu/sum(typees.values()),4)
    fams_fracs = {}
    for fam,abu in fams.items():
        fams_fracs[fam] = round(abu/sum(fams.values()),4)
    
    if 0:
        print(file_.split('.')[0],sorted(fams_fracs.items(),key=lambda x: x[1], reverse=True)[:10])
        print('\n')
    
    if 1:
        print(file_.split('.')[0],typees['LTR'])
        print('\n')
        
        IDE_LTR_counts['reqRO'][file_.split('.')[0]] = typees['LTR']

###/

###
if 0 and 'IDE-dump S1,S2,S3 TEs for "filt fly TEs" from S2':
    Slines_TE_idx_map = {}
    for TEid,samples_data in SAMPLE_TEid_stats.items():
        for sample_raw,data in samples_data.items():
            sampleID = sample_raw.split('.')[0]
            if not sampleID in ('S2_our','SRR497713','SRR497721',): continue
            tmp_save = {}
            tmp_save.update(data)
            tmp_save['sample'] = sampleID
            tmp_save['TEid'] = TEid
            # assign new rcoords
            tmp_save['rcoords_raw'] = tmp_save['rcoords']
            tmp_save['rcoords'] = [None,None]
            if type(tmp_save['rcoords_raw'][0]) == int:     tmp_save['rcoords'][0] = tmp_save['rcoords_raw'][0]
            else:                                           tmp_save['rcoords'][0] = tmp_save['rcoords_raw'][-1]-1
            if type(tmp_save['rcoords_raw'][-1]) == int:    tmp_save['rcoords'][-1] = tmp_save['rcoords_raw'][-1]
            else:                                           tmp_save['rcoords'][-1] = tmp_save['rcoords_raw'][0]+1
            
            assignIdx(tmp_save)
            Slines_TE_idx_map[tmp_save['idx']] = tmp_save
    
    if 0 and 'dump?':
        with open(output_dir+'/'+'Slines_TE_idx_map.pickle','wb') as nf:
            pickle.dump(Slines_TE_idx_map,nf)
###/
    
### Find number of samples per TEid
if 0:
    INFO_numSamples_perTEID = []
    for TEid,samples_data in SAMPLE_TEid_stats.items():
        TEid_numSamples = 0
        
        IDE_save = False
        for sample,data in samples_data.items():
            if data['count_supp'] > 0 and data['ref_TE_ID'] == '-':
                if sample.find('SRR497727') != -1: IDE_save = True # ORE-R FEMALE SAMPLE
                TEid_numSamples += 1
        
        if not IDE_save: continue
        
        INFO_numSamples_perTEID.append(TEid_numSamples)
    
    filt0 = INFO_numSamples_perTEID
    filt1 = filterArr( filt0, 1, 'threshold' )
    filt2 = filterArr( filt1, 200, 'cutoff' )
    fig = plt.figure(figsize=(11,11))
    plt.hist( filt2 ,bins=200)
    plt.title('numVals='+str(len(filt2)))
###/

### Make phylo of TEs
#samples_to_keep = ['S2_our','SRR11076397','SRR497720','SRR497721','SRR612111']
#for i in range(len(samples_to_keep)): samples_to_keep[i] += '.genotypes.txt'

# 1. Parse rnames + coords, 2. Init bins, 3. find TEid ovlps to bins, 4. remove binns without TEids
rname_coords = {}
for TEid,sampleData in SAMPLE_TEid_stats.items():
    for sample,data in sampleData.items():
        rname = data['rname']
        for coord in data['rcoords']:
            if type(coord) == int:
                init(rname_coords,rname,set())
                rname_coords[rname].add(coord)

bin_size = 1000
bin_idx_map = {}
for rname,rcoords in rname_coords.items():
    for i in range(0,max(rcoords),bin_size):
        bstart = i
        bend = i+bin_size
        binn = {'rcoords':[bstart,bend],'rname':rname,'sample_TEids':{}}
        assignIdx(binn)
        bin_idx_map[binn['idx']] = binn
bin_idx_map_ROMD = rangeOverlaps_makeMappingDict(bin_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
for TEid,sampleData in SAMPLE_TEid_stats.items():
    for sample,data in sampleData.items():
        #if not sample in samples_to_keep: continue    
        rcoords = data['rcoords']
        # make rcoord span (for lookup) if other BP of TE is unknown
        if rcoords[0] == '-':       rcoords[0] = rcoords[-1]-1
        if rcoords[-1] == '-':      rcoords[-1] = rcoords[0]+1
        ovlps = rangeOverlaps_lookup(rcoords,bin_idx_map_ROMD,1000,map_lookup_key=rname)
        for oS,idx,oE in ovlps:
            binn = bin_idx_map[idx]
            init(binn['sample_TEids'],sample,[])
            binn['sample_TEids'][sample].append(data)
rm_idxs = []
for binn_idx,binn in bin_idx_map.items():
    if not binn['sample_TEids']:    rm_idxs.append(binn_idx)
for rm_idx in rm_idxs:
    del bin_idx_map[rm_idx]
bin_idx_map_ROMD = rangeOverlaps_makeMappingDict(bin_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
#/
# Make phylip file
sample_phylips = {} # sample -> arr_of_booleans, appending per observation
readsupp_CO = 2
for binn_idx,binn in sorted(bin_idx_map.items(),key=lambda x: (x[1]['rname'],x[1]['rcoords'][0])):
    # pre-check that any sample pass threshold, else we skip this binn
    skip_binn = True
    for sample,data in binn['sample_TEids'].items():
        saveState = 0
        for entry in data:
            if entry['count_supp'] >= readsupp_CO:
                skip_binn = False
    if skip_binn: continue
    #/
    for sample,data in binn['sample_TEids'].items():
        saveState = 0
        for entry in data:
            if entry['count_supp'] >= readsupp_CO and entry['ref_TE_ID'] == '-':
                saveState = 1
                break #only use first dupC call at bin
    
        init(sample_phylips,sample,[])
        sample_phylips[sample].append(saveState)

if 0 and 'dump to file?':
    nf = open(output_dir+'/'+'TEs.req'+str(readsupp_CO)+'.phylip','w')
    # Header
    nf.write(str(len(sample_phylips))+' '+str(len(list(sample_phylips.values())[0]))+'\n')
    # Data
    for sample,obs in sample_phylips.items():
        nf.write(sample+' '+''.join(map(str,obs))+'\n')
    nf.close()
    
    with open(output_dir+'/'+'TEs_per_sample.tsv','w') as nf:
        for sample,phylip in sample_phylips.items():
            sample = sample.replace('.genotypes.txt','')
            nf.write(sample+'\t'+str(phylip.count(1))+'\n')
    #/
#/
###/
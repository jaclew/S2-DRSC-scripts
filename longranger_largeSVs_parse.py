
import os
import pysam
import matplotlib.pyplot as plt
import numpy as np
from statistics import mean,median,stdev
from datetime import datetime,timedelta
import time
import pickle
import vcf as pyvcf
import math
import gzip
from copy import deepcopy
import sys

from functions import *

### INPUTS
chroms = ['2L','2R','3L','3R','4','X']
autosomes = ['2L','2R','3L','3R']
X = ['X']
chroms_to_plot = ('X','2L','2R','3L','3R','4','mitochondrion_genome')

vcf_file = 'longranger_wgs/S2/outs/large_svs.vcf.gz'
bedgraph_file = 'longranger_wgs/S2/outs/10X_LR.bedgraph'

GTF_file = 'dmel-all-chromosome-r6.12.gtf'
masking_file = 'dmel-all-chromosome-r6.12.fasta.out'

bam_file_for_refLens = 'aln.bam'

###/

output_dir = ''

### Define globals
global_vars_modules.global_idx_iterator = 0 # global_vars_modules is a dummy-class from my functions.py
ticToc = None

### Import ref stuff
# Ref lens
refLens = {}
vcf_fo = pyvcf.Reader(filename=vcf_file)
for ctg,ctg_data in vcf_fo.contigs.items():
    refLens[ctg] = ctg_data.length
refLens = importBAMlengths(bam_file_for_refLens)
#/ 
# Import masking
masking = importRepeatMasker(masking_file,rSorted=True)
masking_idx_map = {}
masking_types = {}
for rname in masking:
    for entry in masking[rname]:
        assignIdx(entry)
        masking_idx_map[entry['idx']] = entry
        masking_classi = entry['class']
        init(masking_types,masking_classi,0)
        masking_types[masking_classi] += 1
masking_idx_map_ROMD = rangeOverlaps_makeMappingDict(masking_idx_map,100,coordsKey='rcoords',sortByKey='rname')
#/

# Import GTFs
GTF_data = importGTF(GTF_file,rSorted=True,redundancyReduced=True)
GTF_types = {} #<-- NTS: all GTFs have a gene id!
GTF_idx_map = {}
gene_GTF_idx_map = {} # FBgnXXXXXXXX -> GTF_idx
for rname in GTF_data:
    for entry in GTF_data[rname]:
        assignIdx(entry)
        GTF_idx_map[entry['idx']] = entry
        
        init(GTF_types,entry['type'],0)
        GTF_types[entry['type']] += 1
        
        if entry['type'] == 'gene':
            gene_GTF_idx_map[entry['gene_id']] = entry['idx']
GTF_idx_map_ROMD = rangeOverlaps_makeMappingDict(GTF_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
GTF_feature_prioPurge = {'gene','mRNA','CDS','exon'}

#@ Compile GTF features in convenient format
# Sort by geneXgeneFeatures
gene_featureRanges = {}
for gtf_idx,gtf in GTF_idx_map.items():
    gene_id = gtf['gene_id']
    gene_feature = gtf['type']
    rstart,rend = gtf['rcoords']
    init(gene_featureRanges,gene_id,{})
    init(gene_featureRanges[gene_id],gene_feature,[])
    gene_featureRanges[gene_id][gene_feature].append([rstart,gtf_idx,rend])
#/
# Traverse gene functions, make ranges of those functions.
GFRO_idx_map = {} #geneFeaturesRangeeOvlpss
for gene,gene_features in gene_featureRanges.items():
    # Compile ovlpRanges of features
    features_ovlpRanges = {}
    for feature,ranges in gene_features.items():
        ovlps = computeOvlpRanges_wIdxs(ranges)
        tmp_rangee = []
        for rangee in ovlps:
            init(features_ovlpRanges,feature,[])
            
            # check if chain with previous
            if features_ovlpRanges[feature] and features_ovlpRanges[feature][-1][-1] == rangee[0]:
                features_ovlpRanges[feature][-1][-1] = rangee[-1]
                features_ovlpRanges[feature][-1][1].update(rangee[1])
            # else, save new
            else:
                features_ovlpRanges[feature].append([rangee[0],set(rangee[1]),rangee[-1]])
    #/
    # Save
    for feature,ovlpRanges in features_ovlpRanges.items():
        for ovlpRange in ovlpRanges:
            # Parse rname
            rname = None
            for any_entry in ovlpRanges:
                for any_idx in any_entry[1]:
                    rname = GTF_idx_map[any_idx]['rname']
                    break
            #/
            tmp_entry = {'rname':rname,'rcoords':[ovlpRange[0],ovlpRange[-1]],'gene_id':gene,'type':feature,
                         'members':ovlpRange[1],'len':ovlpRange[-1]-ovlpRange[0]}
            assignIdx(tmp_entry)
            GFRO_idx_map[tmp_entry['idx']] = tmp_entry
    #/
#/

# Compile which GFRO idxs exist for each gene
gene_GFRO_idxs = {} #gene_id -> GFRO idxs
for gfro_idx,gfro in GFRO_idx_map.items():
    gene = gfro['gene_id']
    init(gene_GFRO_idxs,gene,[])
    gene_GFRO_idxs[gene].append(gfro_idx)
#/    

GFRO_idx_map_ROMD = rangeOverlaps_makeMappingDict(GFRO_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
#@/
#/
###/

### Import bedgraph
def parse_masked_covFrac(rname,rstart,rend):
    # parse masked covFrac of bin
    ovlps = rangeOverlaps_lookup([rstart,rend],masking_idx_map_ROMD,100,map_lookup_key=rname)
    ranges = [[rstart,'base',rend]]
    for oS,idx,oE in ovlps:
        if oE==oS: continue #skip single base
        ranges.append([oS,idx,oE])
    ovlps2 = computeOvlpRanges_wIdxs(ranges)
    masked_numBPs = 0
    for rangee in ovlps2:
        if not 'base' in rangee[1]: continue
        if len(rangee[1]) == 1: continue #skip if only base
        masked_numBPs += rangee[-1]-rangee[0]
    masked_covFrac = masked_numBPs / (rend-rstart)
    #/
    return masked_covFrac
"""
rname_bedgraph_bins = {}
bedgraph = {}
bedgraph_binsize = 1000
BGraw_idx_map = {}
with open(bedgraph_file,'r') as f:
    for line in f:
        line = line.strip('\n')
        line = line.split()
        rname,rstart,rend,cov = line
        rstart,rend,cov = int(rstart),int(rend),float(cov)
        ## DONT CARE THAT BINS CAN END UP LARGER THAN X numBP. 0 impact in the end, just make it roughly
        binn = int(rstart/bedgraph_binsize)
        init(rname_bedgraph_bins,rname,{})
        init(rname_bedgraph_bins[rname],binn,[])
        rname_bedgraph_bins[rname][binn].append([rstart,cov,rend])
        
        # parse masked covFrac of bg
        masked_covFrac = parse_masked_covFrac(rname,rstart,rend)           

        tmp_BG = {'rname':rname,'rcoords':[rstart,rend],'cov':cov,'masked_covFrac':masked_covFrac}
        assignIdx(tmp_BG)
        BGraw_idx_map[tmp_BG['idx']] = tmp_BG
BGraw_idx_map_ROMD = rangeOverlaps_makeMappingDict(BGraw_idx_map,100,coordsKey='rcoords',sortByKey='rname')
#/
# parse bins
BG_idx_map = {}
for rname in rname_bedgraph_bins:
    for binn,entries in rname_bedgraph_bins[rname].items():
        binn_start = entries[0][0]
        binn_end = entries[-1][-1]
        binn_size = binn_end-binn_start
        # get numBp per covVal
        vals_parsed = {}
        for entry in entries:
            val = entry[1]
            init(vals_parsed,val,0)
            vals_parsed[val] += (entry[-1]-entry[0])
        #/
        # Calc cov representative for bin
        valMedian = None
        numBP_accu = 0
        numBP_tot = sum(vals_parsed.values())
        for val,numBP in sorted(vals_parsed.items(),key=lambda x: x[1],reverse=True):
            numBP_accu += numBP
            if numBP_accu >= numBP_tot/2:
                valMedian = val
                break
        
        tmp_save = {'rname':rname,'rcoords':[binn_start,binn_end],'cov':valMedian}
        assignIdx(tmp_save)
        BG_idx_map[tmp_save['idx']] = tmp_save
        #/
BG_idx_map_ROMD = rangeOverlaps_makeMappingDict(BG_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
#/
# Plot coverage per chrom
if 0:
    bin_step = 50 #50 for analysis here, 1000 for chromosomeOverview-bedgraph.
    smoothing_add = 0 #0 for analysis here, 4500 for smoothing in genomeOverview coverage track
    rname_BGcovs = {}
    for rname in chroms_to_plot:
        if rname in ('mitochondrion_genome',): continue
        for binnstart in range(0,refLens[rname]+bin_step,bin_step):
            rcoords = [binnstart,binnstart+bin_step]
            rcoords_scan = [rcoords[0]-smoothing_add , rcoords[-1] + smoothing_add]
            
            maskCovFrac = parse_masked_covFrac(rname,rcoords_scan[0],rcoords_scan[-1])           
            
            BG_ovlps = rangeOverlaps_lookup(rcoords_scan,BG_idx_map_ROMD,1000,map_lookup_key=rname)

            BG_ovlps_covs = []
            for oS,idx,oE in BG_ovlps:
                if oS == oE: continue
                masked_covFrac = parse_masked_covFrac(rname,oS,oE)   
                if masked_covFrac > 0: continue #skip repeatmasked
                BG = BG_idx_map[idx]
                BG_ovlps_covs += [BG['cov']]*(oE-oS)
            
            if not len(BG_ovlps_covs) >= bin_step*0.8: continue
        
            binn_cov = median(BG_ovlps_covs)
            
            if maskCovFrac >= 0.3 and binn_cov > 1000: continue #skip if masked AND high cov
            
            init(rname_BGcovs,rname,[])
            rname_BGcovs[rname].append([rcoords[0],binn_cov,rcoords[-1]])
    
    plt.rc('font',size=12)
    fig = plt.figure(figsize=(17,10))
    colmap = plt.cm.gist_rainbow(np.linspace(0,1,5)) #start,step,num_to_generate
    for enum,(rname,BGbins) in enumerate(sorted(rname_BGcovs.items(),key=lambda x: x[0])):
        if rname == '4': continue
        filt0 = [ BGbin[1] for BGbin in BGbins ]
        filt1 = filterArr( filt0, 1, 'threshold' )
        filt2 = filterArr( filt1, 400, 'cutoff' )
        if not filt2: continue
        numBins = 200 #int(max(filt2)-min(filt2))
        #plt.hist( filt2 ,bins=numBins,label=rname,histtype='step',stacked=True,fill=False,linewidth=1.5,alpha=0.5)
        plt.hist( filt2, weights = ((np.zeros_like(filt2)+1)/np.array(filt2).size), bins=numBins,label=rname,histtype='step',stacked=True,fill=False,linewidth=3,alpha=0.5)
    plt.locator_params(nbins=10)
    plt.xlim([0,400])
    plt.legend()
    plt.title(str(rname)+' numVals='+str(len(filt2)))
    if 1 and 'savefig?̈́': plt.savefig(output_dir+'tenX_CNbins.svg')
    
    
    if 0 and 'dump to file?':
        with open(bedgraph_file.replace('.bedgraph','')+'.'+str(bin_step)+'BP_bins.'+str(smoothing_add)+'winMarg.bedgraph','w') as nf:
            for rname in sorted(rname_BGcovs):
                for binn in sorted(rname_BGcovs[rname], key=lambda x: x[0]):
                    writeArr = [rname,binn[0],binn[-1],binn[1]]
                    nf.write('\t'.join(map(str,writeArr))+'\n')
    
    
                    
    
#/
###/
"""

### Import VCF
tenX_SV_events = {}
SV_len_thresh = 0
refLens = {}
INFO_numAlts = []
INFO_lens = {}
INFO_AFs = []
INFO_AFs_hap = []
INFO_pairs = []
INFO_filts = {}
IDE_types = {}
vcf_fo = pyvcf.Reader(filename=vcf_file)
for ctg,ctg_data in vcf_fo.contigs.items():
    refLens[ctg] = ctg_data.length
for entry in vcf_fo:
    chrom1 = entry.CHROM
    pos1 = entry.POS
    qual = entry.QUAL
    filt = entry.FILTER
    for idx,ALT in enumerate(entry.ALT):
        typee = ALT.type
        typee1 = entry.INFO['SVTYPE']
        try:
            typee2 = entry.INFO['SVTYPE2']
        except:
            typee2 = 'None'
        
        entry_id = entry.ID
        AF = entry.INFO['ALLELIC_FRAC']
        AF_hap = entry.INFO['HAP_ALLELIC_FRAC']
        pairs = entry.INFO['PAIRS']
        
        if typee == 'DUP:TANDEM': typee = 'DUP'
        
        ## Save if BND
        if typee == 'BND':
            event_id = entry.INFO['EVENT']
            rname2 = entry.alleles[1].chr
            rcoord2 = entry.alleles[1].pos
            LERI2,LERI1 = 1,1 # 1 mean extension is to left, -1 means to right.
            if not entry.alleles[1].orientation:        LERI2 = -1
            if not entry.alleles[1].remoteOrientation:  LERI1 = -1 # in "vcf alleles" the perspective is relative to pos2 (VA)
            tmp_event = {'rname':chrom1,'rcoords':[pos1,pos1],'AF':AF,'pairs':pairs,'type':typee,'entry_id':entry_id,
                         'qual':qual,'filt':filt,
                         'rname2':rname2,'rcoord2':rcoord2,'LERI1':LERI1,'LERI2':LERI2}
            init(tenX_SV_events,entry_id,[])
            tenX_SV_events[entry_id].append(tmp_event)
        ##/
        else:
        ## Save if SV call (DUP/INV/DEL?)
            pos2 = entry.INFO['END']
            tmp_event = {'rname':chrom1,'rcoords':[pos1,pos2],'AF':AF,'pairs':pairs,'type':typee,'entry_id':entry_id,
                         'qual':qual,'filt':filt,'type2':typee2}
            if entry_id in tenX_SV_events: sys.exit('Had a previous event ID!!')
            if 'EVENT' in entry.INFO and not entry.INFO['EVENT'] == event_id: sys.exit('Had a previous event ID[2]!!')
            tenX_SV_events[entry_id] = tmp_event
        ##/
        
        ## IDE parse
        typess = '||'.join(map(str,[typee,typee1,typee2]))
        init(IDE_types,typess,0)
        IDE_types[typess]+=1
        
        filts = '||'.join(sorted(filt))
        init(INFO_filts,filts,0)
        INFO_filts[filts] += 1
        
        if AF != None:
            INFO_AFs.append(AF)
        if AF_hap != None:
            INFO_AFs_hap.append(AF_hap)
        if pairs != None:
            INFO_pairs.append(pairs)
        ##/
        
if 0:
    filt0 = INFO_pairs
    filt1 = filterArr( filt0, 0, 'threshold' )
    filt2 = filterArr( filt1, 99999, 'cutoff' )
    fig = plt.figure(figsize=(11,11))
    plt.hist( filt2 ,bins=50)
    plt.title('INFO numVals='+str(len(filt2)))
    
## Post-parse events into SVs
SV_idx_map = {}
for event_id,event_entries in tenX_SV_events.items():
    tmp_SV = None
    # Parse if BND
    if type(event_entries) == list:
        if len(event_entries) == 1:
            tmp_SV = deepcopy(event_entries[0])
        elif len(event_entries) > 2:
            print('had more than 2 entries in BND call!!')
            sys.exit()
        else:
            tmp_SV = deepcopy(event_entries[0])            
            tmp_SV['type'] = 'TRA' #set to TRA if we have 2 linked BNDs (keep bnd if only 1)
    # parse if DEL/DUP/INV
    else:
        tmp_SV = deepcopy(event_entries)
    #/
    assignIdx(tmp_SV)
    SV_idx_map[tmp_SV['idx']] = tmp_SV
##/
###/IMPORT VCF

### Write for pygenometracks paint
if 0:
    output_dir_pygenometracks = output_dir+'/tenX_calls_filtFlag/'
    mkdir(output_dir_pygenometracks)
    
    rsqELACP_plotClassis = {'BND','DUP','DEL','INV','UNK'}
    nf_handles = {}
    for classi in rsqELACP_plotClassis:
        nf_handles[classi] = {'sameChrom':open(output_dir_pygenometracks+'/'+classi+'.sameChrom.arcs','w'),
                              'diffChrom':open(output_dir_pygenometracks+'/'+classi+'.diffChrom.arcs','w'),
                              'shortRearr':open(output_dir_pygenometracks+'/'+classi+'.short.arcs','w')}
        
    for SV_idx,SV in SV_idx_map.items():
        if SV['filt']: continue
    
        nf_classi = nf_handles[SV['type']] #assign NF handle
        
        # write for diffChrom
        if 'rname2' in SV:
            nf = nf_classi['diffChrom']
            
            writeArr = [SV['rname'],SV['rcoords'][0]-1000,SV['rcoords'][0]-1000,SV['rname'],SV['rcoords'][0]+1000,SV['rcoords'][0]+1000,1]
            nf.write('\t'.join(map(str,writeArr))+'\n')
            writeArr = [SV['rname2'],SV['rcoord2']-1000,SV['rcoord2']-1000,SV['rname2'],SV['rcoord2']+1000,SV['rcoord2']+1000,1]
            nf.write('\t'.join(map(str,writeArr))+'\n')
        else: #write for sameChrom
            # write as one if not so large
            if SV['rcoords'][-1]-SV['rcoords'][0] <= 100000:
                nf = nf_classi['shortRearr']
                writeArr = [SV['rname'],int(mean(SV['rcoords']))-2,int(mean(SV['rcoords']))-1,SV['rname'],int(mean(SV['rcoords']))+1,int(mean(SV['rcoords']))+2,1]
                nf.write('\t'.join(map(str,writeArr))+'\n')
            else:
                # write with arc if large
                nf = nf_classi['sameChrom']
                writeArr = [SV['rname'],SV['rcoords'][0],SV['rcoords'][0]+1,SV['rname'],SV['rcoords'][-1],SV['rcoords'][-1]+1,1]
                nf.write('\t'.join(map(str,writeArr))+'\n')
        #/
        
    for classi,chromStates_nfs in nf_handles.items():
        for chromState,nf in chromStates_nfs.items():
            nf.close()
###/

### Check SVs
numSVs_by_type = {}
numSVs_pass = 0
numSVs_pass_raw = 0
for SV_idx,SV in SV_idx_map.items():
    numSVs_pass_raw += 1
    #if SV['filt'] and 'LOWQ' in SV['filt']: continue
    numSVs_pass += 1
    
    typee = SV['type']
    
    init(numSVs_by_type,typee,0)
    numSVs_by_type[typee] += 1

if 0 and 'plot?':
    plotDict_barplot(numSVs_by_type)
    print('numTotal (vs pre-filt): '+str(numSVs_pass),numSVs_pass_raw)
    print('Pass filt [LOWQ tag]: '+str(numSVs_pass/numSVs_pass_raw))
    print('Number & percent annotated as BND: ', numSVs_by_type['BND'],(numSVs_by_type['BND']/numSVs_pass))
###/


     
### Check breakpoint repeat content
if 0:
    INFO_repma_stats = {}
    BP_dist_scan = 1000
    for SV_idx,SV in SV_idx_map.items():
        if SV['filt']: continue
    
        if SV['type'] == 'BND':
            if not SV['rname'] in chroms or not SV['rname2'] in chroms: continue
            rcoord1 = SV['rcoords'][0]
            rcoord2 = SV['rcoord2']
            scanCoords1 = [rcoord1,rcoord1+(BP_dist_scan*SV['LERI1'])]
            scanCoords2 = [rcoord2,rcoord2+(BP_dist_scan*SV['LERI2'])]
            
            #if SV['rname'] == '3L' and getRangeOvlp([rcoord1]*2,[3450000]) >= -100000: sys.exit()
            #if getRangeOvlp([rcoord1]*2,[3450000]) >= -100000: sys.exit()
            #if SV['entry_id'] == 'call_267_2': sys.exit()
            #if SV['entry_id'] == 'call_491_2': sys.exit()
            #if SV['entry_id'] == 'call_373_2': sys.exit()
            #if SV['entry_id'] == 'call_426_1': sys.exit()
            #if SV['entry_id'] == 'call_618_1': sys.exit()
            # Check masking
            
            #/
###/
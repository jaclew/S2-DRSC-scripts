
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as colMap
from statistics import mean,median,stdev
from datetime import datetime,timedelta
import pickle
import gzip
from copy import deepcopy
import sys

from functions import *

sys.path.append('libs/pyvenn') #path to pyvenn package
import venn

### INPUTS
chroms = ['2L','2R','3L','3R','4','X']
autosomes = ['2L','2R','3L','3R']
X = ['X']
chroms_to_plot = ('X','2L','2R','3L','3R','4','mitochondrion_genome')

GTF_file = 'dmel-all-chromosome-r6.12.gtf' #gene annotations
masking_file = 'RepeatMasked/dmel-all-chromosome-r6.12.fasta.out' #repeatmasker output file

vsDmel_bam_file = 'aln.bam' # FOR REF LENS
 
freec_file = 'barcoded_vsDmel.rmDups.sorted.bam_ratio.txt' #control-freec software output file
lee_file1 = 'Lee2014_BO.bedgraph' # copy-number from Lee et al (2014, after liftOver from dmel5->dmel6 using UCSC)
lee_file2 = 'Lee2014_DM.bedgraph'

kc_file = 'kc_dm6.bedgraph'
sg4_file = 'sg4_dm6.bedgraph'
###/

### OUTPUTS

###/

### Define globals
global_vars_modules.global_idx_iterator = 0 # global_vars_modules is a dummy-class from my functions.py
ticToc = None
###/

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

########## SCRIPT START

### Import ref stuff
# refLens
refLens = importBAMlengths(vsDmel_bam_file)

# Import masking
rname_maskings = importRepeatMasker(masking_file,rSorted=True)
masking_idx_map = {}
masking_types = {}
for rname in rname_maskings:
    for entry in rname_maskings[rname]:
        assignIdx(entry)
        masking_idx_map[entry['idx']] = entry
        masking_classi = entry['class']
        init(masking_types,masking_classi,0)
        masking_types[masking_classi] += 1
masking_idx_map_ROMD = rangeOverlaps_makeMappingDict(masking_idx_map,100,coordsKey='rcoords',sortByKey='rname')

masking_refInfo_perChrom = {}
masking_refInfo_2RL_centro = {'2L':{},'2R':{}}
maskClassis_lens = {}
maskClassis_blocks_rRanges = {}
mil = 1000000
for masking in masking_idx_map.values():
    rname = masking['rname']
    classi = masking['class']
    numBP = masking['rcoords'][-1] - masking['rcoords'][0]
    init(masking_refInfo_perChrom,rname,{})
    init(masking_refInfo_perChrom[rname],classi,0)
    masking_refInfo_perChrom[rname][classi] += numBP
    
    init(maskClassis_lens,classi,[])
    maskClassis_lens[classi].append(masking['len'])
    
    init(maskClassis_blocks_rRanges,rname,[])
    maskClassis_blocks_rRanges[rname].append([masking['rcoords'][0],masking['class'],masking['rcoords'][-1]])
    
    if rname == '2L' and getRangeOvlp(masking['rcoords'],[22*mil,99*mil]) >= 0:
        init(masking_refInfo_2RL_centro[rname],classi,[0,0])
        masking_refInfo_2RL_centro[rname][classi][0] += numBP
        masking_refInfo_2RL_centro[rname][classi][1] += numBP/mil
    if rname == '2R' and getRangeOvlp(masking['rcoords'],[0,4*mil]) >= 0:
        init(masking_refInfo_2RL_centro[rname],classi,[0,0])
        masking_refInfo_2RL_centro[rname][classi][0] += numBP
        masking_refInfo_2RL_centro[rname][classi][1] += numBP/mil

maskClassis_blockLens = {}
for rname,rranges in maskClassis_blocks_rRanges.items():
    ovlps = computeOvlpRanges_wIdxs(rranges)
    maskClassis_blocks = traverse_ranges_wChaining(ovlps,chain_dist=100)
    for maskClassi,blocks in maskClassis_blocks.items():
        for block in blocks:
            blockLen = block[-1]-block[0]
            init(maskClassis_blockLens,maskClassi,[])
            maskClassis_blockLens[maskClassi].append(blockLen)
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

### Import copynumbers
freecs_idx_map = {}
## S2, our
bin_size = 1000
with open(freec_file,'r') as f:
    for ln,line in enumerate(f): #Use ln to assign IDs
        if ln == 0: continue #skip header
        line = line.strip('\n')
        line = line.split()
        rname,start,ratio,ratio_median,CN = line
        rstart=int(start)
        rend = rstart + (bin_size-1)
        ratio=float(ratio)
        ratio_median=float(ratio_median)
        CN = int(CN)
        
        # Save freec data at LN
        tmp_freecE = {'rname':rname,'rcoords':[rstart,rend],'ratio':ratio,'ratio_median':ratio_median,'CN':CN,'source':'our'}
        assignIdx(tmp_freecE)
        
        freecs_idx_map[tmp_freecE['idx']] = tmp_freecE
##/
## Lee S2's
for enum,lee_file in enumerate( [lee_file1,lee_file2] ):
    with open(lee_file,'r') as f:
        for line in f:
            line = line.strip('\n')
            line = line.split('\t')
            rname = line[0].replace('chr','')
            rstart,rend,CN = list(map(int,line[1:]))
            
            tmp_freecE = {'rname':rname,'rcoords':[rstart,rend],'CN':CN,'source':'lee'+str(enum+1)}
            assignIdx(tmp_freecE)
            
            freecs_idx_map[tmp_freecE['idx']] = tmp_freecE
##/
## Lee sg4 & kc
for enum,lee_file in enumerate( [sg4_file,kc_file] ):
    with open(lee_file,'r') as f:
        for line in f:
            line = line.strip('\n')
            line = line.split('\t')
            rname = line[0].replace('chr','')
            rstart,rend,CN = list(map(int,line[1:]))
            
            tmp_freecE = {'rname':rname,'rcoords':[rstart,rend],'CN':CN,'source':lee_file.split('/')[-1]}
            assignIdx(tmp_freecE)
            
            freecs_idx_map[tmp_freecE['idx']] = tmp_freecE
##/
###/

### make ovlpranges and compile bins of S2 samples CNs (our vs lee1 vs lee2)
rRanges = {}
for freec_idx,freec in freecs_idx_map.items():
    rname = freec['rname']
    if not rname in chroms_to_plot: continue
    rstart,rend = freec['rcoords']
    init(rRanges,rname,[])
    rRanges[rname].append([rstart,freec_idx,rend])

CN_bins = {} #compile CNs between S2 samples
for rname,rangees in rRanges.items():
    ovlps = computeOvlpRanges_wIdxs(rangees)
    for rangee in ovlps:
        if (rangee[-1]-rangee[0]) < 500: continue #skip bullshit flank ovlps
        if len(rangee[1]) < 3: continue
    
        # calc masked covFrac
        masked_covFrac = parse_masked_covFrac(rname,rangee[0],rangee[-1])
        #/
        
        sample_CNs = {}
        for freec_idx in rangee[1]:
            freec = freecs_idx_map[freec_idx]
            sample_CNs[ freec['source'] ] = freec['CN']
            
        tmp_save = {'rname':rname,'rcoords':[rangee[0],rangee[-1]],'sample_CNs':sample_CNs,'masked_covFrac':masked_covFrac}
        assignIdx(tmp_save)
        
        CN_bins[tmp_save['idx']] = tmp_save
###/
### Calc stats about CN across samples
CNbin_stats = {'missingValInAny':0,'allSame':0,'leeDiffs':0,'ourDiff_leeSame':0,'ourSameEither_leeDiff':0,'ourSameLee1_leeDiff':0,'ourSameLee2_leeDiff':0,'allDiff':0}
CNbin_sampleStats_S2 = {} #S2 internal comp
CNbin_sampleStats_oLines = {} #Our S2 vs sg4+kc
TEMPLATE_CNbin_stats = deepcopy(CNbin_stats)
rname_CNbin_stats = {}
for CNbin in CN_bins.values():
    sample_CNs = CNbin['sample_CNs']
    
    our = sample_CNs['our']
    lee1 = sample_CNs['lee1']
    lee2 = sample_CNs['lee2']
    sg4 = sample_CNs['sg4_dm6.bedgraph']
    kc = sample_CNs['kc_dm6.bedgraph']
    
    # Save CNbin_id*CN at sample for pyvenn
    if CNbin['masked_covFrac'] < 0.5:
        for sample,CN in sample_CNs.items():
            # Save for S2 internal comp
            if sample in ('our','lee1','lee2',) and not -1 in (our,lee1,lee2):
                tmp_id = str(CNbin['idx'])+'_'+str(CN)
                init(CNbin_sampleStats_S2,sample,set())
                CNbin_sampleStats_S2[sample].add(tmp_id)
            #/
            # Save for our S2 vs sg4+kc
            if sample in ('our','sg4_dm6.bedgraph','kc_dm6.bedgraph',) and not -1 in (our,sg4,kc): #and not CNbin['rname'] in ('X','Y',):
                tmp_id = str(CNbin['idx'])+'_'+str(CN)
                init(CNbin_sampleStats_oLines,sample,set())
                CNbin_sampleStats_oLines[sample].add(tmp_id)
            #/
    #/

    if -1 in (our,lee1,lee2,):
        CNbin_stats['missingValInAny'] += 1
        continue #skip if either sample has no CN assigned at bin
    
    ## GENOME-WIDE
    # check if all CNs were called the same
    if min(sample_CNs.values()) == max(sample_CNs.values()):
        CNbin_stats['allSame'] += 1
    # if lee's were diff
    if lee1 != lee2:
        CNbin_stats['leeDiffs'] += 1
    # if our diff, but lee's are same
    if lee1 == lee2 and lee1 != our:
        CNbin_stats['ourDiff_leeSame'] += 1
    # if our is similar to any lee, but lee's diff
    if lee1 != lee2 and (our == lee1 or our == lee2):
        CNbin_stats['ourSameEither_leeDiff'] += 1
    ## if our similar to lee1, lee1 diff from lee2 + vice-versa
    if lee1 != lee2:
        if lee1 == our:
            CNbin_stats['ourSameLee1_leeDiff'] += 1
        elif lee2 == our:
            CNbin_stats['ourSameLee2_leeDiff'] += 1
    ##/
    # if all are diff
    if len(set(sample_CNs.values())) == 3:
        CNbin_stats['allDiff'] += 1
    ##/
    ## Chrom-wise
    # init chrom
    rname = CNbin['rname']
    if not rname in rname_CNbin_stats:      rname_CNbin_stats[rname] = deepcopy(TEMPLATE_CNbin_stats)
    #/
    # Run checks
    saveLoc = rname_CNbin_stats[rname]
    # check if all CNs were called the same
    if min(sample_CNs.values()) == max(sample_CNs.values()):
        saveLoc['allSame'] += 1
    # if lee's were diff
    if lee1 != lee2:
        saveLoc['leeDiffs'] += 1
    # if our diff, but lee's are same
    if lee1 == lee2 and lee1 != our:
        saveLoc['ourDiff_leeSame'] += 1
    # if our is similar to any lee, but lee's diff
    if lee1 != lee2 and (our == lee1 or our == lee2):
        saveLoc['ourSameEither_leeDiff'] += 1
    ## if our similar to lee1, lee1 diff from lee2 + vice-versa
    if lee1 != lee2:
        if lee1 == our:
            saveLoc['ourSameLee1_leeDiff'] += 1
        elif lee2 == our:
            saveLoc['ourSameLee2_leeDiff'] += 1
    ##/
    # if all are diff
    if len(set(sample_CNs.values())) == 3:
        saveLoc['allDiff'] += 1
    #/
    ##/
if 0:
    # pre-plot
    rname_CNbin_stats_PLOT = {}
    for rname,data in rname_CNbin_stats.items():
        for key,val in data.items():
            if key == 'missingValInAny': continue
            init(rname_CNbin_stats_PLOT,rname,{})
            rname_CNbin_stats_PLOT[rname][key] = (val / refLens[rname])*10**6
    #/
    plotEventDistribution(rname_CNbin_stats_PLOT,barScale=0.5,xticksrotate_deg=90)
###/
    
### Venn diagram plot shared CN segments
if 0:
    # S2 internal comp
    venn_legend = []
    venn_values = []
    for key,vals in CNbin_sampleStats_S2.items():
        venn_legend.append(key)
        venn_values.append(vals)
        
    venn_labels = venn.get_labels(venn_values,fill=['number']) #percent/number
    fig,ax = venn.venn3(venn_labels,names=venn_legend)
    fig.show()
    #/
    # our S2 vs sg4+kc
    venn_legend = []
    venn_values = []
    for key,vals in CNbin_sampleStats_oLines.items():
        venn_legend.append(key)
        venn_values.append(vals)
    
    venn_labels = venn.get_labels(venn_values,fill=['number'])
    fig,ax = venn.venn3(venn_labels,names=venn_legend)
    fig.show()
    #/
    # Dump to file for online proportional plot
    for key,vals in CNbin_sampleStats_S2.items():
        with open(key+'.list','w') as nf:
            nf.write('\n'.join(sorted(vals))+'\n')
    #/
###/


from functions import *

### INPUTS
chroms = ['2L','2R','3L','3R','4','X']
autosomes = ['2L','2R','3L','3R']
X = ['X']
chroms_to_plot = ('X','2L','2R','3L','3R','4','mitochondrion_genome')

vsDmel_bam_file = 'longread.bam' # FOR REF LENS

CN_file = 'Lee_2014/copy_numbers.txt'
###/

### OUTPUTS
output_dir = 'outputs'
###/

### Define globals
global_vars_modules.global_idx_iterator = 0 # global_vars_modules is a dummy-class from my functions.py
ticToc = None
###/
########## SCRIPT START

### Import copynumbers
sample_enums = {}
freecs_idx_map = {}
bin_size = 1000
with open(CN_file,'r') as f:
    for ln,line in enumerate(f):
        line = line.strip('\n')
        line = line.split('\t')
        
        # parse header
        if ln == 0:
            samples = line[2:]
            for enum,sample in enumerate(samples):
                sample = sample.replace('_CopyNumber','')
                sample = sample #+ '-' + str(enum)
                sample_enums[enum] = sample
            
            continue
        #/
        
        rname = line[0].replace('chr','')
        rstart = int(line[1])-1
        samples_vals = line[2:]
        for enum,sample_val in enumerate(samples_vals):
            sample = sample_enums[enum]
            if sample_val == 'NA': sample_val = -1
            
            sample_val = int(sample_val)
            rend = rstart + (bin_size-1)
            
            # Save freec data at LN
            tmp_freecE = {'rname':rname,'rcoords':[rstart,rend],'CN':sample_val,'source':sample}
            assignIdx(tmp_freecE)
            
            freecs_idx_map[tmp_freecE['idx']] = tmp_freecE
###/

### make ovlpranges and compile bins of samples CNs
rRanges = {}
for freec_idx,freec in freecs_idx_map.items():
    rname = freec['rname']
    if not rname in chroms_to_plot: continue
    rstart,rend = freec['rcoords']
    init(rRanges,rname,[])
    rRanges[rname].append([rstart,freec_idx,rend])

CN_bins = {} #compile CNs between samples
for rname,rangees in rRanges.items():
    ovlps = computeOvlpRanges_wIdxs(rangees)
    for rangee in ovlps:
        if len(rangee[1]) < 3: continue
        
        sample_CNs = {}
        for freec_idx in rangee[1]:
            freec = freecs_idx_map[freec_idx]
            sample_CNs[ freec['source'] ] = freec['CN']
            
        tmp_save = {'rname':rname,'rcoords':[rangee[0],rangee[-1]],'sample_CNs':sample_CNs}
        assignIdx(tmp_save)
        
        CN_bins[tmp_save['idx']] = tmp_save
###/

### Traverse CNbins, check if any is > X in CN, then mark as "is DUPed" (1/0, per sample). if none, then skip that bin.
# Per sample, find rname avg CN
sample_rname_CNs = {}
for freec in freecs_idx_map.values():
    rname = freec['rname']
    CN = freec['CN']
    sample = freec['source']
    if CN < 0: continue #skip if no assignment
    if CN > 8: continue
    init(sample_rname_CNs,sample,{})
    init(sample_rname_CNs[sample],rname,[])
    sample_rname_CNs[sample][rname].append(CN)

sample_rname_CNavgs = {}
for sample in sample_rname_CNs:
    tmp_A = []
    tmp_X = []
    for rname in sample_rname_CNs[sample]:
        if rname in autosomes+['4']:
            tmp_A += sample_rname_CNs[sample][rname]
        elif rname == 'X':
            tmp_X += sample_rname_CNs[sample][rname]
    
    for rname in sample_rname_CNs[sample]:
        CNavg = None
        if rname in autosomes+['4']:
            CNavg = int(median(tmp_A))
        elif rname == 'X':
            CNavg = int(median(tmp_X))
        init(sample_rname_CNavgs,sample,{})
        sample_rname_CNavgs[sample][rname] = CNavg
#/
# Phylip format is COL1 = sample, COL2 = boolean, per observation
CN_markDUP_ratio = 1.5
sample_phylips = {} # sample -> arr_of_booleans, appending per observation
for CNbin_idx,CNbin in CN_bins.items():
    # Skip if any sample does not have data in binn
    if min(CNbin['sample_CNs'].values()) == -1: continue
    for sample,CN in CNbin['sample_CNs'].items():
        init(sample_phylips,sample,[])
        
        # boolean check & save
        CN_boolean = 0
        if CN >= sample_rname_CNavgs[sample][CNbin['rname']]*CN_markDUP_ratio:        CN_boolean = 1
        
        sample_phylips[sample].append(CN_boolean)
###/

### Output
    # Dump phylip for tree
    if 0 and 'dump phylip?':
        nf = open(output_dir+'/'+'freecs_lee.phylip','w')
        # Header
        nf.write(str(len(sample_phylips))+' '+str(len(list(sample_phylips.values())[0]))+'\n')
        # Data
        for sample,obs in sample_phylips.items():
            nf.write(sample+' '+''.join(map(str,obs))+'\n')
        nf.close()
    
    if 0 and 'dump nums?':
        with open(output_dir+'/'+'freecs_per_sample.tsv','w') as nf:
            for sample,phylip in sample_phylips.items():
                nf.write(sample+'\t'+str(phylip.count(1))+'\n')
###/
import pysam
import pickle

from functions import *

masking_file = 'dmel-all-chromosome-r6.12.fasta.out'

vsDmel_bam_file = 'vsRepeatmasked.alenFilt2000.bam'
output_file = 'vsRepeatmasked.alenFilt2000.pickle'

global_vars_modules.global_idx_iterator = 0 # global_vars_modules is a dummy-class from my functions.py


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

## Import alns
qname_alns = {}
mapq_cutoff = 0
bam_fo = pysam.AlignmentFile(vsDmel_bam_file,'rb')
for bam_entry in bam_fo.fetch():
    if bam_entry.reference_id == -1: continue
    if bam_entry.infer_read_length() < 20000: continue
    if bam_entry.mapq <= mapq_cutoff: continue

    entry = parseBAMentry(bam_entry)
    qname = entry['qname']
    init(qname_alns,qname,[])
    
    assignIdx(entry)
    
    # Assign masking
    ovlps = rangeOverlaps_lookup(entry['rcoords'],masking_idx_map_ROMD,100,map_lookup_key=entry['rname'])
    
    # make rRanges, save base
    rRanges = []
    for oS,idx,oE in ovlps:
        if (oE-oS) <= 0: continue #skip low numBP ovlps/single-BP ovlps               
        rRanges.append( [oS,idx,oE] )
    ovlpRanges = computeOvlpRanges_wIdxs(rRanges)
    masked_numBP = 0
    for rangee in ovlpRanges:
        masked_numBP += rangee[-1]-rangee[0]
    masked_covFrac = masked_numBP / (entry['rcoords'][-1]-entry['rcoords'][0])
    entry['masked_covFrac'] = masked_covFrac
    #/
    if masked_covFrac > 1: print(masked_covFrac)
    
    qname_alns[qname].append(entry)
bam_fo.close()
##/

## Dump as rsqEs
with open(output_file,'wb') as nf:
    pickle.dump(qname_alns,nf)
##/
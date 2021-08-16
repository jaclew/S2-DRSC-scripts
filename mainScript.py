from multiprocessing import Pool
import os
import pysam
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as colMap
import numpy as np
from statistics import mean,median,stdev
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

GTF_file = 'dmel-all-chromosome-r6.12.gtf' #gene annotations
masking_file = 'RepeatMasked/dmel-all-chromosome-r6.12.fasta.out' #repeatmasker output file

vsDmel_bam_file = 'aln.nSorted.bam' # n-sorted file (sambamba sort -n <bam_file>)
vsDmel_coordSorted_bam_file = 'aln.sorted.bam' #coordinates/"standard" sort bam file (sambamba sort <bam_file>)
bg_file = 'aln.bedgraph' #bedgraph file of alignments (genomeCoverageBed -ibam <bam_file> -bg)
        

MM2_RAVA_saveFile = 'ava.paf.gz' # read-all-vs-all mappings (minimap2 -x ava-ont -r1000 <reads> <reads>)
    
if 0:
    read_seqs = importReadSeqsFasta('reads.fa') # read sequences, optional import
    ref_seqs = importReadSeqsFasta('dmel-all-chromosome-r6.12.fasta') # reference genome import (used to generate genome graph from variants)
###/

### OUTPUTS
MS_output = 'mainScript_output'
###/

### Define globals
global_vars_modules.global_idx_iterator = 0 # global_vars_modules is a dummy-class from my functions.py
ticToc = None
    
############### SCRIPT START
## Determine contained reads
lenFilt = 2000 #skip read if shorter than this numBP
HT_margin = 750 # require alns to be this close to head/tail of reads in ovlp.
numBP_ovlp_thresh = 3000 # skip alns below this numBP
numBP_banFree_req = 2000 # require this numBP NOT to be banned, e.g. by repeats

rAVA_contained_reads = {} # contained_read -> contained_under (update and select longest only)
    
MM2_write_to_file = False
MM2_readFromFile = True

INFO_rAVA_stats = {'notDovetail':0, 'eitherNameNotIn_qname_rsqRanges':0,'eitherNameNotIn_SVs_reads_callMarked':0,
                   'tooLittleBanFreeNumBP':0,'successful':0,
                   'lowOvlpNumBP':0,'traversed':0}


MM2_read_output_fo  = gzip.open(MM2_RAVA_saveFile,'rb')
for ln,line_raw in enumerate( MM2_read_output_fo ): ### USE THIS IS WE READ WRITTEN FILE!!! ## awk ' { if ((($4-$3)>=3000) && (($9-$8)>=3000)) { print } } '
    if 'gzipped input':        line_raw = line_raw.decode() #Decode line if we read it from file
    
    INFO_rAVA_stats['traversed'] += 1
    # Parse
    line = line_raw.strip('\n')    
    line = line.split()
    qname,qlen,qstart,qend,strand,rname,rlen,rstart,rend,num_matches,aln_len,mapq = line[:12]

    qlen,qstart,qend,rlen,rstart,rend = map(int,[qlen,qstart,qend,rlen,rstart,rend])
    
    if qlen < lenFilt or rlen < lenFilt: continue
    
    qcoords = [qstart,qend]
    rcoords = [rstart,rend]
    q_alnLen = qcoords[-1] - qcoords[0]
    r_alnLen = rcoords[-1] - rcoords[0]
    #/Parse

    # find distance of rAVA-ovlp to H/T of q + r
    ovlp_q_dist_h = -getRangeOvlp([0,0],qcoords)
    ovlp_q_dist_t = -getRangeOvlp([qlen,qlen],qcoords)
    ovlp_q_dists_ht = [ovlp_q_dist_h,ovlp_q_dist_t]
    
    ovlp_r_dist_h = -getRangeOvlp([0,0],rcoords)
    ovlp_r_dist_t = -getRangeOvlp([rlen,rlen],rcoords)
    ovlp_r_dists_ht = [ovlp_r_dist_h,ovlp_r_dist_t]
    #/
    
    # Check if contained
    q_isContained = False
    r_isContained = False

    # q is contained
    if max(ovlp_q_dists_ht) <= HT_margin:
        q_isContained = True
    # r is contained
    if max(ovlp_r_dists_ht) <= HT_margin:
        r_isContained = True
    #/
    
    # Save if either is contained, skip if both are contained in each other (mirror)
    if q_isContained and not r_isContained:
        init(rAVA_contained_reads,qname,[None,0]) #read, readLen
        # check if update containment info
        if rlen > rAVA_contained_reads[qname][1]:       rAVA_contained_reads[qname] = [rname,rlen]
    if r_isContained and not q_isContained:
        init(rAVA_contained_reads,rname,[None,0]) #read, readLen
        # check if update containment info
        if qlen > rAVA_contained_reads[rname][1]:       rAVA_contained_reads[rname] = [qname,qlen]
    #/
    
MM2_read_output_fo.close()
##/

### Import ref stuff
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
masking_refInfo_perChrom_EU = {}
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

### Traverse bam, reshape best alns to reads (sorted by 1:mapq,2:alnScore. distri non-ovlp)
read_length_upper = 100000*999
read_length_lower = 10000
qname_rsqRanges = {} #qname -> reshaped ranges
reads_parsed_PBfix = set()
refLens = importBAMlengths(vsDmel_bam_file)

bam_fo = pysam.AlignmentFile(vsDmel_bam_file,'rb')
bam_ln = -1 #will be set to 0 on first traversal
prev_qname = False
qname_bam_entries = []
print('[vsRef] Importing BAM & reshaping alns...')
print('############ IDE: LENGTH SPAN:',read_length_lower,'-',read_length_upper,'####################')
while True:
    # Parse BAM entry
    try:
        bam_entry = next(bam_fo)
        bam_ln += 1
    except:
        bam_entry = 'EOF'
    #/
    
    # INFOPRINTER
    if bam_ln % 100000 == 0:
        print('\t[vsRef] Processed '+str(bam_ln)+' lines!')
    #/

    if type(bam_entry) != str:
        # skip if unmapped
        if bam_entry.reference_id == -1: continue
        qname = bam_entry.qname
        
        if qname in rAVA_contained_reads: continue #skip if was contained
        if parseReadLenFromCIGAR(bam_entry.cigar) < read_length_lower: continue
        if parseReadLenFromCIGAR(bam_entry.cigar) > read_length_upper: continue
        if len(bam_entry.cigar) == 1: print(bam_entry.cigar)
        if bam_entry.has_tag('CG'):
            print('hey')
            sys.exit()
    else:
        qname = 'EOF' #dummy-set to trigger below..
         
    # Check if qname changed, then handle
    if prev_qname and not prev_qname == qname:
        rn_PBfix = getPacbioMainRead(prev_qname)
        
        # Handle read (keep track of failed pacbio reads...)
        if not rn_PBfix in reads_parsed_PBfix:
            reads_parsed_PBfix.add(rn_PBfix)
    
            # Sort qnames by 1:mapq,2:AS - both reverse. Want maximum mapq, maximum alignmentscore
            #accu_bam_entries_sorted = sorted(qname_bam_entries,key=lambda x: ( x.mapq,x.get_tag('AS') ), reverse=True )                
            accu_bam_entries_sorted = sorted(qname_bam_entries,key=lambda x: ( x.get_tag('AS'),x.mapq ), reverse=True )                
                
            # Parse & fill up alignments on read (take only new alns if it has at least X numbp not covered)
            qDistri = []
            for enum,accu_bam_entry in enumerate(accu_bam_entries_sorted):
                # Parse
                entry = parseBAMentry(accu_bam_entry)
                
                # Calculate which part of aln we should save reference-stuff for (full if not ovlp any existing als..)
                qRange_to_scan = entry['qcoords'][:] #assume we want to keep full aln, recalc with loop below. COPY qcoords to we can tweak it
                # loop will continually chop the qRange_to_scan, if it ovlps with any existing qRange.
                for qDistri_entry in qDistri:
                    exist_qRange = qDistri_entry['qcoords']
                    
                    ## check if ovlp
                    ovlp = getRangeOvlp(exist_qRange,qRange_to_scan)
                    
                    # Check if full ovlp, skip
                    if ovlp >= (qRange_to_scan[-1]-qRange_to_scan[0]):
                        qRange_to_scan = None
                        break #break on first. 
                    
                    # Else calc which region of qRange_to_scan should be chopped
                    elif ovlp > 0:
                        # get range of ovlp
                        ovlpRange = getOvlpRange(exist_qRange,qRange_to_scan)
                        
                        # reshape rcoords outside ovlp:
                        # if qRange_to_scan is left of exist_qRange
                        if qRange_to_scan[0] <= exist_qRange[0] and qRange_to_scan[-1] <= exist_qRange[-1]:
                            qRange_to_scan[-1] = ovlpRange[0]-1 #update end qcoord with start-range of ovlpRange, 1BP margin
                        # if qRange_to_scan is right of exist_qRange
                        if qRange_to_scan[0] >= exist_qRange[0] and qRange_to_scan[-1] >= exist_qRange[-1]:
                            qRange_to_scan[0] = ovlpRange[-1]+1 #update start qcoord with start-range of ovlpRange, 1BP margin
    
                    #/
                    
                    # Check if chopped region became too small to care about, then skip it
                    if (qRange_to_scan[-1]-qRange_to_scan[0]) <= 100:
                        qRange_to_scan = None
                        break
                    #/
                    ##/
                
                # Parse cigar & save (needed for reshape)
                if qRange_to_scan != None:
                    entry['cigar'] = HOLDER(accu_bam_entry.cigar)
                    entry['accu_bam_enum'] = enum
                    
                    # Check if it was reshaped, then flag reshape and do reshape
                    entry['wasReshaped'] = False
                    if entry['qcoords'] != qRange_to_scan:
                        ## do reshape
                        # Try reshape, keep track if some rsqcoord was not reshaped
                        rsrcoords = []
                        qstart_qend_notReshaped = [True,True] # assume they arent, they will prove themselves
                            
                        coords_reshaped = getQueryPoses_fromTargetPoses_byCigar(qRange_to_scan,aln_entry=entry,reshapeToTarget=True)
                        for rsrcoord,rsqcoord in coords_reshaped:
                            rsrcoords.append(rsrcoord)
                            
                            if rsqcoord == qRange_to_scan[0]:
                                qstart_qend_notReshaped[0] = False
                            elif rsqcoord == qRange_to_scan[-1]:
                                qstart_qend_notReshaped[-1] = False
                        #/
                        
                        # Check if missing coords, then assume its close to border range (cigar can be hard to traverse with precision...)
                        for coord_idx,RSfailed in enumerate(qstart_qend_notReshaped):    
                            if RSfailed:
                                if strand_isfv(entry['strand']):
                                    rsrcoords.append(entry['rcoords'][coord_idx])
                                elif strand_isrv(entry['strand']):
                                    # Since reverse, need to switch which rcoord to add...
                                    coord_idx_rvFix = None
                                    if coord_idx == 0:  coord_idx_rvFix = 1
                                    if coord_idx == 1:  coord_idx_rvFix = 0
                                    #/
                                    rsrcoords.append(entry['rcoords'][coord_idx_rvFix])                                    
                        #/
                        ##/
                        
                        # reshape-specifc: save info
                        entry['wasReshaped'] = True
                        entry['orig_qcoords'] = entry['qcoords'][:] #hard-copy
                        entry['qcoords'] = qRange_to_scan[:] #update qcoords at entry
                        entry['orig_rcoords'] = entry['rcoords'][:] #hard-copy
                        entry['rcoords'] = [min(rsrcoords),max(rsrcoords)]
                    
                    entry['bam_ln'] = bam_ln
                    entry['len'] = entry['rcoords'][-1]-entry['rcoords'][0]
                    qDistri.append(entry)
            #/
            
            # Sort by qcoords, assign order
            qDistri_sorted = sorted(qDistri,key=lambda x: x['qcoords'][0])
            for order,entry in enumerate(qDistri_sorted):
                entry['order'] = order
            #/
            
            # Save qDistri
            qname_rsqRanges[prev_qname] = qDistri_sorted
        
        # Reset qname_bam_entries
        qname_bam_entries = []
    #/
    
    # Check if we was EOF, then break
    if type(bam_entry) == str and bam_entry == 'EOF':      break
    
    # Accumulate bam entries of read
    prev_qname = qname
    qname_bam_entries.append(bam_entry)

bam_fo.close()
print('\tDone!')
print('Number of saved reads: '+str(len(qname_rsqRanges)))
###/

### Make IDX-map of qname_rsqRanges entries (ReShaped Query Entries (rsqE))
qname_rsq_idx_map = {}
for qname,qDistris in qname_rsqRanges.items():
    for qDistri in qDistris:
        assignIdx(qDistri)
        qname_rsq_idx_map[qDistri['idx']] = qDistri
qname_rsq_idx_map_ROMD = rangeOverlaps_makeMappingDict(qname_rsq_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
###/
### Parse read lens (also for coverage est)
qnamesLens = {}
for qname,qDistri in qname_rsqRanges.items():
    for anyE in qDistri:
        qnamesLens[qname] = anyE['read_len']
        break
IDE_covTot = 0
for read,readlen in qnamesLens.items(): IDE_covTot += readlen
print('Estimated coverage (given sum of all refLens) = '+str(round(IDE_covTot/sum(refLens.values()),4)))
###/

### Find what exists at qRanges (update/saves entry in-place)
## Multiprocess fill
def fill_qname_rsqEs(rsqEs,fillGFROs=True):
    rsqE_saveData = {}
    for rsqE in rsqEs:
        rcoords = rsqE['rcoords'] #may or may not be reshaped...!!
        rname = rsqE['rname']
        
        ## TWO SCENARIOS;
        # 1: reshape raw ranges for NONMASKED features (i.e. GTF)
        # 2: reshape compiled+selected for masking
        
        # init
        SAVEDATA = {}
        SAVEDATA['gfros'] = HOLDER([])
        SAVEDATA['maskings'] = HOLDER([])
        SAVEDATA['masked_covFrac'] = 0
        
        # 1.
        if 1:
            for RUNVA,feat_id,ROMD,ROMD_step in (
                                            (fillGFROs,'gfros',GFRO_idx_map_ROMD,1000),
                                            ):
                
                if not RUNVA: continue
            
                ovlps = rangeOverlaps_lookup(rcoords,ROMD,ROMD_step,map_lookup_key=rname)
    
                # traverse & save
                for oS,idx,oE in ovlps:
                    if (oE-oS) <= 0: continue #skip low ovlps
                    
                    ## do reshape
                    # Try reshape, keep track if some rsrcoord was not reshaped
                    rsqcoords = []
                    rstart_rend_notReshaped = [True,True] # assume they arent, they will prove themselves
                        
                    coords_reshaped = getQueryPoses_fromTargetPoses_byCigar([oS,oE],aln_entry=rsqE,reshapeToQuery=True)
                    for rsqcoord,rsrcoord in coords_reshaped:
                        rsqcoords.append(rsqcoord)
                        
                        if rsrcoord == oS:
                            rstart_rend_notReshaped[0] = False
                        elif rsrcoord == oE:
                            rstart_rend_notReshaped[-1] = False
                    #/
                    
                    ## IDE CHECK
                    if all(rstart_rend_notReshaped):
                        print('nonMask',oS,oE,rsqE['rcoords'])
                    ##/
                    
                    # Check if missing coords, then assume its close to border range (cigar can be hard to traverse with precision...)
                    for coord_idx,RSfailed in enumerate(qstart_qend_notReshaped):    
                        if RSfailed:
                            if strand_isfv(rsqE['strand']):
                                rsqcoords.append(rsqE['qcoords'][coord_idx])
                            elif strand_isrv(rsqE['strand']):
                                # Since reverse, need to switch which rcoord to add...
                                coord_idx_rvFix = None
                                if coord_idx == 0:  coord_idx_rvFix = 1
                                if coord_idx == 1:  coord_idx_rvFix = 0
                                #/
                                rsqcoords.append(rsqE['qcoords'][coord_idx_rvFix])                                    
                    #/
                    ##/
                    tmp_save = {'qcoords':[min(rsqcoords),idx,max(rsqcoords)],'rcoords':[oS,oE],'gfro_idx':idx}
                    SAVEDATA[feat_id]().append(tmp_save)
            #/
        #/1.
        
        # 2.
        ovlps = rangeOverlaps_lookup(rcoords,masking_idx_map_ROMD,100,map_lookup_key=rname)
        
        ## OvlpRanges: [IDE: dont know if smart... but lets go with computeOvlpRanges....]
        # make rRanges, save base
        rRanges = []
        rRanges.append([rcoords[0],'base',rcoords[-1]])
        for oS,idx,oE in ovlps:
            if (oE-oS) <= 0: continue #skip low numBP ovlps/single-BP ovlps
            feat_entry = None
            if idx in masking_idx_map:            feat_entry = masking_idx_map[idx]
            feat_rstart,feat_rend = feat_entry['rcoords']                
            rRanges.append( [feat_rstart,idx,feat_rend] )
        #/
        
        # Make ovlpsRanges + traverse
        ovlpRanges = computeOvlpRanges_wIdxs(rRanges)
        
        # First traverse and calc numBP for each masking assignemnt.
        # Then, at conflicts, select only the one with most numBP asn.
        maskClasses_numBPs = {}
        masked_numBP_pan = 0
        for rangee in ovlpRanges:
            # skip if not rsqE present
            if not 'base' in rangee[1]: continue
            maskClasses_added = set() #keep track of which maskclasses we addedrs
            for idx in rangee[1]:
                if idx == 'base': continue
                maskClass = masking_idx_map[idx]['class']
                # Check if we already added mask-class
                if maskClass in maskClasses_added: continue
            
                # Add masked numBP
                numBP = rangee[-1]-rangee[0]
                masked_numBP_pan += numBP
                
                # Dont add maskClass if its a low-resolution one (i.e. dont prio it at conflicts)
                if maskClass in ('Simple_repeat','Low_complexity',):
                    numBP = 0
                    
                maskClasses_added.add(maskClass)
                
                # save numBP
                init(maskClasses_numBPs,maskClass,0)
                maskClasses_numBPs[maskClass] += numBP
                
        # Calc overall masked covFrac for RSQ
        covFrac_masked = 0
        if maskClasses_numBPs:
            covFrac_masked = masked_numBP_pan/(rcoords[-1]-rcoords[0])
        SAVEDATA['masked_covFrac'] = covFrac_masked
        #/
    
        # Sort once, for gain in lookup
        maskClasses_numBPs_sorted = sorted(maskClasses_numBPs.items(),key=lambda x: x[1], reverse=True)
        #/
        
        # Traverse again           
        for rangee in ovlpRanges:
            # skip if not rsqE present
            if not 'base' in rangee[1]: continue
            
            # skip if only rsqE present
            if len(rangee[1]) == 1 and 'base' in rangee[1]: continue
        
            # Find which maskClasses exist at rangee
            rangee_maskClasses = {}
            for idx in rangee[1]:
                if idx == 'base': continue
                maskClass = masking_idx_map[idx]['class']
                init(rangee_maskClasses,maskClass,set())
                rangee_maskClasses[maskClass].add(idx)
            #/
            
            # Find numBP of top maskCLass at rangee
            rangee_top_hit_numBP = None
            for maskClass,numBP in maskClasses_numBPs_sorted:
                if maskClass in rangee_maskClasses:
                    rangee_top_hit_numBP = numBP
                    break
            #/
            
            # Add only top masking hits of rangee
            rangee_maskClasses_filt = []
            rangee_maskClasses_filt_wIdx = {}
            for maskClass in rangee_maskClasses:
                if maskClasses_numBPs[maskClass] >= rangee_top_hit_numBP*0.95:
                    rangee_maskClasses_filt.append(maskClass)                  
                    init(rangee_maskClasses_filt_wIdx,maskClass,set())
                    rangee_maskClasses_filt_wIdx[maskClass].update(rangee_maskClasses[maskClass])
            #/
            
            ## do reshape
            # Try reshape, keep track if some rsrcoord was not reshaped
            rsqcoords = []
            rstart_rend_notReshaped = [True,True] # assume they arent, they will prove themselves
                
            coords_reshaped = getQueryPoses_fromTargetPoses_byCigar([rangee[0],rangee[-1]],aln_entry=rsqE,reshapeToQuery=True)
            for rsqcoord,rsrcoord in coords_reshaped:
                rsqcoords.append(rsqcoord)
                
                if rsrcoord == rangee[0]:
                    rstart_rend_notReshaped[0] = False
                elif rsrcoord == rangee[-1]:
                    rstart_rend_notReshaped[-1] = False
            #/

            # Check if missing coords, then assume its close to border range (cigar can be hard to traverse with precision...)
            for coord_idx,RSfailed in enumerate(qstart_qend_notReshaped):    
                if RSfailed:                    
                    if strand_isfv(rsqE['strand']):
                        rsqcoords.append(rsqE['qcoords'][coord_idx])
                    elif strand_isrv(rsqE['strand']):
                        # Since reverse, need to switch which rcoord to add...
                        coord_idx_rvFix = None
                        if coord_idx == 0:  coord_idx_rvFix = 1
                        if coord_idx == 1:  coord_idx_rvFix = 0
                        #/
                        rsqcoords.append(rsqE['qcoords'][coord_idx_rvFix])                                    
            #/
            ##/
            
            # Save
            tmp_save = {'qcoords':[min(rsqcoords),idx,max(rsqcoords)],'rcoords':[rangee[0],rangee[-1]],'maskClassis':rangee_maskClasses_filt,'maskClassis_wIdxs':rangee_maskClasses_filt_wIdx,
                                   'len':max(rsqcoords)-min(rsqcoords)}
            SAVEDATA['maskings']().append(tmp_save)
            #/
        #/
        ##/OvlpRanges
        #/2
        
        # Save to rsqE
        SAVEDATA['maskClassis_numBPs'] = maskClasses_numBPs_sorted
        rsqE_saveData[rsqE['idx']] = SAVEDATA
        #/
        
    return rsqE_saveData
    
if 1 and 'Fill rsqEs?':
    # setup multiprocessing
    num_processes = 10
    batch_size = math.ceil(len(qname_rsq_idx_map)/num_processes)
    print('Setting up multiprocessing jobs! Will use '+str(num_processes)+' threads!')
    jobs_setup = []
    for rsqE in qname_rsq_idx_map.values():
        # check if no job batch started
        if not jobs_setup:                          jobs_setup.append([])
        # Check if job batch is filled
        if len(jobs_setup[-1]) >= batch_size:       jobs_setup.append([])
        
        # Add rsqE to job batch
        jobs_setup[-1].append(rsqE)
    #/
    # start jobs
    pool = Pool(processes=num_processes)
    jobs = {}
    print('\tPreparing and starting jobs...')
    for enum,job_setup_data in enumerate(jobs_setup):
        jobs[enum] = pool.apply_async(fill_qname_rsqEs,args=(job_setup_data,True))
        #break
    pool.close()
    pool.join()
    print('\tRunning '+str(len(jobs))+' jobs!')
    #/
    # parse jobs
    print('\tParsing...')
    SC_oSC_paths_singleExtended = {}
    for enum,job in jobs.items():
        job_result = job.get()
        # Save all data from multiprocessing at rsqE
        for rsqE_idx,data in job_result.items():
            rsqE = qname_rsq_idx_map[rsqE_idx]
            for key,content in data.items():
                rsqE[key] = content
    print('\tDone!')
    #/
##/Multiprocess-fill
###/

### Traverse RSQs, calculate "chromothripsis"-index. NTS: need to distinguish between DUPs and rearrs.
def dumpReadsFromBam(bam_source,bam_nf,qnames=None,repma_frac_CO=0,readLen_CO=0):
    bam_fo = pysam.AlignmentFile(bam_source,'rb')
    bam_nf = pysam.AlignmentFile(bam_nf,'wb',template=bam_fo)
    for entry in bam_fo:
        qname = entry.qname
        if qnames == None or qname in qnames:
            readlen = parseReadLenFromCIGAR(entry.cigar)
            if readlen < readLen_CO: continue
            entry.set_tag('RL',readlen)
            
            # quickly scan repMa
            if repma_frac_CO:
                try:
                    rname = entry.reference_name
                    rcoords = [entry.reference_start,entry.reference_end]
                except:
                    print('skipping',entry.reference_name,entry.reference_start,entry.reference_end)
                    continue
                ovlps = rangeOverlaps_lookup(rcoords,masking_idx_map_ROMD,100,map_lookup_key=rname)
                rangees = [ [rcoords[0],'base',rcoords[-1]] ]
                for oS,idx,oE in ovlps:
                    rangees.append([oS,'masked',oE])
                ovlps2 = computeOvlpRanges_wIdxs(rangees)
                base_numBP = 0
                mask_numBP = 0
                for rangee in ovlps2:
                    if not 'base' in rangee[1]: continue
                    base_numBP += rangee[-1]-rangee[0]
                    if 'masked' in rangee[1]:
                        mask_numBP += rangee[-1]-rangee[0]
                maskCovFrac = mask_numBP / base_numBP
                if maskCovFrac >= repma_frac_CO:
                    continue #skip if repeatmasked too much
            #/
                
            if type(qnames) == dict and qnames[qname]:
                if type(qnames[qname]) == int:
                    entry.set_tag('RN',qnames[qname])
                else:
                    entry.set_tag('RN','||'.join(list(map(str,sorted(qnames[qname])))))
            bam_nf.write(entry)
    bam_fo.close()
    bam_nf.close()

def getRangeBins(rcoords,bin_step = 1000):
    rstart,rend = rcoords
    bins_to_scan = []
    #X BP bin interval. +1 to not miss last one (python rangee open end). -1 on start and +1 on end to not miss coordX = 999 and coordY = 1001
    for i in range(math.floor(rstart/bin_step)-1,math.ceil(rend/bin_step)+1+1):
        if not i in bins_to_scan:
            bins_to_scan.append(i)
    return bins_to_scan

reads_wDup = 0
reads_woDup = 0


rsqE_len_cutoff = 1000 # skip short RSQs
rsqE_masked_covFrac_skip = 0.5 #assume repeat if above this covFrac masked
DUP_len_min = 2000 #skip evidence of DUP if DUP size is lower than this
INSanchor_minSpan = 500 #numBP to shorten refRange anchors to for INScalls
rsqE_adj_clusts_qHT_len_thresh = 3000 # minimum size of rsqE_adjClust which exist at read start/end

qname_rsqE_adj_clusts = {} #ReShaped Query Entry adjacent(alignments) clustered
qname_rsqE_dupClusts = {} 
qname_rsqE_DUPs = {}
INSanchors_pan_idx_map = {}
INScalls_anchors = {} # rname -> binn_id [binned by 1000bp] -> {data_on_idx_map_format}
INScalls_anchors_idx_map = {}
INScall_idx_map = {}
DELcalls_ROUGH = {} # "quick" IDE caller. Rname -> [arr_of_calls] <-- Later cluster them into "same event"
DELcalls_pan_idx_map = {} # call ALL DELs, like with INSes.
panCallINSes = True
DUPclust_idx_map = {}
AC_idx_map = {} #Adjacent Clustered (alignments)
for ln,(qname,qDistri) in enumerate(qname_rsqRanges.items()):
    if ln%100 == 0: print('[WIP1] Processed '+str(ln)+'/'+str(len(qname_rsqRanges)))
    
    # 1. Sort rsqEs
    covFrac_masked_self = 0
    covFrac_masked_nonMaskedRsqEs = 0
    rsqEs_nonMasked = []
    prev_rsqE = None
    for enum,rsqE in enumerate(sorted(qDistri,key=lambda x: x['qcoords'][0])):
        covFrac_masked_self += rsqE['masked_covFrac']*(rsqE['len']/rsqE['read_len'])

        prev_rsqE = rsqE #Update prev
    
        if rsqE['len'] < rsqE_len_cutoff: continue
        if rsqE['masked_covFrac'] < rsqE_masked_covFrac_skip:
            rsqEs_nonMasked.append(rsqE)
            covFrac_masked_nonMaskedRsqEs += rsqE['masked_covFrac']*(rsqE['len']/rsqE['read_len'])
            
    if covFrac_masked_self >= 0.6 or covFrac_masked_nonMaskedRsqEs >= 0.3: continue #skip if large parts of read is just repeat
    #/1
    
    ## 2. Cluster adjacent alns, find if read has indications of large events
    clust_dist = 5000 #merge if within this dist on r.
    clust_dist_max = 25000 #maximum dist to merge, regardless of numBP masked on ref
    
    # Init
    rname_rsqE_adj_clusts = {} # rname -> clustered rsqEs
    for rsqE in rsqEs_nonMasked:
        rname = rsqE['rname']
        rstart,rend = rsqE['rcoords'] #extract out of arr for copy
        
        init(rname_rsqE_adj_clusts,rname,[])
        tmp_save = {'rname':rname,'rcoords':[rstart,rend],'members_raw':[rsqE['idx']],'qname':qname,'len':rend-rstart}
        rname_rsqE_adj_clusts[rname].append(tmp_save)
    
    # Clust
    for rname,entries in rname_rsqE_adj_clusts.items():
        while True:
            resetIter = False
            for enum1,entry1 in enumerate(entries):
                for enum2,entry2 in enumerate(entries):
                    if enum1 == enum2: continue
                    rsqEs_dist = -getRangeOvlp(entry1['rcoords'],entry2['rcoords'])
                    # Calc numBP which was masked between entries
                    numBP_masked_deduct = 0
                    if rsqEs_dist <= clust_dist_max and rsqEs_dist > 0: #check if dist is lower than maximum, and that rsqEs do not ovlp
                        if mean(entry1['rcoords']) < mean(entry2['rcoords']): #entry1 is before entry2
                            scanRange = [entry1['rcoords'][-1],entry2['rcoords'][0]]
                        else: #entry2 is before entry1
                            scanRange = [entry2['rcoords'][-1],entry1['rcoords'][0]]
                        
                        deduct_ranges = []
                        deduct_ranges.append([scanRange[0],'scanRange',scanRange[-1]])
                        for oS,idx,oE in rangeOverlaps_lookup(scanRange,masking_idx_map_ROMD,100,map_lookup_key=rname):
                            deduct_ranges.append([oS,idx,oE])
                        ovlps = computeOvlpRanges_wIdxs(deduct_ranges)
                        for rangee in ovlps:
                            if not 'scanRange' in rangee[1]: continue #skip if not scanRange exist
                            if len(rangee[1]) == 1: continue #skip if only scanRange exist
                            numBP_masked_deduct += rangee[-1]-rangee[0] #add numBP to deduct
                        
                    #/
                    if (rsqEs_dist-numBP_masked_deduct) <= clust_dist:
                        entry1['members_raw'] += [[rsqEs_dist,'dist']] + entry2['members_raw']
                        entry1['rcoords'][0] = min([entry1['rcoords'][0],entry2['rcoords'][0]])
                        entry1['rcoords'][-1] = max([entry1['rcoords'][-1],entry2['rcoords'][-1]])
                        del entries[enum2]
                        resetIter = True
                        break
                        
                if resetIter: break
            if not resetIter: break
    
    # Extract clusts
    rsqE_adj_clusts = []
    for rname in rname_rsqE_adj_clusts:
        for clust in rname_rsqE_adj_clusts[rname]:
            rsqE_adj_clusts.append(clust)
            
    # Add rsqE idxs as members + parse qSpan
    for enum,clust in enumerate(rsqE_adj_clusts):
        clust_rsqE_idxs = []
        for member in clust['members_raw']:
            if type(member) == int and member in qname_rsq_idx_map:
                clust_rsqE_idxs.append(member)
        clust['members'] = clust_rsqE_idxs
        
        clust_qcoords = []
        for member in clust['members']:
            if type(member) == int and member in qname_rsq_idx_map:
                rsqE = qname_rsq_idx_map[member]
                clust_qcoords += rsqE['qcoords']
        clust['qcoords'] = [min(clust_qcoords),max(clust_qcoords)]
        
    # Remove HT clusts which are too short
    for scanKey in (0,-1):
        if rsqE_adj_clusts: #check we still have clusts
            clust = rsqE_adj_clusts[scanKey]
            if (clust['qcoords'][-1]-clust['qcoords'][0]) < rsqE_adj_clusts_qHT_len_thresh:     del rsqE_adj_clusts[scanKey]
    
    
    ## Go through clusts, find which are contained in another on q (e.g. INS. If multiple clusts contain each other, means rearr)
    # Find clust ovlps on q
    for enum1,clust1 in enumerate(rsqE_adj_clusts):
        for enum2,clust2 in enumerate(rsqE_adj_clusts):
            if enum1 == enum2: continue
            numBP_ovlp = getRangeOvlp(clust1['qcoords'],clust2['qcoords'])
         
            clust1_fracContained = numBP_ovlp / (clust1['qcoords'][-1]-clust1['qcoords'][0])
            clust2_fracContained = numBP_ovlp / (clust2['qcoords'][-1]-clust2['qcoords'][0])            
            
            # Only consider clust1, since clust2 will become clust1 later ...
            if clust1_fracContained >= 1:
                SK = 'fullyContained'
                init(clust1,SK,[])
                clust1[SK].append([clust1_fracContained,enum2])
            # If not 100% contained, means partial ovlp and potential REARR
            elif clust1_fracContained > 0 and clust2_fracContained > 0:
                SK = 'partialContainment'
                init(clust1,SK,[])
                clust1[SK].append([clust1_fracContained,enum2])
    #/
    # Find ref BP + LERI
    for enum,clust in enumerate(rsqE_adj_clusts):
        # grab rsqEs, sort by q
        clust_rsqEs = []
        for member in clust['members']:
            if type(member) == int and member in qname_rsq_idx_map:
                clust_rsqEs.append(qname_rsq_idx_map[member])
        clust_rsqEs = sorted(clust_rsqEs,key=lambda x: x['qcoords'][0])
        clust_rsqE_first = clust_rsqEs[0]
        clust_rsqE_last = clust_rsqEs[-1]
        #/
        
        # Calc ref BP + LERI for clust BPs. BPL = BP-left, BPR = BP-right (on q)
        clust_qBPL_refBPLERI = [None,None,'notReadEnd'] # refBP,refLERI,isReadEnd
        clust_qBPR_refBPLERI = [None,None,'notReadEnd'] # refBP,refLERI,isReadEnd
        if strand_isfv(clust_rsqE_first):
            clust_qBPL_refBPLERI[0] = clust_rsqE_first['rcoords'][0]
            clust_qBPL_refBPLERI[1] = 'left'
        elif strand_isrv(clust_rsqE_first):
            clust_qBPL_refBPLERI[0] = clust_rsqE_first['rcoords'][-1]
            clust_qBPL_refBPLERI[1] = 'right'
        if strand_isfv(clust_rsqE_last):
            clust_qBPR_refBPLERI[0] = clust_rsqE_last['rcoords'][-1]
            clust_qBPR_refBPLERI[1] = 'right'
        elif strand_isrv(clust_rsqE_last):
            clust_qBPR_refBPLERI[0] = clust_rsqE_last['rcoords'][0]
            clust_qBPR_refBPLERI[1] = 'left'
        if enum == 0:                            clust_qBPL_refBPLERI[2] == 'isReadEnd'
        if enum == len(rsqE_adj_clusts)-1:       clust_qBPR_refBPLERI[2] == 'isReadEnd'
        #/
        # Save
        clust['qBPL_refBPLERI'] = clust_qBPL_refBPLERI
        clust['qBPR_refBPLERI'] = clust_qBPR_refBPLERI
        #/
    #/
    # Assign enum + save mark at rsqE
    rsqE_adj_clusts_enumed = {}
    for enum,clust in enumerate(rsqE_adj_clusts):
        clust['enum'] = enum
        for member in clust['members']:
            if type(member) == int and member in qname_rsq_idx_map:
                rsqE = qname_rsq_idx_map[member]
                if 'IDE-check' and 'memberOfClust' in rsqE and rsqE['memberOfClust'] != enum: print('uahisuuuuuuuuuuuuuuuuuduahsiudhasd')
                rsqE['memberOfClust'] = clust['enum']
        rsqE_adj_clusts_enumed[enum] = clust
    #/
    ##/
    
    # save
    if rsqE_adj_clusts_enumed:
        qname_rsqE_adj_clusts[qname] = rsqE_adj_clusts_enumed
        # Update 27march Save at idx-map
        for AC in rsqE_adj_clusts_enumed.values():
            assignIdx(AC)
            AC_idx_map[AC['idx']] = AC
        #/
    #/
    ##/2
    
    ## 3. Check if read has evidence of DUP (find ovlping rsqE's on ref)
    rsqE_rclusts = []
    for clust in rsqE_adj_clusts:
        # Parse rsqEs of adjClust
        rRanges = []
        for rsqE_idx in clust['members']:
            rsqE = qname_rsq_idx_map[rsqE_idx]
            rstart,rend = rsqE['rcoords']
            rRanges.append([rstart,rsqE_idx,rend])
    
        # Cluster rsqEs into DUPranges
        ovlps = computeOvlpRanges_wIdxs(rRanges)

        for rangee in ovlps:
            if (rangee[-1]-rangee[0]) < 100: continue #skip potential double-aln ovlp peaks
            
            # Only select rangees where at least x2 exist (dont extend DUPclust in anchors!)
            if len(rangee[1]) < 2: continue
        
            # Check if adjacent to existing clust
            initClust = True #assume we will add, prove us wrong
            for existClust in rsqE_rclusts:
                if existClust[1] == rname and getRangeOvlp(existClust,rangee) >= 0:
                    existClust[-1] = rangee[-1] #update clust end coord
                    existClust[2].update(rangee[1]) #add rsqE members
                    initClust = False
            #/
            if initClust:
                rname = clust['rname']
                rsqE_rclusts.append([rangee[0],rname,set(rangee[1]),clust['idx'],rangee[-1]])
    #/
    
    # Filter small clusts
    rsqE_rclusts_filt1 = []
    for clust in rsqE_rclusts:
        clust_size = clust[-1]-clust[0]
        if clust_size < DUP_len_min: continue
        rsqE_rclusts_filt1.append(clust)
    #/
    
    # Estimate DUP CN at clusts
    DUP_clusts_parsed = {}
    for enum,clust in enumerate(rsqE_rclusts_filt1):
        DUP_clusts_parsed[enum] = {'clust_raw':clust,'CN_est':len(clust[2])}
    #/

    # Chain adjacent clusts using shared alns or ClustContainmentOnQ ( find panspan DUPs with intra-DEL/DUPs/INS-splitters [ example: bantam ] )
    rsqE_rclusts_chained = [ entry[:] for entry in rsqE_rclusts_filt1 ]
    while True:
        resetIter = False
        for enum1,clust1 in enumerate(rsqE_rclusts_chained):
            for enum2,clust2 in enumerate(rsqE_rclusts_chained):
                clust1_rcoords = [clust1[0],clust1[-1]]
                clust2_rcoords = [clust2[0],clust2[-1]]
                if enum1 == enum2: continue
                clusts_dist = -getRangeOvlp(clust1_rcoords,clust2_rcoords)
                # Calc numBP which was masked between entries
                numBP_masked_deduct = 0
                if clusts_dist <= clust_dist_max and clusts_dist > 0: #check if dist is lower than maximum, and that rsqEs do not ovlp
                    if mean(clust1_rcoords) < mean(clust2_rcoords): #entry1 is before entry2
                        scanRange = [clust1_rcoords[-1],clust2_rcoords[0]]
                    else: #entry2 is before entry1
                        scanRange = [clust2_rcoords[-1],clust1_rcoords[0]]
                    
                    deduct_ranges = []
                    deduct_ranges.append([scanRange[0],'scanRange',scanRange[-1]])
                    for oS,idx,oE in rangeOverlaps_lookup(scanRange,masking_idx_map_ROMD,100,map_lookup_key=rname):
                        deduct_ranges.append([oS,idx,oE])
                    ovlps = computeOvlpRanges_wIdxs(deduct_ranges)
                    for rangee in ovlps:
                        if not 'scanRange' in rangee[1]: continue #skip if not scanRange exist
                        if len(rangee[1]) == 1: continue #skip if only scanRange exist
                        numBP_masked_deduct += rangee[-1]-rangee[0] #add numBP to deduct
                    
                #/
                if (clusts_dist-numBP_masked_deduct) <= clust_dist:
                    clusts_shared_rsqEs = clust1[2].intersection(clust2[2])
                    
                    # check if contained within each other on q
                    clusts_qcoordsspan = {}
                    for tmp_enum,tmp_clust in enumerate( [clust1,clust2] ):
                        for rsqE_idx in tmp_clust[2]:
                            rsqE = qname_rsq_idx_map[rsqE_idx]
                            init(clusts_qcoordsspan,tmp_enum,[])
                            clusts_qcoordsspan[tmp_enum] += rsqE['qcoords']
                    clusts_qspans = {}
                    for enum,coords in clusts_qcoordsspan.items():
                        clusts_qspans[enum] = [min(coords),max(coords)]
                    clusts_qspans_ovlp = getRangeOvlp(clusts_qspans[0],clusts_qspans[1])
                    #/
                    
                    if clusts_shared_rsqEs or clusts_qspans_ovlp >= 100:
                        clust1[2].update(clust2[2])
                        clust1[0] = min([clust1[0],clust2[0]])
                        clust1[-1] = max([clust1[-1],clust2[-1]])
                        del rsqE_rclusts_chained[enum2]
                        resetIter = True
                        break
            if resetIter: break
        if not resetIter: break   
    #/
    
    # Compile DUPs from chained clusts, assign CN ests of DUP parts
    DUPs = {}
    for DUPenum,chainedClusts in enumerate(rsqE_rclusts_chained):
        # Init
        tmp_save = {'rname':chainedClusts[1],'rcoords':[chainedClusts[0],chainedClusts[-1]],'len':chainedClusts[-1]-chainedClusts[0],
                    'raw':chainedClusts,'AC_idx':chainedClusts[3],'qname':qname,'members':chainedClusts[2]}
        
        
        # Fill with subClusts, calc numBP per CNest
        for enum,dupClust in DUP_clusts_parsed.items():
            if dupClust['clust_raw'][2].intersection(chainedClusts[2]) == tmp_save['members']:
                """
                SK = 'clusts'
                init(tmp_save,SK,[])
                tmp_save[SK].append(dupClust)
                """
                SK = 'CNest'
                SK2 = dupClust['CN_est']
                init(tmp_save,SK,{})
                init(tmp_save[SK],SK2,0)
                tmp_save[SK][SK2] += dupClust['clust_raw'][-1]-dupClust['clust_raw'][0]
        
        
        # Distribute all rsqEs on Q, sort, check left-/right-most if they extend outside DUP. If so, mark as anchor + reshape where DUP exist on q.
        DUP_rsqEs = []
        rsqE_idxs_added = set() #keep track of which we added already
        #for dupClust in tmp_save['clusts']:
        for rsqE_idx in tmp_save['members']:
            if rsqE_idx in rsqE_idxs_added: continue
            rsqE_idxs_added.add(rsqE_idx)
            DUP_rsqEs.append(qname_rsq_idx_map[rsqE_idx])
        DUP_rsqEs = sorted(DUP_rsqEs,key=lambda x: x['qcoords'][0])
        
        rsqE_qLeftRight = { 'left':{'rsqE':DUP_rsqEs[0]},'right':{'rsqE':DUP_rsqEs[-1]} }
        tmp_DUP_qcoords = []
        for LERI,data in rsqE_qLeftRight.items():
            LERIrsqE_rcoords = data['rsqE']['rcoords']
            
            # Find BP on r where LERIrsqE and DUP stop/start ovlping each other
            rRanges = []
            rRanges.append([LERIrsqE_rcoords[0],'LERIrsqE',LERIrsqE_rcoords[-1]])
            rRanges.append([tmp_save['rcoords'][0],'DUP',tmp_save['rcoords'][-1]])
            
            ovlps = computeOvlpRanges_wIdxs(rRanges)
            prevRangee = None
            LERIrsqE_rBP = None
            DUP_anchor_rRange = []
            for rangee in ovlps:
                # check if prevRangee was set and if we changed
                if prevRangee != None and len(rangee[1]) != len(prevRangee[1]):
                    LERIrsqE_rBP = rangee[0]
                    break
                
                prevRangee = rangee
                
                # Check if add to DUP anchor rRange
                if len(rangee[1]) == 1 and 'LERIrsqE' in rangee[1]:
                    DUP_anchor_rRange = [rangee[0],rangee[-1]]
            #/
            
            # Check if rBP was found, else means LERIrsqE is 100% contained.
            DUP_qBP = None
            if LERIrsqE_rBP != None:
                coords_reshaped = getQueryPoses_fromTargetPoses_byCigar([LERIrsqE_rBP],aln_entry=data['rsqE'],reshapeToQuery=True)
                for rsqcoord,rsrcoord in coords_reshaped:
                    DUP_qBP = rsqcoord
            else:
                if LERI == 'left':              DUP_qBP = data['rsqE']['qcoords'][0]
                elif LERI == 'right':           DUP_qBP = data['rsqE']['qcoords'][-1]
            #/
            
            # Check LERI on ref for the anchor_rRange
            rLERI = None
            if DUP_anchor_rRange:
                if   mean(DUP_anchor_rRange) < mean(tmp_save['rcoords']):  rLERI = 'left'
                elif mean(DUP_anchor_rRange) > mean(tmp_save['rcoords']):  rLERI = 'right'    
            #/

            # Save
            tmp_DUP_qcoords.append(DUP_qBP)
            SK = LERI
            tmp_save['r_anchor_q'+LERI] = {'rRange':DUP_anchor_rRange,'rLERI':rLERI,'rsqE_idx':data['rsqE']['idx']}
        
        tmp_save['qcoords'] = [min(tmp_DUP_qcoords),max(tmp_DUP_qcoords)] #save qcoords which we parsed from anchor-alns
        #/
        
        ## Traverse DUP rsqEs, track when we "recoil" on ref
        # Call recoils
        DUP_recoils = []
        for rsqEenum,_ in enumerate(DUP_rsqEs):
            if rsqEenum == 0: continue #backwards checking
            prev_rsqE = DUP_rsqEs[rsqEenum-1]
            cur_rsqE = DUP_rsqEs[rsqEenum]
            
            # check if recoiled: 1/ they ovlp on r, or 2/ prev is before cur
            comp_is_recoil = False
            if getRangeOvlp(prev_rsqE['rcoords'],cur_rsqE['rcoords']) >= 100: comp_is_recoil = True #req more than 0 numBP to not have bullshit ovlps giving false calls
            #if mean(prev_rsqE['rcoords']) < mean(cur_rsqE['rcoords']):        comp_is_recoil = True
            if strand_isfv(prev_rsqE): # aln is fv, then recoil is left on ref
                if mean(prev_rsqE['rcoords']) > mean(cur_rsqE['rcoords']):        comp_is_recoil = True
            if strand_isrv(prev_rsqE): # aln is rv, then recoil is right on ref
                if mean(prev_rsqE['rcoords']) < mean(cur_rsqE['rcoords']):        comp_is_recoil = True
            if comp_is_recoil:
                recoil_qSpan = [prev_rsqE['qcoords'][-1],cur_rsqE['qcoords'][0]]
                recoil_rcoords_pan = prev_rsqE['rcoords'] + cur_rsqE['rcoords']
                # Calc span on r for the recoil. Correct it for out-of-DUP-boundary
                recoil_rSpan = [max([min(recoil_rcoords_pan),tmp_save['rcoords'][0]]),min([max(recoil_rcoords_pan),tmp_save['rcoords'][-1]])]
                tmp_recoilSave = {'qLeft_rsqE_idx':prev_rsqE['idx'],'qRight_rsqE_idx':cur_rsqE['idx'],'qcoords':recoil_qSpan,
                                  'rcoords':recoil_rSpan,'len':recoil_rSpan[-1]-recoil_rSpan[0]}
                # Save
                #DUP_recoils.append(tmp_recoilSave)
                DUP_recoils.append([recoil_rSpan[0],[tmp_recoilSave],recoil_rSpan[-1]])
        #/
        # Cluster recoils
        while True:
            breakIter = False
            for e1,recoil1 in enumerate(DUP_recoils):
                for e2,recoil2 in enumerate(DUP_recoils):
                    if e1 == e2: continue
                    if getRangeOvlp(recoil1,recoil2) >= max( [recoil1[-1]-recoil1[0],recoil2[-1]-recoil2[0]] )*0.9:
                        # merge members
                        recoil1[1] += recoil2[1]
                        # Update coords
                        recoil1[0] = int(mean([recoil1[0],recoil2[0]]))
                        recoil1[-1] = int(mean([recoil1[-1],recoil2[-1]]))
                        del DUP_recoils[e2]
                        breakIter = True
                        break
                if breakIter: break
            if not breakIter: break
                        
        #/
        # Compile and save
        tmp_recoils = []
        for recoil in DUP_recoils:
            tmp_recoils.append( {'rcoords':[recoil[0],recoil[-1]],'entries':recoil[1]} )
        #/
        
        tmp_save['recoils'] = DUP_recoils
        ##/
        
        # Save
        DUPs[DUPenum] = tmp_save
        
        # Assign rsqEs to DUP + mark them at rsqE_self
        for rsqE in DUP_rsqEs:
            # Save rsqEs at DUP
            SK = 'members'
            init(tmp_save,SK,set())
            tmp_save[SK].add(rsqE['idx'])
            
            # Save DUP at rsqE (can be multiple, due to anchors)
            SK = 'memberOfDUPs'
            init(rsqE,SK,set())
            rsqE[SK].add(DUPenum)
        #/
        
        # Mark DUP at adjClust
        SK = 'containsDUPs'
        AC = AC_idx_map[tmp_save['AC_idx']]
        init(AC,SK,set())
        AC[SK].add(DUPenum)
        
    #/

    if DUP_clusts_parsed:
        qname_rsqE_dupClusts[qname] = DUP_clusts_parsed
        qname_rsqE_DUPs[qname] = DUPs
        # Update 27march Save at idx-map
        for DUP in DUPs.values():
            assignIdx(DUP)
            DUPclust_idx_map[DUP['idx']] = DUP
        #/
    #/
    ##/3
    
DUPclust_idx_map_ROMD = rangeOverlaps_makeMappingDict(DUPclust_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
AC_idx_map_ROMD = rsqE_adjClusts_idx_map_ROMD = rangeOverlaps_makeMappingDict(AC_idx_map,1000,coordsKey='rcoords',sortByKey='rname')

# Sort ACs by qname
qname_ACs = {}
for AC_idx,AC in AC_idx_map.items():
    qname = AC['qname']
    init(qname_ACs,qname,set())
    qname_ACs[qname].add(AC_idx)
#/
# Mark AC-containedment at rsqE-level (could be moved to loop above)
for AC in AC_idx_map.values():
    if 'fullyContained' in AC:
        for member in AC['members']:
            rsqE = qname_rsq_idx_map[member]
            rsqE['AC_contained'] = AC['idx']
#/
# Save ACs qspan including potential repeatmasked rsqE in AC-switch
for AC in AC_idx_map.values():
    if 'fullyContained' in AC:
        # Find rsqE orders of AC
        AC_members_order = []
        for member in AC['members']:
            rsqE = qname_rsq_idx_map[member]
            AC_members_order.append(rsqE['order'])
        #/
        # Get left+right member which is part of another AC
        left_oAC_rsqE = None
        right_oAC_rsqE = None
        for rsqE in qname_rsqRanges[AC['qname']]:
            if not 'memberOfClust' in rsqE: continue
            if rsqE['order'] < min(AC_members_order):
                if left_oAC_rsqE == None or rsqE['order'] < left_oAC_rsqE['order']:
                    left_oAC_rsqE = rsqE
            if rsqE['order'] > max(AC_members_order):
                if right_oAC_rsqE == None or rsqE['order'] > right_oAC_rsqE['order']:
                    right_oAC_rsqE = rsqE
        #/
        # Update AC with qspan including gap in AC-switch
        if not (left_oAC_rsqE and right_oAC_rsqE): continue
        if left_oAC_rsqE['qcoords'][-1] > right_oAC_rsqE['qcoords'][0]: sys.exit('not supposed to happen: ID=uasdh1879213')
        qspan = [ left_oAC_rsqE['qcoords'][-1] , right_oAC_rsqE['qcoords'][0] ]
        AC['qspan_wGap'] = qspan
        AC['left_oAC_rsqE'] = left_oAC_rsqE
        AC['right_oAC_rsqE'] = right_oAC_rsqE
        #/
#/
    
## Check number of reads used for SV calling vs. skipped
qnames_usedOrNot = {'used':0,'notUsed':0,'tot':0}
for qname,rsqEs in qname_rsqRanges.items():
    qlen = rsqEs[0]['read_len']
    if qlen >= 50000:
        qnames_usedOrNot['tot'] += 1
        if qname in qname_rsqE_adj_clusts:
            qnames_usedOrNot['used'] += 1
        else:
            qnames_usedOrNot['notUsed'] += 1
if 1: print('Fraction of dataset >50kb used=',qnames_usedOrNot['used']/qnames_usedOrNot['tot'])
##/

### Calc rAVA/refBPspan/refBPCovSelf-supps of rsqE BPs
if 1 and 'rAVA-supp qname_rsqE_clusts?':
    
    ## Pre-flight: Cluster rBPcovs/rBPsupps
    # Make idx-map of rsqE BPs at q. Pre-flight-make idx-map of rsqE BPs at r (clustering rBPs that are very close, for mem efficiency)
    qname_rsqE_BPs_idx_map = {}
    rname_rsqE_rBPs_merge = {} #will cluster rsqE rBPs: rname -> coord_bin -> rBP + LERI. Then will link this at rsqE. Must not be perfect, its just to prevent double-parsing data from BAM later. Must not even be rsqEs involved in same event. Its really just the ref range to scan.
    rsqE_lenThresh = 1000
    rsqE_HT_margin = rsqE_lenThresh-1 # Dont call BP supp if BP is this close to read HT
    readSupp_margin = 500 #require this numBP +- BP for supp (for rAVA,refSupp,refBPcov)
    
    rAVA_supp_onlyACrsqEs = True #only choose rsqEs that have AC assigned?
    
    for qname,rsqEs in qname_rsqRanges.items():
        if len(rsqEs) < 2: continue #skip if just single rsqE
        prev_rsqE = None
        rsqE_idxs_done = set()
        for cur_rsqE in sorted(rsqEs,key=lambda x: x['qcoords'][0]):
            # Skip if too short
            if cur_rsqE['len'] < rsqE_lenThresh: continue
            # If we only want those rsqEs which have an AC clust, check for it
            if rAVA_supp_onlyACrsqEs and not 'memberOfClust' in cur_rsqE: continue

            # Only save rsqE coords if NOT BOTH current and previous are masked. (heuristic, repetitive regions hard)
            if prev_rsqE != None and (prev_rsqE['masked_covFrac'] <= 0.7 or cur_rsqE['masked_covFrac'] <= 0.7):
                for tmp_rsqE,coord_idx in ( [prev_rsqE,-1],[cur_rsqE,0] ):                    
                    # Define the BPs around each rsqE. Round it and hash, for quicker lookup in rAVA below (skip if it is read head/tail)
                    qcoord = tmp_rsqE['qcoords'][coord_idx]
                    # Skip if coord is close to read HT
                    if qcoord < rsqE_HT_margin: continue
                    if qcoord >= (tmp_rsqE['read_len']-rsqE_HT_margin): continue
                    #/
                    
                    rcoord = None
                    rLERI_self_numToggle = None
                    if qcoord == tmp_rsqE['qcoords'][0]: #we are dealing with qstart
                        if strand_isfv(tmp_rsqE):
                            rcoord = tmp_rsqE['rcoords'][0]
                            rLERI_self_numToggle = 1
                        if strand_isrv(tmp_rsqE):
                            rcoord = tmp_rsqE['rcoords'][-1]
                            rLERI_self_numToggle = -1
                    if qcoord == tmp_rsqE['qcoords'][-1]: #we are dealing with qend
                        if strand_isfv(tmp_rsqE):
                            rcoord = tmp_rsqE['rcoords'][-1]
                            rLERI_self_numToggle = -1
                        if strand_isrv(tmp_rsqE):
                            rcoord = tmp_rsqE['rcoords'][0]
                            rLERI_self_numToggle = 1
                    
                    # Define range to scan on q
                    qBPsupp_range = [qcoord-readSupp_margin,qcoord+readSupp_margin]
                    #/
                    
                    ## Check if merge rBPcov/rBPsupp to previous exist
                    # Define bin from rBP
                    rname = tmp_rsqE['rname']
                    rBP_bin = int(rcoord/1000)
                    rBP_bins_to_scan = []
                    for i in range(rBP_bin-1,rBP_bin+1): # Check adjacent bins, to handle e.g. rcoord=1950 vs rcoord=2001 cases.
                        rBP_bins_to_scan.append(i)
                    
                    # Check if previous exist
                    rBP_entry = None
                    if rname in rname_rsqE_rBPs_merge:
                        for rBP_bin_to_scan in rBP_bins_to_scan:
                            if rBP_bin_to_scan in rname_rsqE_rBPs_merge[rname]:
                                for existEntry in rname_rsqE_rBPs_merge[rname][rBP_bin_to_scan]:
                                    # Check if the rBP is close and they have same LERI assigned
                                    if abs(existEntry['rBP'] - rcoord) <= 10 and existEntry['rLERI_self_numToggle'] == rLERI_self_numToggle:
                                        rBP_entry = existEntry['idx'] #save the idx. gonna compile rname->bin->entries structure to idx_map structure post-loop.
                                        break #break on first. This clustering need not be perfect! its just to prevent double-parsing the BAM. see comment at dict init.
                            # Break if we found rBP_entry
                            if rBP_entry != None: break
                    
                    # Init new rBP entry if we didnt find an existing one
                    if rBP_entry == None:
                        # Define ranges to scan on r
                        rBPcov_range_coords = [rcoord,rcoord+(((readSupp_margin*2)+500)*rLERI_self_numToggle)] #add X numBP [+500] for leniency
                        rBPcov_range = [min(rBPcov_range_coords),max(rBPcov_range_coords)]
                        rBPcov_oRange_coords = [rcoord,rcoord+(((readSupp_margin*2)+500)*rLERI_self_numToggle*-1)]
                        rBPcov_oRange = [min(rBPcov_oRange_coords),max(rBPcov_oRange_coords)]
                        rBPsupp_range = [rcoord-readSupp_margin,rcoord+readSupp_margin]
                        #/
                        # Check if rsqE is repeatmasked, then mark it, so we wont look for rBPcov here (efficiency)
                        isMasked = False
                        if tmp_rsqE['masked_covFrac'] >= 0.7:
                            isMasked = True
                        #/
                        # Save
                        tmp_rBP_entry = {'rname':rname,'rBP':rcoord,'rLERI_self_numToggle':rLERI_self_numToggle,
                                         'rBPcov_range':rBPcov_range,'rBPsupp_range':rBPsupp_range,'isMasked':isMasked,'rBPcov_oRange':rBPcov_oRange,
                                         'supp_rBPcov':HOLDER({}),'supp_rBPspan':HOLDER({})}
                        assignIdx(tmp_rBP_entry)
                        init(rname_rsqE_rBPs_merge,rname,{})
                        init(rname_rsqE_rBPs_merge[rname],rBP_bin,[])
                        rname_rsqE_rBPs_merge[rname][rBP_bin].append(tmp_rBP_entry)
                        #/
                        # Update variable for save at q
                        rBP_entry = tmp_rBP_entry['idx']
                        #/
                        ##/
                        
                    # Save at q
                    tmp_save = {'qname':qname,'qcoord':qcoord,'qBPsupp_range':qBPsupp_range,
                                'rname':tmp_rsqE['rname'],'rcoord':rcoord,
                                'rBP_idx':rBP_entry,'supp_rAVA':HOLDER({}),
                                'masked_covFrac':tmp_rsqE['masked_covFrac'],'rsqE_idx':tmp_rsqE['idx']}
                    
                    assignIdx(tmp_save)
                    qname_rsqE_BPs_idx_map[tmp_save['idx']] = tmp_save
                    #/
                    
                    # Save at rsqE
                    SK = 'qBP_idxs'
                    init(tmp_rsqE,SK,{})
                    tmp_rsqE[SK][qcoord] = tmp_save['idx']
                    #/
                #/
                
            # Update prev rsqE
            prev_rsqE = cur_rsqE
    #/
    # Compile idx-map of rsqE BPs at r
    rname_rsqE_rBPs_merge_idx_map = {}
    for rname in rname_rsqE_rBPs_merge:
        for rBP_bin in rname_rsqE_rBPs_merge[rname]:
            for rBP_entry in rname_rsqE_rBPs_merge[rname][rBP_bin]:
                rname_rsqE_rBPs_merge_idx_map[rBP_entry['idx']] = rBP_entry
    #/
    # Traverse qBPs -> rBPs. Add rBPs->qBPs info
    SK = 'qBP_idxs'
    for qBP in qname_rsqE_BPs_idx_map.values():
        rBP_idx = qBP['rBP_idx']
        rBP_entry = rname_rsqE_rBPs_merge_idx_map[rBP_idx]
        init(rBP_entry,SK,HOLDER(set()))
        rBP_entry[SK]().add(qBP['idx'])
    #/

    # ROMD BPs on q, on refBPcov, refSpan
    qname_rsqE_BPs_idx_map_qBPsupp_ROMD = rangeOverlaps_makeMappingDict(qname_rsqE_BPs_idx_map,2000,coordsKey='qBPsupp_range',sortByKey='qname')
    qname_rsqE_BPs_idx_map_rBPsupp_ROMD = rangeOverlaps_makeMappingDict(rname_rsqE_rBPs_merge_idx_map,2000,coordsKey='rBPsupp_range',sortByKey='rname')
    qname_rsqE_BPs_idx_map_rBPcov_ROMD_noMasked = rangeOverlaps_makeMappingDict(rname_rsqE_rBPs_merge_idx_map,2000,coordsKey='rBPcov_range',sortByKey='rname',ignoreIfHasKeys=['isMasked'])
    #/
    ##/
    
    ## IMPORT rAVA
    HT_margin = 750 # require alns to be this close to head/tail of reads in ovlp.
    numBP_ovlp_thresh = 2000 # skip alns below this numBP
    numBP_banFree_req = 2000 # require this numBP NOT to be banned, e.g. by repeats
    rAVA_min_numBP_ovlp = (readSupp_margin*2) -1
    
    rAVA_idx_map = {}
    
    MM2_write_to_file = False
    MM2_readFromFile = True
    
    INFO_rAVA_stats = {'notDovetail':0, 'eitherNameNotIn_qname_rsqRanges':0,'eitherNameNotIn_SVs_reads_callMarked':0,
                       'tooLittleBanFreeNumBP':0,'successful':0,
                       'lowOvlpNumBP':0,'traversed':0,'rname==qname':0}
    
    MM2_read_output_fo  = gzip.open(MM2_RAVA_saveFile,'rb')
    for ln,line_raw in enumerate( MM2_read_output_fo ): ### USE THIS IS WE READ WRITTEN FILE!!! ## awk ' { if ((($4-$3)>=3000) && (($9-$8)>=3000)) { print } } '
        if ln%100000000 == 0: print('[rAVA_import] Processed '+str(ln)+' entries. (rAVA np.3kb has 4.000.000.000 lines)')
        if 'gzipped input':        line_raw = line_raw.decode() #Decode line if we read it from file
        
        INFO_rAVA_stats['traversed'] += 1
        # Parse
        line = line_raw.strip('\n')    
        line = line.split()
        qname,qlen,qstart,qend,strand,rname,rlen,rstart,rend,num_matches,aln_len,mapq = line[:12]
        
        if qname == rname:
            INFO_rAVA_stats['rname==qname'] += 1
            continue #skip self mappings
        
        # Insta-skip if niether qname nor rname had rsqEs which we make BP-entries from
        if not (qname in qname_rsqE_BPs_idx_map_qBPsupp_ROMD or rname in qname_rsqE_BPs_idx_map_qBPsupp_ROMD): continue
        
        # Insta-skip aln if ovlp is too little
        aln_len = int(aln_len)
        if not aln_len >= numBP_ovlp_thresh:
            INFO_rAVA_stats['lowOvlpNumBP'] += 1
            continue
        #/
        
        qlen,qstart,qend,rlen,rstart,rend = map(int,[qlen,qstart,qend,rlen,rstart,rend])
        
        qcoords = [qstart,qend]
        rcoords = [rstart,rend]
        q_alnLen = qcoords[-1] - qcoords[0]
        r_alnLen = rcoords[-1] - rcoords[0]
        #/Parse
        
        ### RUN CHECKS OF MINIMAP2 ENTRY
        # 2. Ovlp between reads must be within range to head/tail on both reads, OR one read must be contained
        # 3. Ovlp between reads must be either end-to-end or either is contained in the other one.
        # 4. Check so ovlp occur in specified region
        
        ## 2&3: Check rAVA-ovlp dist to H/T
        # find distance of rAVA-ovlp to H/T of q + r
        ovlp_q_dist_h = -getRangeOvlp([0,0],qcoords)
        ovlp_q_dist_t = -getRangeOvlp([qlen,qlen],qcoords)
        ovlp_q_dists_ht = [ovlp_q_dist_h,ovlp_q_dist_t]
        
        ovlp_r_dist_h = -getRangeOvlp([0,0],rcoords)
        ovlp_r_dist_t = -getRangeOvlp([rlen,rlen],rcoords)
        ovlp_r_dists_ht = [ovlp_r_dist_h,ovlp_r_dist_t]
        #/
        
        # Check which regions of aln is covered by either read:
        # can have either, End<->End, or contained.
        # In other words: Require either <H/T at read> + <H/T at oRead> OR <H+T> at either
        valid_rAVA_ovlp = False
        q_isContained = False
        r_isContained = False
        isNotDoveTail = False
        # End <-> End 
        if ( min(ovlp_q_dists_ht) <= HT_margin and min(ovlp_r_dists_ht) <= HT_margin ):
            valid_rAVA_ovlp = True
        # q is contained
        if max(ovlp_q_dists_ht) <= HT_margin:
            valid_rAVA_ovlp = True
            q_isContained = True
        # r is contained
        if max(ovlp_r_dists_ht) <= HT_margin:
            valid_rAVA_ovlp = True
            r_isContained = True
            
        # Check if ovlp was not dovetail
        if not valid_rAVA_ovlp:
            isNotDoveTail = True
        #/
        ##/2&3

        
        # Pre-flight: setup rAVA
        tmp_data = {'qcoords':qcoords,'qname':qname,'qlen':qlen,
                    'rcoords':rcoords,'rname':rname,'rlen':rlen,
                    'isNotDoveTail':isNotDoveTail,'q_isContained':q_isContained,'r_isContained':r_isContained,
                    'strand':strand,'NM':int(num_matches),
                    'ovlp_q_dist_h':ovlp_q_dist_h,'ovlp_q_dist_t':ovlp_q_dist_t,'ovlp_r_dist_h':ovlp_r_dist_h,'ovlp_r_dist_t':ovlp_r_dist_t
                    }
        #/
        
        # Check for supp around BP
        BPsupp = False
        for qr in ( 'q','r' ):
            qrname = tmp_data[qr+'name']
            qrcoords = tmp_data[qr+'coords']
            
            oqr = None
            if qr == 'q':           oqr = 'r'
            elif qr == 'r':         oqr = 'q'
            oqrname = tmp_data[oqr+'name']
            
            # Check ovlps, only select those with full ovlp to the range we scan (set via skipNumBPOvlpBelow)
            ovlps = rangeOverlaps_lookup(qrcoords,qname_rsqE_BPs_idx_map_qBPsupp_ROMD,2000,map_lookup_key=qrname,skipNumBPOvlpBelow=rAVA_min_numBP_ovlp)
            if ovlps:
                BPsupp = True
        #/
        # Save if it supp'ed a BP
        if BPsupp:
            assignIdx(tmp_data)
            rAVA_idx_map[tmp_data['idx']] = tmp_data
        #/
        
    MM2_read_output_fo.close()
    
    # Map rname/qname -> rAVA_idx
    rAVA_read_mappings = {} #read -> comp_id
    for rAVA_idx,rAVA in rAVA_idx_map.items():
        for qr in ('q','r'):
            qrname = rAVA[qr+'name']
            init(rAVA_read_mappings,qrname,set())
            rAVA_read_mappings[qrname].add(rAVA_idx)
    #/
    
    if 0 and 'dump as pickle?':
        pickleSaveLoc = MM2_RAVA_saveFile+'.rAVA_idx_map.update09apr.rsqEinAConly.BPsuppMargin='+str(readSupp_margin)+'.PICKLE'
        if not os.path.exists(pickleSaveLoc): ### SAVE
            with open(pickleSaveLoc,'wb') as nf:
                pickle.dump(rAVA_idx_map,nf)
        with open(pickleSaveLoc,'rb') as f: ### LOAD
            rAVA_idx_map = pickle.load(f)
    ##/IMPORT rAVA
    
    ### Plot supps between rsqEs
    refNumBP_req = 1000 -1 #numBP from left to right end around BP. (this value divided by 2 is the numBP US or DS...)
    masked_covFrac_thresh = 0.7
    
    ## REFERENCE
    # Traverse BAM
    print('[rsqE_BP_supps] Traversing BAM to find refSupps (span BP at ref) and BPcovs (spans at self aln)')
    bam_fo = pysam.AlignmentFile(vsDmel_coordSorted_bam_file,'rb')
    for bam_ln,bam_entry in enumerate(bam_fo):
        # INFOPRINTER
        if bam_ln % 100000 == 0:
            print('\t[rsqE_BP_supps] Processed '+str(bam_ln)+' lines! (LF2000 has 7,806,000)')
        #/    
        if bam_entry.reference_id == -1: continue
        rname = bam_entry.reference_name    
        rcoords = [bam_entry.reference_start,bam_entry.reference_end]
        if (rcoords[-1]-rcoords[0]) < (rAVA_min_numBP_ovlp*0.7): continue #skip if short numBP ovlp
        
        # Get ovlps to rBPsupps
        ovlps = rangeOverlaps_lookup(rcoords,qname_rsqE_BPs_idx_map_rBPsupp_ROMD,2000,map_lookup_key=rname,skipNumBPOvlpBelow=rAVA_min_numBP_ovlp)
        for oS,rsqE_BP_idx,oE in ovlps:
            rsqE_BP_entry = rname_rsqE_rBPs_merge_idx_map[rsqE_BP_idx]
            rBP = rsqE_BP_entry['rBP']
            # calc spanRange DS of BP + US of BP
            DSsupp_range = [rcoords[0],rBP]
            USsupp_range = [rBP,rcoords[-1]]
            init(rsqE_BP_entry,'supp_rBPspan',HOLDER({}))
            saveLoc = rsqE_BP_entry['supp_rBPspan']()
            alnHashID = alnIDhasher(bam_entry)
            qname = bam_entry.query_name
            init(saveLoc,qname,[])
            saveLoc[qname].append( {'alnHashID':alnHashID,'DSsupp_range':DSsupp_range,'USsupp_range':USsupp_range} )
        
        # Get ovlps to rBPcovs. Take only alns which extend beyond the BP
        ovlps = rangeOverlaps_lookup(rcoords,qname_rsqE_BPs_idx_map_rBPcov_ROMD_noMasked,2000,map_lookup_key=rname,skipNumBPOvlpBelow=rAVA_min_numBP_ovlp)
        for oS,rsqE_BP_idx,oE in ovlps:
            rsqE_BP_entry = rname_rsqE_rBPs_merge_idx_map[rsqE_BP_idx]
            rBP = rsqE_BP_entry['rBP']
            alnHashID,readlen = alnIDhasher(bam_entry,returnReadLen=True)
            rID,rstart,rend,qstart,qend,strand,NM = alnHashID.split('-')
            rstart,rend,qstart,qend = int(rstart),int(rend),int(qstart),int(qend)
            
            save_aln_entry = False
            # Check if BPcovLERI is to left, then we expect bam aln to have at least X numBP right of BP
            if rsqE_BP_entry['rLERI_self_numToggle'] == -1:
                # Calc offset we want to find on qcoord
                qOffset_req = (rBP+1000) - rend #will be huge negative if rend spans far ahead. will be over 0 if we must have headroom on q
                if strand_isfv(strand) and (readlen-qend) >= qOffset_req:   save_aln_entry = True
                if strand_isrv(strand) and qstart >= qOffset_req:           save_aln_entry = True
            # Check if BPcovLERI is to right, then we expect bam aln to have at least X numBP left of BP
            elif rsqE_BP_entry['rLERI_self_numToggle'] == 1:
                # Calc offset we want to find on qcoord
                qOffset_req = rstart - (rBP-1000) #will be huge negative if rend spans far ahead. will be over 0 if we must have headroom on q
                if strand_isfv(strand) and qstart >= qOffset_req:           save_aln_entry = True
                if strand_isrv(strand) and (readlen-qend) >= qOffset_req:   save_aln_entry = True
                    
            if save_aln_entry:
                # calc spanRange DS of BP + US of BP
                DSsupp_range = [rcoords[0],rBP]
                USsupp_range = [rBP,rcoords[-1]]
                init(rsqE_BP_entry,'supp_rBPcov',HOLDER({}))
                saveLoc = rsqE_BP_entry['supp_rBPcov']()
                qname = bam_entry.query_name
                init(saveLoc,qname,[])
                saveLoc[qname].append( {'alnHashID':alnHashID,'DSsupp_range':DSsupp_range,'USsupp_range':USsupp_range} )
        
    bam_fo.close()
    print('\tDone!')
    #/
###/
    
#### FUNCTIONS2
def get_qcoords_rAVAsupp(qname,qcoords,return_rAVA_entries=False):
    """
    Will check rAVA_read_mappings -> rAVA_idx_map support between provided qcoords.
    Returns dict of rAVA-reads giving supp with array entries of True/False if it is with/without broken aln. Doing this way to better recapture DUPs
    """
    tmp_rAVA_supp = {}
    if qname in rAVA_read_mappings:
        for rAVA_idx in rAVA_read_mappings[qname]:
            rAVA = rAVA_idx_map[rAVA_idx]
            
            # Parse rAVA
            qname_qr = None
            rAVAread_qr = None
            if qname == rAVA['qname']:
                qname_qr = 'q'
                rAVAread_qr = 'r'
            if qname == rAVA['rname']:
                qname_qr = 'r'
                rAVAread_qr = 'q'
            
            qname_coords = rAVA[qname_qr+'coords']
            rAVAread = rAVA[rAVAread_qr+'name']
            #/
            # Check ovlp with input coords
            if getRangeOvlp(qcoords,qname_coords) >= (qcoords[-1]-qcoords[0]):
                init(tmp_rAVA_supp,rAVAread,[])
                if return_rAVA_entries:
                    tmp_rAVA_supp[rAVAread].append(rAVA)
                else:
                    tmp_rAVA_supp[rAVAread].append(not rAVA['isNotDoveTail']) #want TRUE for no haplo conflict
    
    return tmp_rAVA_supp

def get_qcoords_rAVA_bedgraph(qname,qcoords=None,bin_size=500,onlyWithoutHapConflict=None):
    if qcoords == None:
        qstart,qend = [0,qname_rsqRanges[qname][0]['read_len']]
    else:
        qstart,qend = qcoords
        
    bins = []
    for i in range(qstart,qend,bin_size):
        bins.append([i,0,i+bin_size-1])
        
    if qname in rAVA_read_mappings:
        for rAVA_idx in rAVA_read_mappings[qname]:
            rAVA = rAVA_idx_map[rAVA_idx]
            
            if onlyWithoutHapConflict != None and rAVA['isNotDoveTail']: continue
            
            # Parse rAVA
            qname_qr = None
            rAVAread_qr = None
            if qname == rAVA['qname']:
                qname_qr = 'q'
                rAVAread_qr = 'r'
            if qname == rAVA['rname']:
                qname_qr = 'r'
                rAVAread_qr = 'q'
            
            qname_coords = rAVA[qname_qr+'coords']
            rAVAread = rAVA[rAVAread_qr+'name']
            #/
            # Check ovlp with input coords
            if qcoords == None or getRangeOvlp(qcoords,qname_coords) >= (qcoords[-1]-qcoords[0]):
                for binn in bins:
                    if getRangeOvlp(binn,qname_coords) >= (bin_size-1):
                        binn[1] += 1
    return bins

def parse_rspanGTF(rname,rcoords,returnGenes=False,covFrac_cutoff=0):
    """
    Super simple function atm. Just return what types we hit. Will be used to parse tRNA-adj INSes, frac of INSes at rnas, etc.
    """
    ovlps_gfro = rangeOverlaps_lookup(rcoords,GFRO_idx_map_ROMD,1000,map_lookup_key=rname)
    
    gene_ids_parses = {}
    gene_feats_parses = {}
    for oS,idx,oE in ovlps_gfro:
        gfro = GFRO_idx_map[idx]
        gfro_covered = (oE-oS)/gfro['len']
        if gfro_covered < covFrac_cutoff: continue
        gene_id = gfro['gene_id']
        typee = gfro['type']
        init(gene_ids_parses,gene_id,{})
        gene_ids_parses[gene_id][typee] = gfro_covered
        init(gene_feats_parses,typee,{})
        gene_feats_parses[typee][gene_id] = gfro_covered
        
    if returnGenes:     return gene_ids_parses
    
    return gene_feats_parses

def parse_rspanContent(rname,rcoords,return_maskedCovFrac=False,return_subTypes=False,scanGFRO=False,returnFullCovGene=False):
    ovlps_masking = rangeOverlaps_lookup(rcoords,masking_idx_map_ROMD,100,map_lookup_key=rname)
    
    mask_content = {}
    mask_content_subtypes = {}
    maskRanges = [ [rcoords[0],'base',rcoords[-1]] ]
    for oS,idx,oE in ovlps_masking:
        numBP = getRangeOvlp(rcoords,[oS,oE])
        if numBP < 0: continue
        maskE = masking_idx_map[idx]
        maskClass = maskE['class']
        maskSubClass = maskE['repeat_id']
        
        init(mask_content,maskClass,0)
        mask_content[maskClass] += numBP
        
        init(mask_content_subtypes,maskSubClass,0)
        mask_content_subtypes[maskSubClass] += numBP
        
        # Check if return_maskedCovFrac
        if return_maskedCovFrac:
            if oE == oS: oE += 1 #fix single-base ovlps for computeOvlpRanges
            maskRanges.append([oS,'mask',oE])
    
    # Check if we want to add GFRO stats too
    if scanGFRO:
        ## Parse gfro of BP
        gfros = parse_rspanGTF(rname,rcoords)
        fullGenes = []
        #check if we want to parse out fully covered genes
        if returnFullCovGene:
            fullGenes = []
            if 'gene' in gfros:
                for gene_id,covFrac in gfros['gene'].items():
                    if covFrac >= 0.99:
                        fullGenes.append(gene_id)
            for fullGene in fullGenes:
                for key in gfros:
                    if key == 'gene': continue
                    if fullGene in gfros[key]:      del gfros[key][fullGene]
            rmEmptyy = []
            for checkEmpty,checkEmptyy in gfros.items():
                if not checkEmptyy: rmEmptyy.append(checkEmpty)
            for rmEmpty in rmEmptyy:
                del gfros[rmEmpty]
        #/
        if 'gene' in gfros:     del gfros['gene']
        if 'mRNA' in gfros and len(gfros) >= 2: del gfros['mRNA'] #remove mRNA if we have informative subentries
        if mask_content and 'mRNA' in gfros: del gfros['mRNA'] #remove mRNA if we have maskings found
        
        # Add back fullyparsed genes
        for fullGene in fullGenes:
            init(gfros,'gene',{})
            gfros['gene'][fullGene] = 1
        #/
        
        for entry in gfros:
            mask_content[entry] = 1
        # check if tRNA within X numBP
        tRNA_check = parse_rspanGTF(rname,[rcoords[0]-2000,rcoords[1]+2000])
        if 'tRNA' in tRNA_check:
            mask_content['tRNA_withinXkb'] = 1
        #/
        ##/
    #/
    
    returnContent = [mask_content]
    
    # Check if return mask subtypes
    if return_subTypes:     returnContent.append(mask_content_subtypes)
    
    # Check if return_maskedCovFrac
    if return_maskedCovFrac:
        maskOvlps = computeOvlpRanges_wIdxs(maskRanges)
        numBP_masked = 0
        for rangee in maskOvlps:
            if not 'base' in rangee[1]: continue
            if not 'mask' in rangee[1]: continue
            numBP_masked += rangee[-1]-rangee[0]
        
        returnContent.append(numBP_masked/(rcoords[-1]-rcoords[0]))
    
    # Check if return just 1 val, then parse it out of returnArr
    if len(returnContent) == 1:     return returnContent[0]
    
    return returnContent

def parse_qspanGTF(qname,qcoords,returnGenes=False,covFrac_cutoff=0,returnCovFrac=False,requireGFRO_numBPorCovFrac=[0,0]):
    """
    Super simple function atm. Just return what GTF-types we hit. Will be used to parse tRNA-adj INSes, frac of INSes at rnas, etc.
    """
    
    gene_ids_parses = {}
    gene_feats_parses = {}
    gfroRanges = [ [qcoords[0],'base',qcoords[-1]] ]
    for rsqE in qname_rsqRanges[qname]:
        if getRangeOvlp(rsqE['qcoords'],qcoords) >= 0:
            # parse gfros
            for gfroE in rsqE['gfros']():
                numBP_ovlp = getRangeOvlp(gfroE['qcoords'],qcoords)
                if numBP_ovlp < 0: continue #skip if they are apart...
                gfro_idx = gfroE['gfro_idx']
                gfro = GFRO_idx_map[gfro_idx]
                gfro_covered = numBP_ovlp/gfro['len']
                if gfro_covered < covFrac_cutoff: continue
                if not gfro_covered > requireGFRO_numBPorCovFrac[1] and numBP_ovlp < requireGFRO_numBPorCovFrac[0]: continue #skip if too little ovlp
                gene_id = gfro['gene_id']
                typee = gfro['type']
                init(gene_ids_parses,gene_id,{})
                gene_ids_parses[gene_id][typee] = gfro_covered
                init(gene_feats_parses,typee,{})
                gene_feats_parses[typee][gene_id] = gfro_covered
                
                if returnCovFrac:    
                    gfroRanges.append([gfroE['qcoords'][0],'gfro',gfroE['qcoords'][-1]])
    
    returnContent = [gene_feats_parses]
    if returnGenes:     returnContent.append(gene_ids_parses)
    
    # Parse covFrac
    if returnCovFrac:
        maskOvlps = computeOvlpRanges_wIdxs(gfroRanges)
        numBP_covered = 0
        for rangee in maskOvlps:
            if not 'base' in rangee[1]: continue
            if not 'gfro' in rangee[1]: continue
            numBP_covered += rangee[-1]-rangee[0]
        
        returnContent.append(numBP_covered/(qcoords[-1]-qcoords[0]))
    #/
    
    return returnContent

def parse_qspanContent(qname,qcoords,return_maskedCovFrac=False,scanGFRO=False,returnGFRO_covFrac=False,requireGFRO_numBPorCovFrac=[0,0],skipNumBPBelow=0):
    tmp_content = {}
    tmp_content_detailed = {}
    maskRanges = [ [qcoords[0],'base',qcoords[-1]] ]
    if qname in qname_rsqRanges:
        for rsqE in qname_rsqRanges[qname]:
            if getRangeOvlp(rsqE['qcoords'],qcoords) >= 0:
                # parse maskings
                for maskEntry in rsqE['maskings']():
                    numBP_ovlp = getRangeOvlp(maskEntry['qcoords'],qcoords)
                    if numBP_ovlp < skipNumBPBelow: continue #skip if they are apart...
                    for maskClassi in maskEntry['maskClassis']:
                        init(tmp_content,maskClassi,0)
                        tmp_content[maskClassi] += numBP_ovlp
                    
                    # detailed
                    for maskClassi,idxs in maskEntry['maskClassis_wIdxs'].items():
                        for idx in idxs:
                            maskE = masking_idx_map[idx]
                            subClass = maskE['repeat_id']
                            init(tmp_content_detailed,subClass,0)
                            tmp_content_detailed[subClass] += numBP_ovlp
                    #/
                    # Add to maskRanges. USED FOR COVFRAC-masked + parse GFROs where no maskings
                    if return_maskedCovFrac or scanGFRO:    
                        maskRanges.append([maskEntry['qcoords'][0],'mask',maskEntry['qcoords'][-1]])
    
    returnContent = [tmp_content,tmp_content_detailed]
    
    # Calc masked covfrac / gfro
    if return_maskedCovFrac or scanGFRO:
        # Parse covFrac, pre-flight gfro
        maskOvlps = computeOvlpRanges_wIdxs(maskRanges)
        numBP_masked = 0
        nonMasked_ranges = []
        for rangee in maskOvlps:
            if not 'base' in rangee[1]: continue
            if not 'mask' in rangee[1]:
                nonMasked_ranges.append([rangee[0],rangee[-1]])
            else:
                numBP_masked += rangee[-1]-rangee[0]
        
        masking_covFrac = 0
        if qcoords[0] != qcoords[-1]:       masking_covFrac = numBP_masked/(qcoords[-1]-qcoords[0])
        returnContent.append(masking_covFrac)
        #/
        # Parse & sum gfro content
        if scanGFRO:
            gene_feats_parses = {}
            gfro_numBP_cum = 0 #keep track of cumulative numBP which were gfro_covfrac'ed
            for nonMasked_rangee in nonMasked_ranges:
                rangee_gfros,gfro_covFrac = parse_qspanGTF(qname,nonMasked_rangee,returnCovFrac=True,requireGFRO_numBPorCovFrac=requireGFRO_numBPorCovFrac)
                for gfro_type,geneIDs_covFracs in rangee_gfros.items():
                    for geneID,covFrac in geneIDs_covFracs.items():
                        init(gene_feats_parses,gfro_type,{})
                        init(gene_feats_parses[gfro_type],geneID,0)
                        gene_feats_parses[gfro_type][geneID] += covFrac
                        
                # Save numBP gfro'ed
                gfro_numBP_cum += (nonMasked_rangee[-1]-nonMasked_rangee[0])*gfro_covFrac
                        
            # purge entries for genes which are >= X numBP covered
            genes_fully_covered = []
            if 'gene' in gene_feats_parses:
                for geneID in gene_feats_parses['gene']:
                    gene_has_start_and_end = False
                    try: #check if start codon + stop codon is covered at gene. Also check so we have some covFrac in gene itself. Cant be too stringent here, since short genes may be affected alot of repeat.
                        if gene_feats_parses['start_codon'][geneID] >= 1 and gene_feats_parses['stop_codon'][geneID] >= 1 and gene_feats_parses['gene'][geneID] >= 0.5:
                            gene_has_start_and_end = True
                    except:     pass
                    if gene_has_start_and_end:
                        genes_fully_covered.append(geneID)
                # Remove oEntries for genes which are fully covered
                for geneID in genes_fully_covered:
                    for feature in gene_feats_parses:
                        if feature == 'gene': continue
                        if geneID in gene_feats_parses[feature]:      del gene_feats_parses[feature][geneID]
            # Remove gene entries for genes which have oFeats covered
            if 'gene' in gene_feats_parses:
                for feature in gene_feats_parses:
                    if feature == 'gene': continue
                    for geneID in gene_feats_parses[feature]:
                        if geneID in gene_feats_parses['gene'] and not geneID in genes_fully_covered:     del gene_feats_parses['gene'][geneID]
            # Remove all oEntries for those that are "xRNA"
            xxRNAs = []
            for feature in gene_feats_parses:
                if feature.find('RNA') != -1 and not feature == 'mRNA':
                    for geneID in gene_feats_parses[feature]:
                        xxRNAs.append(geneID)
            for xxRNA in xxRNAs:
                for feature in gene_feats_parses:
                    if feature.find('RNA') == -1 and not feature == 'mRNA' and xxRNA in gene_feats_parses[feature]:
                        del gene_feats_parses[feature][xxRNA]
            # Clearup empty keys
            empty_keys = []
            for feature in gene_feats_parses:
                if not gene_feats_parses[feature]:
                    empty_keys.append(feature)
            for empty_key in empty_keys:
                del gene_feats_parses[empty_key]
                        
            returnContent += [gene_feats_parses]
            
            gfro_covFrac = 0
            if qcoords[0] != qcoords[-1]:   gfro_covFrac = gfro_numBP_cum/(qcoords[-1]-qcoords[0]) #handle single-BP fuckers
            if returnGFRO_covFrac:      returnContent.append(gfro_covFrac)
            #/
        #/
    
    return returnContent

def get_rsqElinks_BPs_data(rsqElink,reqMargin=1000):
    qname = rsqElink['qname']
    # Parse rBPcov/span/rAVA supps
    rBPs_data = {}
    for intstr in ('1','2',):
        rBP = rname_rsqE_rBPs_merge_idx_map[ rsqElink['rBP_idx'+intstr] ]
        rBPspan = rBP['supp_rBPspan']()
        rBPcov = rBP['supp_rBPcov']()
        
        BP = rsqElink['qcoord'+intstr]
        BP_range = [BP-reqMargin,BP+reqMargin]
        
        rAVAsupp_key = 'rAVA_supp'+intstr+'_margin='+str(reqMargin)
        if not rAVAsupp_key in rsqElink:
            rAVA_supp = get_qcoords_rAVAsupp(qname,BP_range)
            rsqElink[rAVAsupp_key] = HOLDER(rAVA_supp)
        else:
            rAVA_supp = rsqElink[rAVAsupp_key]()
            
        # parse numSupps by entries
        rBPspan_counted = 0
        for readname,entries in rBPspan.items():
            rBPspan_counted += len(entries)
        rBPcov_counted = 0
        for readname,entries in rBPcov.items():
            rBPcov_counted += len(entries)
        rAVA_counted = 0
        for readname,entries in rAVA_supp.items():
            rAVA_counted += len(entries)
            
        rBPs_data[intstr] = {'rBPspan':rBPspan,'rBPcov':rBPcov,'rAVA':rAVA_supp,
                             'rBPspan_counted':rBPspan_counted,'rBPcov_counted':rBPcov_counted,'rAVA_counted':rAVA_counted}
    return rBPs_data

def parse_rBGcov(rname,rcoords,skipMaskedCovFrac=None):
    """
    returns: cov_med,cov_mean,cov_numBPs
    """
    cov_numBPs = {}
    numBP_masked = 0
    for oS,idx,oE in rangeOverlaps_lookup(rcoords,BGraw_idx_map_ROMD,100,map_lookup_key=rname):
        BG = BGraw_idx_map[idx]
        if skipMaskedCovFrac != None and BG['masked_covFrac'] >= skipMaskedCovFrac:
            numBP_masked += (oE-oS)
            continue
        bg_cov = BG['cov']
        init(cov_numBPs,bg_cov,0)
        cov_numBPs[bg_cov] += (oE-oS)
        #/
    # Add upp with 0's for non-found bases (handle repeat-skipped, dont wanna add 0's for these.)
    init(cov_numBPs,0,0)
    cov_numBPs[0] += (rcoords[-1]-rcoords[0]) - (sum(cov_numBPs.values()) + numBP_masked)
    #/
    # Calc statistics
    cov_med = -1
    cov_mean = -1
    
    numBP_traversed = 0
    rcoords_cov_stacked = 0 #for mean
    for cov,numBP in sorted(cov_numBPs.items(),key=lambda x: x[1], reverse=True):
        # Parse median
        numBP_traversed += numBP
        if cov_med == -1 and numBP_traversed >= (rcoords[-1]-rcoords[0]-numBP_masked)*0.5: # must deduct numBP masked
            cov_med = cov
        #/
        # Save for mean calc
        rcoords_cov_stacked += numBP*cov
        #/
    
    cov_mean = rcoords_cov_stacked / (rcoords[-1]-rcoords[0])
    #/
    return cov_med,cov_mean,cov_numBPs

def assignHaplotype(inp_val,haplotypeE,haplotypeEmarg,ignoreError=False):
    hapClassi = None
    if inp_val == 0:
        hapClassi = '0'
    elif getRangeOvlp([0,(haplotypeE*1)-haplotypeEmarg],[inp_val]*2) >= 0:
        hapClassi = '<1'
    elif getRangeOvlp([(haplotypeE*1)-haplotypeEmarg,(haplotypeE*1)+haplotypeEmarg],[inp_val]*2) >= 0:  
        hapClassi = '1'
    elif getRangeOvlp([(haplotypeE*2)-haplotypeEmarg,(haplotypeE*2)+haplotypeEmarg],[inp_val]*2) >= 0:  
        hapClassi = '2'
    elif getRangeOvlp([(haplotypeE*3)-haplotypeEmarg,(haplotypeE*3)+haplotypeEmarg],[inp_val]*2) >= 0:  
        hapClassi = '3'
    elif getRangeOvlp([(haplotypeE*4)-haplotypeEmarg,(haplotypeE*4)+haplotypeEmarg],[inp_val]*2) >= 0:  
        hapClassi = '4'
    elif getRangeOvlp([(haplotypeE*4)+haplotypeEmarg,9999],[inp_val]*2) >= 0:  
        hapClassi = '>4'
    else:
        if not ignoreError:
            print('wtf?')
    return hapClassi
######/FUNCTIONS

### rsqElinks (ReShaped Query Entries "link" [link=grouped rsqEs indicative of same event]): Find which rsqE-links ovlp on ref + in rAVA
# Traverse rsqEs which are not masked, define links between nonmasked rsqEs
if 1 and 'IDE-clear':
    clearKey(qname_rsqE_BPs_idx_map,'qBP_link')

IDE_links_made = 0
print('Traversing rsqEs and adding information about links between rsqEs (only between non-masked and non-short)')
for ln,(qname,rsqEs) in enumerate(qname_rsqRanges.items()):
    qBP_pairs_done = set()
    prev_rsqE = None
    for cur_rsqE in sorted(rsqEs,key=lambda x: x['qcoords'][0]):
        # Check if short or masked rsqE. then continue: Dont use as prev, dont use as cur. Will want to link two rsqE(n-4) & rsqE(n), if (n-3),(n-2),(n-1) are short/masked, etc.
        if cur_rsqE['len'] < 1000: continue #skip short rsqEs, and dont consider them for prev_rsqE
        if cur_rsqE['masked_covFrac'] >= 0.7: continue #skip masked rsqEs, and dont consider them for prev_rsqE
        if rAVA_supp_onlyACrsqEs and not 'memberOfClust' in cur_rsqE: continue #skip rsqE if we only want to use those in AC clusts
        #/
        if prev_rsqE != None:
            # Check so qBP exist at both entries
            qBP_prev = None
            qBP_cur = None
            try:        qBP_prev = prev_rsqE['qBP_idxs'][prev_rsqE['qcoords'][-1]]
            except:
                if not (getRangeOvlp(prev_rsqE['qcoords'],[0]*2) < 1000 and getRangeOvlp(prev_rsqE['qcoords'],[prev_rsqE['read_len']]*2) < 1000): sys.exit('BUGCHECK:98vu98ajsd')
            try:        qBP_cur = cur_rsqE['qBP_idxs'][cur_rsqE['qcoords'][0]]
            except:
                if not (getRangeOvlp(cur_rsqE['qcoords'],[0]*2) < 1000 and getRangeOvlp(cur_rsqE['qcoords'],[cur_rsqE['read_len']]*2) < 1000): sys.exit('BUGCHEC:as9di012')
            if qBP_prev == None or qBP_cur == None: continue
            #/
            
            # Init link
            qBP_idx_prev = prev_rsqE['qBP_idxs'][prev_rsqE['qcoords'][-1]]
            qBP_idx_cur = cur_rsqE['qBP_idxs'][cur_rsqE['qcoords'][0]]
            qBP_prev = qname_rsqE_BPs_idx_map[qBP_idx_prev]
            qBP_cur = qname_rsqE_BPs_idx_map[qBP_idx_cur]
            
            saveKey = 'qBP_link'
            # save cur_qBP at prev_qBP
            init(qBP_prev,saveKey,{})
            qBP_prev[saveKey] = qBP_cur['idx']
            
            # save prev_qBP at cur_qBP
            init(qBP_cur,saveKey,{})
            qBP_cur[saveKey] = qBP_prev['idx']
            #/
            IDE_links_made += 1
            
        # Update prev
        prev_rsqE = cur_rsqE
print('\tDone! \tMade '+str(IDE_links_made)+' links!')
#/

rsqElinks_idx_map = {}
rsqElinks_anchors_idx_map = {}
rsqElinks_anchors = {} # rname -> bins -> anchor_idxs

IDE_num_merged = 0
min_HT_dist_for_init = 4000 #require this numBP in HT for read to init a rsqElink call
## Go through all qnames with rsqEs that have links
for ln,(qname,rsqEs) in enumerate(qname_rsqRanges.items()):
    if ln%1000 == 0: print('[clust_rsqE_links] Processed '+str(ln)+'/'+str(len(qname_rsqRanges))+'. Number of rsqElinks: '+str(len(rsqElinks_idx_map))+'. Number of merges: '+str(IDE_num_merged))
    qBP_pairs_done = set()
    #indent
    for rsqE in rsqEs:
        #indent
        CK1 = 'qBP_idxs'
        CK2 = 'qBP_link'
        if not CK1 in rsqE: continue
        for qcoord,qBP_entry_idx in rsqE[CK1].items():
            # Check if it ovlps existing rsqElink
            exist_rsqElinks = {} # rsqElink_idx -> anchor_idx -> rBP_idx -> numBP
            qBP1 = qname_rsqE_BPs_idx_map[qBP_entry_idx]
            if not CK2 in qBP1: continue #skip if no link in qBP entry
            qBP2 = qname_rsqE_BPs_idx_map[ qBP1['qBP_link'] ]
            
            # Check if already done qBP pair
            qBPs_ID = '||'.join(map(str,sorted([qBP1['idx'],qBP2['idx']])))
            if qBPs_ID in qBP_pairs_done: continue
            qBP_pairs_done.add(qBPs_ID)
            #/
            
            # Parse spans on ref
            rBP1 = rname_rsqE_rBPs_merge_idx_map[ qBP1['rBP_idx'] ]
            rBP2 = rname_rsqE_rBPs_merge_idx_map[ qBP2['rBP_idx'] ]
            
            # Find if ovlps to existing
            for tmp_rBP in ( rBP1,rBP2 ):
                rname = tmp_rBP['rname']
                rBPcov_range = tmp_rBP['rBPcov_range']
                bins_to_scan = getRangeBins(rBPcov_range)
                if rname in rsqElinks_anchors:
                    for binn in bins_to_scan:
                        if binn in rsqElinks_anchors[rname]:
                            for existAnchor_idx in rsqElinks_anchors[rname][binn]:
                                existAnchor = rsqElinks_anchors_idx_map[existAnchor_idx]
                                if existAnchor['rLERI_self_numToggle'] == tmp_rBP['rLERI_self_numToggle']:
                                    numBP_ovlp = getRangeOvlp(rBPcov_range,existAnchor['rcoords'])
                                    if numBP_ovlp >= 100:
                                        SK1 = existAnchor['link_idx']
                                        SK2 = existAnchor['idx']
                                        init(exist_rsqElinks,SK1,{})
                                        init(exist_rsqElinks[SK1],SK2,{})
                                        exist_rsqElinks[SK1][SK2][tmp_rBP['idx']] = numBP_ovlp
            #/
            # Pre-check found rsqElink ovlps to see if current ready was already added somewhere. Then add current link too.
            qname_appended_to_exist = False
            for link_idx in exist_rsqElinks:
                if len(exist_rsqElinks[link_idx]) >= 2:
                    # Check so both current anchors ovlp different (or both) existing rsqElink anchors
                    curAnchors_existAnchors_ovlps = {}
                    existAnchors_curAnchors_ovlps = {}
                    for anchor_idx in exist_rsqElinks[link_idx]:
                        for rBP_idx in exist_rsqElinks[link_idx][anchor_idx]:
                            init(curAnchors_existAnchors_ovlps,rBP_idx,{})
                            curAnchors_existAnchors_ovlps[rBP_idx][anchor_idx] = exist_rsqElinks[link_idx][anchor_idx][rBP_idx] #copy numBP
                    
                    # Check if existing rBP_anchors and current rBP_anchors ovlp at both locations
                    if len(curAnchors_existAnchors_ovlps) >= 2 and len(exist_rsqElinks[link_idx]) >= 2:
                        qBPs_span = [ min([qBP1['qcoord'],qBP2['qcoord']])-1000 , max([qBP1['qcoord'],qBP2['qcoord']])+1000 ]
                        exist_rsqElink = rsqElinks_idx_map[link_idx]
                        
                        if qname in exist_rsqElink['reads']():
                            init_rsqElink = False
                            qname_appended_to_exist = True
                            
                            tmp_save = { 'qBPs_span_wMargin':qBPs_span }
                            saveLoc = exist_rsqElink['reads']()
                            init(saveLoc,qname,[])
                            exist_rsqElink['reads']()[qname].append( tmp_save )
                            exist_rsqElink['qBP_pairs']().append([qBP1['idx'],qBP2['idx']])
                            
                            IDE_num_merged += 1
                            break
            if qname_appended_to_exist: continue
            #/                
            # Check which rsqElinks we ovlped at both anchors
            init_rsqElink = True
            for link_idx in exist_rsqElinks:
                if len(exist_rsqElinks[link_idx]) >= 2:
                    # Check so both current anchors ovlp different (or both) existing rsqElink anchors
                    curAnchors_existAnchors_ovlps = {}
                    for anchor_idx in exist_rsqElinks[link_idx]:
                        for rBP_idx in exist_rsqElinks[link_idx][anchor_idx]:
                            init(curAnchors_existAnchors_ovlps,rBP_idx,{})
                            curAnchors_existAnchors_ovlps[rBP_idx][anchor_idx] = exist_rsqElinks[link_idx][anchor_idx][rBP_idx] #copy numBP
                    
                    # Check if existing rBP_anchors and current rBP_anchors ovlp at both locations
                    if len(curAnchors_existAnchors_ovlps) >= 2 and len(exist_rsqElinks[link_idx]) >= 2:
                        qBPs_span = [ min([qBP1['qcoord'],qBP2['qcoord']])-1000 , max([qBP1['qcoord'],qBP2['qcoord']])+1000 ]
                        
                        exist_rsqElink = rsqElinks_idx_map[link_idx]
                                                    
                        # Check if it ovlps in rAVA
                        for existRead,existRead_entries in exist_rsqElink['reads']().items():
                            # Parse rAVA idxs between reads
                            rAVA_idxs = []
                            if qname in rAVA_read_mappings and existRead in rAVA_read_mappings:     rAVA_idxs = rAVA_read_mappings[qname].intersection(rAVA_read_mappings[existRead])
                            #/
                            for existRead_entry in existRead_entries:
                                for rAVA_idx in rAVA_idxs:
                                    rAVA_entry = rAVA_idx_map[rAVA_idx]
                                    # Parse rAVA
                                    qname_qr = None
                                    existRead_qr = None
                                    if qname == rAVA_entry['qname']:
                                        qname_qr = 'q'
                                        existRead_qr = 'r'
                                    if qname == rAVA_entry['rname']:
                                        qname_qr = 'r'
                                        existRead_qr = 'q'
                                    
                                    qname_rAVAcoords = rAVA_entry[qname_qr+'coords']
                                    existRead_rAVAcoords = rAVA_entry[existRead_qr+'coords']
                                    #/
                                    # Check so rsqElink spans ovlp on both reads in rAVA  
                                    existRead_qBPs_span = existRead_entry['qBPs_span_wMargin']
                                    if ( (getRangeOvlp(qBPs_span,qname_rAVAcoords) >= qBPs_span[-1]-qBPs_span[0]) and 
                                         (getRangeOvlp(existRead_qBPs_span,existRead_rAVAcoords) >= existRead_qBPs_span[-1]-existRead_qBPs_span[0]) ) :

                                        init_rsqElink = False
                                        
                                        tmp_save = { 'qBPs_span_wMargin':qBPs_span }
                                        saveLoc = exist_rsqElink['reads']()
                                        init(saveLoc,qname,[])
                                        exist_rsqElink['reads']()[qname].append( tmp_save )
                                        exist_rsqElink['qBP_pairs']().append([qBP1['idx'],qBP2['idx']])
                                        
                                        IDE_num_merged += 1
                                        break
                                    #/
                            #/
                            if not init_rsqElink: break # Break all outer if we added
                        if not init_rsqElink: break # Break all outer if we added
                        #/
                if not init_rsqElink: break # Break all outer if we added
                #/
            #/
            
            # Check if we did not append to existing rsqElink, then init new
            if init_rsqElink:
                # Check if HT margins are too low, then skip init
                qspan = [min([qBP1['qcoord'],qBP2['qcoord']]),max([qBP1['qcoord'],qBP2['qcoord']])]
                if qspan[0] < min_HT_dist_for_init or (qname_rsqRanges[qname][0]['read_len']-qspan[-1] < min_HT_dist_for_init):
                    
                    continue
                
                # Init rsqElink anchors
                tmp_anchorIdxs = set() #keep track of such idxs, for rsqElink save below, so we can add link-idx into them
                for tmp_rBP,tmp_qBP in ( [rBP1,qBP1],[rBP2,qBP2] ):
                    rname = tmp_rBP['rname']
                    rBPcov_range = tmp_rBP['rBPcov_range']
                    
                    tmp_anchor = {'rname':rname,'rcoords':rBPcov_range,'rBP_idx':tmp_rBP['idx'],'rLERI_self_numToggle':tmp_rBP['rLERI_self_numToggle'],
                                  'qBP_idx':tmp_qBP['idx']}
                    assignIdx(tmp_anchor)
                    
                    # Save at IDX map
                    rsqElinks_anchors_idx_map[tmp_anchor['idx']] = tmp_anchor
                    
                    # Save at rname->bin dict
                    init(rsqElinks_anchors,rname,{})
                    for binn_id in getRangeBins(tmp_anchor['rcoords']):
                        init(rsqElinks_anchors[rname],binn_id,[])
                        rsqElinks_anchors[rname][binn_id].append(tmp_anchor['idx'])
                        
                    # Save for track to rsqElink entry
                    tmp_anchorIdxs.add(tmp_anchor['idx'])
                #/
                
                # Init rsqElink (double-save len for clarity + compliance)
                tmp_rsqElink = {'anchor_idxs':tmp_anchorIdxs,'reads':HOLDER({}),'qBP_pairs':HOLDER([])}
                
                tmp_rsqElink_readsupp = { 'qBPs_span_wMargin':[ min([qBP1['qcoord'],qBP2['qcoord']])-1000 , max([qBP1['qcoord'],qBP2['qcoord']])+1000 ] }
                tmp_rsqElink['reads']()[qname] = [tmp_rsqElink_readsupp]
                tmp_rsqElink['qBP_pairs']().append([qBP1['idx'],qBP2['idx']])
                
                assignIdx(tmp_rsqElink)                    
                #/
                # Save link-idx to anchors
                for anchor_idx in tmp_anchorIdxs:
                    rsqElinks_anchors_idx_map[anchor_idx]['link_idx'] = tmp_rsqElink['idx']
                #/
                # Save
                rsqElinks_idx_map[tmp_rsqElink['idx']] = tmp_rsqElink
                #/
            
            #/
        #/
    #/
##/

## Save qBP/rBP coords from best qname span at rsqElink (best US/DS span)
for rsqElink_idx,rsqElink in rsqElinks_idx_map.items():
    # parse qnames US/DS spans
    qname_USDS_stats = []
    for qBP_pair in rsqElink['qBP_pairs']():
        qBP1,qBP2 = qname_rsqE_BPs_idx_map[qBP_pair[0]],qname_rsqE_BPs_idx_map[qBP_pair[-1]]
        qname = qBP1['qname']
        rsqElink_qSpan = [min([qBP1['qcoord'],qBP2['qcoord']]),max([qBP1['qcoord'],qBP2['qcoord']])]
        qname_len = qname_rsq_idx_map[qBP1['rsqE_idx']]['read_len'] #parse read len from any rsqE
        spanH,spanT = -getRangeOvlp([0]*2,[min(rsqElink_qSpan)]*2) , -getRangeOvlp([qname_len]*2,[max(rsqElink_qSpan)]*2)
        spanHT_min = min([spanH,spanT])
        qname_USDS_stats.append( [[qname,qBP1,qBP2],[spanHT_min,spanH,spanT,qBP_pair]] )
    #/
    # select qname with highest minHTspan
    qname_entry_select = sorted(qname_USDS_stats,key=lambda x: x[1][0],reverse=True)[0]
    #/
    #
    qname,qBP1,qBP2 = qname_entry_select[0]
    rBP1,rBP2 = rname_rsqE_rBPs_merge_idx_map[qBP1['rBP_idx']],rname_rsqE_rBPs_merge_idx_map[qBP2['rBP_idx']]
    
    rsqElink['len'] = abs(rBP1['rBP']-rBP2['rBP'])
    rsqElink['len_r'] = abs(rBP1['rBP']-rBP2['rBP'])
    rsqElink['len_q'] = abs(qBP1['qcoord']-qBP2['qcoord'])
    rsqElink['qname'] = qname
    rsqElink['qcoords'] = [min([qBP1['qcoord'],qBP2['qcoord']]),max([qBP1['qcoord'],qBP2['qcoord']])]
    rsqElink['qBPs_span_wMargin'] = [ min([qBP1['qcoord'],qBP2['qcoord']])-1000 , max([qBP1['qcoord'],qBP2['qcoord']])+1000 ]
    
    saveKeys_rBP = ['rname','rBP','rBPcov_range','rLERI_self_numToggle']
    saveKeys_qBP = ['rsqE_idx','qcoord']
    for enum,(rBP,qBP) in enumerate([(rBP1,qBP1),(rBP2,qBP2)]):
        enum += 1 #have it on 1,2 format (not 0,1)
        rsqElink['rBP_idx'+str(enum)] = rBP['idx']
        rsqElink['qBP_idx'+str(enum)] = qBP['idx']
        
        # Parse content from rBPs
        for saveKey in saveKeys_rBP:
            rsqElink[saveKey+str(enum)] = rBP[saveKey]
            
        # Parse content from qBPs
        for saveKey in saveKeys_qBP:
            rsqElink[saveKey+str(enum)] = qBP[saveKey]
##/

rsqElinks_idx_map_ROMD1 = rangeOverlaps_makeMappingDict(rsqElinks_idx_map,1000,coordsKey='rBPcov_range1',sortByKey='rname1')
rsqElinks_idx_map_ROMD2 = rangeOverlaps_makeMappingDict(rsqElinks_idx_map,1000,coordsKey='rBPcov_range2',sortByKey='rname2')

## Map rsqE_idx -> rsqElink_idx
rsqE_rsqElink_map = {}
for rsqElink in rsqElinks_idx_map.values():
    for strint in ('1','2',):
        rsqE_idx = rsqElink['rsqE_idx'+strint]
        rsqE_rsqElink_map[rsqE_idx] = rsqElink['idx']
##/

## Find which rsqElinks ovlp each other at rBPs (goal: Find 1/ non-merged rsqElinks [ovlp at x2 rBPs], 2/ BPs with multiple rearrs )
clearKey(rsqElinks_idx_map,'rBPs_orsqElinks_ovlps')
IDE_num_orsqElinks = []
IDE_skip_rsqElinks = []
IDE_num_rBP_clusters = []
IDE_sizes_rBP_clusters = []
IDE_num_reasn = 0
for rsqElink in rsqElinks_idx_map.values():
    #if rsqElink['idx'] in banned_rsqElinks: continue
    # Traverse + find
    orsqElink_rBP_ovlps = {} # orsqElink_idx -> rBP idx -> self_rBP_idx
    for strint in ('1','2',):
        rname = rsqElink['rname'+strint]
        rcoords = rsqElink['rBPcov_range'+strint]
        #rcoords = [rsqElink['rBP'+strint]-1000,rsqElink['rBP'+strint]+1000]
        rLERI_self_numToggle = rsqElink['rLERI_self_numToggle'+strint]
        
        ovlps1 = rangeOverlaps_lookup(rcoords,rsqElinks_idx_map_ROMD1,1000,map_lookup_key=rname)
        ovlps2 = rangeOverlaps_lookup(rcoords,rsqElinks_idx_map_ROMD2,1000,map_lookup_key=rname)
        
        idxs_done = set()
        for oS,oIdx,oE in ovlps1+ovlps2:
            #if oIdx in banned_rsqElinks: continue #skip if banned
            if oIdx == rsqElink['idx']: continue #skip if self
            # check if already cone comparison to this orsqElink
            if oIdx in idxs_done: continue
            idxs_done.add(oIdx)
            
            # Check which rBPs we ovlped
            orsqElink = rsqElinks_idx_map[oIdx]
            for ostrint in ('1','2',):
                orname = orsqElink['rname'+ostrint]
                orcoords = orsqElink['rBPcov_range'+ostrint]
                #orcoords = [orsqElink['rBP'+ostrint]-1000,orsqElink['rBP'+ostrint]+1000]
                orLERI_self_numToggle = orsqElink['rLERI_self_numToggle'+ostrint]
                
                # check so rname is same, LERI is same, and coords ovlp
                if orname == rname and orLERI_self_numToggle == rLERI_self_numToggle and getRangeOvlp(orcoords,rcoords) >= 0:
                    init(orsqElink_rBP_ovlps,oIdx,{})
                    init(orsqElink_rBP_ovlps[oIdx],ostrint,{})
                    orsqElink_rBP_ovlps[oIdx][ostrint][strint] = getRangeOvlp(orcoords,rcoords)
                    if 0 and 'remove?':
                        ## Save at rsqElinks
                        SK1 = 'rBPs_orsqElinks_ovlps'
                        # Save at rsqElink
                        SK2 = 'rBP'+strint
                        init(rsqElink,SK1,{})
                        init(rsqElink[SK1],SK2,{})
                        rsqElink[SK1][SK2][oIdx] = 'rBP'+ostrint
                        
                        # Save at orsqElink
                        SK2 = 'rBP'+ostrint
                        init(orsqElink,SK1,{})
                        init(orsqElink[SK1],SK2,{})
                        orsqElink[SK1][SK2][rsqElink['idx']] = 'rBP'+strint
                        ##/
    IDE_num_orsqElinks.append(len(orsqElink_rBP_ovlps))
    
    ## Calc "target platforms" based on orsqElink rcoords
    rRanges = {}
    # Add self rBPs
    for strint in ('1','2',):
        rname = rsqElink['rname'+strint]
        rcoords = rsqElink['rBPcov_range'+strint]
        init(rRanges,rname,[])
        rRanges[rname].append([rcoords[0],str(rsqElink['idx'])+'||'+strint,rcoords[-1]])
    #/
    # Add orsqElink rBPs
    for orsqElink_idx in orsqElink_rBP_ovlps:
        orsqElink = rsqElinks_idx_map[orsqElink_idx]
        for ostrint in ('1','2',):
            orname = orsqElink['rname'+ostrint]
            orcoords = orsqElink['rBPcov_range'+ostrint]
            init(rRanges,orname,[])
            rRanges[rname].append([orcoords[0],str(orsqElink['idx'])+'||'+ostrint,orcoords[-1]])
    #/
    # computeOvlpRanges + make clusters
    clustersArr = [] # clusters of "platforms" in the genome
    clustersArr_rBPs_ovlp = {'numBP':0,'ovlpRanges':[]} #special handle for cases where rBPs ovlp. Can be an INV or small DUP, must handle separately
    for rname,ranges in rRanges.items():
        ovlps = computeOvlpRanges_wIdxs(ranges)
        for rangee in ovlps:
            # Sort rangee rsqElinks + rBPs
            rangee_rsqElink_rBPs = {}
            for member in rangee[1]:
                rsqElink_idx,rBPnum = member.split('||')
                rsqElink_idx = int(rsqElink_idx)
                init(rangee_rsqElink_rBPs,rsqElink_idx,set())
                rangee_rsqElink_rBPs[rsqElink_idx].add(rBPnum)
            
            # Check if rangee contained both rBPs of rsqElink_self. This can happen in INV/small DUP scenarios. Lets skip such cases now, and handle them separately later.
            if rsqElink['idx'] in rangee_rsqElink_rBPs and len(rangee_rsqElink_rBPs[rsqElink['idx']]) == 2:
                rangee_copy = deepcopy(rangee)
                rangee_copy[1][rname] = True
                clustersArr_rBPs_ovlp['ovlpRanges'].append(rangee_copy) #gonna be a WIP
                clustersArr_rBPs_ovlp['numBP'] += rangee[-1]-rangee[0]
                continue
        
            # Add to (or init) cluster
            if clustersArr and clustersArr[-1]['rname'] == rname and (clustersArr[-1]['rcoords'][-1] - rangee[0]) >= 0:
                # Update
                clustersArr[-1]['rcoords'][-1] = rangee[-1]
                for rsqElink_idx in rangee_rsqElink_rBPs:
                    for rBPnum in rangee_rsqElink_rBPs[rsqElink_idx]:
                        init(clustersArr[-1]['rsqElink_idxs'],rsqElink_idx,{})
                        init(clustersArr[-1]['rsqElink_idxs'][rsqElink_idx],rBPnum,[rangee[0],rangee[-1]]) #init with rangee. Will update end coords if traverse more
                        # update end coord of rsqElink -> rBP
                        clustersArr[-1]['rsqElink_idxs'][rsqElink_idx][rBPnum][-1] = rangee[-1]
                    
            else:
                # Init
                clustersArr.append({'rname':rname,'rcoords':[rangee[0],rangee[-1]],'rsqElink_idxs':{}})
                for rsqElink_idx in rangee_rsqElink_rBPs:
                    for rBPnum in rangee_rsqElink_rBPs[rsqElink_idx]:
                        init(clustersArr[-1]['rsqElink_idxs'],rsqElink_idx,{})
                        init(clustersArr[-1]['rsqElink_idxs'][rsqElink_idx],rBPnum,[rangee[0],rangee[-1]]) #init with rangee. Will update end coords if traverse more
    #/
    # Dict clusters
    clusters = {}
    for enum,clust in enumerate(clustersArr):
        clusters[enum] = clust
        clust['enum'] = enum
        IDE_sizes_rBP_clusters.append(clust['rcoords'][-1]-clust['rcoords'][0])
    IDE_num_rBP_clusters.append(len(clusters))
    #/
    # Sort clusts to each rsqElink + rBP
    rsqElink_rBPs_clusts = {}
    for enum,clust in clusters.items():
        for rsqElink_idx in clust['rsqElink_idxs']:
            for rBPnum in clust['rsqElink_idxs'][rsqElink_idx]:
                init(rsqElink_rBPs_clusts,rsqElink_idx,{})
                # Check if rBPnum was already assigned to clust. Then check if it has more numBP to cur clust + update
                if rBPnum in rsqElink_rBPs_clusts[rsqElink_idx]:
                    IDE_num_reasn += 1
                    curAsn_clust = clusters[rsqElink_rBPs_clusts[rsqElink_idx][rBPnum]]
                    curAsn_numBP = curAsn_clust['rsqElink_idxs'][rsqElink_idx][rBPnum][-1]-curAsn_clust['rsqElink_idxs'][rsqElink_idx][rBPnum][0]
                    newAsn_numBP = clust['rsqElink_idxs'][rsqElink_idx][rBPnum][-1]-clust['rsqElink_idxs'][rsqElink_idx][rBPnum][0]
                    if curAsn_numBP > newAsn_numBP:
                        continue
                rsqElink_rBPs_clusts[rsqElink_idx][rBPnum] = enum
    #/
    #@ -- Checkpoint --
    # If not rsqElink in rsqElink_rBPs_clusts means both rsqElink rBPs ovlp --> Flag for INV/short DUP.
    # If we still have clusters, means this INV/DUP breakpoint have oRearrs
    rsqElink_rBPs_ovlp = False
    if clustersArr_rBPs_ovlp['numBP'] >= 500:
        # Handle special scenario: Make a dummy-cluster
        IDE_skip_rsqElinks.append('rsqElink rBPs ovlp='+str(clustersArr_rBPs_ovlp['numBP'])+' ,can be INV/small DUP. skipping for now.')
        # Save empty (for compliance)
        for strint in ('1','2',):
            rsqElink['rBPcov'+strint+'_clusts'] = HOLDER({'dummy':True,'oClusts':{},'rBPnum':strint})
        continue
    #@/    
    # For each rBP "platform"-cluster in rsqElink_self, find which orsqElinks exist and then which oClusters/"oPlatforms" it links to via such orsqElinks
    rBPs_clust_oclusts = {}
    for rBPnum,clustEnum in rsqElink_rBPs_clusts[rsqElink['idx']].items():
        clust = clusters[clustEnum]
        # Parse orsqElinks
        orsqElinks_orBPs = clust['rsqElink_idxs']
        # Parse oclusts per each orsqElink
        oClust_enums = {} #oClust -> orsqElink making oClust link
        for orsqElink_idx in orsqElinks_orBPs:
            if orsqElink_idx == rsqElink['idx']: continue #skip self
            for orBPnum,oClustEnum in rsqElink_rBPs_clusts[orsqElink_idx].items():
                if orBPnum in orsqElinks_orBPs[orsqElink_idx]: continue # skip if orsqEidx orBP_self
                init(oClust_enums,oClustEnum,{})
                init(oClust_enums[oClustEnum],orsqElink_idx,set())
                oClust_enums[oClustEnum][orsqElink_idx].add(orBPnum)  
        # Save
        rBPs_clust_oclusts[rBPnum] = oClust_enums
    #/
    # Compile clust per rBP & save
    for strint in ('1','2',):
        rBP_clust = deepcopy(clusters[ rsqElink_rBPs_clusts[rsqElink['idx']][strint] ]) #hard-copy for modify
        # Remove rsqElink_self from cluster rsqElink_idxs
        del rBP_clust['rsqElink_idxs'][rsqElink['idx']]
        
        rBP_clust['rBPnum'] = strint
        rBP_clust['oClusts'] = {}
        for oClust_enum in rBPs_clust_oclusts[strint]:
            oClust = deepcopy(clusters[oClust_enum]) #hard-copy
            # Remove rsqElink_self from cluster rsqElink_idxs
            if rsqElink['idx'] in oClust['rsqElink_idxs']:       del oClust['rsqElink_idxs'][rsqElink['idx']]
            rBP_clust['oClusts'][oClust_enum] = oClust
            
        # Save at rsqElink
        rsqElink['rBPcov'+strint+'_clusts'] = HOLDER(rBP_clust)
    #/
    ##/
if IDE_skip_rsqElinks:
    print('\t[IDE] Had num='+str(len(IDE_skip_rsqElinks))+' rsqLink self rBPs with large overlaps. Print variable "IDE_skip_rsqElinks" for messages.')
##/

###/

### Analyze rsqElinks
rsqElinks_numStats = {'num_raw':0,'num_filt':0,
                      'qdists':{'<100bp':0,'<1000bp':0,'<10kb':0,'<30kb':0,'>30kb':0},
                      'insert_1kb_maskedCovFracs':{'<0.25':0,'<0.5':0,'<0.75':0,'>0.75':0},
                      'rdists':{'TRA':0,'<1000bp':0,'<10kb':0,'<100kb':0,'<1mbp':0,'>1mbp':0},
                      'rname_breakpoints':{'2L':0,'2R':0,'3L':0,'3R':0,'X':0},
                      'insert_1kb_90percClassi':{'LTR/Gypsy':0,'LTR/Copia':0,'LTR/Pao':0,'LINE':0,'Simple_repeat':0,'other':0},
                      'insert_1kb_topClassi':{'LTR/Gypsy':0,'LTR/Copia':0,'LTR/Pao':0,'LINE':0,'Simple_repeat':0,'other':0},
                      'insert_1kb_topClassi_element_lengths':{'LTR/Gypsy':HOLDER([]),'LTR/Copia':HOLDER([]),'LTR/Pao':HOLDER([]),'LINE':HOLDER([]),'Simple_repeat':HOLDER([]),'other':HOLDER([])},
                      'insert_1kb_70masked_numBP':{'LTR/Gypsy':{},'LTR/Copia':{},'LTR/Pao':{},'LINE':{},'Simple_repeat':{},'other':{}},
                      'asd':0,'fullyContained_AC':0,'isDup':0,'isRepeat_andNotDup':0,'isOther':0,'isOther_local':0,
                      'strand_maintain':0,'strand_change':0
                      }

rAVA_supp_req = 3
fullyContained_ACs_done = [] # keep track of which contained ACs were parsed
reqMargin_1 = 3000
reqMargin_2 = 1000 # will use later, just init the readcount here
for ln,rsqElink in enumerate(rsqElinks_idx_map.values()):
    rsqElinks_numStats['num_raw'] += 1
    rsqElink['isUsed'] = False
    rsqElink['isSwitchOfContainedAC'] = False
    qname = rsqElink['qname']
    
    rsqE1 = qname_rsq_idx_map[ rsqElink['rsqE_idx1'] ]
    rsqE2 = qname_rsq_idx_map[ rsqElink['rsqE_idx2'] ]

    # Check rAVA supp at BPs. require some semi-stringent to get rid of false calls and low frequency calls
    BP_passed_req = []
    BP_supps = []
    BP_rnames = set()
    for intstr in ('1','2',):
        BP_rnames.add(rsqElink['rname'+intstr])
        BP = rsqElink['qcoord'+intstr]
        BP_range = [BP-reqMargin_1,BP+reqMargin_1]
        rAVAsupp_key = 'rAVA_supp'+intstr+'_margin='+str(reqMargin_1)
        
        # Compute supp1
        if not rAVAsupp_key in rsqElink:
            rAVA_supp = get_qcoords_rAVAsupp(qname,BP_range)
            rsqElink[rAVAsupp_key] = HOLDER(rAVA_supp)
        else:
            rAVA_supp = rsqElink[rAVAsupp_key]()
        #/
        # Compute supp2 (used later)
        rAVAsupp_key = 'rAVA_supp'+intstr+'_margin='+str(reqMargin_2)
        if not rAVAsupp_key in rsqElink:
            BP_range = [BP-reqMargin_2,BP+reqMargin_2]
            rsqElink[rAVAsupp_key] = HOLDER(get_qcoords_rAVAsupp(qname,BP_range))
        #/  
        BP_passed_req.append( len(rAVA_supp)>=rAVA_supp_req )
        BP_supps.append(len(rAVA_supp))
    if not all(BP_passed_req):
        continue
    
    if not BP_rnames.intersection(set(chroms_to_plot)):
        continue
    
    if len(rsqElink['reads']()) <= 1:
        continue
    #/
    # Check if rsqElink involves an AC_switch of which the AC is fully contained (I.e. big INS)
    isSwitchOfContainedAC = False
    if ('AC_contained' in rsqE1 or 'AC_contained' in rsqE2) and not ('AC_contained' in rsqE1 and 'AC_contained' in rsqE2):
        isSwitchOfContainedAC = True
    #@ mark AC-switch status at rsqElink
    if isSwitchOfContainedAC:  rsqElink['isSwitchOfContainedAC'] = True
    #/@
    #/
    # If rsqElink passed to here, mark it as "used in analysis". <-- e.g. Sniffles DEL analysis
    rsqElink['isUsed'] = True
    #/
    
    #@ parse stats for MS text
    qSpan = [min([rsqElink['qcoord1'],rsqElink['qcoord2']]),max([rsqElink['qcoord1'],rsqElink['qcoord2']])]
    len_q = rsqElink['len_q']
    len_r = rsqElink['len_r']
    rname1 = rsqElink['rname1']
    rname2 = rsqElink['rname2']
    
    # @@ Check if contained AC, then update variables with AC values
    if isSwitchOfContainedAC:
        # Get the contained AC
        AC_idx = None
        for rsqE in (rsqE1,rsqE2,):
            if 'AC_contained' in rsqE:
                AC_idx = rsqE['AC_contained']
        AC = AC_idx_map[AC_idx]     
        #/
        # Check if it was already parsed (>90% similar to an already called fullyContained AC)
        AC_was_done = False
        for exist_AC in fullyContained_ACs_done:
            if exist_AC['rname'] == AC['rname'] and getRangeOvlp(exist_AC['rcoords'],AC['rcoords']) > (AC['rcoords'][-1]-AC['rcoords'][0])*0.9:
                AC_was_done = True
                break
        if AC_was_done: continue
        fullyContained_ACs_done.append(AC)
        #/
        
        rsqElinks_numStats['fullyContained_AC'] += 1
    
        qSpan = AC['qspan_wGap']
        len_q = qSpan[-1]-qSpan[0]
        len_r = getRangeOvlp(AC['left_oAC_rsqE']['rcoords'],AC['right_oAC_rsqE']['rcoords'])
        rname1 = AC['left_oAC_rsqE']['rname']
        rname2 = AC['right_oAC_rsqE']['rname']
    #/@@
    
    rsqElinks_numStats['num_filt'] += 1

    # count rsqElink by DUP / repetitive (r or q) / Other_local / other
    rspan = [ min([rsqElink['rBP1'],rsqElink['rBP2']]),max([rsqElink['rBP1'],rsqElink['rBP2']]) ]
    if 'memberOfDUPs' in rsqE1 and 'memberOfDUPs' in rsqE2:
        rsqElinks_numStats['isDup'] += 1
    elif ( 
        (len_q >= 1000 and parse_qspanContent(qname,qSpan,return_maskedCovFrac=True,scanGFRO=True,returnGFRO_covFrac=True)[2] >= 0.7) or
        (len_r >= 1000 and rsqElink['rname1'] == rsqElink['rname2'] and parse_rspanContent(rsqElink['rname1'],rspan,return_maskedCovFrac=True)[1] >= 0.7)
        ):
        rsqElinks_numStats['isRepeat_andNotDup'] += 1
    elif len_r < 1000 and len_q < 1000:
        rsqElinks_numStats['isOther_local'] += 1
    else:
        rsqElinks_numStats['isOther'] += 1
    #/
    
    # Check strandedness (maintained vs. changed)
    if rsqE1['strand'] == rsqE2['strand']:      rsqElinks_numStats['strand_maintain'] += 1
    else:                                       rsqElinks_numStats['strand_change'] += 1
    #/
    
    # qdists
    if len_q <= 100:         rsqElinks_numStats['qdists']['<100bp'] += 1
    elif len_q <= 1000:      rsqElinks_numStats['qdists']['<1000bp'] += 1
    elif len_q <= 10000:     rsqElinks_numStats['qdists']['<10kb'] += 1
    elif len_q <= 30000:     rsqElinks_numStats['qdists']['<30kb'] += 1
    else:                                rsqElinks_numStats['qdists']['>30kb'] += 1
    #/
    # rdists
    if rname1 != rname2:         rsqElinks_numStats['rdists']['TRA'] += 1
    elif len_r <= 1000:                      rsqElinks_numStats['rdists']['<1000bp'] += 1
    elif len_r <= 10000:                     rsqElinks_numStats['rdists']['<10kb'] += 1
    elif len_r <= 100000:                    rsqElinks_numStats['rdists']['<100kb'] += 1
    elif len_r <= 1000000:                   rsqElinks_numStats['rdists']['<1mbp'] += 1
    else:                                                rsqElinks_numStats['rdists']['>1mbp'] += 1
    #/
    # ins_masked_classis
    if len_q >= 1000:
        INScontent,INScontent_detailed,masked_covFrac,gfros_notInMask,gfros_notInMaskCovFrac = parse_qspanContent(qname,qSpan,return_maskedCovFrac=True,scanGFRO=True,returnGFRO_covFrac=True)
        if masked_covFrac < 0.25:       rsqElinks_numStats['insert_1kb_maskedCovFracs']['<0.25'] += 1
        elif masked_covFrac < 0.5:      rsqElinks_numStats['insert_1kb_maskedCovFracs']['<0.5'] += 1
        elif masked_covFrac < 0.75:     rsqElinks_numStats['insert_1kb_maskedCovFracs']['<0.75'] += 1
        else:                           rsqElinks_numStats['insert_1kb_maskedCovFracs']['>0.75'] += 1
    #/
    # number of breakpoints per rname
    for intstr in ('1','2',):
        rname = rsqElink['rname'+intstr]
        if rname in rsqElinks_numStats['rname_breakpoints']:
            rsqElinks_numStats['rname_breakpoints'][rname] += 1
    #/
    # A: Number of ins_masked >1kb: 1/ top_mask_classi, 2/ numBP inserted by each annotation, 3/ >90% annotated to one repeat-class
    # B: Element length distribution
    if len_q >= 1000:
        INScontent,INScontent_detailed,masked_covFrac,gfros_notInMask,gfros_notInMaskCovFrac = parse_qspanContent(qname,qSpan,return_maskedCovFrac=True,scanGFRO=True,returnGFRO_covFrac=True)
        if INScontent:
            top_masking,top_numBP = sorted(INScontent.items(),key=lambda x: x[1],reverse=True)[0]
            if top_masking.find('LINE') != -1:      top_masking = 'LINE'
            # A1
            if masked_covFrac >= 0.7:
                if not top_masking in rsqElinks_numStats['insert_1kb_topClassi']:     top_masking = 'other'
                rsqElinks_numStats['insert_1kb_topClassi'][top_masking] += 1
                
                # B
                rsqElinks_numStats['insert_1kb_topClassi_element_lengths'][top_masking]().append(len_q)
                #/
            #/
            # A2
            for mask,numBP in INScontent.items():
                if mask.find('LINE') != -1:      mask = 'LINE'
                if not mask in rsqElinks_numStats['insert_1kb_70masked_numBP']:     mask = 'other'
                
                for intstr in ('1','2',):
                    rname = rsqElink['rname'+intstr]
                    init(rsqElinks_numStats['insert_1kb_70masked_numBP'][mask],rname,0)
                    rsqElinks_numStats['insert_1kb_70masked_numBP'][mask][rname] += numBP/2 #add half on insert per breakpoint rname
            #/
            # A3
            if top_numBP/len_q >= 0.9:
                if not top_masking in rsqElinks_numStats['insert_1kb_90percClassi']:     top_masking = 'other'
                rsqElinks_numStats['insert_1kb_90percClassi'][top_masking] += 1
            #/
    #/
    #@/

if 0 and 'dump to file':
    with open(MS_output+'/'+'rsqElinks_numStats.pickle','wb') as nf:
        pickle.dump(MS_output,nf)
        
if 0 and 'plot element length distribution?':
    ####### VIOLIN OF INS LENDISTRI, GYPSY VS "OTHER"
    preFilt = {}
    for classification,vals in rsqElinks_numStats['insert_1kb_topClassi_element_lengths'].items():
        vals = vals()
        if classification.find('LTR') == -1:      classification = 'Other'
        filt1 = filterArr( vals, 1000, 'threshold' )
        filt2 = filterArr( filt1, 10000, 'cutoff')
        init(preFilt,classification,[])
        preFilt[classification] += filt2
    maxVal = None
    panVals = []
    for i,j in preFilt.items():
        if maxVal == None or len(j) > maxVal: maxVal = len(j)
        panVals += j
    
    elements_at_full_length = {'LTR/Gypsy':0,'LTR/Pao':0,'LTR/Copia':0}
    element_full_lengths = {'LTR/Gypsy':7469,'LTR/Pao':8832,'LTR/Copia':5183}
    violins = []
    xlabels = []
    widths = []
    plt.figure(figsize=(11,11))
    cols = ['black','blue','red','green'] # Gypsy,Pao,Copia,Other
    for enum,classification in enumerate(['LTR/Gypsy','LTR/Pao','LTR/Copia','Other']):
        vals = preFilt[classification]
        
        ## number of elements at full length
        if classification in element_full_lengths:
            for val in vals:
                if getRangeOvlp([val]*2,[element_full_lengths[classification]*0.9,element_full_lengths[classification]*1.1]) >= 0:
                    elements_at_full_length[classification] += 1
        ##/
        
        violins.append(vals)
        xlabels.append(classification)
        widths.append(len(vals)/len(panVals))
        
        violin = plt.violinplot(vals,widths=len(vals)/maxVal,showmeans=False,showmedians=False,showextrema=False)
        for pc in violin['bodies']:
            #pc.set_facecolor(cols[enum])
            pc.set_facecolor('None')
            pc.set_edgecolor(cols[enum])
            pc.set_alpha(1)
    
    # Add lines for gypsy full length
    plt.plot([0.5,1.5],[element_full_lengths['LTR/Gypsy']*0.9,element_full_lengths['LTR/Gypsy']*0.9],'black')
    plt.plot([0.5,1.5],[element_full_lengths['LTR/Gypsy']*1.1,element_full_lengths['LTR/Gypsy']*1.1],'black')
    #/
    
    # adjust the fill
    for pc in violin['bodies']:
        pc.set_facecolor('None')
        pc.set_edgecolor('black')
    
    if 1 and 'savefig':
        plt.savefig(MS_output+'/'+'element_length_distri.svg')
    #######/
###/

### COVBINS: Bedgraph
# Import bedgraph
rname_bedgraph_bins = {}
bedgraph = {}
bedgraph_binsize = 1000
BGraw_idx_map = {}
with open(bg_file,'r') as f:
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
        INS,masked_covFrac = parse_rspanContent(rname,[rstart,rend],return_maskedCovFrac=True)
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
        
        tmp_save = {'rname':rname,'rcoords':[binn_start,binn_end],'cov':valMedian,'INScall_rsqElinks':set()}
        assignIdx(tmp_save)
        BG_idx_map[tmp_save['idx']] = tmp_save
        #/
BG_idx_map_ROMD = rangeOverlaps_makeMappingDict(BG_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
#/
#
bin_step = 50 #50 for analysis here, 1000 for chromosomeOverview-bedgraph.
smoothing_add = 0 #0 for analysis here, 4500 for smoothing in genomeOverview coverage track
rname_BGcovs = {}
genomeAverageCovBins_A = []
genomeAverageCovBins_X = []
tmp_for_median_cov_est_X = []
for rname in chroms_to_plot:
    if rname in ('mitochondrion_genome',): continue
    for binnstart in range(0,refLens[rname]+bin_step,bin_step):
        rcoords = [binnstart,binnstart+bin_step]
        rcoords_scan = [rcoords[0]-smoothing_add , rcoords[-1] + smoothing_add]
        
        BG_ovlps = rangeOverlaps_lookup(rcoords_scan,BG_idx_map_ROMD,1000,map_lookup_key=rname)

        BG_ovlps_covs = []
        for oS,idx,oE in BG_ovlps:
            if oS == oE: continue
            #INS,masked_covFrac = parse_rspanContent(rname,[oS,oE],return_maskedCovFrac=True)
            #if masked_covFrac >= 0.3: continue #skip repeatmasked
            BG = BG_idx_map[idx]
            BG_ovlps_covs += [BG['cov']]*(oE-oS)
        
        if not len(BG_ovlps_covs) >= bin_step*0.8: continue
    
        binn_cov = median(BG_ovlps_covs)
        
        mask,maskCovFrac = parse_rspanContent(rname,rcoords_scan,return_maskedCovFrac=True)
        #if maskCovFrac >= 0.3 and binn_cov > 1000: continue #skip if masked AND high cov
        if binn_cov > 600 or maskCovFrac > 0.3: continue
        
        if rname == 'X':
            if getRangeOvlp([binn_cov]*2,[90,140]) >= 0:    genomeAverageCovBins_X.append(binn_cov)
        else:
            if getRangeOvlp([binn_cov]*2,[180,250]) >= 0:   genomeAverageCovBins_A.append(binn_cov)
        
        init(rname_BGcovs,rname,[])
        rname_BGcovs[rname].append([rcoords[0],binn_cov,rcoords[-1]])

genomeAverageCov_A = int(median(genomeAverageCovBins_A))
genomeAverageCov_X = int(median(genomeAverageCovBins_X))
hap_eCov = round(mean( [genomeAverageCov_A/4,genomeAverageCov_X/2] ))
#/
# Plot coverage per chrom
if 0:
    plt.rc('font',size=12)
    fig = plt.figure(figsize=(17,10))
    colmap = plt.cm.gist_rainbow(np.linspace(0,1,5)) #start,step,num_to_generate
    # For regions to calc area-under-curve for (e.g. haplo regions). Need estimated coverage of haplotype. Will calculcate AOC at each chrom and haplotype conformation coverage
    chrom_hap_AOCs = {}
    
    #/
    for enum,(rname,BGbins) in enumerate(sorted(rname_BGcovs.items(),key=lambda x: x[0])):
        if rname == '4': continue
        filt0 = [ BGbin[1] for BGbin in BGbins ]
        filt1 = filterArr( filt0, 1, 'threshold' )
        filt2 = filterArr( filt1, 1000, 'cutoff' )
        if not filt2: continue
        numBins = 200 #int(max(filt2)-min(filt2))
        #plt.hist( filt2 ,bins=numBins,label=rname,histtype='step',stacked=True,fill=False,linewidth=3,alpha=0.5)
        # Plot hist with relative frequency per chromosome (not absolute numbers)
        bin_vals,bins,patches = plt.hist( filt2, weights = ((np.zeros_like(filt2)+1)/np.array(filt2).size), bins=numBins,label=rname,histtype='step',stacked=True,fill=False,linewidth=3,alpha=0.5)
        
        # Calculate area under curve regions 
        hist_AOC = {} #region -> area-under-curve
        for bin_idx,_ in enumerate(bins):
            if bin_idx == 0: continue #backwards checking
            bin_range = [bins[bin_idx-1],bins[bin_idx]]
            bin_val = bin_vals[bin_idx-1]
            
            # check which hap range overlaps
            nearest_hap = round(mean(bin_range)/hap_eCov)
            init(chrom_hap_AOCs,rname,{})
            init(chrom_hap_AOCs[rname],nearest_hap,0)
            chrom_hap_AOCs[rname][nearest_hap] += bin_val
        #/
    # Plot haplotype shaded regions
    ax = fig.gca()
    ylim = ax.get_ylim()
    for i in range(0,8):
        plt.plot([i*hap_eCov-(0.5*hap_eCov)]*2,ylim,color='black')
        plt.plot([i*hap_eCov+(0.5*hap_eCov)]*2,ylim,color='black')
    #/
    plt.locator_params(nbins=10)
    plt.xlim([0,600])
    plt.legend()
    plt.title(str(rname)+' numVals='+str(len(filt2)))
    if 1 and 'savefig?̈́': plt.savefig(MS_output+'/'+'CNbins.svg')
    
    if 0 and 'dump to file?':
        with open(bg_file.replace('.bedgraph','')+'.'+str(bin_step)+'BP_bins.'+str(smoothing_add)+'winMarg.bedgraph','w') as nf:
            for rname in sorted(rname_BGcovs):
                for binn in sorted(rname_BGcovs[rname], key=lambda x: x[0]):
                    writeArr = [rname,binn[0],binn[-1],binn[1]]
                    nf.write('\t'.join(map(str,writeArr))+'\n')
#/
###/COVINBS:Bedgraph

### Check supps between local INSes
localINS_stats = {'rAVA_supps':{},'rBPspans':{},'ratios':{}}
for k,v in localINS_stats.items():
    for chrom in chroms:
        v[chrom] = []

MS_selectINS = [] # array of [rsqElink,qBP strint, rAVA_supp]. Gonna take haplotype number from used INSes and find S1/S3 absence/presence to filter reference-error
for ln,rsqElink in enumerate(rsqElinks_idx_map.values()):
    if not rsqElink['isUsed']: continue

    qname = rsqElink['qname']
    qSpan = [min([rsqElink['qcoord1'],rsqElink['qcoord2']]),max([rsqElink['qcoord1'],rsqElink['qcoord2']])]
    
    rsqE1 = qname_rsq_idx_map[ rsqElink['rsqE_idx1'] ]
    rsqE2 = qname_rsq_idx_map[ rsqElink['rsqE_idx2'] ]
    rBPclust1 = rsqElink['rBPcov1_clusts']()
    rBPclust2 = rsqElink['rBPcov2_clusts']() 
    if 'dummy' in rBPclust1 or 'dummy' in rBPclust2: continue #skip if reference-BP was not parsed properly

    # Run filters
    CK = 'memberOfDUPs'
    if CK in rsqE1 or CK in rsqE2: continue
    if not rsqElink['rname1'] == rsqElink['rname2']: continue
    if not rsqElink['len_q'] >= 4000: continue #skip if "short" insertion
    if not (rsqE1['rname'] in chroms_to_plot and rsqE2['rname'] in chroms_to_plot): continue # Skip if not both breakpoints are on major chromosomes
    if not getRangeOvlp(rBPclust1['rcoords'],rBPclust2['rcoords']) >= -1000: continue #expect them to be close
    if getRangeOvlp(rBPclust1['rcoords'],rBPclust2['rcoords']) > 250: continue #expect them not to have too much ovlp (mapping error)
    
    INScontent,INScontent_detailed,masked_covFrac = parse_qspanContent(qname,qSpan,return_maskedCovFrac=True)
    if not masked_covFrac >= 0.7: continue # skip if not repeatmasked at least 70%
    #/
    
    # Get coverage at breakpoints
    rBPcovs_passed = []
    for strint in ('1','2',):
        cov_med,cov_mean,cov_data = parse_rBGcov(rsqElink['rname'+strint],rsqElink['rBPcov_range'+strint],skipMaskedCovFrac=0.7)
        if rsqElink['rname'+strint] == 'X':
            cov_target = [genomeAverageCov_X*0.875,genomeAverageCov_X*1.125]
        else:
            cov_target = [genomeAverageCov_A*0.875,genomeAverageCov_A*1.125]
        if getRangeOvlp([cov_med]*2,cov_target) >= 0:
            rBPcovs_passed.append(strint)
    if not len(rBPcovs_passed) >= 2: continue #skip if not both breakpoints had expected tetraploid/diploid coverage
    #/

    # Parse rBPcov/span/rAVA supps
    rBPs_data = get_rsqElinks_BPs_data(rsqElink)
    #/
    # Parse orsqElinks supp data
    orsqElink_rBPs_data = {}
    for orsqElink_idx in set(rBPclust1['rsqElink_idxs']).union(rBPclust2['rsqElink_idxs']):
        orsqElink_rBPs_data[orsqElink_idx] = get_rsqElinks_BPs_data( rsqElinks_idx_map[orsqElink_idx] )
    #/    
    # For each rBP, sort its rBPcov reads into: 1/ rBPspan, 2/ "self_rsqElink_target", 3/ "orsqElink_oTargets"
    rBPcov1_supp = rBPs_data['1']['rBPcov']
    rBPcov2_supp = rBPs_data['2']['rBPcov']
    
    rBPspan1_supp = set(rBPcov1_supp).intersection(rBPs_data['1']['rBPspan'])
    rBPspan2_supp = set(rBPcov2_supp).intersection(rBPs_data['2']['rBPspan'])
    
    ## Pool + count
    # Pool all reads, save counts to where they come from (rBPspan/rBPcov/rAVA+sameTarget/diffTarget)
    rBPclusts_readPool = {'1':{},'2':{}}
    for strint in ('1','2',):
        saveLoc = rBPclusts_readPool[strint]
        for read,read_entries in rBPs_data[strint]['rAVA'].items():
            init(saveLoc,read, {'rBPspan':0,'rBPcov':0,'rAVA':0,'diffTarget':0} )
            saveLoc[read]['rAVA'] += len(read_entries)
            
        for read,read_entries in rBPs_data[strint]['rBPcov'].items():
            init(saveLoc,read, {'rBPspan':0,'rBPcov':0,'rAVA':0,'diffTarget':0} )
            saveLoc[read]['rBPcov'] += len(read_entries)
            
        for read,read_entries in rBPs_data[strint]['rBPspan'].items():
            init(saveLoc,read, {'rBPspan':0,'rBPcov':0,'rAVA':0,'diffTarget':0} )
            saveLoc[read]['rBPspan'] += len(read_entries)
        
        for orsqElink_idx in rsqElink['rBPcov'+strint+'_clusts']()['rsqElink_idxs']:
            if orsqElink_idx in rBPclust2['rsqElink_idxs']:
                continue #skip for now
            else:
                for orBP in rsqElink['rBPcov'+strint+'_clusts']()['rsqElink_idxs'][orsqElink_idx]:
                    for read,read_entries in orsqElink_rBPs_data[orsqElink_idx][orBP]['rAVA'].items():
                        init(saveLoc,read, {'rBPspan':0,'rBPcov':0,'rAVA':0,'diffTarget':0} )
                        saveLoc[read]['diffTarget'] += len(read_entries)
                        
    # distribute counts
    rBPclusts_counted = {'1':{'rBPspan':0,'rBPcov':0,'rAVA':0,'diffTarget':0,'rBPcov_used':0},'2':{'rBPspan':0,'rBPcov':0,'rAVA':0,'diffTarget':0,'rBPcov_used':0}}
    for strint in ('1','2',):
        saveLoc = rBPclusts_counted[strint]
        for read,counts in rBPclusts_readPool[strint].items():
            # Take only reads at rBPcov?
            if counts['rBPcov'] == 0: continue
        
            #diffTarget_count = counts['diffTarget']
            diffTarget_count = 0
            rBPcov_count = counts['rBPcov'] - diffTarget_count
            rBPspan_count = counts['rBPspan'] - diffTarget_count
            rAVA_count = counts['rAVA'] - diffTarget_count - rBPspan_count
            
            if rBPcov_count > 0:        saveLoc['rBPcov'] += rBPcov_count
            if rBPspan_count > 0:       saveLoc['rBPspan'] += rBPspan_count
            if rAVA_count > 0:          saveLoc['rAVA'] += rAVA_count
            if diffTarget_count > 0:    saveLoc['diffTarget'] += rAVA_count
            if rBPspan_count > 0 or rAVA_count > 0 or diffTarget_count > 0: saveLoc['rBPcov_used'] += max([rBPspan_count,rAVA_count,diffTarget_count])
    #/
    # Make ratios
    rBPclusts_ratios = {'1':{'rAVA_vs_rBPcov':-2,'rBPsupp_vs_rBPcov':-2,'rAVA_vs_rBPspan':-2},'2':{'rAVA_vs_rBPcov':-2,'rBPsupp_vs_rBPcov':-2,'rAVA_vs_rBPspan':-2}}
    for strint in ('1','2',):
        saveLoc = rBPclusts_ratios[strint]
        try:        saveLoc['rAVA_vs_rBPcov'] = rBPclusts_counted[strint]['rAVA'] / rBPclusts_counted[strint]['rBPcov_used']
        except:     saveLoc['rAVA_vs_rBPcov'] = -1
        try:        saveLoc['rBPsupp_vs_rBPcov'] = rBPclusts_counted[strint]['rBPspan'] / rBPclusts_counted[strint]['rBPcov_used']
        except:     saveLoc['rBPsupp_vs_rBPcov'] = -1
        try:        saveLoc['rAVA_vs_rBPspan'] = rBPclusts_counted[strint]['rAVA'] / rBPclusts_counted[strint]['rBPspan']
        except:     saveLoc['rAVA_vs_rBPspan'] = -1
    ##/
    
    # Sort orsqElinks: those that have same targets and those that have diff targets
    rBPclust1_rAVA_sameTarget = set(rBPs_data['1']['rAVA'])
    rBPclust1_rAVA_diffTargets = set()
    for orsqElink_idx in rBPclust1['rsqElink_idxs']:
        if orsqElink_idx in rBPclust2['rsqElink_idxs']:
            saveLoc = rBPclust1_rAVA_sameTarget
        else:
            saveLoc = rBPclust1_rAVA_diffTargets
        for orBP in rBPclust1['rsqElink_idxs'][orsqElink_idx]:
            saveLoc.update(orsqElink_rBPs_data[orsqElink_idx][orBP]['rAVA'])
    
    rBPclust2_rAVA_sameTarget = set(rBPs_data['2']['rAVA'])
    rBPclust2_rAVA_diffTargets = set()
    for orsqElink_idx in rBPclust2['rsqElink_idxs']:
        if orsqElink_idx in rBPclust1['rsqElink_idxs']:
            saveLoc = rBPclust2_rAVA_sameTarget
        else:
            saveLoc = rBPclust2_rAVA_diffTargets
        for orBP in rBPclust2['rsqElink_idxs'][orsqElink_idx]:
            saveLoc.update(orsqElink_rBPs_data[orsqElink_idx][orBP]['rAVA'])
    #/rsq
    
    # Save for plotting
    for strint in ('1','2',):
        localINS_stats['rAVA_supps'][rsqElink['rname1']].append(rBPclusts_counted[strint]['rAVA'])
        localINS_stats['rBPspans'][rsqElink['rname1']].append(rBPclusts_counted[strint]['rBPspan'])
        localINS_stats['ratios'][rsqElink['rname1']].append(rBPclusts_ratios[strint]['rAVA_vs_rBPcov'])
    #/
    
    # Check so reported rAVA supp correspond between breakpoints
    #if not (rBPclusts_counted['1']['rAVA'] > 1 and rBPclusts_counted['2']['rAVA'] > 1): continue #skip if not both had rAVA countes
    #if ( min([rBPclusts_counted['1']['rAVA'],rBPclusts_counted['2']['rAVA']]) / max([rBPclusts_counted['2']['rAVA'],rBPclusts_counted['2']['rAVA']]) ) <= 0.9: continue
    #/
    
    # Assign haplotype-number to INS
    AOs = {}
    for strint in ('1','2',):
        hapE = int(mean([40,40,47,51,44])) # mean of "medians per chromosome"
        hapE_marg = int(hapE/2)
        hapClassi = assignHaplotype(rBPclusts_counted[strint]['rAVA'],hapE,hapE_marg)
        AOs[strint] = hapClassi
        
    if not AOs['1'] == AOs['2']: continue #skip if not same estimated AO at both INS breakpoints
    #/
    for strint in ('1','2'):
        # Update 3march: include the qLen of insertion at same index as the rAVA supp
        MS_selectINS.append( [rsqElink,strint,rBPclusts_counted[strint],AOs[strint]] )

if 0 and 'plot?':
    for rname in chroms:
        rAVA_supp = localINS_stats['rAVA_supps'][rname]
        rAVA_supp_ratio = filterArr( localINS_stats['ratios'][rname],0,'threshold')
        rsupp = localINS_stats['rBPspans'][rname]
        fig,axes = plt.subplots(1,3,figsize=(13,7))
        fig.suptitle(rname)
        nbins=100
        axes[0].hist(rAVA_supp,bins=nbins)
        axes[0].set_title('rAVA_supp')
        axes[0].set_xlabel('n='+str(len(rAVA_supp)))
        axes[1].hist(rAVA_supp_ratio,bins=nbins)
        axes[1].set_title('rAVA_supp_ratio')
        axes[2].hist(rsupp,bins=nbins)
        axes[2].set_title('ref_supp')
        plt.tight_layout()
        plt.show()

    if 0 and 'get AO numbers, plot qlen distri per maskType, per hapClassi':
        INS_rnameLocs_by_AO = {} #rname -> coords+AO

        x4HapINS_AO_supps = {}
        for chunk in MS_selectINS:
    		# Parse & calc
            rsqElink = chunk[0]
            qBP_strint = chunk[1]
            AO_hapClassi = chunk[3]
            rname = rsqElink['rname'+qBP_strint]
            
            init(x4HapINS_AO_supps,rname,{})
            init(x4HapINS_AO_supps[rname],AO_hapClassi,0)
            x4HapINS_AO_supps[rname][AO_hapClassi] += 1
            
            RO = chunk[2]['rBPspan']
            qname = rsqElink['qname']
            qSpan = [min(rsqElink['qcoords']) , max(rsqElink['qcoords'])]
            qlen = rsqElink['len_q']
            INScontent,INScontent_detailed,masked_covFrac = parse_qspanContent(qname,qSpan,return_maskedCovFrac=True,skipNumBPBelow=1)
            if INScontent == None or not INScontent:
                INScontent = {'_none_':rsqElink['len_q']}
            maskClassi = sorted(INScontent.items(), key=lambda x: x[1], reverse=True)[0][0]

            #
            init(INS_rnameLocs_by_AO,rsqElink['rname1'],[])
            INS_rnameLocs_by_AO[rsqElink['rname1']].append([rsqElink['rname1'],rsqElink['rBP1'],rsqElink['rname2'],rsqElink['rBP2'],AO_hapClassi])
            #/
            
        if 0 and 'dump for INS + SNP AO overview':
            with open('INSes.pickle','wb') as nf:
                pickle.dump(INS_rnameLocs_by_AO,nf)

## MS 23april2021: longread INS calls vs shortread INS calls
if 0:
    # import 
    with open(MS_output+'/../cell_line_similarity/outputs'+'/'+'Slines_TE_idx_map.pickle','rb') as f:
        Slines_TE_idx_map = pickle.load(f)
    Slines_TE_idx_map_ROMD = rangeOverlaps_makeMappingDict(Slines_TE_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
    #/
    # Traverse used INS-rsqElinks, check if they 1. ovlp with S2 teflon call, 2. do not have S1/S3
    found_in_teflon = 0
    x4hapReg_INS_SlinesControlled = {}
    for INS_data in MS_selectINS:
        rsqElink = INS_data[0]
        rname = rsqElink['rname1']
        rcoords_raw = [rsqElink['rBP1'],rsqElink['rBP2']]
        rcoords = [min(rcoords_raw),max(rcoords_raw)]
        rcoords_wmarg = [rcoords[0]-100,rcoords[-1]+100]
        
        # Check ovlp to teflon calls
        tmp_samples = ('S2_our','SRR497713','SRR497721','SRR497720','SRR11076397','SRR612111',)
        tmp_keys = ('supp','span','dist',)
        sample_stats = {}
        for sample in tmp_samples:
            init(sample_stats,sample,{})
            for key in tmp_keys:
                init(sample_stats[sample],key,[])
        
        teflons = []
        for oS,idx,oE in rangeOverlaps_lookup(rcoords_wmarg,Slines_TE_idx_map_ROMD,1000,map_lookup_key=rname):
            teflon = Slines_TE_idx_map[idx]
            teflons.append(teflon)
            
            sample = teflon['sample']
            sample_stats[sample]['supp'].append(teflon['count_supp'])
            sample_stats[sample]['span'].append(teflon['count_absence'])
            
            dist = abs(getRangeOvlp([oS,oE],rcoords))
            sample_stats[sample]['dist'].append(dist)
        
        if not teflons: continue # skip if we didnt find a teflon here
        
        if 'found in S2' and max(sample_stats['S2_our']['supp']) >= 1:
        #if 'found in S2, not in S1 nor S3, refSpan in S1/S3' and (
        #    max(sample_stats['S2_our']['supp']) >= 1 and 
        #    (min(sample_stats['SRR497713']['span'])>=1 or min(sample_stats['SRR497721']['span'])>=1) and
        #    (max(sample_stats['SRR497713']['supp']) == 0 and max(sample_stats['SRR497721']['supp']) == 0)
        #        ):
            found_in_teflon += 1
            
            #init(x4hapReg_INS_SlinesControlled,rname,[])
            #x4hapReg_INS_SlinesControlled[rname].append(len(rsqElink['rAVA_supp1_margin=3000']()))
            #x4hapReg_INS_SlinesControlled[rname].append(len(rsqElink['rAVA_supp2_margin=3000']()))
            init(x4hapReg_INS_SlinesControlled,rname,{})
            hapClassi = INS_data[3]['1']
            init(x4hapReg_INS_SlinesControlled[rname],hapClassi,0)
            x4hapReg_INS_SlinesControlled[rname][hapClassi] += 1
            
            if 0 and 'printinfo?':
                print('###',rname,rcoords_wmarg,len(rsqElink['rAVA_supp1_margin=1000']()))
                for enum,teflon in enumerate(teflons):
                    print(teflon['sample'],teflon['count_supp'],teflon['count_absence'],teflon['genotype'])
            
        #/
    #/
    # plot
    if 1 and 'plot? NOTE: just uses breakpoint 1':
        for rname,vals_dict in sorted(x4hapReg_INS_SlinesControlled.items(),key=lambda x: x[0]):
            for key,val in sorted(vals_dict.items()):
                print(rname,key,val)
            #plotDict_barplot(vals_dict,title=rname,file_dest=rname+'.png')
    #/
##/
###/

### Duplications
## Add supp of DUPclusts (check supp in recoils)
for qname,DUPclusts_enumed in qname_rsqE_DUPs.items():
    for DUPenum,DUP in DUPclusts_enumed.items():
        DUP_supps = []
        for recoil in DUP['recoils']:
            for entry in recoil[1]:
                prev_rsqE = qname_rsq_idx_map[entry['qLeft_rsqE_idx']]
                cur_rsqE = qname_rsq_idx_map[entry['qRight_rsqE_idx']]
                
                recoil_qSpan = [ entry['qcoords'][0]-1000,entry['qcoords'][-1]+1000 ]
                rsqEs_dist = -getRangeOvlp(prev_rsqE['qcoords'],cur_rsqE['qcoords'])
                if rsqEs_dist < 0: continue #bullshit, should never happen
                prev_qBPrange = [ prev_rsqE['qcoords'][-1]-500,prev_rsqE['qcoords'][-1]+500 ]
                
                cur_qBPrange = [ cur_rsqE['qcoords'][0]-500,cur_rsqE['qcoords'][-1]+500 ]
                recoil_qBPspan_rAVA_supp = get_qcoords_rAVAsupp(qname,recoil_qSpan)
                prev_rsqE_qBP_rAVA_supp = get_qcoords_rAVAsupp(qname,prev_qBPrange)
                cur_rsqE_qBP_rAVA_supp = get_qcoords_rAVAsupp(qname,cur_qBPrange)
                BPs_intersection = prev_rsqE_qBP_rAVA_supp.keys() & cur_rsqE_qBP_rAVA_supp.keys()
                
                # check what exists in-between
                recoil_gapRange_maskings = {}
                recoil_gapRange_masked_covFrac = 0
                if rsqEs_dist:
                    recoil_gapRange_maskings,detailed = parse_qspanContent(qname,entry['qcoords'])
                    recoil_gapRange_masked_covFrac = sum(recoil_gapRange_maskings.values()) / rsqEs_dist
                #/
                
                # Find refBPspan at BPs
                try:
                    prev_qBP = qname_rsqE_BPs_idx_map[ prev_rsqE['qBP_idxs'][prev_rsqE['qcoords'][-1]] ]
                    prev_rBP = rname_rsqE_rBPs_merge_idx_map[prev_qBP['rBP_idx']]
                    cur_qBP = qname_rsqE_BPs_idx_map [ cur_rsqE['qBP_idxs'][cur_rsqE['qcoords'][0]] ]
                    cur_rBP = rname_rsqE_rBPs_merge_idx_map[cur_qBP['rBP_idx']]
                except:
                    print('No qBP saved for entry! (readLen,rsqEqcoords[prev+cur])=',prev_rsqE['read_len'],prev_rsqE['qcoords'],cur_rsqE['qcoords'])
                #/
                
                tmp_save = {'dist':rsqEs_dist,'supp_span':HOLDER(recoil_qBPspan_rAVA_supp),'gapRange_maskings':recoil_gapRange_maskings,'gapRange_masked_covFrac':recoil_gapRange_masked_covFrac,
                            'supp_prev':HOLDER(prev_rsqE_qBP_rAVA_supp),'supp_cur':HOLDER(cur_rsqE_qBP_rAVA_supp),
                            'supp_prev_rBPspan':prev_rBP['supp_rBPspan'],'supp_cur_rBPspan':cur_rBP['supp_rBPspan'],
                            'supp_prev_rBPcov':prev_rBP['supp_rBPcov'],'supp_cur_rBPcov':cur_rBP['supp_rBPcov']}
                entry.update(tmp_save)
                DUP_supps.append(tmp_save)
##/
## Merge

# Sort DUPclusts by qname
qname_DUPclusts = {}
for DUPclust_idx,DUPclust in DUPclust_idx_map.items():
    qname = DUPclust['qname']
    init(qname_DUPclusts,qname,set())
    qname_DUPclusts[qname].add(DUPclust_idx)

# Cluster DUPcalls by ovlp in rAVA and anchors + find supps
DUPclusts_clustered = {}
DUPclusts_cluster_pointers = {}
for DUPclust_idx,DUPclust in DUPclust_idx_map.items():
    qname = DUPclust['qname']
    if not qname in qname_rsqE_BPs_idx_map_qBPsupp_ROMD: continue #skip if no entries for this qname

    # Get rAVA supp between anchors
    rAVA_supp = get_qcoords_rAVAsupp(qname,[DUPclust['qcoords'][0]-1000,DUPclust['qcoords'][-1]+1000])
    
    # Go through rAVA supp, check which of them also call INSes. Check such inses for anchor ovlps
    tmp_DUPclust_idxs = set()
    tmp_DUPclust_idxs.add(DUPclust_idx) # Add self
    if rAVA_supp:
        for rAVAread,rAVA_data in rAVA_supp.items():
            if not all(rAVA_data): continue
            if rAVAread in qname_DUPclusts:
                for oDUPclust_idx in qname_DUPclusts[rAVAread]:
                    oDUPclust = DUPclust_idx_map[oDUPclust_idx]
                    
                    # Check so they ovlp on r and either is contained in the other (+ some leniency)
                    if DUPclust['rname'] == oDUPclust['rname'] and getRangeOvlp(DUPclust['rcoords'],oDUPclust['rcoords']) >= min([DUPclust['len'],oDUPclust['len']])*0.8:
                        tmp_DUPclust_idxs.add(oDUPclust_idx)
                    
    # Save to outer
    if tmp_DUPclust_idxs:
        # Check if we append to previously add
        existClusts = set()
        for tmp_DUPclust_idx in tmp_DUPclust_idxs:
            if tmp_DUPclust_idx in DUPclusts_cluster_pointers:
                existClusts.add(DUPclusts_cluster_pointers[tmp_DUPclust_idx])
        
        # Check if just one previous found, then append to it. Else, merge all (including previous found) and init new
        initNew = True
        if existClusts:
            if len(existClusts) == 1:
                DUPclust_clustIdx = list(existClusts)[0]
                DUPclusts_clustered[DUPclust_clustIdx]['DUPclust_idxs'].update(tmp_DUPclust_idxs)
                initNew = False
                for tmp_DUPclust_idx in tmp_DUPclust_idxs:
                    DUPclusts_cluster_pointers[tmp_DUPclust_idx] = DUPclust_clustIdx
            else:
                for existClust in existClusts:
                    tmp_DUPclust_idxs.update(DUPclusts_clustered[existClust]['DUPclust_idxs'])
                    del DUPclusts_clustered[existClust]
        if initNew:
            tmp_init = {'DUPclust_idxs':tmp_DUPclust_idxs}
            assignIdx(tmp_init)
            DUPclusts_clustered[tmp_init['idx']] = tmp_init
            for tmp_DUPclust_idx in tmp_DUPclust_idxs:
                DUPclusts_cluster_pointers[tmp_DUPclust_idx] = tmp_init['idx']
            
    #/
print('\tDone!')
#/
##/
###/

### Compute ovlpranges of DUPclusts
DC_fatten = 0
DC_rRanges = {}
for DC_idx,DC in DUPclust_idx_map.items():
    # Check so DC has supp in a recoil. If not, skip add
    add_DC = False
    for recoil_chunk in DC['recoils']:
        for recoil in recoil_chunk[1]:
            qBPs_pass = []
            for qBP in recoil['qcoords']:
                rAVA_supp = get_qcoords_rAVAsupp(DC['qname'],[qBP-2000,qBP+2000])
                if len(rAVA_supp) >= 3:
                    qBPs_pass.append(True)
                else:
                    qBPs_pass.append(False)
            if all(qBPs_pass):
                add_DC = True
    if not add_DC: continue
    #/
    
    rname = DC['rname']
    rcoords = DC['rcoords']
    init(DC_rRanges,rname,[])
    DC_rRanges[rname].append([rcoords[0]-DC_fatten,DC['idx'],rcoords[-1]+DC_fatten])

DCro_idx_map = {}
for rname,rRanges in DC_rRanges.items():
    ovlps = computeOvlpRanges_wIdxs(rRanges)
    # Compile DCro
    DCranges = []
    for rangee in ovlps:
        # check if append to old (rcoords are adjacent) or add new
        if DCranges and rangee[0] == DCranges[-1][-1]:
            existEntry = DCranges[-1]
            existEntry[1].update(set(rangee[1]))
            existEntry[-1] = rangee[-1] #update rcoord
        else:
            DCranges.append([rangee[0],set(rangee[1]),rangee[-1]])
    #/
    ## Traverse DCro's and save
    for DCrange in DCranges:
        # Calc intra-DCro individual DUP borders
        DCro_rRanges = []
        for DC_idx in DCrange[1]:
            DC = DUPclust_idx_map[DC_idx]
            DCro_rRanges.append([DC['rcoords'][0],DC['idx'],DC['rcoords'][-1]])
        DCro_ovlps = computeOvlpRanges_wIdxs(DCro_rRanges)
        
        DCro_columns = []
        for dupColumn in DCro_ovlps:
            col_size = dupColumn[-1]-dupColumn[0]
            if col_size < 250: continue
            
            # Parse BG cov in DCro
            dupColumn_cov_med,dupColumn_cov_mean,cov_data = parse_rBGcov(rname,[dupColumn[0],dupColumn[-1]],skipMaskedCovFrac=0.7)
            #/
            # Per DC_idx: Parse qname, find all alns of qnames in DCro. Calc CNest
            CNests = {}
            CNests_fvrv = {}
            for DC_idx in DCrange[1]:
                DC = DUPclust_idx_map[DC_idx]
                qname = DC['qname']
                qname_DCro_aln_ovlps = []
                qname_DCro_aln_ovlps_fvrv = {}
                for rsqE in qname_rsqRanges[qname]:
                    if getRangeOvlp(rsqE['rcoords'],dupColumn) >= (dupColumn[-1]-dupColumn[0]) or getRangeOvlp(rsqE['rcoords'],dupColumn) >= 1000:
                        qname_DCro_aln_ovlps.append([rsqE['idx'],rsqE['strand']])
                        
                        init(qname_DCro_aln_ovlps_fvrv,rsqE['strand'],set())
                        qname_DCro_aln_ovlps_fvrv[rsqE['strand']].add(rsqE['idx'])
                if not qname_DCro_aln_ovlps: continue #skip if no alns of qname ovlped this DCro perfectly
                CNest = len(qname_DCro_aln_ovlps)
                if not qname_DCro_aln_ovlps: sys.exit()
                init(CNests,CNest,{})
                CNests[CNest][qname] = qname_DCro_aln_ovlps
                for fvrv,rsqE_idxs in qname_DCro_aln_ovlps_fvrv.items():
                    CNest_fvrv = len(rsqE_idxs)
                    init(CNests_fvrv,fvrv,{})
                    init(CNests_fvrv[fvrv],CNest_fvrv,{})
                    CNests_fvrv[fvrv][CNest_fvrv][qname] = rsqE_idxs
            #/
            # Determine CNest: Take highest of count "all same direction" vs. "fv+rv CNests summarized, where each CNest must have >=2 reads" <-- to prevent false NP-dups
            fvrv_CNest_filt = {}
            for fvrv,CNests_readsData in CNests_fvrv.items():
                for CNest,readsData in CNests_readsData.items():
                    init(fvrv_CNest_filt,fvrv,{})
                    if len(readsData) < 2: continue #require at least 2 reads to validate fvrv CNest
                    # Check if this CNest has higher CNest than previous add
                    if not fvrv_CNest_filt[fvrv] or list(fvrv_CNest_filt[fvrv])[0] < CNest: #parse CN of previously added. will be only key. Check if higher
                        fvrv_CNest_filt[fvrv][CNest] = readsData
            
            fvrv_CNest_sum = 0
            for fvrv,bestCNest_data in fvrv_CNest_filt.items():
                # check if no data. this happens e.g. when only one read is called in DCro
                if not bestCNest_data: continue
                bestCNest = list(bestCNest_data)[0] # best CN will be the only key in dict, parse it
                fvrv_CNest_sum += bestCNest
            
            #@ compare CNest by fv/rv only vs. sum of fv+rv
            fvPLUSrv_CNest = 0
            if CNests:      fvPLUSrv_CNest = max(CNests.keys())
            DCro_CNest = max([fvrv_CNest_sum,fvPLUSrv_CNest])
            #@/
            #/
            # Save dupColumn
            DCro_columns.append( {'rcoords':[dupColumn[0],dupColumn[-1]],'rname':DC['rname'],'DC_idxs':set(dupColumn[1]),
                                  'cov':dupColumn_cov_med,'cov_mean':dupColumn_cov_mean,
                                  'CNest':DCro_CNest,'CNests_fvrv':HOLDER(qname_DCro_aln_ovlps_fvrv),'CNests_alns':HOLDER(CNests) })
            #/
        #/
        # Save DCro
        DCro_save = {'rname':rname,'rcoords':[DCrange[0],DCrange[-1]],'DC_idxs':set(DCrange[1]),'columns':DCro_columns}
        assignIdx(DCro_save)
        DCro_idx_map[DCro_save['idx']] = DCro_save
        #/
    ##/

# Compute DCro anchors
DCro_anchor_margin = 1000
for DCro in DCro_idx_map.values():
    # Pre-flight: Calc DCro qspans on all its members
    qnames_DCro_qspan = {}
    for DC_idx in DCro['DC_idxs']:
        DC = DUPclust_idx_map[DC_idx]
        # Parse rsqE which is left/right of DC on Q. Include "left/right" q info (so we can parse correct rBP, etc, later)
        qname = DC['qname']
        qspan = DC['qcoords']
        init(qnames_DCro_qspan,qname,[])
        tmp_qcoords = qnames_DCro_qspan[qname] + qspan
        qnames_DCro_qspan[qname] = [ min(tmp_qcoords),max(tmp_qcoords) ]
    #/
    # Parse rsqEs adjacent to DCro, for each qname DCs at DCro
    anchor_rRanges = {}
    for qname,DCro_qspan in qnames_DCro_qspan.items():
        # Parse rsqE which is left/right of DC on Q. Include "left/right" q info (so we can parse correct rBP, etc, later)
        anchor_rsqEs = []
        toggle = False
        for rsqE in sorted(qname_rsqRanges[qname],key=lambda x: x['qcoords'][0]):
            if not 'memberOfClust' in rsqE: continue #skip if not member of an AC
            
            # Check so current rsqE is NOT "fully" contained in DCro
            if not (getRangeOvlp(rsqE['qcoords'],DCro_qspan) >= (rsqE['qcoords'][-1]-rsqE['qcoords'][0])*0.9 or 
                    getRangeOvlp(rsqE['rcoords'],DCro['rcoords']) >= (rsqE['rcoords'][-1]-rsqE['rcoords'][0])*0.9):
                # Check if rsqE ovlps DCro
                rsqE_ovlps_DCro = False
                if getRangeOvlp(rsqE['qcoords'],DCro_qspan) >= 1:
                    rsqE_ovlps_DCro = True
                #/
                
                SK = 'qLeft'
                if toggle:              SK = 'qRight'
                if rsqE_ovlps_DCro:     SK = 'ovlp'                
                # Check if update previous (we didnt bump into DC rsqEs)
                if not toggle:
                    if not anchor_rsqEs:
                        anchor_rsqEs.append([rsqE,SK])
                    else:
                        anchor_rsqEs[-1] = [rsqE,SK] #update previous add
                # Check if toggle was switched, then add and break (take only first rsqE to right of DC on q)
                if toggle:
                    anchor_rsqEs.append([rsqE,SK])
                    break #break on first
                
            # Check if we "bumped" into DCro, then toggle. Must have it here, since a rsqE can be member of DCro but have region outside of it.
            if getRangeOvlp(rsqE['qcoords'],DCro_qspan) >= 1 or getRangeOvlp(rsqE['rcoords'],DCro['rcoords']) >= 1:
                toggle = True
        #/
        # Add found (if any) rsqEs to anchor_rRanges
        for rsqE,qLERI in anchor_rsqEs:
            rname = rsqE['rname']
            rstart,rend = rsqE['rcoords']
            init(anchor_rRanges,rname,[])
            anchor_rRanges[rname].append([rstart,rsqE['idx'],qLERI,rend])
        #/
    #/
    # Add DCro to ovlpranges. make sure DCro split a "merged-anchor-due-to-short-DCro" into multiple anchors
    if DCro['rname'] in anchor_rRanges:     anchor_rRanges[ DCro['rname'] ].append([ DCro['rcoords'][0],'DCro',DCro['rcoords'][-1] ])
    #/
    # Traverse anchor rRanges -> Cluster into DCro anchors
    DCro_anchors = []
    for rname,rRanges in anchor_rRanges.items():
        ovlps = computeOvlpRanges_wIdxs(rRanges)
        for rangee in ovlps:
            if 'DCro' in rangee[1]: continue #skip rangee if DCro is present
            # check if chain to previous
            if DCro_anchors and DCro_anchors[-1]['rname'] == rname and DCro_anchors[-1]['rcoords'][-1] == rangee[0]:
                DCro_anchors[-1]['rcoords'][-1] = rangee[-1] #update end
                DCro_anchors[-1]['rsqE_idxs'].update(rangee[1]) #update members
            # else, start new
            else:
                DCro_anchors.append( { 'rname':rname,'rcoords':[rangee[0],rangee[-1]],'rsqE_idxs':deepcopy(rangee[1]) } )
    #/
    # Compile DCro anchors
    DCro_anchors_dict = {}
    for anchor in DCro_anchors:
        ## Calc cov at anchor X bp from rBP
        # Get rBP + rSelfLERI from any anchor rsqE
        anchor_rBP = None
        anchor_selfLERI_numToggle = None
        for rsqE_idx,qLERIs in anchor['rsqE_idxs'].items():
            # Check if have both qLERIs. happens when DCro is short and has aln which spans over it, making one rsqE become anchors on both qLERIs.
            if len(qLERIs) >= 2 or 'ovlp' in qLERIs:
                # Check if anchor is left or right of DCro
                if mean(anchor['rcoords']) < mean(DCro['rcoords']): # LEFT
                    anchor_rBP = anchor['rcoords'][-1]
                    anchor_selfLERI_numToggle = -1
                else: # ELSE: RIGHT
                    anchor_rBP = anchor['rcoords'][0]
                    anchor_selfLERI_numToggle = 1
            # Else, parse qLERI from rsqE and calc rBP + rLERI
            else:
                qLERI = list(qLERIs)[0]
                rsqE = qname_rsq_idx_map[rsqE_idx]
                if qLERI == 'qLeft' and strand_isfv(rsqE):
                    anchor_rBP = rsqE['rcoords'][-1]
                    anchor_selfLERI_numToggle = -1
                if qLERI == 'qLeft' and strand_isrv(rsqE):
                    anchor_rBP = rsqE['rcoords'][0]
                    anchor_selfLERI_numToggle = 1
                if qLERI == 'qRight' and strand_isfv(rsqE):
                    anchor_rBP = rsqE['rcoords'][0]
                    anchor_selfLERI_numToggle = 1
                if qLERI == 'qRight' and strand_isrv(rsqE):
                    anchor_rBP = rsqE['rcoords'][-1]
                    anchor_selfLERI_numToggle = -1
            
            break #break on first
        #/
        # Grab cov
        cov_data = None
        scanRange_iter = 1
        while (cov_data == None or sum(cov_data.values()) < DCro_anchor_margin*0.7): #req X numBP non-masked
            anchor_rSpan_coords = [anchor_rBP,anchor_rBP+((DCro_anchor_margin*scanRange_iter)*anchor_selfLERI_numToggle)]
            anchor_rSpan = [min(anchor_rSpan_coords),max(anchor_rSpan_coords)]
            cov_med,cov_mean,cov_data = parse_rBGcov(anchor['rname'],anchor_rSpan,skipMaskedCovFrac=0.7)
            scanRange_iter += 1
        anchor['anchor_rSpan'] = anchor_rSpan
        anchor['cov_med'] = cov_med
        anchor['cov_mean'] = cov_mean
        #/
        ##/
        assignIdx(anchor)
        DCro_anchors_dict[anchor['idx']] = anchor
    #/
    # Save
    DCro['anchors'] = DCro_anchors_dict
    #/
#/
DCro_idx_map_ROMD = rangeOverlaps_makeMappingDict(DCro_idx_map,1000,coordsKey='rcoords',sortByKey='rname')

### Parse rSpan and qSpan content of DCros
for DCro in DCro_idx_map.values():
    # Parse rsqElinks (to parse qSpan content) + rsqE_idxs (to calc ovlpRanges on ref, traverse, and call rSpancontent + DEL) of DUPclust 
    DCro_rsqElinks = []
    DCro_rsqEs = []
    for DC_idx in DCro['DC_idxs']:
        DC = DUPclust_idx_map[DC_idx]
        for rsqE_idx in DC['members']:
            if rsqE_idx in rsqE_rsqElink_map:
                rsqElink = rsqElinks_idx_map[rsqE_rsqElink_map[rsqE_idx]]
                DCro_rsqElinks.append(rsqElink)
            DCro_rsqEs.append(qname_rsq_idx_map[rsqE_idx])
    #/
    ## Parse rSpan content
    # Setup ranges
    rRanges = []
    rRanges.append([DCro['rcoords'][0],'base',DCro['rcoords'][-1]])
    for rsqE in DCro_rsqEs:
        rRanges.append([rsqE['rcoords'][0],rsqE['idx'],rsqE['rcoords'][-1]])
    #/
    # Calc ovlps + traverse
    ovlps = computeOvlpRanges_wIdxs(rRanges)
    anch_ranges = []
    DEL_ranges = []
    refScan_ranges = []
    for rangee in ovlps:
        if not 'base' in rangee[1]:
            # check if extend previous range
            if anch_ranges and anch_ranges[-1][-1] == rangee[0]:
                anch_ranges[-1][-1] = rangee[-1]
                anch_ranges[-1][1].update(rangee[1])
            else:
                anch_ranges.append(rangee)
        if 'base' in rangee[1] and len(rangee[1]) == 1:
            DEL_ranges.append(rangee)
        if 'base' in rangee[1] and len(rangee[1]) >= 2:
            # check if extend previous range
            if refScan_ranges and refScan_ranges[-1][-1] == rangee[0]:
                refScan_ranges[-1][-1] = rangee[-1]
                refScan_ranges[-1][1].update(rangee[1])
            else:
                refScan_ranges.append(rangee)
    #/
    # Parse rSpan content
    DCro_rcontent = {}
    DCro_rcontent_detailed = {}
    masked_numBP = 0
    rNotDel_numBP = 0
    for rangee in refScan_ranges:
        rcoords = [rangee[0],rangee[-1]]
        rcontent = parse_rspanContent(DCro['rname'],rcoords,return_maskedCovFrac=True,return_subTypes=True,scanGFRO=True,returnFullCovGene=True)
        
        masked_numBP += (rangee[-1]-rangee[0])*rcontent[2] ## calc numBP which were masked. Later we will calc overall masked covFrac
        rNotDel_numBP += rangee[-1]-rangee[0]
        
        for key,value in rcontent[0].items(): ## Parse "normal"
            init(DCro_rcontent,key,0)
            DCro_rcontent[key] += value
        for key,value in rcontent[1].items(): ## Parse detailed
            init(DCro_rcontent_detailed,key,0)
            DCro_rcontent_detailed[key] += value
    #/
    # Save at DCro
    DCro['rcontent'] = DCro_rcontent
    DCro['rcontent_detailed'] = DCro_rcontent_detailed
    DCro['rmasked_covFrac'] = masked_numBP / rNotDel_numBP
    #/
    ##/
    ## Parse qSpan content
    DCro_INScontent = {}
    DCro_INScontent_detailed = {}
    DCro_INScontent_notInMask = {}
    DCro_masked_numBP = 0
    DCro_INS_numBP_tot = 0
    if 0 and 'call q content from rsqElinks':
        for rsqElink in DCro_rsqElinks:
            qname = rsqElink['qname']
            qcoordss = [rsqElink['qcoord1'],rsqElink['qcoord2']]
            qSpan = [min(qcoordss),max(qcoordss)]
            INScontent,INScontent_detailed,masked_covFrac,gfros_notInMask,gfros_covFrac = parse_qspanContent(qname,qSpan,return_maskedCovFrac=True,scanGFRO=True,returnGFRO_covFrac=True,requireGFRO_numBPorCovFrac=[100,0.6])
            DCro_masked_numBP += (qSpan[-1]-qSpan[0])*masked_covFrac
            DCro_INS_numBP_tot += (qSpan[-1]-qSpan[0])
            for classi,numBP in INScontent.items():
                init(DCro_INScontent,classi,0)
                DCro_INScontent[classi] += numBP
            for classi,numBP in INScontent_detailed.items():
                init(DCro_INScontent_detailed,classi,0)
                DCro_INScontent_detailed[classi] += numBP
            for gfro,gene_covFracs in gfros_notInMask.items():
                for gene,covFrac in gene_covFracs.items():
                    init(DCro_INScontent_notInMask,gfro,0)
                    DCro_INScontent_notInMask[gfro] += 1
    if 1 and 'call q content from rsqEs':
        # Pre-flight: parse qnames and find qspans on them
        qname_DCro_qcoordss = {}
        qname_DCro_rsqE_idxs = {}
        for rsqE in DCro_rsqEs:
            init(qname_DCro_qcoordss,rsqE['qname'],[])
            qname_DCro_qcoordss[rsqE['qname']] += rsqE['qcoords']
            init(qname_DCro_rsqE_idxs,rsqE['qname'],set())
            qname_DCro_rsqE_idxs[rsqE['qname']].add(rsqE['idx'])
        qname_DCro_qspans = {}
        for qname,qcoordss in qname_DCro_qcoordss.items():
            qname_DCro_qspans[qname] = [min(qcoordss),max(qcoordss)]
        qname_DCro_INSrsqEs = {}
        for qname,DCro_qspan in qname_DCro_qspans.items():
            for rsqE in qname_rsqRanges[qname]:
                if getRangeOvlp(rsqE['qcoords'],DCro_qspan) >= (rsqE['qcoords'][-1]-rsqE['qcoords'][0]) and not rsqE['idx'] in qname_DCro_rsqE_idxs:
                    init(qname_DCro_INSrsqEs,qname,[])
                    qname_DCro_INSrsqEs[qname].append(rsqE)
        #/
        for qname in qname_DCro_INSrsqEs:
            for rsqE in qname_DCro_INSrsqEs[qname]:
                INScontent,INScontent_detailed,masked_covFrac,gfros_notInMask,gfros_covFrac = parse_qspanContent(qname,rsqE['qcoords'],return_maskedCovFrac=True,scanGFRO=True,returnGFRO_covFrac=True,requireGFRO_numBPorCovFrac=[100,0.6])
                DCro_masked_numBP += (qSpan[-1]-qSpan[0])*masked_covFrac
                DCro_INS_numBP_tot += (qSpan[-1]-qSpan[0])
                for classi,numBP in INScontent.items():
                    init(DCro_INScontent,classi,[])
                    DCro_INScontent[classi].append(numBP)
                for classi,numBP in INScontent_detailed.items():
                    init(DCro_INScontent_detailed,classi,[])
                    DCro_INScontent_detailed[classi].append(numBP)
                for gfro,gene_covFracs in gfros_notInMask.items():
                    for gene,covFrac in gene_covFracs.items():
                        init(DCro_INScontent_notInMask,gfro,[])
                        DCro_INScontent_notInMask[gfro].append(covFrac)
    DCro['INScontent'] = DCro_INScontent
    DCro['INScontent_detailed'] = DCro_INScontent_detailed
    DCro['INScontent_gfro'] = DCro_INScontent_notInMask
    DCro['INSmasked_covFrac'] = DCro_masked_numBP / max([DCro_INS_numBP_tot,1]) #sloppy-handle of division by 0...
    ##/
###/

### INFO
# Per DCro, run checks and calc fraction of DCro pool which pass criteria
# 1/ Have at least Xkb INSed. 2/ Have an LTR/Gypsy, 3/ Have at least 1 full gene, 4/ Have at least one xxRNA,... might come more!
INFO_DCro_stats = {'num':0,'hasINS':0,'hasGypsy':0,'hasHeli':0,'hasLTR':0,'hasFullGene':0,'hasRNA':0,'has_refLTR':0,'has_refRepeat':0,'has_tRNA_withinXkb':0,'passed_any':0,'ovlps_coding':0}
INFO_DCro_sets = {'hasINS':set(),'hasGypsy':set(),'hasHeli':set(),'ovlps_coding':set()}
DCro_INSes = {}
DCro_INSes_detailed = {}
DCro_rContents = {}
DCro_rContents_detailed = {}
DCro_lens = []
DCro_lens_perChrom = {}
for DCro in DCro_idx_map.values():
    # Filter
    #if not len(DCro['columns']) >= 2: continue
    if not len(DCro['DC_idxs']) >= 2: continue
    #/
    INFO_DCro_stats['num'] += 1
    DCro_lens.append(DCro['rcoords'][-1]-DCro['rcoords'][0])
    init(DCro_lens_perChrom,DCro['rname'],[])
    DCro_lens_perChrom[DCro['rname']].append(DCro['rcoords'][-1]-DCro['rcoords'][0])
    
    # 1 + 2
    num_passed = 0
    if DCro['INScontent']:
        # check if X kb INS
        for maskClassi,numBPs in DCro['INScontent'].items():
            if max(numBPs) >= 1000:
                INFO_DCro_stats['hasINS'] += 1
                INFO_DCro_sets['hasINS'].add(DCro['idx'])
                num_passed += 1
                break
        #/
        if 'LTR/Gypsy' in DCro['INScontent'] and max(DCro['INScontent']['LTR/Gypsy']) >= 4000:
            INFO_DCro_stats['hasGypsy'] += 1
            INFO_DCro_sets['hasGypsy'].add(DCro['idx'])
            num_passed += 1
        LTR_numBP = 0
        for maskClassi,numBPs in DCro['INScontent'].items():
            if maskClassi.find('LTR/') != -1:
                LTR_numBP += max(numBPs)
        if LTR_numBP >= 1000:
            INFO_DCro_stats['hasLTR'] += 1
        
        for maskClassi,numBPs in DCro['INScontent'].items():
            if maskClassi.lower().find('helitron') == -1: continue
            if max(numBPs) >= 300:
                INFO_DCro_stats['hasHeli'] += 1
                INFO_DCro_sets['hasHeli'].add(DCro['idx'])
                num_passed += 1
                break
            
    #/      
    # 3 + 4
    if DCro['rcontent']:
        if set(['gene','exon','CDS','pseudogene','snoRNA','ncRNA','rRNA','pre_miRNA','miRNA','snRNA']).intersection(set(DCro['rcontent'])):
            INFO_DCro_stats['ovlps_coding'] += 1
            INFO_DCro_sets['ovlps_coding'].add(DCro['idx'])
            num_passed += 1
        if 'gene' in DCro['rcontent']:
            INFO_DCro_stats['hasFullGene'] += 1
            num_passed += 1
        if set(['ncRNA','snoRNA','miRNA','pre_miRNA','snRNA','rRNA']).intersection(DCro['rcontent']):
            INFO_DCro_stats['hasRNA'] += 1
            num_passed += 1
        if 'tRNA_withinXkb' in DCro['rcontent']:
            INFO_DCro_stats['has_tRNA_withinXkb'] += 1
        for maskClassi,numBP in DCro['rcontent'].items():
            if maskClassi.find('LTR') != -1:
                if numBP >= 500:
                    INFO_DCro_stats['has_refLTR'] += 1
                    num_passed += 1
                    break
            if not maskClassi in ('Simple_repeat','Low_complexity',):
                if numBP >= 500:
                    INFO_DCro_stats['has_refRepeat'] += 1
            
    #/
    if num_passed:  INFO_DCro_stats['passed_any'] += 1
    # Save pan q INSes
    for maskClassi,numBPs in DCro['INScontent'].items():
        init(DCro_INSes,maskClassi,{'num':0,'numBP':0})
        DCro_INSes[maskClassi]['num'] += 1
        DCro_INSes[maskClassi]['numBP'] += max(numBPs)
    for maskClassi,numBPs in DCro['INScontent_detailed'].items():
        init(DCro_INSes_detailed,maskClassi,{'num':0,'numBP':0})
        DCro_INSes_detailed[maskClassi]['num'] += 1
        DCro_INSes_detailed[maskClassi]['numBP'] += max(numBPs)
    #/
    # Save rcontent covered by DCro
    for maskClassi,numBP in DCro['rcontent'].items():
        init(DCro_rContents,maskClassi,{'num':0,'numBP':0})
        DCro_rContents[maskClassi]['num'] += 1
        DCro_rContents[maskClassi]['numBP'] += numBP
    for maskClassi,numBP in DCro['rcontent_detailed'].items():
        init(DCro_rContents_detailed,maskClassi,{'num':0,'numBP':0})
        DCro_rContents_detailed[maskClassi]['num'] += 1
        DCro_rContents_detailed[maskClassi]['numBP'] += numBP
    #/
    
INFO_DCro_stats_PERC = {}
for key,val in INFO_DCro_stats.items():
    if key in ('num',): continue
    INFO_DCro_stats_PERC[key]=val/INFO_DCro_stats['num']
    
if 0 and 'plot':
    plotDict_barplot(DCro_INSes,xticksrotate_deg=90,selectArrElement='num')
    plotDict_barplot(DCro_rContents,xticksrotate_deg=90,selectArrElement='num')
    plotDict_barplot(INFO_DCro_stats,xticksrotate_deg=90)
    plotDict_barplot(INFO_DCro_stats_PERC,xticksrotate_deg=90,ytickdens=21,ylim=[0,1.05],printValues=True)
    
    ##############
    filt0 = DCro_lens
    filt1 = filterArr( filt0, 0, 'threshold' )
    filt2 = filterArr( filt1, 200000, 'cap' )
    fig = plt.figure(figsize=(11,11))
    plt.hist( filt2 ,bins=300)
    plt.title('DCro_lens'+' numVals='+str(len(filt2)))
    #/////////////
    ##############
    xlabels = []
    violins = []
    for rname,DCro_lens in sorted(DCro_lens_perChrom.items(),key=lambda x: x[0]):
        filt0 = DCro_lens
        filt1 = filterArr( filt0, 0, 'threshold' )
        filt2 = filterArr( filt1, 200000, 'cap' )
        violins.append(filt2)
        xlabels.append(rname+'\nN='+str(len(filt2)))
    fig = plt.figure(figsize=(11,11))
    plt.violinplot(violins,showmedians=True)
    ax = plt.gca()
    ax.set_xticks(np.arange(len(xlabels))+1)
    ax.set_xticklabels(xlabels)
    plt.title('DCro_lens'+' numVals='+str(len(filt2)))

    #/////////////
###/INFO

### Check copy-number where rsqElinks change rname (hypothesis: TRA occur in DUPed regions)
TRA_copyNumDeviation = {} # deviation of copynumber by TRA from base ploidy
IDE_TRA_qnames = set()
for rsqElink_idx,rsqElink in rsqElinks_idx_map.items():
    if rsqElink['rname1'] == rsqElink['rname2']: continue
    if not (rsqElink['rname1'] in chroms and rsqElink['rname2'] in chroms): continue
    if not rsqElink['isUsed']: continue
    if rsqElink['isSwitchOfContainedAC']: continue

    rsqE1,rsqE2 = qname_rsq_idx_map[ rsqElink['rsqE_idx1'] ],qname_rsq_idx_map[ rsqElink['rsqE_idx2'] ]
    if not ('memberOfClust' in rsqE1 and 'memberOfClust' in rsqE2): continue

    rBPcov1 = parse_rBGcov(rsqElink['rname1'],rsqElink['rBPcov_range1'])[0] # returns: cov_med,cov_mean,cov_numBPs
    rBPcov2 = parse_rBGcov(rsqElink['rname2'],rsqElink['rBPcov_range2'])[0]
    hapNum1 = round(rBPcov1/hap_eCov)
    hapNum2 = round(rBPcov2/hap_eCov)
    dev1 = None
    dev2 = None
    if rsqElink['rname1'] in autosomes:     dev1 = hapNum1-4
    if rsqElink['rname1'] in X:             dev1 = hapNum1-2
    if rsqElink['rname2'] in autosomes:     dev2 = hapNum2-4
    if rsqElink['rname2'] in X:             dev2 = hapNum2-2
    
    try: topDev = max([dev1,dev2])
    except: topDev = 9999
    init(TRA_copyNumDeviation,topDev,0)
    TRA_copyNumDeviation[topDev] += 1
    
    IDE_TRA_qnames.add(rsqElink['qname'])

if 0 and 'dump reads to paint':
    with open('TRA_CNdev_qnamePaints.pickle','wb') as nf:
        pickle.dump(IDE_TRA_qnames,nf)
###/
"""
### Rearr hotspots from rsqElinks involved in AC switches
## Traverse rsqElinks, find rsqElinks involved in an AC switch. Take rSpans from all such, calc ovlpRanges, find those with x1 + x2 "ro ovlps". e.g. same rearr vs. diff rearr but share x1 BP
# Make ranges per rname
keepACcontainment_numBP = 10000
rsqELAC_platforms_rRanges = {} #rsqEL = rsqElink AC
for rsqElink in rsqElinks_idx_map.values():
    qname = rsqElink['qname']
    # Check if filter
    tmp_reqMargin = 2000
    BP_passed_req = []
    for intstr in ('1','2',):
        BP_rnames.add(rsqElink['rname'+intstr])
        BP = rsqElink['qcoord'+intstr]
        BP_range = [BP-tmp_reqMargin,BP+tmp_reqMargin]
        rAVAsupp_key = 'rAVA_supp'+intstr+'_margin='+str(tmp_reqMargin)
        if not rAVAsupp_key in rsqElink:
            rAVA_supp = get_qcoords_rAVAsupp(qname,BP_range)
            rsqElink[rAVAsupp_key] = HOLDER(rAVA_supp)
        else:
            rAVA_supp = rsqElink[rAVAsupp_key]()
        BP_passed_req.append( len(rAVA_supp)>=3 )
    if not all(BP_passed_req):
        continue
    if len(rsqElink['reads']()) <= 1:
        continue
    #/
      
    rsqE1 = qname_rsq_idx_map[rsqElink['rsqE_idx1']]
    rsqE2 = qname_rsq_idx_map[rsqElink['rsqE_idx2']]
    CK = 'memberOfClust'
    if not (CK in rsqE1 and CK in rsqE2): continue
    AC1 = qname_rsqE_adj_clusts[rsqE1['qname']][rsqE1['memberOfClust']]
    AC2 = qname_rsqE_adj_clusts[rsqE2['qname']][rsqE2['memberOfClust']]
    if AC1['idx'] == AC2['idx']: continue
        
    # Check if either AC is contained. Skip if it is not a large contained AC.
    if 'fullyContained' in AC1 or 'fullyContained' in AC2:
        # check how many numBP the contained AC has. Example: I'm missing 3R<->MITO for variant graph since its fullContained. Should I keep fullyContained which are large? e.g. >10kb or such. Goal is to get rid of the "small" containments
        numBP_contained_q = 0
        numBP_contained_r = 0
        for tmp_AC in ( AC1,AC2 ):
            if 'fullyContained' in tmp_AC:
                numBP_contained_q = tmp_AC['qcoords'][-1]-tmp_AC['qcoords'][0]
                numBP_contained_r = tmp_AC['rcoords'][-1]-tmp_AC['rcoords'][0]
        if keepACcontainment_numBP < 0 or (numBP_contained_q < keepACcontainment_numBP):
            continue
    #/

    for rsqE,AC in ( [rsqE1,AC1],[rsqE2,AC2] ):
        rname = AC['rname']
        rstart,rend = AC['rcoords']
        init(rsqELAC_platforms_rRanges,rname,[])
        rsqELAC_platforms_rRanges[rname].append([rstart,rsqE['idx'],rsqElink['idx'],rend])
#/
# Calc overlaps on rnames, compile into "platforms"
clearKey(rsqElinks_idx_map,'rsqELACP_idx1')
clearKey(rsqElinks_idx_map,'rsqELACP_idx2')
rsqELAC_platforms_chain_dist = 2000 # chain platforms within this dist on R, after deducting ref maskings
rsqELAC_platforms_idx_map = {}
for rname,rRanges in rsqELAC_platforms_rRanges.items():
    ovlps = computeOvlpRanges_wIdxs(rRanges)
    # parse "platforms" of rsqElink's ACs in genome
    AC_platforms = []
    for rangee in ovlps:
        ## check if append to old (rcoords are adjacent) or add new
        # Check numBP masked at ref between platforms
        rangee_existRangee_numBP_refMasked = 0
        if AC_platforms:
            prevEntry = AC_platforms[-1]
            rspan_to_scan = [prevEntry[-1],rangee[0]]
            rspan_len = rspan_to_scan[-1]-rspan_to_scan[0]
            # Only calc numBP mask at ref if within reasonable dist
            if rspan_len > rsqELAC_platforms_chain_dist and rspan_len <= 20000:
                _,masked_covFrac = parse_rspanContent(rname,rspan_to_scan,return_maskedCovFrac=True)
                rangee_existRangee_numBP_refMasked = int(masked_covFrac * rspan_len)
        #/
        if AC_platforms and ((rangee[0] - AC_platforms[-1][-1] - rangee_existRangee_numBP_refMasked) < rsqELAC_platforms_chain_dist):
            existEntry = AC_platforms[-1]
            existEntry[1].update(rangee[1])
            existEntry[-1] = rangee[-1] #update rcoord
        else:
            AC_platforms.append([rangee[0],rangee[1],rangee[-1]])
    #/
    # Compile into idx-map entry, save data into entry
    for AC_platform in AC_platforms:
        rsqElinks_data = {} #compile data from rsqElink at AC platform, per BP strint. E.g. rBP+rBP_self_LERI. If two strints, then we know both BPs of rsqElink is contained in platform
        for rsqE_idx,rsqElink_idxs in AC_platform[1].items():
            for rsqElink_idx in rsqElink_idxs:            
                # Identify strint of rsqE at rsqElink
                rsqElink = rsqElinks_idx_map[rsqElink_idx]
                strint = None
                if rsqE_idx == rsqElink['rsqE_idx1']:       strint = '1'
                if rsqE_idx == rsqElink['rsqE_idx2']:       strint = '2'
                #/
                # Parse & save data
                init(rsqElinks_data,rsqElink_idx,{})
                rsqElinks_data[rsqElink_idx][strint] = {'strint':strint,'rLERI_self_numToggle':rsqElink['rLERI_self_numToggle'+strint],
                                                        'rBP':rsqElink['rBP'+strint]}

        # save rsqELAC platforms
        tmp_rsqELAC_platform = {'rname':rname,'rcoords':[AC_platform[0],AC_platform[-1]],'rsqE_idxs':rangee[1],
                                'rsqElinks_data':rsqElinks_data}
        assignIdx(tmp_rsqELAC_platform)
        rsqELAC_platforms_idx_map[tmp_rsqELAC_platform['idx']] = tmp_rsqELAC_platform
        #/
        # Mark at rsqElink which rsqELAC platforms it belongs to
        for rsqElink_idx in rsqElinks_data:
            for strint in rsqElinks_data[rsqElink_idx]:
                rsqElink = rsqElinks_idx_map[rsqElink_idx]
                SK = 'rsqELACP_idx'+strint
                rsqElink[SK] = tmp_rsqELAC_platform['idx']
        #/
#/
# Traverse rsqElinks, check those that have "rsqELACP" assignments. Take those rsqELACPs and assign orsqELACP connections
for rsqElink in rsqElinks_idx_map.values():
    if 'rsqELACP_idx1' in rsqElink: #check for any of strints, if have either, should have both - else bug?!
        rsqELACP1 = rsqELAC_platforms_idx_map[ rsqElink['rsqELACP_idx1'] ]
        rsqELACP2 = rsqELAC_platforms_idx_map[ rsqElink['rsqELACP_idx2'] ]
        
        init(rsqELACP1,'orsqELACP_idxs',set())
        init(rsqELACP2,'orsqELACP_idxs',set())
        
        rsqELACP1['orsqELACP_idxs'].add(rsqELACP2['idx'])
        rsqELACP2['orsqELACP_idxs'].add(rsqELACP1['idx'])
#/
# Assign breakpoint characteristics at RrsqELACP: 1/ "normal" fuse, 2/ overlapping regions (LERI point inwards), 3/ back-to-back (LERI point outwards)
for rsqELACP in rsqELAC_platforms_idx_map.values():
    rname = rsqELACP['rname']
    # Preflight: Find "leaps" at rsqELACP ==> rsqElink report of (rBP + rBPselfLERI + orsqELACP.)
    leaps = {}
    for rsqElink_idx,strints_data in rsqELACP['rsqElinks_data'].items():
        if len(strints_data) >= 2: continue #skip if rsqElink has multiple assignments (i.e. BP1 and BP2 of rsqElink exist at rsqELACP)
        for strint,data in strints_data.items():
            rsqElink = rsqElinks_idx_map[rsqElink_idx]
            rsqE1 = qname_rsq_idx_map[rsqElink['rsqE_idx1']]
            rsqE2 = qname_rsq_idx_map[rsqElink['rsqE_idx2']]
            any_rsqE_has_DUP = False
            if 'memberOfDUPs' in rsqE1 or 'memberOfDUPs' in rsqE2:
                any_rsqE_has_DUP = True
            # Parse orsqELACP
            orsqELACP = None
            if rsqELACP['idx'] == rsqElink['rsqELACP_idx1']:        orsqELACP = rsqELAC_platforms_idx_map[ rsqElink['rsqELACP_idx2'] ]
            if rsqELACP['idx'] == rsqElink['rsqELACP_idx2']:        orsqELACP = rsqELAC_platforms_idx_map[ rsqElink['rsqELACP_idx1'] ]
            #/
            # Check cov in ovlp vs. outside rsqELACP
            rsqELACP_scanRange = [min(data['rBP'],data['rBP']+(data['rLERI_self_numToggle']*2000)),max(data['rBP'],data['rBP']+(data['rLERI_self_numToggle']*2000))]
            scanRange_left = [data['rBP']-2000,data['rBP']]
            scanRange_right = [data['rBP'],data['rBP']+2000]
            scanRange_cov = parse_rBGcov(rname,rsqELACP_scanRange) # returns: cov_med,cov_mean,cov_numBPs
            scanRange_left_cov = parse_rBGcov(rname,scanRange_left)
            scanRange_right_cov = parse_rBGcov(rname,scanRange_right)
            #/
            init(leaps,orsqELACP['idx'],[])
            leaps[orsqELACP['idx']].append({'rBP':data['rBP'],'rLERI_self_numToggle':data['rLERI_self_numToggle'],'rsqElink_idx':rsqElink_idx,
                                            'any_rsqE_has_DUP':any_rsqE_has_DUP,'cov_med':scanRange_cov[0],
                                            'covLeft_med':scanRange_left_cov[0],'covRight_med':scanRange_right_cov[0]})
    
    rsqELACP['leaps'] = leaps
    #/
    ## Traverse leaps, assign BP characteristics as header (1+2+3)
    leaps_classis = {} # classi -> orsqELACP-pairs ("comp_id")
    # Check if 0 leaps passed strints thresh. Then means its probably fullyContained AC
    if len(leaps) == 0:
        classi = 'INS'
        leaps_classis[classi] = leaps
    # Check if just x1 orsqELACP
    elif len(leaps) == 1:
        classi = 'single'
        leaps_classis[classi] = leaps
    else:
        # Check 2 + 3. Will later assign 1 (="normal") to every orsqELACP not assigned with 2 or 3.
        orsqELACP_idxs_assigned = {}
        for orsqELACP_idx1,leaps1 in leaps.items():
            for orsqELACP_idx2,leaps2 in leaps.items():
                # Check if multiple leaps from rsqELACP to orsqELACP. This is case for large-scale INV.
                if orsqELACP_idx1 == orsqELACP_idx2:
                    # @@@ Expect to have "back-to-back" with diffDir. @@@
                    # Parse both leap entries
                    leaps_strintSorted_LERIsorted = {}
                    for strint,leaps12 in ( ['1',leaps1],['2',leaps2] ):
                        for leap in leaps12:
                            leap_LERI = leap['rLERI_self_numToggle']
                            leap_rsqElink_idx = leap['rsqElink_idx']
                            init(leaps_strintSorted_LERIsorted,strint,{})
                            init(leaps_strintSorted_LERIsorted[strint],leap_LERI,set())
                            leaps_strintSorted_LERIsorted[strint][leap_LERI].add(leap_rsqElink_idx)

                    IDE = 'IDE: this must be done in another loop. Gonna need to find leap_LERIs in both directions at two platforms. That is the INV.'
                    WIP = 'WIP'
                    continue
                # Else, check 1+2
                else:
                    comp_id = '||'.join(map(str,sorted([orsqELACP_idx1,orsqELACP_idx2])))
                    if comp_id in orsqELACP_idxs_assigned: continue
                
                    # Parse first leap if multiple leaps (if multiple, means haplotype specific rsqElink, but same path. Here we are not interested in haplotypes. so just take any of entries.)
                    leap1 = leaps1[0]
                    leap2 = leaps2[0]
                    #/
                    # Sort by rBP pos
                    leaps_sorted = sorted([leap1,leap2],key=lambda x: x['rBP'])
                    first = leaps_sorted[0]
                    last = leaps_sorted[-1]
                    #/
                    # Parse leaps ovlpRange + dist
                    leaps_ovlpRange = [first['rBP'],last['rBP']]
                    # tweak ovlpRange if they perfectly overlap (fucks up algorithm)
                    if leaps_ovlpRange[0] == leaps_ovlpRange[-1]:       leaps_ovlpRange[-1] += 1
                    orsqELACPs_dist = -getRangeOvlp([leap1['rBP']]*2,[leap2['rBP']]*2)
                    #/
                    # Parse rsqEs of leaps
                    leaps1_rsqEs = []
                    for leap in leaps1:
                        rsqElink = rsqElinks_idx_map[leap['rsqElink_idx']]
                        for strint in ('1','2',):
                            rsqE = qname_rsq_idx_map[rsqElink['rsqE_idx'+strint]]
                            leaps1_rsqEs.append(rsqE)
                    leaps2_rsqEs = []
                    for leap in leaps2:
                        rsqElink = rsqElinks_idx_map[leap['rsqElink_idx']]
                        for strint in ('1','2',):
                            rsqE = qname_rsq_idx_map[rsqElink['rsqE_idx'+strint]]
                            leaps2_rsqEs.append(rsqE)
                    #/
                    ## To get ACs at rsqELACP (we need this for TRAblock classi), we must traverse all reads at each leaps rsqElink
                    leaps_ACs = {'1':[],'2':[]}
                    for strenum,tmp_leaps in ( ['1',leaps1],['2',leaps2] ):
                        for tmp_leap in tmp_leaps:
                            rsqElink = rsqElinks_idx_map[tmp_leap['rsqElink_idx']]
                            for read in rsqElink['reads']():
                                if not read in qname_ACs: continue #skip if read has no ACs
                                for AC_idx in qname_ACs[read]:
                                    AC = AC_idx_map[AC_idx]
                                    if AC['rname'] == rsqELACP['rname'] and getRangeOvlp(AC['rcoords'],rsqELACP['rcoords']) >= 0:
                                        leaps_ACs[strenum].append(AC)
                    ##/
                    # Run checks
                    classi = 'ovlp_sameDir' # "normal" but this implies that BPs ovlp at rsqELACP and could be haplotype diffs
                    if first['rLERI_self_numToggle'] == 1 and last['rLERI_self_numToggle'] == -1: # Check 2
                        # Check if both have extension into rsqELACP outside leaps_ovlpRange (if not, then its more an INS-TRA-BLOCK)
                        #if (rsqELACP['rcoords'][0] - leaps_ovlpRange[0]) > -500 and (rsqELACP['rcoords'][-1] - leaps_ovlpRange[-1]) < 500:
                        # Update: check if either have AC which is 100% in leaps_ovlpRange. Then its TRA-block.
                        leaps_AC_notContained = False
                        for strenum,ACs in leaps_ACs.items():
                            for AC in ACs:
                                # check if AC exist outside ovlp
                                if (AC['rcoords'][-1]-AC['rcoords'][0]) > getRangeOvlp(AC['rcoords'],leaps_ovlpRange):
                                    leaps_AC_notContained = True
                                    break #break on first
                        #/
                        # make classi
                        if not leaps_AC_notContained:
                            classi = 'TRAblock'
                        else:
                            classi = 'ovlp_diffDir'
                        #/
                    if first['rLERI_self_numToggle'] == -1 and last['rLERI_self_numToggle'] == 1: # Check 3
                        classi = 'backToBack'
                    #/
                    
                    # Check cov in ovlp vs. outside range
                    scanRange_left = [leaps_ovlpRange[0]-2000,leaps_ovlpRange[0]]
                    scanRange_right = [leaps_ovlpRange[-1],leaps_ovlpRange[-1]+2000]
                    leaps_ovlpRange_cov = parse_rBGcov(rname,leaps_ovlpRange) # returns: cov_med,cov_mean,cov_numBPs
                    scanRange_left_cov = parse_rBGcov(rname,scanRange_left)
                    scanRange_right_cov = parse_rBGcov(rname,scanRange_right)
                    #/
                    # Save
                    init(leaps_classis,classi,{})
                    leaps_classis[classi][comp_id] = {'rBPs_dist':orsqELACPs_dist,'r_first':first,'r_last':last,
                                                      'leaps1_rsqEs':HOLDER(leaps1_rsqEs),'leaps2_rsqEs':HOLDER(leaps2_rsqEs),
                                                      'leaps1_ACs':HOLDER(leaps_ACs['1']),'leaps2_ACs':HOLDER(leaps_ACs['2']),
                                                      'leaps_range':leaps_ovlpRange,'leaps_range_cov_med':leaps_ovlpRange_cov[0],
                                                      'covLeft_med':scanRange_left_cov[0],'covRight_med':scanRange_right_cov[0]}
                    orsqELACP_idxs_assigned[comp_id] = classi
                    #/
    #/
    # Save at rsqELACP
    rsqELACP['leaps_classis'] = leaps_classis
    #/
    ##/
#/

### INFO: plot number of orsqELACPs per orsqELACPs
if 0:
    rsqELACP_stats = {'num_orsqELACPs':[],'LERI_states':[]}
    for rsqELACP in rsqELAC_platforms_idx_map.values():
        num_orsqELACPs = len(rsqELACP['orsqELACP_idxs'])
        rsqELACP_stats['num_orsqELACPs'].append(num_orsqELACPs)        
        
        if num_orsqELACPs >= 2:
            rBP_LERIs = {}
            for rsqElink_idx,strints_data in rsqELACP['rsqElinks_data'].items():
                if len(strints_data) >= 2: continue #skip if rsqElink has multiple assignments (i.e. BP1 and BP2 of rsqElink exist at rsqELACP)
                for strint,data in strints_data.items():
                    rLERI = data['rLERI_self_numToggle']
                    init(rBP_LERIs,rLERI,0)
                    rBP_LERIs[rLERI] += 1
            if len(rBP_LERIs) >= 2:
                rsqELACP_stats['LERI_states'].append(1)
                print(num_orsqELACPs,rsqELACP['rname'],rsqELACP['rcoords'])
                print()
            else:
                rsqELACP_stats['LERI_states'].append(0)

    filt0 = rsqELACP_stats['num_orsqELACPs']
    filt1 = filterArr( filt0, 0, 'threshold' )
    filt2 = filterArr( filt1, 500, 'cap' )
    fig = plt.figure(figsize=(11,11))
    plt.hist( filt2 ,bins=300)
    plt.title('rsqELACP_stats'+' numVals='+str(len(filt2)))
    
    # rsqELACP classis
    rsqELACP_leapClassis_abu = {}
    rsqELACP_leapClassiss_abu = {}
    rsqELACP_leapClassis_lens = []
    rsqELACP_leapClassis_juncLens = {}
    rsqELACP_leapClassis_readSupps_scatter = {}
    rsqELACP_leapClassis_readSupps_ratio = {}
    rsqELACP_leapClassis_hasINSinLink = []
    rsqELACP_leapClassis_hasDUPs = {}
    rsqELACP_leapClassis_hasHigherCov = {}
    rsqELACP_leapClassis_hasHigherCov_vals = []
    rsqELACP_diffDirs = []
    rsqELACP_TRAblocks = []
    rsqELACP_diffDirs_bothPlatforms = set()
    rsqELACP_simple_num_rsqElinks = []
    rsqELACP_leapClassis_rdists = {}
    rsqELACP_leapClassis_qdists = {}
    rsqELACP_leapClassis_dists_compsDone = set()
    rsqELACP_comps_scatterSuppPlot = {} # rsqElink_idx -> plot data
    rsqELACP_comps_scatterSuppPlot2 = {} # rsqELACP_comp_id -> rsqELACP_idx -> rsqElink_idx -> data
    rsqELACP_scatterSuppPlot = {} # rsqELACP_idx -> rsqElink_idx -> data
    rsqELACP_comps_classified_wCov = {}
    rsqELACP_PWcomps_wCov = {}
    IDE_extractReads = set()
    
    capCovAt = 600
    for rsqELACP_idx,rsqELACP in rsqELAC_platforms_idx_map.items():
        # Parse rDist + qDist for leapClassis
        for classi in rsqELACP['leaps_classis']:
            for orsqELACP_idxs,data in rsqELACP['leaps_classis'][classi].items():
                # check if its a "simple" BND, then reformat for loop convenience
                if type(orsqELACP_idxs) == int:
                    orsqELACP_idxs = [orsqELACP_idxs]
                else:
                    orsqELACP_idxs = list(map(int,orsqELACP_idxs.split('||')))
                #/
                for orsqELACP_idx in orsqELACP_idxs:
                    # check if done this comp already
                    comp_id = '||'.join(sorted(map(str,[rsqELACP_idx,orsqELACP_idx])))
                    if comp_id in rsqELACP_leapClassis_dists_compsDone: continue
                    rsqELACP_leapClassis_dists_compsDone.add(comp_id)
                    #/
                    rsqElink = rsqElinks_idx_map[ rsqELACP['leaps'][orsqELACP_idx][0]['rsqElink_idx'] ]
                    # parse rDist
                    rname = 'TRA'
                    dist_r = -1
                    if rsqElink['rname1'] == rsqElink['rname2']:
                        rname = rsqElink['rname1']
                        dist_r = abs(rsqElink['rBP1'] - rsqElink['rBP2'])
                    init(rsqELACP_leapClassis_rdists,classi,{})
                    init(rsqELACP_leapClassis_rdists[classi],rname,[])
                    rsqELACP_leapClassis_rdists[classi][rname].append(dist_r)
                    #/
                    # Parse qDist
                    dist_q = abs(rsqElink['qcoord1'] - rsqElink['qcoord2'])
                    init(rsqELACP_leapClassis_qdists,classi,{})
                    init(rsqELACP_leapClassis_qdists[classi],rname,[])
                    rsqELACP_leapClassis_qdists[classi][rname].append(dist_q)
                    #/
                #/
        #/
        
        # Parse connections between rsqELACPs for scatterplot (single-mode)
        for orsqELACP_idx,leapDatas in rsqELACP['leaps'].items():
            for leapData in leapDatas:
                rsqElink = rsqElinks_idx_map[ leapData['rsqElink_idx'] ]
                rsqElink_supp_min = min( [ len(rsqElink['rAVA_supp1_margin=3000']()) , len(rsqElink['rAVA_supp2_margin=3000']()) ] )
                cov_med = leapData['cov_med']
                
                init(rsqELACP_PWcomps_wCov,rsqELACP_idx,{})
                init(rsqELACP_PWcomps_wCov[rsqELACP_idx],orsqELACP_idx,[])
                rsqELACP_PWcomps_wCov[rsqELACP_idx][orsqELACP_idx].append( [rsqElink['idx'],rsqElink_supp_min,cov_med] )
        #/
        
        classis = []
        for classi in rsqELACP['leaps_classis']:
            classis.append(classi)
            init(rsqELACP_leapClassis_abu,classi,0)
            rsqELACP_leapClassis_abu[classi] += 1
            if classi in ('single','INS',):
                for orsqELACP_idx in rsqELACP['leaps_classis'][classi]:
                    for rsqElink_entry in rsqELACP['leaps_classis'][classi][orsqELACP_idx]:
                        if rsqElink_entry['any_rsqE_has_DUP']:
                            init(rsqELACP_leapClassis_hasDUPs,classi,0)
                            rsqELACP_leapClassis_hasDUPs[classi] += 1
                            
                if classi == 'single':
                    rsqElinks = set()
                    for orsqELACP_idx in rsqELACP['leaps_classis'][classi]:
                        for rsqElink_entry in rsqELACP['leaps_classis'][classi][orsqELACP_idx]:
                            rsqElinks.add(rsqElink_entry['rsqElink_idx'])
                    rsqELACP_simple_num_rsqElinks.append(len(rsqElinks))
                    if 0 and rsqELACP_simple_num_rsqElinks[-1] == 2:
                        for rsqElink_idx in rsqElinks:
                            print(len(rsqElinks_idx_map[rsqElink_idx]['reads']()))
                        print()
                        
                    ## rAVA supp vs rsqELACP cov
                    for orsqELACP_idx,data_entries in rsqELACP['leaps_classis'][classi].items():
                        for data in data_entries:
                            #data = data_entries[0] # parse data from first rsqElink
                            rsqElink = rsqElinks_idx_map[ data['rsqElink_idx'] ]
                            
                            rsqElink_supp_min = min( [ len(rsqElink['rAVA_supp1_margin=3000']()) , len(rsqElink['rAVA_supp2_margin=3000']()) ] )
                                
                            rsqElink_platform_cov = data['cov_med']
                            
                            # save for detailed scatter plot
                            init(rsqELACP_scatterSuppPlot,rsqELACP_idx,{})
                            rsqELACP_scatterSuppPlot[rsqELACP_idx][rsqElink['idx']] = [rsqElink_supp_min,rsqElink_platform_cov,classi]
                            
                            rsqELACPs_comp_id = '||'.join(map(str,sorted([rsqELACP_idx,orsqELACP_idx])))
                            init(rsqELACP_comps_classified_wCov,rsqELACPs_comp_id,{})
                            init(rsqELACP_comps_classified_wCov[rsqELACPs_comp_id],classi,[])
                            rsqELACP_comps_classified_wCov[rsqELACPs_comp_id][classi].append( [rsqElink['idx'],rsqElink_supp_min,rsqElink_platform_cov] )
                            #/
                            
                            init(rsqELACP_leapClassis_readSupps_ratio,classi,[])
                            #rsqElink_leapsOvlp_cov = 1
                            rsqELACP_leapClassis_readSupps_ratio[classi].append(rsqElink_supp_min/rsqElink_platform_cov)
                            # make supp on negative axis if high cov
                            if rsqElink_platform_cov >= capCovAt:   rsqELACP_leapClassis_readSupps_ratio[classi][-1] = -rsqELACP_leapClassis_readSupps_ratio[classi][-1]
                            
                            if rsqElink_platform_cov >= capCovAt:   rsqElink_platform_cov = -1
                            if rsqElink_supp_min >= capCovAt:        rsqElink_supp_min = -1
                            init(rsqELACP_leapClassis_readSupps_scatter,classi,{'x':[],'y':[]})
                            rsqELACP_leapClassis_readSupps_scatter[classi]['x'].append(rsqElink_supp_min)
                            rsqELACP_leapClassis_readSupps_scatter[classi]['y'].append(rsqElink_platform_cov)
                    ##/
            else:
                rBPs_dist = list(rsqELACP['leaps_classis'][classi].values())[0]['rBPs_dist']
                rsqELACP_leapClassis_lens.append(rBPs_dist)
                
                # Check if has DUPclust
                rsqEs_w_DUP = []
                for CK in ('leaps1_rsqEs','leaps2_rsqEs',):
                    for comp_id in rsqELACP['leaps_classis'][classi]:
                        for rsqE in rsqELACP['leaps_classis'][classi][comp_id][CK]():
                            if 'memberOfDUPs' in rsqE:
                                rsqEs_w_DUP.append(rsqE)
                if rsqEs_w_DUP:
                    init(rsqELACP_leapClassis_hasDUPs,classi,0)
                    rsqELACP_leapClassis_hasDUPs[classi] += 1
                else:
                    if classi == 'ovlp_diffDir':
                        pass
                #/
                # Check if leap ovlpRange has higher cov than left/right regions of it
                leaps_range_has_higher_cov = False
                covDiffs = []
                rsqElink_qlens = []
                for comp_id in rsqELACP['leaps_classis'][classi]:
                    comp_entry = rsqELACP['leaps_classis'][classi][comp_id]
                    if comp_entry['leaps_range_cov_med'] > comp_entry['covLeft_med'] and comp_entry['leaps_range_cov_med'] > comp_entry['covRight_med']:
                        cov_diff = comp_entry['leaps_range_cov_med'] - max(comp_entry['covLeft_med'],comp_entry['covRight_med'])
                        if cov_diff >= 10:
                            leaps_range_has_higher_cov = True
                            covDiffs.append(cov_diff)
                            
                    rsqElink1 = rsqElinks_idx_map[comp_entry['r_first']['rsqElink_idx']]
                    rsqElink2 = rsqElinks_idx_map[comp_entry['r_last']['rsqElink_idx']]
                    len_q1 = rsqElink1['len_q']
                    len_q2 = rsqElink2['len_q']
                    rsqElink_qlens.append(len_q1)
                    rsqElink_qlens.append(len_q2)
                
                init(rsqELACP_leapClassis_juncLens,classi,[])
                rsqELACP_leapClassis_juncLens[classi].append(max(rsqElink_qlens))

                if leaps_range_has_higher_cov:
                    init(rsqELACP_leapClassis_hasHigherCov,classi,0)
                    rsqELACP_leapClassis_hasHigherCov[classi] += 1
                    rsqELACP_leapClassis_hasHigherCov_vals.append(covDiffs)
                
                # Extract rsqELACPs for classis
                if classi == 'TRAblock':
                    rsqELACP_TRAblocks.append(rsqELACP)
                if classi == 'ovlp_diffDir':
                    rsqELACP_diffDirs.append(rsqELACP)
                #/
                # Extract read supp for classis
                for comp_id in rsqELACP['leaps_classis'][classi]:
                    comp_entry = rsqELACP['leaps_classis'][classi][comp_id]
                    rsqElink_rFirst = rsqElinks_idx_map[ comp_entry['r_first']['rsqElink_idx'] ]
                    rsqElink_rLast = rsqElinks_idx_map[ comp_entry['r_last']['rsqElink_idx'] ]
                    for rsqElink in ( rsqElink_rFirst,rsqElink_rLast ):
                        rsqElink_supp_min = min( [ len(rsqElink['rAVA_supp1_margin=3000']()) , len(rsqElink['rAVA_supp2_margin=3000']()) ] )
                        
                        rsqElink_leapsOvlp_cov = comp_entry['leaps_range_cov_med']
                        # correct cov if back-to-back
                        if classi.find('back') != -1:
                            if rsqElink['idx'] == rsqElink_rFirst['idx']:       rsqElink_leapsOvlp_cov = comp_entry['covLeft_med']
                            if rsqElink['idx'] == rsqElink_rLast['idx']:        rsqElink_leapsOvlp_cov = comp_entry['covRight_med']
                        
                        ## save for detailed scatter plot
                        init(rsqELACP_scatterSuppPlot,rsqELACP_idx,{})
                        rsqELACP_scatterSuppPlot[rsqELACP_idx][rsqElink['idx']] = [rsqElink_supp_min,rsqElink_leapsOvlp_cov,classi]
                        
                        # Parse rsqElink for rsqELACP <-> orsqELACP comp
                        for orsqELACP_idx in list(map(int,comp_id.split('||'))):
                            orsqELACP = rsqELAC_platforms_idx_map[orsqELACP_idx]
                            leaps_data = rsqELACP['leaps'][orsqELACP_idx]
                            for leap_data in leaps_data:
                                if rsqElink['idx'] == leap_data['rsqElink_idx']:
                                    rsqELACPs_comp_id = '||'.join(map(str,sorted([rsqELACP_idx,orsqELACP_idx])))
                                    init(rsqELACP_comps_classified_wCov,rsqELACPs_comp_id,{})
                                    init(rsqELACP_comps_classified_wCov[rsqELACPs_comp_id],classi,[])
                                    rsqELACP_comps_classified_wCov[rsqELACPs_comp_id][classi].append( [rsqElink['idx'],rsqElink_supp_min,rsqElink_leapsOvlp_cov] )
                        #/
                        ##/
                        
                        init(rsqELACP_leapClassis_readSupps_ratio,classi,[])
                        #rsqElink_leapsOvlp_cov = 1
                        rsqELACP_leapClassis_readSupps_ratio[classi].append(rsqElink_supp_min/rsqElink_leapsOvlp_cov)
                        # make supp on negative axis if high cov
                        if rsqElink_leapsOvlp_cov >= capCovAt:   rsqELACP_leapClassis_readSupps_ratio[classi][-1] = -rsqELACP_leapClassis_readSupps_ratio[classi][-1]
                        
                        if rsqElink_leapsOvlp_cov >= capCovAt:   rsqElink_leapsOvlp_cov = -1
                        if rsqElink_supp_min >= capCovAt:        rsqElink_supp_min = -1
                        init(rsqELACP_leapClassis_readSupps_scatter,classi,{'x':[],'y':[]})
                        rsqELACP_leapClassis_readSupps_scatter[classi]['x'].append(rsqElink_supp_min)
                        rsqELACP_leapClassis_readSupps_scatter[classi]['y'].append(rsqElink_leapsOvlp_cov)
                #/
        # Save all classis
        classii = '||'.join(sorted(classis))
        init(rsqELACP_leapClassiss_abu,classii,0)
        rsqELACP_leapClassiss_abu[classii] += 1
        #/
        # Check diffDir at both platforms (e.g. 3R12.9<->3R21)
        if 'ovlp_diffDir' in rsqELACP['leaps_classis']:
            for orsqELACP_idx in rsqELACP['leaps']:
                orsqELACP = rsqELAC_platforms_idx_map[orsqELACP_idx]
                if 'ovlp_diffDir' in orsqELACP['leaps_classis']:
                    comp_id = '||'.join(sorted(map(str,[rsqELACP_idx,orsqELACP_idx])))
                    
                    # IDE
                    if not comp_id in rsqELACP_diffDirs_bothPlatforms:
                        print(rsqELACP['rname'],rsqELACP['rcoords'],orsqELACP['rname'],orsqELACP['rcoords'])
                    #/
                    
                    rsqELACP_diffDirs_bothPlatforms.add(comp_id)
                    
        #/
        #if 'TRAblock' in rsqELACP['leaps_classis']: print(rsqELACP['rname'],rsqELACP['rcoords'])


if 0 and 'plot?':
    fig = plt.figure(figsize=(11,11))
    for key,vals in rsqELACP_leapClassis_juncLens.items():
        #fig = plt.figure(figsize=(11,11))
        filt0 = vals
        filt1 = filterArr( filt0, 0, 'threshold' )
        filt2 = filterArr( filt1, 1500, 'cutoff' )
        if not filt2: continue
        plt.hist( filt2 ,label=key,bins=300)
        #plt.title(key+' || numVals='+str(len(filt2)))
    plt.title('rsqELACP_leapClassis_juncLens')
    plt.legend()
    
    ## Scatter of classis
    for key,xy in rsqELACP_leapClassis_readSupps_scatter.items():
        fig = plt.figure(figsize=(11,11))
        plt.scatter(xy['x'],xy['y'])
        plt.plot([0,600],[genomeAverageCov_A,genomeAverageCov_A],linestyle='dashed')
        plt.plot([40,40],[0,600],linestyle='dashed')
        plt.plot([80,80],[0,600],linestyle='dashed')
        plt.title(key+' || numVals='+str(len(xy['x'])))
        plt.xlim([-5,max(xy['x'])+5])
        plt.ylim([-5,max(xy['y'])+5])
    ##/
    ## Scatter of comp_ids, with markers for e.g. MITO
    fig = plt.figure(figsize=(11,11))
    
    # Points
    plotLines_comps = {}
    rsqELACP_compID_pointers = {} #for marking those with platforms to multiple breakpoints
    for rsqELACP_idx,orsqELACP_idxs_data in rsqELACP_PWcomps_wCov.items():
        for orsqELACP_idx,datas in orsqELACP_idxs_data.items():
            comp_id = '||'.join(map(str,sorted([rsqELACP_idx,orsqELACP_idx])))
            # select top supp between rsqELACPs
            top_entry = sorted(datas,key=lambda x: x[1],reverse=True)[0]
            #/
            rsqELACP = rsqELAC_platforms_idx_map[rsqELACP_idx]
            orsqELACP = rsqELAC_platforms_idx_map[orsqELACP_idx]
            rsqElink_idx,rsqElink_supp_min,rsqElink_platform_cov = top_entry
            rsqElink = rsqElinks_idx_map[rsqElink_idx]
            
            plotLineBetweenPoints = False
            #if rsqElink['rname1'] == 'mitochondrion_genome' and rsqElink['rname2'] == 'mitochondrion_genome': continue
            if 0 and 'mitochondrion_genome' in ( rsqElink['rname1'],rsqElink['rname2'] ):
                pointCol = 'red'
                plt.text(rsqElink_supp_min,rsqElink_platform_cov,rsqELACP['rname'],fontsize=9)
                plotLineBetweenPoints = True
                
            elif rsqELACP['rname'] != orsqELACP['rname']:
                pointCol = 'black'
                rcoords_mil = [round(rsqELACP['rcoords'][0],-5)/10**6,round(rsqELACP['rcoords'][-1],-5)/10**6]
                plt.text(rsqElink_supp_min,rsqElink_platform_cov,rsqELACP['rname']+':'+'-'.join(map(str,rcoords_mil)),fontsize=9)
                plotLineBetweenPoints = True
            
            if plotLineBetweenPoints:
                plt.plot([rsqElink_supp_min],[rsqElink_platform_cov],marker='o',markersize=4,color=pointCol) # x,y
                
                init(plotLines_comps,comp_id,[])
                plotLines_comps[comp_id].append([rsqElink_supp_min,rsqElink_platform_cov])
                
                init(rsqELACP_compID_pointers,rsqELACP_idx,[])
                rsqELACP_compID_pointers[rsqELACP_idx].append([rsqElink_supp_min,rsqElink_platform_cov])
    #/
    # Lines
    for comp_id,lineData in plotLines_comps.items():
        rsqElink_supp_min1,rsqElink_platform_cov1 = lineData[0]
        rsqElink_supp_min2,rsqElink_platform_cov2 = lineData[1]
        plt.plot([rsqElink_supp_min1,rsqElink_supp_min2],[rsqElink_platform_cov1,rsqElink_platform_cov2],linestyle='solid',color='black',linewidth=0.5)    
        #/
    #/
    # Lines, multi-BP-platforms
    for rsqELACP_idx,lineDatas in rsqELACP_compID_pointers.items():
        for enum,_ in enumerate(lineDatas):
            if enum == 0: continue
            rsqElink_supp_min1,rsqElink_platform_cov1 = lineDatas[enum-1]
            rsqElink_supp_min2,rsqElink_platform_cov2 = lineDatas[enum]
            plt.plot([rsqElink_supp_min1,rsqElink_supp_min2],[rsqElink_platform_cov1,rsqElink_platform_cov2],linestyle='dashed',color='green',linewidth=0.5)    
        #/
    #/
    
    ax = plt.gca()
    xlim = ax.get_xlim()[1]
    ylim = ax.get_ylim()[1]
    plt.plot([0,xlim],[genomeAverageCov_A*2,genomeAverageCov_A*0.5],linestyle='dashed',color='dodgerblue',linewidth=1)
    plt.plot([0,xlim],[genomeAverageCov_A,genomeAverageCov_A],linestyle='dashed',color='dodgerblue',linewidth=1)
    plt.plot([0,xlim],[genomeAverageCov_A*2,genomeAverageCov_A*2],linestyle='dashed',color='dodgerblue',linewidth=1)
    plt.plot([40,40],[0,ylim],linestyle='dashed',color='dodgerblue',linewidth=1)
    plt.plot([80,80],[0,ylim],linestyle='dashed',color='dodgerblue',linewidth=1)
    plt.plot([120,120],[0,ylim],linestyle='dashed',color='dodgerblue',linewidth=1)
    plt.plot([160,160],[0,ylim],linestyle='dashed',color='dodgerblue',linewidth=1)
    plt.yscale('log')
    plt.xscale('log')
    #plt.xlim([10,600])
    #plt.ylim([50,1000])
    ##/
    
    for key,vals in rsqELACP_leapClassis_readSupps_ratio.items():
        filt0 = rsqELACP_simple_num_rsqElinks
        filt1 = filterArr( filt0, -100, 'threshold' )
        filt2 = filterArr( filt1, 500, 'cap' )
        fig = plt.figure(figsize=(11,11))
        plt.hist( filt2 ,bins=300)
        plt.title(key+' || numVals='+str(len(filt2)))
    
    ## Violin-plot rDists of classis
    mergedClassis_violins = {} #plotMerged, below
    for classi,rname_rdists in rsqELACP_leapClassis_rdists.items():
        xlabels = []
        violins = []
        for rname,data in sorted(rname_rdists.items(),key=lambda x: x[0]):
            if not rname in chroms: continue
            filt0 = data
            filt1 = filterArr( filt0, 1, 'raise' )
            filt2 = filterArr( filt1, 99000000, 'cap' )
            if not filt2: continue
            if 1 and 'log10?':      filt2 = [math.log(x,10) for x in filt2]
            violins.append(filt2)
            xlabels.append(rname+'\nN='+str(len(filt2)))
            #plotMerged
            init(mergedClassis_violins,rname,[])
            mergedClassis_violins[rname] += filt2
            #/
        fig = plt.figure(figsize=(11,11))
        plt.violinplot(violins,showmedians=True)
        ax = plt.gca()
        ax.set_xticks(np.arange(len(xlabels))+1)
        ax.set_xticklabels(xlabels,rotation=0)
        plt.title(classi+' || numVals='+str(sum(len(x) for x in violins)))
    
    # ALL IN SAME (+hist)
    xlabels = []
    violins = []
    for rname,data in sorted(mergedClassis_violins.items(),key=lambda x: x[0]):
        violins.append(data)
        xlabels.append(rname+'\nN='+str(len(data)))
    fig = plt.figure(figsize=(11,11))
    plt.violinplot(violins,showmedians=True)
    ax = plt.gca()
    ax.set_xticks(np.arange(len(xlabels))+1)
    ax.set_xticklabels(xlabels,rotation=0)
    plt.title('rsqELACPs_all'+' || numVals='+str(sum(len(x) for x in violins)))
    
    vals = []
    for classi,rname_rdists in rsqELACP_leapClassis_qdists.items():
        for rname,dists in rname_rdists.items():
            vals += dists
    filt0 = vals
    filt1 = filterArr( filt0, 0, 'threshold' )
    filt2 = filterArr( filt1, 99*20000, 'cap' )
    fig = plt.figure(figsize=(11,11))
    plt.hist( filt2 ,bins=50)
    plt.title('rsqELACP qdists || numVals='+str(len(filt2)))
    #/
    ##/
###/
###/

### INFO calc haplo supps for SIMPLE rsqELACPs + check ovlps with DCros
rsqELACP_DCro_ovlps = {'yes':0,'no':0}
reads_3LmitoX = set()
for rsqELACP in rsqELAC_platforms_idx_map.values():    
    
    # IDE fetch reads
    for rsqElink_idx in rsqELACP['rsqElinks_data']:
        rsqElink = rsqElinks_idx_map[rsqElink_idx]
        ovlps_selection = []
        for strint in ('1','2',):
            if rsqElink['rname'+strint] == '3L' and getRangeOvlp(rsqElink['rBPcov_range'+strint],[3352186,3486147]) >= 0:
                ovlps_selection.append(strint)
            #if rsqElink['rname'+strint] == 'X' and getRangeOvlp(rsqElink['rBPcov_range'+strint],[14908014,14932317]) >= 0:
            #    ovlps_selection.append(strint)
            if rsqElink['rname'+strint] == 'mitochondrion_genome' and getRangeOvlp(rsqElink['rBPcov_range'+strint],[0,18000]) >= 0:
                ovlps_selection.append(strint)
        if len(ovlps_selection) >= 2:
            reads_3LmitoX.add(rsqElink['qname'])
    #/
        
    for classi in rsqELACP['leaps_classis']:
        if classi != 'single': continue
        ## DCro ovlps
        # Parse rBP + rBP_selfLERI from any rsqElink
        rBP = None
        rBP_selfLERI = None
        for orsqELACP_idx,orsqELACP_data in rsqELACP['leaps_classis'][classi].items():
            for rsqElink_entry in orsqELACP_data:
                rsqElink_idx = rsqElink_entry['rsqElink_idx']
                
                for strint,rsqElink_rsqELACPdata in rsqELACP['rsqElinks_data'][rsqElink_idx].items():
                    rBP = rsqElink_rsqELACPdata['rBP']
                    rBP_selfLERI = rsqElink_rsqELACPdata['rLERI_self_numToggle']
                    break
                if rBP != None: break
            if rBP != None: break
        #/
        # Check ovlp
        rcoords_PRE = [rBP,rBP+(1000*rBP_selfLERI)]
        rcoords_to_scan = [ min(rcoords_PRE),max(rcoords_PRE) ]
        DCro_ovlps = rangeOverlaps_lookup(rcoords_to_scan,DCro_idx_map_ROMD,1000,map_lookup_key=rsqELACP['rname'],skipNumBPOvlpBelow=500)
        if DCro_ovlps:
            rsqELACP_DCro_ovlps['yes'] += 1
        else:
            rsqELACP_DCro_ovlps['no'] += 1
        #/
        ##/
###/
"""
### VGSEGS (variant graph genome segments)
readSupp_req = 3
## make rRanges of rBPs, merge + make gsegs
# init rBPs
rname_rBPs = {}
for rsqElink in rsqElinks_idx_map.values():
    ## Check if filter
    if not (rsqElink['rname1'] in chroms+['mitochondrion_genome'] and rsqElink['rname1'] in chroms+['mitochondrion_genome']): continue
    if rsqElink['rname1'] == rsqElink['rname2'] and not rsqElink['len_r'] >= 100000: continue
    if rsqElink['isSwitchOfContainedAC'] and not (rsqElink['rname1'] == 'mitochondrion_genome' or rsqElink['rname2'] == 'mitochondrion_genome'): continue #skip if a "big-ins"
    rsqE1,rsqE2 = qname_rsq_idx_map[rsqElink['rsqE_idx1']],qname_rsq_idx_map[rsqElink['rsqE_idx2']]
    if 'memberOfDUPs' in rsqE1 and 'memberOfDUPs' in rsqE2 and rsqE1['rname'] == rsqE2['rname']: continue #skip if local duplication
    #/
    if not rsqElink['isUsed']: continue

    # Check read supp for rsqElink
    rAVA_supps = []
    for strint in ('1','2',):
        rAVA_supps.append( len(rsqElink['rAVA_supp'+strint+'_margin=3000']()) >= readSupp_req )
    if not all(rAVA_supps): continue
    #/
    
    for strint in ('1','2',):
        rname = rsqElink['rname'+strint]
        rBP = rsqElink['rBP'+strint]
        init(rname_rBPs,rname,[])
        rname_rBPs[rname].append([rBP,rsqElink['rLERI_self_numToggle'+strint],strint,rsqElink['idx']])
#/
# Cluster
def init_vgSeg(): # Run inside funct for easy "last-add". Remove DEF-line and dedent to have vanilla back
    # Check if we had a previous rBP clust, then sort that clust selfLERIs:
    # *selfLERIs with 1 should go to prev gseg tail
    # *selfLERIs with -1 should go to cur gseg head
    rsqElinks_selfLERIplus = []
    rsqElinks_selfLERIneg = []
    if prev_rBPclust:
        for rsqElink_entry in prev_rBPclust['source_rsqElink_entries']:
            if rsqElink_entry[1] == 1:      rsqElinks_selfLERIplus.append(rsqElink_entry)
            if rsqElink_entry[1] == -1:     rsqElinks_selfLERIneg.append(rsqElink_entry)
    #/
    # If prev gseg, add negative self LERIs to that one (if it doesnt exist, then add it to current. Probably its a bullshit at chrom end, but w/e)
    if rGsegs:
        rGsegs[-1]['tail_rsqElink_entries'] = rsqElinks_selfLERIneg
    else:
        rsqElinks_selfLERIplus += rsqElinks_selfLERIneg
        print('[Warning] This is not expected to happens more than once per chromosome',rname)
    #/
    
    # init gseg
    maskClassi,maskCovFrac = parse_rspanContent(rname,[rangee[0],rangee[-1]],return_maskedCovFrac=True)
    
    tmp_gseg = {'rname':rname,'rcoords':[rangee[0],rangee[-1]],'len':rangee[-1]-rangee[0],
                'head_rsqElink_entries':rsqElinks_selfLERIplus,'tail_rsqElink_entries':[],
                'head_refGseg':None,'tail_refGseg':None,'masked_covFrac':maskCovFrac,'maskClassi':HOLDER(maskClassi)}
    assignIdx(tmp_gseg)
    #/
    # Get+Set ref-link info: From "prev gseg to cur gseg" + "from cur gseg to prev gseg"
    if rGsegs:
        tmp_gseg['head_refGseg'] = rGsegs[-1]['idx']
        rGsegs[-1]['tail_refGseg'] = tmp_gseg['idx']
    #/
    # Save
    rGsegs.append(tmp_gseg)
    #/
    
chain_dist = 1000 # chain rBPs within this dist into same vgseg
IDE_clustLens = []
vgseg_idx_map = {}
for rname,rBPs in rname_rBPs.items():
    # Chain clusters
    rBP_clusters = []
    for rBP_entry in sorted(rBPs,key=lambda x: x[0]):
        rBP = rBP_entry[0]
        # check if extend existing cluster
        if rBP_clusters and (rBP - rBP_clusters[-1]['rBPs'][-1]) < chain_dist:
                rBP_clusters[-1]['rBPs'].append(rBP)
                rBP_clusters[-1]['source_rsqElink_entries'].append(rBP_entry)
        # else init new cluster
        else:
            rBP_clusters.append( {'rBPs':[rBP],'source_rsqElink_entries':[rBP_entry] } )
    #/
    # Make gsegs
    rRanges = []
    rRanges.append([0,'base',refLens[rname]])
    for enum,clust in enumerate(rBP_clusters):
        rstart,rend = min(clust['rBPs']),max(clust['rBPs'])
        if rend == rstart:  rend += 1
        rRanges.append([rstart,enum,rend])
    ovlps = computeOvlpRanges_wIdxs(rRanges)
   
    rGsegs = []
    prev_rBPclust = None
    for rangeeEnum,rangee in enumerate(ovlps):
        # check if only base, then we want to init gseg
        if len(rangee[1]) == 1 and 'base' in rangee[1]:
            init_vgSeg()
        #/
        # else, update tracker variable "prev_rBPclust"
        else:
            if len(rangee[1]) >= 3: sys.exit('didnt expect this?!?!?!')
            for entry in rangee[1]:
                if not entry == 'base':
                    prev_rBPclust = rBP_clusters[entry] #entry = enum of cluster
        #/
    #@ TO DO: Handle last-add, if it wasnt a 'base', then we will have missed an entry.
    if not len(ovlps[-1][1]) == 1:
        init_vgSeg()
        #print('IDE:occ, missed last. Implement later if strategy is viable',rname)
    #@
    #/
    # Save to outer
    for enum,gseg in enumerate(rGsegs):
        vgseg_idx_map[gseg['idx']] = gseg
    #/
vgseg_idx_map_ROMD = rangeOverlaps_makeMappingDict(vgseg_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
#/
# add cov to vgseg
for vgseg in vgseg_idx_map.values():
    rname=vgseg['rname']
    rstart,rend=vgseg['rcoords']
    cov_med,cov_mean = 1,1
    #cov_med,cov_mean,cov_data = parse_rBGcov(rname,[rstart,rend],skipMaskedCovFrac=0.01)
    vgseg['cov'] = cov_med
    vgseg['cov_mean'] = cov_mean
#/

# REFSUPP: Check which refGseg connections exist in S2
for vgseg in vgseg_idx_map.values():
    for HT in ('head','tail',):
        # init save
        oVgseg_ref_SK = HT+'oVgsegs_ref'
        vgseg[oVgseg_ref_SK] = {}
        #/
        
        oVgseg_idx = vgseg[HT+'_refGseg']
        if oVgseg_idx == None: continue #skip if no ref connection was established (happens in e.g. chrom ends)
        oVgseg = vgseg_idx_map[oVgseg_idx]
        
        # Check if either vgseg is repeatMasked, then skip
        #if vgseg['masked_covFrac'] >= 0.7 or oVgseg['masked_covFrac'] >= 0.7: continue
        #/
        
        # parse read ACs at each vgseg. then find which have same AC at both vgsegs
        vgsegs_qnameAC_supps_pre = {}
        for key,vgseg_var in ( ['vgseg',vgseg] , ['oVgseg',oVgseg] ):
            for oS,idx,oE in rangeOverlaps_lookup(vgseg_var['rcoords'],qname_rsq_idx_map_ROMD,1000,map_lookup_key=vgseg_var['rname']):
                if (oE-oS) < 500: continue #require X numBP (handle ref maskings which were deleted)
                rsqE = qname_rsq_idx_map[idx]
                qname = rsqE['qname']
                if 'memberOfClust' in rsqE:
                    ACenum = rsqE['memberOfClust']
                    init(vgsegs_qnameAC_supps_pre,qname,{})
                    init(vgsegs_qnameAC_supps_pre[qname],ACenum,set())
                    vgsegs_qnameAC_supps_pre[qname][ACenum].add(key)
        #/
        # Check which qnames had same AC at both vgsegs
        vgsegs_qnameAC_supps = {}
        for qname,ACenum_vgsegOvlps in vgsegs_qnameAC_supps_pre.items():
            for ACenum,vgsegOvlps in ACenum_vgsegOvlps.items():
                if len(vgsegOvlps) == 2:
                    vgsegs_qnameAC_supps[qname] = ACenum
        #/

        # Save
        if vgsegs_qnameAC_supps:
            ### IDE: now, save on same format as for vgsegs inferred by variant calls <-- I.e. with strand info. Default-save and dont care about local inv's here.
            oVgseg_ref_HT = None
            if HT == 'head':    oVgseg_ref_HT = 'tail'
            if HT == 'tail':    oVgseg_ref_HT = 'head'
            init(vgseg[oVgseg_ref_SK],oVgseg_idx,{})
            vgseg[oVgseg_ref_SK][oVgseg_idx][oVgseg_ref_HT] = HOLDER(vgsegs_qnameAC_supps)
        #/
#/

## VARSUPP: Link clusters via rsqElink_idxs
# find vGsegs to each rsqElink
rsqElink_idxs_vGsegs = {}
for vgseg in vgseg_idx_map.values():
    for HT in ('head','tail',):
        for rsqElink_gsegEntry in vgseg[HT+'_rsqElink_entries']:
            rsqElink_idx = rsqElink_gsegEntry[3]
            init(rsqElink_idxs_vGsegs,rsqElink_idx,[])
            rsqElink_idxs_vGsegs[rsqElink_idx].append( [vgseg['idx'],HT,rsqElink_gsegEntry] )
#/
# Traverse rsqElinks, save the link between vGsegs
#@ preflight: init/reset links
for vgseg in vgseg_idx_map.values():
    vgseg['head_oVgsegs'] = {}
    vgseg['tail_oVgsegs'] = {}
#@/
IDE_stats = {'bothRsqEsNotHaveACs':0,'qname_hadNoACs':0,'rsqElink_isACcontained':0,'singleReadInRsqElink':0}
for rsqElink_idx,linkData in rsqElink_idxs_vGsegs.items():
    rsqElink = rsqElinks_idx_map[rsqElink_idx]
    if not rsqElink['isUsed']: continue
    # Check read supp for rsqElink
    rAVA_supps = []
    for strint in ('1','2',):
        rAVA_supps.append( len(rsqElink['rAVA_supp'+strint+'_margin=3000']()) >= readSupp_req )
    if not all(rAVA_supps): continue
    #/
    if len(linkData) == 2:
        vgseg_idx1,HT1,_ = linkData[0]
        vgseg_idx2,HT2,_ = linkData[1]
        vgseg1 = vgseg_idx_map[vgseg_idx1]
        vgseg2 = vgseg_idx_map[vgseg_idx2]
        SK1 = HT1+'_oVgsegs'
        SK2 = HT2+'_oVgsegs'
        
        init(vgseg1,SK1,{})
        init(vgseg1[SK1],vgseg_idx2,{})
        init(vgseg1[SK1][vgseg_idx2],HT2,{})
        vgseg1[SK1][vgseg_idx2][HT2][rsqElink_idx] = linkData
        init(vgseg2,SK2,{})
        init(vgseg2[SK2],vgseg_idx1,{})
        init(vgseg2[SK2][vgseg_idx1],HT1,{})
        vgseg2[SK2][vgseg_idx1][HT1][rsqElink_idx] = linkData
#/
##/


## dump in GFA
if 0 and 'dump GFA?':
    with open('varGraph.ACs.noACcontain.noDups.gfa','w') as nf:
        # Write SEGMENT section: col1=S,col2=name,col3=SEQ,colN="tag"[dp:i:50 is depth/cov 50 at ctg. Google GFA spec for more tags]
        for vgseg in vgseg_idx_map.values():
            rname=vgseg['rname']
            rstart,rend=vgseg['rcoords']
            vgseg_id = rname+':'+str(rstart)+'-'+str(rend)
            writeArr = ['S',vgseg_id,ref_seqs[rname][rstart:rend],'dp:i:'+str(int(vgseg['cov']))]
            nf.write('\t'.join(writeArr)+'\n')
        #/
        
        # Write variant connections
        num_HTs = {}
        for vgseg in vgseg_idx_map.values():
            rname=vgseg['rname']
            rstart,rend=vgseg['rcoords']
            vgseg_id = rname+':'+str(rstart)+'-'+str(rend)
            writeEntries = []
            for HT in ('head','tail',):
                HT_ = None
                if HT == 'head':        HT_ = '-' #gfa spec. super unintuitive..
                if HT == 'tail':        HT_ = '+'
                for ovgseg_idx,conn_data in vgseg[HT+'_oVgsegs'].items():
                    num_oHTs = len(conn_data)
                    init(num_HTs,num_oHTs,0)
                    num_HTs[num_oHTs] += 1
                    for oHT,data in conn_data.items():
                        oHT_ = None
                        if oHT == 'head':        oHT_ = '+'
                        if oHT == 'tail':        oHT_ = '-'
                        
                        ovgseg = vgseg_idx_map[ovgseg_idx]
                        ovgseg_id = ovgseg['rname']+':'+str(ovgseg['rcoords'][0])+'-'+str(ovgseg['rcoords'][-1])
                        
                        conn_id = '||'.join(list(map(str,[vgseg['idx'],ovgseg_idx])))
                        
                        writeArr = ['L',vgseg_id,HT_,ovgseg_id,oHT_,'0M','ID:Z:'+conn_id]
                        nf.write('\t'.join(writeArr)+'\n')
        #/
        
        # Write refSupp connections
        num_HTs = {}
        for vgseg in vgseg_idx_map.values():
            rname=vgseg['rname']
            rstart,rend=vgseg['rcoords']
            vgseg_id = rname+':'+str(rstart)+'-'+str(rend)
            writeEntries = []
            for HT in ('head','tail',):
                HT_ = None
                if HT == 'head':        HT_ = '-' #gfa spec. super unintuitive..
                if HT == 'tail':        HT_ = '+'
                for ovgseg_idx,conn_data in vgseg[HT+'oVgsegs_ref'].items():
                    for oHT,data in conn_data.items():
                        oHT_ = None
                        if oHT == 'head':        oHT_ = '+'
                        if oHT == 'tail':        oHT_ = '-'
                        
                        ovgseg = vgseg_idx_map[ovgseg_idx]
                        ovgseg_id = ovgseg['rname']+':'+str(ovgseg['rcoords'][0])+'-'+str(ovgseg['rcoords'][-1])
                        
                        writeArr = ['L',vgseg_id,HT_,ovgseg_id,oHT_,'0M']
                        nf.write('\t'.join(writeArr)+'\n')
        #/
##/IDE
###/VGSEGS

### Calc chrom-densities of stuff, like LTR insertion rate.
# save INS-associated rBPs as idx-map
INSBPs_idx_map = {}
rname_BP_numBP_perClassi = {}
for rsqElink_idx,rsqElink in rsqElinks_idx_map.items():
    if not rsqElink['isUsed']: continue
    qname = rsqElink['qname']
    qSpan = [min([rsqElink['qcoord1'],rsqElink['qcoord2']]),max([rsqElink['qcoord1'],rsqElink['qcoord2']])]
    
    rsqE1 = qname_rsq_idx_map[ rsqElink['rsqE_idx1'] ]
    rsqE2 = qname_rsq_idx_map[ rsqElink['rsqE_idx2'] ]
    rBPclust1 = rsqElink['rBPcov1_clusts']()
    rBPclust2 = rsqElink['rBPcov2_clusts']()
    
    # Parse rsqElink INS content
    INScontent,INScontent_detailed,masked_covFrac = parse_qspanContent(qname,qSpan,return_maskedCovFrac=True,returnGFRO_covFrac=True,requireGFRO_numBPorCovFrac=[100,0.6],skipNumBPBelow=1)
    if not INScontent: continue
    if max(INScontent.values()) < 1000: continue
    #/
    # Add half of INS at each BP
    for intstr in ('1','2',):
        rname = rsqElink['rname'+intstr]
        if not rname in chroms_to_plot: continue
        rBP = rsqElink['rBP'+intstr]
        
        tmp_INSBP = {'rname':rname,'rcoords':[rBP,rBP],'rsqElink_idx':rsqElink_idx,'intstr':intstr}
        assignIdx(tmp_INSBP)
        INSBPs_idx_map[tmp_INSBP['idx']] = tmp_INSBP
    #@ Save stats for maskClassi numBasepair abu + numBreakpoints at rnames
    for intstr in ('1','2',):
        rname = rsqElink['rname'+intstr]
        if not rname in chroms: continue
        for maskClassi,numBP in INScontent.items():
            if numBP >= 1000:
                init(rname_BP_numBP_perClassi,rname,{})
                init(rname_BP_numBP_perClassi[rname],maskClassi,[])
                rname_BP_numBP_perClassi[rname][maskClassi].append(numBP/2) #divide by 2, since save numBP per breakpoint
    #/@
    #/
INSBPs_idx_map_ROMD = rangeOverlaps_makeMappingDict(INSBPs_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
#@ normalize numBP to INS breakpoints by 1/ chromLen, 2/ chromRefAbu

#@@ find refMask which is within range of rsqElinks (e.g. calling-range)
masking_refInfo_perChrom_callingRangees = {}
for rsqElink_idx,rsqElink in rsqElinks_idx_map.items():
    if not rsqElink['rname1'] == rsqElink['rname2']: continue
    #if not abs(rsqElink['rBP1']-rsqElink['rBP2']) < 1000000: continue
    search_margin = 100000
    rcoords = [rsqElink['rBP1']-search_margin,rsqElink['rBP2']+search_margin]
    mask_ovlps = rangeOverlaps_lookup(rcoords,masking_idx_map_ROMD,100,map_lookup_key=rsqElink['rname1'])
    for oS,maskIdx,oE in mask_ovlps:
        masking = masking_idx_map[maskIdx]
        Mclassi = masking['class']
        Mrname = masking['rname']
        Mrcoords = masking['rcoords']
        init(masking_refInfo_perChrom_callingRangees,Mrname,{})
        init(masking_refInfo_perChrom_callingRangees[Mrname],Mclassi,[])
        masking_refInfo_perChrom_callingRangees[Mrname][Mclassi].append([Mrcoords[0],str(maskIdx)+'||'+str(rsqElink_idx),Mrcoords[-1]])

masking_refInfo_perChrom_callingRanges = {}
for rname,classis_rangees in masking_refInfo_perChrom_callingRangees.items():
    for classi,rangees in classis_rangees.items():
        ovlps = computeOvlpRanges_wIdxs(rangees)
        for rangee in rangees:
            init(masking_refInfo_perChrom_callingRanges,rname,{})
            init(masking_refInfo_perChrom_callingRanges[rname],classi,0)
            masking_refInfo_perChrom_callingRanges[rname][classi] += rangee[-1]-rangee[0]
#@@/

rname_BP_numBP_perClassi_num = {}
rname_BP_numBP_perClassi_normChromLen = {}
rname_BP_numBP_perClassi_normChromRefAbu = {}
for rname,classi_numBPs in rname_BP_numBP_perClassi.items():
    for classi,numBPs in classi_numBPs.items():
        if (classi.find('LTR')==-1 and classi.find('Jockey')==-1): continue
        init(rname_BP_numBP_perClassi_num,rname,{})
        init(rname_BP_numBP_perClassi_normChromLen,rname,{})
        init(rname_BP_numBP_perClassi_normChromRefAbu,rname,{})
        rname_BP_numBP_perClassi_num[rname][classi] = len(numBPs)
        rname_BP_numBP_perClassi_normChromLen[rname][classi] = (sum(numBPs) / refLens[rname])*10**6
        try:
            if 0 and 'old':
                rname_BP_numBP_perClassi_normChromRefAbu[rname][classi] = sum(numBPs) / masking_refInfo_perChrom[rname][classi]
            if 1 and 'new':
                rname_BP_numBP_perClassi_normChromRefAbu[rname][classi] = sum(numBPs) / masking_refInfo_perChrom_callingRanges[rname][classi]
            
        except:
            rname_BP_numBP_perClassi_normChromRefAbu[rname][classi] = -(sum(numBPs) / refLens[rname]*0.1)
rname_BP_numBP_perClassi_num['.4'] = rname_BP_numBP_perClassi_num['4']
del (rname_BP_numBP_perClassi_num['4'])
rname_BP_numBP_perClassi_normChromRefAbu['.4'] = rname_BP_numBP_perClassi_normChromRefAbu['4']
del (rname_BP_numBP_perClassi_normChromRefAbu['4'])
rname_BP_numBP_perClassi_normChromLen['.4'] = rname_BP_numBP_perClassi_normChromLen['4']
del (rname_BP_numBP_perClassi_normChromLen['4'])
if 0:
    plotEventDistribution(rname_BP_numBP_perClassi_num,figsize=(10,10),xticksrotate_deg=90,barScale=0.15,tight_layout=True)
    plotEventDistribution(rname_BP_numBP_perClassi_normChromLen,figsize=(10,10),xticksrotate_deg=90,barScale=0.25,tight_layout=True)
    plotEventDistribution(rname_BP_numBP_perClassi_normChromRefAbu,figsize=(10,10),xticksrotate_deg=90,barScale=0.15,tight_layout=True)
#@/
#/
# Define bins at chroms, scan over rRanges and check dens
binwin_step = 100000
binwin_size = binwin_step
rname_INSbinwins = {}
rsqElinks_parsed = set()
mask_types_inserted = {}
mask_types_ref = {}
for rname in chroms_to_plot:
    for binwin_start in range(0,refLens[rname],binwin_step):
        binwin_rcoords = [binwin_start,binwin_start+binwin_size]
        # parse INSes in binwin
        """
        INSBP_ovlps = rangeOverlaps_lookup(binwin_rcoords,INSBPs_idx_map_ROMD,1000,map_lookup_key=rname)
        binwin_numINSBPs = len(INSBP_ovlps)
        tmp_INSbinwin = {'rcoords':binwin_rcoords,'numINSBPs':binwin_numINSBPs}
        """
        ## Parse number of repeats+Gypsy
        numRepeats = 0
        numGypsy = 0
        ovlps1 = rangeOverlaps_lookup(binwin_rcoords,rsqElinks_idx_map_ROMD1,1000,map_lookup_key=rname)
        ovlps2 = rangeOverlaps_lookup(binwin_rcoords,rsqElinks_idx_map_ROMD2,1000,map_lookup_key=rname)
        for oS,idx,oE in ovlps1+ovlps2:
            if idx in rsqElinks_parsed: continue
            rsqElinks_parsed.add(idx)
            rsqElink = rsqElinks_idx_map[idx]
            if rsqElink['len_q'] < 1000: continue
            if rsqElink['rname1'] != rsqElink['rname2'] or abs(rsqElink['rBP1']-rsqElink['rBP2']) >= 1000: continue
            INScontent,INScontent_detailed,masked_covFrac = parse_qspanContent(rsqElink['qname'],rsqElink['qcoords'],return_maskedCovFrac=True,skipNumBPBelow=1)
            if INScontent == None or not INScontent:
                INScontent = {'_none_':rsqElink['len_q']}
            topMaskClassi,topMaskClassiFrac = sorted(INScontent.items(), key=lambda x: x[1], reverse=True)[0]
            
            if masked_covFrac >= 0.7 and topMaskClassiFrac >= 0.9:
                init(mask_types_inserted,topMaskClassi,0)
                mask_types_inserted[topMaskClassi] += 1
                numRepeats += 1
                if topMaskClassi == 'LTR/Gypsy':
                    numGypsy += 1
        ##/
        ## Parse reference-annotation of repeats+Gypsy
        numRepeats_ref = 0
        numGypsy_ref = 0
        for oS,idx,oE in rangeOverlaps_lookup(binwin_rcoords,masking_idx_map_ROMD,100,map_lookup_key=rname):
            masking = masking_idx_map[idx]
            maskClassi = masking['class']
            if not maskClassi in mask_types_inserted: continue
            if masking['len'] >= 1000:
                init(mask_types_ref,maskClassi,0)
                mask_types_ref[maskClassi] += 1
                numRepeats_ref += 1
                if maskClassi == 'LTR/Gypsy':
                    numGypsy_ref += 1
        ##/
        
        tmp_INSbinwin = {'rcoords':binwin_rcoords,'numRepeats':numRepeats,'numGypsy':numGypsy,'numRepeats_ref':numRepeats_ref,'numGypsy_ref':numGypsy_ref}
        init(rname_INSbinwins,rname,[])
        rname_INSbinwins[rname].append(tmp_INSbinwin)
#/
## NTS: should dump these and plot using pyGenomeTracks. Can plot easily there the INS dens, arcs for TRAs, DUPclusts, etc.
if 0 and 'dump for pygenometracks':
    output_dir_pygenometracks = 'pygenometracks'
    mkdir(output_dir_pygenometracks)
    # INSbins
    for key in ('numRepeats','numGypsy','numRepeats_ref','numGypsy_ref'):
        nf = open(output_dir_pygenometracks+'/'+'INSbins.'+key+'.100kb.bedgraph','w')
        for rname in sorted(rname_INSbinwins):
            for INSbin in sorted(rname_INSbinwins[rname],key=lambda x: x['rcoords'][0]):
                writeArr = [ rname,INSbin['rcoords'][0],INSbin['rcoords'][-1],float(INSbin[key]) ]
                nf.write('\t'.join(map(str,writeArr))+'\n')
        nf.close()
    #/
    # DUPclusts
    nf = open(output_dir_pygenometracks+'/'+'DCros.np_pbRF.bed','w')
    tmp_written = 0
    for DCro in DCro_idx_map.values():
        # Filter
        #if not len(DCro['anchors']) >= 1: continue
        #if not len(DCro['columns']) >= 2: continue
        if not len(DCro['DC_idxs']) >= 2: continue
        direction1,direction2 = '+','+'# DUMMY WRITE AS FW AT THIS POINT. such a hassle to parse direction from what I have.
        top_col_cov = sorted(DCro['columns'],key=lambda x: x['cov'], reverse=True)[0]['cov']
        # "chrom1","start1","end1","chrom2","start2","end2","name","score","strand1","strand2","samplenum","supp"
        writeArr = [DCro['rname'],DCro['rcoords'][0],DCro['rcoords'][-1],top_col_cov]

        nf.write('\t'.join(map(str,writeArr))+'\n')
        #/
        tmp_written += 1
    nf.close()
    #/
    # Large-scale rearrs (replace RSQELACPs)
    nf_sameChrom = open(output_dir_pygenometracks+'/'+'large.np_pbRF.'+classi+'.sameChrom.arcs','w')
    nf_diffChrom = open(output_dir_pygenometracks+'/'+'large.np_pbRF.'+classi+'.diffChrom.arcs','w')
    for rsqElink in rsqElinks_idx_map.values():
        # Filters
        if not rsqElink['isUsed']: continue
        rsqE1 = qname_rsq_idx_map[ rsqElink['rsqE_idx1'] ]
        rsqE2 = qname_rsq_idx_map[ rsqElink['rsqE_idx2'] ]
        if 'memberOfDUPs' in rsqE1 and 'memberOfDUPs' in rsqE2 and rsqE1['rname'] == rsqE2['rname']: continue
        if not (rsqElink['rname1'] != rsqElink['rname2'] or rsqElink['len_r'] >= 100000): continue
        if rsqElink['isSwitchOfContainedAC']: continue
        #/
        
        # write for diffChrom
        if rsqElink['rname1'] != rsqElink['rname2']:
            writeArr = [rsqElink['rname1'],rsqElink['rBP1']-1000,rsqElink['rBP1']-1000,rsqElink['rname1'],rsqElink['rBP1']+1000,rsqElink['rBP1']+1000,1]
            nf_diffChrom.write('\t'.join(map(str,writeArr))+'\n')
            writeArr = [rsqElink['rname2'],rsqElink['rBP2']-1000,rsqElink['rBP2']-1000,rsqElink['rname2'],rsqElink['rBP2']+1000,rsqElink['rBP2']+1000,1]
            nf_diffChrom.write('\t'.join(map(str,writeArr))+'\n')
        else: #write for sameChrom
            writeArr = [rsqElink['rname1'],rsqElink['rBP1'],rsqElink['rBP1']+1,rsqElink['rname2'],rsqElink['rBP2'],rsqElink['rBP2']+1,1]
            nf_sameChrom.write('\t'.join(map(str,writeArr))+'\n')
            #/
    
    nf_sameChrom.close()
    nf_diffChrom.close()
    #/
#/
###/

### GENE FUNCTION ANALYSIS (22 sep 2020)
## Find genes which are not disrupted
INFO_rBP_numSupps = []
INFO_qBP_numSupps = []
INFO_disruptions_count = {}
IDE_disrupted_genes = {}
IDE_covDel_GTFs = {}

refLacking_req = 20 # if site has  less than this number of reads supportive of reference, determine it disrupted
varSupp_req = 10 # if site has this number of reads to variant, consider it rearranged
for gene_id,GTF_idxs in gene_GFRO_idxs.items():
    # parse disruptions per gene GTF entry
    GTF_types_disruptions = {}
    for GTF_idx in GTF_idxs:
        GTF = GFRO_idx_map[GTF_idx]
        typee = GTF['type']
        rname = GTF['rname']
        rcoords = GTF['rcoords']
        
        """
        # Calc covFrac & numBP missing of feature by rsqE's
        ovlps = rangeOverlaps_lookup(rcoords,qname_rsq_idx_map_ROMD,1000,map_lookup_key=rname)
        rRanges = []
        for oS,rsqE_idx,oE in ovlps:
            if oE==oS:  oE += 1 #handle fuckup ovlpRanges
            rRanges.append([oS,rsqE_idx,oE])
        ovlps2 = computeOvlpRanges_wIdxs(rRanges)
        numBP_ovlp = 0
        for rangee in ovlps2:
            numBP_ovlp += rangee[-1]-rangee[0]
        ovlp_covFrac = numBP_ovlp / (rcoords[-1]-rcoords[0])
        numBP_missing = (rcoords[-1]-rcoords[0])-numBP_ovlp
        #/
        """
        
        # Check for disruptions (rsqElinks)
        supps_and_disrupts = {} #rsqElink_idx -> BP_strint -> refSupp+rAVAsupp data
        for strint,ROMD in [ ('1',rsqElinks_idx_map_ROMD1,) , ('2',rsqElinks_idx_map_ROMD2,) ]:
            rsqElink_ovlps = rangeOverlaps_lookup(rcoords,ROMD,1000,map_lookup_key=rname)
            for oS,rsqElink_idx,oE in rsqElink_ovlps:
                rsqElink = rsqElinks_idx_map[rsqElink_idx]
                if not rsqElink['isUsed']: continue
            
                # parse q1+r1/q2+r2 breakpoint data
                rBP = rname_rsqE_rBPs_merge_idx_map[rsqElink['rBP_idx'+strint]]
                qBP = qname_rsqE_BPs_idx_map[rsqElink['qBP_idx'+strint]]
                
                rBP_numSupp = len(rBP['supp_rBPspan']())
                #qBP_numSupp = len(qBP['supp_rAVA']())
                qBP_numSupp = len(rsqElink['rAVA_supp'+strint+'_margin=1000']())
                
                init(supps_and_disrupts,rsqElink_idx,{})
                supps_and_disrupts[rsqElink_idx][strint] = [rBP_numSupp,qBP_numSupp]
                
                
                if 'IDE-plot':
                    INFO_qBP_numSupps.append(qBP_numSupp)
                    if qBP_numSupp > varSupp_req:
                        INFO_rBP_numSupps.append(rBP_numSupp)
                    
                    
                if rBP_numSupp <= refLacking_req and qBP_numSupp > varSupp_req:
                    init(GTF_types_disruptions,typee,set())
        #/
        
        # Check if TSS, then scan for upstream rearr
        if typee == 'start_codon':
            gene_GTF = GTF_idx_map[gene_GTF_idx_map[GTF['gene_id']]]
            
            US_rangee = [None,None]
            if strand_isfv(gene_GTF['strand']):     US_rangee = [rcoords[0]-1000,rcoords[0]]
            elif strand_isrv(gene_GTF['strand']):   US_rangee = [rcoords[-1],rcoords[-1]+1000]
            
            for strint,ROMD in [ ('1',rsqElinks_idx_map_ROMD1,) , ('2',rsqElinks_idx_map_ROMD2,) ]:
                rsqElink_ovlps = rangeOverlaps_lookup(US_rangee,ROMD,1000,map_lookup_key=rname)
                for oS,rsqElink_idx,oE in rsqElink_ovlps:
                    rsqElink = rsqElinks_idx_map[rsqElink_idx]
                    if not rsqElink['isUsed']: continue
                    
                    # parse q1+r1/q2+r2 breakpoint data
                    rBP = rname_rsqE_rBPs_merge_idx_map[rsqElink['rBP_idx'+strint]]
                    qBP = qname_rsqE_BPs_idx_map[rsqElink['qBP_idx'+strint]]
                    
                    rBP_numSupp = len(rBP['supp_rBPspan']())
                    #qBP_numSupp = len(qBP['supp_rAVA']())
                    qBP_numSupp = len(rsqElink['rAVA_supp'+strint+'_margin=1000']())
                    
                    init(supps_and_disrupts,rsqElink_idx,{})
                    supps_and_disrupts[rsqElink_idx][strint] = [rBP_numSupp,qBP_numSupp]
                    
                    
                    if 'IDE-plot':
                        INFO_qBP_numSupps.append(qBP_numSupp)
                        if qBP_numSupp > varSupp_req:
                            INFO_rBP_numSupps.append(rBP_numSupp)
                        
                        
                    if rBP_numSupp <= refLacking_req and qBP_numSupp > varSupp_req:
                        init(GTF_types_disruptions,'gene_US',set())
        #/
        
        # Check if disruption by coverage
        if not GTF_types_disruptions and not typee in ('mRNA','gene'):
            # ugly-parse covered coords of GTF feature
            coords_cov = set()
            for oS,idx,oE in rangeOverlaps_lookup(rcoords,BGraw_idx_map_ROMD,100,map_lookup_key=rname):
                BG = BGraw_idx_map[idx]
                if BG['cov'] >= 1:
                    for coord in range(oS,oE+1):
                        coords_cov.add(coord)
            #/
            # ugly-parse non-covered coords of GTF feature
            coords_noCov = []
            for coord in range(rcoords[0],rcoords[-1]+1):
                if not coord in coords_cov:
                    coords_noCov.append(coord)
            #/
            # chain coords, flag if >X numBP without cov
            coords_chained = []
            for coord_enum,coord in enumerate(coords_noCov):
                if coord_enum == 0: continue
                if coords_chained and coords_chained[-1][-1] == (coord-1):
                    coords_chained[-1][-1] = coord
                else:
                    coords_chained.append([coord,coord])
                    
            for rangee in coords_chained:
                if rangee[-1]-rangee[0] >= 10:
                    init(IDE_covDel_GTFs,gene_id,[])
                    IDE_covDel_GTFs[gene_id].append([rname,rcoords,typee])
            #/
        #/
        # Check if disruption by coverage, TSS upstream
        if not GTF_types_disruptions and typee == 'start_codon':
            gene_GTF = GTF_idx_map[gene_GTF_idx_map[GTF['gene_id']]]
            
            US_rangee = [None,None]
            if strand_isfv(gene_GTF['strand']):     US_rangee = [rcoords[0]-1000,rcoords[0]]
            elif strand_isrv(gene_GTF['strand']):   US_rangee = [rcoords[-1],rcoords[-1]+1000]
            
            # ugly-parse covered coords of GTF feature
            coords_cov = set()
            for oS,idx,oE in rangeOverlaps_lookup(US_rangee,BGraw_idx_map_ROMD,100,map_lookup_key=rname):
                BG = BGraw_idx_map[idx]
                if BG['cov'] >= 1 or BG['masked_covFrac'] >= 0.7: #add bases as covered if its repeat. Dont care to call in repeats?
                    for coord in range(oS,oE+1):
                        coords_cov.add(coord)
            #/
            # ugly-parse non-covered coords of GTF feature
            coords_noCov = []
            for coord in range(US_rangee[0],US_rangee[-1]+1):
                if not coord in coords_cov:
                    coords_noCov.append(coord)
            #/
            # chain coords, flag if >X numBP without cov
            coords_chained = []
            for coord_enum,coord in enumerate(coords_noCov):
                if coord_enum == 0: continue
                if coords_chained and coords_chained[-1][-1] == (coord-1):
                    coords_chained[-1][-1] = coord
                else:
                    coords_chained.append([coord,coord])
                    
            for rangee in coords_chained:
                if rangee[-1]-rangee[0] >= 10:
                    init(IDE_covDel_GTFs,gene_id,[])
                    IDE_covDel_GTFs[gene_id].append([rname,rcoords,typee])
            #/
        #/
        
    disrupted_feats_id = '||'.join(GTF_types_disruptions)
    init(INFO_disruptions_count,disrupted_feats_id,0)
    INFO_disruptions_count[disrupted_feats_id] += 1
    
    # Finally, label disrupted status
    for key in ('exon','start_codon','stop_codon','ncRNA','snoRNA','snRNA','miRNA','gene_US',):
        if disrupted_feats_id.find(key) != -1:
            IDE_disrupted_genes[gene_id] = disrupted_feats_id
    #/
    #/  
##/

## Find copynumber per gene (03/12/2020)
gene_covs_CNest = {}
for gene_id,GTF_idx in gene_GTF_idx_map.items():
    GTF = GTF_idx_map[GTF_idx]
    rname = GTF['rname']
    rcoords = GTF['rcoords']
    
    cov_med,cov_mean,covs = parse_rBGcov(rname,rcoords,skipMaskedCovFrac=0.3)
    
    CNest = None
    if rname in autosomes+['4']:
        CNest = round(cov_med / (genomeAverageCov_A/4),2)
    elif rname in X:
        CNest = round(cov_med / (genomeAverageCov_X/2),2)
    
    gene_covs_CNest[gene_id] = [cov_med,CNest]
##/

## Compile list:
if 0:
    with open(MS_output+'/'+'gene_CNest_disruptStatus.list','w') as nf:
        nf.write('\t'.join(map(str,['gene_id','rname','CNest_raw','CNest','geneDisrupted','gainLoss']))+'\n')
        for gene_id,data in gene_covs_CNest.items():
            # Check if gene was disrupted
            geneDisrupted = 0
            if gene_id in set(IDE_disrupted_genes).union(set(IDE_covDel_GTFs)):
                geneDisrupted = 1
        
            rname = GTF_idx_map[gene_GTF_idx_map[gene_id]]['rname']
            CNest = data[1]
            CNest_rounded = None
            if CNest != None:     CNest_rounded = round(CNest)
            
            gainLoss = None
            if rname in autosomes+['4']:
                if CNest_rounded < 4:        gainLoss = 'loss'
                if CNest_rounded > 4:        gainLoss = 'gain'
            elif rname in X:
                if CNest_rounded < 2:       gainLoss = 'loss'
                if CNest_rounded > 2:       gainLoss = 'gain'
            
            writeArr = [gene_id,rname,CNest,CNest_rounded,geneDisrupted,gainLoss]
            nf.write('\t'.join(map(str,writeArr))+'\n')
##/
###/
sys.exit()



### CIRCOS_genome
if 0 and 'pycircos':
    import pandas
    import pycircos

    tmp_chroms = []
    tmp_lens = []
    for chrom,chromLen in sorted(refLens.items(),key=lambda x: x[0]):
        if not chrom in chroms_to_plot: continue
        tmp_chroms.append(chrom)
        tmp_lens.append(chromLen)
        
    chroms_df = pandas.DataFrame(tmp_lens,columns=['length'],index=tmp_chroms)

    figSize = 15
    cg = pycircos.Circos(chroms_df,gap=2,figsize=(figSize,figSize))
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    BASE_rad = figSize/2
    OFFSET_rad = 0
    # ------------------------------------ TRAs
    minJumpDist = 100000
    for qname,ACs_enumed in qname_rsqE_adj_clusts.items():
        comps_done = set()
        # Go through ACs and check supp
        for ACenum,AC in ACs_enumed.items():
            if not 'switchSupp_rAVA' in AC: continue
            for comp_id,comp_data in AC['switchSupp_rAVA'].items():
                if comp_id in comps_done: continue
                comps_done.add(comp_id)
                # Parse ACs in comp, check their dist & skip if not a "large rearr"
                enum1,enum2 = list(map(int,comp_id.split('||')))
                AC1,AC2 = ACs_enumed[enum1],ACs_enumed[enum2]
                if not ((AC1['rname'] != AC2['rname']) or getRangeOvlp(AC1['rcoords'],AC2['rcoords']) < -500000): continue
                if not (AC1['rname'] in chroms_to_plot and AC2['rname'] in chroms_to_plot): continue
                #/
                numSupp_qBPs_span = len(comp_data['supp_span']())
                numSupp_qBP_prev = len(comp_data['supp_prev']())
                numSupp_qBP_cur = len(comp_data['supp_cur']())
                
                if not (numSupp_qBP_prev >= 10 and numSupp_qBP_cur >= 10): continue
                
                AC1_rstart,AC1_rend= AC1['rcoords']
                AC1_rname = AC1['rname']
                AC2_rstart,AC2_rend= AC2['rcoords']
                AC2_rname = AC2['rname']
                if (getRangeOvlp(AC1['rcoords'],AC2['rcoords']) < -minJumpDist) or (AC1_rname != AC2_rname):
                    # draw links [ [rname1,rname2] , [start1,start2] , [end1,end2] ]
                    if AC1_rname != AC2_rname:
                        cg.draw_link(BASE_rad+OFFSET_rad,[AC1_rname,AC2_rname],[AC1_rstart,AC2_rstart],[AC1_rend,AC2_rend],color='purple',alpha=0.5)
                    else:
                        cg.draw_link(BASE_rad+OFFSET_rad,[AC1_rname,AC2_rname],[AC1_rstart,AC2_rstart],[AC1_rend,AC2_rend],color='orange',alpha=0.5)
    # ------------------------------------/TRAs
    # ------------------------------------ CHROM
    cg.draw_scaffold(BASE_rad+OFFSET_rad+0.1, 0.2)
    cg.draw_ticks(BASE_rad+OFFSET_rad+0.3, 0.4, inside=False, tick_gap=5000000)
    cg.draw_scaffold_ids(BASE_rad+OFFSET_rad+2.5, inside=False, fontsize=10)
    #/------------------------------------ CHROM
    # ------------------------------------ DUPs
    OFFSET_rad += 1.5
    cg.draw_scaffold(BASE_rad + OFFSET_rad, 0.01)
    dupNorm = 7 # divide by this num to obtain fraction
    clust_min_size = 5000
    tmp_DUPs = {'chrom':[],'rstarts':[],'rends':[],'fracs':[]}
    IDE_qname_DUPclusts = {}
    for qname,dupClusts in qname_rsqE_dupClusts.items():
        for DC_enum,DC in dupClusts.items():
            DC_rangee = DC['clust_raw']
            DC_rname = DC_rangee[1]
            if not DC_rname in chroms_to_plot: continue
            rstart,rend = DC_rangee[0],DC_rangee[-1]
            DC_CNe = DC['CN_est']
            DC_size = rend-rstart
            if DC_size >= clust_min_size:
                IDE_qname_DUPclusts[qname] = dupClusts
                tmp_DUPs['chrom'].append(DC_rname)
                tmp_DUPs['rstarts'].append(rstart)
                tmp_DUPs['rends'].append(rend)
                tmp_DUPs['fracs'].append( DC_CNe/dupNorm )

    DUPs_df = pandas.DataFrame(tmp_DUPs)
    cg.fill_between(BASE_rad +OFFSET_rad, DUPs_df,start='rstarts',end='rends',score='fracs',scale=2,facecolor='red',alpha=0.5,cutoff=0)
    #/------------------------------------ DUPs
    
    # ------------------------------------ INScallDens
    OFFSET_rad += 2
    cg.draw_scaffold(BASE_rad +OFFSET_rad, 0.001)
    eINSavg = 3
    eINSavg = [eINSavg-eINSavg*0.15 , eINSavg+eINSavg*0.15]
    tmp_INS_pan = {'chrom':[],'rstarts':[],'rends':[],'fracs':[]}
    tmp_INS_GYP = {'chrom':[],'rstarts':[],'rends':[],'fracs':[]}
    tmp_INS_stats = {}
    for BG in BG_idx_map.values():
        if not BG['rname'] in tmp_chroms: continue #IDE, for dev
        INSdens = len(BG['INScall_rsqElinks'])
        if INSdens <= 0: continue
        
        tmp_INS_pan['chrom'].append(BG['rname'])
        tmp_INS_pan['rstarts'].append(BG['rcoords'][0])
        tmp_INS_pan['rends'].append(BG['rcoords'][-1])
        tmp_INS_pan['fracs'].append(INSdens)

    INS_pan_df = pandas.DataFrame(tmp_INS_pan)
    cg.fill_between(BASE_rad +OFFSET_rad,INS_pan_df,start='rstarts',end='rends',score='fracs',scale=0.3,facecolor='black',alpha=1,cutoff=0)
    #/------------------------------------ INScallDens
    """
    # ------------------------------------ COVERAGE
    OFFSET_rad += 3
    cg.draw_scaffold(BASE_rad +OFFSET_rad, 0.001)
    eCovAvg = 140
    eCovAvg_rangee = [eCovAvg-eCovAvg*0.15 , eCovAvg+eCovAvg*0.15]
    tmp_CN_abo = {'chrom':[],'rstarts':[],'rends':[],'fracs':[]}
    tmp_CN_below = {'chrom':[],'rstarts':[],'rends':[],'fracs':[]}
    tmp_covStats = {}
    for BG in BG_idx_map.values():
        if not BG['rname'] in tmp_chroms: continue #IDE, for dev
        BG_cov = BG['cov']
        if BG_cov > 1000: BG_cov = 1000 #cap it
        BG_cov = (BG['cov']/eCovAvg) - 1
        
        # Check if lower than eCovAvg vs. above eCovAvg
        if BG_cov < 0:
            saveLoc = tmp_CN_below
        else:
            saveLoc = tmp_CN_abo
        
        saveLoc['chrom'].append(BG['rname'])
        saveLoc['rstarts'].append(BG['rcoords'][0])
        saveLoc['rends'].append(BG['rcoords'][-1])
        saveLoc['fracs'].append(BG_cov)
        
        init(tmp_covStats,BG['rname'],[])
        tmp_covStats[BG['rname']].append(BG['cov'])
        
    CN_abo_df = pandas.DataFrame(tmp_CN_abo)
    CN_below_df = pandas.DataFrame(tmp_CN_below)
    cg.fill_between(BASE_rad +OFFSET_rad,CN_abo_df,start='rstarts',end='rends',score='fracs',scale=0.5,facecolor='red',alpha=0.5,cutoff=0)
    cg.fill_between(BASE_rad +OFFSET_rad,CN_below_df,start='rstarts',end='rends',score='fracs',scale=0.5,facecolor='blue',alpha=0.5,cutoff=0)
    #/------------------------------------ COVERAGE
    """
    plt.show()
    
    if 0 and 'IDE, set eCovAvg':
        for rname,covStat in tmp_covStats.items():
            print(rname,median(covStat))
###/

### CIRCOS, read perspective. Plot read-all-vs-all coverage on read at breakpoints. Plot coverage at read alignments to reference genome.
# Idea: For interesting reads, plot ref in CIRCOS. (to visualize rearrs)
import pandas
import pycircos
repeats_colors = {'DNA/P':'azure','LTR':'lime','LINE1':'neon purple','other':'black'}
colorCycler = ['red','blue','black','gray','purple','orange','green']
numberColors = plt.cm.nipy_spectral(np.linspace(0,1,10)) #start,step,num_to_generate
numberColors_dict = {}
for i in range(len(numberColors)):
    numberColors_dict[i] = numberColors[i]

if 0:
    fig = plt.figure(figsize=(3,1))
    for i,col in numberColors_dict.items():
        plt.bar(i,1,color=col)
        plt.xticks(ticks=list(range(len(numberColors_dict))))
    
parse_regions = True #if we want to zoom in to regions where rearrs occur on ref, instead of looking chrom-wise
parse_regions_fattenRange = 50000 #fatten rname regions by this numBP

CIRCOS_entries = []
#CIRCOS_entries.append(['4b9da816-c65b-4e3a-b05e-7eda6987631f']) #3L_X
#CIRCOS_entries.append(['234040fe-d115-4bd9-ae81-275598af6a1e','99a18323-bc9f-422b-bb01-6e5f1f342f28']) # banX+ban3L
#CIRCOS_entries.append(['7f6aa397-4cc2-4c5a-9c0f-fa43bdf49cc7','4b9da816-c65b-4e3a-b05e-7eda6987631f']) # 3L_mito+3L_X
#CIRCOS_entries.append(['4027726f-d731-4eac-abc7-262e71842bb5','05e5a681-f3e5-4bf8-875b-6d5645cdf227']) # 3L3.3_X+3L0.8_X
#CIRCOS_entries.append(['5ea58da1-4372-461e-b7ae-f5223f67ebdb_Basecall_1D_template']) # 2L->3R->2R. This this one is an orphan event??
#CIRCOS_entries.append(['8ba72d7b-0656-436d-bc07-1d8c10e09c88']),CIRCOS_entries.append(['e7deb45b-1d36-438a-b969-7b3b0461047f']) # c88=[weird mapper. Might be sequencing fuckup. from 2017.] 47f=[from 2020]
#CIRCOS_entries.append(['598f04c7-3361-4f3e-86d5-b10e6b72045a','c4981883-4806-46b8-b118-5ae380f30fae'])
#CIRCOS_entries.append(['8ba72d7b-0656-436d-bc07-1d8c10e09c88'])
CIRCOS_entries.append(['234040fe-d115-4bd9-ae81-275598af6a1e'])
#circosPlot(CIRCOS_entries)
def circosPlot(CIRCOS_entries,qonly=False,cleanPlot=False):
    # Check if input is single read, then convert to arr 
    if type(CIRCOS_entries) == str:                 CIRCOS_entries = [[CIRCOS_entries]]
    # Check if not on chunk->entry format (arr-in-arr)
    if not type(CIRCOS_entries[0]) == list:         CIRCOS_entries = [CIRCOS_entries]
    
    SCALE_plot = 1
    CHROMTICKS_interval,CHROMTICKS_unitLabel,CHROMTICKS_unitSize = 20000,'K',1000
    BASE_rad = 6.5*SCALE_plot
    refCoords_scale = 1 #2/100
    maskCovFrac = 0.5
    cov_norm_to = 20
    for qnames in CIRCOS_entries:
        OFFSET_rad = 0
        
        ## INIT rnames + reads
        # Add chroms, in log10
        qnames_rnames = {} #rname -> qname
        for qname in qnames:
            prev = None
            for enum,rsqE in enumerate(sorted(qname_rsqRanges[qname],key=lambda x: x['qcoords'][0])):
                if rsqE['masked_covFrac'] <= maskCovFrac or rsqE['rname'] == 'rDNA':
                    qnames_rnames[rsqE['rname']] = True
        
        tmp_chroms = []
        tmp_lens = []
        for chrom,chromLen in sorted(refLens.items(),key=lambda x: x[0]):
            if not chrom in chroms_to_plot: continue
            if not chrom in qnames_rnames: continue
            tmp_chroms.append(chrom)
            #tmp_lens.append(chromLen*refCoords_scale)
            tmp_lens.append(chromLen*refCoords_scale)
        
        #--- Check if we wanted to do with chrom extractions, then scrap and redo
        if parse_regions:
            #refCoords_scale = 1 #nullify
    
            # Setup rangees per rname with BASE and rsqE's
            rname_regions_to_parse = {}
            for qname in qnames:
                for enum,rsqE in enumerate(sorted(qname_rsqRanges[qname],key=lambda x: x['qcoords'][0])):
                    if rsqE['masked_covFrac'] <= maskCovFrac or rsqE['rname'] == 'rDNA':
                        rname = rsqE['rname']
                        rstart,rend = rsqE['rcoords']
                        init(rname_regions_to_parse,rname,[ [0,'base',refLens[rname]] ])
                        rstart_fattened = rstart-parse_regions_fattenRange
                        if rstart_fattened < 0: rstart_fattened = 0
                        rend_fattened = rend+parse_regions_fattenRange
                        if rend_fattened > refLens[rname]: rend_fattened = refLens[rname]
                        rname_regions_to_parse[rname].append([rstart_fattened,rsqE['idx'],rend_fattened])
            #/
            # Traverse rnames, compile ranges where rsqE's exist
            rname_regions_PRE = [] # [rname, [rposes], set(rsqIdxs)] to plot at this "artificial chrom"
            for rname, rangees in rname_regions_to_parse.items():
                ovlpRangees = computeOvlpRanges_wIdxs(rangees)
                for rangee in ovlpRangees:
                    tmp_save = {'rname':rname,'rposes':[],'rposes_wFatten':[],'rsqE_idxs':set()}
                    
                    # check if base only, then init new artificialChrom and do nothing more...
                    if len(rangee[1]) == 1 and 'base' in rangee[1]: 
                        if not rname_regions_PRE or len(rname_regions_PRE[-1]) >= 1:
                            rname_regions_PRE.append(tmp_save)
                    else:
                        # Append coords, rsqE_idxs to existing artificialChrom
                        if not rname_regions_PRE:                           rname_regions_PRE.append(tmp_save)
                        if not rname_regions_PRE[-1]['rname'] == rname:     rname_regions_PRE.append(tmp_save) # check if rname was changed
                        
                        for member in rangee[1]:
                            if member in qname_rsq_idx_map:
                                rsqE = qname_rsq_idx_map[member]
                                rname_regions_PRE[-1]['rposes_wFatten'].append(rangee[0])
                                rname_regions_PRE[-1]['rposes_wFatten'].append(rangee[-1])
                                rname_regions_PRE[-1]['rposes'] += rsqE['rcoords']
                                rname_regions_PRE[-1]['rsqE_idxs'].add(member)
            #/
            # Compile ranges into artificial chroms
            rname_regions = {} # rname:rstart-rend -> rsqEs
            rsqE_rnameRegions_map = {} # Map beloning of rsqEs: rsqE -> rname:rstart-rend
            for rname_region_PRE in rname_regions_PRE:
                if not rname_region_PRE['rposes']: continue #skip if no data was added here...
                rname = rname_region_PRE['rname']
                rstart, rend = [min(rname_region_PRE['rposes_wFatten']),max(rname_region_PRE['rposes_wFatten'])]
                rstart_rsqEs = min(rname_region_PRE['rposes'])
                rcoords_offset = rstart_rsqEs+(rstart-rstart_rsqEs)
                rnameRegion_ID = rname+':'+ str(f"{rstart:,}") +'-'+str(f"{rend:,}") 
                rname_regions[rnameRegion_ID] = {'rcoords_offset':rcoords_offset,'rname':rname,'rcoords':[rstart,rend],'rsqE_idxs':rname_region_PRE['rsqE_idxs']}
                for rsqE_idx in rname_regions[rnameRegion_ID]['rsqE_idxs']:
                    rsqE_rnameRegions_map[rsqE_idx] = rnameRegion_ID
            #/
            # Reset plot chromosomes and add ariticialChroms of regions instead..
            tmp_chroms = []
            tmp_lens = []
            for rnameReg,data in rname_regions.items():
                tmp_chroms.append(rnameReg)
                tmp_lens.append(data['rcoords'][-1] - data['rcoords'][0])
            #/
        #/--
        
        # Add reads
        for qname in qnames:
            tmp_chroms.append(qname)
            tmp_lens.append(qname_rsqRanges[qname][0]['read_len'])
            
        # Make DF + add to plot
        chroms_df = pandas.DataFrame(tmp_lens,columns=['length'],index=tmp_chroms)
        cg = pycircos.Circos(chroms_df,gap=2,figsize=(12,12))
        # ------------------------------------ CHROM
        cg.draw_scaffold(BASE_rad+OFFSET_rad+0.1*SCALE_plot, 0.3*SCALE_plot)
        cg.draw_ticks(BASE_rad+OFFSET_rad+0.3*SCALE_plot, 1.5*SCALE_plot, inside=False, tick_gap=CHROMTICKS_interval,unit=CHROMTICKS_unitSize,unit_label=CHROMTICKS_unitLabel)
        cg.draw_scaffold_ids(BASE_rad+OFFSET_rad+4*SCALE_plot, inside=False, fontsize=10, rotation=90)
        #/------------------------------------ CHROM
        ##/
        
        ## Add stuff
        qnames_qcoords_to_rAVA = {}
        for qname in qnames:
            # read rsqEs + qBPs + rBPs
            prev = None
            for enum,rsqE in enumerate(sorted(qname_rsqRanges[qname],key=lambda x: x['qcoords'][0])):
                for coord in rsqE['qcoords']:
                    if coord < 1000: continue
                    if coord > rsqE['read_len']-1000: continue
                    init(qnames_qcoords_to_rAVA,qname,[])
                    qnames_qcoords_to_rAVA[qname].append(coord)
                if rsqE['masked_covFrac'] <= maskCovFrac or rsqE['rname'] == 'rDNA':
                    rname = rsqE['rname']
                    rstart,rend = None,None
                    if strand_isfv(rsqE):   rstart,rend = rsqE['rcoords']
                    if strand_isrv(rsqE):   rend,rstart = rsqE['rcoords']
                    
                    if parse_regions: # Correct coords if we extract!!
                        rname = rsqE_rnameRegions_map[rsqE['idx']] # Update the rname!
                        rstart = rstart - rname_regions[rname]['rcoords_offset']
                        rend = rend - rname_regions[rname]['rcoords_offset']
    
                    cg.draw_link(BASE_rad+OFFSET_rad,[qname,rname],[rsqE['qcoords'][0],rstart*refCoords_scale],[rsqE['qcoords'][-1],rend*refCoords_scale],
                                 color=colorCycler[enum%len(colorCycler)],alpha=0.3)
                    
                    # qBP + rBP
                    if 'qBP_idxs' in rsqE:
                        for qBP,qBP_idx in rsqE['qBP_idxs'].items():
                            qBP = qname_rsqE_BPs_idx_map [qBP_idx]
                            rBP = rname_rsqE_rBPs_merge_idx_map[qBP['rBP_idx']]
                            
                            tmp_qBP = {'qnames':[qname],'qstarts':[qBP['qcoord']-500],'qends':[qBP['qcoord']+500],'fracs':[len(qBP['supp_rAVA']())/150]}
                            cg.fill_between(BASE_rad +0.15*SCALE_plot, tmp_qBP,gid='qnames',start='qstarts',end='qends',score='fracs',scale=0.3,
                                        facecolor='red',alpha=0.8,cutoff=0)

                            if parse_regions and not cleanPlot:
                                rname = rsqE_rnameRegions_map[rsqE['idx']] # Update the rname!
                                rstart = rBP['rBP'] - 500 - rname_regions[rname]['rcoords_offset']
                                rend = rBP['rBP'] + 500 - rname_regions[rname]['rcoords_offset']
                                tmp_rBP = {'rnames':[rname],'rstarts':[rstart],'rends':[rend],'fracs':[len(rBP['supp_rBPspan']())/150]}
                                cg.fill_between(BASE_rad +1.55*SCALE_plot, tmp_rBP,gid='rnames',start='rstarts',end='rends',score='fracs',scale=0.3,
                                                facecolor='green',alpha=1,cutoff=0)
                                
                                ## Try plot larger ref supps
                                ovlps = rangeOverlaps_lookup([ rBP['rBP']-1000,rBP['rBP']+1000 ] ,qname_rsqE_BPs_idx_map_rBPcov_ROMD_noMasked,2000,map_lookup_key=rBP['rname'])
                                reads_per_rBPcov = {}
                                for oS,rsqE_BP_idx,oE in ovlps:
                                    if rsqE_BP_idx == rBP['idx']: continue #skip if self
                                    rsqE_BP_entry = rname_rsqE_rBPs_merge_idx_map[rsqE_BP_idx]
                                    init(reads_per_rBPcov,rsqE_BP_idx,set())
                                    reads_per_rBPcov[rsqE_BP_idx].update(rsqE_BP_entry['supp_rBPcov']())
                                
                                # Plot pairwise cov supps between rBP and orBPs
                                rspans_suppsIntersect = []
                                for orBP_idx,reads in reads_per_rBPcov.items():
                                    orBP = rname_rsqE_rBPs_merge_idx_map[orBP_idx]
                                    suppsIntersect = reads.intersection(set(rBP['supp_rBPcov']()))
                                    if suppsIntersect:
                                        rspans_suppsIntersect.append([ [min([rBP['rBP'],orBP['rBP']])-1000 , max([rBP['rBP'],orBP['rBP']])+1000],len(suppsIntersect)])
                                
                                for tmp_rspan,numSuppsIntersect in rspans_suppsIntersect:
                                    tmp_rBP = {'rnames':[rname],'rstarts':[tmp_rspan[0]-rname_regions[rname]['rcoords_offset']],
                                               'rends':[tmp_rspan[-1]-rname_regions[rname]['rcoords_offset']],'fracs':[numSuppsIntersect/150]}
                                    cg.fill_between(BASE_rad +1.55*SCALE_plot, tmp_rBP,gid='rnames',start='rstarts',end='rends',score='fracs',scale=-0.3,
                                                    facecolor='red',alpha=1,cutoff=0)
                                ##/
            #/
            
            # read annotated masks
            prev = None
            for enum,rsqE in enumerate(sorted(qname_rsqRanges[qname],key=lambda x: x['qcoords'][0])):
                if rsqE['masked_covFrac'] > maskCovFrac:
                    tmp_clusts = {'qnames':[],'qstarts':[],'qends':[],'fracs':[]}
                    tmp_clusts['qnames'].append(qname) #plot on self...
                    tmp_clusts['qstarts'].append(rsqE['qcoords'][0])
                    tmp_clusts['qends'].append(rsqE['qcoords'][-1])
                    tmp_clusts['fracs'].append( 1 )
                    topMaskClassi,numBP = rsqE['maskClassis_numBPs'][0] #is already sorted
                    tmp_color = repeats_colors['other']
                    for key,color in repeats_colors.items():
                        if topMaskClassi.find(key) != -1:
                            tmp_color = color
                            break
                        
                    cg.fill_between(BASE_rad -0.075*SCALE_plot, tmp_clusts,gid='qnames',start='qstarts',end='qends',score='fracs',scale=0.3,
                                facecolor=tmp_color,alpha=0.8,cutoff=0)
            #/
            
            # Clusts
            if qname in qname_rsqE_adj_clusts and not cleanPlot:
                #clust_min_size = 5000
                for clustEnum,clust in qname_rsqE_adj_clusts[qname].items():
                    # plot on q
                    clust_qstart,clust_qend = clust['qcoords']
                    clust_rname = clust['rname']
                    #if not clust_rname in chroms_to_plot: continue
                    #if not clust_qend-clust_qstart >= clust_min_size: continue
                    tmp_clusts = {'qnames':[],'qstarts':[],'qends':[],'fracs':[]}
                    tmp_clusts['qnames'].append(qname) #plot on self...
                    tmp_clusts['qstarts'].append(clust_qstart)
                    tmp_clusts['qends'].append(clust_qend)
                    tmp_clusts['fracs'].append( (clustEnum+1)/len(qname_rsqE_adj_clusts[qname]) )
                    cg.fill_between(BASE_rad -0.33*SCALE_plot, tmp_clusts,gid='qnames',start='qstarts',end='qends',score='fracs',scale=0.3,
                                facecolor=colorCycler[clustEnum],alpha=1,cutoff=0)
                    #/
                    # plot on r
                    tmp_clusts = {'rnames':[clust['rname']],'rstarts':[clust['rcoords'][0]],'rcoords':[clust['rcoords'][-1]],'fracs':[]}
                    if parse_regions: # Correct coords if we extract!!
                        for any_rsqE_idx in clust['members']:
                            rname = rsqE_rnameRegions_map[any_rsqE_idx] # Update the rname!
                            rcoords_offset = rname_regions[rname]['rcoords_offset']
                            tmp_clusts = {'rnames':[rname],'rstarts':[clust['rcoords'][0]-rcoords_offset],'rends':[clust['rcoords'][-1]-rcoords_offset],'fracs':[]}
                            break
    
                    tmp_clusts['fracs'].append( (clustEnum+1)/len(qname_rsqE_adj_clusts[qname]) )
                    cg.fill_between(BASE_rad -0.33*SCALE_plot, tmp_clusts,gid='rnames',start='rstarts',end='rends',score='fracs',scale=0.3,
                                facecolor=colorCycler[clustEnum],alpha=1,cutoff=0)
                    #/
                    # Plot rAVAsupp
                    if not cleanPlot:
                        if not 'switchSupp_rAVA' in clust: continue
                        for clusts_comp_ID,rAVA_supp_entry in clust['switchSupp_rAVA'].items():
                            for qcoord,CP in ( [clust['qcoords'][0],'prev'] , [clust['qcoords'][-1],'cur'] ): 
                                tmp_q = {'qnames':[qname],'qstarts':[qcoord-500],'qends':[qcoord+500],'fracs':[len(rAVA_supp_entry['supp_'+CP]())/150]}
                                cg.fill_between(BASE_rad +0.15*SCALE_plot, tmp_q,gid='qnames',start='qstarts',end='qends',score='fracs',scale=0.3,
                                    facecolor='red',alpha=0.8,cutoff=0)
                    #/
            #/
            
            # dupClusts
            if qname in qname_rsqE_DUPs and not cleanPlot:
                clust_min_size = 5000
                for clustEnum,clust in qname_rsqE_DUPs[qname].items():
                    # plot on q
                    clust_qcoordss = []
                    for rsqE_idx in clust['members']:
                        clust_qcoordss += qname_rsq_idx_map[rsqE_idx]['qcoords']
                    clust_qstart,clust_qend = min(clust_qcoordss),max(clust_qcoordss)
                    clust_rname = clust['rname']
                    if not clust_rname in chroms_to_plot: continue
                    if not clust_qend-clust_qstart >= clust_min_size: continue
                    CNest = qname_rsqE_dupClusts[qname][clustEnum]['CN_est']
                    tmp_clusts = {'qnames':[],'qstarts':[],'qends':[],'fracs':[]}
                    tmp_clusts['qnames'].append(qname) #plot on self...
                    tmp_clusts['qstarts'].append(clust_qstart)
                    tmp_clusts['qends'].append(clust_qend)
                    tmp_clusts['fracs'].append( (clustEnum+1)/len(qname_rsqE_dupClusts[qname]) )
                    try:        tmp_col = numberColors[CNest]
                    except:     tmp_col = 'cyan'
                    cg.fill_between(BASE_rad -1*SCALE_plot, tmp_clusts,gid='qnames',start='qstarts',end='qends',score='fracs',scale=0.3,
                                facecolor=tmp_col,alpha=0.3,cutoff=0,linewidth=1,edgecolor='black',linestyle='-.')
                    #/
                    # plot on r
                    tmp_clusts = {'rnames':[clust['rname']],'rstarts':[clust['rcoords'][0]],'rends':[clust['rcoords'][-1]],'fracs':[]}
                    
                    if parse_regions: # Correct coords if we extract!!
                        for any_rsqE_idx in clust['members']:
                            rname = rsqE_rnameRegions_map[any_rsqE_idx] # Update the rname!
                            rcoords_offset = rname_regions[rname]['rcoords_offset']
                            tmp_clusts = {'rnames':[rname],'rstarts':[clust['rcoords'][0]-rcoords_offset],'rends':[clust['rcoords'][-1]-rcoords_offset],'fracs':[]}
                            break
    
                    tmp_clusts['fracs'].append( (clustEnum+1)/len(qname_rsqE_dupClusts[qname]) )
                    cg.fill_between(BASE_rad -1*SCALE_plot, tmp_clusts,gid='rnames',start='rstarts',end='rends',score='fracs',scale=0.3,
                                facecolor=tmp_col,alpha=0.3,cutoff=0,linewidth=1,edgecolor='black',linestyle='-.')
                    #/
            #/
            
            # other rearrs in range
            if not cleanPlot:
                for rnameRegion,regionData in rname_regions.items():
                    rname = regionData['rname']
                    rcoords = regionData['rcoords']
                    offset_val = regionData['rcoords_offset']
                    
                    rsqElinks_ovlps1 = rangeOverlaps_lookup(rcoords,rsqElinks_idx_map_ROMD1,1000,map_lookup_key=rname)
                    rsqElinks_ovlps2 = rangeOverlaps_lookup(rcoords,rsqElinks_idx_map_ROMD2,1000,map_lookup_key=rname)
    
                    rsqElinks_done = set()
                    for oS,idx,oE in rsqElinks_ovlps1+rsqElinks_ovlps2:
                        if idx in rsqElinks_done: continue
                        rsqElinks_done.add(idx)
                        rsqElink = rsqElinks_idx_map[idx]
                        if qname in rsqElink['reads'](): continue
                        # Check if both rBPs are within range
                        entries = []
                        for intstr in ('1','2',):
                            tmp_rname = rsqElink['rname'+intstr]
                            tmp_rcoord = rsqElink['rBP'+intstr]
                            if tmp_rname == rname and getRangeOvlp([tmp_rcoord]*2,rcoords) >= 0:
                                entries.append([rname,tmp_rcoord])
                        if not entries: continue
                        tmp_rcoords = []
                        for entry in entries:
                            tmp_rcoords.append(entry[1])
                        tmp_rcoords2 = [min(tmp_rcoords),max(tmp_rcoords)]
                        if tmp_rcoords[-1]-tmp_rcoords[0] < 1000:
                            tmp_rcoords2[0] = mean(tmp_rcoords)-500
                            tmp_rcoords2[-1] = mean(tmp_rcoords)+500
                            
                        if len(entries) == 1: 
                            col = 'red'
                            scale = 0.5
                            extraOffset = 0.5
                        elif len(entries) == 2:
                            col = 'blue'
                            scale = 0.1
                            extraOffset = 0
                        else:
                            col = 'yellow'
                            scale = 0.3
                            extraOffset = -0.5
                        
                        # Count how many of the rsqEs are DUPclusts
                        num_isDUPclust = 0
                        for intstr in ('1','2',):
                            for tmp_rsqE in qname_rsq_idx_map[rsqElink['rsqE_idx'+intstr]]:
                                if 'memberOfDUPs' in tmp_rsqE:
                                    num_isDUPclust += 1
                        if num_isDUPclust == 2:
                            col = 'lime'
                        #/
    
                        tmp_BG = {'chrom':[],'rstarts':[],'rends':[],'fracs':[]}
                        tmp_BG['chrom'].append(rnameRegion)
                        tmp_BG['rstarts'].append(min(tmp_rcoords)-offset_val)
                        tmp_BG['rends'].append(max(tmp_rcoords)-offset_val)
                        tmp_BG['fracs'].append(1)
                            
                        BG_df = pandas.DataFrame(tmp_BG)
                        cg.fill_between(BASE_rad +OFFSET_rad+0.3+extraOffset,BG_df,start='rstarts',end='rends',score='fracs',scale=scale,facecolor=col,alpha=0.5,cutoff=0)
                #/
        ##/
       
        # Check if multiple reads, then add pairwise rAVAs
        if 0 and 'rAVA-plot':
            for read1 in qnames:
                for read2 in qnames:
                    if read1 == read2: continue
                    for rAVA_comp_id in rAVA_read_mappings[read1]:
                        if rAVA_comp_id.find(read2) != -1:
                            for entry in rAVA_mappings[rAVA_comp_id]:
                                cg.draw_link(BASE_rad+OFFSET_rad,[entry['qname'],entry['rname']],[entry['qcoords'][0],entry['rcoords'][0]],[entry['qcoords'][-1],entry['rcoords'][-1]],
                                             color='yellow',alpha=0.1)
                            
        #/
        
        ## Add covs
        OFFSET_rad += 3
        cg.draw_scaffold(BASE_rad +OFFSET_rad, 0.001)
        
        # Qnames
        for qname in qnames:
            # Plot pan cov
            qname_BG = get_qcoords_rAVA_bedgraph(qname,bin_size=1000)
            
            tmp_BG = {'chrom':[],'rstarts':[],'rends':[],'fracs':[]}
            for binn in qname_BG:
                binn_cov = binn[1]/200
                
                tmp_BG['chrom'].append(qname)
                tmp_BG['rstarts'].append(binn[0])
                tmp_BG['rends'].append(binn[-1])
                tmp_BG['fracs'].append(binn_cov)
                
            BG_df = pandas.DataFrame(tmp_BG)
            cg.fill_between(BASE_rad +OFFSET_rad,BG_df,start='rstarts',end='rends',score='fracs',scale=0.5,facecolor='blue',alpha=0.5,cutoff=0)
            #/
            
            # Plot specifically around qBPs
            if not cleanPlot:
                # large span
                tmp_cov = {'chrom':[],'rstarts':[],'rends':[],'fracs':[]}
                for qcoord in qnames_qcoords_to_rAVA[qname]:
                    qcoord_range = [qcoord-5000,qcoord+5000]
                    cov = len(get_qcoords_rAVAsupp(qname,qcoord_range)) / 200
                    
                    tmp_cov['chrom'].append(qname)
                    tmp_cov['rstarts'].append(qcoord_range[0])
                    tmp_cov['rends'].append(qcoord_range[-1])
                    tmp_cov['fracs'].append(cov)
                    
                cov_df = pandas.DataFrame(tmp_cov)
                cg.fill_between(BASE_rad +OFFSET_rad,cov_df,start='rstarts',end='rends',score='fracs',scale=0.5,facecolor='red',alpha=0.3,cutoff=0)
                #/
                # short span
                tmp_cov = {'chrom':[],'rstarts':[],'rends':[],'fracs':[]}
                for qcoord in qnames_qcoords_to_rAVA[qname]:
                    qcoord_range = [qcoord-2500,qcoord+2500]
                    cov = len(get_qcoords_rAVAsupp(qname,qcoord_range))
                    if cov > 5: cov = 5
                    tmp_cov['chrom'].append(qname)
                    tmp_cov['rstarts'].append(qcoord_range[0])
                    tmp_cov['rends'].append(qcoord_range[-1])
                    tmp_cov['fracs'].append(-cov/25) #adjust this for the scale
                    
                cov_df = pandas.DataFrame(tmp_cov)
                cg.fill_between(BASE_rad +OFFSET_rad,cov_df,start='rstarts',end='rends',score='fracs',scale=1,facecolor='purple',alpha=1,cutoff=0)
                #/
        #/
        # Rnames
        for rnameRegion,regionData in rname_regions.items():
            rname = regionData['rname']
            rcoords = regionData['rcoords']
            offset_val = regionData['rcoords_offset']
            BG_ovlps = rangeOverlaps_lookup(rcoords,BG_idx_map_ROMD,1000,map_lookup_key=rname)
            
            tmp_BG = {'chrom':[],'rstarts':[],'rends':[],'fracs':[]}
            for oS,idx,oE in BG_ovlps:
                BG = BG_idx_map[idx]
                cov = BG['cov']
                tmp_BG['chrom'].append(rnameRegion)
                tmp_BG['rstarts'].append(oS-offset_val)
                tmp_BG['rends'].append(oE-offset_val)
                tmp_BG['fracs'].append(cov/200)
                
            BG_df = pandas.DataFrame(tmp_BG)
            cg.fill_between(BASE_rad +OFFSET_rad,BG_df,start='rstarts',end='rends',score='fracs',scale=0.5,facecolor='blue',alpha=0.5,cutoff=0)
                
        #/
        
        # Call plot
        plt.show()
###/
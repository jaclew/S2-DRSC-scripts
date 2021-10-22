import os
import pysam
import matplotlib.pyplot as plt
import numpy as np
from statistics import mean,median,stdev
import pickle
import sys

from functions import *

### INPUTS
chroms = ['2L','2R','3L','3R','4','X']
autosomes = ['2L','2R','3L','3R']
X = ['X']

vsDmel_bam_file = 'barcoded.rmDups.sorted.bam' 
longread_nSorted_bam_file = 'longreads.nSorted.bam'

## SUBSMPL RUN
masking_file = 'teflonReference.fa.out' #repeatmasker output file
bam_file_dirOfSamples = 'bams_subsampled' #Subsampled bams from Teflon
bam_file_basename = 'RAW'
req_rsqE_ovlp = False
##/
### OUTPUTS
output_dir = 'outputs'
###/

### Define globals
global_vars_modules.global_idx_iterator = 0 # global_vars_modules is a dummy-class from my functions.py
ticToc = None

########### SCRIPT START

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
#/

# Import longread rsqEs (define search-space of short-reads by long-read alignments which pass filtering thresholds to alleviate mappings including problematic regions)
if 0:
    print('Trying to import qname_rsq_idx_map!')
    with open('qname_rsq_idx_map.pickle','rb') as f:
        qname_rsq_idx_map = pickle.load(f)
    qname_rsq_idx_map_ROMD = rangeOverlaps_makeMappingDict(qname_rsq_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
    qname_rsqEs_map = {} #longread qname -> rsqE_idxs
    for rsqE in qname_rsq_idx_map.values():
        init(qname_rsqEs_map,rsqE['qname'],set())
        qname_rsqEs_map[rsqE['qname']].add(rsqE['idx'])
        
else:
    print('Running without check against qname_rsq_idx_map!')
    qname_rsq_idx_map = {}
    qname_rsq_idx_map_ROMD = {}
    qname_rsqEs_map = {}
#/

def readOrientation_dupC_check(read1,read2):
    pair_is_dupC = False
    #        vanilla case (read1 is fw)
    #    --read1-> <--read2--
    #    ----ref------>------dupC---------->ref-----
    #                                   --> read1
    #                   <--read2
    #
    #        vanilla case (read1 is rv)
    #    --read2-> <--read1--
    #    ----ref------>------dupC---------->ref-----
    #                                   --> read2
    #                   <--read1
    #
    
    # If mate1 is fw, then expect rv-read to be at a lower rcoords on ref 
    if strand_isfv(read1['strand']) and strand_isrv(read2['strand']):
        if mean(read1['rcoords']) > mean(read2['rcoords']):     pair_is_dupC = True
    
    # If mate1 is rv, then expect rv-read to be at a higher rcoords on ref
    elif strand_isrv(read1['strand']) and strand_isfv(read2['strand']):
        if mean(read1['rcoords']) < mean(read2['rcoords']):     pair_is_dupC = True
        
    return pair_is_dupC

def readOrientation_dupC_check_NEXTERA(read1,read2):
    pair_is_dupC = False
    #        vanilla case (read1 is fw and expected to be to right of read2 in vanilla case)
    #    <-read2-- --read1-->
    #    ----ref------>------dupC---------->ref-----
    #                              <-- read2 
    #                  --> read1
    #
    #        vanilla case (read1 is rv and expected to be left of read2 in vanilla case)
    #    <-read1-- --read2-->
    #    ----ref------>------dupC---------->ref-----
    #                                   <-- read1
    #                   --> read2
    #
    # If mate1 is fw, then expect rv-read to be at a lower rcoords on ref 
    if strand_isfv(read1['strand']) and strand_isrv(read2['strand']):
        if mean(read1['rcoords']) < mean(read2['rcoords']):     pair_is_dupC = True
    
    # If mate1 is rv, then expect rv-read to be at a higher rcoords on ref
    elif strand_isrv(read1['strand']) and strand_isfv(read2['strand']):
        if mean(read1['rcoords']) > mean(read2['rcoords']):     pair_is_dupC = True

    return pair_is_dupC

### Run loop per sample
SAMPLE_CALLS = {}
SAMPLE_READLENS = {}

# Construct samples to analyze
sample_paths = {} #sample_id -> bam_file
for sample in os.listdir(bam_file_dirOfSamples):
    tmp_path = bam_file_dirOfSamples + '/' + sample
    
    if bam_file_basename.find('<ID>') != -1:
        tmp_path2 = tmp_path + '/' + bam_file_basename.replace('<ID>',sample)
    elif bam_file_basename.find('RAW') != -1:
        tmp_path2 = tmp_path
        sample = sample.split('.')[0]
        sample_paths[sample] = tmp_path2
    else:
        tmp_path2 = tmp_path + '/' + bam_file_basename
        
    if os.path.isdir(tmp_path) and os.path.exists(tmp_path2):
        sample_paths[sample] = tmp_path2
#/

# Run loop
for sample,vsDmel_bam_file in sample_paths.items():        
    ## Parse BAM, find pairs which have distance on r and are not masked
    INFO_isizes = []
    INFO_masked_covFracs = []
    INFO_NMs = []
    INFO_alens = []
    INFO_NM_alen_ratios = []
    
    far_pairs = {} # readname -> <pair_entry_data>
    mates_readLens = {} #mate -> maxLen
    bam_fo = pysam.AlignmentFile(vsDmel_bam_file,'rb')
    for entry in bam_fo:
        #if entry.is_proper_pair: continue #skip perfect pairs
        if abs(entry.isize) < 1000: continue
        
        qname = entry.qname
        rname = entry.reference_name
        rcoords = entry.reference_start,entry.reference_end
        
        # Check numBP ovlp to maskings
        ovlps = rangeOverlaps_lookup(rcoords,masking_idx_map_ROMD,100,map_lookup_key=rname)
        maskRanges = []
        for oS,maskIdx,oE in ovlps:
            if oE==oS:  oE += 1 #handle fuckup ovlpRanges
            maskRanges.append([oS,maskIdx,oE])
        ovlps2 = computeOvlpRanges_wIdxs(maskRanges)
        masked_numBP = 0
        for rangee in ovlps2:
            masked_numBP += rangee[-1]-rangee[0]
        masked_covFrac = masked_numBP / entry.qlen
        INFO_masked_covFracs.append(masked_covFrac)
        if masked_covFrac >= 0.7: continue
        #/
        
        cigar = entry.cigar
        mate = None
        if entry.is_read1: mate = '1'
        if entry.is_read2: mate = '2'
        strand = 'fv'
        if entry.is_reverse:    strand = 'rv'
        
        pair_id = entry.qname
        
        NM = entry.get_tag('NM')
        alen = entry.alen
        
        NM_alen_ratio = NM/alen
        
        tmp_save = {'rname':entry.reference_name,'rcoords':rcoords,'mate':mate,'strand':strand,'isize':entry.isize,
                    'masked_covFrac':masked_covFrac,'AS':entry.get_tag('AS'),'NM':NM,'alen':alen,'NM_alen_ratio':NM_alen_ratio}
        
        init(mates_readLens,mate,0)
        if entry.qlen > mates_readLens[mate]:      mates_readLens[mate] = entry.qlen
        
        init(far_pairs,pair_id,[])
        far_pairs[pair_id].append(tmp_save)
        
        INFO_isizes.append(entry.isize)
        INFO_NMs.append(NM)
        INFO_alens.append(alen)
        INFO_NM_alen_ratios.append(NM_alen_ratio)
    bam_fo.close()
    
    if 0:
        filt0 = INFO_isizes
        filt1 = filterArr( filt0, 0, 'threshold' )
        filt2 = filterArr( filt1, 100000, 'cutoff' )
        fig = plt.figure(figsize=(11,11))
        plt.hist( filt2 ,bins=200)
        plt.title('numVals='+str(len(filt2)))
    ##/
    
    ### Call rearrs
    ## Hash pair read-aln locations; Post-traverse found pairs, take only those that have both reads fetched
    INFO_numEntries = []
    read_pairs_filtered = {} #save best mates read_entires
    
    rearr_pairs_readEntries_idxs = {} # pair_id -> arr of idx of read_entries
    rearr_readEntries_idx_map = {}
    for pair_id,read_entries in far_pairs.items():
        if not len(read_entries) >= 2: continue #skip if not even 2 aln entries for pair
    
        # Fetch best AS per mate & save
        mates_best_entries = {}
        for read_entry in read_entries:
            mate = read_entry['mate']
            if not mate in mates_best_entries or mates_best_entries[mate]['AS'] < read_entry['AS']:
                mates_best_entries[mate] = read_entry
                
        if not len(mates_best_entries) == 2: continue #skip if we didnt have both reads in pair
        
        read_pairs_filtered[pair_id] = mates_best_entries
        #/
        # Infer rearr
        if mates_best_entries['1']['rname'] == mates_best_entries['2']['rname']:
            mates_idx_track = [] #keep track of idx assigned to each mate in pair, add post-loop
            for mate_enum,read_entry in mates_best_entries.items():
                read_entry['pair_id'] = pair_id
                assignIdx(read_entry)
                rearr_readEntries_idx_map[read_entry['idx']] = read_entry
                mates_idx_track.append(read_entry['idx'])
            
            # save mates idx at each other
            rearr_readEntries_idx_map[ mates_idx_track[0] ]['mate_idx'] = mates_idx_track[1]
            rearr_readEntries_idx_map[ mates_idx_track[1] ]['mate_idx'] = mates_idx_track[0]
            rearr_pairs_readEntries_idxs[pair_id] = mates_idx_track
            #/
        #/
    rearr_readEntries_idx_map_ROMD = rangeOverlaps_makeMappingDict(rearr_readEntries_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
    
    if 0 and 'plot pair distances':
        INFO_isizes = []
        for pair_id,readEntryIdxs in rearr_pairs_readEntries_idxs.items():
            read_entry = rearr_readEntries_idx_map[ readEntryIdxs[0] ]
            INFO_isizes.append(abs(read_entry['isize']))
            
        filt0 = INFO_isizes
        filt1 = filterArr( filt0, 0, 'threshold' )
        filt2 = filterArr( filt1, 100000, 'cutoff' )
        fig = plt.figure(figsize=(11,11))
        plt.hist( filt2 ,bins=200)
        plt.title('numVals='+str(len(filt2)))
                
    ##/
    ## Chain pairs which have reads that ovlp at both read locations
    
    chain_dist = 3000 #attempt to chain pair read ovlps to oPair reads within this numBP
    pair_oPair_ovlps = set()
    rearrs_byPairs = {} # any_pair -> oPairs supping it
    pairs_used_as_supps = {} # keep track of which pairs were used as supp for pairs which made a call
    INFO_numEntryIdxs = []
    INFO_dupCindicator_stats = {'pairs_tot':0,'pairs_wOpair':0,'pairs_wMultipleOpairs':0,'pairs_wLRconfirm':0,'pairs_wLRconfirm_wOpair':0,'pairs_asSupp':0,
                                'pairs_wMultipleLRconfirm':0,'pairs_wMultipleLRconfirm_wOpair':0,'pairs_woLRconfirm_wMultipleOpair':0,
                                'pairs_woLRconfirm':0,'pairs_wOpair_woLRconfirm':0}
    
    IDE_pairs_woLRconfirm = {} # rname -> arr_of_entries
    IDE_pairs_woLRconfirm2 = {} # tenX readname -> rname -> arr_of_entries
    IDE_calls = {}
    for enum,(pair_id,read_entryIdxs) in enumerate(rearr_pairs_readEntries_idxs.items()):
        if pair_id in pairs_used_as_supps:
            INFO_dupCindicator_stats['pairs_asSupp'] += 1
            continue #skip if already used pair in supp
        
        #
        # Dont attempt to call rearr if too much masking (NOTE: it can still be used to support the call from another pair!)
        mates_passed = 0
        for entryIdx in read_entryIdxs:
            readEntry = rearr_readEntries_idx_map[entryIdx]
            if readEntry['masked_covFrac'] >= 0.2: continue
            # Skip entry if too many mismatches
            if readEntry['mate'] == '1' and (readEntry['alen'] < round(mates_readLens['1']*0.95) or readEntry['NM'] >= round(mates_readLens['1']*0.95)): continue
            if readEntry['mate'] == '2' and (readEntry['alen'] < round(mates_readLens['2']*0.95) or readEntry['NM'] >= round(mates_readLens['2']*0.95)): continue
            #/
            mates_passed += 1
        if not mates_passed >= 2: continue
        #/
        
        ## PRE IDE-analysis
        read1 = read_pairs_filtered[pair_id]['1']
        read2 = read_pairs_filtered[pair_id]['2']
        
        #@ RUN FILTERS
        
        # expect them to be mapped on diff strands
        if read1['strand'] == read2['strand']: continue
    
        # Check if ovlp rsqEs at both mates
        if req_rsqE_ovlp:
            mates_passed = 0
            for entryIdx in read_entryIdxs:
                readEntry = rearr_readEntries_idx_map[entryIdx]
                rname = readEntry['rname']
                rcoords = readEntry['rcoords']
                ovlps = rangeOverlaps_lookup(rcoords,qname_rsq_idx_map_ROMD,1000,map_lookup_key=rname)
                if ovlps:
                    for oS,idx,oE in ovlps:
                        rsqE = qname_rsq_idx_map[idx]
                        if rsqE['masked_covFrac'] >= 0.7: continue
                        if oE-oS >= readEntry['alen']*0.5:
                            mates_passed += 1
                            break
            if not mates_passed >= 2: continue
        #/
        
        rname = read1['rname']
        rspan = [min(read1['rcoords']+read2['rcoords']),max(read1['rcoords']+read2['rcoords'])]
        rspan_len = rspan[-1]-rspan[0]
    
        # Run filters
        if not rname in chroms: continue
        if rspan_len < 5000: continue
        if rspan_len > 50000: continue
    
        # parse masking of inferred dupC (to filter bullshitters)
        ovlps = rangeOverlaps_lookup(rspan,masking_idx_map_ROMD,100,map_lookup_key=rname)
        maskRanges = []
        for oS,maskIdx,oE in ovlps:
            if oE==oS:  oE += 1 #handle fuckup ovlpRanges
            maskRanges.append([oS,maskIdx,oE])
        ovlps2 = computeOvlpRanges_wIdxs(maskRanges)
        masked_numBP = 0
        for rangee in ovlps2:
            masked_numBP += rangee[-1]-rangee[0]
        masked_covFrac = masked_numBP / rspan_len
        if masked_covFrac >= 0.5: continue
        #/
        
        pair_is_dupC = False
        pair_is_dupC = readOrientation_dupC_check(read1,read2)
        if pair_is_dupC:
            INFO_dupCindicator_stats['pairs_tot'] += 1
        else:
            continue
        #@/
        ##/
        
        ## Cluster and check supp with other readpairs
        # check ovlps to oPairReads
        oPair_entryIdxs_ovlp = {} #oPair_id -> pairReadEntryIdx -> data
        for entryIdx in read_entryIdxs:
            readEntry = rearr_readEntries_idx_map[entryIdx]
            rname = readEntry['rname']
            #if rname == '2L': continue
            rstart,rend = readEntry['rcoords']
    
            ovlps = rangeOverlaps_lookup([rstart-chain_dist,rend+chain_dist],rearr_readEntries_idx_map_ROMD,1000,map_lookup_key=rname)
            for oS,oEntryIdx,oE in ovlps:
                oEntry = rearr_readEntries_idx_map[oEntryIdx]
                oPair_id = oEntry['pair_id']
                if oPair_id == pair_id: continue #skip self
                oRead1,oRead2 = read_pairs_filtered[oPair_id]['1'],read_pairs_filtered[oPair_id]['2']
                
                # Check if expected dupC orientation
                pair_is_dupC = False
                pair_is_dupC = readOrientation_dupC_check(oRead1,oRead2)
                if not pair_is_dupC: continue
            
                # Check if too much mismatch in oPair reads aln (use read1 & read2 to relate and handle potential different readlens between datasets)
                if oRead1['AS'] < read1['AS']*0.8: continue
                if oRead2['AS'] < read2['AS']*0.8: continue
                #/
                # check so call region of oPair is similar to current pair
                oPair_rspan = [min(oRead1['rcoords']+oRead2['rcoords']),max(oRead1['rcoords']+oRead2['rcoords'])]
                if not getRangeOvlp(oPair_rspan,rspan) >= rspan_len*0.8: continue
                #/
                
                init(oPair_entryIdxs_ovlp,oPair_id,{})
                init(oPair_entryIdxs_ovlp[oPair_id],entryIdx,[])
                oPair_entryIdxs_ovlp[oPair_id][entryIdx].append(oEntryIdx)
        #/
        # Filter oPairs that mapped at both pair reads
        oPair_entryIdxs_filt = {}
        for oPair_id,entryIdxs in oPair_entryIdxs_ovlp.items():
            INFO_numEntryIdxs.append(len(entryIdxs))
            if len(entryIdxs) == 2:
                oPair_entryIdxs_filt[oPair_id] = entryIdxs
                
        #/
        # Infer rearr
        if oPair_entryIdxs_filt:
            for oPair_id in oPair_entryIdxs_filt:
                if oPair_id in pairs_used_as_supps: continue #dont add pair as supp to current if already supping something else
                pairs_used_as_supps[oPair_id] = pair_id
                
                init(rearrs_byPairs,pair_id,set())
                rearrs_byPairs[pair_id].add(oPair_id)
            
            for oPair_id in oPair_entryIdxs_filt:
                comp_id = '||'.join(sorted(map(str,[pair_id,oPair_id])))
                pair_oPair_ovlps.add(comp_id)
        #/
        # IDE INFO
        if oPair_entryIdxs_filt:
            INFO_dupCindicator_stats['pairs_wOpair'] += 1
            if len(oPair_entryIdxs_filt) >= 5:
                INFO_dupCindicator_stats['pairs_wMultipleOpairs'] += 1
        #/
        
        # Save call
        num_oPairs = len(oPair_entryIdxs_filt)
        tmp_rcoords = [min(read1['rcoords']+read2['rcoords']),max(read1['rcoords']+read2['rcoords'])]
        init(IDE_calls,read1['rname'],[])
        IDE_calls[read1['rname']].append([read1['pair_id'],tmp_rcoords,num_oPairs,round(masked_covFrac,4)])
        ##/
        
        ## IDE-INFO: Check if ovlp with DUP called by longreads
        # Find rsqE ovlps at both mates, sorted by longread qname
        mates_qnameRsqEIdxs_ovlp = {} # mate_idx -> longread_qname
        for mate_enum,entryIdx in enumerate(read_entryIdxs):
            readEntry = rearr_readEntries_idx_map[entryIdx]
            rname = readEntry['rname']
            rcoords = readEntry['rcoords']
            ovlps = rangeOverlaps_lookup(rcoords,qname_rsq_idx_map_ROMD,1000,map_lookup_key=rname)
            for oS,idx,oE in ovlps:
                if oE-oS >= readEntry['alen']*0.5:
                    rsqE = qname_rsq_idx_map[idx]
                    longread_qname = rsqE['qname']
                    init(mates_qnameRsqEIdxs_ovlp,mate_enum,{})
                    init(mates_qnameRsqEIdxs_ovlp[mate_enum],longread_qname,set())
                    mates_qnameRsqEIdxs_ovlp[mate_enum][longread_qname].add(idx)
        #/
        mate_longRead_dupSupp = {} # mate_enum -> longread_qname -> rsqE_idxs
        if mates_qnameRsqEIdxs_ovlp:
            # Traverse longread_qnames which ovlp'ed both mates & check which mate longread_qname can confirm DUP
            for longread_qname in mates_qnameRsqEIdxs_ovlp[0].keys() & mates_qnameRsqEIdxs_ovlp[1].keys():
                # Grab rsqEs
                rsqEs = []
                for rsqE_idx in qname_rsqEs_map[longread_qname]:
                    rsqEs.append(qname_rsq_idx_map[rsqE_idx])
                #/
                # Make ovlpRanges with mates and check for DUP-call support by the longread
                rRanges = []
                rRanges.append([read1['rcoords'][0],'mate1',read1['rcoords'][-1]])
                rRanges.append([read2['rcoords'][0],'mate2',read2['rcoords'][-1]])
                for rsqE in rsqEs:
                    rname = rsqE['rname']
                    rstart,rend = rsqE['rcoords']
                    rRanges.append([rstart,rsqE['idx'],rend])
                
                ovlps = computeOvlpRanges_wIdxs(rRanges)
                mates_numBP_DUPconfirmed = {} # mate -> numBP_with_dupConfirm
                for rangee in ovlps:
                    if 'IDE-bugcheck' and  'mate1' in rangee[1] and 'mate2' in rangee[1]: sys.exit('IDE BUG!!!!ahduiasduhasuid41')
                    # Find rangee mate & rsqE idxs
                    rangee_mate = None
                    rangee_rsqE_idxs = []
                    for entry in rangee[1]:
                        if entry in ('mate1','mate2',):
                            rangee_mate = entry
                        else:
                            rangee_rsqE_idxs.append(entry)
                    #/
                    
                    if not rangee_mate: continue #skip if not a rangee where pair mate exist (it was only longreads)
                    if not len(rangee_rsqE_idxs) >= 2: continue # skip if not at least x2 ovlp'ign rsqE at rangee (expected if DUP)
                    
                    # Check so LR alns are not differential fw/rv only (prevent bullshit inverted reads)
                    rangee_rsqE_orientations = {}
                    for rsqE_idx in rangee_rsqE_idxs:
                        rsqE = qname_rsq_idx_map[rsqE_idx]
                        strand = rsqE['strand']
                        init(rangee_rsqE_orientations,strand,0)
                        rangee_rsqE_orientations[strand] += 1
                    
                    if not max(rangee_rsqE_orientations.values()) >= 2: continue #skip if we did not have at least x1 occ of fw+fw or rv+rv aln from longread
                    #/
                    
                    init(mates_numBP_DUPconfirmed,rangee_mate,0)
                    mates_numBP_DUPconfirmed[rangee_mate] += rangee[-1]-rangee[0]
                #/
                # Check if save to outer
                for mate,numBP_dupConfirmed in mates_numBP_DUPconfirmed.items():
                    if numBP_dupConfirmed >= min(read1['alen'],read2['alen'])*0.5:
                        init(mate_longRead_dupSupp,mate,set())
                        mate_longRead_dupSupp[mate].add(longread_qname)
            #/
            # Check longread DUP confirm
            if mate_longRead_dupSupp:
                bothMates_numLRs_confirm = 0
                if 'mate1' in mate_longRead_dupSupp and 'mate2' in mate_longRead_dupSupp:
                    bothMates_numLRs_confirm = len(mate_longRead_dupSupp['mate1'].intersection(mate_longRead_dupSupp['mate2']))
                    
                mate_numLRs_confirm = {}
                for mate,LRs in mate_longRead_dupSupp.items():
                    mate_numLRs_confirm[mate] = len(LRs)
                
                ### SAVE INFO-STATS
                if mate_numLRs_confirm and max(mate_numLRs_confirm.values()) >= 1:
                    INFO_dupCindicator_stats['pairs_wLRconfirm'] += 1
                    
                    if oPair_entryIdxs_filt:
                        INFO_dupCindicator_stats['pairs_wLRconfirm_wOpair'] += 1
                        
                    if mate_numLRs_confirm and max(mate_numLRs_confirm.values()) >= 3:
                        INFO_dupCindicator_stats['pairs_wMultipleLRconfirm'] += 1
                        
                        if oPair_entryIdxs_filt:
                            INFO_dupCindicator_stats['pairs_wMultipleLRconfirm_wOpair'] += 1
                ###/
        #/
        
        ### INFO-STATS
        if not mate_longRead_dupSupp:
            INFO_dupCindicator_stats['pairs_woLRconfirm'] += 1
            
            num_oPairs = len(oPair_entryIdxs_filt)
            tmp_rcoords = [min(read1['rcoords']+read2['rcoords']),max(read1['rcoords']+read2['rcoords'])]
            
            init(IDE_pairs_woLRconfirm,read1['rname'],[])
            IDE_pairs_woLRconfirm[read1['rname']].append([read1['pair_id'],tmp_rcoords,num_oPairs,round(masked_covFrac,4)])
            
            init(IDE_pairs_woLRconfirm2,read1['pair_id'],[])
            IDE_pairs_woLRconfirm2[read1['pair_id']].append([read1['rname'],tmp_rcoords,num_oPairs,round(masked_covFrac,4)])
            
        if oPair_entryIdxs_filt and not mate_longRead_dupSupp:
            INFO_dupCindicator_stats['pairs_wOpair_woLRconfirm'] += 1
            
            if len(oPair_entryIdxs_filt) >= 4 and not mate_longRead_dupSupp:
                INFO_dupCindicator_stats['pairs_woLRconfirm_wMultipleOpair'] += 1
        ###/
        ##/
        
        if enum%10000 == 0:
            print(len(pair_oPair_ovlps),enum,'/',len(rearr_pairs_readEntries_idxs))
            
    for rname,entries in IDE_calls.items():
        IDE_pairs_woLRconfirm[rname] = sorted(entries,key=lambda x: x[1][0])
    for rname,entries in IDE_pairs_woLRconfirm.items():
        IDE_pairs_woLRconfirm[rname] = sorted(entries,key=lambda x: x[1][0])
        
    ## Save at sample
    SAMPLE_CALLS[sample] = IDE_calls
    SAMPLE_READLENS[sample] = mates_readLens
    ##/
#/

if 0 and 'make phylo from dupC calls':
    # 1. Parse rnames + coords, 2. Init bins, 3. find dupC ovlps to bins, 4. remove binns without dupCs
    rname_coords = {}
    for sample,IDE_calls in SAMPLE_CALLS.items():
        for rname,entries in IDE_calls.items():
            for entry in entries:
                init(rname_coords,rname,[])
                rname_coords[rname] += entry[1]
    bin_size = 1000
    bin_idx_map = {}
    for rname,rcoords in rname_coords.items():
        for i in range(0,max(rcoords),bin_size):
            bstart = i
            bend = i+bin_size
            binn = {'rcoords':[bstart,bend],'rname':rname,'sample_dupCs':{}}
            assignIdx(binn)
            bin_idx_map[binn['idx']] = binn
    bin_idx_map_ROMD = rangeOverlaps_makeMappingDict(bin_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
    for sample,IDE_calls in SAMPLE_CALLS.items():
        for rname,entries in IDE_calls.items():
            for entry in entries:
                rcoords = entry[1]
                ovlps = rangeOverlaps_lookup(rcoords,bin_idx_map_ROMD,1000,map_lookup_key=rname)
                for oS,idx,oE in ovlps:
                    binn = bin_idx_map[idx]
                    init(binn['sample_dupCs'],sample,[])
                    binn['sample_dupCs'][sample].append(entry)
    rm_idxs = []
    for binn_idx,binn in bin_idx_map.items():
        if not binn['sample_dupCs']:    rm_idxs.append(binn_idx)
    for rm_idx in rm_idxs:
        del bin_idx_map[rm_idx]
    bin_idx_map_ROMD = rangeOverlaps_makeMappingDict(bin_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
    #/
    # Make phylip file
    sample_phylips = {} # sample -> arr_of_booleans, appending per observation
    readsupp_CO = 2
    for binn_idx,binn in sorted(bin_idx_map.items(),key=lambda x: (x[1]['rname'],x[1]['rcoords'][0])):
        for sample in SAMPLE_CALLS:
            saveState = 0
            if sample in binn['sample_dupCs']:
                for entry in binn['sample_dupCs'][sample]:
                    if entry[2] >= readsupp_CO:
                        saveState = 1
                        break #only use first dupC call at bin
            
            init(sample_phylips,sample,[])
            sample_phylips[sample].append(saveState)
    
    if 1 and 'dump to file?':
        nf = open(output_dir+'/'+'dupCs.req'+str(readsupp_CO)+'.phylip','w')
        # Header
        nf.write(str(len(sample_phylips))+' '+str(len(list(sample_phylips.values())[0]))+'\n')
        # Data
        for sample,obs in sample_phylips.items():
            nf.write(sample+' '+''.join(map(str,obs))+'\n')
        nf.close()
        
    if 1 and 'dump num found?':
        with open(output_dir+'/'+'dupCs_per_sample.tsv','w') as nf:
            for sample,obs in sample_phylips.items():
                nf.write(sample+'\t'+str(obs.count(1))+'\n')
    #/
    
if 0 and 'IDE-dump pairs which had no reference supp':
    bam_fo = pysam.AlignmentFile(vsDmel_bam_file,'rb')
    bam_nf = pysam.AlignmentFile('noLRconfirm.bam','wb',template=bam_fo)
    bam_nf2 = pysam.AlignmentFile('noLRconfirm_suppPairs.bam','wb',template=bam_fo)
    for entry in bam_fo:
        #if entry.is_proper_pair: continue #skip perfect pairs
        if abs(entry.isize) < 1000: continue
        qname = entry.qname
        rname = entry.reference_name
        rcoords = entry.reference_start,entry.reference_end        
        
        if qname in IDE_pairs_woLRconfirm2:
            if IDE_pairs_woLRconfirm2[qname][0][0] == rname and getRangeOvlp(IDE_pairs_woLRconfirm2[qname][0][1],rcoords) >= 50:
                bam_nf.write(entry)
        elif qname in pairs_used_as_supps and pairs_used_as_supps[qname] in IDE_pairs_woLRconfirm2:
            bam_nf2.write(entry)
    
    bam_fo.close()
    bam_nf.close()
    bam_nf2.close()
###/
    
sys.exit('PRE-longread')

### Check frequency of longreads DUPclust in hetero, including shorter LRs as well
#   (traverse nSorted BAM)
IDE_DUPC_LONGREADS = {}

read_length_upper = 100000*999
read_length_lower = 2500
reads_parsed_PBfix = set()
refLens = importBAMlengths(longread_nSorted_bam_file)
bam_fo = pysam.AlignmentFile(longread_nSorted_bam_file,'rb')
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
              
                
            ## Traverse alns:
            # 1. make ovlpRanges
            # 2. parse ovlpRanges with DUP-indi
            # 3. Assign maskings to dupC-indis
            # 4. Filter1 & find dupC-indi rcoords (e.g. recoil rcoords)
            
            # 1.
            rRanges = {}
            for enum,entry in enumerate(qname_bam_entries):
                # Parse
                entry['enum'] = enum #save enum
                rname = entry['rname']
                rstart,rend = entry['rcoords']
                init(rRanges,rname,[])
                rRanges[rname].append([rstart,enum,rend])
            #/
            # 2.
            rDupCindis = {}
            for rname,rangees in rRanges.items():
                ovlps = computeOvlpRanges_wIdxs(rangees)
                for rangee in ovlps:
                    if len(rangee[1]) < 2: continue
                    # check so rangee has at least x2 alns of same direction (expected in dupC + prevent bullshit reads)
                    rangee_strands = {}
                    for entry_enum in rangee[1]:
                        entry = qname_bam_entries[entry_enum]
                        strand = entry['strand']
                        init(rangee_strands,strand,[])
                        rangee_strands[strand].append(entry)
                    
                    for strand,entries in rangee_strands.items():
                        if len(entries) < 2: continue # skip if not at least x2 entries at strand
                        
                        entries_enums = []
                        for entry in entries:
                            entries_enums.append(entry['enum'])
                            
                        # save rangee (chain with previous if possible)
                        init(rDupCindis,rname,[])
                        
                        if rDupCindis[rname] and rDupCindis[rname][-1]['rcoords'][-1] == rangee[0]:
                            rDupCindis[rname][-1]['rcoords'][-1] = rangee[-1]
                            rDupCindis[rname][-1]['bam_entries_enums'].update(set(entries_enums))
                        else:
                            rDupCindis[rname].append( {'rname':rname,'rcoords':[rangee[0],rangee[-1]],'bam_entries_enums':set(entries_enums) } )
            #/
            # 3.
            for rname,dupCindis in rDupCindis.items():
                for dupCindi in dupCindis:
                    rcoords = dupCindi['rcoords']
                    # parse masking of inferred dupC (to filter bullshitters)
                    ovlps = rangeOverlaps_lookup(rcoords,masking_idx_map_ROMD,100,map_lookup_key=rname)
                    maskRanges = []
                    for oS,maskIdx,oE in ovlps:
                        if oE==oS:  oE += 1 #handle fuckup ovlpRanges
                        maskRanges.append([oS,maskIdx,oE])
                    ovlps2 = computeOvlpRanges_wIdxs(maskRanges)
                    masked_numBP = 0
                    for rangee in ovlps2:
                        masked_numBP += rangee[-1]-rangee[0]
                    masked_covFrac = masked_numBP / (rcoords[-1]-rcoords[0])
                    dupCindi['masked_covFrac'] = masked_covFrac
                    #/
            #/
            # 4. 
            rDupCindis_filt1 = {}
            for rname,dupCindis in rDupCindis.items():
                for dupCindi in dupCindis:
                    rcoords = dupCindi['rcoords']
                    rlen = rcoords[-1]-rcoords[0]
                    dupCindi['rlen'] = rlen
                    if dupCindi['masked_covFrac'] >= 0.7: continue # Skip if too much masking
                    if rlen < 500: continue # skip if overlapping region at recoil is too short
                    
                    # find qspan of bam_entries that give dupC indi. Then parse all bam_entries within qspan
                    qcoordss = []
                    for entry_enum in dupCindi['bam_entries_enums']:
                        entry = qname_bam_entries[entry_enum]
                        qcoordss += entry['qcoords']
                    qspan = [min(qcoordss),max(qcoordss)]
                    
                    entries_in_span = []
                    for entry in qname_bam_entries:
                        entry_rame = entry['rname']
                        entry_rcoords = entry['rcoords']
                        entry_qcoords = entry['qcoords']
                        if entry_rame == rname and getRangeOvlp(entry_qcoords,qspan) >= (entry_qcoords[-1]-entry_qcoords[0]):
                            # assign masking to entries at this stage
                            ovlps = rangeOverlaps_lookup(entry_rcoords,masking_idx_map_ROMD,100,map_lookup_key=entry_rame)
                            maskRanges = []
                            for oS,maskIdx,oE in ovlps:
                                if oE==oS:  oE += 1 #handle fuckup ovlpRanges
                                maskRanges.append([oS,maskIdx,oE])
                            ovlps2 = computeOvlpRanges_wIdxs(maskRanges)
                            masked_numBP = 0
                            for rangee in ovlps2:
                                masked_numBP += rangee[-1]-rangee[0]
                            masked_covFrac = masked_numBP / (entry_rcoords[-1]-entry_rcoords[0])
                            entry['masked_covFrac'] = masked_covFrac
                            #/
                            if masked_covFrac < 0.7:
                                entries_in_span.append(entry)
                    #/
                    # Traverse entries in span to calc rSpan and infer recoil rcoords:
                    # Sort by qcoords and select entries adjacent to dupC-indi bam_entry (depending on strand)
                    # FW-case (expect aln2 at qcoords to have rstart2 < rend1[or rend3])
                    #  --------1------->[--opt3-->]
                    #             -------2---------->
                    # RV-case (expect aln2 at qcoords to have rend2[or rend3] > rstart1)
                    #   <------2--------------
                    #        [<--opt3--]<---1---------
                    
                    entries_in_span.sort(key=lambda x: x['qcoords'])
                    
                    dupCindi_entries_strandSorted = {}
                    for entry_enum in dupCindi['bam_entries_enums']:
                        entry = qname_bam_entries[entry_enum]
                        init(dupCindi_entries_strandSorted,entry['strand'],[])
                        dupCindi_entries_strandSorted[entry['strand']].append(entry)
                    for strand,entries in dupCindi_entries_strandSorted.items():
                        entries.sort(key=lambda x: x['qcoords'])
                    
                    # traverse dupC-indi spans in-order and pairwise, then grab the other recoil rcoord from "entries_in_span"
                    tmp_dupCs = []
                    for strand,entries in dupCindi_entries_strandSorted.items():
                        for idx,entry in enumerate(entries):
                            if idx == 0: continue #backwards checking, see schemata above
                            
                            # find corresponding recoil entry
                            for tmp_idx,tmp_entry in enumerate(entries_in_span):
                                if tmp_entry['enum'] == entry['enum']:
                                    recoil_entry = entries_in_span[tmp_idx-1] #grab entry which is previous to "aln2" in schemata above
                                    # Check if rnames are the same and within calling-distance (in IDE it was 50kb, e.g. we only attempt to call dupClusts with maxRspan of 50kb)
                                    if recoil_entry['rname'] == entry['rname'] and getRangeOvlp(recoil_entry['rcoords'],entry['rcoords']) > -50000:
                                        # Check so rend of recoil alignment is in correct relation to dupC-indi-aln (see schemata; expect recoil rend to be greater than aln2 rstart)
                                        if recoil_entry['rcoords'][-1] > entry['rcoords'][0]: 
                                            tmp_dupCs.append([entry,recoil_entry])
                                        
                    if tmp_dupCs:
                        IDE_DUPC_LONGREADS[rn_PBfix] = tmp_dupCs
            #/
            

        # Reset qname_bam_entries
        qname_bam_entries = []
    #/
    
    # Check if we was EOF, then break
    if type(bam_entry) == str and bam_entry == 'EOF':      break
    
    # Accumulate bam entries of read
    prev_qname = qname
    qname_bam_entries.append(parseBAMentry(bam_entry))

bam_fo.close()
print('\tDone!')


if 0 and 'IDE-bamDUMP dupC-indi reads':
    bam_fo = pysam.AlignmentFile(longread_nSorted_bam_file,'rb')
    bam_nf = pysam.AlignmentFile('LRs_dupCindi.bam','wb',template=bam_fo)
    for entry in bam_fo:
        qname = entry.qname
        rname = entry.reference_name
        rcoords = entry.reference_start,entry.reference_end
        
        if qname in IDE_DUPC_LONGREADS:
            # Skip if aln is repeatMasked
            ovlps = rangeOverlaps_lookup(rcoords,masking_idx_map_ROMD,100,map_lookup_key=rname)
            maskRanges = []
            for oS,maskIdx,oE in ovlps:
                if oE==oS:  oE += 1 #handle fuckup ovlpRanges
                maskRanges.append([oS,maskIdx,oE])
            ovlps2 = computeOvlpRanges_wIdxs(maskRanges)
            masked_numBP = 0
            for rangee in ovlps2:
                masked_numBP += rangee[-1]-rangee[0]
            masked_covFrac = masked_numBP / (rcoords[-1]-rcoords[0])
            if masked_covFrac >= 0.7: continue
            #/
            
            bam_nf.write(entry)
    
    bam_fo.close()
    bam_nf.close()
###/
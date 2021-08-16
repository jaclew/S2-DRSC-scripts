import os
import matplotlib.pyplot as plt
import numpy as np
from statistics import mean,median,stdev
import pickle
import vcf as pyvcf
import sys

from functions import *

### INPUTS
chroms = ['2L','2R','3L','3R','4','X']
autosomes = ['2L','2R','3L','3R']
X = ['X']

GTF_file = 'dmel-all-chromosome-r6.12.gtf' #gene annotations
masking_file = 'RepeatMasked/dmel-all-chromosome-r6.12.fasta.out' #repeatmasker output file

vcf_file = 'MNS15.genotype.cluster.s10.vcf' #sniffles output variant calls
###/

### Define globals
global_vars_modules.global_idx_iterator = 0 # global_vars_modules is a dummy-class from my functions.py
ticToc = None

######## SCRIPT START
### Import repma
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
###/

### Import GTFs
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
###/
  
### Import Sniffles VCF
SV_qual_thresh = 0 #arbitrarily picked
INFO_SV_types = {}
INFO_filter_types = {}
IDE_AOs = {}
IDE_ROs = {}
IDE_AFs = {}
SV_idx_map = {}
vcf_fo = pyvcf.Reader(filename=vcf_file)
print('Importing SNIFFLES VCF...')
for entry in vcf_fo:
    # Format of VCF:
    # chrom1,pos1 -> [chrom2+pos2,...,chromN+posN]
    
    #if entry.INFO['SVTYPE'] == 'DUP': sys.exit()
    SV_id = entry.ID
    chrom1 = entry.CHROM
    pos1 = entry.affected_start
    reads = entry.INFO['RNAMES']
    qual = entry.QUAL
    # Parse type of SV
    typee = entry.INFO['SVTYPE']
    
    if typee == 'DEL/INV':      typee = 'INVDEL'
    
    for filter_ in entry.FILTER:
        init(INFO_filter_types,filter_,{})
        init(INFO_filter_types[filter_],typee,0)
        INFO_filter_types[filter_][typee] += 1
    
    if qual != None and qual < SV_qual_thresh: continue
    #if 'not_fully_covered' in entry.FILTER: continue #skip bullshitters
    
    reads_PBfix = set()
    SV_num_NP_reads = 0
    for read in reads:
        if read[0] == 'm':
            reads_PBfix.add(getPacbioMainRead(read))
        else:
            reads_PBfix.add(read)
            SV_num_NP_reads += 1
    
    # NOTE: SVIM only reports one entry in the ALT
    if len(entry.ALT) >= 2: sys.exit('IDE, we thought we always had 1 alt...')
    for idx,ALT in enumerate(entry.ALT):
        # If SV is INS, then parse the seq
        seq = ''
        if typee == 'INS':
            seq = ALT.sequence
        
        ## INFO
        init(INFO_SV_types,typee,0)
        INFO_SV_types[typee]+=1
        ##/INFO
        #/
        
        ## Parse variant/ref supps
        sample_data = entry.samples[0]
        AO = sample_data.data.DV
        RO = sample_data.data.DR
        AF = entry.INFO['AF'][idx]
            
        LERI1 = 1 # 1 means alns extend to right, -1 means extend to left
        LERI2 = 1 
        if entry.INFO['STRANDS'][idx][0] == '+':    LERI1 = -1
        if entry.INFO['STRANDS'][idx][1] == '+':    LERI2 = -1
        ##/
        
        ## INFO
        if not typee == 'BND':
            init(IDE_AOs,chrom1,[])
            IDE_AOs[chrom1].append(AO)
            init(IDE_ROs,chrom1,[])
            IDE_ROs[chrom1].append(RO)
            init(IDE_AFs,chrom1,[])
            IDE_AFs[chrom1].append(AF)
        ##/INFO
        
        # Try parsing chrom & end coord of SV type
        chrom2 = None
        if 'CHR2' in entry.INFO:        chrom2= entry.INFO['CHR2']
        elif 'chr' in ALT.__dict__:     chrom2= ALT.chr
        
        pos2 = None
        pos2 = entry.affected_end
        
        if 'END' in entry.INFO:
            pos2 = entry.INFO['END']
        #/
        
        # Try parsing len of SV
        SV_len = None
        if 'SVLEN' in entry.INFO:       SV_len = abs(entry.INFO['SVLEN'])
        else:
            # Len is not reported for INV. calc it myself
            if typee == 'INV':
                if not chrom1 == chrom2: sys.exit('IDE: INV had diff chroms!!!')
                SV_len = pos2-pos1

        refseqlen = len(entry.REF)
        
        if not type(pos2) == int:
            print('IDE: Wasnt able to parse end coord?! handle')
            sys.exit()
        #/
        
        # Save
        if typee == 'BND':
            rcoords = [pos1,pos1]
        else:
            rcoords = [pos1,pos2]
        
        SV_entry = {'type':typee,'rname':chrom1,'rcoords':rcoords,
                    'len':SV_len,'refseqlen':refseqlen,
                    'LERI1':LERI1,'LERI2':LERI2,
                    'reads':HOLDER(reads),'reads_PBfix':HOLDER(reads_PBfix),
                    'qual':qual,'AF':AF,'RO':RO,'AO':AO,'reads_numNP':SV_num_NP_reads,
                    'seq':HOLDER(seq),'reads_num':len(reads),'name':SV_id,'flags':set()}
        
        if typee == 'BND':
            SV_entry['rname2'] = chrom2
            SV_entry['rcoord2'] = ALT.pos #stored here for sniffles
        
        # Try parse SVIM-specific stuff
        if 'COPYNUMBER_EST' in entry.INFO:
            SV_entry['CN_EST'] = entry.INFO['COPYNUMBER_EST']
            
        SV_entry['filter_flags'] = entry.FILTER
        
        if 'SOURCECHROM' in entry.INFO:
            SV_entry['sourceChrom'] = entry.INFO['SOURCECHROM'][0]
            sourcePos1 = int(entry.INFO['SOURCESTART'][0])
            sourcePos2 = int(entry.INFO['SOURCEEND'][0])
            # Tweak pos2 if they are identical
            if sourcePos1 == sourcePos2:
                sourcePos2 += 1
                
            SV_entry['sourceRcoords'] = [ sourcePos1,sourcePos2 ]
            SV_entry['sourceChrom1'] = entry.INFO['SOURCECHROM'][0]
            SV_entry['sourceChrom2'] = entry.INFO['SOURCECHROM'][0]
            SV_entry['sourcePos1'] = sourcePos1
            SV_entry['sourcePos2'] = sourcePos2
        #/
        assignIdx(SV_entry)
        SV_idx_map[SV_entry['idx']] = SV_entry

if 0:
    filt0 = IDE_AFs['2L']
    filt1 = filterArr( filt0, 0, 'threshold' )
    filt2 = filterArr( filt1, 0.99, 'cutoff' )
    fig = plt.figure(figsize=(11,11))
    plt.hist( filt2 ,bins=400)
    plt.title('INFO numVals='+str(len(filt2)))
###/

### INS check
if 0:
    INSrcoords_lens = []
    SVlens = []
    refLens = []
    refCalls = []
    for SV_idx,SV in SV_idx_map.items():
        #SVlens.append(SV['len'])
        if not SV['type'] == 'INS': continue
        rcoordslen = SV['rcoords'][-1]-SV['rcoords'][0]
        INSrcoords_lens.append(rcoordslen)
        
        SVlen = SV['len']
        refLen = SV['refseqlen']
        
        refLens.append(refLen)
        
        if min(SVlen,refLen)/max(SVlen,refLen) >= 0.9:
            refCalls.append(1)
        else:
            refCalls.append(0)
    
    filt0 = refLens
    filt1 = filterArr( filt0, 0, 'threshold' )
    filt2 = filterArr( filt1, 10000, 'cutoff' )
    fig = plt.figure(figsize=(11,11))
    plt.hist( filt2 ,bins=400)
    plt.title('INFO numVals='+str(len(filt2)))
###/
    
### Output for PyGenomeTracks
if 0 and 'output for pygenometracks?':
    output_dir_pygenometracks = 'sniffles_calls'
    mkdir(output_dir_pygenometracks)
    
    rsqELACP_plotClassis = {'BND','INS','DUP','DEL','INV','INVDEL','INVDUP'}
    nf_handles = {}
    for classi in rsqELACP_plotClassis:
        nf_handles[classi] = {'sameChrom':open(output_dir_pygenometracks+'/'+classi+'.sameChrom.arcs','w'),
                              'diffChrom':open(output_dir_pygenometracks+'/'+classi+'.diffChrom.arcs','w')}
        
    for SV_idx,SV in SV_idx_map.items():
        if SV['type'] != 'BND' and SV['len'] < 5000: continue
        
        nf_classi = nf_handles[SV['type']] #assign NF handle
        
        # write for BNDs
        if 'rname2' in SV:
            # diffchroms
            if not SV['rname2'] == SV['rname']:
                nf = nf_classi['diffChrom']
                
                writeArr = [SV['rname'],SV['rcoords'][0]-1000,SV['rcoords'][0]-1000,SV['rname'],SV['rcoords'][0]+1000,SV['rcoords'][0]+1000,1]
                nf.write('\t'.join(map(str,writeArr))+'\n')
                writeArr = [SV['rname2'],SV['rcoord2']-1000,SV['rcoord2']-1000,SV['rname2'],SV['rcoord2']+1000,SV['rcoord2']+1000,1]
                nf.write('\t'.join(map(str,writeArr))+'\n')
            # samechroms (never occurrs...)
            else:
                nf = nf_classi['sameChrom']
                writeArr = [SV['rname'],SV['rcoords'][0],SV['rcoords'][0]+1,SV['rname2'],SV['rcoord2'],SV['rcoord2']+1,1]
                nf.write('\t'.join(map(str,writeArr))+'\n')
                
        else: #write for sameChrom
            nf = nf_classi['sameChrom']
            writeArr = [SV['rname'],SV['rcoords'][0],SV['rcoords'][0]+1,SV['rname'],SV['rcoords'][-1],SV['rcoords'][-1]+1,1]
            nf.write('\t'.join(map(str,writeArr))+'\n')
        #/
        
    for classi,chromStates_nfs in nf_handles.items():
        for chromState,nf in chromStates_nfs.items():
            nf.close()
###/
            
### Output for comparison of used reads vs. our algorithm
if 0 and 'output for comparison vs. our algorithm':
    tmp_readsUsedInSniffles = {}
    for SV in SV_idx_map.values():
        if SV['type'] != 'BND' and (SV['len'] < 1000 or SV['len'] > 500000): continue
        for read in SV['reads']():
            init(tmp_readsUsedInSniffles,read,{})
            init(tmp_readsUsedInSniffles[read],SV['type'],set())
            tmp_readsUsedInSniffles[read][SV['type']].add(SV_idx)
    with open('qnames.sniffles.pickle','wb') as nf:
        pickle.dump(tmp_readsUsedInSniffles,nf)
###/

### Check breakpoint repeat content
INFO_repma_stats = {}
BP_dist_scan = 1000
for SV_idx,SV in SV_idx_map.items():
    if SV['type'] == 'BND':
        rcoord1 = SV['rcoords'][0]
        rcoord2 = SV['rcoord2']
        scanCoords1 = [rcoord1,rcoord1+(BP_dist_scan*SV['LERI1'])]
        scanCoords2 = [rcoord2,rcoord2+(BP_dist_scan*SV['LERI2'])]
        
        # Check masking
        
        #/
###/

### Grab genes subject to dup+del
rname_DUPDELranges = {}
for SV_idx,SV in SV_idx_map.items():
    if SV['type'] in ('DUP','DEL',):
        if SV['len'] >= 500000: continue #skip if too long call
        rname = SV['rname']
        rstart,rend = SV['rcoords']
        
        init(rname_DUPDELranges,rname,[])
        rname_DUPDELranges[rname].append([rstart,SV_idx,rend])

genes_DUPDEL = {} #fbgn -> SV_idxs
for rname,rangees in rname_DUPDELranges.items():
    ovlps = computeOvlpRanges_wIdxs(rangees)
    
    for rangee in ovlps:
        rangee_len = rangee[-1] - rangee[0]
        if rangee_len < 1000: continue #skip small overlaps
        rangee_SVtypes = set()
        for SV_idx in rangee[1]:
            rangee_SVtypes.add(SV_idx_map[SV_idx]['type'])
        
        if len(rangee_SVtypes.intersection(set(['DUP','DEL']))) == 2:
            # check if ovlp full gene
            for oS,gtf_idx,oE in rangeOverlaps_lookup([rangee[0],rangee[-1]],GTF_idx_map_ROMD,1000,map_lookup_key=rname):
                GTF = GTF_idx_map[gtf_idx]
                if GTF['type'] == 'gene':
                    if (oE-oS) >= GTF['len']:
                        fbgn = GTF['gene_id']
                        init(genes_DUPDEL,fbgn,[])
                        genes_DUPDEL[fbgn].append(rangee)
            #/
if 0 and 'dump gene list':
    with open('DUPDEL.list','w') as nf:
        nf.write('\n'.join(list(genes_DUPDEL)))
###/
        
### MS 20210601: Find reads to plot path
## Import read lengths
read_lens = {}
for read,seq in importReadSeqsFasta('longreads.fa').items():
    read_lens[read] = len(seq)
##/
        
## Compile read->SVs
read_SVs_byType = {}
for SV in SV_idx_map.values():
    if SV['len'] < 1000 or SV['len'] > 500000: continue
    for read in SV['reads']():
        init(read_SVs_byType,read,{})
        init(read_SVs_byType[read],SV['type'],set())
        read_SVs_byType[read][SV['type']].add(SV['idx'])
##/
## Assign reads with location by SV
locs_to_parse = {'Xreg':'X:12,649,792-12,787,997',
                 '3Lban':'3L:575,723-683,018','3Lmito':'3L:3,317,803-3,490,354',
                 '2Lhotspot6.4':'2L:6,411,900-6,459,698'
                 }
loc_read_SVs = {}
for read,SVtypes_idxs in read_SVs_byType.items():
    for SVtype,idxs in SVtypes_idxs.items():
        for idx in idxs:
            SV = SV_idx_map[idx]
            
            for reg,loc in locs_to_parse.items():
                rname,rstart,rend = loc.split(':')[0],int(loc.split(':')[1].split('-')[0].replace(',','')),int(loc.split(':')[1].split('-')[1].replace(',',''))
                if SV['rname'] == rname and getRangeOvlp(SV['rcoords'],[rstart,rend]) >= 0:
                    init(loc_read_SVs,reg,{})
                    init(loc_read_SVs[reg],read,{})
                    init(loc_read_SVs[reg][read],SV['type'],set())
                    loc_read_SVs[reg][read][SV['type']].add(SV['idx'])
##/
### INFO
if 0:
    def output_pygenometracks(output_name,read_SVs):
        output_dir_PRE = 'read_paintings/'+output_name
        output_dir = output_dir_PRE+'/'+'arcs'
        mkdir(output_dir_PRE)
        mkdir(output_dir)
        
        rsqELACP_plotClassis = {'BND','INS','DUP','DEL','INV','INVDEL','INVDUP'}
        nf_handles = {}
        for classi in rsqELACP_plotClassis:
            nf_handles[classi] = {'sameChrom':open(output_dir+'/'+classi+'.sameChrom.arcs','w'),
                                  'diffChrom':open(output_dir+'/'+classi+'.diffChrom.arcs','w')}
            
        for SV_idx,SV in read_SVs.items():
            #if SV['type'] != 'BND' and SV['len'] < 5000: continue
            
            nf_classi = nf_handles[SV['type']] #assign NF handle
            
            # write for BNDs
            if 'rname2' in SV:
                # diffchroms
                if not SV['rname2'] == SV['rname']:
                    nf = nf_classi['diffChrom']
                    
                    writeArr = [SV['rname'],SV['rcoords'][0]-1,SV['rcoords'][0]-1,SV['rname'],SV['rcoords'][0]+1,SV['rcoords'][0]+1,1]
                    nf.write('\t'.join(map(str,writeArr))+'\n')
                    writeArr = [SV['rname2'],SV['rcoord2']-1,SV['rcoord2']-1,SV['rname2'],SV['rcoord2']+1,SV['rcoord2']+1,1]
                    nf.write('\t'.join(map(str,writeArr))+'\n')
                # samechroms (never occurrs...)
                else:
                    nf = nf_classi['sameChrom']
                    writeArr = [SV['rname'],SV['rcoords'][0],SV['rcoords'][0]+1,SV['rname2'],SV['rcoord2'],SV['rcoord2']+1,1]
                    nf.write('\t'.join(map(str,writeArr))+'\n')
                    
            else: #write for sameChrom
                nf = nf_classi['sameChrom']
                writeArr = [SV['rname'],SV['rcoords'][0],SV['rcoords'][0]+1,SV['rname'],SV['rcoords'][-1],SV['rcoords'][-1]+1,1]
                nf.write('\t'.join(map(str,writeArr))+'\n')
            #/
            
        for classi,chromStates_nfs in nf_handles.items():
            for chromState,nf in chromStates_nfs.items():
                nf.close()
                
        return True
    
    ## check distribution of readlengths at 2Lhotspot
    readLens_2Lhotspot = []
    for read in loc_read_SVs['2Lhotspot6.4']:
        if not read in read_lens: continue
        readLens_2Lhotspot.append(read_lens[read])
    plotHist(readLens_2Lhotspot)
    ##/
    ## get longest reads from Xreg+3Lban+3Lmito
    reads_written = set()
    for loc,locdata in loc_read_SVs.items():
        for read,readdata in locdata.items():
            if not read in read_lens: continue
            read_len = read_lens[read]
            if read_len >= 150000:
                if not read in reads_written:
                    read_SVs = {}
                    for SV in SV_idx_map.values():
                        if read in SV['reads']():
                            read_SVs[SV['idx']] = SV
                    output_pygenometracks(loc+'_'+read,read_SVs)
                    reads_written.add(read)
                if 1 and 'IDE':
                    for SV in SV_idx_map.values():
                        if read in SV['reads']():
                            print(loc,read,read_len,SV['rname'],SV['rcoords'],SV['len'],SV['type'],len(SV['reads']()))
    ##/
###/
###/
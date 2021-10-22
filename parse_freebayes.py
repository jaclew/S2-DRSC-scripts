
import os
import matplotlib.pyplot as plt
import numpy as np
from statistics import mean,median,stdev
import pickle
import vcf as pyvcf
from copy import deepcopy
import sys

from functions import *

### INPUTS
chroms = ['2L','2R','3L','3R','4','X']
autosomes = ['2L','2R','3L','3R']
X = ['X']
chroms_to_plot = ('X','2L','2R','3L','3R','4','mitochondrion_genome')

vcf_file = 'freebayes.vcf.gz' # variant calls from freebayes

GTF_file = 'dmel-all-chromosome-r6.12.gtf' # gene annotations
masking_file = 'dmel-all-chromosome-r6.12.fasta.out' # repeatmasker output

#samples_selection = set(['SRR497713','SRR497721','S2'])
samples_selection = None
if samples_selection:
    print('WILL ONLY IMPORT SAMPLES\n'+'\n'.join(samples_selection))

### Define globals
global_vars_modules.global_idx_iterator = 0 # global_vars_modules is a dummy-class from my functions.py
ticToc = None

SV_len_thresh = 0
refLens = {}
SV_idx_map = {}
INFO_numAlts = []
INFO_types = {}
INFO_lens = {}
INFO_AFs = []
INFO_ABs = []
vcf_fo = pyvcf.Reader(filename=vcf_file)
for ctg,ctg_data in vcf_fo.contigs.items():
    refLens[ctg] = ctg_data.length
for entry in vcf_fo:
    chrom1 = entry.CHROM
    if not chrom1 in chroms: continue # skip if not on major chromosome
    pos1 = entry.start
    pos2 = entry.end
    tmp_pos_partners = set()
    
    tmp_SV = {'rname':chrom1,'rcoords':[pos1,pos2],'samples':{}}
    for sample_data in entry.samples:
        sample = sample_data.sample
    
        if samples_selection and not sample in samples_selection: continue
    
        DP = sample_data.data.DP
        AOs = sample_data.data.AO
        RO = sample_data.data.RO
        
        if DP == None:  DP = -1
        if RO == None:  RO = -1
        if AOs == None: AOs = -1
        if type(AOs) != list:    AOs = [AOs] # have on same format always - if multiple variants at pos and not
        
        AFs = []
        for AO in AOs:
            if AO == None: continue
            AFs.append(AO/DP)
        
        tmp_SV['samples'][sample] = {'DP':DP,'AOs':AOs,'AFs':AFs,'RO':RO}
    
    # Skip if not at least X supp at variant in any sample
    skip_variant = True
    for sample,data in tmp_SV['samples'].items():
        for AO in data['AOs']:
            if AO >= 3:
                skip_variant = False
                break
        if not skip_variant:    break
    if skip_variant: continue
    #/

    # Save
    assignIdx(tmp_SV)
    
    SV_idx_map[tmp_SV['idx']] = tmp_SV

#SV_idx_map_ROMD = rangeOverlaps_makeMappingDict(SV_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
###/IMPORT VCF

### INFO
# Plot AFs for sample X
if 0:
    samples_data = {}
    for SV_idx,SV in SV_idx_map.items():
        for sample,data in SV['samples'].items():
            if sample != 'S2': continue
            for key in ('AOs','RO','DP','AFs'):
                value = data[key]
                init(samples_data,sample,{})
                init(samples_data[sample],key,[])
                
                if type(value) != list:
                    samples_data[sample][key].append(value)
                else:
                    samples_data[sample][key] += value
                
                
    filt0 = samples_data['S2']['AFs']
    filt1 = filterArr( filt0, 0.01, 'threshold' )
    filt2 = filterArr( filt1, 30, 'cutoff' )
    fig = plt.figure(figsize=(11,11))
    plt.hist( filt2 ,bins=400)
    plt.title('INFO numVals='+str(len(filt2)))
###/
    
### Make phylo from SNPs in genomic bins
if 0:
    sample_phylips = {}
    for SV_idx,SV in sorted(SV_idx_map.items(),key=lambda x: (x[1]['rname'],x[1]['rcoords'][0])):
        SV_sample_AFmax = {}
        for sample,data in SV['samples'].items():
            init(sample_phylips,sample,[])
    
            # boolean check & save        
            SNP_boolean = 0
            if sum(data['AFs']) >= 0.5:     SNP_boolean = 1
            
            sample_phylips[sample].append(SNP_boolean)
    
    # Dump phylip for tree
    if 0 and 'dump?':
        nf = open('snps.phylip','w')
        # Header
        nf.write(str(len(sample_phylips))+' '+str(len(list(sample_phylips.values())[0]))+'\n')
        # Data
        for sample,obs in sample_phylips.items():
            nf.write(sample+' '+''.join(map(str,obs))+'\n')
        nf.close()
    #/
###/
    
sys.exit()

### S2-related: Filter SNPs which do not have AF = 1 in ALL S2-relative lines "sLines" = schneider-lines
sLines_SV_idx_map = {}
#sLines_IDs = ['S2','SRR11076397','SRR497720','SRR497721','SRR497713','SRR497728','SRR612111'] # 20=Sg4,21=S3,13=S1,97=S2R+,28=Mbn2,111=S2_2
sLines_IDs = ['S2','SRR497721','SRR497713']
#sLines_IDs = ['S2','SRR11076397','SRR497720','SRR497721','SRR497713','SRR497728','SRR612111','SRR497710','SRR497711','SRR497725','SRR497729']
INFO_sLines_stats = {'inAll':0,'inNone':0,'inHetero_any':0}
for SV_idx,SV in SV_idx_map.items():
    # INFO: check if AF is 1 in all sLines, then we want to skip.
    AFs_sum = [] # if used, can miss the few haplotype-specific SNPs. OPPOSITE: can be used to identify potential source SNPs from diploid ancestor
    AFs_max = [] #
    for sample in sLines_IDs:
        AFs_sum.append(sum(SV['samples'][sample]['AFs']))
        AFs_max.append(max(SV['samples'][sample]['AFs']))
    
    if min(AFs_max) >= 1 or min(AFs_sum) >= 1: # skip if minimum AF found across samples is X. Then we cant separate from reference-error
        INFO_sLines_stats['inAll'] += 1
    #/INFO
    
    if 0 and 'try1':
        if max(AFs_sum) == 0: # add to INFO if none of sLines had call at SV
            INFO_sLines_stats['inNone'] += 1
            continue
        
        if max(AFs_sum) <= 0.1: # Skip if this site appears to be heterogeneous
            continue
        #/
        # Check if any has >X AO, then we want to keep. Dont wanna analyze heterogeneous bullshitters
        AOs_max = []
        for sample in sLines_IDs:
            AOs_max.append(max(SV['samples'][sample]['AOs']))
            
        if max(AOs_max) < 4: #skip if no sample has at least X amount of reads in call
            INFO_sLines_stats['inHetero_any'] += 1
            continue 
        #/
    
    AOs_summed = []
    AFs_summed = []
    ROs = []
    isRO = []
    for sample in sLines_IDs:
        RO = SV['samples'][sample]['RO']
        AO_sum = sum(SV['samples'][sample]['AOs'])
        AF_sum = sum(SV['samples'][sample]['AFs'])
        AOs_summed.append(AO_sum)
        AFs_summed.append(AF_sum)
        ROs.append(RO)
        if RO == SV['samples'][sample]['DP']:
            isRO.append(sample)
        
    
    if len(isRO) == 0 or max(AOs_summed) < 4 or max(AFs_summed) < 0.125: continue
    
    # Filter oSamples (for IDE development)
    SV_cp = deepcopy(SV)
    for sample in list(SV_cp['samples']):
        if not sample in sLines_IDs:
            del SV_cp['samples'][sample]
    #/
    
    # Save
    sLines_SV_idx_map[SV_idx] = SV_cp
    #/
    
# Find ovlpRanges of SVs
def assignHaplotype(inp_val,haplotypeE,haplotypeEmarg):
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
        print('wtf?')
    return hapClassi

# Find SVs which does not exist in any line, but exist in another line
IDE_SVs = []
IDE_genes = set()
IDE_genes_perSample = {}
IDE_numAOhas = []
IDE_multiAOhas_genes = {}
IDE_gene_sampleHas = {}
IDE_variant_sampleHas_wGeneData = {}
IDE_sample_variants = {}

hapClassis_perSample = {}
AFs_perSample = {}
AOs_perSample = {}
DPs_perSample = {}
for SV_idx,SV in sLines_SV_idx_map.items():
    RO_only = []
    AO_has = []
    AO_has_AFpass = []
    AF_sums = []
    DPs = []
    for sample,data in SV['samples'].items():
        # Check if sample has DP within range X, else skip. By not adding the DP to arr, we will not even consider it in below IF
        if 1:
            if sample == 'SRR497713' and getRangeOvlp([data['DP']]*2,[8,12]) < 0: continue
            if sample == 'S2' and getRangeOvlp([data['DP']]*2,[8,12]) < 0: continue
            if sample == 'SRR497721' and getRangeOvlp([data['DP']]*2,[7,11]) < 0: continue
        if 0:
            if getRangeOvlp([data['DP']*2],[7,20]) < 0: continue
        #/
        AF_sums.append(sum(data['AFs']))
        DPs.append(data['DP'])
        if data['RO'] >= 4 and data['RO'] == data['DP']:
            RO_only.append(sample)
        elif max(data['AOs']) >= 4:
            AO_has.append(sample)
            if max(data['AFs']) >= 0.125:
                AO_has_AFpass.append(sample)
    
    if len(DPs) < 3: continue #for AF plot 3 samples
    if not SV['rname'] in autosomes: continue #for AF plot 3 samples
    if len(RO_only) >= 1 and len(AO_has_AFpass) >= 1 and SV['rcoords'][-1]-SV['rcoords'][0] <= 1: #and min(DPs) >= 8 and max(DPs) <= 12:
        # Stats
        for sample,data in SV['samples'].items():
            #if data['DP'] < 7 or data['DP'] > 15: continue
            init(AFs_perSample,sample,[])
            init(AOs_perSample,sample,[])
            init(DPs_perSample,sample,[])
            init(hapClassis_perSample,sample,{})
            #AFs_perSample[sample].append(sum(data['AFs']))
            AFs_perSample[sample] += data['AFs']
            AOs_perSample[sample].append(sum(data['AOs']))
            DPs_perSample[sample].append(data['DP'])
            
            for AF in data['AFs']:
                if AF == 1:
                    hapClassi = 'all'
                else:
                    hapClassi = assignHaplotype(AF,0.25,0.125)

                init(hapClassis_perSample[sample],hapClassi,0)
                hapClassis_perSample[sample][hapClassi] += 1
        #/
        
        IDE_numAOhas.append(len(AO_has))
        # check if ovlp gfro
        ovlps = rangeOverlaps_lookup(SV['rcoords'],gene_features_alns_idx_map,1000,map_lookup_key=SV['rname'])
        ovlp_genes = []
        for oS,idx,oE in ovlps:
            geneAln = gene_features_alns[idx]
            gene_id = geneAln['qname']
        
            IDE_genes.add(gene_id)
            
            for sample in AO_has_AFpass:
                init(IDE_genes_perSample,sample,set())
                IDE_genes_perSample[sample].add(gene_id)
                
                init(IDE_gene_sampleHas,gene_id,set())
                IDE_gene_sampleHas[gene_id].add(sample)
            
            # Add all, regardless of supp
            for sample,data in SV['samples'].items():
                if min(data['AFs']) >= 0.5:
                    init(IDE_variant_sampleHas_wGeneData,SV_idx,{'samples':set(),'gene_id':gene_id})
                    IDE_variant_sampleHas_wGeneData[SV_idx]['samples'].add(sample)
            #/
            
        # "number of ovlp between samples"-analysis
        samples_hasAF = set()
        for sample,data in SV['samples'].items():
            if min(data['AFs']) >= 0.5:        
                samples_hasAF.add(sample)
        if samples_hasAF:
            ID = '||'.join(sorted(samples_hasAF))
            init(IDE_sample_variants,ID,[])
            IDE_sample_variants[ID].append(SV_idx)
        #/
            
        if ovlps and len(AO_has_AFpass) >= 1:
            ID = '||'.join(sorted(AO_has))
            init(IDE_multiAOhas_genes,ID,set())
            IDE_multiAOhas_genes[ID].add(gene_id)
            IDE_SVs.append(SV)

if 0:
    for sample,genes in IDE_genes_perSample.items():
        with open('sample.'+sample+'.list','w') as nf:
            nf.write('\n'.join(genes))
    for sampleID,genes in IDE_multiAOhas_genes.items():
        with open('sampleID.'+sampleID+'.list','w') as nf:
            nf.write('\n'.join(genes))
            
    with open('abus.panSNP.S1S2S3.list','w') as nf:
        for sampleID,SV_idxs in IDE_sample_variants.items():
            nf.write(sampleID+'\t'+str(len(SV_idxs))+'\n')

    """
    with open('geneSampleHas.list','w') as nf:
        nf.write('gene\t'+'\t'.join(sLines_IDs)+'\n')
        for geneID,sampleHas in sorted(IDE_gene_sampleHas.items(),key=lambda x: len(x[1]), reverse=True):
            writeArr = []
            for sample in sLines_IDs:
                if sample in sampleHas:
                    writeArr.append('1')
                else:
                    writeArr.append('0')
            nf.write(geneID+'\t'+'\t'.join(writeArr)+'\n')
    """
    
    with open('variant_wGeneID_sampleHas.len1.numRO2.AO2.list','w') as nf:
        nf.write('SV_idx\trname\trstart[TEFLONref]\trend[TEFLONref]\tgene\t'+'\t'.join(sLines_IDs)+'\n')
        for SV_idx,data in sorted(IDE_variant_sampleHas_wGeneData.items(),key=lambda x: len(x[1]), reverse=True):
            SV = SV_idx_map[SV_idx]
            writeArr = [SV_idx,SV['rname'],SV['rcoords'][0],SV['rcoords'][-1],data['gene_id']]
            for sample in sLines_IDs:
                if sample in data['samples']:
                    writeArr.append('1')
                else:
                    writeArr.append('0')
            nf.write('\t'.join(map(str,writeArr))+'\n')
        
#/

if 0 and 'plot?':
    plotDict = AFs_perSample
    for sample,vals in plotDict.items():
        filt0 = vals
        filt1 = filterArr( filt0, 0.01, 'threshold' )
        filt2 = filterArr( filt1, 20, 'cutoff' )
        fig = plt.figure(figsize=(11,11))
        plt.hist( filt2 ,bins=100)
        plt.title(sample+' || INFO numVals='+str(len(filt2)))
###/

# Find SNPs which exist in multiple S1+S2+S3 (evidence against cycling)
AFs_perSample = {}
for SV_idx,SV in sLines_SV_idx_map.items():
    RO_only = []
    AO_has = []
    AO_has_AFpass = []
    AF_sums = []
    DPs = []
    for sample,data in SV['samples'].items():
        #if sample == 'SRR497713' and getRangeOvlp([data['DP']]*2,[8,12]) < 0: continue
        #if sample == 'S2' and getRangeOvlp([data['DP']]*2,[8,12]) < 0: continue
        #if sample == 'SRR497721' and getRangeOvlp([data['DP']]*2,[7,11]) < 0: continue
    
        AF_sums.append(sum(data['AFs']))
    
    if len(AF_sums) != 3: continue
        
    if min(AF_sums) >= 0.2:
        for sample,data in SV['samples'].items():
            #if data['DP'] < 7 or data['DP'] > 15: continue
            init(AFs_perSample,sample,[])
            #AFs_perSample[sample].append(sum(data['AFs']))
            AFs_perSample[sample].append(sum(data['AFs']))

if 0 and 'plot?':
    for sample,AFs in AFs_perSample.items():
        plotHist(AFs,title=sample)
#/
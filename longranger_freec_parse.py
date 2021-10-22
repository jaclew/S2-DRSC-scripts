
import matplotlib.pyplot as plt
from statistics import mean,median,stdev
import pickle
import vcf as pyvcf
import sys

from functions import *

### INPUTS
chroms = ['2L','2R','3L','3R','4','X']
autosomes = ['2L','2R','3L','3R']
X = ['X']
chroms_to_plot = ('X','2L','2R','3L','3R','4','mitochondrion_genome')

vcf_file = 'phased_variants.vcf.gz' # Longranger output variant calls
freec_file = 'phased_possorted_bam.bam_ratio.txt' #CONTROL-FREEC output file
bedgraph_file = '10X_LR.bedgraph' # bedgraph coverage from short-read alignments

GTF_file = 'dmel-all-chromosome-r6.12.gtf' #gene annotations
masking_file = 'RepeatMasked/dmel-all-chromosome-r6.12.fasta.out' #repeatmasker output file

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
    bin_step = 50
    rname_BGcovs = {}
    for rname in chroms_to_plot:
        if rname in ('mitochondrion_genome',): continue
        for binnstart in range(0,refLens[rname]+bin_step,bin_step):
            rcoords = [binnstart,binnstart+bin_step]
            rstart,rend = rcoords
            
            maskCovFrac = parse_masked_covFrac(rname,rstart,rend)           
            
            BG_ovlps = rangeOverlaps_lookup(rcoords,BG_idx_map_ROMD,1000,map_lookup_key=rname)

            BG_ovlps_covs = []
            for oS,idx,oE in BG_ovlps:
                if oS == oE: continue
                masked_covFrac = parse_masked_covFrac(rname,oS,oE)   
                #if masked_covFrac >= 0.3: continue #skip repeatmasked
                BG = BG_idx_map[idx]
                BG_ovlps_covs += [BG['cov']]*(oE-oS)
            
            if not len(BG_ovlps_covs) >= bin_step*0.8: continue
        
            binn_cov = median(BG_ovlps_covs)
            
            if maskCovFrac >= 0.3 and binn_cov > 1000: continue #skip if masked AND high cov
            
            init(rname_BGcovs,rname,[])
            rname_BGcovs[rname].append([rcoords[0],binn_cov,rcoords[-1]])
    
    plt.rc('font',size=12)
    fig = plt.figure(figsize=(25*0.7,20*0.7))
    colmap = plt.cm.gist_rainbow(np.linspace(0,1,5)) #start,step,num_to_generate
    for enum,(rname,BGbins) in enumerate(sorted(rname_BGcovs.items(),key=lambda x: len(x[1]),reverse=True)):
        if rname == '4': continue
        filt0 = [ BGbin[1] for BGbin in BGbins ]
        filt1 = filterArr( filt0, 1, 'threshold' )
        filt2 = filterArr( filt1, 400, 'cutoff' )
        if not filt2: continue
        numBins = 200 #int(max(filt2)-min(filt2))
        plt.hist( filt2 ,bins=numBins,label=rname,histtype='step',stacked=True,fill=False,linewidth=1.5,alpha=0.5)
    plt.locator_params(nbins=50)
    plt.xlim([0,400])
    plt.legend()
    plt.title(str(rname)+' numVals='+str(len(filt2)))
    
    
    if 0 and 'dump to file?':
        with open(bedgraph_file.replace('.bedgraph','')+'.'+str(bin_step)+'BP_bins.bedgraph','w') as nf:
            for rname in sorted(rname_BGcovs):
                for binn in sorted(rname_BGcovs[rname], key=lambda x: x[0]):
                    writeArr = [rname,binn[0],binn[-1],binn[1]]
                    nf.write('\t'.join(map(str,writeArr))+'\n')
#/
###/
"""
### Import FREEC
freec_idx_map = {}
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
        
        freec_idx_map[tmp_freecE['idx']] = tmp_freecE

freec_idx_map_ROMD = rangeOverlaps_makeMappingDict(freec_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
###/

### Import VCF
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
    pos1 = entry.POS
    tmp_pos_partners = set()
    
    bxRO = len(entry.samples[0]['BX'][0].split(';'))
    bxROs = HOLDER(set(entry.samples[0]['BX'][0].split(';')))
    
    for idx,ALT in enumerate(entry.ALT):
        try:
            typee = entry.INFO['TYPE'][idx]
        except:
            typee = 'NA'
        try:
            SV_len = entry.INFO['LEN'][idx]
        except:
            SV_len = -1
        try:
            DP = entry.INFO['DP']
        except:
            DP = 0
        
        AO = entry.INFO['AO'][idx]
        AF = entry.INFO['AF'][idx]
        RO = entry.INFO['RO']
        
        bxAO = len(entry.samples[0]['BX'][idx+1].split(';')) #idx=0 is ref at this field
        bxAOs = HOLDER(set(entry.samples[0]['BX'][idx+1].split(';')))
        
        # Check if current entry + alt is a snp. then set length to 0.
        if typee == 'snp': SV_len = 0 
        
        if 1 and 'IDE':
            INFO_AFs.append(entry.INFO['AF'][idx])
            
            init(INFO_types,typee,0)
            INFO_types[typee] += 1
            
            init(INFO_lens,SV_len,0)
            INFO_lens[SV_len] += 1

        pos2 = pos1 + SV_len
        if pos2 <= pos1:         pos2 = pos1+1
        
        if 'IDE' and typee == 'complex':
            nucl_count = 0
            for nucl in ('A','T','G','C'):
                nucl_count += ALT.sequence.count(nucl)
            if not nucl_count == len(ALT.sequence):
                print(ALT)
        
        # Check numBP ovlp to maskings
        ovlps = rangeOverlaps_lookup([pos1,pos2],masking_idx_map_ROMD,100,map_lookup_key=rname)
        maskRanges = []
        for oS,maskIdx,oE in ovlps:
            if oE==oS: oE += 1 #handle ovlpranges fuckup
            maskRanges.append([oS,maskIdx,oE])
        ovlps2 = computeOvlpRanges_wIdxs(maskRanges)
        maskClassis = {}
        maskClassis_detailed = {}
        masked_numBP = 0
        for rangee in ovlps2:
            masked_numBP += rangee[-1]-rangee[0]
            for maskIdx in rangee[1]:
                maskE = masking_idx_map[maskIdx]
                classi = maskE['class']
                subclassi = maskE['repeat_id']
                init(maskClassis,classi,0)
                maskClassis[classi] += rangee[-1]-rangee[0]
                init(maskClassis_detailed,subclassi,0)
                maskClassis_detailed[subclassi] += rangee[-1]-rangee[0]
        masked_covFrac = masked_numBP / (pos2-pos1)
        #/
        
        # Save
        #SV_entry = {'vcf_obj':entry,'rname':chrom1,'rcoords':[pos1,pos2]}
        SV_entry = {'type':typee,'rname':chrom1,'rcoords':[pos1,pos2],'len':SV_len,'ALT':ALT.sequence,'REF':entry.REF,
                    'AO':AO,'DP':DP,'AF':AF,'RO':RO,'masked_covFrac':masked_covFrac,'masked_numBP':masked_numBP,
                    'maskings':maskClassis,'maskings_detailed':maskClassis_detailed,
                    'bxRO':bxRO,'bxAO':bxAO,'bxROs':bxROs,'bxAOs':bxAOs}

        assignIdx(SV_entry)
        
        SV_idx_map[SV_entry['idx']] = SV_entry
        tmp_pos_partners.add(SV_entry['idx'])
        
    if len(tmp_pos_partners) >= 2:
        for tmp_SV_idx in tmp_pos_partners:
            SV_idx_map[tmp_SV_idx]['pos_partners'] = tmp_pos_partners.difference(set([tmp_SV_idx]))

SV_idx_map_ROMD = rangeOverlaps_makeMappingDict(SV_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
###/IMPORT VCF

### Assign SV ovlps to masking, gtf, freecs
for SV_idx,SV in SV_idx_map.items():
    rname = SV['rname']
    rcoords = SV['rcoords']
    
    # Masking
    ovlps_masking = rangeOverlaps_lookup(rcoords,masking_idx_map_ROMD,100,map_lookup_key=rname)
    if ovlps_masking:
        SV['isMasked'] = True
    else:
        SV['isMasked'] = False
    #/
    # FREEC
    ovlps_freec = rangeOverlaps_lookup(rcoords,freec_idx_map_ROMD,1000,map_lookup_key=rname)
    if ovlps_freec:
        top_freec = sorted(ovlps_freec,key=lambda x: (x[-1]-x[0]),reverse=True)[0]
        freecE = freec_idx_map[top_freec[1]]
        SV['freec'] = freecE['CN']
    else:
        SV['freec'] = -1
    #/
        
###/

### INFO
if 0:
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
    
    INFO_SV_types = {}
    INFO_INS_SEQS_byLen = {}
    INFO_ROs = []
    INFO_AFs = []
    INFO_AOs = []
    INFO_DPs = []
    INFO_AFros = []
    INFO_bxAOs = []
    INFO_bxROs = []
    INFO_bxAFs = []
    INFO_bxDPs = []
    INFO_freecs = []
    INFO_SVlens = []
    INFO_seq_counts = {}
    AFs_hapClassi = {}
    num_RO_0 = 0
    SNP_rnameLocs_by_AO = {}
    
    read_cutoff = 10
    for SV_idx,SV in SV_idx_map.items():
        if not SV['rname'] in ('2R','2L','3R','3L'): continue
        if SV['rname'] == 'X': continue
        
        """
        # Ide: check if ovlp with x2-hap SVs called from lrsvhq_chromothripsis script
        if not SV['rname'] in x2hap_rRanges: continue
        skip_SV = True
        hit_entry = None
        for rRange in x2hap_rRanges[SV['rname']]:
            if getRangeOvlp(rRange,SV['rcoords']) >= 0:
                skip_SV = False
                hit_entry = rRange
                break
        if skip_SV: continue
        #/
        """
    
        #if 'pos_partners' in SV: continue
        #if not SV['rname'] in chroms_to_plot: continue
        SV_len = SV['len']
        SV_type = SV['type']
        #if not SV_type == 'snp': continue
        freec = SV['freec']
        if not freec == 4: continue
        INFO_freecs.append(freec)
    
        AF = SV['AF']
        DP = SV['DP']
        AO = SV['AO']
        RO = SV['RO']
        bxAO = SV['bxAO']
        bxRO = SV['bxRO']
        bxDP = bxAO+bxRO
        
        INFO_bxDPs.append(bxDP)
        
        #if DP < 102-10 or DP > 102+10: continue
        #if bxDP < 105-int(105/8) or bxDP > 105+int(105/8): continue
        cov_avg = 101
        if bxDP < cov_avg-int(cov_avg/8) or bxDP > cov_avg+int(cov_avg/8): continue
        
        if SV['DP']:
            AFro = SV['RO'] / SV['DP']
        else:
            AFro = -1
        if SV['RO'] or SV['AO']:
            AFmy = SV['AO'] / (SV['AO']+SV['RO'])
        else:
            AFmy = -1
    
        #if AF == 1:
        if RO == 0:
            bxRO = 0 #correct bullshit when bxRO has a count but AF is 1. We get 0 bxAF = 0.
            num_RO_0 += 1
            continue
        
        if (bxRO < read_cutoff or bxAO < read_cutoff): continue
        
        bxAF = bxAO / (bxRO+bxAO)
        
        INFO_bxAOs.append(bxAO)
        INFO_bxROs.append(bxRO)
        INFO_bxAFs.append(bxAF)
    
        #if AO < 30: continue
        INFO_DPs.append(DP)
        #INFO_AOs.append(AO)
        INFO_ROs.append(RO)
        INFO_AFs.append(AFmy)
        INFO_AFros.append(AFro)
        
        init(INFO_SV_types,SV_type,0)
        INFO_SV_types[SV_type] += 1
        
        INFO_SVlens.append(SV_len)
        
        hapClassi = assignHaplotype(bxAF,0.25,0.125)
        init(AFs_hapClassi,hapClassi,0)
        AFs_hapClassi[hapClassi] += 1

    filt0 = INFO_bxAFs
    filt1 = filterArr( filt0, 0, 'threshold' )
    filt2 = filterArr( filt1, 1, 'cutoff' )
    fig = plt.figure(figsize=(11,11))
    plt.hist( filt2 ,bins=100)
    plt.title('INFO numVals='+str(len(filt2)))
    plt.ylim([0,30000])
    if 0 and 'save plot?':  plt.savefig(output_dir+'bxAFs_'+str(read_cutoff)+'.svg')
    
    filt0 = INFO_bxAOs
    filt1 = filterArr( filt0, 0, 'threshold' )
    filt2 = filterArr( filt1, 3, 'cap' )
    fig = plt.figure(figsize=(11,11))
    plt.hist( filt2 ,bins=100)
    plt.title('INFO numVals='+str(len(filt2)))
    ###/INFO
    
    if 0 and 'dump for SNP + INS AO overview?':
        SNP_rnameLocs_by_AO = {}
        cov_avg = 101
        for SV_idx,SV in SV_idx_map.items():
            if SV['RO'] == 0:
                continue
            
            freec = SV['freec']
            bxAO = SV['bxAO']
            bxRO = SV['bxRO']
            bxDP = bxAO+bxRO
            bxAF = bxAO / (bxRO+bxAO)
            
            hapClassi = None
            if SV['rname'] in ('2R','2L','3R','3L'):
                if not freec == 4: continue
                if bxDP < cov_avg-int(cov_avg/8) or bxDP > cov_avg+int(cov_avg/8): continue
                hapClassi = assignHaplotype(bxAF,0.25,0.125)
                
            if SV['rname'] == 'X':
                # Mult DP by 2, since X chromosome
                if bxDP*2 < cov_avg-int(cov_avg/8) or bxDP*2 > cov_avg+int(cov_avg/8): continue
                hapClassi = assignHaplotype(bxAF,0.25,0.125)
            
            
            if hapClassi != None:
                init(SNP_rnameLocs_by_AO,SV['rname'],[])
                SNP_rnameLocs_by_AO[SV['rname']].append([SV['rname'],SV['rcoords'],hapClassi])
        
        with open('SNPs.pickle','wb') as nf:
            pickle.dump(SNP_rnameLocs_by_AO,nf)
    ### ?
    IDE_traversed = 0
    IDE_repeatSeqs = {}
    IDE_repeatSeqs['gypsy'] = {'sigs':['TACA','TGTA'],'SV_idxs':set()}
    IDE_repeatSeqs['copia'] = {'sigs':['TATA','TATA'],'SV_idxs':set()}
    for SV_idx,SV in SV_idx_map.items():
        if not SV['rname'] in chroms_to_plot: continue
        SV_len = SV['len']
        SV_type = SV['type']
        
        init(INFO_SV_types,SV_type,0)
        INFO_SV_types[SV_type] += 1
        
        INFO_SVlens.append(SV_len)
        
        if SV_len >= 4:
            IDE_traversed += 1
            seq = SV['seq']
            for repeat,data in IDE_repeatSeqs.items():
                for sig in data['sigs']:
                    if seq.find(sig) != -1 and SV['REF'].find(sig) == -1:
                        data['SV_idxs'].add(SV_idx)
###/

### Link SNPs via BCIDS
if 0:
    # 1. Create BCID -> SV_idx map
    BCID_SVs = {}
    SVs_BCIDs = {}
    for SV_idx,SV in SV_idx_map.items():
        # Filter
        if SV['RO'] == 0:     continue
        AF = SV['AO'] / (SV['AO']+SV['RO'])
        if AF >= 0.875: continue # skip if in "all"
        if AF <= 0.125: continue # skip if in "none"
        #/
        
        init(SVs_BCIDs,SV_idx,set())
        for BCID in SV['bxROs']():
            init(BCID_SVs,BCID,set())
            BCID_SVs[BCID].add(SV_idx)
            SVs_BCIDs[SV_idx].add(BCID)
        for BCID in SV['bxAOs']():
            init(BCID_SVs,BCID,set())
            BCID_SVs[BCID].add(SV_idx)
            SVs_BCIDs[SV_idx].add(BCID)
    #/1
    # 2. Cluster SNPs by BCID
    SVs_comp = {}
    for BCID1,SV_idxs in BCID_SVs.items():
        for SV_idx1 in SV_idxs:
            for SV_idx2 in SV_idxs:
                if SV_idx1 == SV_idx2: continue #skip self
                
                # Check distance between SVs
                SV1 = SV_idx_map[SV_idx1]
                SV2 = SV_idx_map[SV_idx2]
                if SV1['rname'] == SV2['rname']:
                    SVs_dist = getRangeOvlp(SV1['rcoords'],SV2['rcoords'])     
                    if SVs_dist <= 100: continue
                    if SVs_dist >= 75000: continue
                else:
                    SVs_dist = -1
                #/
                
                comp_ID = '||'.join(map(str,sorted([SV_idx1,SV_idx2])))
                if comp_ID in SVs_comp: continue #skip if already made comp
                
                shared_BCIDs = SVs_BCIDs[SV_idx1].intersection(SVs_BCIDs[SV_idx2])
                if len(shared_BCIDs) <= 3: continue #skip if no shared SV_idxs
                
                SVs_comp[comp_ID] = {'BCIDs':shared_BCIDs,'dist':SVs_dist}
    #/
    if 0 and 'plot':
        numSharedBCIDs = []
        rname_AFs = {}
        for comp_ID,comp_data in SVs_comp.items():
            numSharedBCIDs.append(len(comp_data['BCIDs']))
            SV_idx1,SV_idx2 = list(map(int,comp_ID.split('||')))
            for SV_idx in (SV_idx1,SV_idx2,):
                SV = SV_idx_map[SV_idx]
                init(rname_AFs,SV['rname'],[])
                rname_AFs[SV['rname']].append(SV['AO'])
        plotHist(numSharedBCIDs,threshold=10,cutoff=200)
        plotHist(rname_AFs['2L'])
    # 3. Run homologous recombination check
    num_hit = 0
    for comp_ID,comp_data in SVs_comp.items():
        SV_idx1,SV_idx2 = list(map(int,comp_ID.split('||')))
        SV1 = SV_idx_map[SV_idx1]
        SV2 = SV_idx_map[SV_idx2]
        
        if not (SV1['rname'] in autosomes and SV2['rname'] in autosomes): continue
    
        # Select by depth
        SVs_passed = []
        for SV in (SV1,SV2,):
            bxAO = SV['bxAO']
            bxRO = SV['bxRO']
            bxDP = bxAO+bxRO
            
            #if bxDP < 105-int(105/8) or bxDP > 105+int(105/8): continue
            if not SV['freec'] == 4: continue
            SVs_passed.append(SV['idx'])
        if not len(SVs_passed) == 2: continue
        #/
        
        # Sort reads as span/supp at SVs
        bcid_status = {}
        for bcid in SV1['bxROs']():
            init(bcid_status,bcid,set())
            bcid_status[bcid].add('span1')
        for bcid in SV1['bxAOs']():
            init(bcid_status,bcid,set())
            bcid_status[bcid].add('supp1')
        for bcid in SV2['bxROs']():
            init(bcid_status,bcid,set())
            bcid_status[bcid].add('span2')
        for bcid in SV2['bxAOs']():
            init(bcid_status,bcid,set())
            bcid_status[bcid].add('supp2')
        #/
        # Purge bcids which both span+supp at same SV
        rm_bcids = set()
        for bcid,statuses in bcid_status.items():
            for intstr in ('1','2',):
                if 'span'+intstr in statuses and 'supp'+intstr in statuses:
                    rm_bcids.add(bcid)
        for rm_bcid in rm_bcids:
            del bcid_status[rm_bcid]
        #/
        # Compile SVs span/supps bcids
        status_bcids = {}
        for bcid,statuses in bcid_status.items():
            for status in statuses:
                init(status_bcids,status,set())
                status_bcids[status].add(bcid)
        #/
        # Run SV span/supps check
        if not len(status_bcids) == 4: continue
    
        su1su2 = status_bcids['supp1'].intersection(status_bcids['supp2'])
        su1sp2 = status_bcids['supp1'].intersection(status_bcids['span2'])
        su2sp1 = status_bcids['supp2'].intersection(status_bcids['span1'])
        sp1sp2 = status_bcids['span1'].intersection(status_bcids['span2'])
        #/
        # Test?
        if su1su2 and su1sp2 and su2sp1:
            num_hit += 1
            
            print(SV1['rname'],SV2['rname'],SV1['rcoords'],SV2['rcoords'])
            #sys.exit()
        #/
    #/
###/
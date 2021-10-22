from matplotlib import pyplot as plt
import math
import numpy as np
import scipy.stats
import sys
import pysam
from statistics import mean,median,stdev


# Module-globals
class dummyClass():
    """
    Stackoverflow: Global variables in python are intra-module only. So when we import functions here,
                    but assign the global variable in main script, it will not be used in imported functions.
                    Solution will be to have a class in which I store variables I wish to be used across modules.
                    
    In scripts, example:
        global_vars_modules.global_idx_map = 0
    """
    def __repr__(self):
        return str(self.__dict__)
    def __repr__str(self):
        return str(self.__dict__)
global_vars_modules = dummyClass()
#/

####### FUNCTIONS

class fitData_tTest():
    """
    Inputs: arrays of X,Y values.
    Currently:
        Throw dataset at it, make regression.
        -Can then test datapoints from dataset to see if they fit in CI / PI, etc.
    To improve:
        Make predictions available for points outside of dataset scope (e.g. outside dataset boundaries, or within boundaries but not in model generation)
    """
    def __init__(self,x_raw,y_raw,polynomial_degree,t_test_probability):
        self.x_raw = x_raw
        self.y_raw = y_raw
        self.polynomial_degree = polynomial_degree
        self.t_test_probability = t_test_probability
        
        # Make regression on init
        self.sortInput()
        self.makeModel()

    def sortInput(self):
        # sort input, based on X, needed for the regression
        x,y = [],[]
        for X,Y in sorted(zip(self.x_raw,self.y_raw),key=lambda x: x[0]):
            x.append(X)
            y.append(Y)
        self.x = x
        self.y = y
            
    def makeModel(self):
        x = self.x
        y = self.y
        eq_poly = np.polyfit(x,y,self.polynomial_degree)
        yfit = np.polyval(eq_poly,x)
        
        eq_fit_size = eq_poly.size
        num_obs = len(yfit)
        deg_of_freedom = num_obs - eq_fit_size
        # make t-test for CI/PI
        t = scipy.stats.t.ppf(self.t_test_probability, deg_of_freedom)
        
        # calc residuals from input x,y to their polynomial of fit
        residuals = y - yfit
    
        stdev_of_err = np.sqrt(np.sum(residuals**2) / deg_of_freedom)
    
        # Calc CI,PI
        CI = t * stdev_of_err * np.sqrt( (1/num_obs) + ( (x-np.mean(x))**2 / np.sum((x-np.mean(x))**2) ) )
        PI = t * stdev_of_err * np.sqrt(1 + (1/num_obs) + ( (x - np.mean(x))**2 / np.sum((x - np.mean(x))**2) ) )
        
        # save
        self.yfit = yfit
        self.CI_lower = yfit-CI
        self.CI_upper = yfit+CI
        self.PI_lower = yfit-PI
        self.PI_upper = yfit+PI
        
        self.plotPoints = [] #will be set in internal functions for testing predictions <x,y,col>
    
    def plot(self,showCI=False,showPI=False):
        fig = plt.figure(figsize=(11,11))
        
        plt.scatter( self.x_raw,self.y_raw ) #x,y - plot raws. see if sorting fucks up? shouldnt :O
        plt.plot(self.x,self.yfit,'black')
        if showCI:
            plt.plot(self.x,self.CI_lower,'blue')
            plt.plot(self.x,self.CI_upper,'blue')
        if showPI:
            plt.plot(self.x,self.PI_lower,'red')
            plt.plot(self.x,self.PI_upper,'red')
            
        for x,y,col in self.plotPoints:
            plt.scatter([x],[y],color=col)
    
        plt.title('fitData_tTest')
        
    def clearMarkedPoints(self):
        self.plotPoints=[] #reset
        
    def __findDatasetPoints__(self,inp_x,inp_y):
        # Two scenarios:
        # 1. inp_x exist in model_x: then return that X and Y. (or those, if multiple)
        # 2. inp_x does not exist in model_x: then find <model_x1,inp_x,model_x2> + <model_y1,model_y2>. <--- TODO
        
        x_idxs = []
        y_idxs = []
        # Check 1/
        for enum,xval in enumerate(self.x):
            if xval == inp_x:
                x_idxs.append(enum)
                y_idxs.append(enum)
                
        return x_idxs,y_idxs
    
    def checkFit_CI(self,inp_x,inp_y,showPlot=False,loadPlot=False):
        x_idxs,y_idxs = self.__findDatasetPoints__(inp_x,inp_y)
            
        CI_upper = mean( [ self.CI_upper[x_idx] for x_idx in x_idxs ] )
        CI_lower = mean( [ self.CI_lower[x_idx] for x_idx in x_idxs ] )
        
        # check if input Y fit within CI
        yval_fit = False
        if inp_y >= CI_lower and inp_y <= CI_upper:
            yval_fit = True
            
        if showPlot or loadPlot:
            if yval_fit:        col='lime'
            if not yval_fit:    col='m'
            self.plotPoints.append([inp_x,inp_y,col])
            if showPlot:
                self.plot(showCI=True)
        return yval_fit
            
    def checkFit_PI(self,inp_x,inp_y,showPlot=False,loadPlot=False):
        x_idxs,y_idxs = self.__findDatasetPoints__(inp_x,inp_y)
            
        PI_upper = mean( [ self.PI_upper[x_idx] for x_idx in x_idxs ] )
        PI_lower = mean( [ self.PI_lower[x_idx] for x_idx in x_idxs ] )
        
        # check if input Y fit within CI
        yval_fit = False
        if inp_y >= PI_lower and inp_y <= PI_upper:
            yval_fit = True
            
        if showPlot or loadPlot:
            if yval_fit:        col='lime'
            if not yval_fit:    col='m'
            self.plotPoints.append([inp_x,inp_y,col])
            if showPlot:
                self.plot(showPI=True)
        return yval_fit
    
def alnIDhasher(pysamAlnObj,returnReadLen=False):
    """
    WAS MADE TO CONNECT SVIM (LONG READ SV CALLER) INTO BAM ALIGNMENT PARSING. Connect SV-called reads to their alignments.
    After investigations: Cannot use fixed BAM ID in a tag, since caller takes secondary alns from "SA"-tag.
    Cigar in SA doesnt include detailed info, so we cannot use it in hash creation.
    It allows for passing of: rcoords,qcoords,strand,mapq,NM-tag. Lets use all coords + Nm tag. that should make it pretty unique.
    
    Update2: But after investigation, coords dont match up. lets instead just add a messy rname,rcoords,qcoords,NM. then find ovlps
    for inferrment.
    
    Update3: Seems it fucks up due to Qcoords!!! we must parse that from cigar!!!! Simply parse number of clippings from end!    
    """
    
    cigar = pysamAlnObj.cigar
    
    ## Cig example: 100S10M1D2I500S
    # reverse:      S005I2D1M01S001
    
    ## Parse start until hard/soft clip
    ##OBS: Must find qcoords on our own!! pysam cant handle hardclipped alns.
    ## find qstart,qend from cigar
    tmp_qstart = 0
    qend = None
    # find start
    if cigar[0][0] in (4,5): #check if hard (5) or soft-clipped (4), else we start at coords 0
        tmp_qstart = cigar[0][1]
        
    # find end, by summing cigarstring - ignoring soft/hard-clips and DEL (2) made on query by reference
    cigar_sum = 0
    cigar_add = 0 #add soft/hard-clipping to get read length
    for i,j in cigar:
        if not i in (4,5,2):
            cigar_sum += j
        elif not i == 2:
            cigar_add += j

    # Compute read length from cigar. needed for fixing reverse mapped qstart
    read_len = cigar_sum + cigar_add
    
    # Check if entry is reverse and fix poses accordingly...
    strand = None
    if pysamAlnObj.is_reverse:
        strand = 'rv'
        qend = read_len - tmp_qstart
        qstart = qend - cigar_sum
    else:
        strand = 'fw'
        qstart = tmp_qstart
        qend = tmp_qstart + cigar_sum
    ##/
    
    if len(pysamAlnObj.tags) == 0:
        print(pysamAlnObj.cigar,read_len)
    
    tmp_arr = [pysamAlnObj.reference_id,pysamAlnObj.reference_start,pysamAlnObj.reference_end,qstart,qend,strand,pysamAlnObj.get_tag('NM')]
    
    aln_id = '-'.join(map(str,tmp_arr))
    
    if returnReadLen: return aln_id,read_len
    
    return aln_id

def mkdir(path):
    import os
    if not os.path.exists(path): os.makedirs(path)
    
def tic():
    global ticToc
    ticToc = datetime.now()
    
def toc():
    global ticToc
    time_elapsed = datetime.now() - ticToc
    returnTime = str(timedelta(seconds=time_elapsed.seconds))
    return returnTime
 
def init(inp_dict,key,value):
    if not key in inp_dict:     inp_dict[key] = value
    
def existIn(inp_dict,key,value):
    if key in inp_dict and (value == 'any' or value in inp_dict[key]): return True
    else: return False
    
def sumRanges(inp_ranges):
    sum_ = 0
    for rangee in inp_ranges:
        sum_ += rangee[-1] - rangee[0]
    return sum_

def strand_isfv(strand):
    if type(strand) == dict and 'strand' in strand: strand = strand['strand'] #parse out if input is a "aln-dict"
    if strand in ('fw','fv','+'): return True
    else: return False

def strand_isrv(strand):
    if type(strand) == dict and 'strand' in strand: strand = strand['strand'] #parse out if input is a "aln-dict"
    if strand in ('rv','-'): return True
    else: return False
    
def clearKey(inp_dict,key):
    for i in inp_dict:
        if key in inp_dict[i]:
            del inp_dict[i][key]

def filterArr(input_arr,input_val,treatment):
    """ Idea of function:
        1. Take in an array.
        2. Take in treatment
            cutoff = cut all values above VAL
            threshold = cut all values below VAL
            cap = cap all values above VAL to VAL
            raise = raise all values below VAL to VAL
    """
    
    tmp_arr = []
    for val in input_arr:
        if treatment == 'cutoff' and val > input_val: continue
        elif treatment == 'threshold' and val < input_val: continue
        elif treatment == 'cap':
            if val > input_val: val = input_val
        elif treatment == 'raise':
            if val < input_val: val = input_val
        
        tmp_arr.append(val)
        
    return tmp_arr

def assignValueToBin(input_val,bin_rules_dict):
    """
    Structure of bin_rules_dict:
        {'id' : [start,end],
         'id2': [start,end],
         'id3': [coord_ovlp]}
    Will check input coord and see which rules it OVERLAPS. so must have sharp borders, not open borders.
    e.g. rules with [5,15] [15,25], then if input is 15, the IDs of both those ranges will be returned
    """
    val_bins_assigned = []
    for bin_id,bin_rule in bin_rules_dict.items():
        # Check if None in LOWER boundary, that means "inf". Then check if our val is lower than UPPER boundary.
        if bin_rule[0] == None and not bin_rule[-1] == None and input_val <= bin_rule[-1]:
            val_bins_assigned.append(bin_id)
        # Check if None in UPPER boundary, that means "inf". Then check if our val is greater than LOWER boundary.
        if bin_rule[-1] == None and not bin_rule[0] == None and input_val >= bin_rule[0]:
            val_bins_assigned.append(bin_id)
            
        # Check cases where no input rule is with none
        if bin_rule[0] != None and bin_rule[-1] != None:
            if input_val >= bin_rule[0] and input_val <= bin_rule[-1]:
                val_bins_assigned.append(bin_id)
    return val_bins_assigned

def getPacbioMainRead(qname):
    if qname[0] == 'm':
        return '/'.join(qname.split('/')[:-1])
    else:
        return qname

def assignIdx(input_entry):
    if not 'idx' in input_entry:
        input_entry['idx'] = global_vars_modules.global_idx_iterator
        global_vars_modules.global_idx_iterator += 1
        
def sort_idx_map_dict_by_entry(idx_mapping_dict,key,shuffle_idxs=False):
    idx_mapping_dict_sorted = {}
    dict_idxs = list(idx_mapping_dict.keys())
    
    if shuffle_idxs: #Goal of randomization: Have fewer overlapRanges when doing processAlignments. Dont know if this really made a diff or not.
        random.shuffle(dict_idxs)
        
    for idx in dict_idxs:
        entry = idx_mapping_dict[idx]
        init(idx_mapping_dict_sorted,entry[key],[])
        idx_mapping_dict_sorted[entry[key]].append(idx)
    return idx_mapping_dict_sorted

def importRepeatMasker(path,rSorted=False):
    # check if no input, then return empty
    if not path: return {}
    print('Import RepeatMasking ...\n\t'+path)
    
    ## Terminology:         SCs will be consiedered "REFERENCE", RepMa's will be considered "QUERY"
    masking_arr = []
    with open(path,'r') as f:
        for line in f:
            line = line.strip('\n')
            line = line.split()
            
            # Identify entry line by having a number at first entry
            if not line: continue
            if not line[0].isdigit():
                continue
            
            # Strand is '+' or 'C'
            # qstart1 and qstart2 is dependent on strand
            
            # Example lines
            # OUT-file
            #   319    4.3  0.0  0.0  211000022278031                                   1      336      (685) + (TATAA)n          Simple_repeat           1    336     (0)      1  
            #   997    4.9  0.0  0.0  211000022278031                                 345      467      (554) C HETA              LINE/I-Jockey      (4357)   1724    1602      2  
            # Corresponding GFF-file
            # 211000022278031	RepeatMasker	similarity	1	336	 4.3	 +	.	Target "Motif:(TATAA)n" 1 336
            # 211000022278031	RepeatMasker	similarity	345	467	 4.9	 -	.	Target "Motif:HETA" 1602 1724
            
            # Sometimes the line contains a '*' in end.
            if len(line) == 15: line.append('.') #append a dummy
            SW_score,percDic,percDel,percIns,rname,rstart,rend,rstartLeftMost,strand,qname,repeat_class,qstart_rightMost,qend,qstart_leftMost,ID,data = line
            rstart = int(rstart)
            rend = int(rend)
            
            # assign qstart
            if strand == '+':
                qstart = int(qstart_rightMost)
            else:
                qstart = int(qstart_leftMost)
                strand = '-'
            qend = int(qend)
            
            # IDE CHECK IF STRAND NEEDS FIXING
            if rstart > rend: print('Handle coords, REFERENCE')
            if qstart > qend: print('Handle coords, QUERY')
            #/
            
            # Save repeat to SC
            tmp_data = {}
            tmp_data['qname'] = qname+':'+rname+':'+str(rstart)+':'+str(rend)
            tmp_data['repeat_id'] = qname
            tmp_data['qcoords'] = [qstart,qend]
            tmp_data['rname'] = rname
            tmp_data['rcoords'] = [rstart,rend]
            tmp_data['strand'] = strand
            tmp_data['type'] = 'repeat'
            tmp_data['len'] = rend-rstart
            tmp_data['strand'] = strand
            tmp_data['class'] = repeat_class
            
            #if not rname in masking:        masking[rname] = []
            #masking[rname].append(dict(tmp_data))
            masking_arr.append(tmp_data)

    if rSorted:
        masking = {}
        for entry in masking_arr:
            rname = entry['rname']
            if not rname in masking:        masking[rname] = []
            masking[rname].append(entry)
    else:
        masking =  masking_arr
    print('\tDone!')
    return masking

def importGTF(path,rSorted=False,redundancyReduced=False):
    # Redundancy reduced: filter overlapping GTFs of same type + coords. (many exons have stacked variations)
    if not path: return {}
    print('Import GTF features...\n\t'+path)
    GTF_data = []
    GTF_types = set()
    with open(path,'r') as f:
        for line in f:
            line = line.strip('\n')
            line = line.split('\t')
            chrom,typee, start, end, strand, data = line[0],line[2],int(line[3]),int(line[4]), line[6], line[-1]
            gene_id = data.split('gene_id "')[1].split('";')[0]
            GTF_types.add(typee)
            
            # correct coords if 0, do not allow that
            if start == end:
                end += 1
            #/
            
            tmp_data = {'rname':chrom,'rcoords':[start,end], 'len':end-start, 'strand':strand,
                        'type':typee, 'gene_id':gene_id, 'data':data,'qname':gene_id+':'+typee } 
            
            GTF_data.append(tmp_data)
            
    ##
    if redundancyReduced:
        # Calc IDs which make multiple feature types of same gene at identical poses redundant.
        GTF_type_coord_tracker = {}
        rm_idxs = []
        for order_idx,entry in enumerate(GTF_data):
            entry_id = '_'.join( [entry['gene_id'],entry['type'],entry['rname'],'..'.join(map(str,entry['rcoords']))] )
            init(GTF_type_coord_tracker,entry_id,[])
            GTF_type_coord_tracker[entry_id].append(order_idx)
            
            # Check if we have multiple entries, then add current entry to be removed
            if len(GTF_type_coord_tracker[entry_id]) >= 2:
                rm_idxs.append(order_idx)
            
        # Remove IDXs which were marked as redundant
        for rm_idx in rm_idxs[::-1]:
            del GTF_data[rm_idx]
    ##
    
    if rSorted:
        GTF_data2 = {}
        for entry in GTF_data:
            rname = entry['rname']
            if not rname in GTF_data2:      GTF_data2[rname] = []
            GTF_data2[rname].append(entry)
        GTF_data = GTF_data2
        
    print('\tDone!')
    return GTF_data

def importGFF(path,rSorted=False):
    print('I WAS MADE FOR REPEATMASKER, CORRECT ME!')
    # check if no input, then return empty
    if not path: return {}
    print('Import GFF features...\n\t'+path)
    
    ## Terminology:         SCs will be consiedered "REFERENCE", RepMa's will be considered "QUERY"
    masking_arr = []
    with open(path,'r') as f:
        for line in f:
            if line[0] == '#': continue
            line = line.strip('\n')
            line = line.split()
            
            rname,software,feature,rstart,rend,score,strand,frame = line[:8]
            rstart = int(rstart)
            rend = int(rend)
            data = line[8:]
            qname = data[1] #TERMINOLOGY:       repMa will be "qname / qcoords"
            qstart = int(data[2])
            qend = int(data[3])
            
            # IDE CHECK IF STRAND NEEDS FIXING
            if rstart > rend: print('Handle coords, REFERENCE')
            if qstart > qend: print('Handle coords, QUERY')
            #/
            
            # Save repeat to SC
            tmp_data = {}
            tmp_data['qname'] = qname+':'+rname+':'+str(rstart)+':'+str(rend)
            tmp_data['repeat_id'] = qname
            tmp_data['qcoords'] = [qstart,qend]
            tmp_data['rname'] = rname
            tmp_data['rcoords'] = [rstart,rend]
            tmp_data['strand'] = strand
            tmp_data['type'] = 'repeat'
            tmp_data['len'] = rend-rstart
            
            #if not rname in masking:        masking[rname] = []
            #masking[rname].append(dict(tmp_data))
            masking_arr.append(tmp_data)

    if rSorted:
        masking = {}
        for entry in masking_arr:
            rname = entry['rname']
            if not rname in masking:        masking[rname] = []
            masking[rname].append(entry)
    else:
        masking =  masking_arr
    print('\tDone!')
    return masking

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

def parseSVIMentry(entry):
    chrom1 = entry.CHROM
    pos1 = entry.POS
    reads = entry.INFO['READS']
    qual = entry.QUAL
    
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
        if 'withinMainAssembly' in ALT.__dict__ and not ALT.withinMainAssembly: print('IDE: OCC')
        
        # Parse type of SV
        typee = entry.INFO['SVTYPE']
        #/
        
        # Try parsing chrom & end coord of SV type
        chrom2 = None
        if 'chr' in ALT.__dict__:   chrom2 = ALT.chr
        elif not typee == 'BND':    chrom2 = chrom1
        else:                       print('IDE: handle!!')
        
        pos2 = None
        if 'END' in entry.INFO:     pos2 = entry.INFO['END']
        elif 'pos' in ALT.__dict__: pos2 = ALT.pos
        #/
        
        # Try parsing len of SV
        SV_len = None
        if 'SVLEN' in entry.INFO:       SV_len = entry.INFO['SVLEN']
        else:
            # Len is not reported for INV. calc it myself
            if typee == 'INV':
                if not chrom1 == chrom2: sys.exit('IDE: INV had diff chroms!!!')
                SV_len = pos2-pos1

        if not type(pos2) == int:
            print('IDE: Wasnt able to parse end coord?! handle')
            sys.exit()
        #/
            
        # Try parsing AAF of SV
        aaf = 0
        try:
            aaf = entry.aaf[0]
            if len(entry.aaf) >= 2: sys.exit('IDE-check, multiple AAF?!')
        except ZeroDivisionError:
            pass
        ref_aaf = 1-aaf
        #/
        
        # Check if CUTPASTE key is in INFO
        SV_is_cutPaste = False
        if 'CUTPASTE' in entry.INFO:        SV_is_cutPaste = True
        #/
            
        
        # Save
        #SV_entry = {'vcf_obj':entry,'rname':chrom1,'rcoords':[pos1,pos2]}
        # "isDefined" states whether feature had reported length and thus is "proven". We have BND's and long INVs which
        # do not have lens reported, thus we need to treat them separately and with caution.
        SV_entry = {'type':typee,'chrom1':chrom1,'pos1':pos1,'chrom2':chrom2,'pos2':pos2,'len':SV_len,
                    'isDefined':False,'isCutPaste':SV_is_cutPaste,'reads':reads,'reads_PBfix':reads_PBfix,
                    'qual':qual,'gt':aaf,'ref_gt':ref_aaf,'reads_numNP':SV_num_NP_reads}
        
        # If SV has length assigned to it, means we have variant called by SVIM which is "safe"?!
        # Then add it as entry
        if SV_len != None:
            if pos1 > pos2: print('IDE:occ')
            #if pos1 == pos2: pos2 += 1 #adjust if start is end
            poses = [pos1,pos2]
            SV_entry['rcoords'] = [min(poses),max(poses)]
            SV_entry['rname'] = chrom1
            SV_entry['isDefined'] = True #state that it is defined
        
        # Else, means its a long SV or BND. then add breakpoints as individual SVs
        else:
            for rname,pos in ( [chrom1,pos1],[chrom2,pos2] ):
                SV_entry['rname'] = rname
                SV_entry['rcoords'] = [pos,pos] #dummy range for computeOvlpRanges
                    
    return SV_entry

def importReadSeqsFasta(path,noSequence=False,selection=set()):
    read_seqs = {}
    with open(path,'r') as f:
        #add reads
        rn=None
        seq=''
        for line in f:
            #keep storing seq to current RN
            if not line[0] == '>':
                seq += line.strip('\n')
            #if we find new rn, store old and write it down along with folded seq
            else:
                if rn:
                    rn_old = rn.strip('\n').split()[0][1:]
                    if noSequence:      seq = None
                    
                    # check if we have selection for reads
                    if not selection or rn_old in selection:
                        read_seqs[rn_old] = seq

                rn = line.strip('\n')
                seq = '' #reset
                

        #/add reads
        #add last
        if seq:
            # Check if we had only seq in input
            if rn == None:
                rn = '>NoFastaHeader' #make dummy
            rn_old = rn.strip('\n').split()[0][1:]
            if noSequence:      seq = None
            
            # check if we have selection for reads
            if not selection or rn_old in selection:
                read_seqs[rn_old] = seq
                
            seq = '' #reset
        #/add last
    return read_seqs

def getRangeOvlp(range1,range2):
    num_bp_ovlp =-( ( max([range1[-1],range2[-1]]) - min(range1[0],range2[0]) ) - ( (range1[-1]-range1[0]) + (range2[-1]-range2[0]) ) )
    return num_bp_ovlp #Return: positive vals = numBPovlp. Neg vals: Distance ranges are apart

def getOvlpRange(range1,range2):
    ovlp_start = max([ range1[0],range2[0] ])
    ovlp_end = min([ range1[-1],range2[-1] ])
    return [ovlp_start,ovlp_end]

def computeOvlpRanges_wIdxs(hit_ranges,attribute=None):
    """
    inputs: array of ranges, input of form [start,identifier,end]
    outputs: array of ranges and the active hit_ranges index for that range
    UPDATE summer2019, now support with objects too
    """   
    events = {}
    for hit_range in hit_ranges:
        start = hit_range[0]
        end = hit_range[-1]
        query = hit_range[1]
        entry_id = hit_range[2]
        
        # Check if we have attribute we should parse from hit_range[1]
        if attribute != None: query = query.attribute
        
        # init pos
        if not start in events:     events[start] = {'starts':{},'ends':{}}
        if not end in events:     events[end] = {'starts':{},'ends':{}}

        # init query
        if not query in events[start]['starts']:     events[start]['starts'][query] = set()
        if not query in events[end]['ends']:         events[end]['ends'][query] = set()
        
        # add entry
        events[start]['starts'][query].add(entry_id)
        events[end]['ends'][query].add(entry_id)

    ranges = []
    active = {}
    old_active = {}
    for pos in sorted(events):
        old_active = dict(active)
        
        ## update active
        # add query+entry_id from START
        for query in events[pos]['starts']:
            # check if query it doesnt exist, then start
            if not query in active:      active[query] = set()
            # add all query entries
            active[query].update(events[pos]['starts'][query])
          
        # remove query+entry_id from END
        remove_list = set()
        for query in events[pos]['ends']:
            # remove end count
            active[query].difference_update(events[pos]['ends'][query])
            # check if empty, then add to remove list
            if not active[query]:
                remove_list.add(query)
                
        # remove empty entries
        for query in remove_list:
            del active[query]
        ##/
        
        ## handle range
        #check if start new range
        if not old_active:
            ranges.append([pos,{},None])
            for query in active:
                for entry_id in active[query]:
                    if not query in ranges[-1][1]:      ranges[-1][1][query] = set()
                    ranges[-1][1][query].add(entry_id)
            
        #check if close old range
        if (not active) and old_active:
            ranges[-1][-1] = pos

        #check if close old + start new (tts overlap change)          
        if active and old_active and set(old_active) != set(active):
            ranges[-1][-1] = pos #close old
            ranges.append([pos,{},None])
            for query in active:
                for entry_id in active[query]:
                    if not query in ranges[-1][1]:      ranges[-1][1][query] = set()
                    ranges[-1][1][query].add(entry_id)
        ##/
        
    return ranges

def parseFreebayesDebug(file_path,SV_idx_dict):
    print('[functions.py->parseFreebayesDebug] Warning, I\'m in IDE/MVP-mode!')
    # First traverse SV_dict, compile it to format chrom->pos->ALT=>SV_idx.
    # --> WE DO IT with SV_dict, so that we can easily map out "reference" reads across alts.
    SV_dict_mapped = {}
    for SV_idx,SV in SV_idx_dict.items():
        chrom = SV['rname']
        pos = SV['rcoords'][0]
        alt = SV['alt']
        
        init(SV_dict_mapped,chrom,{})
        init(SV_dict_mapped[chrom],pos,{})
        SV_dict_mapped[chrom][pos][alt] = SV_idx
    #/
    
    ### Parse debug file
    fo = open(file_path,'r')
    chrom_heartbeat = None
    #NEEDED? tmp_reference_supports = {} #store reference calls here, since VCF is... chrom->pos->ALT->
    for ln,line in enumerate(fo):
        line = line.split()
        # Keep update tracker of chromosome
        if line[0] == 'position:':
            chrom,pos = line[1].split(':')
            chrom_heartbeat = chrom
            continue
            pos = int(pos)+1 #make +1, since debug is 0based. VCF is 1based
            cov = int(line[-1].strip('\n'))
            sys.exit()
            
        if not line[0] == 'haplo_obs': continue
        
        #if not chrom_heartbeat == '2L': continue
        pos = int(line[1])+1 #make +1, since debug is 0based. VCF is 1based
        #inquality = float(line[2])
        
        # Check if pos exist in vcf calls
        if not (chrom_heartbeat in SV_dict_mapped and pos in SV_dict_mapped[chrom_heartbeat]): continue

        call = line[3]        
        info = line[-1].split(':')
        #RGID = info[0]
        #freebayes_score1 = info[-1]
        #freebayes_score2 = info[-2]
        #variant_cigar = info[-3]
        #read_rstart = info[-4]
        #read_call = info[-5] #same ass "call" in line[3]
        #read_call_strand = info[-6]
        #read_call_pos = info[-7]
        #unknown1 = info[-8]
        #unknown2 = info[-9]
        call_type = info[-10]
        
        # Check if call exist at pos, or is reference
        if not (call in SV_dict_mapped[chrom_heartbeat][pos] or call_type == 'reference'): continue

        readname = line[-1].split(':',1)[1] # remove readgroupID
        readname = readname.split(':'+call_type)[0] #split by call type, which appears after readname
        
        ## Save: Two ways, save at ALT, or save for ALTs if current call is reference. Modify SV entry in-place
        # Call is ref-case
        if call_type == 'reference':
            for alt in SV_dict_mapped[chrom_heartbeat][pos]:
                SV_idx = SV_dict_mapped[chrom_heartbeat][pos][alt]
                SV_entry = SV_idx_dict[SV_idx] #grab entry for clarity
                init(SV_entry,'reads_ref',HOLDER(set()))
                SV_entry['reads_ref']().add(readname)
        
        # Call is alt-case
        else:
            SV_idx = SV_dict_mapped[chrom_heartbeat][pos][call]
            SV_entry = SV_idx_dict[SV_idx] #grab entry for clarity
            init(SV_entry,'reads',HOLDER(set()))
            SV_entry['reads']().add(readname)
        ##/Save
    fo.close()
    
def traverse_ranges_wChaining(ovlps,chain_dist=0,requireMarker=False):
    """ Function is made to traverse ranges, such as computeOvlpRanges. RequireMarker requires stuff at rangee[1]
    Optional inputs:
        chain_dist = chain ranges appearing within distance for same entry
        requireMarker = Skip all ranges which do not have this marker present
    """
    input_ranges_chained = {}
    for rangee in ovlps:
        if requireMarker and not requireMarker in rangee[1]: continue
    
        # Check if we have [1] stuff (computeOvlpRanges), otherwise assume its a simple range and dummy convert it to ovlpRanges format
        if not len(rangee) >= 3:
            rangee = [rangee[0],{'DUMMY:simpleRange'},rangee[-1]]
        for entry in rangee[1]:
            if requireMarker and entry == requireMarker: continue

            # init entry if non started
            if not entry in input_ranges_chained:     input_ranges_chained[entry] = []
            # init chain if none started, or init new chain if current range is not adjacent to previous
            if (not input_ranges_chained[entry]) or ((rangee[0] - input_ranges_chained[entry][-1][-1]) > chain_dist):
                input_ranges_chained[entry].append([rangee[0],rangee[-1]])
            # check if extend previous chain, then update end coord
            if ((rangee[0] - input_ranges_chained[entry][-1][-1]) <= chain_dist):
                input_ranges_chained[entry][-1][-1] = rangee[-1]
                
    # Check if we had simpleRanges as input, then return at same format
    if 'DUMMY:simpleRange' in input_ranges_chained:
        input_ranges_chained = input_ranges_chained['DUMMY:simpleRange']
        
    return input_ranges_chained

def rangeOvlp_reshapeCoords(ref_coords,rcoords,qcoords,relative_strand):
    """
    Function will compute overlap of something at reference + query, then reshape query qcoords accordingly.
    
    ref_coords = "Coordinates where something exists"
    rcoords = "Coordinates of the query, on reference"
    qcoords = "Coordinates of the query, on self"
    relative_strand = "Determines in which end on q we chop coordinates"
    """
    
    refStart,refEnd = ref_coords
    rStart,rEnd = rcoords
    qStart,qEnd = qcoords
    
    # Find where on ref they ovlp
    ref_ranges = []
    ref_ranges.append([refStart,'r',refEnd])
    ref_ranges.append([rStart,'q',rEnd])
    ref_ovlps = computeOvlpRanges_wIdxs(ref_ranges)
    ref_coords = []
    for rangee in ref_ovlps:
        if 'r' in rangee[1] and 'q' in rangee[1]:
            ref_coords.append(rangee[0])
            ref_coords.append(rangee[-1])
    
    # Check if no overlap at all, then return empty
    if not ref_coords: return []
    
    refOvlp_start,refOvlp_end = min(ref_coords),max(ref_coords)
    
    ## Find offset this causes on aln, adjust qcoords accordingly
    qStart_offset = refOvlp_start-rStart # think about it this way. "shrink start"
    qEnd_offset = rEnd-refOvlp_end #shrink end               
    # Adjust offsets based on overall differences between aln_idx qcoords rcoords discrepancy (EMPRICALLY DETERMINED WE NEED TO DO THIS, logical too)
    qr_ratio = (qcoords[-1] - qcoords[0]) / float((rcoords[-1] - rcoords[0]))
    qStart_offset = int(qStart_offset*qr_ratio)
    qEnd_offset = int(qEnd_offset*qr_ratio)
    ##/

    # Reshape rcoords and qcoords according to strand of alignment
    qcoords_reshaped = None
    if relative_strand in ('fw','+','fv'):
        qcoords_reshaped = [qcoords[0]+qStart_offset , qcoords[-1]-qEnd_offset]
    elif relative_strand in ('rv','-'):
        qcoords_reshaped = [qcoords[0]+qEnd_offset , qcoords[-1]-qStart_offset]
        
    rcoords_reshaped = [refOvlp_start,refOvlp_end]
    #/
    
    return qcoords_reshaped,rcoords_reshaped

def importPAF(path,rselection=set(),qselection=set(),saveTags=False):
    alns = {}
    if rselection: print('parsing PAF: looking for '+ str(len(rselection)) + ' targets!')
    if qselection: print('parsing PAF: looking for '+ str(len(qselection)) + ' queries!')

    # Check if input is gzipped or not
    gzipped = False
    if path.endswith('.gz'):
        import gzip
        gzipped = True
    if gzipped:     fo = gzip.open(path,'rb')
    else:           fo = open(path,'r')
        
    for ln,line in enumerate(fo):
        if gzipped:     line = line.decode('utf-8')
        line=line.strip('\n')
        line = line.split()
        
        qname,qlen,qstart,qend,strand,rname,rlen,rstart,rend,num_matches,aln_len,mapq = line[:12]
        
        if saveTags:
            tags = line[12:]
        
        # check if selection and not readname in selection
        if qselection and not qname in qselection: continue
        if rselection and not rname in rselection: continue
        
        tmp_data = { 'rname':rname, 'rlen':int(rlen),'qlen':int(qlen),'qcoords':[int(qstart),int(qend)],
                              'qname':qname,'rcoords':[int(rstart),int(rend)],
                              'num_matches':int(num_matches), 'strand':strand, 'ln':ln }    
        
        if saveTags: tmp_data['tags'] = tags
    
        if not qname in alns:    alns[qname] = []
        alns[qname].append( tmp_data ) #'mapq':int(mapq),'aln_len':int(aln_len),'num_matches':int(num_matches)
    
    return alns

def importBAMlengths(path):
    bam_fo = pysam.AlignmentFile(path,'rb')
    ref_lens = {}
    for entry in bam_fo.header['SQ']:
        ref,length = entry['SN'],entry['LN']
        ref_lens[ref] = length
    bam_fo.close()
    return ref_lens

def parseReadLenFromCIGAR(cigarTuple):
    # find end, by summing cigarstring - ignoring soft/hard-clips and DEL (2) made on query by reference
    cigar_sum = 0
    cigar_add = 0 #add soft/hard-clipping to get read length
    for i,j in cigarTuple:
        if not i in (4,5,2):
            cigar_sum += j
        elif not i == 2:
            cigar_add += j

    # Compute read length from cigar. needed for fixing reverse mapped qstart
    read_len = cigar_sum + cigar_add
    return read_len

def parseBAMentry(entry, skipLenFromCigar=False):
    rname,rstart,rend = entry.reference_name,entry.reference_start,entry.reference_end   
    qname,cigar = entry.qname,entry.cigar
    
    if entry.has_tag('NM'):
        mismatches = entry.get_tag('NM')
    else:
        mismatches = -1
    
    
    ##OBS: Must find qcoords on our own!! pysam cant handle hardclipped alns.
    if not skipLenFromCigar:
        ## find qstart,qend from cigar
        tmp_qstart = 0
        qend = None
        # find start
        if cigar[0][0] in (4,5): #check if hard (5) or soft-clipped (4), else we start at coords 0
            tmp_qstart = cigar[0][1]
            
        # find end, by summing cigarstring - ignoring soft/hard-clips and DEL (2) made on query by reference
        cigar_sum = 0
        cigar_add = 0 #add soft/hard-clipping to get read length
        for i,j in cigar:
            if not i in (4,5,2):
                cigar_sum += j
            elif not i == 2:
                cigar_add += j
    
        # Compute read length from cigar. needed for fixing reverse mapped qstart
        read_len = cigar_sum + cigar_add
    else:
        tmp_qstart = -1
        qend = -1
        read_len = 0
        cigar_sum = 0
        
    
    # Check if entry is reverse and fix poses accordingly...
    strand = None
    if entry.is_reverse:
        strand = 'rv'
        qend = read_len - tmp_qstart
        qstart = qend - cigar_sum
    else:
        strand = 'fw'
        qstart = tmp_qstart
        qend = tmp_qstart + cigar_sum
    ##/
    
    data = {'qcoords':[qstart,qend],'rcoords':[rstart,rend],
            'rname':rname,'qname':qname,'read_len':read_len,
            'strand':strand, 'MM':mismatches}
        
    return data

def printVariableMemoryUsage():
    print('Note: seems to be very inaccurate? have cases where it reports total usage of <1GB but mem in top is >15%?!?!?')
    def sizeof_fmt(num, suffix='B'):
        ''' by Fred Cirera,  https://stackoverflow.com/a/1094933/1870254, modified'''
        for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
            if abs(num) < 1024.0:
                return "%3.1f %s%s" % (num, unit, suffix)
            num /= 1024.0
        return "%.1f %s%s" % (num, 'Yi', suffix)
    
    tot_size = 0
    iter_ = 0
    for name, size in sorted(((name, sys.getsizeof(value)) for name, value in globals().items()),
                             key= lambda x: -x[1]):
        tot_size += size
        iter_ += 1
        if iter_ < 20:    
            print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))
        
    print("{:>30}: {:>8}".format('Tot size of all vars:', sizeof_fmt(tot_size)))

class HOLDER(): #Dummy-class to hold stuff, not to spam it on IDE-work
    # NOTE: MUST CALL IT, to retrieve its holdings. E.g.: 
    # saveDataInHolder = HOLDER(data)
    # retrieveDataFromHolder = saveDataInHolder()
    def __init__(self,input_):
        self.input_ = input_
    def __call__(self):
        return self.input_
    
def getQueryPoses_fromTargetPoses_byCigar(target_coords,aln_entry=None,aln_rstart=None,aln_qcoords=None,
                                          aln_strand=None,aln_cigar=None,traverse_reverse=False,
                                          reshapeToQuery=False,reshapeToTarget=False):
    
    if reshapeToQuery and reshapeToTarget:
        print('IDE: cannot reshape to both query and target at the same time!!')
        sys.exit()
    if not (reshapeToQuery or reshapeToTarget):
        #print('IDE: must reshape to query or target!! function was updated on nov2019, perhaps old code must add the "reshapeToQuery=True" in call=')
        #sys.exit()
        reshapeToQuery = True #for compilance with old code... 
    
    # Reshape-reverse, will try to get targetPoses_fromQueryPoses...
    
    # Traverse-reverse: Should we flip cigar before traversal? Useful if there are INS-seqs and we want to find
    # an INDEL which have rcoords X,X+50 - if we scan for those coords in read, then we may up with stack around Y, missing Y+50.
    # ^May be hard to understand. WHOLE TRAVERSE_REVERSE IS ONLY TESTED IN FW-mode!!! was implemented later, cant see why it would
    # be wrong though...
    
    # Check if input was on "aln_entry"-format
    if type(aln_entry) == dict:
        rstart,rend = aln_entry['rcoords']
        qstart,qend = aln_entry['qcoords']
        cig = aln_entry['cigar']()
        strand = aln_entry['strand']
    # Else, parse function args
    else:
        # Add support for rcoords/rangees input
        if type(aln_rstart) != int and len(aln_rstart) >= 2:
            rstart = aln_rstart[0]
            rend = aln_rstart[-1]
        else:
            rstart = aln_rstart
        qstart,qend = aln_qcoords
        cig = aln_cigar
        strand = aln_strand
        
    # Traverse cigar
    q_iter_direction = 1
    cig_qpos = qstart
    cig_rpos = rstart
    target_qposes = [] #must correct for qstart/qend according to aln strand
    # check if reverse, then flip cigar and traverse from "end-to-start", from reference perspective
    if strand_isrv(strand):
        cig = cig[::-1]
        cig_rpos = rstart
        cig_qpos = qend
        q_iter_direction = -1
    
    # Check if flip cigar + stuff
    iter_direction = 1
    if traverse_reverse:
        cig = cig[::-1]
        cig_rpos = rend
        cig_qpos = qend
        iter_direction = -1
        
    for cig_e,num in cig:
        # Skip clippings
        if cig_e in (4,5): continue
    
        ## Check if hit targets at ends of aln
        # [reshapeToQuery] Check if hit target pos
        if reshapeToQuery and cig_rpos in target_coords:
            tmp_save = [cig_qpos,cig_rpos]
            if not tmp_save in target_qposes:       target_qposes.append(tmp_save)
            
        # [reshapeToTarget] Check if hit target pos
        if reshapeToTarget and cig_qpos in target_coords:
            tmp_save = [cig_rpos,cig_qpos]
            if not tmp_save in target_qposes:       target_qposes.append(tmp_save)
        ##/
        
        # if match or mismatch or N, append q and r
        if cig_e in (0,3,7):
            for i in range(num):
                cig_qpos += 1 *iter_direction*q_iter_direction
                cig_rpos += 1 *iter_direction
                
                # [reshapeToQuery] Check if hit target pos
                if reshapeToQuery and cig_rpos in target_coords:
                    tmp_save = [cig_qpos,cig_rpos]
                    if not tmp_save in target_qposes:       target_qposes.append(tmp_save)
                    
                # [reshapeToTarget] Check if hit target pos
                if reshapeToTarget and cig_qpos in target_coords:
                    tmp_save = [cig_rpos,cig_qpos]
                    if not tmp_save in target_qposes:       target_qposes.append(tmp_save)
    
        # if ins, append q
        if cig_e in (1,):
            for i in range(num):
                cig_qpos += 1*iter_direction*q_iter_direction
                # [reshapeToTarget] Check if hit target pos
                if reshapeToTarget and cig_qpos in target_coords:
                    tmp_save = [cig_rpos,cig_qpos]
                    if not tmp_save in target_qposes:       target_qposes.append(tmp_save)
        
        # if del, append r
        if cig_e in (2,):
            for i in range(num):
                cig_rpos += 1*iter_direction
                
                # [reshapeToQuery] Check if hit target pos
                if reshapeToQuery and cig_rpos in target_coords:
                    tmp_save = [cig_qpos,cig_rpos]
                    if not tmp_save in target_qposes:       target_qposes.append(tmp_save)
    
    return target_qposes #RETURN IS [<reshaped target>,<target>]

def getQueryPoses_fromTargetPoses_byCigar_BKBEFOREgettargetposesfromqueryposes(target_coords,aln_entry=None,aln_rstart=None,aln_qcoords=None,
                                          aln_strand=None,aln_cigar=None,traverse_reverse=False):
    
    # Reshape-reverse, will try to get targetPoses_fromQueryPoses...
    
    # Check if input was on "aln_entry"-format
    if type(aln_entry) == dict:
        rstart,rend = aln_entry['rcoords']
        qstart,qend = aln_entry['qcoords']
        cig = aln_entry['cigar']()
        strand = aln_entry['strand']
    # Else, parse function args
    else:
        # Add support for rcoords/rangees input
        if type(aln_rstart) != int and len(aln_rstart) >= 2:
            rstart = aln_rstart[0]
        else:
            rstart = aln_rstart
        qstart,qend = aln_qcoords
        cig = aln_cigar
        strand = aln_strand
        
    # Traverse cigar
    cig_qpos = 0 #will correct later for true qstart. depends on aln span (qcoords) and strand (start from qstart or qend)
    cig_rpos = rstart
    target_qposes = [] #must correct for qstart/qend according to aln strand
    
    # Check if flip cigar + stuff
    iter_direction = 1
    if traverse_reverse:
        cig = cig[::-1]
        cig_rpos = rend
        cig_qpos = qend
        iter_direction = -1
        
    for cig_e,num in cig:
        # Skip clippings
        if cig_e in (4,5): continue
    
        # if match or mismatch or N, append q and r
        if cig_e in (0,3,7):
            for i in range(num):
                cig_qpos += 1 *iter_direction
                cig_rpos += 1 *iter_direction
                
                # Check if hit target pos
                if cig_rpos in target_coords:
                    target_qposes.append([cig_qpos,cig_rpos])
    
        # if ins, append q
        if cig_e in (1,):
            cig_qpos += num *iter_direction
        
        # if del, append r
        if cig_e in (2,):
            for i in range(num):
                cig_rpos += 1*iter_direction
                
                # Check if hit target pos
                if cig_rpos in target_coords:
                    target_qposes.append([cig_qpos,cig_rpos])
            
    target_qposes_corrected = []
    for target_qpos,target_rpos in target_qposes:
        # Check if flip due to reverse strand
        if strand in ('fw','fv','+'):
            if not iter_direction == -1:
                target_qpos_corrected = qstart + target_qpos
            else:
                target_qpos_corrected = target_qpos
        if strand in ('rv','-'):
            if not iter_direction == -1:
                target_qpos_corrected = qend - target_qpos
            else:
                target_qpos_corrected = target_qpos
            
        target_qposes_corrected.append([target_qpos_corrected,target_rpos])
        
    return target_qposes_corrected # RETURN IS the corresponding qpos for input target_rpos

def plotEventDistribution(inp_data,file_dest='',xlabel='lvl2_keys',figsize=(11,11),xticksrotate_deg=0,ylim=None,selectArrElement=None,selectArrLen=False,barScale=1,tight_layout=False):
    """ How input data should be structured:
        data[chrom][Sv_type]->number, will plot the number of each Sv type on each chrom, grouped by Sv type
    """
    data = {}
    
    # find which categories we have, fill data with 0's if needed
    tmp_categs = set()
    for chrom in inp_data:
        for sv_type in inp_data[chrom]:
            tmp_categs.add(sv_type)
    categs = sorted(tmp_categs)
    
    for chrom in inp_data:
        data[chrom] = {}
        for categ in categs:
            if not categ in inp_data[chrom]:
                data[chrom][categ] = 0
            else:
                tmp_vals = inp_data[chrom][categ]
                # Check if parse input vals
                if selectArrElement != None:
                    tmp_vals = tmp_vals[selectArrElement]
                if selectArrLen != False:
                    tmp_vals = len(tmp_vals)
                #/
                data[chrom][categ] = tmp_vals
    #/
    
    num_categs = len(categs)
    categs_idx = np.arange(num_categs)
    
    bar_width = (1/(float(num_categs)*0.5))*barScale
    x_bars = []
    colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
    fig,ax = plt.subplots(figsize=figsize)
    for idx,chrom in enumerate(sorted(data)):
        chrom_data = []
        for SV_type in sorted(data[chrom]):
            chrom_data.append(data[chrom][SV_type])
            
            
        x_bars.append(ax.bar(categs_idx+(idx*bar_width),chrom_data,bar_width,color=colors[idx]))
        
      
    ax.set_xticks(categs_idx + bar_width)
    ax.set_xticklabels(categs,rotation=xticksrotate_deg)
    x_bars_legend = []
    for x_bar in x_bars: x_bars_legend.append(x_bar[0])
    ax.legend(x_bars_legend,sorted(data))
    
    ax.set_ylabel('abundance')
    ax.set_xlabel(xlabel)
    if ylim:        ax.set_ylim(ylim)
    if tight_layout:    plt.tight_layout()
    
    if file_dest:   fig.savefig(file_dest)
    
def plotDictKeyAbundance(inp_data,file_dest='',title='',xlabel='Event types',xticksrotate_deg=0,logFrequency=False,ylim=None):
    """ How input data should be structured:
        dict[key_lvl1] -> keys_lvl2 [--> any values of those keys]
        Function will summarize all lvl2 keys on X axis, plot abundance of lvl2 keys at lvl1 keys.
    """
    
    data = {}
    for lvl1 in inp_data:
        for lvl2 in inp_data[lvl1]:
            if not inp_data[lvl1][lvl2]: continue #skip if empty
            
            # init
            init(data,lvl2,0)
            data[lvl2] += 1
            
    yvals = []
    xvals = []
    for categ,abundance in sorted(data.items()):
        if logFrequency: abundance = math.log10(abundance)
        yvals.append(abundance)
        xvals.append(categ)
    
    fig,ax = plt.subplots(figsize=(15,15))
    ax.bar(xvals,yvals)
    ax.set_xticklabels(xvals,rotation=xticksrotate_deg)   
    ax.set_title(title)
    ax.set_ylabel('abundance')
    ax.set_xlabel(xlabel)
    if ylim:        ax.set_ylim(ylim)
    
    if file_dest:   fig.savefig(file_dest)
    
def plotDict_barplot(inp_data,file_dest='',title='',xlabel='lvl2_keys',figsize=(11,11),xticksrotate_deg=0,ylim=None,selectArrElement=None,selectArrLen=False,barScale=1,ytickdens=None,printValues=False):
    """ How input data should be structured:
        data[chrom][Sv_type]->number, will plot the number of each Sv type on each chrom, grouped by Sv type
    """
    
    # find which categories we have, fill data with 0's if needed
    tmp_categs = []
    for categ,val in inp_data.items():
        tmp_categs.append(categ)
    categs = sorted(tmp_categs)
    
    yarr = []
    xarr = []
    for categ in categs:
        val = inp_data[categ]
        if selectArrElement != None:    val=val[selectArrElement]
        yarr.append(val)
        xarr.append(categ)
    #/
    
    fig,ax = plt.subplots(figsize=figsize)
    
    ax.bar(range(len(xarr)),yarr)
    ax.set_xticks(range(len(xarr)))
    ax.set_xticklabels(xarr,rotation=xticksrotate_deg)
    ax.set_title(title)
    ax.set_ylabel('abundance')
    ax.set_xlabel(xlabel)
    if ylim:        ax.set_ylim(ylim)
    if ytickdens:   ax.locator_params(nbins=ytickdens,axis='y')
    if printValues:
        for patch in ax.patches:
            ax.annotate("%.3f" % patch.get_height(), [(patch.get_x()+(patch.get_width()/float(2))),patch.get_height()*0.05] ,ha='center',va='center' )
        
    if file_dest:   fig.savefig(file_dest)
    
def plotHist(inp_arr,title=None,cutoff=None,cap=None,threshold=None,raisee=None,nbins=100):
    filt = inp_arr
    if cutoff != None:          filt = filterArr( filt, cutoff, 'cutoff' )
    if cap != None:             filt = filterArr( filt, cap, 'cap' )
    if threshold!= None:        filt = filterArr( filt, threshold, 'threshold' )
    if raisee != None:          filt = filterArr( filt, raisee, 'raise' )
    
    fig = plt.figure(figsize=(11,11))
    plt.hist( filt,bins=nbins)
    
    title_str = ''
    if title != None:       title_str = title
    plt.title(title_str+' || numVals='+str(len(filt)))

def invertScfLoc(scfLoc):
    link_inv = None
    if scfLoc.find('::tail') != -1: link_inv = scfLoc.replace('::tail','::head')
    if scfLoc.find('::head') != -1: link_inv = scfLoc.replace('::head','::tail')
    return link_inv

def reverseComplement(seq):
    # First invert input seq
    seq = seq[::-1]
    
    # then reverse complement it
    comp_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    seq_comp = ''
    for letter in seq:
        seq_comp += comp_dict[letter]

    return seq_comp

### RANGEOVERLAPS, 2 functions
def rangeOverlaps_makeMappingDict(input_dict,bin_step,coordsKey=None,startKey='rstart',endKey='rend',sortByKey=None,ignoreIfHasKeys=None):
    ## Function1: make map
    """"
    # Make functions for making overlap-queries to a mapped dict, to semi-replace computeOvlpRanges
    # For example, map GTF/SV, query alignments at them. Or, map alignments, query alignments, and so on.
    # Two functions: 1/ make map, 2/ look for overlap
    # It will procude bins of poses, per rname, containing which IDXs exist in each bin.
    # Then we can dump GTF_idx_map into it, recieve back GTF_idxs_to_check -> Bin_pos -> GTF_idxs.
    # Then have another function to look in "RANGES", take start+end coord, calc start/end bin, scan all bins in-between
    # take out all all "GTF idxs" and corresponding GTF ranges, then look for overlap between RANGEE and GTF ranges
    
    ## EXAMPLE USAGE, 1/ map GTF file with bin size 1000. 2/ Query a range on X.
    GTF_idx_map_ROMD = rangeOverlaps_makeMappingDict(GTF_idx_map,1000,coordsKey='rcoords',sortByKey='rname') #ROMD = "range ovlp mapping dict"
    test_ovlp = rangeOverlaps_lookup([132000,138000] , GTF_idx_map_ROMD , 1000 , map_lookup_key='X')

    """
    coordBins = {} # "sortByKey" -> coordBins -> Idxs present
    for idx in input_dict:
        entry = input_dict[idx]
        
        # Check if we specified to filter
        skipEntry = False
        if ignoreIfHasKeys != None:
            for key in ignoreIfHasKeys:
                if key in entry and entry[key]:
                    skipEntry = True
                    break
        if skipEntry: continue
        #/
        
        # Parse start/end
        if coordsKey != None:
            start,end = entry[coordsKey]
        else:
            start,end = entry[startKey],entry[endKey]
            
        sortBy = None
        if sortByKey != None:           sortBy = entry[sortByKey] #update sortBy if we provided key
        
        # Round up/down coords
        start_roundD = int(start/float(bin_step)) *bin_step
        end_roundU = int(end/float(bin_step)) *bin_step + bin_step
            
        # Save at map, add idx to all bins in-between start/end rounded coords
        init(coordBins,sortBy,{})
        for bin_to_save in range(start_roundD,end_roundU+bin_step,bin_step): #rangee is open end
            init(coordBins[sortBy],bin_to_save,[])
            coordBins[sortBy][bin_to_save].append([start,idx,end])
        #/Save
        
    # If no sortByKey is provided, return dict of bin->idxs
    if sortByKey == None:       return coordBins[sortByKey]
    # Else, return full ("default")
    return coordBins
    ##/Function1

def rangeOverlaps_lookup(inp_coords,mapping_dict,bin_step,map_lookup_key=None,haltOnNumBPOvlp=None,skipNumBPOvlpBelow=None):
    ## Function2: overlap scan, using map (Query-function for input vs mapping dict)
    # Might need to improve this to parse mapping dict from GLOBAL variable space, in case we do multiprocessing
        
    # pre-flight, if we provided key that the mapping dict is sorted by, parse it out
    if map_lookup_key != None:
        if not map_lookup_key in mapping_dict: return []
        mapping_dict = mapping_dict[map_lookup_key]
    
    # round input coords
    start,end = inp_coords
    start_roundD = int(start/float(bin_step)) *bin_step
    end_roundU = int(end/float(bin_step)) *bin_step + bin_step
    
    ## Get ovlp ranges
    ovlps = {} #run per target idx
    haltOnNumBPOvlp_toggle = False #for haltOnNumBPOvlp input
    for bin_to_scan in range(start_roundD,end_roundU+bin_step,bin_step): #rangee is open end
        if not bin_to_scan in mapping_dict: continue #skip binn if it doesnt exist
        # Traverse bin entries
        for (target_start,target_idx,target_end) in mapping_dict[bin_to_scan]:
            # check if we did target_idx already
            if target_idx in ovlps: continue
        
            # check if ovlp at all
            if not getRangeOvlp(inp_coords,[target_start,target_end]) >= 0: continue #0 means border ovlp/coord ovlp
            
            # calc ovlp
            ovlp_start = max([start,target_start])
            ovlp_end = min([end,target_end])
            
            numBP_ovlp = (ovlp_end-ovlp_start)
            
            # Check if we have requirements for minimal numBP ovlp for save
            if skipNumBPOvlpBelow != None and numBP_ovlp <= skipNumBPOvlpBelow:
                continue
            
            # save ovlp
            ovlps[target_idx] = [ovlp_start,ovlp_end]
            
            # Check if we want to break on first found occurrance, for instance on repeatmasking.
            if haltOnNumBPOvlp != None and (numBP_ovlp >= haltOnNumBPOvlp):
                haltOnNumBPOvlp_toggle = True
                break
        
        # Check if break outer as well
        if haltOnNumBPOvlp_toggle:
            break
    ##/
        
    # Return on "RANGES" format
    ovlps_arr = []
    for target_idx,(ovlp_start,ovlp_end) in ovlps.items():
        ovlps_arr.append([ovlp_start,target_idx,ovlp_end])
    
    return ovlps_arr
    ##/Function2
###/RANGEOVERLAPS, 2 functions
    
####/FUNCTIONS

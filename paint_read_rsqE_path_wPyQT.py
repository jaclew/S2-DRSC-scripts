from PyQt5.QtGui import * #QPainter,QPainterPath,QBrush,QPen,QFont
from PyQt5.QtWidgets import * #QApplication,QMainWindow,QWidget,QLabel
from PyQt5.QtCore import * #Qt,QCoreApplication
from PyQt5.QtSvg import *

import os
import pickle
import time

import sys
from functions import *

## INPUTS
qnames_dirOfDirs = 'reads_to_paint'
qname_rsqEs_file = 'rsqEs.pickle'
painting_output_dir = 'paints' # generated with sniffles_pygenometracks.py
##/


# Import read rsqEs
if not 'qname_rsqRanges' in globals(): # for "spammable" script. Do not import this if it is already in memory.
    with open(qname_rsqEs_file,'rb') as f:
        qname_rsqRanges = pickle.load(f)
#/
#### Function
def calc_aChroms_from_rsqEs(rsqEs,lenReq=0,show_all_for_qnames={}):
    ### Define regions to parse
    # Setup rangees per rname with BASE and rsqE's
    parse_regions_fattenRange = 10000
    rname_regions_to_parse = {}
    for enum,rsqE in enumerate(rsqEs):
        if rsqE['masked_covFrac'] <= 0.7 and ((rsqE['qcoords'][-1]-rsqE['qcoords'][0] >= lenReq) or (show_all_for_qnames and rsqE['rname'] in show_all_for_qnames)):
            rname = rsqE['rname']
            rstart,rend = rsqE['rcoords']
            init(rname_regions_to_parse,rname,[])
            rstart_fattened = rstart-parse_regions_fattenRange
            if rstart_fattened < 0: rstart_fattened = 0
            rend_fattened = rend+parse_regions_fattenRange
            rname_regions_to_parse[rname].append([rstart_fattened,enum,rend_fattened])
    for rname,rangees in rname_regions_to_parse.items():
        rname_regions_to_parse[rname].append([0,'base',sorted(rangees,key=lambda x: x[-1])[-1][-1]])
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
                    if member == 'base': continue
                    rsqE = rsqEs[member]
                    rname_regions_PRE[-1]['rposes_wFatten'].append(rangee[0])
                    rname_regions_PRE[-1]['rposes_wFatten'].append(rangee[-1])
                    rname_regions_PRE[-1]['rposes'] += rsqE['rcoords']
                    rname_regions_PRE[-1]['rsqE_idxs'].add(rsqE['idx'])
    #/
    # Extract "artifical chroms"
    aChroms = {}
    rsqE_aChromEnums = {}
    for enum,aChrom_PRE in enumerate(sorted(rname_regions_PRE,key=lambda x: x['rname'])):
        rstart,rend = min(aChrom_PRE['rposes']),max(aChrom_PRE['rposes'])
        aChroms[enum] = {'rcoords':[rstart,rend],'rname':aChrom_PRE['rname'],'len':rend-rstart}
        for rsqE_idx in aChrom_PRE['rsqE_idxs']:
            rsqE_aChromEnums[rsqE_idx] = enum
    #/
    # Return
    return aChroms,rsqE_aChromEnums
    #/
####/Function

# Paint read paths
class QTwindow(QMainWindow):
    def __init__(self,qname='',saveLoc='',rsqE_len_req=0,rsqE_show_all_for_qnames={},offset_add=33,ref_flanks=0,ref_margins=10):
        super().__init__()
        self.qname = qname
        self.saveLoc = saveLoc
        self.rsqE_len_req = rsqE_len_req
        self.rsqE_show_all_for_qnames = rsqE_show_all_for_qnames
        self.offset_add = offset_add # Y-axis distance between paints
        self.ref_margins = ref_margins # Adds X-axis distance between reference regions
        self.ref_flanks = ref_flanks # Adds margin around reference region paints
        self.title = 'Hey'
        self.top=800
        self.left=800
        self.width=800
        self.height=800
        self.InitWindow()
    
    def InitWindow(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.top,self.left,self.width,self.height)
    
    def paintEvent(self,event):
        min_x,max_x = 100,700
        min_y,max_y = 100,700
        
        if not self.qname:
            qname = '6ab7f4d2-201f-4ba4-a4c7-1af7da6529d9' # Banreg

            #qname = '5acbe08c-56b9-48ab-9334-59ba285fa966' # Xreg
            #qname = '0002e070-d7d2-4d02-9e4b-f13e39050afe'
        else:
            qname = self.qname
        rsqEs = qname_rsqRanges[qname]

        
        ## Compile offset coords for aChroms
        # Comopute aChroms from rsqEs
        aChroms,rsqE_aChromEnums = calc_aChroms_from_rsqEs(rsqEs,lenReq=self.rsqE_len_req,show_all_for_qnames=self.rsqE_show_all_for_qnames)
        #/
        
        # Shrink max_x depending on (1) how many aChroms we plot and (2) if we want margin around them
        # (painter will paint outside boundaries, so we need to pre-shrink it to fit after boundary violation)
        if self.ref_flanks > 0:
            # get all bases which are supposed to fit within coordinates
            x_bases_tot = 0
            for aChrom in aChroms.values():
                x_bases_tot += aChrom['len']
            #/
            # deduct all X-margins which we will add during plot
            x_margins = self.ref_margins*len(aChroms)
            #/
            # shrink the x_max
            pre_coord_pix_ratio = ((max_x-min_x) - x_margins) / x_bases_tot
            max_x -= pre_coord_pix_ratio * (self.ref_flanks * len(aChroms) * 1)
            #/
        #/
        
        #1. get total length of all achroms to calc coord_pix_ratio
        regions_tot_lens = 0
        for enum,achrom in aChroms.items():
            regions_tot_lens += achrom['len']
            
        #2. calc coord_pix_ratio
        achroms_paint_margin = self.ref_margins #number of pixels in-between plotted regions
        try:
            coord_pix_ratio = ((max_x-min_x) - (len(aChroms)*achroms_paint_margin)) / regions_tot_lens
        except:
            print('unable to print for '+qname)
        #3. Save offset ranges for achroms
        offset_x = achroms_paint_margin
        if self.ref_flanks > 0:         offset_x += (self.ref_flanks*coord_pix_ratio) # check if we paint references with margin
        for achrom in aChroms.values():
            achrom['x_offset'] = offset_x
            offset_x += coord_pix_ratio*achrom['len'] + achroms_paint_margin
            if self.ref_flanks > 0:         offset_x += (self.ref_flanks*coord_pix_ratio)*2 # check if we paint references with margin
            
        #4. Function to obtain x-coord from input rsqE
        def map_rsqE_rcoords_to_pix(rsqE):
            achrom = aChroms[ rsqE_aChromEnums[rsqE['idx']] ]
            # To get both rstart+rend pix, the distance to achrom start is scaled by coord_pix_ratio
            rstart_pix = achrom['x_offset'] + (rsqE['rcoords'][0] - achrom['rcoords'][0])*coord_pix_ratio
            rend_pix = achrom['x_offset'] + (rsqE['rcoords'][-1] - achrom['rcoords'][0])*coord_pix_ratio
            return rstart_pix,rend_pix
        ##/
        ###/

        ## DRAW
        painter = QPainter(self)
        
        label_dist = 15
        y_offset = 25 + label_dist * len(aChroms)
        # Draw refRegion
        #painter.setPen(QPen(Qt.blue,1,Qt.SolidLine))
        #painter.drawLine(mapCoordToPix(rcoords_reg[0]),y_offset,mapCoordToPix(rcoords_reg[-1]),y_offset)
        #/
        # Draw refRegions
        painter.setPen(QPen(Qt.blue,3,Qt.SolidLine))
        for enum,achrom in aChroms.items():
            x_start,x_end = achrom['x_offset'],achrom['x_offset']+(achrom['len']*coord_pix_ratio)
            # check if we paint references with margin
            if self.ref_flanks > 0:
                x_start -= self.ref_flanks*coord_pix_ratio
                x_end += self.ref_flanks*coord_pix_ratio
            #/
            painter.drawLine(x_start,y_offset,x_end,y_offset)
            if self.ref_flanks > 0:
                text = achrom['rname']+':'+str(round((achrom['rcoords'][0]-self.ref_flanks)/1000,2))+'-'+str(round((achrom['rcoords'][-1]+self.ref_flanks)/1000,2))+' kbp'
            else:
                text = achrom['rname']+':'+str(round(achrom['rcoords'][0]/1000,2))+'-'+str(round(achrom['rcoords'][-1]/1000,2))+' kbp'
            refLabel = QLabel(self)
            refLabel.setText(text)
            refLabel.setFont(QFont('Helvetica',12))
            refLabel.adjustSize() # for font size to update label size
            refLabel.move(x_start,y_offset-25-(label_dist*enum))
            refLabel.show()
        #/
        
        ## Draw alns & junctions
        # alns
        painter.setPen(QPen(Qt.black,2,Qt.SolidLine))
        y_offset += self.offset_add
        rsqEs_drawCoords = {} #rsqE_idx -> drawData
        rsqEs_INSes = [ {'rsqE_left':None,'rsqE_right':None,'rsqEs':[]} ] #init first
        for rsqE in sorted(rsqEs,key=lambda x: x['qcoords'][0]):
            # Check if write rsqE as "INS" (its outside plotting region)
            #if getRangeOvlp(rcoords_reg,rsqE['rcoords']) < min((rsqE['rcoords'][-1]-rsqE['rcoords'][0])*0.7,1000): #allow some liberality, we can write some outside span
            if not rsqE['idx'] in rsqE_aChromEnums:
                rsqEs_INSes[-1]['rsqEs'].append(rsqE)
                continue
            #/
            # Make new INS-section each time we write an rsqE
            rsqEs_INSes[-1]['rsqE_right'] = rsqE #close old
            rsqEs_INSes.append( {'rsqE_left':rsqE,'rsqE_right':None,'rsqEs':[]} )
            #/

            #rstart,rend = rsqE['rcoords']
            #pix_start_x,pix_end_x = mapCoordToPix(rstart),mapCoordToPix(rend)
            pix_start_x,pix_end_x = map_rsqE_rcoords_to_pix(rsqE)
            pix_start_y,pix_end_y = y_offset,y_offset
            painter.drawLine(pix_start_x,pix_start_y,pix_end_x,pix_end_y)
            rsqEs_drawCoords[rsqE['idx']] = [pix_start_x,pix_start_y,pix_end_x,pix_end_y]
            
            # draw aln strand indicator
            if strand_isfv(rsqE):           painter.drawLine(pix_end_x,pix_end_y,pix_end_x-10,pix_end_y-5)
            if strand_isrv(rsqE):           painter.drawLine(pix_start_x,pix_start_y,pix_start_x+10,pix_start_y-5)
            #/
            y_offset += self.offset_add
        
        # junctions
        junc_margin = int(self.offset_add/4)
        for junction in rsqEs_INSes:
            # parse drawn rsqEs
            rsqE_left = junction['rsqE_left']
            rsqE_right = junction['rsqE_right']
            if rsqE_left == None or rsqE_right == None: continue #for now, skip HT INSes
            
            pix_start_x,pix_start_y = rsqEs_drawCoords[rsqE_left['idx']][:2]
            pix_end_x,pix_end_y = rsqEs_drawCoords[rsqE_right['idx']][2:]

            # Draw junction (select pen color depending on size)
            penColor = Qt.red #default
            
            if junction['rsqEs']:
                juncLen = 0
                for rsqE in junction['rsqEs']:      juncLen += rsqE['qcoords'][-1]-rsqE['qcoords'][0]
                
                # Parse junction content
                maskings_numBPs = {}
                for rsqE in junction['rsqEs']:
                    if 'maskClassis_numBPs' in rsqE:
                        for repeat,numBP in rsqE['maskClassis_numBPs']:
                            init(maskings_numBPs,repeat,0)
                            maskings_numBPs[repeat] += numBP
                top_masking = 'NA'
                if maskings_numBPs:
                    #top_masking = str(sorted(maskings_numBPs.items(),key=lambda x: x[1], reverse=True))
                    top_masking = sorted(maskings_numBPs.items(),key=lambda x: x[1], reverse=True)[0][0]
                #/
                
                INSlabel = QLabel(self)
                INSlabel.setText(str(round(juncLen,-2)) + ' ' + top_masking)
                INSlabel.setFont(QFont('Helvetica',10))
                INSlabel.move(((pix_start_x+pix_end_x)/2) + junc_margin,pix_start_y)
                INSlabel.resize(600,20) #dirty way of making the label wide to fit long texts
                INSlabel.show()
                
                penColor = Qt.green #update to green if LTR
            
            painter.setPen(QPen(penColor,2,Qt.DotLine))
            painter.drawLine(pix_start_x,pix_start_y+junc_margin,pix_end_x,pix_end_y-junc_margin)
            #/        
        ##/
        #/

RUNINTERACTIVE = False
dumpRegions = True # if true, will dump regions in output (for pygenometracks)

App = QCoreApplication.instance()
if App is None:
    App = QApplication(sys.argv)

mito_reads = ['f77d70ed-28c2-4d02-9bd0-15b3c614fdd8', #3L
              '9b04600c-8446-4eb5-8e17-32d1a24270ce', #3R
              '63d7681f-d9cc-4e8f-93c3-51fe0edb2fbd'] #4

SV_trouble_reads = ['1d1e5ed0-a54d-4c08-b1f0-c8b200f8fca8',
                    '2b66788b-1981-4ae7-b812-99f56987aead',
                    '2d039b29-7a8c-4e1b-b491-c7da47ae5671',
                    '2f47cb16-cf01-4c76-b44d-7db1340ce5d0',
                    '3da36f4b-6863-4eb5-99a7-6a1e1ece9dd9',
                    '4be6f1bd-8028-4eb8-862a-b7a466998f9a',
                    '9d385047-8d01-4f76-acc9-664b1d7fae70',
                    '68c90335-d432-4874-a2f2-5660f5e096a7',
                    '201f58ea-8743-45d2-ba17-c90c8edfc089',
                    '407a58a6-d2f0-41c7-9638-e84f9161fba8',
                    'a284e11f-bbd8-4297-9c06-e6787105cb66',# "cleaved-DUP"
                    'd08527a7-2b07-4280-8f6b-d1c657b97046',# "DUP+DEL" at 3R
                    'e7b979ae-8374-4f65-8fb6-08f0eba14247',
                    ]

#qnames_to_paint = ['e7deb45b-1d36-438a-b969-7b3b0461047f', # or 'ce042b0c-11d7-48bb-a95b-20a4833ca352'
#                   'd08527a7-2b07-4280-8f6b-d1c657b97046']

qnames_to_paint = SV_trouble_reads
#for file_ in os.listdir(qnames_dirOfDirs):
for file_ in qnames_to_paint:
    #loc,qname = file_.split('_')
    qname = file_
    if qname in qname_rsqRanges:
        tmp_out = painting_output_dir+'/'+file_
        mkdir(tmp_out)
        ## Paint read path
        runPainter = QTwindow()
        runPainter.qname=qname
        
        # SETTINGS
        runPainter.rsqE_len_req = 1000 # len cutoff for rsqEs to plot
        runPainter.rsqE_show_all_for_qnames = {'mitochondrion_genome'}
        runPainter.offset_add = 20
        runPainter.ref_flanks = 10000
        runPainter.ref_margins= 10
        #/
        
        runPainter.show()
        # Save image
        picture = runPainter.grab()
        #picture.save(tmp_out+'/'+'path.'+qname+'.png') # save in subfolder
        picture.save(tmp_out+'.png') # save in same folder
        #/
        
        if not RUNINTERACTIVE:
            runPainter.close()
        else:
            break
        ##/
        ## Dump coordinates to parse pygenometrack sniffles calls
        rsqEs = qname_rsqRanges[qname]
        aChroms,rsqE_aChromEnums = calc_aChroms_from_rsqEs(rsqEs)
        if dumpRegions and 'qnames_as_subfolders':
            #print(len(aChroms),loc,qname)
            with open(tmp_out+'/'+'regions.list','w') as nf:
                for enum,achrom in aChroms.items():
                    nf.write(str(enum)+'\t'+achrom['rname']+':'+str(achrom['rcoords'][0]-1000)+'-'+str(achrom['rcoords'][-1]+1000)+'\n')
            
        ##/

if RUNINTERACTIVE:
    sys.exit(App.exec())
###/
if 0 and 'qnames_as_subfolders': print('Now run CMD.generate_pyGenomeTracks')
#/
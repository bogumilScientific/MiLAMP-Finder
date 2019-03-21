# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 02:33:51 2018

@author: michael bogumil
"""

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import SeqUtils
from Bio import Restriction as rstrn
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from io import StringIO
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from time import sleep
import primer3
import time
from random import choice


class MiLAMPFinder:
    def __init__(self):
        self.numMethSitesUsed = -2
        self.resEnzWithRange = []
        self.checkSites = []
        
        self.foldername = ''
        self.foldernames = []
        self.filename = ''
        
        self.F1cBest = 0
        self.B1cBest = 0
        self.F2Best = 0
        self.B2Best = 0
        self.F3Best = 0
        self.B3Best = 0
        self.FLPBest = 0
        self.BLPBest = 0
        
        self.customLoopPrimers = False
    
    def loadSeqFile(self, filepath, fileformat = "fasta"):
        for record in SeqIO.parse(filepath, fileformat):
            self.gene = record
            
    def loadSeqText(self, seq):
        self.gene = SeqRecord(Seq(seq, IUPAC.ambiguous_dna))
        
    def getSeqText(self):
        return str(self.gene.seq)
        
    def loadSettings(self, path):
        pass
    
    def saveSettings(self, path):
        pass
    
    def setStart(self, strtPos = 1):
        self.strtPos = strtPos
    
    def setMethSite(self, methSite = []):
        self.methSite = np.array(methSite)
        
        
    def methSiteCheck(self):
        # Perform a sanity check on the methylation site, asuring that the addresses correspond to CpG islands
        # Returns False if any addresses don't correspond to a CpG island otherwise it returns True
        noCpGFlag = True
        
        for i in (np.nditer(self.methSite-self.strtPos)):
            if (self.gene.seq[i] == 'C') and (self.gene.seq[i+1] == 'G'):
                self.checkSites.append((i + self.strtPos, True))
            else:
                self.checkSites.append((i + self.strtPos, False))
                noCpGFlag = False
                
        return noCpGFlag
    
    def setResEnz(self):
        self.rb_supp = rstrn.RestrictionBatch(first=[], suppliers=['C','B','E','I','K','J','M','O','N','Q','S','R','V','Y','X'])
    
    def numOfMethSitesUsed(self, methSitesUsed = 1):
        self.numMethSitesUsed = methSitesUsed
    
    def findResEnz(self):
        
        # Finder cuts sites on loaded gene sequence
        reEnzymeSites = self.rb_supp.search(self.gene.seq)
        
        # Loop over each restriction enzyme
        for key, value in reEnzymeSites.items() :
            
            # Use only the restriction enzymes that intersect with methylation sites
            match = set(np.array([np.arange(2*len(key)-1)-(len(key) - 1) + i for i in (self.methSite - self.strtPos)]).flatten()) & set(value)
            if (len(match) > 0):
                # Return all results regardless of number of methylation sites used (to be implemeted)
                if (self.numMethSitesUsed == -1):
                    pass
                # Return results using best methylaion site regardless of number of methylation sites used 
                if (self.numMethSitesUsed == -2):
                    bestCutDistance = -1
                    # Go thru all the methylation site
                    for i in (np.array(list(match)) + self.strtPos):
                        # Calculate which methylation site has the biggest free range(no nonmethylation cut sites) around it
                        
                        # Calculate first cut site above methylation site
                        if (len(np.array(value)[np.array(value) > i]) > 0):
                            topCut = np.min(np.array(value)[np.array(value) > i])
                        else:
                            topCut = self.strtPos + len(self.gene) - 1
                        
                        # Calculate first cut site below methylation site
                        if (len(np.array(value)[np.array(value) < i]) > 0 ):
                            bottomCut = np.max(np.array(value)[np.array(value) < i])
                        else:
                            bottomCut = self.strtPos
                        
                        # Find the distance between both non-methylation cut sites
                        cutDistance = topCut-bottomCut
                        
                        # Use the methylation site that has the largest range
                        if (cutDistance > bestCutDistance):
                            bestCutDistance = cutDistance
                            bestTopCut = topCut
                            bestBottomCut = bottomCut
                            bestMethCutSite = i
                        
                # Return results of the use of one methylation sites (to be implemeted)
                if (self.numMethSitesUsed == 1):
                    pass
                # Return results of the use of a specific number of methylation sites (to be implemeted)
                elif (self.numMethSitesUsed > 0):
                    pass
                
                # Store the restriction enzymes with the best mythylation site it has to offer
                self.resEnzWithRange.append((key, bestMethCutSite,bestCutDistance,bestTopCut,bestBottomCut,bestMethCutSite))
                
    def getCurrentEnz(self):
        return self.resEnzWithRange
    
    def setEnz(self, enzUsed = 0):
        self.useEnz = enzUsed
     
        
    def setDefaultPrimerParam(self):
        # ============ Generate Standard Loop Primers ============
        self.customLoopPrimers = False
        
        # ============ Primer3 variables ============
        self.innerLengthMax = 30
        self.innerLengthMin = 18
        self.innerLengthOpt = 20
        
        self.innerTempMax = 61
        self.innerTempMin = 59
        self.innerTempOpt = 60
        
        self.innerGCMax = 80
        self.innerGCMin = 40
        
        self.complementLengthMax = 30
        self.complementLengthMin = 18
        self.complementLengthOpt = 20
        
        self.complementTempMax = 67
        self.complementTempMin = 55
        self.complementTempOpt = 60
        
        self.complementGCMax = 90
        self.complementGCMin = 20
        
        self.loopLengthMax = 30
        self.loopLengthMin = 18
        self.loopLengthOpt = 20
        
        self.loopTempMax = 65
        self.loopTempMin = 59
        self.loopTempOpt = 60
        
        self.loopGCMax = 90
        self.loopGCMin = 20
        
        self.outerLengthMax = 30
        self.outerLengthMin = 18
        self.outerLengthOpt = 20
        
        self.outerTempMax = 65
        self.outerTempMin = 57
        self.outerTempOpt = 60
        
        self.outerGCMax = 70
        self.outerGCMin = 40
        
        self.maxHairpinTm = 35
        
        self.numToReturn = 10
        
        # ============ Sequence spacing variables ============
        
        # get methylation cut site for chosen enzyme
        self.methSite = self.resEnzWithRange[self.useEnz][1] - self.strtPos
        
        self.methSiteBuf = 10
        
        self.maxFrameSize = 350
        
        if(self.customLoopPrimers == False):
            self.stdBuf = int(np.ceil((self.maxFrameSize - 2*self.methSiteBuf)/8))
        else:
            self.stdBuf = int(np.ceil((self.maxFrameSize - 2*self.methSiteBuf)/6))
        
        self.B1start = self.methSite - self.methSiteBuf - self.stdBuf
        self.B1end = self.B1start + self.stdBuf
        
        self.F1start = self.methSite + self.methSiteBuf
        self.F1end = self.F1start + self.stdBuf
        
        if(self.customLoopPrimers == False):
            self.BLPstart = self.B1start - self.stdBuf
            self.BLPend = self.B1start
            
            self.FLPstart = self.F1end
            self.FLPend = self.FLPstart + self.stdBuf
            
        else:
            # Need to check these
            self.BLPstart = self.B1start
            self.FLPstart = self.F1end
            
            self.BLPend = self.B1start
            self.FLPend = self.F2start 
    
        
        self.B2start = self.BLPstart - self.stdBuf
        self.B2end = self.BLPstart
        
        self.F2start = self.FLPend
        self.F2end = self.F2start + self.stdBuf
        
        self.B3start = self.B2start - self.stdBuf
        self.B3end = self.B2start
        
        self.F3start = self.F2end
        self.F3end = self.F3start + self.stdBuf
        
        self.polyTinterConnect = 0
        
    def setPrimerParam(self, param):
        # ============ Generate Standard Loop Primers ============
        if 'customLoopPrimers' in param.keys():
            self.customLoopPrimers = param['customLoopPrimers']
        
        # ============ Primer3 variables ============
        if 'innerLengthMax' in param.keys():
            self.innerLengthMax = param['innerLengthMax']
        if 'innerLengthMin' in param.keys():
            self.innerLengthMin = param['innerLengthMin']            
        if 'innerLengthOpt' in param.keys():
            self.innerLengthOpt = param['innerLengthOpt']
        
        if 'innerTempMax' in param.keys():
            self.innerTempMax = param['innerTempMax']
        if 'innerTempMin' in param.keys():
            self.innerTempMin = param['innerTempMin']            
        if 'innerTempOpt' in param.keys():
            self.innerTempOpt = param['innerTempOpt']
        
        if 'innerGCMax' in param.keys():
            self.innerGCMax = param['innerGCMax']            
        if 'innerGCMin' in param.keys():
            self.innerGCMin = param['innerGCMin']
        
        if 'complementLengthMax' in param.keys():
            self.complementLengthMax = param['complementLengthMax']
        if 'complementLengthMin' in param.keys():
            self.complementLengthMin = param['complementLengthMin']            
        if 'complementLengthOpt' in param.keys():
            self.complementLengthOpt = param['complementLengthOpt']

        if 'complementTempMax' in param.keys():
            self.complementTempMax = param['complementTempMax']
        if 'complementTempMin' in param.keys():
            self.complementTempMin = param['complementTempMin']            
        if 'complementTempOpt' in param.keys():
            self.complementTempOpt = param['complementTempOpt']
        
        if 'complementGCMax' in param.keys():
            self.complementGCMax = param['complementGCMax']            
        if 'complementGCMin' in param.keys():
            self.complementGCMin = param['complementGCMin']
        
        if 'loopLengthMax' in param.keys():
            self.loopLengthMax = param['loopLengthMax']
        if 'loopLengthMin' in param.keys():
            self.loopLengthMin = param['loopLengthMin']            
        if 'loopLengthOpt' in param.keys():
            self.loopLengthOpt = param['loopLengthOpt']
        
        if 'loopTempMax' in param.keys():
            self.loopTempMax = param['loopTempMax']
        if 'loopTempMin' in param.keys():
            self.loopTempMin = param['loopTempMin']            
        if 'loopTempOpt' in param.keys():
            self.loopTempOpt = param['loopTempOpt']
        
        if 'loopGCMax' in param.keys():
            self.loopGCMax = param['loopGCMax']            
        if 'loopGCMin' in param.keys():
            self.loopGCMin = param['loopGCMin']
        
        if 'outerLengthMax' in param.keys():
            self.outerLengthMax = param['outerLengthMax']
        if 'outerLengthMin' in param.keys():
            self.outerLengthMin = param['outerLengthMin']            
        if 'outerLengthOpt' in param.keys():
            self.outerLengthOpt = param['outerLengthOpt']
        
        if 'outerTempMax' in param.keys():
            self.outerTempMax = param['outerTempMax']
        if 'outerTempMin' in param.keys():
            self.outerTempMin = param['outerTempMin']            
        if 'outerTempOpt' in param.keys():
            self.outerTempOpt = param['outerTempOpt']
        
        if 'outerGCMax' in param.keys():
            self.outerGCMax = param['outerGCMax']            
        if 'outerGCMin' in param.keys():
            self.outerGCMin = param['outerGCMin']
        
        
        if 'maxHairpinTm' in param.keys():
            self.maxHairpinTm = param['maxHairpinTm']
        
        if 'numToReturn' in param.keys():
            self.numToReturn = param['numToReturn']
        
        if 'strtPos' in param.keys():
            self.strtPos = param['strtPos']
        
        if 'methSiteBuf' in param.keys():
            self.methSiteBuf = param['methSiteBuf']
        
        if 'maxFrameSize' in param.keys():
            self.maxFrameSize = param['maxFrameSize']
        
        if 'gene' in param.keys():
            self.gene = param['gene']
            
        if 'methSite' in param.keys():
            self.methSite = param['methSite']
            
        if 'useEnz' in param.keys():
            self.useEnz = param['useEnz']
         
        if 'polyTinterConnect' in param.keys():
            self.polyTinterConnect = param['polyTinterConnect']
        # ============ Sequence spacing variables ============
        
        # get methylation cut site for chosen enzyme
        self.methSite = self.resEnzWithRange[self.useEnz][1] - self.strtPos
        

        
        if(self.customLoopPrimers == False):
            self.stdBuf = int(np.ceil((self.maxFrameSize - 2*self.methSiteBuf)/8))
        else:
            self.stdBuf = int(np.ceil((self.maxFrameSize - 2*self.methSiteBuf)/6))
        
        self.B1start = self.methSite - self.methSiteBuf - self.stdBuf
        self.B1end = self.B1start + self.stdBuf
        
        self.F1start = self.methSite + self.methSiteBuf
        self.F1end = self.F1start + self.stdBuf
        
        if(self.customLoopPrimers == False):
            self.BLPstart = self.B1start - self.stdBuf
            self.BLPend = self.B1start
            
            self.FLPstart = self.F1end
            self.FLPend = self.FLPstart + self.stdBuf
            
        else:
            # Need to check these
            self.BLPstart = self.B1start
            self.FLPstart = self.F1end
            
            self.BLPend = self.B1start
            self.FLPend = self.F2start 
    
        
        self.B2start = self.BLPstart - self.stdBuf
        self.B2end = self.BLPstart
        
        self.F2start = self.FLPend
        self.F2end = self.F2start + self.stdBuf
        
        self.B3start = self.B2start - self.stdBuf
        self.B3end = self.B2start
        
        self.F3start = self.F2end
        self.F3end = self.F3start + self.stdBuf
        
    def recalculateCutSites(self):
        if(self.customLoopPrimers == False):
            self.stdBuf = int(np.ceil((self.maxFrameSize - 2*self.methSiteBuf)/8))
        else:
            self.stdBuf = int(np.ceil((self.maxFrameSize - 2*self.methSiteBuf)/6))
        
        self.B1start = self.methSite - self.methSiteBuf - self.stdBuf
        self.B1end = self.B1start + self.stdBuf
        
        self.F1start = self.methSite + self.methSiteBuf
        self.F1end = self.F1start + self.stdBuf
        
        if(self.customLoopPrimers == False):
            self.BLPstart = self.B1start - self.stdBuf
            self.BLPend = self.B1start
            
            self.FLPstart = self.F1end
            self.FLPend = self.FLPstart + self.stdBuf
            
        else:
            # Need to check these
            self.BLPstart = self.B1start
            self.FLPstart = self.F1end
            
            self.BLPend = self.B1start
            self.FLPend = self.F2start 
    
        
        self.B2start = self.BLPstart - self.stdBuf
        self.B2end = self.BLPstart
        
        self.F2start = self.FLPend
        self.F2end = self.F2start + self.stdBuf
        
        self.B3start = self.B2start - self.stdBuf
        self.B3end = self.B2start
        
        self.F3start = self.F2end
        self.F3end = self.F3start + self.stdBuf
    
    def generateRandomPrimers(self, n=1000):
        self.cutSiteList = np.zeros((16,n))
        self.primerRating = np.zeros(n)
        
        self.cutSiteShiftSTD = 10     
        
        self.cutSiteList[0,:] = np.random.normal(loc=self.B3start,scale=self.cutSiteShiftSTD,size=n).astype(int)
        self.cutSiteList[1,:] = np.random.normal(loc=self.B3end,scale=self.cutSiteShiftSTD,size=n).astype(int)
        self.cutSiteList[2,:] = np.copy(self.cutSiteList[1,:])
        self.cutSiteList[3,:] = np.random.normal(loc=self.B2end,scale=self.cutSiteShiftSTD,size=n).astype(int)
        self.cutSiteList[4,:] = np.copy(self.cutSiteList[3,:])
        self.cutSiteList[5,:] = np.random.normal(loc=self.BLPend,scale=self.cutSiteShiftSTD,size=n).astype(int)
        self.cutSiteList[6,:] = np.copy(self.cutSiteList[5,:])
        self.cutSiteList[7,:] = np.full(n, self.B1end)
        self.cutSiteList[8,:] = np.full(n, self.F1start)
        self.cutSiteList[9,:] = np.random.normal(loc=self.F1end,scale=self.cutSiteShiftSTD,size=n).astype(int)
        self.cutSiteList[10,:] = np.copy(self.cutSiteList[9,:])
        self.cutSiteList[11,:] = np.random.normal(loc=self.FLPend,scale=self.cutSiteShiftSTD,size=n).astype(int)
        self.cutSiteList[12,:] = np.copy(self.cutSiteList[11,:])
        self.cutSiteList[13,:] = np.random.normal(loc=self.F2end,scale=self.cutSiteShiftSTD,size=n).astype(int)
        self.cutSiteList[14,:] = np.copy(self.cutSiteList[13,:])
        self.cutSiteList[15,:] = np.random.normal(loc=self.F3end,scale=self.cutSiteShiftSTD,size=n).astype(int)        
        
        for i in range(n):
            print(i)         
            self.B3start = int(self.cutSiteList[0,i])
            self.B3end = int(self.cutSiteList[1,i])
            self.B2start = int(self.cutSiteList[2,i])
            self.B2end = int(self.cutSiteList[3,i])
            self.BLPstart = int(self.cutSiteList[4,i])
            self.BLPend = int(self.cutSiteList[5,i])
            self.B1start = int(self.cutSiteList[6,i])
            self.B1end = int(self.cutSiteList[7,i])
            self.F1start = int(self.cutSiteList[8,i])
            self.F1end = int(self.cutSiteList[9,i])
            self.FLPstart = int(self.cutSiteList[10,i])
            self.FLPend = int(self.cutSiteList[11,i])
            self.F2start = int(self.cutSiteList[12,i])
            self.F2end = int(self.cutSiteList[13,i])
            self.F3start = int(self.cutSiteList[14,i])
            self.F3end = int(self.cutSiteList[15,i])
            
#            self.printAddress()
            
            if (self.areStartEndLargeEnough()==False):
                self.primerRating[i] = np.nan
                continue
            try:
                myFinder.genPrimers()
                myFinder.createXIP()
                pass
            except:
                self.primerRating[i] = np.nan
                continue
#             if a particular cut site doesn't make primers set the rating to 
#             -1, otherwise calculate the primer set rating
            if(self.arePrimersEmpty()==False):
                self.primerRating[i] = self.rateCurrentPrimerSet()
            else:
                self.primerRating[i] = np.nan
#        print(self.primerRating)
        
    def areStartEndLargeEnough(self, thres=3):
        addr = np.array([self.B3start,
                        self.B3end,
                        self.B2end,
                        self.BLPend,
                        self.B1end,
                        self.F1start,
                        self.F1end,
                        self.FLPend,
                        self.F2end,
                        self.F3end])
        diff = np.ediff1d(addr)
#        print(diff)
        return all(x > thres for x in diff)
        

    def arePrimersEmpty(self):
        if(self.primerF2['PRIMER_RIGHT_NUM_RETURNED']==0):
            return True
        elif(self.primerB2['PRIMER_LEFT_NUM_RETURNED']==0):
            return True
        elif(self.primerF3['PRIMER_RIGHT_NUM_RETURNED']==0):
            return True
        elif(self.primerB3['PRIMER_LEFT_NUM_RETURNED']==0):
            return True
        elif(self.primerF1c['PRIMER_LEFT_NUM_RETURNED']==0):
            return True
        elif(self.primerB1c['PRIMER_RIGHT_NUM_RETURNED']==0):
            return True
        elif(self.primerFLP['PRIMER_RIGHT_NUM_RETURNED']==0):
            return True
        elif(self.primerBLP['PRIMER_LEFT_NUM_RETURNED']==0):
            return True
        else:
            return False

        


    def rateCurrentPrimerSet(self):
        tempWeight = 1/20000
        hairpinWeight = 1/700
        GCPWeight = 1/20000
        endStabWeight = 1/100
        homoWeight = 1/50000000
        heteroWeight = 1/200000000
        
        tempOpt = 60
        hairpinOpt = 40
        CGPOpt = 50
        endStabOpt = 3
        homoOpt = -2500
        heteroOpt = -2500
        
        
        def lstneqZero(value, level):
            for idx, val in enumerate(value):
                if val < level:
                    value[idx] = 0
                else:
                    value[idx] = val
            return value
        
        
        tempSum = np.array([])
        hairpinSum = np.array([])
        CGPSum = np.array([])
        endStabSum = np.array([])
        homoSum = np.array([])
        
        primer = self.primerF2
        indexNum = self.F2Best
        
        tempSum = np.append(tempSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_TM'])
        hairpinSum = np.append(hairpinSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_HAIRPIN_TH'])
        CGPSum = np.append(CGPSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_GC_PERCENT'])
        endStabSum = np.append(endStabSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_END_STABILITY'])
        homoSum = np.append(homoSum, primer3.calcHomodimer(primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE']).dg)
        
        primer = self.primerB2
        indexNum = self.B2Best
        
        tempSum = np.append(tempSum, primer['PRIMER_LEFT_' + str(indexNum) + '_TM'])
        hairpinSum = np.append(hairpinSum, primer['PRIMER_LEFT_' + str(indexNum) + '_HAIRPIN_TH'])
        CGPSum = np.append(CGPSum, primer['PRIMER_LEFT_' + str(indexNum) + '_GC_PERCENT'])
        endStabSum = np.append(endStabSum, primer['PRIMER_LEFT_' + str(indexNum) + '_END_STABILITY'])
        homoSum = np.append(homoSum, primer3.calcHomodimer(primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE']).dg)
        
        
        primer = self.primerF3
        indexNum = self.F3Best
        
        tempSum = np.append(tempSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_TM'])
        hairpinSum = np.append(hairpinSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_HAIRPIN_TH'])
        CGPSum = np.append(CGPSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_GC_PERCENT'])
        endStabSum = np.append(endStabSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_END_STABILITY'])
        homoSum = np.append(homoSum, primer3.calcHomodimer(primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE']).dg)
        
        
        
        primer = self.primerB3
        indexNum = self.B3Best
        
        tempSum = np.append(tempSum, primer['PRIMER_LEFT_' + str(indexNum) + '_TM'])
        hairpinSum = np.append(hairpinSum, primer['PRIMER_LEFT_' + str(indexNum) + '_HAIRPIN_TH'])
        CGPSum = np.append(CGPSum, primer['PRIMER_LEFT_' + str(indexNum) + '_GC_PERCENT'])
        endStabSum = np.append(endStabSum, primer['PRIMER_LEFT_' + str(indexNum) + '_END_STABILITY'])
        homoSum = np.append(homoSum, primer3.calcHomodimer(primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE']).dg)
        
        
        primer = self.primerF1c
        indexNum = self.F1cBest
        
        tempSum = np.append(tempSum, primer['PRIMER_LEFT_' + str(indexNum) + '_TM'])
        hairpinSum = np.append(hairpinSum, primer['PRIMER_LEFT_' + str(indexNum) + '_HAIRPIN_TH'])
        CGPSum = np.append(CGPSum, primer['PRIMER_LEFT_' + str(indexNum) + '_GC_PERCENT'])
        endStabSum = np.append(endStabSum, primer['PRIMER_LEFT_' + str(indexNum) + '_END_STABILITY'])
        homoSum = np.append(homoSum, primer3.calcHomodimer(primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE']).dg)
        
        
        primer = self.primerB1c
        indexNum = self.B1cBest
        
        tempSum = np.append(tempSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_TM'])
        hairpinSum = np.append(hairpinSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_HAIRPIN_TH'])
        CGPSum = np.append(CGPSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_GC_PERCENT'])
        endStabSum = np.append(endStabSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_END_STABILITY'])
        homoSum = np.append(homoSum, primer3.calcHomodimer(primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE']).dg)
        
        if(self.customLoopPrimers == True):
            # Needs implemteations
            pass
        else:
            primer = self.primerFLP
            indexNum = self.FLPBest
            
            tempSum = np.append(tempSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_TM'])
            hairpinSum = np.append(hairpinSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_HAIRPIN_TH'])
            CGPSum = np.append(CGPSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_GC_PERCENT'])
            endStabSum = np.append(endStabSum, primer['PRIMER_RIGHT_' + str(indexNum) + '_END_STABILITY'])
            homoSum = np.append(homoSum, primer3.calcHomodimer(primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE']).dg)
        
            primer = self.primerBLP
            indexNum = self.BLPBest
            
            tempSum = np.append(tempSum, primer['PRIMER_LEFT_' + str(indexNum) + '_TM'])
            hairpinSum = np.append(hairpinSum, primer['PRIMER_LEFT_' + str(indexNum) + '_HAIRPIN_TH'])
            CGPSum = np.append(CGPSum, primer['PRIMER_LEFT_' + str(indexNum) + '_GC_PERCENT'])
            endStabSum = np.append(endStabSum, primer['PRIMER_LEFT_' + str(indexNum) + '_END_STABILITY'])
            homoSum = np.append(homoSum, primer3.calcHomodimer(primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE']).dg)
            

        heteroSum = np.array([primer3.calcHeterodimer(str(self.FIP), str(self.BIP)).dg,
                              primer3.calcHeterodimer(str(self.FIP), self.primerF3['PRIMER_RIGHT_' + str(self.F3Best) + '_SEQUENCE']).dg,
                              primer3.calcHeterodimer(str(self.FIP), self.primerB3['PRIMER_LEFT_' + str(self.B3Best) + '_SEQUENCE']).dg,
                              primer3.calcHeterodimer(str(self.FIP), self.primerFLP['PRIMER_RIGHT_' + str(self.FLPBest) + '_SEQUENCE']).dg,
                              primer3.calcHeterodimer(str(self.FIP), self.primerBLP['PRIMER_LEFT_' + str(self.BLPBest) + '_SEQUENCE']).dg,
                              primer3.calcHeterodimer(str(self.BIP), self.primerF3['PRIMER_RIGHT_' + str(self.F3Best) + '_SEQUENCE']).dg,
                              primer3.calcHeterodimer(str(self.BIP), self.primerB3['PRIMER_LEFT_' + str(self.B3Best) + '_SEQUENCE']).dg,
                              primer3.calcHeterodimer(str(self.BIP), self.primerFLP['PRIMER_RIGHT_' + str(self.FLPBest) + '_SEQUENCE']).dg,
                              primer3.calcHeterodimer(str(self.BIP), self.primerBLP['PRIMER_LEFT_' + str(self.BLPBest) + '_SEQUENCE']).dg,
                              primer3.calcHeterodimer(self.primerFLP['PRIMER_RIGHT_' + str(self.FLPBest) + '_SEQUENCE'], self.primerBLP['PRIMER_LEFT_' + str(self.BLPBest) + '_SEQUENCE']).dg,
                              primer3.calcHeterodimer(self.primerFLP['PRIMER_RIGHT_' + str(self.FLPBest) + '_SEQUENCE'], self.primerF3['PRIMER_RIGHT_' + str(self.F3Best) + '_SEQUENCE']).dg,
                              primer3.calcHeterodimer(self.primerFLP['PRIMER_RIGHT_' + str(self.FLPBest) + '_SEQUENCE'], self.primerB3['PRIMER_LEFT_' + str(self.B3Best) + '_SEQUENCE']).dg,
                              primer3.calcHeterodimer(self.primerBLP['PRIMER_LEFT_' + str(self.BLPBest) + '_SEQUENCE'], self.primerF3['PRIMER_RIGHT_' + str(self.F3Best) + '_SEQUENCE']).dg,
                              primer3.calcHeterodimer(self.primerBLP['PRIMER_LEFT_' + str(self.BLPBest) + '_SEQUENCE'], self.primerB3['PRIMER_LEFT_' + str(self.B3Best) + '_SEQUENCE']).dg,
                              primer3.calcHeterodimer(self.primerF3['PRIMER_RIGHT_' + str(self.F3Best) + '_SEQUENCE'], self.primerB3['PRIMER_LEFT_' + str(self.B3Best) + '_SEQUENCE']).dg])

                
        total =  np.sum([tempWeight*np.sum((tempSum-tempOpt)**2),
                         hairpinWeight*np.sum((lstneqZero(hairpinSum,hairpinOpt)-hairpinOpt)**2),
                         GCPWeight*np.sum((CGPSum-CGPOpt)**2),
                         endStabWeight*np.sum((lstneqZero(endStabSum, endStabOpt)-endStabOpt)**2),
                         homoWeight*np.sum((lstneqZero(homoSum,homoOpt)-homoOpt)**2),
                         heteroWeight*np.sum((lstneqZero(heteroSum,heteroOpt)-heteroOpt)**2)])
    
#        print("tempSum: " + str(tempWeight*np.sum(tempSum**2)))
#        print("hairpinSum: " + str(hairpinWeight*np.sum(hairpinSum**2)))
#        print("CGPSum: " + str(GCPWeight*np.sum(CGPSum**2)))
#        print("endStabSum: " + str(endStabWeight*np.sum(endStabSum**2)))
#        print("homoSum: " + str(homoWeight*np.sum(homoSum**2)))
#        print("heteroSum: " + str(heteroWeight*np.sum(heteroSum**2)))
        
        #print("Total: " + str(total))
        
        return total
        
    def printAddress(self):
        print("Standard Frame Size: " + str(self.stdBuf))
        
        print("B3start Address: " + str(self.B3start))
        print("B3end Address: " + str(self.B3end))
        
        print("B2start Address: " + str(self.B2start))
        print("B2end Address: " + str(self.B2end))
        
        print("BLPstart Address: " + str(self.BLPstart))
        print("BLPend Address: " + str(self.BLPend))
        
        print("B1start Address: " + str(self.B1start))
        print("B1end Address: " + str(self.B1end))
        
        print("F1start Address: " + str(self.F1start))
        print("F1end Address: " + str(self.F1end))
        
        print("FLPstart Address: " + str(self.FLPstart))
        print("FLPend Address: " + str(self.FLPend))
        
        print("F2start Address: " + str(self.F2start))
        print("F2end Address: " + str(self.F2end))
        
        print("F3start Address: " + str(self.F3start))
        print("F3end Address: " + str(self.F3end))
        
    def genSinglePrimer(self, primerType, startAdd, endAdd):
        if(primerType==''):
            # ============ Find F2 ============
            pass
        elif(primerType=='B2'):
            # ============ Find B2 ============
            self.primerB2 = (primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_ID': 'B2',
                        'SEQUENCE_TEMPLATE': str(self.gene.seq[startAdd:endAdd]),
                        'SEQUENCE_EXCLUDED_REGION': [0, 0] 
                    },
                    {
                        'PRIMER_TASK': 'generic',
                        'PRIMER_PICK_LEFT_PRIMER': 1,
                        'PRIMER_PICK_INTERNAL_OLIGO': 0,
                        'PRIMER_PICK_RIGHT_PRIMER': 0,
                        'PRIMER_OPT_SIZE': self.innerLengthOpt,
                        'PRIMER_MIN_SIZE': self.innerLengthMin,
                        'PRIMER_MAX_SIZE': self.innerLengthMax,
                        'PRIMER_OPT_TM': self.innerTempOpt,
                        'PRIMER_MIN_TM': self.innerTempMin,
                        'PRIMER_MAX_TM': self.innerTempMax,
                        'PRIMER_MIN_GC': self.innerGCMin,
                        'PRIMER_MAX_GC': self.innerGCMax,
                        'PRIMER_MAX_POLY_X': 5,
                        'PRIMER_SALT_MONOVALENT': 50.0,
                        'PRIMER_DNA_CONC': 50.0,
                        'PRIMER_MAX_NS_ACCEPTED': 0,
                        'PRIMER_MAX_SELF_ANY': 12,
                        'PRIMER_MAX_SELF_END': 8,
                        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                        'PRIMER_PAIR_MAX_COMPL_END': 8,
                        'PRIMER_MAX_HAIRPIN_TH': self.maxHairpinTm,
                        'PRIMER_NUM_RETURN': self.numToReturn,}))
        
        elif(primerType==''):
            # ============ Find F3 ============
            pass
        elif(primerType==''):
            # ============ Find B3 ============
            pass
        elif(primerType==''):
            # ============ Find F1c ============
            pass
        elif(primerType=='B1c'):
            # ============ Find B1c ============
            self.primerB1c = (primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_ID': 'B1c',
                        'SEQUENCE_TEMPLATE': str(self.gene.seq[startAdd:endAdd]),
                        'SEQUENCE_EXCLUDED_REGION': [0, 0] 
                    },
                    {
                        'PRIMER_TASK': 'generic',
                        'PRIMER_PICK_LEFT_PRIMER': 0,
                        'PRIMER_PICK_INTERNAL_OLIGO': 0,
                        'PRIMER_PICK_RIGHT_PRIMER': 1,
                        'PRIMER_OPT_SIZE': self.complementLengthOpt,
                        'PRIMER_MIN_SIZE': self.complementLengthMin,
                        'PRIMER_MAX_SIZE': self.complementLengthMax,
                        'PRIMER_OPT_TM': self.complementTempOpt,
                        'PRIMER_MIN_TM': self.complementTempMin,
                        'PRIMER_MAX_TM': self.complementTempMax,
                        'PRIMER_MIN_GC': self.complementGCMin,
                        'PRIMER_MAX_GC': self.complementGCMax,
                        'PRIMER_MAX_POLY_X': 5,
                        'PRIMER_SALT_MONOVALENT': 50.0,
                        'PRIMER_DNA_CONC': 50.0,
                        'PRIMER_MAX_NS_ACCEPTED': 0,
                        'PRIMER_MAX_SELF_ANY': 12,
                        'PRIMER_MAX_SELF_END': 8,
                        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                        'PRIMER_PAIR_MAX_COMPL_END': 8,
                        'PRIMER_MAX_HAIRPIN_TH': self.maxHairpinTm,
                        'PRIMER_NUM_RETURN': self.numToReturn,}))
                    
        elif(primerType==''):
            # ============ Find FLP ============
            pass
        elif(primerType==''):
            # ============ Find BLP ============
            pass
               
    def something(self):
        primer = self.primerB2
        hairpinTm = [[i,primer['PRIMER_LEFT_' + str(i) + '_SELF_ANY_TH']] for i in range(0,primer['PRIMER_LEFT_NUM_RETURNED'])]
        hairpinTm = sorted(hairpinTm, key=lambda l:l[1])
        self.B2Best = hairpinTm[0][0]
        
    def genPrimers(self):   
        # ============ Find F2 ============
        self.primerF2 = (primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_ID': 'F2',
                        'SEQUENCE_TEMPLATE': str(self.gene.seq[self.F2start:self.F2end]),
                        'SEQUENCE_EXCLUDED_REGION': [0, 0] 
                    },
                    {
                        'PRIMER_TASK': 'generic',
                        'PRIMER_PICK_LEFT_PRIMER': 0,
                        'PRIMER_PICK_INTERNAL_OLIGO': 0,
                        'PRIMER_PICK_RIGHT_PRIMER': 1,
                        'PRIMER_OPT_SIZE': self.innerLengthOpt,
                        'PRIMER_MIN_SIZE': self.innerLengthMin,
                        'PRIMER_MAX_SIZE': self.innerLengthMax,
                        'PRIMER_OPT_TM': self.innerTempOpt,
                        'PRIMER_MIN_TM': self.innerTempMin,
                        'PRIMER_MAX_TM': self.innerTempMax,
                        'PRIMER_MIN_GC': self.innerGCMin,
                        'PRIMER_MAX_GC': self.innerGCMax,
                        'PRIMER_MAX_POLY_X': 5,
                        'PRIMER_SALT_MONOVALENT': 50.0,
                        'PRIMER_DNA_CONC': 50.0,
                        'PRIMER_MAX_NS_ACCEPTED': 0,
                        'PRIMER_MAX_SELF_ANY': 12,
                        'PRIMER_MAX_SELF_END': 8,
                        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                        'PRIMER_PAIR_MAX_COMPL_END': 8,
                        'PRIMER_MAX_HAIRPIN_TH': self.maxHairpinTm,
                        'PRIMER_NUM_RETURN': self.numToReturn,}))
        
        
        # ============ Find B2 ============
        self.primerB2 = (primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_ID': 'B2',
                        'SEQUENCE_TEMPLATE': str(self.gene.seq[self.B2start:self.B2end]),
                        'SEQUENCE_EXCLUDED_REGION': [0, 0] 
                    },
                    {
                        'PRIMER_TASK': 'generic',
                        'PRIMER_PICK_LEFT_PRIMER': 1,
                        'PRIMER_PICK_INTERNAL_OLIGO': 0,
                        'PRIMER_PICK_RIGHT_PRIMER': 0,
                        'PRIMER_OPT_SIZE': self.innerLengthOpt,
                        'PRIMER_MIN_SIZE': self.innerLengthMin,
                        'PRIMER_MAX_SIZE': self.innerLengthMax,
                        'PRIMER_OPT_TM': self.innerTempOpt,
                        'PRIMER_MIN_TM': self.innerTempMin,
                        'PRIMER_MAX_TM': self.innerTempMax,
                        'PRIMER_MIN_GC': self.innerGCMin,
                        'PRIMER_MAX_GC': self.innerGCMax,
                        'PRIMER_MAX_POLY_X': 5,
                        'PRIMER_SALT_MONOVALENT': 50.0,
                        'PRIMER_DNA_CONC': 50.0,
                        'PRIMER_MAX_NS_ACCEPTED': 0,
                        'PRIMER_MAX_SELF_ANY': 12,
                        'PRIMER_MAX_SELF_END': 8,
                        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                        'PRIMER_PAIR_MAX_COMPL_END': 8,
                        'PRIMER_MAX_HAIRPIN_TH': self.maxHairpinTm,
                        'PRIMER_NUM_RETURN': self.numToReturn,}))
        
        
        # ============ Find F3 ============
        self.primerF3 = (primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_ID': 'F3',
                        'SEQUENCE_TEMPLATE': str(self.gene.seq[self.F3start:self.F3end]),
                        'SEQUENCE_EXCLUDED_REGION': [0, 0] 
                    },
                    {
                        'PRIMER_TASK': 'generic',
                        'PRIMER_PICK_LEFT_PRIMER': 0,
                        'PRIMER_PICK_INTERNAL_OLIGO': 0,
                        'PRIMER_PICK_RIGHT_PRIMER': 1,
                        'PRIMER_OPT_SIZE': self.outerLengthOpt,
                        'PRIMER_MIN_SIZE': self.outerLengthMin,
                        'PRIMER_MAX_SIZE': self.outerLengthMax,
                        'PRIMER_OPT_TM': self.outerTempOpt,
                        'PRIMER_MIN_TM': self.outerTempMin,
                        'PRIMER_MAX_TM': self.outerTempMax,
                        'PRIMER_MIN_GC': self.outerGCMin,
                        'PRIMER_MAX_GC': self.outerGCMax,
                        'PRIMER_MAX_POLY_X': 5,
                        'PRIMER_SALT_MONOVALENT': 50.0,
                        'PRIMER_DNA_CONC': 50.0,
                        'PRIMER_MAX_NS_ACCEPTED': 0,
                        'PRIMER_MAX_SELF_ANY': 12,
                        'PRIMER_MAX_SELF_END': 8,
                        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                        'PRIMER_PAIR_MAX_COMPL_END': 8,
                        'PRIMER_MAX_HAIRPIN_TH': self.maxHairpinTm,
                        'PRIMER_NUM_RETURN': self.numToReturn,}))
        
        
        # ============ Find B3 ============
        self.primerB3 = (primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_ID': 'B3',
                        'SEQUENCE_TEMPLATE': str(self.gene.seq[self.B3start:self.B3end]),
                        'SEQUENCE_EXCLUDED_REGION': [0, 0] 
                    },
                    {
                        'PRIMER_TASK': 'generic',
                        'PRIMER_PICK_LEFT_PRIMER': 1,
                        'PRIMER_PICK_INTERNAL_OLIGO': 0,
                        'PRIMER_PICK_RIGHT_PRIMER': 0,
                        'PRIMER_OPT_SIZE': self.outerLengthOpt,
                        'PRIMER_MIN_SIZE': self.outerLengthMin,
                        'PRIMER_MAX_SIZE': self.outerLengthMax,
                        'PRIMER_OPT_TM': self.outerTempOpt,
                        'PRIMER_MIN_TM': self.outerTempMin,
                        'PRIMER_MAX_TM': self.outerTempMax,
                        'PRIMER_MIN_GC': self.outerGCMin,
                        'PRIMER_MAX_GC': self.outerGCMax,
                        'PRIMER_MAX_POLY_X': 5,
                        'PRIMER_SALT_MONOVALENT': 50.0,
                        'PRIMER_DNA_CONC': 50.0,
                        'PRIMER_MAX_NS_ACCEPTED': 0,
                        'PRIMER_MAX_SELF_ANY': 12,
                        'PRIMER_MAX_SELF_END': 8,
                        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                        'PRIMER_PAIR_MAX_COMPL_END': 8,
                        'PRIMER_MAX_HAIRPIN_TH': self.maxHairpinTm,
                        'PRIMER_NUM_RETURN': self.numToReturn,}))
        
        
        # ============ Find F1c ============
        self.primerF1c = (primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_ID': 'F1c',
                        'SEQUENCE_TEMPLATE': str(self.gene.seq[self.F1start:self.F1end]),
                        'SEQUENCE_EXCLUDED_REGION': [0, 0] 
                    },
                    {
                        'PRIMER_TASK': 'generic',
                        'PRIMER_PICK_LEFT_PRIMER': 1,
                        'PRIMER_PICK_INTERNAL_OLIGO': 0,
                        'PRIMER_PICK_RIGHT_PRIMER': 0,
                        'PRIMER_OPT_SIZE': self.complementLengthOpt,
                        'PRIMER_MIN_SIZE': self.complementLengthMin,
                        'PRIMER_MAX_SIZE': self.complementLengthMax,
                        'PRIMER_OPT_TM': self.complementTempOpt,
                        'PRIMER_MIN_TM': self.complementTempMin,
                        'PRIMER_MAX_TM': self.complementTempMax,
                        'PRIMER_MIN_GC': self.complementGCMin,
                        'PRIMER_MAX_GC': self.complementGCMax,
                        'PRIMER_MAX_POLY_X': 5,
                        'PRIMER_SALT_MONOVALENT': 50.0,
                        'PRIMER_DNA_CONC': 50.0,
                        'PRIMER_MAX_NS_ACCEPTED': 0,
                        'PRIMER_MAX_SELF_ANY': 12,
                        'PRIMER_MAX_SELF_END': 8,
                        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                        'PRIMER_PAIR_MAX_COMPL_END': 8,
                        'PRIMER_MAX_HAIRPIN_TH': self.maxHairpinTm,
                        'PRIMER_NUM_RETURN': self.numToReturn,}))
        
        
        # ============ Find B1c ============
        self.primerB1c = (primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_ID': 'B1c',
                        'SEQUENCE_TEMPLATE': str(self.gene.seq[self.B1start:self.B1end]),
                        'SEQUENCE_EXCLUDED_REGION': [0, 0] 
                    },
                    {
                        'PRIMER_TASK': 'generic',
                        'PRIMER_PICK_LEFT_PRIMER': 0,
                        'PRIMER_PICK_INTERNAL_OLIGO': 0,
                        'PRIMER_PICK_RIGHT_PRIMER': 1,
                        'PRIMER_OPT_SIZE': self.complementLengthOpt,
                        'PRIMER_MIN_SIZE': self.complementLengthMin,
                        'PRIMER_MAX_SIZE': self.complementLengthMax,
                        'PRIMER_OPT_TM': self.complementTempOpt,
                        'PRIMER_MIN_TM': self.complementTempMin,
                        'PRIMER_MAX_TM': self.complementTempMax,
                        'PRIMER_MIN_GC': self.complementGCMin,
                        'PRIMER_MAX_GC': self.complementGCMax,
                        'PRIMER_MAX_POLY_X': 5,
                        'PRIMER_SALT_MONOVALENT': 50.0,
                        'PRIMER_DNA_CONC': 50.0,
                        'PRIMER_MAX_NS_ACCEPTED': 0,
                        'PRIMER_MAX_SELF_ANY': 12,
                        'PRIMER_MAX_SELF_END': 8,
                        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                        'PRIMER_PAIR_MAX_COMPL_END': 8,
                        'PRIMER_MAX_HAIRPIN_TH': self.maxHairpinTm,
                        'PRIMER_NUM_RETURN': self.numToReturn,}))
        
        if(self.customLoopPrimers == False):
            # ============ Find Forward Loop Primer ============
            self.primerFLP = (primer3.bindings.designPrimers(
                        {
                            'SEQUENCE_ID': 'F3',
                            'SEQUENCE_TEMPLATE': str(self.gene.seq[self.FLPstart:self.FLPend]),
                            'SEQUENCE_EXCLUDED_REGION': [0, 0] 
                        },
                        {
                            'PRIMER_TASK': 'generic',
                            'PRIMER_PICK_LEFT_PRIMER': 0,
                            'PRIMER_PICK_INTERNAL_OLIGO': 0,
                            'PRIMER_PICK_RIGHT_PRIMER': 1,
                            'PRIMER_OPT_SIZE': self.loopLengthOpt,
                            'PRIMER_MIN_SIZE': self.loopLengthMin,
                            'PRIMER_MAX_SIZE': self.loopLengthMax,
                            'PRIMER_OPT_TM': self.loopTempOpt,
                            'PRIMER_MIN_TM': self.loopTempMin,
                            'PRIMER_MAX_TM': self.loopTempMax,
                            'PRIMER_MIN_GC': self.loopGCMin,
                            'PRIMER_MAX_GC': self.loopGCMax,
                            'PRIMER_MAX_POLY_X': 5,
                            'PRIMER_SALT_MONOVALENT': 50.0,
                            'PRIMER_DNA_CONC': 50.0,
                            'PRIMER_MAX_NS_ACCEPTED': 0,
                            'PRIMER_MAX_SELF_ANY': 12,
                            'PRIMER_MAX_SELF_END': 8,
                            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                            'PRIMER_PAIR_MAX_COMPL_END': 8,
                            'PRIMER_MAX_HAIRPIN_TH': self.maxHairpinTm,
                            'PRIMER_NUM_RETURN': self.numToReturn,}))
            
            # ============ Find Backward Loop Primer ============
            self.primerBLP = (primer3.bindings.designPrimers(
                        {
                            'SEQUENCE_ID': 'F3',
                            'SEQUENCE_TEMPLATE': str(self.gene.seq[self.BLPstart:self.BLPend]),
                            'SEQUENCE_EXCLUDED_REGION': [0, 0] 
                        },
                        {
                            'PRIMER_TASK': 'generic',
                            'PRIMER_PICK_LEFT_PRIMER': 1,
                            'PRIMER_PICK_INTERNAL_OLIGO': 0,
                            'PRIMER_PICK_RIGHT_PRIMER': 0,
                            'PRIMER_OPT_SIZE': self.loopLengthOpt,
                            'PRIMER_MIN_SIZE': self.loopLengthMin,
                            'PRIMER_MAX_SIZE': self.loopLengthMax,
                            'PRIMER_OPT_TM': self.loopTempOpt,
                            'PRIMER_MIN_TM': self.loopTempMin,
                            'PRIMER_MAX_TM': self.loopTempMax,
                            'PRIMER_MIN_GC': self.loopGCMin,
                            'PRIMER_MAX_GC': self.loopGCMax,
                            'PRIMER_MAX_POLY_X': 5,
                            'PRIMER_SALT_MONOVALENT': 50.0,
                            'PRIMER_DNA_CONC': 50.0,
                            'PRIMER_MAX_NS_ACCEPTED': 0,
                            'PRIMER_MAX_SELF_ANY': 12,
                            'PRIMER_MAX_SELF_END': 8,
                            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                            'PRIMER_PAIR_MAX_COMPL_END': 8,
                            'PRIMER_MAX_HAIRPIN_TH': self.maxHairpinTm,
                            'PRIMER_NUM_RETURN': self.numToReturn,}))
            
        
        
            
        
        
        
    def setCustomLoopPrimer(self, sequenceF, sequenceB):
        self.customLoopPrimerF = Seq(sequenceF)
        self.customLoopPrimerB = Seq(sequenceB)
    
    
    def createXIP(self):       
        if(self.customLoopPrimers == True):
            # ============ Create FIP ============
            self.FIP = Seq(self.primerF1c['PRIMER_LEFT_' + str(self.F1cBest) + '_SEQUENCE']) + Seq(''.join([ 'T' for i in range(0,self.polyTinterConnect)])) + self.customLoopPrimerF + Seq(''.join([ 'T' for i in range(0,self.polyTinterConnect)])) + Seq(self.primerF2['PRIMER_RIGHT_' + str(self.F2Best) + '_SEQUENCE'])
            
            # ============ Create BIP ============
            self.BIP = Seq(self.primerB1c['PRIMER_RIGHT_' + str(self.B1cBest) + '_SEQUENCE']) + Seq(''.join([ 'T' for i in range(0,self.polyTinterConnect)])) + self.customLoopPrimerB + Seq(''.join([ 'T' for i in range(0,self.polyTinterConnect)])) + Seq(self.primerB2['PRIMER_LEFT_' + str(self.B2Best) + '_SEQUENCE'])
        else:
            # ============ Create FIP ============
            self.FIP = Seq(self.primerF1c['PRIMER_LEFT_' + str(self.F1cBest) + '_SEQUENCE']) + Seq(''.join([ 'T' for i in range(0,self.polyTinterConnect)])) + Seq(self.primerF2['PRIMER_RIGHT_' + str(self.F2Best) + '_SEQUENCE'])
            
            # ============ Create BIP ============
            self.BIP = Seq(self.primerB1c['PRIMER_RIGHT_' + str(self.B1cBest) + '_SEQUENCE']) + Seq(''.join([ 'T' for i in range(0,self.polyTinterConnect)])) + Seq(self.primerB2['PRIMER_LEFT_' + str(self.B2Best) + '_SEQUENCE'])
            
            
    def printPrimerComponentStats(self):
        primer = self.primerF2
        indexNum = self.F2Best
        
        print('F2 primer: ' + primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE'])
        print('     |---> Tm = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_TM']) + 'C')
        print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_HAIRPIN_TH']) + 'C')
        print('     |---> CG% = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_GC_PERCENT']) + '%')
        print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + str(primer['PRIMER_RIGHT_' + str(indexNum) + '_END_STABILITY']) + ' kcal/mol')
        print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE']).dg) + ' kcal/mol')
        
        primer = self.primerB2
        indexNum = self.B2Best
        
        print('B2 primer: ' + primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE'])
        print('     |---> Tm = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_TM']) + 'C')
        print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_HAIRPIN_TH']) + 'C')
        print('     |---> CG% = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_GC_PERCENT']) + '%')
        print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + str(primer['PRIMER_LEFT_' + str(indexNum) + '_END_STABILITY']) + ' kcal/mol')
        print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE']).dg) + ' kcal/mol')
        
        
        primer = self.primerF3
        indexNum = self.F3Best
        
        print('F3 primer: ' + primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE'])
        print('     |---> Tm = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_TM']) + 'C')
        print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_HAIRPIN_TH']) + 'C')
        print('     |---> CG% = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_GC_PERCENT']) + '%')
        print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + str(primer['PRIMER_RIGHT_' + str(indexNum) + '_END_STABILITY']) + ' kcal/mol')
        print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE']).dg) + ' kcal/mol')
        
        
        
        primer = self.primerB3
        indexNum = self.B3Best
        
        print('B3 primer: ' + primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE'])
        print('     |---> Tm = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_TM']) + 'C')
        print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_HAIRPIN_TH']) + 'C')
        print('     |---> CG% = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_GC_PERCENT']) + '%')
        print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + str(primer['PRIMER_LEFT_' + str(indexNum) + '_END_STABILITY']) + ' kcal/mol')
        print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE']).dg) + ' kcal/mol')
        
        
        primer = self.primerF1c
        indexNum = self.F1cBest
        
        print('F1c primer: ' + primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE'])
        print('     |---> Tm = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_TM']) + 'C')
        print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_HAIRPIN_TH']) + 'C')
        print('     |---> CG% = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_GC_PERCENT']) + '%')
        print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + str(primer['PRIMER_LEFT_' + str(indexNum) + '_END_STABILITY']) + ' kcal/mol')
        print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE']).dg) + ' kcal/mol')
        
        
        primer = self.primerB1c
        indexNum = self.B1cBest
        
        print('B1c primer: ' + primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE'])
        print('     |---> Tm = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_TM']) + 'C')
        print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_HAIRPIN_TH']) + 'C')
        print('     |---> CG% = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_GC_PERCENT']) + '%')
        print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + str(primer['PRIMER_RIGHT_' + str(indexNum) + '_END_STABILITY']) + ' kcal/mol')
        print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE']).dg) + ' kcal/mol')
        
        if(self.customLoopPrimers == True):
            primer = self.customLoopPrimerF.complement()
            print('F Loop Complement primer: ' + str(primer))
            print('     |---> Tm = ' + '{:02.2f}'.format(primer3.calcTm(str(primer))) + 'C')
            print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer3.calcHairpin(str(primer)).dg) + 'C')
            print('     |---> CG% = ' + '{:02.2f}'.format(SeqUtils.GC(primer)) + '%')
            #print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + '{:02.2f}'.format(primer3.calcEndStability(str(primer))) + ' kcal/mol')
            print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(str(primer)).dg) + ' kcal/mol')
            
            primer = self.customLoopPrimerB.complement()
            print('B Loop Complement primer: ' + str(primer))
            print('     |---> Tm = ' + '{:02.2f}'.format(primer3.calcTm(str(primer))) + 'C')
            print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer3.calcHairpin(str(primer)).dg) + 'C')
            print('     |---> CG% = ' + '{:02.2f}'.format(SeqUtils.GC(primer)) + '%')
            #print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + '{:02.2f}'.format(primer3.calcEndStability(str(primer))) + ' kcal/mol')
            print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(str(primer)).dg) + ' kcal/mol')
            
            primer = self.customLoopPrimerF
            print('F Loop primer: ' + str(primer))
            print('     |---> Tm = ' + '{:02.2f}'.format(primer3.calcTm(str(primer))) + 'C')
            print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer3.calcHairpin(str(primer)).dg) + 'C')
            print('     |---> CG% = ' + '{:02.2f}'.format(SeqUtils.GC(primer)) + '%')
            #print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + '{:02.2f}'.format(primer3.calcEndStability(str(primer))) + ' kcal/mol')
            print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(str(primer)).dg) + ' kcal/mol')
            
            primer = self.customLoopPrimerB
            print('B Loop primer: ' + str(primer))
            print('     |---> Tm = ' + '{:02.2f}'.format(primer3.calcTm(str(primer))) + 'C')
            print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer3.calcHairpin(str(primer)).dg) + 'C')
            print('     |---> CG% = ' + '{:02.2f}'.format(SeqUtils.GC(primer)) + '%')
            #print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + '{:02.2f}'.format(primer3.calcEndStability(str(primer))) + ' kcal/mol')
            print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(str(primer)).dg) + ' kcal/mol')
        else:
            primer = self.primerFLP
            indexNum = self.FLPBest
            
            print('FLP primer: ' + primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE'])
            print('     |---> Tm = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_TM']) + 'C')
            print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_HAIRPIN_TH']) + 'C')
            print('     |---> CG% = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_GC_PERCENT']) + '%')
            print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + str(primer['PRIMER_RIGHT_' + str(indexNum) + '_END_STABILITY']) + ' kcal/mol')
            print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE']).dg) + ' kcal/mol')
        
            primer = self.primerBLP
            indexNum = self.BLPBest
            
            print('BLP primer: ' + primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE'])
            print('     |---> Tm = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_TM']) + 'C')
            print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_HAIRPIN_TH']) + 'C')
            print('     |---> CG% = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_GC_PERCENT']) + '%')
            print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + str(primer['PRIMER_LEFT_' + str(indexNum) + '_END_STABILITY']) + ' kcal/mol')
            print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE']).dg) + ' kcal/mol')
            
            
            

    
    def printPrimerStats(self):
        primer = self.primerF3
        indexNum = self.F3Best
        
        print('F3 primer: ' + primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE'])
        print('     |---> Tm = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_TM']) + 'C')
        print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_HAIRPIN_TH']) + 'C')
        print('     |---> CG% = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_GC_PERCENT']) + '%')
        print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + str(primer['PRIMER_RIGHT_' + str(indexNum) + '_END_STABILITY']) + ' kcal/mol')
        print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE']).dg) + ' kcal/mol')
        
        primer = self.primerB3
        indexNum = self.B3Best
        
        print('B3 primer: ' + primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE'])
        print('     |---> Tm = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_TM']) + 'C')
        print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_HAIRPIN_TH']) + 'C')
        print('     |---> CG% = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_GC_PERCENT']) + '%')
        print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + str(primer['PRIMER_LEFT_' + str(indexNum) + '_END_STABILITY']) + ' kcal/mol')
        print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE']).dg) + ' kcal/mol')
        
        if(self.customLoopPrimers == True):
            primer = self.customLoopPrimerF.complement()
            print('F Loop Complement primer: ' + str(primer))
            print('     |---> Tm = ' + '{:02.2f}'.format(primer3.calcTm(str(primer))) + 'C')
            print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer3.calcHairpin(str(primer)).dg) + 'C')
            print('     |---> CG% = ' + '{:02.2f}'.format(SeqUtils.GC(primer)) + '%')
            #print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + '{:02.2f}'.format(primer3.calcEndStability(str(primer))) + ' kcal/mol')
            print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(str(primer)).dg) + ' kcal/mol')
            
            primer = self.customLoopPrimerB.complement()
            print('B Loop Complement primer: ' + str(primer))
            print('     |---> Tm = ' + '{:02.2f}'.format(primer3.calcTm(str(primer))) + 'C')
            print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer3.calcHairpin(str(primer)).dg) + 'C')
            print('     |---> CG% = ' + '{:02.2f}'.format(SeqUtils.GC(primer)) + '%')
            #print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + '{:02.2f}'.format(primer3.calcEndStability(str(primer))) + ' kcal/mol')
            print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(str(primer)).dg) + ' kcal/mol')
            
            primer = self.FIP
            print('FIP primer: ' + str(primer))
            print('     |---> Tm = ' + '{:02.2f}'.format(primer3.calcTm(str(primer))) + 'C')
            #print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer3.calcHairpin(str(primer)).dg) + 'C')
            print('     |---> CG% = ' + '{:02.2f}'.format(SeqUtils.GC(primer)) + '%')
            #print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + '{:02.2f}'.format(primer3.calcEndStability(str(primer))) + ' kcal/mol')
            #print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(str(primer)).dg) + ' kcal/mol')
            
            primer = self.BIP
            print('BIP primer: ' + str(primer))
            print('     |---> Tm = ' + '{:02.2f}'.format(primer3.calcTm(str(primer))) + 'C')
            #print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer3.calcHairpin(str(primer)).dg) + 'C')
            print('     |---> CG% = ' + '{:02.2f}'.format(SeqUtils.GC(primer)) + '%')
            #print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + '{:02.2f}'.format(primer3.calcEndStability(str(primer))) + ' kcal/mol')
            #print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(str(primer)).dg) + ' kcal/mol')
        else:
            primer = self.primerFLP
            indexNum = self.FLPBest
            
            print('FLP primer: ' + primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE'])
            print('     |---> Tm = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_TM']) + 'C')
            print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_HAIRPIN_TH']) + 'C')
            print('     |---> CG% = ' + '{:02.2f}'.format(primer['PRIMER_RIGHT_' + str(indexNum) + '_GC_PERCENT']) + '%')
            print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + str(primer['PRIMER_RIGHT_' + str(indexNum) + '_END_STABILITY']) + ' kcal/mol')
            print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(primer['PRIMER_RIGHT_' + str(indexNum) + '_SEQUENCE']).dg) + ' kcal/mol')
            
            primer = self.primerBLP
            indexNum = self.BLPBest
            
            print('BLP primer: ' + primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE'])
            print('     |---> Tm = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_TM']) + 'C')
            print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_HAIRPIN_TH']) + 'C')
            print('     |---> CG% = ' + '{:02.2f}'.format(primer['PRIMER_LEFT_' + str(indexNum) + '_GC_PERCENT']) + '%')
            print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + str(primer['PRIMER_LEFT_' + str(indexNum) + '_END_STABILITY']) + ' kcal/mol')
            print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(primer['PRIMER_LEFT_' + str(indexNum) + '_SEQUENCE']).dg) + ' kcal/mol')
        
            
            primer = self.FIP
            print('FIP primer: ' + str(primer))
            print('     |---> Tm = ' + '{:02.2f}'.format(primer3.calcTm(str(primer))) + 'C')
            print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer3.calcHairpin(str(primer)).dg) + 'C')
            print('     |---> CG% = ' + '{:02.2f}'.format(SeqUtils.GC(primer)) + '%')
            #print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + '{:02.2f}'.format(primer3.calcEndStability(str(primer))) + ' kcal/mol')
            print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(str(primer)).dg) + ' kcal/mol')
            
            primer = self.BIP
            print('BIP primer: ' + str(primer))
            print('     |---> Tm = ' + '{:02.2f}'.format(primer3.calcTm(str(primer))) + 'C')
            print('     |---> Hairpin Tm = ' + '{:02.2f}'.format(primer3.calcHairpin(str(primer)).dg) + 'C')
            print('     |---> CG% = ' + '{:02.2f}'.format(SeqUtils.GC(primer)) + '%')
            #print('     |---> Stablity of last 5 nucleotides of 3\' end = ' + '{:02.2f}'.format(primer3.calcEndStability(str(primer))) + ' kcal/mol')
            print('     |---> Homodimer Formation Energy = ' + '{:02.2f}'.format(primer3.calcHomodimer(str(primer)).dg) + ' kcal/mol')
            
        
    
    def printPrimerPotStats(self): 
        print('=============================================')
        print('Homodimer Formation Energy Analysis for FIP/BIP:')
        print('FIP: ' + '{:02.2f}'.format(primer3.calcHomodimer(str(self.FIP)).dg) + ' kcal/mol')
        print('BIP: ' + '{:02.2f}'.format(primer3.calcHomodimer(str(self.BIP)).dg) + ' kcal/mol')
        
        print('=============================================')
        print('Heterodimer Formation Energy Analysis:')
        print('FIP-BIP: ' + '{:02.2f}'.format(primer3.calcHeterodimer(str(self.FIP), str(self.BIP)).dg) + ' kcal/mol')
        print('FIP-F3: ' + '{:02.2f}'.format(primer3.calcHeterodimer(str(self.FIP), self.primerF3['PRIMER_RIGHT_' + str(self.F3Best) + '_SEQUENCE']).dg ) + ' kcal/mol')
        print('FIP-B3: ' + '{:02.2f}'.format(primer3.calcHeterodimer(str(self.FIP), self.primerB3['PRIMER_LEFT_' + str(self.B3Best) + '_SEQUENCE']).dg ) + ' kcal/mol')
        print('FIP-FLP: ' + '{:02.2f}'.format(primer3.calcHeterodimer(str(self.FIP), self.primerFLP['PRIMER_RIGHT_' + str(self.FLPBest) + '_SEQUENCE']).dg ) + ' kcal/mol')
        print('FIP-BLP: ' + '{:02.2f}'.format(primer3.calcHeterodimer(str(self.FIP), self.primerBLP['PRIMER_LEFT_' + str(self.BLPBest) + '_SEQUENCE']).dg ) + ' kcal/mol')
        print('BIP-F3: ' + '{:02.2f}'.format(primer3.calcHeterodimer(str(self.BIP), self.primerF3['PRIMER_RIGHT_' + str(self.F3Best) + '_SEQUENCE']).dg ) + ' kcal/mol')
        print('BIP-B3: ' + '{:02.2f}'.format(primer3.calcHeterodimer(str(self.BIP), self.primerB3['PRIMER_LEFT_' + str(self.B3Best) + '_SEQUENCE']).dg ) + ' kcal/mol')
        print('BIP-FLP: ' + '{:02.2f}'.format(primer3.calcHeterodimer(str(self.BIP), self.primerFLP['PRIMER_RIGHT_' + str(self.FLPBest) + '_SEQUENCE']).dg ) + ' kcal/mol')
        print('BIP-BLP: ' + '{:02.2f}'.format(primer3.calcHeterodimer(str(self.BIP), self.primerBLP['PRIMER_LEFT_' + str(self.BLPBest) + '_SEQUENCE']).dg ) + ' kcal/mol')
        print('FLP-BLP: ' + '{:02.2f}'.format(primer3.calcHeterodimer(self.primerFLP['PRIMER_RIGHT_' + str(self.FLPBest) + '_SEQUENCE'], self.primerBLP['PRIMER_LEFT_' + str(self.BLPBest) + '_SEQUENCE']).dg ) + ' kcal/mol')
        print('FLP-F3: ' + '{:02.2f}'.format(primer3.calcHeterodimer(self.primerFLP['PRIMER_RIGHT_' + str(self.FLPBest) + '_SEQUENCE'], self.primerF3['PRIMER_RIGHT_' + str(self.F3Best) + '_SEQUENCE']).dg ) + ' kcal/mol')
        print('FLP-B3: ' + '{:02.2f}'.format(primer3.calcHeterodimer(self.primerFLP['PRIMER_RIGHT_' + str(self.FLPBest) + '_SEQUENCE'], self.primerB3['PRIMER_LEFT_' + str(self.B3Best) + '_SEQUENCE']).dg ) + ' kcal/mol')
        print('BLP-F3: ' + '{:02.2f}'.format(primer3.calcHeterodimer(self.primerBLP['PRIMER_LEFT_' + str(self.BLPBest) + '_SEQUENCE'], self.primerF3['PRIMER_RIGHT_' + str(self.F3Best) + '_SEQUENCE']).dg ) + ' kcal/mol')
        print('BLP-B3: ' + '{:02.2f}'.format(primer3.calcHeterodimer(self.primerBLP['PRIMER_LEFT_' + str(self.BLPBest) + '_SEQUENCE'], self.primerB3['PRIMER_LEFT_' + str(self.B3Best) + '_SEQUENCE']).dg ) + ' kcal/mol')
        print('F3-B3: ' + '{:02.2f}'.format(primer3.calcHeterodimer(self.primerF3['PRIMER_RIGHT_' + str(self.F3Best) + '_SEQUENCE'], self.primerB3['PRIMER_LEFT_' + str(self.B3Best) + '_SEQUENCE']).dg ) + ' kcal/mol')
        
        
        print('=============================================')

if __name__ == '__main__':
    myFinder = MiLAMPFinder()
    
    myFinder.loadSeqFile("BRCA1Thomas.txt")
    
    myFinder.setStart(43124000)
    
    myFinder.setMethSite([43125347, 43125364, 43125372, 43125377, 43125383, 43125409, 43125411, 43125419])
##    
##    check = myFinder.methSiteCheck()
##    
##    checkList = myFinder.checkSites
##    
    myFinder.setResEnz()
    
    myFinder.numOfMethSitesUsed(-2)
    
    myFinder.findResEnz()
    
    rezEnz = myFinder.getCurrentEnz()
    
    myFinder.setEnz(0)
    
    myFinder.setDefaultPrimerParam()
    
    mySettings = dict()
    
    mySettings['polyTinterConnect'] = 3
#    
#    myFinder.setPrimerParam(mySettings)
#    
#    print(myFinder.polyTinterConnect)
#    
#    myFinder.genPrimers()
#    
#    #myFinder.setCustomLoopPrimer('CCCCCCCCCCCCCCCCCCCCC', 'GGGGGGGGGGGGGGGGGGGGGG')
#    
#    myFinder.createXIP()
    
#    print("============================================================================")
#    #myFinder.printPrimerComponentStats()
#    print("============================================================================")
#    #myFinder.printPrimerStats()
#    print("============================================================================")
#    #myFinder.printPrimerPotStats()
#    print("============================================================================")
#    t = time.time()
#    myFinder.generateRandomPrimers(n=50)
#    elapsed = time.time() - t
#    print(elapsed)
##    print(myFinder.rateCurrentPrimerSet())
#    #primerRating  cutSiteList
#    fmt='%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %f'
#    
#    np.savetxt("ratedPrimers2.csv", np.transpose(np.append(myFinder.cutSiteList, [myFinder.primerRating], axis=0)), delimiter=",",fmt=fmt)
#    
    # =========== data creatation ===========
    def randomGene(length):
       DNA=""
       for count in range(length):
          DNA+=choice("CGTA")
       return DNA
   
    
    myFinder.setStart = 0
    # Only use DNA sequence of of 500 basepairs
    seqLen = 600
    
    # Force methylation site to be in center
    myFinder.methSite = seqLen/2
    

    
    # Number of sequence used in the study
    numOfSeq = 750
    
    # Number of primer used for each sequence in the study
    numOfPrimers = 100
    
    seqList = []
       
    for i in range(numOfSeq):
        seqList.append(randomGene(seqLen))
    
    print('Processing...')
    
    cutHeaders = ['B3start',
                    'B3end',
                    'B2start',
                    'B2end',
                    'BLPstart',
                    'BLPend',
                    'B1start',
                    'B1end',
                    'F1start',
                    'F1end',
                    'FLPstart',
                    'FLPend',
                    'F2start',
                    'F2end',
                    'F3start',
                    'F3end',
                    'Rating']
    
#    data = pd.DataFrame(columns = cutHeaders.append('Sequence'))
    dataList = []
    t = time.time()
    for randInx in range(numOfSeq):
        myFinder.loadSeqText(seqList[randInx])
        myFinder.recalculateCutSites()
        print('Sequence: ' + str(randInx))
        myFinder.generateRandomPrimers(n=numOfPrimers)
        dataBuf = pd.DataFrame(np.transpose(np.append(myFinder.cutSiteList, [myFinder.primerRating], axis=0)), columns=cutHeaders)
        dataBuf['Sequence'] = seqList[randInx]*numOfPrimers
        dataList.append(dataBuf)
#        print(dataBuf)
        #myFinder.cutSiteList
        #seqList[randInx]
        
#        print(elapsed)#    data.to_csv(path_or_buf='ratedPrimersData.csv',index=False)
    elapsed = time.time() - t
    print(elapsed)
    
    
    data = pd.DataFrame(columns = cutHeaders.append('Sequence'))
    data = data.append(dataList, ignore_index=True)
    
    data['BoolRating'] = data['Rating']
    
    data.loc[data.notnull()['BoolRating'],'BoolRating'] = True
    
    data.loc[data.isnull()['BoolRating'],'BoolRating'] = False
    
    data['BoolRating'] = data['BoolRating'].astype('bool')
    
    data.to_pickle('../../UCLA Feed Thru/MiLAMPFinderData/ratedPrimersData9')
    
    print('Successfully created primers: {a:6.2f}%'.format(a=(data['Rating'].shape[0]-data.isnull().sum()['Rating'])/data['Rating'].shape[0]*100))
    
#    data = pd.read_pickle('ratedPrimersData3')
    
#    print(data2)





























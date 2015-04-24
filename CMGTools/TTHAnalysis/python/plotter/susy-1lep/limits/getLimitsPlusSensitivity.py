#!/usr/bin/env python

import shutil
import subprocess
import os


# check CMSSW version
cmssw_vers = os.environ['CMSSW_VERSION']

if 'CMSSW_7_1' not in cmssw_vers:
    print "You've got the wrong CMSSW env:"
    print "Need CMSSW_7_1_X with limit tool!"
    exit(0)

cardDirectory="susy_cards_1l_4fb"
cardDirectory = os.path.abspath(cardDirectory)

limitDir = cardDirectory+"/limits/"
if not os.path.isdir(limitDir):
    os.mkdir(limitDir)

print 'Entering out dir', limitDir
os.chdir(limitDir)
print

Samples = []
Samples.append("T1tttt_HL_1500_100")
#Samples.append("T1ttbbWW_1300_300_290")
#Samples.append("T1ttbbWW_1300_300_295")
Samples.append("T1tttt_HM_1200_800")
#Samples.append("T5qqqqWW_Gl1200_Chi1000_LSP800")
#Samples.append("T5ttttDeg_1300_300_280")



#VariantSnippet = ["standardnJ_HTTTYesNo","standardnJ_HT1","finenJ_HTTTYesNo","finenJ_HT1"]
#VariantSnippet = ["standardnJ_HTLowLepPt","standardnJ_HTLowLepPtDPhi","standardnJ_HT1","finenJ_HTLowLepPt","finenJ_HTLowLepPtDPhi","finenJ_HT1"]
#VariantSnippet = ["standardnJ_HTLowLepPt","standardnJ_HT1","standardnJ_HighLowLepPt","finenJ_HTLowLepPt","finenJ_HT1","finenJ_HighLowLepPt"]
#VariantSnippet = ["standardnJ_HighLowLepPt","finenJ_HighLowLepPt"]
#VariantSnippet = ["standardnJ", "finenJ","finenJ_HT1","standardnJ_HT1"]
#VariantSnippet = ["finenJ_HT1","standardnJ_HT1","finenJ_HTTop","standardnJ_HTTop"]
VariantSnippet = ["standardnJ", "finenJ"]


#standard: baseline
#finenJ: baseline plus highnJ-bin
#finenJ_HT*: only high HT-bin
#"finenJ_HT1": baselin, only high HT
#"finenJ_HT1TT": baseline+ toptag-counting
#"finenJ_HTTop": baseline+ dphi>.05 single topness >1.25, topness > 5
#"finenJ_HTTopTT: baseline+ toptag-counting + dphi>.05 single topness >1.25, topness > 5
#"finenJ_HighLowLepPt" combination of hard/soft leptons

limitdict = {}
sigdict = {}
for s_i, s in enumerate(Samples):
    for v_i, v in enumerate(VariantSnippet):

        print 80*'#'
        print  "Limits_"+s+"_"+v+".txt"

        os.system("combine -M Asymptotic "+cardDirectory+"/"+s+"/CnC2015X_"+v+".card.txt -n Limits_"+s+"_"+v+" &> Limits_"+s+"_"+v+".txt")
        os.system("combine -M ProfileLikelihood --significance "+cardDirectory+"/"+s+"/CnC2015X_"+v+".card.txt -t -1 --expectSignal=1 -n Significance_"+s+"_"+v+" &> Significance_"+s+"_"+v+".txt")

        lf=open("Limits_"+s+"_"+v+".txt")
        llines=lf.readlines()

        sf=open("Significance_"+s+"_"+v+".txt")
        slines=sf.readlines()

        limitline = llines[6]
        print limitline#,"bla0"
#        print limitline.split("<")[1],"bla1"
#        snipp=  limitline.split("<")[1]
#        print snipp, "bla11"
#        print (limitline.split("<")[1]).split()[0],"bla2"
        limitdict[(s,v)]= (limitline.split("<")[1]).split()[0]
        sigline= slines[2]
        print sigline





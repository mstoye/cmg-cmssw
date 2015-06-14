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
#Samples.append("T2tt_425_325")
#Samples.append("T2tt_650_325")
#Samples.append("T2tt_850_100")



#VariantSnippet = ["standardnJ_HTTTYesNo","standardnJ_HT1","finenJ_HTTTYesNo","finenJ_HT1"]
#VariantSnippet = ["standardnJ_HTLowLepPt","standardnJ_HTLowLepPtDPhi","standardnJ_HT1","finenJ_HTLowLepPt","finenJ_HTLowLepPtDPhi","finenJ_HT1"]
#VariantSnippet = ["standardnJ_HTLowLepPt","standardnJ_HT1","standardnJ_HighLowLepPt","finenJ_HTLowLepPt","finenJ_HT1","finenJ_HighLowLepPt"]
#VariantSnippet = ["standardnJ_HighLowLepPt","finenJ_HighLowLepPt"]
#VariantSnippet = ["standardnJ", "finenJ","finenJ_HT1","standardnJ_HT1"]
#VariantSnippet = ["finenJ_HT1","standardnJ_HT1","finenJ_HTTop","standardnJ_HTTop"]
#VariantSnippet = ["standardnJ", "finenJ","allnJ", "allnJTopnesses"]

VariantSnippet = ["NonExtreme_standardnJ", "NonExtreme_68nJ", "NonExtreme_68nJTopness", "NonExtreme_68nJDPhi05", "NonExtreme_68nJSingleTopness"]


VariantSnippet = []

#for ST in ["ST4"]:
for ST in ["ST0", "ST1", "ST2", "ST3", "ST4"]:
#    for nJ in ["68j", "6Infj", "9Infj"]:
    for nJ in ["45j"]:
        for nB in ["2B", "3p"]:
#        for nB in ["2B"]:
            for HT in ["HT0", "HT1"]:
#            for HT in ["HT0"]:
                for RD in ["DPhi00", "DPhi05", "DPhi10", "Stop", "Top"]:
                    VariantSnippet.append(nB+"_"+ST+"_"+nJ+"_"+HT+"_"+RD)

#VariantSnippet = ["2B_ST0_9Infj_HT0_Top"]


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
        print cardDirectory+"/"+s+"/CnC2015X_"+v+".card.txt"
        print  "Limits_"+s+"_"+v+".txt"

        os.system("combine -M Asymptotic "+cardDirectory+"/"+s+"/CnC2015X_"+v+".card.txt -n Limits_"+s+"_"+v+" &> Limits_"+s+"_"+v+".txt")
        os.system("combine -M ProfileLikelihood --significance "+cardDirectory+"/"+s+"/CnC2015X_"+v+".card.txt -t -1 --expectSignal=1 -n Significance_"+s+"_"+v+" &> Significance_"+s+"_"+v+".txt")

        lf=open("Limits_"+s+"_"+v+".txt")
        llines=lf.readlines()

        sf=open("Significance_"+s+"_"+v+".txt")
        slines=sf.readlines()

        RunSuccessfully = False
        if len(llines)>5 and len(slines)>1:
            RunSuccessfully = True
            if "You must have at least one signal process (id <= 0)" in llines[8]:
                print llines[8]
                RunSuccessfully= False
        if RunSuccessfully:
            limitline = llines[6]
            limitdict[(s,v)]= (limitline.split("<")[1]).split()[0]
            print limitdict[(s,v)]
            sigline= slines[2]
            sigdict[(s,v)]= (sigline.split(":")[1]).split()[0]
            print sigdict[(s,v)]
        else: 
            print "WARNING: Output files not in the correct format..."
            for line in llines: print line,
            for line in slines: print line,
            open("Failed_"+s+"_"+v+".txt", 'a').close()




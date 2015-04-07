from CMGTools.TTHAnalysis.treeReAnalyzer import *
from ROOT import TLorentzVector, TVector2
import ROOT

def mt_2(p4one, p4two):
    return sqrt(2*p4one.Pt()*p4two.Pt()*(1-cos(p4one.Phi()-p4two.Phi())))


class EventVars1LGenLevel:
    def __init__(self):
        self.branches = ["GenDeltaPhiLepWSum", "GenDeltaPhiLepWDirect", "GenWSumMass", "GenWDirectMass",
                         "nidxGenWs", "GenmTLepNu"
                         ]
    
    def listBranches(self):
        return self.branches[:]
    
    def __call__(self,event):
        # The following variables still need to be double-checked for validity
        genleps = [l for l in Collection(event,"genLep","ngenLep")] 
        genParts = [l for l in Collection(event,"GenPart","nGenPart")] 
        ngenleps = len(genleps); ngenParts = len (genParts)

        GenDeltaPhiLepWSum=-999
        GenDeltaPhiLepWDirect=-999
        GenWSumMass=-999
        GenWDirectMass=-999
        GenmTLepNu=-999
        idx_genWs=[]
        idx_genLeps=[]
        idx_genNus=[]
        # find gen-level neutrinos (status 23), calculate deltaPhi (lep, nu), and genW-masses m(lep+nu)
        # for this: start from genleps (status 23)
        for i_lep, genlep in enumerate(genleps):
            if genlep.status == 23 and abs(genlep.motherId) == 24: # genlep is outgoing and has W as mother
                W_idx = genlep.motherIndex
                idx_genWs.append(W_idx)
                idx_genLeps.append(i_lep)
                for i_nu, genPart in enumerate(genParts):
                    if genPart.motherIndex==W_idx and genPart.status == 23: # find W as mother
                        if abs(genPart.pdgId) == 12 or abs(genPart.pdgId) == 14 or abs(genPart.pdgId) == 16: #check whether it is a neutrino
                            idx_genNus.append(i_nu)
  
    
        if(len(idx_genLeps)>=1):
            genlepP4 = genleps[idx_genLeps[0]].p4()
            genNuP4 = genParts[idx_genNus[0]].p4()
            genWSumP4 = genlepP4 + genNuP4
            genWDirectP4 = genParts[genleps[idx_genLeps[0]].motherIndex].p4()
            GenDeltaPhiLepWSum = genlepP4.DeltaPhi(genWSumP4)
            GenDeltaPhiLepWDirect = genlepP4.DeltaPhi(genWDirectP4)
            GenWSumMass = genWSumP4.M()
            GenWDirectMass = genWDirectP4.M()
            GenmTLepNu = mt_2(genlepP4,genNuP4)


        ret    =  { 'GenDeltaPhiLepWSum'   : GenDeltaPhiLepWSum } #initialize the dictionary with a first entry   
        ret["GenDeltaPhiLepWDirect"] = GenDeltaPhiLepWDirect
        ret["GenWSumMass"]           = GenWSumMass
        ret["GenWDirectMass"]        = GenWDirectMass
        ret["GenmTLepNu"]            = GenmTLepNu

        ret["nidxGenWs"]= len(idx_genWs)
        
        return ret

if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("tree") 
    class Tester(Module):
        def __init__(self, name):
            Module.__init__(self,name,None)
            self.sf = EventVars1LGenLevel()
        def analyze(self,ev):
            print "\nrun %6d lumi %4d event %d: leps %d" % (ev.run, ev.lumi, ev.evt, ev.nLepGood)
            print self.sf(ev)
    el = EventLoop([ Tester("tester") ])
    el.loop([tree], maxEvents = 50)

        


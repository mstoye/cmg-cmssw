from CMGTools.TTHAnalysis.treeReAnalyzer import *
from ROOT import TLorentzVector, TVector2, std
import ROOT
import time
import itertools
import PhysicsTools.Heppy.loadlibs
import array
import operator

ROOT.gInterpreter.GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector") #need this to be able to use topness code

mt2wSNT = ROOT.heppy.mt2w_bisect.mt2w()
topness = ROOT.Topness.Topness()

def getPhysObjectArray(j): # https://github.com/HephySusySW/Workspace/blob/72X-master/RA4Analysis/python/mt2w.py
    px = j.pt*cos(j.phi )
    py = j.pt*sin(j.phi )
    pz = j.pt*sinh(j.eta )
    E = sqrt(px*px+py*py+pz*pz) #assuming massless particles...
    return array.array('d', [E, px, py,pz])

class EventVars1L_extra:
    def __init__(self):
        self.branches = [ "MT2W", "Topness" ]

    def listBranches(self):
        return self.branches[:]

    def __call__(self,event,base):

        # output dict:
        ret = {}

        # get some collections from initial tree
        leps = [l for l in Collection(event,"LepGood","nLepGood")]
        jets = [j for j in Collection(event,"Jet","nJet")]

        njet = len(jets); nlep = len(leps)

        # MET
        metp4 = ROOT.TLorentzVector(0,0,0,0)
        metp4.SetPtEtaPhiM(event.met_pt,event.met_eta,event.met_phi,event.met_mass)
        pmiss  =array.array('d',[event.met_pt * cos(event.met_phi), event.met_pt * sin(event.met_phi)] )

        ####################################
        # import output from previous step #
        #base = keyvals
        ####################################

        # get selected leptons
        tightLeps = []
        tightLepsIdx = base['tightLepsIdx']
        tightLeps = [leps[idx] for idx in tightLepsIdx]
        nTightLeps = len(tightLeps)

        # get selected jets
        centralJet30 = []
        centralJet30idx = base['centralJet30idx']
        centralJet30 = [jets[idx] for idx in centralJet30idx]
        nCentralJet30 = len(centralJet30)

        # B jets
        BJetCMVAMedium30 = []
        BJetCMVAMedium30idx = base['BJetCMVAMedium30idx']
        NonBJetCMVAMedium30 = []
        nBJetCMVAMedium30 = base['nBJetCMVAMedium30']

        for idx,jet in enumerate(jets):
            if idx in BJetCMVAMedium30idx:
                BJetCMVAMedium30.append(jet)
            else:
                NonBJetCMVAMedium30.append(jet)

        ##################################################################
        # The following variables need to be double-checked for validity #
        ##################################################################

        #add topness and mt2W-variable (timing issue with topness: slows down the friend tree production by a factor of ~3)
        mt2w_values=[]

        if nTightLeps>=1: #topness and mt2w only make sense for
            lep = getPhysObjectArray(tightLeps[0])
            if nBJetCMVAMedium30==0 and nCentralJet30>=3: #All combinations from the highest three light (or b-) jets
                consideredJets = [ getPhysObjectArray(jet) for jet in NonBJetCMVAMedium30[:3] ] # only throw arrays into the permutation business
                ftPerms = itertools.permutations(consideredJets, 2)
                for perm in ftPerms:
                    mt2wSNT.set_momenta(lep, perm[0], perm[1], pmiss)
                    mt2w_values.append(mt2wSNT.get_mt2w())
            elif nBJetCMVAMedium30==1 and nCentralJet30>=2: #All combinations from one b and the highest two light jets
                consideredJets = [ getPhysObjectArray(jet) for jet in NonBJetCMVAMedium30[:2] ] # only throw arrays into the permutation business
                consideredJets.append(getPhysObjectArray(BJetCMVAMedium30[0]))
                ftPerms = itertools.permutations(consideredJets, 2)
                for perm in ftPerms:
                    mt2wSNT.set_momenta(lep, perm[0], perm[1], pmiss)
                    mt2w_values.append(mt2wSNT.get_mt2w())
            elif nBJetCMVAMedium30==2:
                consideredJets = [ getPhysObjectArray(jet) for jet in BJetCMVAMedium30[:2] ] # only throw arrays into the permutation business
                ftPerms = itertools.permutations(consideredJets, 2)
                for perm in ftPerms:
                    mt2wSNT.set_momenta(lep, perm[0], perm[1], pmiss)
                    mt2w_values.append(mt2wSNT.get_mt2w())
            elif nBJetCMVAMedium30>=3: #All combinations from the highest three b jets
                consideredJets = [ getPhysObjectArray(jet) for jet in BJetCMVAMedium30[:3] ] # only throw arrays into the permutation business
                ftPerms = itertools.permutations(consideredJets, 2)
                for perm in ftPerms:
                    mt2wSNT.set_momenta(lep, perm[0], perm[1], pmiss)
                    mt2w_values.append(mt2wSNT.get_mt2w())

        if len(mt2w_values)>0:
            ret["MT2W"]=min(mt2w_values)
        #else:
        #    ret["MT2W"]=-999

        # topness
        #ret['Topness']=-999
        if nTightLeps>=1: #topness and mt2w only make sense for
            # does not seem to work for njet =3 ??! # need to edit btag working point in the code...!! did not quickly find a twiki with official phys14 cmva working points
            if (nCentralJet30>=3) and (nBJetCMVAMedium30>=1) :

                p4_jets = std.vector(TLorentzVector)();
                bdisc_jets = std.vector('float')();

                for jet in centralJet30:
                    jetTLorentz = ROOT.TLorentzVector(0,0,0,0)
                    jetTLorentz.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass)
                    p4_jets.push_back(jetTLorentz)
                    bdisc_jets.push_back(jet.btagCMVA)

                lepTLorentz = ROOT.TLorentzVector(0,0,0,0)
                lepTLorentz.SetPtEtaPhiM(tightLeps[0].pt, tightLeps[0].eta, tightLeps[0].phi, tightLeps[0].mass)

                # calc topness
                tempTopness = topness.GetTopness(p4_jets,bdisc_jets,lepTLorentz,metp4) #this is really slow!
                if tempTopness <=0:
                    print tempTopness, "this will fail"
                else:
                    ret['Topness'] = log(tempTopness) #this is really slow!

        # return branches
        return ret

if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("tree")
    class Tester(Module):
        def __init__(self, name):
            Module.__init__(self,name,None)
            self.sf = EventVars1L()
        def analyze(self,ev):
            print "\nrun %6d lumi %4d event %d: leps %d" % (ev.run, ev.lumi, ev.evt, ev.nLepGood)
            print self.sf(ev)
    el = EventLoop([ Tester("tester") ])
    el.loop([tree], maxEvents = 50)

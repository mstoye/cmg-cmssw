AlCaRecoMatrix = {'ExpressCosmics' : 'SiStripCalZeroBias+TkAlCosmics0T',
                  'StreamExpress'  : 'SiStripCalZeroBias+TkAlMinBias+MuAlCalIsolatedMu+DtCalib',
                  'MinimumBias'    : 'SiStripCalMinBias+TkAlMinBias',
                  'Commissioning'  : 'HcalCalIsoTrk',
                  'SingleMu'       : 'MuAlCalIsolatedMu+MuAlOverlaps+TkAlMuonIsolated+DtCalib',
                  'DoubleMu'       : 'MuAlCalIsolatedMu+MuAlOverlaps+DtCalib+TkAlZMuMu',
                  'MuOnia'         : 'TkAlJpsiMuMu+TkAlUpsilonMuMu',
                  'SingleElectron' : 'EcalCalElectron',
                  'DoubleElectron' : 'EcalCalElectron',
                  'Cosmics'        : 'TkAlCosmics0T+MuAlGlobalCosmics+HcalCalHOCosmics+DtCalibCosmics',
                  # These two cannot run on RAW, they are just meant to run on the dedicated AlcaRAW so they do not enter the allForPrompt list
                  'AlCaP0'         : 'EcalCalPi0Calib+EcalCalEtaCalib',
                  # ---------------------------------------------------------------------------------------------------------------------------
                  'HcalNZS'        : 'HcalCalMinBias'
                  # This is in the AlCaRecoMatrix, but no RelVals are produced
                  # 'TestEnablesTracker' : 'TkAlLAS'
                  }

def buildList(pdList, matrix):
    """Takes a list of primary datasets (PDs) and the AlCaRecoMatrix (a dictinary) and returns a string with all the AlCaRecos for the selected PDs separated by the '+' character without duplicates."""
    alCaRecoList = []
    for pd in pdList:
        alCaRecoList.extend(matrix[pd].split("+"))
    # remove duplicates converting to a set
    alCaRecoList = set(alCaRecoList)
    stringList = ''
    for alCaReco in alCaRecoList:
        if stringList == '':
            stringList += alCaReco
        else:
            stringList += '+'+alCaReco
    return stringList

# Update the lists anytime a new PD is added to the matrix
autoAlca = { 'allForPrompt'         : buildList(['MinimumBias', 'Commissioning', 'SingleMu', 'DoubleMu', 'MuOnia', 'SingleElectron', 'DoubleElectron', 'HcalNZS'], AlCaRecoMatrix),
             'allForExpress'        : buildList(['StreamExpress'], AlCaRecoMatrix),
             'allForPromptCosmics'  : buildList(['Cosmics'], AlCaRecoMatrix),
             'allForExpressCosmics' : buildList(['ExpressCosmics'], AlCaRecoMatrix) }
autoAlca.update(AlCaRecoMatrix)

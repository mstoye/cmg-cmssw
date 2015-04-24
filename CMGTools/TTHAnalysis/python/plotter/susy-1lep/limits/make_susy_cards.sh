#!/bin/bash

if [[ "$1" == "SingleLepAFS" ]]; then
    shift # shift register
    T="/afs/cern.ch/work/k/kirschen/public/PlotExampleSamples/V3";
    FT="/afs/cern.ch/work/k/kirschen/public/PlotExampleSamples/V3/PHYS14_V3_Friend_Batool"
    J=4;
elif [[ "$HOSTNAME" == *"lxplus"* ]] ; then
    T="/afs/cern.ch/work/k/kirschen/public/PlotExampleSamples/V3";
    FT="/afs/cern.ch/work/k/kirschen/public/PlotExampleSamples/V3/PHYS14_V3_Friend_Batool"
    J=4;
elif [[ "$1" == "DESYV3" ]] ; then
    shift # shift register
    T="/nfs/dust/cms/group/susy-desy/Run2/MC/CMGtuples/Phys14_v3/ForCMGplot";
    FT="/nfs/dust/cms/group/susy-desy/Run2/MC/CMGtuples/Phys14_v3/Phys14_V3_Friend_CSVbtag"
    J=8;
elif [[ "$HOSTNAME" == *"naf"* ]] ; then
    T="/nfs/dust/cms/group/susy-desy/Run2/MC/CMGtuples/Phys14_v3/ForCMGplot";
    FT="/nfs/dust/cms/group/susy-desy/Run2/MC/CMGtuples/Phys14_v3/Phys14_V3_Friend_CSVbtag"
    J=8;
else
    echo "Didn't specify location!"
    echo "Usage: ./make_cards.sh location analysis "
    exit 0
fi

LUMI=4.0
OUTDIR="susy_cards_1l_4fb"
OPTIONS=" -P $T -j $J -l $LUMI -f --s2v --tree treeProducerSusySingleLepton --od $OUTDIR --asimov "

# Get current plotter dir
#PLOTDIR="$CMSSW_BASE/src/CMGTools/TTHAnalysis/python/plotter/"
PLOTDIR=$(pwd -P)
PLOTDIR=${PLOTDIR/plotter/plotterX}
PLOTDIR=$(echo $PLOTDIR |  cut -d 'X' -f 1 )

# Append FriendTree dir
OPTIONS=" $OPTIONS -F sf/t $FT/evVarFriend_{cname}.root "

function makeCard_1l {
    local EXPR=$1; local BINS=$2; local SYSTS=$3; local OUT=$4; local GO=$5

    CutFlowCard="1l_CardsFullCutFlow.txt"

    # b-jet cuts
    case $nB in
        0B)  GO="${GO} -R 1nB 0nB nBJetMedium30==0 " ;;
        1B)  GO="${GO} -R 1nB 1nB nBJetMedium30==1 " ;;
        2B)  GO="${GO} -R 1nB 2nB nBJetMedium30==2 " ;;
        2Btop)  GO="${GO} -R 1nB 2nB nBJetMedium30==2&&Topness>5 " ;;
        3p)  GO="${GO} -R 1nB 3nBp nBJetMedium30>=3 " ;;
    esac;

    # ST categories
    case $ST in
        ST0)  GO="${GO} -R st200 st200250 ST>200&&ST<250 " ;;
        ST1)  GO="${GO} -R st200 st250350 ST>250&&ST<350 " ;;
        ST2)  GO="${GO} -R st200 st350450 ST>350&&ST<450 " ;;
        ST3)  GO="${GO} -R st200 st450600 ST>450&&ST<600 " ;;
        ST4)  GO="${GO} -R st200 st600Inf ST>600 " ;;
    esac;

    # jet multiplicities
    case $nJ in
        45j)  GO="${GO} -R geq6j 45j nCentralJet30>=4&&nCentralJet30<=5"  ;;
        68j)  GO="${GO} -R geq6j 67j nCentralJet30>=6&&nCentralJet30<=8"  ;;
        6Infj)  GO="${GO} -R geq6j geq6j nCentralJet30>=6"  ;;
        9Infj)  GO="${GO} -R geq6j geq8j nCentralJet30>=9"  ;;
        68TTj)  GO="${GO} -R geq6j 68TTj nCentralJet30+2*nHighPtTopTagPlusTau23>=6&&nCentralJet30+2*nHighPtTopTagPlusTau23<9"  ;;
        9InfTTj)  GO="${GO} -R geq6j 9InfTTj nCentralJet30+2*nHighPtTopTagPlusTau23>=9"  ;;
    esac;

    # HT and "R&D" categories
    case $HT in
        HT0) GO="${GO} -R ht500 ht5001000 HT>500&&HT<1000"  ;;
        HT1) GO="${GO} -R ht500 ht1000Inf HT>=1000"  ;;
        HTDPhi) GO="${GO} -R ht500 ht1000Inf HT>=1000 -R dp1 dp05 fabs(DeltaPhiLepW)>0.5 "  ;;
        HTStop) GO="${GO} -R ht500 ht1000Inf HT>=1000 -R dp1 dp05 fabs(DeltaPhiLepW)>0.5 -A dp1 stopness (TopVarsMETovTopMin[0]-0.5)/0.5+(TopVarsMtopMin[0]-175)/175>1.25"  ;;
        HTTop) GO="${GO} -R ht500 ht1000Inf HT>=1000 -R dp1 dp05 fabs(DeltaPhiLepW)>0.5 -A dp1 stopness (TopVarsMETovTopMin[0]-0.5)/0.5+(TopVarsMtopMin[0]-175)/175>1.25&&Topness>5"  ;;
        HTLowLepPt) GO="${GO} -R ht500 ht1000Inf HT>=1000 -R 1tl 1tllowpt nTightLeps==1&&LepGood1_pt<=25  -R dp1 dp00 fabs(DeltaPhiLepW)>0.0 -A dp1 stopness (TopVarsMETovTopMin[0]-0.5)/0.5+(TopVarsMtopMin[0]-175)/175>1.25&&Topness>5"  ;;
        HTLowLepPtDPhi) GO="${GO} -R ht500 ht1000Inf HT>=1000 -R 1tl 1tllowpt nTightLeps==1&&LepGood1_pt<=25"  ;;
        HTTTYes) GO="${GO} -R ht500 ht1000Inf HT>=1000&&nHighPtTopTagPlusTau23>=1"  ;;
        HTTTNo) GO="${GO} -R ht500 ht1000Inf HT>=1000&&nHighPtTopTagPlusTau23==0"  ;;
    esac;

    if [[ "$PRETEND" == "1" ]]; then
        echo "making datacard $OUT from makeShapeCardsSusy.py mca-Phys14_1l.txt $CutFlowCard \"$EXPR\" \"$BINS\" $SYSTS $GO --dummyYieldsForZeroBkg;"
    else
        echo "making datacard $OUT from makeShapeCardsSusy.py mca-Phys14_1l.txt $CutFlowCard \"$EXPR\" \"$BINS\" $SYSTS $GO --dummyYieldsForZeroBkg;"
        python $PLOTDIR/makeShapeCardsSusy.py $PLOTDIR/mca-Phys14_1l.txt $PLOTDIR/susy-1lep/$CutFlowCard "$EXPR" "$BINS" $SYSTS -o $OUT $GO --dummyYieldsForZeroBkg;
        echo "  -- done at $(date)";
    fi;
}

function combineCardsSmart {

    DummyC=0
    AllC=0
    CMD=""
    for C in $*; do
        # missing datacards
        test -f $C || continue;

        if grep -q "DummyContent" $C; then
            echo "empty bin ${C}" >&2
            DummyC=$(($DummyC+1))
            if grep -q "observation 0.0$" $C; then
                echo "this is not the way it was intended..."
            fi
        fi

        grep -q "observation 0.0$" $C && continue
        AllC=$((AllC+1))
        CMD="${CMD} $(basename $C .card.txt)=$C ";
    done
    if [[ "$CMD" == "" ]]; then
        echo "Not any card found in $*" 1>&2 ;
    else
#       echo "combineCards.py $CMD" >&2
        combineCards.py $CMD
    fi
    if [[ "$DummyC" != "0" ]]; then
        echo "In total $DummyC out of $AllC are empty, but taken into account by adding DummyContent." >&2
    fi
    #echo "In total $DummyC out of $AllC are empty, but taken into account by adding DummyContent." >&2
}

if [[ "$1" == "--pretend" ]]; then
    PRETEND=1; shift;
    echo "# Pretending to run"
fi;

if [[ "$1" == "1l-makeCards" ]]; then

    SYSTS="syst/susyDummy.txt"
    CnC_expr="1" #not used as of now
    CnC_bins="[0.5,1.5]"


    echo "Making individual datacards"
    for ST in ST0 ST1 ST2 ST3 ST4; do for nJ in 68j 6Infj 9Infj; do for nB in 2B 3p; do for HT in HT0 HT1; do
#    for ST in ST0 ST1 ST2 ST3 ST4; do
        echo " --- CnC2015X_${nB}_${ST}_${nJ}_${HT} ---"
        makeCard_1l $CnC_expr $CnC_bins $SYSTS CnC2015X_${nB}_${ST}_${nJ}_${HT} "$OPTIONS";
                done; done; done; done
#    done;
    #exit

fi

if [[ "$1" == "1l-combine" ]]; then

    if [[ ! $CMSSW_VERSION == *"CMSSW_7_1_"* ]] ;then
	echo "You don't have the correct CMSSW environment!"
	echo "Found: $CMSSW_VERSION, need CMSSW_7_1_X"
	exit 0
    fi

    echo "Making combined datacards"

    for D in $OUTDIR/T[0-9]*; do
        test -f $D/CnC2015X_2B_ST0_68j_HT0.card.txt || continue
        (cd $D && echo "    $D";
            for nB in 2B 3p; do
                combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_6Infj_{HT0,HT1}.card.txt >  CnC2015X_${nB}_standardnJ.card.txt
                combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_{68j,9Infj}_{HT0,HT1}.card.txt >  CnC2015X_${nB}_finenJ.card.txt

            done
            combineCardsSmart CnC2015X_{2B,3p}_standardnJ.card.txt >  CnC2015X_standardnJ.card.txt # standard nJet-binning; HT-binning
            combineCardsSmart CnC2015X_{2B,3p}_finenJ.card.txt >  CnC2015X_finenJ.card.txt #fine nJet-binning; HT-binning

        )
    done
fi

echo "Done at $(date)";

source /vols/grid/cms/setup.sh

# 2HDMC
#thdmccmssw="10_2_19"
#scramv1 project CMSSW CMSSW_${thdmccmssw}
#pushd CMSSW_${thdmccmssw}/src
#eval `scram runtime -sh`
#curl -o 2HDMC-1.8.0.tar.gz  https://2hdmc.hepforge.org/downloads/?f=2HDMC-1.8.0.tar.gz
#tar -xvf 2HDMC-1.8.0.tar.gz
#popd
#
#cp setup_tools/TestPointVaryingm12_2.cpp CMSSW_${thdmccmssw}/src/2HDMC-1.8.0/src/ 
#sed -i '/PROG/s/$/ TestPointVaryingm12_2/' CMSSW_${thdmccmssw}/src/2HDMC-1.8.0/Makefile
#
#pushd CMSSW_${thdmccmssw}/src/2HDMC-1.8.0
#eval `scram runtime -sh`
#make
#popd

## 2HDECAY
#hdecaycmssw="10_2_19"
#scramv1 project CMSSW CMSSW_${hdecaycmssw}
#pushd CMSSW_${hdecaycmssw}/src
#eval `scram runtime -sh`
#source /vols/grid/cms/setup.sh
#git clone git@github.com:marcel-krause/2HDECAY.git
#popd
#
#pushd CMSSW_${hdecaycmssw}/src/2HDECAY
#eval `scram runtime -sh`
#curl -O https://feynarts.de/looptools/LoopTools-2.14.tar.gz
#tar -xvf LoopTools-2.14.tar.gz
#python setup.py
#popd

# higgstools
htoolscmssw="10_2_19"
pushd CMSSW_${htoolscmssw}/src
eval `scram runtime -sh`
git clone https://gitlab.com/higgsbounds/higgstools.git
git clone https://gitlab.com/higgsbounds/hbdataset.git
git clone https://gitlab.com/higgsbounds/hsdataset.git
popd

pushd CMSSW_${htoolscmssw}/src/higgstools
pip3 install .
popd

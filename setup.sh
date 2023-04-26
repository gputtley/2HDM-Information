source /vols/grid/cms/setup.sh

# 2HDMC
scramv1 project CMSSW CMSSW_12_4_8
pushd CMSSW_12_4_8/src
eval `scram runtime -sh`
curl -o 2HDMC-1.8.0.tar.gz  https://2hdmc.hepforge.org/downloads/?f=2HDMC-1.8.0.tar.gz
tar -xvf 2HDMC-1.8.0.tar.gz
popd

cp setup_tools/TestPointVaryingm12_2.cpp CMSSW_12_4_8/src/2HDMC-1.8.0/src/ 
sed -i '/PROG/s/$/ TestPointVaryingm12_2/' CMSSW_12_4_8/src/2HDMC-1.8.0/Makefile

pushd CMSSW_12_4_8/src/2HDMC-1.8.0
eval `scram runtime -sh`
make
popd

# 2HDECAY
scramv1 project CMSSW CMSSW_10_2_19
pushd CMSSW_10_2_19/src
eval `scram runtime -sh`
source /vols/grid/cms/setup.sh
git clone git@github.com:marcel-krause/2HDECAY.git
popd

pushd CMSSW_10_2_19/src/2HDECAY
eval `scram runtime -sh`
curl -O https://feynarts.de/looptools/LoopTools-2.14.tar.gz
tar -xvf LoopTools-2.14.tar.gz
python setup.py
popd

# higgstools
pushd CMSSW_12_4_8/src
eval `scram runtime -sh`
git clone https://gitlab.com/higgsbounds/higgstools.git
git clone https://gitlab.com/higgsbounds/hbdataset.git
git clone https://gitlab.com/higgsbounds/hsdataset.git
popd

pushd CMSSW_12_4_8/src/higgstools
pip3 install .
popd

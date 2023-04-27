import os
import subprocess
import itertools
import math as m
import sys

def IterateAndCheckTheoreticalConstraintsAndFindm12_2(type_2hdms=["X"],mhs=[125.0],mHs=[200.0],mAs=[160.0],mHcs=[200.0],tanbs=[10.0,20.0,30.0,40.0,50.0],sinbmas=[1.0]):
  cmd_list = [
              "pushd CMSSW_12_4_8/src/2HDMC-1.8.0/ > /dev/null",
              "source /vols/grid/cms/setup.sh",
              "eval `scram runtime -sh`"]
  for (type_2hdm,mh,mH,mA,mHc,tanb,sinbma) in itertools.product(type_2hdms,mhs,mHs,mAs,mHcs,tanbs,sinbmas):
    name = MakeName(type_2hdm=type_2hdm,mh=mh,mH=mH,mA=mA,mHc=mHc,tanb=tanb,sinbma=sinbma)
    cmd_list.append("./TestPointVaryingm12_2 --mh=%(mh)s --mH=%(mH)s --mA=%(mA)s --mC=%(mHc)s --tb=%(tanb)s --sba=%(sinbma)s &> ../../../output/2HDMC_output_%(name)s.txt" % vars())
  cmd_list.append("popd > /dev/null"),
  os.system(";".join(cmd_list))

  valid_dict = {}
  m12_2_dict = {}
  for (type_2hdm,mh,mH,mA,mHc,tanb,sinbma) in itertools.product(type_2hdms,mhs,mHs,mAs,mHcs,tanbs,sinbmas):
    name = MakeName(type_2hdm=type_2hdm,mh=mh,mH=mH,mA=mA,mHc=mHc,tanb=tanb,sinbma=sinbma)
    valid_dict[name] = str(subprocess.check_output(['tail', '-1', "output/2HDMC_output_%(name)s.txt" % vars()]).decode("utf-8").rstrip()) == "Valid"
    m12_2_dict[name] = float(subprocess.check_output(['tail', '-8', "output/2HDMC_output_%(name)s.txt" % vars()]).decode("utf-8").split()[1]) if valid_dict[name] else None

  return valid_dict, m12_2_dict

def MakeName(type_2hdm="X",mh=125.0,mH=200.0,mA=160.0,mHc=200.0,tanb=20.0,sinbma=1.0):
  return "t"+type_2hdm+"_mh"+str(mh)+"_mH"+str(mH)+"_mA"+str(mA)+"_mHc"+str(mHc)+"_tanb"+str(tanb)+"_sinbma"+str(sinbma)

def WriteListToFile(lst,output_name):
  textfile = open(output_name, "w")
  for i in lst:
    textfile.write(i + "\n")
  textfile.close()

def IterateAndFindWidthsAndBranchingRatios(type_2hdms=["X"],mhs=[125.0],mHs=[200.0],mAs=[160.0],mHcs=[200.0],tanbs=[10.0,20.0,30.0,40.0,50.0],sinbmas=[1.0],m12_2_dict={},renscheme=7,input_name=""):


  if os.path.isdir("CMSSW_10_2_19/src/2HDECAY/Input{}".format(input_name)):
    os.system("rm -r CMSSW_10_2_19/src/2HDECAY/Input{}".format(input_name))
  os.system("mkdir CMSSW_10_2_19/src/2HDECAY/Input{}".format(input_name))

  width_dict = {}
  br_dicts = {}

  for (type_2hdm,mh,mH,mA,mHc,tanb,sinbma) in itertools.product(type_2hdms,mhs,mHs,mAs,mHcs,tanbs,sinbmas):
    name = MakeName(type_2hdm=type_2hdm,mh=mh,mH=mH,mA=mA,mHc=mHc,tanb=tanb,sinbma=sinbma)

    if name not in m12_2_dict: continue

    if type_2hdm == "X":
      num_2hdm = "3"

    file_list = [
                 "SLHAIN   = 1",
                 "SLHAOUT  = 1",
                 "COUPVAR  = 1",
                 "HIGGS    = 5",
                 "OMIT ELW = 1",
                 "OMIT ELW2= 0",
                 "SM4      = 0",
                 "FERMPHOB = 0",
                 "2HDM     = 1",
                 "MODEL    = 1",
                 "TGBET    = 5.07403e+01",
                 "MABEG    = 4.67967e+02",
                 "MAEND    = 4.67967e+02",
                 "NMA      = 1",
                 "********************* hMSSM (MODEL = 10) *********************************",
                 "MHL      = 125.D0",
                 "**************************************************************************",
                 "ALS(MZ)  = 1.18000e-01",
                 "MSBAR(2) = 9.50000e-02",
                 "MCBAR(3) = 0.98600e+00",
                 "MBBAR(MB)= 4.18000e+00",
                 "MT       = 1.73200e+02",
                 "MTAU     = 1.77682e+00",
                 "MMUON    = 1.056583715e-01",
                 "1/ALPHA  = 1.37036e+02",
                 "ALPHAMZ  = 7.754222173973729e-03",
                 "GF       = 1.1663787e-05",
                 "GFCALC   = 0.000000000",
                 "GAMW     = 2.08500e+00",
                 "GAMZ     = 2.49520e+00",
                 "MZ       = 9.11876e+01",
                 "MW       = 8.0385e+01",
                 "VTB      = 9.9910e-01",
                 "VTS      = 4.040e-02",
                 "VTD      = 8.67e-03",
                 "VCB      = 4.12e-02",
                 "VCS      = 9.7344e-01",
                 "VCD      = 2.252e-01",
                 "VUB      = 3.51e-03",
                 "VUS      = 2.2534e-01",
                 "VUD      = 9.7427e-01",
                 "********************* 4TH GENERATION *************************************",
                 "  SCENARIO FOR ELW. CORRECTIONS TO H -> GG (EVERYTHING IN GEV):",
                 "  GG_ELW = 1: MTP = 500    MBP = 450    MNUP = 375    MEP = 450",
                 "  GG_ELW = 2: MBP = MNUP = MEP = 600    MTP = MBP+50*(1+LOG(M_H/115)/5)",
                 "",
                 "GG_ELW   = 1",
                 "MTP      = 500.D0",
                 "MBP      = 450.D0",
                 "MNUP     = 375.D0",
                 "MEP      = 450.D0",
                 "************************** 2 Higgs Doublet Model *************************",
                 "  TYPE: 1 (I), 2 (II), 3 (lepton-specific), 4 (flipped)",
                 "  PARAM: 1 (masses), 2 (lambda_i)",
                 "",
                 "PARAM    = 1",
                 "TYPE     = {}".format(num_2hdm),
                 "RENSCHEM = {}".format(renscheme),
                 "REFSCHEM = 5",
                 "********************",
                 "TGBET2HDM= {}D0".format(tanb),
                 "M_12^2   = {}D0".format(m12_2_dict[name]),
                 "INSCALE  = 125.0D0",
                 "OUTSCALE = MIN",
                 "******************** PARAM=1:",
                 "ALPHA_H  = {}D0".format(m.atan(tanb) - m.asin(sinbma)),
                 "MHL      = {}D0".format(mh),
                 "MHH      = {}D0".format(mH),
                 "MHA      = {}D0".format(mA),
                 "MH+-     = {}D0".format(mHc),
                 "******************** PARAM=2:",
                 "LAMBDA1  = 6.368674377530086700D0",
                 "LAMBDA2  = 0.235570240072350970D0",
                 "LAMBDA3  = 1.780416490847621700D0",
                 "LAMBDA4  = -1.52623758540479430D0",
                 "LAMBDA5  = 0.074592764717552856D0",
                 "**************************************************************************",
                 "SUSYSCALE= 2.22449e+03",
                 "MU       = -1.86701e+03",
                 "M2       = -2.39071e+02",
                 "MGLUINO  = 7.32754e+02",
                 "MSL1     = 1.49552e+03",
                 "MER1     = 1.62210e+03",
                 "MQL1     = 9.30379e+01",
                 "MUR1     = 2.77029e+03",
                 "MDR1     = 1.76481e+03",
                 "MSL      = 1.97714e+03",
                 "MER      = 9.29678e+02",
                 "MSQ      = 2.68124e+03",
                 "MUR      = 1.85939e+03",
                 "MDR      = 2.28235e+03",
                 "AL       = -4.62984e+03",
                 "AU       = 5.31164e+03",
                 "AD       = 2.54430e+03",
                 "ON-SHELL = 0",
                 "ON-SH-WZ = 0",
                 "IPOLE    = 0",
                 "OFF-SUSY = 0",
                 "INDIDEC  = 0",
                 "NF-GG    = 5",
                 "IGOLD    = 0",
                 "MPLANCK  = 2.40000e+18",
                 "MGOLD    = 1.00000e-13",
                 "******************* VARIATION OF HIGGS COUPLINGS *************************",
                 "ELWK     = 1",
                 "CW       = 1.D0",
                 "CZ       = 1.D0",
                 "Ctau     = 1.D0",
                 "Cmu      = 1.D0",
                 "Ct       = 1.D0",
                 "Cb       = 1.D0",
                 "Cc       = 1.D0",
                 "Cs       = 1.D0",
                 "Cgaga    = 0.D0",
                 "Cgg      = 0.D0",
                 "CZga     = 0.D0",
                 "********************* 4TH GENERATION *************************************",
                 "Ctp      = 0.D0",
                 "Cbp      = 0.D0",
                 "Cnup     = 0.D0",
                 "Cep      = 0.D0",
                 ]

    WriteListToFile(file_list,"CMSSW_10_2_19/src/2HDECAY/Input{}/2hdecay_{}.in".format(input_name,name))
   
  cmd_list = [
              "pushd CMSSW_10_2_19/src/2HDECAY/ > /dev/null",
              "source /vols/grid/cms/setup.sh",
              "eval `scram runtime -sh`",
              "python 2HDECAY.py --input=%(input_name)s &> /dev/null" % vars(),
              "popd > /dev/null"
               ]
  os.system(";".join(cmd_list))


  for (type_2hdm,mh,mH,mA,mHc,tanb,sinbma) in itertools.product(type_2hdms,mhs,mHs,mAs,mHcs,tanbs,sinbmas):
    name = MakeName(type_2hdm=type_2hdm,mh=mh,mH=mH,mA=mA,mHc=mHc,tanb=tanb,sinbma=sinbma)

    if name not in m12_2_dict: continue

    with open('CMSSW_10_2_19/src/2HDECAY/Results{}/2hdecay_{}_BR.out'.format(input_name,name)) as f:
      use_lines = False
      for line in f:
        if "QCD and EW" in line: use_lines = True
        if "QCD Only" in line: use_lines = False

        if use_lines and "DECAY QCD&EW" in line:
          width = float(line.split()[3].replace("E","e"))
          part = line.split("# ")[1].split()[0].replace("H+","Hc") 
          if name in width_dict.keys():
            width_dict[name][part] = width
          else:
            width_dict[name] = {part:width}

        if use_lines and "BR(" in line and ("h ->" in line or "H ->" in line or "A ->" in line or "H+ ->" in line):
          val = float(line.split()[0].replace("E","e"))
          proc_name = str("".join(line.split("(")[1].split(")")[0].rstrip().replace("->","_to_").replace("H+","Hc").replace("+","").replace("-","").split()))
          if name in br_dicts.keys():
            br_dicts[name][proc_name] = val
          else:
            br_dicts[name] = {proc_name:val}

  return width_dict, br_dicts

def calc_effective_couplings(al, be, yuktype):
    if yuktype == 1:
        dcs = calc_effective_cpls_type1(al, be)
    elif yuktype == 2:
        dcs = calc_effective_cpls_type2(al, be)
    elif yuktype == 3:
        dcs = calc_effective_cpls_type3(al, be)
    elif yuktype == 4:
        dcs = calc_effective_cpls_type4(al, be)
    else:
        raise RuntimeError
    return dcs

# calculate effective couplings in type 1
def calc_effective_cpls_type1(a, b):
    uu = cs(a) / ss(b)
    dd = uu
    ll = uu
    vv = ss(b - a)
    cplh = {
        'tautau': ll,
        'mumu': ll,
        'tt': uu,
        'bb': dd,
        'ZZ': vv,
        'WW': vv}
    uu = ss(a) / ss(b)
    dd = uu
    ll = uu
    vv = cs(b - a)
    cplH = {
        'tautau': ll,
        'mumu': ll,
        'tt': uu,
        'bb': dd,
        'ZZ': vv,
        'WW': vv}
    uu = 1 / np.tan(b)
    dd = 1 / np.tan(b)
    ll = 1 / np.tan(b)
    cplA = {
        'tautau': ll,
        'mumu': ll,
        'tt': uu,
        'bb': dd}
    return [cplh, cplH, cplA]

# calculate effective couplings in type 2
def calc_effective_cpls_type2(a, b):
    uu = cs(a) / ss(b)
    dd = -ss(a) / cs(b)
    ll = dd
    vv = ss(b - a)
    cplh = {
        'tautau': ll,
        'mumu': ll,
        'tt': uu,
        'bb': dd,
        'ZZ': vv,
        'WW': vv}
    uu = ss(a) / ss(b)
    dd = cs(a) / cs(b)
    ll = dd
    vv = cs(b - a)
    cplH = {
        'tautau': ll,
        'mumu': ll,
        'tt': uu,
        'bb': dd,
        'ZZ': vv,
        'WW': vv}
    uu = 1 / np.tan(b)
    dd = np.tan(b)
    ll = dd
    cplA = {
        'tautau': ll,
        'mumu': ll,
        'tt': uu,
        'bb': dd}
    return [cplh, cplH, cplA]

# calculate effective couplings in type 3
def calc_effective_cpls_type3(a, b):
    uu = cs(a) / ss(b)
    dd = uu
    ll = -ss(a) / cs(b)
    vv = ss(b - a)
    cplh = {
        'tautau': ll,
        'mumu': ll,
        'tt': uu,
        'bb': dd,
        'ZZ': vv,
        'WW': vv}
    uu = ss(a) / ss(b)
    dd = uu
    ll = cs(a) / cs(b)
    vv = cs(b - a)
    cplH = {
        'tautau': ll,
        'mumu': ll,
        'tt': uu,
        'bb': dd,
        'ZZ': vv,
        'WW': vv}
    uu = 1 / np.tan(b)
    dd = -1 / np.tan(b)
    ll = np.tan(b)
    cplA = {
        'tautau': ll,
        'mumu': ll,
        'tt': uu,
        'bb': dd}
    return [cplh, cplH, cplA]

# calculate effective couplings in type 4
def calc_effective_cpls_type4(a, b):
    uu = cs(a) / ss(b)
    dd = -ss(a) / cs(b)
    ll = uu
    vv = ss(b - a)
    cplh = {
        'tt': uu,
        'bb': dd,
        'tautau': ll,
        'mumu': ll,
        'ZZ': vv,
        'WW': vv}
    uu = ss(a) / ss(b)
    dd = cs(a) / cs(b)
    ll = uu
    vv = cs(b - a)
    cplH = {
        'tt': uu,
        'bb': dd,
        'tautau': ll,
        'mumu': ll,
        'ZZ': vv,
        'WW': vv}
    uu = 1 / np.tan(b)
    dd = np.tan(b)
    ll = uu
    cplA = {
        'tautau': ll,
        'mumu': ll,
        'tt': uu,
        'bb': dd}
    return [cplh, cplH, cplA]

def IterateAndCheckExperimentalConstraints(type_2hdms=["X"],mhs=[125.0],mHs=[200.0],mAs=[160.0],mHcs=[200.0],tanbs=[10.0,20.0,30.0,40.0,50.0],sinbmas=[1.0],widths_dict={},br_dicts={}):

  import Higgs.predictions as HP
  import Higgs.bounds as HB
  import Higgs.signals as HS

  pred = HP.Predictions()
  bounds = HB.Bounds('CMSSW_12_4_8/src/hbdataset') # load HB dataset
  signals = HS.Signals('CMSSW_12_4_8/src/hsdataset') # load HS dataset
  
  cs = np.cos
  ss = np.sin

  for (type_2hdm,mh,mH,mA,mHc,tanb,sinbma) in itertools.product(type_2hdms,mhs,mHs,mAs,mHcs,tanbs,sinbmas):
    name = MakeName(type_2hdm=type_2hdm,mh=mh,mH=mH,mA=mA,mHc=mHc,tanb=tanb,sinbma=sinbma)

    if name not in widths_dict: continue

    h = pred.addParticle(HP.BsmParticle("h", "neutral", "even"))
    h.setMass(mh)
    if widths_dict[name]["h"] > 0:
      h.setTotalWidth(widths_dict[name]["h"])
    else:
      print "Warning: h has negative width for {}. Setting to 1.".format(name)
      h.setTotalWidth(1.0)
 
    H = pred.addParticle(HP.BsmParticle("H", "neutral", "even"))
    H.setMass(mH)
    if widths_dict[name]["H"] > 0:
      H.setTotalWidth(widths_dict[name]["H"])
    else:
      print "Warning: H has negative width for {}. Setting to 1.".format(name)
      H.setTotalWidth(1.0)
    
    A = pred.addParticle(HP.BsmParticle("A", "neutral", "odd"))
    A.setMass(mA)
    if widths_dict[name]["A"] > 0:
      h.setTotalWidth(widths_dict[name]["A"])
    else:
      print "Warning:A has negative width for {}. Setting to 1.".format(name)
      A.setTotalWidth(1.0)

    X = pred.addParticle(HP.BsmParticle("X", "single"))
    X.setMass(mHc)
    if widths_dict[name]["Hc"] > 0:
      h.setTotalWidth(widths_dict[name]["X"])
    else:
      print "Warning: Hc has negative width for {}. Setting to 1.".format(name)
      X.setTotalWidth(1.0)

    for proc, br in br_dicts[name].items():

      if br > 0.0 and br < 1.0:
       
        final_state = i.split("_to_")[1] 
        final_state = final_state.replace("munu_mu","munu").replace("taunu_tau","taunu")
        if len(final_state) == 3 and final_state[-1] == "b": final_state = final_state[:-1]
 
        if "h_to" in proc:
          if not ("A" in final_state or "H" in final_state or "h" in final_state):
            h.setBr(final_state, val)
          else:
            h.setBr(final_state[0], final_state[1], val)

        if "H_to" in proc:
          if not ("A" in final_state or "H" in final_state or "h" in final_state):
            H.setBr(final_state, val)
          else:
            H.setBr(final_state[0], final_state[1], val)

        if "A_to" in proc:
          if not ("A" in final_state or "H" in final_state or "h" in final_state):
            A.setBr(final_state, val)
          else:
            A.setBr(final_state[0], final_state[1], val)

        if "Hc_to" in proc:
          if not ("A" in final_state or "H" in final_state or "h" in final_state):
            X.setBr(final_state, val)
          else:
            X.setBr(final_state[0], final_state[1], val)

    b = m.atan(tanb)
    a = b - a.sin(sinbma)
    cpl = calc_effective_couplings(a, b, thdm)

    cpls = HP.NeutralEffectiveCouplings()
    cpls.tt = cpl[0]['tt'] if 'tt' in cpl[0] else 0.0
    cpls.bb = cpl[0]['bb'] if 'bb' in cpl[0] else 0.0
    cpls.ZZ = cpl[0]['ZZ'] if 'ZZ' in cpl[0] else 0.0
    cpls.WW = cpl[0]['WW'] if 'WW' in cpl[0] else 0.0
    cpls.tautau = cpl[0]['tautau'] if 'tautau' in cpl[0] else 0.0
    cpls.mumu = cpl[0]['mumu'] if 'mumu' in cpl[0] else 0.0
    HP.effectiveCouplingInput(h,cpls,reference=HP.ReferenceModel.SMHiggsEW)

    cpls = HP.NeutralEffectiveCouplings()
    cpls.tt = cpl[1]['tt'] if 'tt' in cpl[1] else 0.0
    cpls.bb = cpl[1]['bb'] if 'bb' in cpl[1] else 0.0
    cpls.ZZ = cpl[1]['ZZ'] if 'ZZ' in cpl[1] else 0.0
    cpls.WW = cpl[1]['WW'] if 'WW' in cpl[1] else 0.0
    cpls.tautau = cpl[1]['tautau'] if 'tautau' in cpl[1] else 0.0
    cpls.mumu = cpl[1]['mumu'] if 'mumu' in cpl[1] else 0.0
    HP.effectiveCouplingInput(H,cpls,reference=HP.ReferenceModel.SMHiggsEW)

    cpls = HP.NeutralEffectiveCouplings()
    cpls.tt = cpl[2]['tt'] if 'tt' in cpl[2] else 0.0
    cpls.bb = cpl[2]['bb'] if 'bb' in cpl[2] else 0.0
    cpls.ZZ = cpl[2]['ZZ'] if 'ZZ' in cpl[2] else 0.0
    cpls.WW = cpl[2]['WW'] if 'WW' in cpl[2] else 0.0
    cpls.tautau = cpl[2]['tautau'] if 'tautau' in cpl[2] else 0.0
    cpls.mumu = cpl[2]['mumu'] if 'mumu' in cpl[2] else 0.0
    HP.effectiveCouplingInput(A,cpls,reference=HP.ReferenceModel.SMHiggsEW)

    res = bounds(pred)
    print res


valid, m12_2 = IterateAndCheckTheoreticalConstraintsAndFindm12_2()
width, brs = IterateAndFindWidthsAndBranchingRatios(m12_2_dict=m12_2)
print(width)
print(brs)
IterateAndCheckExperimentalConstraints(widths_dict=width,br_dicts=brs)

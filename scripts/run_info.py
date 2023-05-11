import functions
import numpy as np
import math as m
import os
import argparse
import glob
import json
import copy
import itertools
import ROOT
from array import array

parser = argparse.ArgumentParser()
parser.add_argument('--collect', help= 'Collect results and make root file of outputs',  action='store_true')
args = parser.parse_args()

num = 100
#mphis = [100,110,125,140,160,180,200,250,300]
#mAs = [60,70,80,90,100,125,140,160]
mphis = [200.0]
mAs = [160.0]
tanbs = [round(tanb,4) for tanb in list(np.logspace(np.log10(1.0),np.log10(100.0),num))]
csbmas = list(np.linspace(-1.0,1.0,num))

tanbs_uncentered = copy.deepcopy(tanbs)
csbmas_uncentered = copy.deepcopy(csbmas)

tanbs = functions.GetBinCenters(tanbs)
csbmas = functions.GetBinCenters(csbmas)

if not args.collect:
  ind = 0
  for mphi in mphis:
  
    if mphi > 125:
      lmH = [mphi]
      lmh = [125.0]
      sinbmas = [round(m.sin(m.acos(csbma)),4) for csbma in csbmas]
    else:
      lmH = [125.0]
      lmh = [mphi]
      sinbmas = [round(csbma,4) for csbma in csbmas]
  
    sinbmas = sorted(list(set(sinbmas)))
  
    for mA in mAs:
      print("mphi = {}, mA = {}".format(mphi,mA))
      for tanb in tanbs:
        ind += 1
  
        lmA = [mA]
        lmHc = [mphi]
        ltanb = [tanb]
  
        cmds = [
                "import sys",
                "sys.path.insert(0, '/vols/cms/gu18/4tau_v3/CMSSW_12_4_8/src/UserCode/2HDM-Information/scripts')",
                "import functions as f",
                "import json",
#                "m12_2, theory_excluded, exp_excluded, widths, brs = f.GetAllInformation(mhs={},mHs={},mAs={},mHcs={},tanbs={},sinbmas={},input_name={})".format(lmh,lmH,lmA,lmHc,ltanb,sinbmas,ind),
                "m12_2, theory_excluded, exp_excluded_hb, exp_excluded_hs, widths, brs = f.GetAllInformation(mhs={},mHs={},mAs={},mHcs={},tanbs={},sinbmas={})".format(lmh,lmH,lmA,lmHc,ltanb,sinbmas),
                "with open('output/thdm_info_theory_exclusions_{}.json', 'w') as outfile: json.dump(theory_excluded, outfile)".format(ind),
                "with open('output/thdm_info_experimental_exclusions_hb_{}.json', 'w') as outfile: json.dump(exp_excluded_hb, outfile)".format(ind),
                "with open('output/thdm_info_experimental_chisq_hs_{}.json', 'w') as outfile: json.dump(exp_excluded_hs, outfile)".format(ind),
                "with open('output/thdm_info_widths_{}.json', 'w') as outfile: json.dump(widths, outfile)".format(ind),
                "with open('output/thdm_info_brs_{}.json', 'w') as outfile: json.dump(brs, outfile)".format(ind),
                ]
  
        functions.CreateJob("jobs/thdm_info_{}.py".format(ind),cmds)
        os.system("python3 jobs/thdm_info_{}.py".format(ind))
        #cmssw_base = os.getcwd().replace('src/UserCode/2HDM-Information','')
        #functions.CreateBatchJob("jobs/thdm_info_{}.sh".format(ind),cmssw_base,["python3 jobs/thdm_info_{}.py".format(ind)])
        #functions.SubmitBatchJob("jobs/thdm_info_{}.sh".format(ind))

else:

  # collect json into one dictionary for each element
  br_fs = glob.glob("output/thdm_info_brs*.json")
  widths_fs = glob.glob("output/thdm_info_widths*.json")
  exp_exc_hb_fs = glob.glob("output/thdm_info_experimental_exclusions_hb*.json")
  exp_chisq_hs_fs = glob.glob("output/thdm_info_experimental_chisq_hs*.json")
  the_exc_fs = glob.glob("output/thdm_info_theory_exclusions*.json")

  br_d = {}
  widths_d = {}
  the_exc_d = {}
  exp_exc_hb_d = {}
  exp_chisq_hs_d = {}
  exp_exc_hs_d = {}

  for i in br_fs:
    f = open(i)
    data = json.load(f)
    for k, v in data.items():
      br_d[k] = v

  for i in widths_fs:
    f = open(i)
    data = json.load(f)
    for k, v in data.items():
      widths_d[k] = v

  for i in the_exc_fs:
    f = open(i)
    data = json.load(f)
    for k, v in data.items():
      the_exc_d[k] = v

  for i in exp_exc_hb_fs:
    f = open(i)
    data = json.load(f)
    for k, v in data.items():
      exp_exc_hb_d[k] = v

  for i in exp_chisq_hs_fs:
    f = open(i)
    data = json.load(f)
    for k, v in data.items():
      exp_chisq_hs_d[k] = v


  min_chisq = min(exp_chisq_hs_d.values())
  for k, v in exp_chisq_hs_d.items():
    exp_exc_hs_d[k] = ((v - min_chisq) > 5.99)

  # set up template histogram to fill
  bins1 = array('f', map(float,csbmas_uncentered))
  bins2 = array('f', map(float,tanbs_uncentered))
  hout = ROOT.TH2D('hout','',len(bins1)-1, bins1,len(bins2)-1, bins2)

  hists = {}
 
  for (mphi,mA,tanb,csbma) in itertools.product(mphis,mAs,tanbs,csbmas): 

    if mphi > 125:
      mH = mphi
      mh = 125.0
      sinbma = round(m.sin(m.acos(csbma)),4) 
    else:
      mH = 125.0
      mh = mphi
      sinbma = round(csbma,4) 

    entry_name = functions.MakeName(mh=mh,mH=mH,mA=mA,mHc=mphi,tanb=tanb,sinbma=sinbma)
    hist_name = "mh"+str(mh)+"_mH"+str(mH)+"_mA"+str(mA)

    x_bin = hout.GetXaxis().FindBin(csbma)
    y_bin = hout.GetYaxis().FindBin(tanb)

    # do exclusions

    exp_exc_hb_hist_name = "exp_exc_hb_"+hist_name
    exp_exc_hs_hist_name = "exp_exc_hs_"+hist_name
    the_exc_hist_name = "the_exc_"+hist_name

    if exp_exc_hb_hist_name not in hists.keys(): 
      hists[exp_exc_hb_hist_name] = copy.deepcopy(hout)
      hists[exp_exc_hb_hist_name].SetName(exp_exc_hb_hist_name.replace(".","p"))
    if exp_exc_hs_hist_name not in hists.keys():
      hists[exp_exc_hs_hist_name] = copy.deepcopy(hout)
      hists[exp_exc_hs_hist_name].SetName(exp_exc_hs_hist_name.replace(".","p"))
    if the_exc_hist_name not in hists.keys(): 
      hists[the_exc_hist_name] = copy.deepcopy(hout)
      hists[the_exc_hist_name].SetName(the_exc_hist_name.replace(".","p"))

    if entry_name in exp_exc_hb_d.keys(): hists[exp_exc_hb_hist_name].SetBinContent(x_bin,y_bin,exp_exc_hb_d[entry_name])
    if entry_name in exp_exc_hs_d.keys(): hists[exp_exc_hs_hist_name].SetBinContent(x_bin,y_bin,exp_exc_hs_d[entry_name])
    if entry_name in the_exc_d.keys(): hists[the_exc_hist_name].SetBinContent(x_bin,y_bin,the_exc_d[entry_name])

    # do branching ratios

    if entry_name in br_d.keys():
      for k,v in br_d[entry_name].items():            
        br_hist_name = "br_"+hist_name+"_"+k
        if br_hist_name not in hists.keys(): 
          hists[br_hist_name] = copy.deepcopy(hout) 
          hists[br_hist_name].SetName(br_hist_name.replace(".","p"))
        hists[br_hist_name].SetBinContent(x_bin,y_bin,v)

    # do widths

    if entry_name in widths_d.keys():
      for k,v in widths_d[entry_name].items():
        widths_hist_name = "widths_"+hist_name+"_"+k
        if widths_hist_name not in hists.keys(): 
          hists[widths_hist_name] = copy.deepcopy(hout)
          hists[widths_hist_name].SetName(widths_hist_name.replace(".","p"))
        hists[widths_hist_name].SetBinContent(x_bin,y_bin,v)


  hists["total_exp_exc_"+hist_name] = copy.deepcopy(hout)
  hists["total_exp_exc_"+hist_name].SetName("total_exp_exc_"+hist_name.replace(".","p"))
  hists["total_exc_"+hist_name] = copy.deepcopy(hout)
  hists["total_exc_"+hist_name].SetName("total_exc_"+hist_name.replace(".","p"))

  for k in hists.keys():
    if "total" in k: continue
    if "exp_exc_" in k:
      hists["total_exp_exc_"+hist_name].Add(hists[k])
    if "exc_" in k:
      hists["total_exc_"+hist_name].Add(hists[k])

  # make contours of limits
  conts = {}
  for k in hists.keys():
    if "_exc_" in k:
      cont_name = "contour_"+k.replace(".","p")
      clone = copy.deepcopy(hists[k].Clone())
      clone.Scale(-1)
      for xbin in range(1, clone.GetNbinsX() + 1):
        for ybin in range(1, clone.GetNbinsY() + 1):
          bin_content = clone.GetBinContent(xbin, ybin)
          clone.SetBinContent(xbin, ybin, bin_content + 1)
      conts[cont_name] = functions.contourFromTH2(clone,0.5,name=cont_name)

  # make contours of limits in the alignment scenario
  # WILL HAVE TO LOOP THROUGH ALL MASSES
  #bin_number = h.GetXaxis().FindBin(0)
  #h1 = h.ProjectionY("h1", bin_number, bin_number)
     
  out_file = "typeX_info.root"
  fout = ROOT.TFile(out_file, 'RECREATE')
  for k,v in hists.items(): 
    v.Write()
  for k,v in conts.items():
    print(k,v) 
    v.Write()
  fout.Close()


#for k, v in hists.items():
#  print(k)
#  v.Print("all")

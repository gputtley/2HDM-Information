import functions
import numpy as np
import math as m
import argparse
import json
import itertools
import ROOT
import glob
import time
import copy
import sys
import os
from array import array

parser = argparse.ArgumentParser()
parser.add_argument('--step', help= 'Step to run', type=str, default=None, choices=["2HDMC","2HDECAY","HiggsTools","Collect"])
parser.add_argument('--submit', help= 'Submit to the IC batch if possible',  action='store_true')
parser.add_argument('--n-per-job', help= 'Number of points per jobs', type=int, default=1000)
parser.add_argument('--output-name', help= 'Offset to begin from',  default="thdm_information")
parser.add_argument('--specific-ind', help= 'Specific ind to run', type=int, default=-1)
args = parser.parse_args()

# Input parameters
type_2hdm = "X"
num = 100
inputs = {
  "mphi" : [40.0,50.0,60.0,70.0,80.0,90.0,100.0,110.0,125.0,140.0,160.0,180.0,200.0,250.0,300.0,400.0,600.0,800.0],
  "mA" : [mA for mA in list(np.linspace(20, 1000, num))],
  "tanb" : [tanb for tanb in list(np.logspace(np.log10(1.0),np.log10(100.0),num))],
  "csbma" : [0.0],
}

# Histogramming parameters
x_hist_variable = "mA"
y_hist_variable = "tanb"

# center histogramming bins
x_hist_bins = copy.deepcopy(inputs[x_hist_variable])
y_hist_bins = copy.deepcopy(inputs[y_hist_variable])
x_bin_centers = [round((x_hist_bins[i] + x_hist_bins[i+1])/2,4) for i in range(len(x_hist_bins)-1)]
y_bin_centers = [round((y_hist_bins[i] + y_hist_bins[i+1])/2,4) for i in range(len(y_hist_bins)-1)]
inputs[x_hist_variable] = x_bin_centers
inputs[y_hist_variable] = y_bin_centers

# Useful parameters
job_ind = 0
type_2hdms = []
mhs = []
mHs = []
mHcs = []
mAs = []
tanbs = []
sinbmas = []
loop_through = list(itertools.product(inputs["mphi"],inputs["mA"],inputs["tanb"],inputs["csbma"]))

# Collect results
if args.step == "Collect":

  br_fs = glob.glob("output/br_dict_*.json")
  widths_fs = glob.glob("output/width_dict_*.json")
  exp_exc_hb_fs = glob.glob("output/hb_excluded_*.json")
  exp_chisq_hs_fs = glob.glob("output/hs_excluded_*.json")
  the_exc_fs = glob.glob("output/th_excluded_*.json")

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

  bins1 = array('f', map(float,x_hist_bins))
  bins2 = array('f', map(float,y_hist_bins))
  hout = ROOT.TH2D('hout','',len(bins1)-1, bins1,len(bins2)-1, bins2)

  hists = {}
  conts = {}


# Begin loop
for ind, (mphi,mA,tanb,csbma) in enumerate(loop_through): 

  # Set useful information
  if mphi > 125:
    mH = mphi
    mh = 125.0
    sinbma = round(m.sin(m.acos(csbma)),4) 
  else:
    mH = 125.0
    mh = mphi
    sinbma = round(csbma,4)

  entry_name = functions.MakeName(mh=mh,mH=mH,mA=mA,mHc=mphi,tanb=tanb,sinbma=sinbma)
  naming = {
    "mphi" : mphi,
    "mA" : mA,
    "tanb" : tanb,
    "csbma" : csbma,
  }

  # make job variables lists
  type_2hdms.append(type_2hdm)
  mhs.append(mh)
  mHs.append(mH)
  mHcs.append(mH)
  mAs.append(mA)
  tanbs.append(tanb)
  sinbmas.append(sinbma)

  if ((ind + 1) % args.n_per_job == 0) or (ind + 1 == len(loop_through)):
    
    print(f"- Job ind {job_ind}")
    start_time = time.time()

    if (args.specific_ind == -1 or args.specific_ind == job_ind):

      if args.submit and args.step != "2HDECAY":

        cmd = f"python3 {' '.join([i for i in sys.argv if '--submit' not in i and '--specific-ind' not in i])} --specific-ind={job_ind}"
        job_name = f"jobs/thdm_info_{args.step}_{job_ind}.sh"
        functions.CreateBatchJob(job_name, os.getcwd(), [cmd])
        functions.SubmitBatchJob(job_name)

      else:

        if args.step == "2HDMC":
          # Run 2HDMC
          print(" - Running 2HDMC")
          valid_dict, m12_2_dict = functions.CheckTheoreticalConstraintsAndFindm12_2(
            type_2hdms=type_2hdms,
            mhs=mhs,
            mHs=mHs,
            mAs=mAs,
            mHcs=mHcs,
            tanbs=tanbs,
            sinbmas=sinbmas
            )
      
          with open(f'output/valid_dict_{job_ind}.json', 'w') as json_file: json.dump(valid_dict, json_file)
          with open(f'output/m12_2_dict_{job_ind}.json', 'w') as json_file: json.dump(m12_2_dict, json_file)

        elif args.step == "2HDECAY":
          # Run 2HDECAY
          print(" - Running 2HDECAY")

          with open(f'output/valid_dict_{job_ind}.json', 'r') as json_file: valid_dict = json.load(json_file)
          with open(f'output/m12_2_dict_{job_ind}.json', 'r') as json_file: m12_2_dict = json.load(json_file)

          width_dict, br_dict, th_excluded = functions.FindWidthsAndBranchingRatios(
            type_2hdms=type_2hdms,
            mhs=mhs,
            mHs=mHs,
            mAs=mAs,
            mHcs=mHcs,
            tanbs=tanbs,
            sinbmas=sinbmas,
            m12_2_dict=m12_2_dict,
            excluded=valid_dict,
            renscheme=7,
            input_name=""
            )

          with open(f'output/width_dict_{job_ind}.json', 'w') as json_file: json.dump(width_dict, json_file)
          with open(f'output/br_dict_{job_ind}.json', 'w') as json_file: json.dump(br_dict, json_file)
          with open(f'output/th_excluded_{job_ind}.json', 'w') as json_file: json.dump(th_excluded, json_file)

        elif args.step == "HiggsTools":    
          # Run HiggsTools
          print(" - Running HiggsTools")

          with open(f'output/width_dict_{job_ind}.json', 'r') as json_file: width_dict = json.load(json_file)
          with open(f'output/br_dict_{job_ind}.json', 'r') as json_file: br_dict = json.load(json_file)

          hb_excluded, hs_excluded = functions.CheckExperimentalConstraints(
            type_2hdms=type_2hdms,
            mhs=mhs,
            mHs=mHs,
            mAs=mAs,
            mHcs=mHcs,
            tanbs=tanbs,
            sinbmas=sinbmas,
            widths_dict=width_dict,
            br_dicts=br_dict
            )

          with open(f'output/hb_excluded_{job_ind}.json', 'w') as json_file: json.dump(hb_excluded, json_file)
          with open(f'output/hs_excluded_{job_ind}.json', 'w') as json_file: json.dump(hs_excluded, json_file)

    # Update parameters
    job_ind += 1
    type_2hdm = []
    mhs = []
    mHs = []
    mHcs = []
    mAs = []
    tanbs = []
    sinbmas = []  

    end_time = time.time()
    print(f" - Time elapsed: {round(end_time-start_time,2)} seconds")

  if args.step == "Collect":

    # fill histograms
    entry_name = functions.MakeName(mh=mh,mH=mH,mA=mA,mHc=mphi,tanb=tanb,sinbma=sinbma)
    hist_name = "mh"+str(mh)+"_mH"+str(mH)+"_csbma"+str(csbma)
    hist_name = ""
    for k, v in naming.items():
      if k not in [x_hist_variable,y_hist_variable]:
        if hist_name != "": hist_name += "_"
        hist_name += k+str(v)

    x_bin = hout.GetXaxis().FindBin(naming[x_hist_variable])
    y_bin = hout.GetYaxis().FindBin(naming[y_hist_variable])

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

    if entry_name in exp_exc_hb_d.keys():
      if exp_exc_hb_d[entry_name] is not None:
        hists[exp_exc_hb_hist_name].SetBinContent(x_bin,y_bin,exp_exc_hb_d[entry_name])

    if entry_name in exp_exc_hs_d.keys():
      if exp_exc_hs_d[entry_name] is not None:
        hists[exp_exc_hs_hist_name].SetBinContent(x_bin,y_bin,exp_exc_hs_d[entry_name])

    if entry_name in the_exc_d.keys():
      if the_exc_d[entry_name] is not None:
        hists[the_exc_hist_name].SetBinContent(x_bin,y_bin,the_exc_d[entry_name])

    # do branching ratios
    if entry_name in br_d.keys():
      for k,v in br_d[entry_name].items():            
        br_hist_name = "br_"+hist_name+"_"+k
        if br_hist_name not in hists.keys(): 
          hists[br_hist_name] = copy.deepcopy(hout) 
          hists[br_hist_name].SetName(br_hist_name.replace(".","p"))
        if v is not None and not np.isnan(v):
          hists[br_hist_name].SetBinContent(x_bin,y_bin,v)

    # do widths
    if entry_name in widths_d.keys():
      for k,v in widths_d[entry_name].items():
        widths_hist_name = "widths_"+hist_name+"_"+k
        if widths_hist_name not in hists.keys(): 
          hists[widths_hist_name] = copy.deepcopy(hout)
          hists[widths_hist_name].SetName(widths_hist_name.replace(".","p"))
        if v is not None and not np.isnan(v):
          hists[widths_hist_name].SetBinContent(x_bin,y_bin,v)

    # make contours of limits
    for k in hists.keys():
      if "_exc_" in k and hist_name in k:
        cont_name = "contour_"+k.replace(".","p")
        clone = copy.deepcopy(hists[k].Clone())
        clone.Scale(-1)
        for xbin in range(1, clone.GetNbinsX() + 1):
          for ybin in range(1, clone.GetNbinsY() + 1):
            bin_content = clone.GetBinContent(xbin, ybin)
            clone.SetBinContent(xbin, ybin, bin_content + 1)
        conts[cont_name] = functions.contourFromTH2(clone,0.5,name=cont_name)


if args.step == "Collect":
  out_file = f"{args.output_name}.root"
  fout = ROOT.TFile(out_file, 'RECREATE')
  for k,v in hists.items(): 
    v.Write()
    #if "tautau" in k:
    #  print(k)
    #  v.Print("all")
  for k,v in conts.items():
    v.Write()
  fout.Close()


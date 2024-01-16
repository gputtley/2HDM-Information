# 2HDM-Information

This is a repository that combines 3 specific Higgs packages to provide information about the 2HDM:

 - 2HDMC - Places theoretical constraints on a theory.
 - 2HDECAY - Calculates branching fractions and widths of an allowed theory.
 - HiggsTools - Finds experimental exclusion limits from the inputs of the widths and branching fractions.

The 2HDMC code is set up to try and find a value of $m_{12}^{2}$ for given values of $m_{A}$, $m_{h}$, $m_{H}$, $tan\beta$ and $cos(\beta - \alpha)$, as well as the type of 2HDM. If it cannot find one it will return back information saying the theory is not theoretically allowed, if it can it will return the value of $m_{12}^{2}$. The 2HDECAY code then takes in this value of $m_{12}^{2}$ as well as the other parameters and calculates and returns the widths and branching fractions of each additional Higgs boson. This information is then passed to HiggsTools, which uses both the HiggsBounds and HiggsSignal databases, to return 95% CL exclusion points from the experimental measurements. This is all wrappered within the code to return many 2D root $m_{A}$-$tan\beta$ histograms and contours.

The code is run through the `run_step.py` code in `scripts`. This calls many functions in the `functions.py` file. It also contains functionality to parallelise on the IC batch. To run the code you must run the following commands:

```
python3 scripts/run_step.py --step=2HDMC --submit
```
```
python3 scripts/run_step.py --step=2HDECAY
```
```
python3 scripts/run_step.py --step=HiggsTools --submit
```
```
python3 scripts/run_step.py --step=Collect
```
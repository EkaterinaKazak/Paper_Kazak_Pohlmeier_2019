Readme file for the paper

Testing Out-of-Sample Portfolio Performance

Kazak and Pohlmeier, 2019

https://doi.org/10.1016/j.ijforecast.2018.09.010

The following folders contain MATLAB codes

1. size and power simulations:

- size_sim.m generates null rejection probabilites for size/power analysis
- 30.mat data used for empirical Monte Carlo
- deltalw_ce.m - funciton which computes standard errors for CE difference based on Delta Method from Ledoit and Wolf (2004)
- findx.m and rsimul.m are complementary functions used in size_sim.m

2. ROC
-ROC.m generates ROC curves Figure 3 in the paper.

3. pretest - empirical application of the pretest strategy from Section 4
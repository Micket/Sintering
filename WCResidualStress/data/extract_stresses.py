#!/usr/bin/python3

import numpy as np


# 0 "IST_StressTensor:0","IST_StressTensor:1","IST_StressTensor:2","IST_StressTensor:3","IST_StressTensor:4","IST_StressTensor:5","IST_StressTensor:6","IST_StressTensor:7","IST_StressTensor:8",
# 9 "IST_StrainTensor:0","IST_StrainTensor:1","IST_StrainTensor:2","IST_StrainTensor:3","IST_StrainTensor:4","IST_StrainTensor:5","IST_StrainTensor:6","IST_StrainTensor:7","IST_StrainTensor:8",
# 18 "IST_vonMisesStress",
# 19 "IST_CrossSectionNumber",
# 20 "IST_PlasticStrainTensor:0","IST_PlasticStrainTensor:1","IST_PlasticStrainTensor:2","IST_PlasticStrainTensor:3","IST_PlasticStrainTensor:4","IST_PlasticStrainTensor:5","IST_PlasticStrainTensor:6","IST_PlasticStrainTensor:7","IST_PlasticStrainTensor:8",
# 29 "IST_StressHyd",
# 30 "Cell Type"

files = [
        'wc_phase_elast_johansson',
        'co_phase_elast_johansson',
        'wc_phase_elast',
        'co_phase_elast',
        #'wc_phase_klow_plast',
        #'co_phase_klow_plast',
        'wc_phase_klow_plast2',
        'co_phase_klow_plast2',
        #'wc_phase_plast',
        #'co_phase_plast',
        'wc_phase_plast2',
        'co_phase_plast2',
        #'wc_phase_therm_plast',
        #'co_phase_therm_plast',
        'wc_phase_therm_plast2',
        'co_phase_therm_plast2'
        ]

for f in files:
    d = np.loadtxt('/mnt/dump/residual_stress/'+f+'.csv',  skiprows=1, delimiter=',')
    hyd = np.sum(d[:, [0,4,8]], 1)/3
    np.savetxt(f + '_hyd.txt', hyd[::] * 1e-9)
    np.savetxt(f + '_vm.txt', d[::,18] * 1e-9)
    print(f + " eq:", np.mean(d[:,18]) * 1e-6, "\pm", np.std(d[:,18]) * 1e-6)
    print(f + " Hyd:", np.mean(hyd) * 1e-6, "\pm", np.std(hyd) * 1e-6)


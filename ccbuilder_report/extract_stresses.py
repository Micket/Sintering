#!/usr/bin/python3

import numpy as np
d = np.loadtxt('../../CCBuilder/wc_phase.csv',       skiprows=1, delimiter=','); np.savetxt('wc_vol.txt',       d[::8,-2]); np.savetxt('wc_hyd.txt',       np.sum(d[::8,[0,5,8]], 1)/3 * 1e-9); np.savetxt('wc_vm.txt',       d[::8,18] * 1e-9);
print("WC Hyd:", np.mean(d[:,[0,5,8]]) * 1e-6, "\pm", np.std(d[:,[0,5,8]]) * 1e-6); print("WC eq:", np.mean(d[:,18]) * 1e-6, "\pm", np.std(d[:,18]) * 1e-6)

d = np.loadtxt('../../CCBuilder/wc_phase_elast.csv', skiprows=1, delimiter=','); np.savetxt('wc_vol_elast.txt', d[::8,-2]); np.savetxt('wc_hyd_elast.txt', np.sum(d[::8,[0,5,8]], 1)/3 * 1e-9); np.savetxt('wc_vm_elast.txt', d[::8,18] * 1e-9);
print("WC elastic Hyd:", np.mean(d[:,[0,5,8]]) * 1e-6, "\pm", np.std(d[:,[0,5,8]]) * 1e-6); print("WC elastic eq:", np.mean(d[:,18]) * 1e-6, "\pm", np.std(d[:,18]) * 1e-6)

d = np.loadtxt('../../CCBuilder/co_phase.csv',       skiprows=1, delimiter=','); np.savetxt('co_vol.txt',       d[::2,-2]); np.savetxt('co_hyd.txt',       np.sum(d[::2,[0,5,8]], 1)/3 * 1e-9); np.savetxt('co_vm.txt',       d[::2,18] * 1e-9);
print("Co Hyd:", np.mean(d[:,[0,5,8]]) * 1e-6, "\pm", np.std(d[:,[0,5,8]]) * 1e-6); print("Co eq:", np.mean(d[:,18]) * 1e-6, "\pm", np.std(d[:,18]) * 1e-6)

d = np.loadtxt('../../CCBuilder/co_phase_elast.csv', skiprows=1, delimiter=','); np.savetxt('co_vol_elast.txt', d[::2,-2]); np.savetxt('co_hyd_elast.txt', np.sum(d[::2,[0,5,8]], 1)/3 * 1e-9); np.savetxt('co_vm_elast.txt', d[::2,18] * 1e-9);
print("Co elastic Hyd:", np.mean(d[:,[0,5,8]]) * 1e-6, "\pm", np.std(d[:,[0,5,8]]) * 1e-6); print("Co elastic eq:", np.mean(d[:,18]) * 1e-6, "\pm", np.std(d[:,18]) * 1e-6)

d = np.loadtxt('../../CCBuilder/wc_only.csv', skiprows=1, delimiter=','); np.savetxt('wc_vol_only.txt', d[::2,-2]); np.savetxt('wc_hyd_only.txt', np.sum(d[::2,[0,5,8]], 1)/3 * 1e-9); np.savetxt('wc_vm_only.txt', d[::2,18] * 1e-9);
print("WC only Hyd:", np.mean(d[:,[0,5,8]]) * 1e-6, "\pm", np.std(d[:,[0,5,8]]) * 1e-6); print("WC only eq:", np.mean(d[:,18]) * 1e-6, "\pm", np.std(d[:,18]) * 1e-6)

#!/usr/bin/env python

import scipy as sp
import matplotlib.pyplot as plt
from cmath import phase, sqrt
from math import pi
from optparse import OptionParser

num_access_codes = 3
num_streams = 2

parser = OptionParser()
parser.add_option("-i", "--sync-index", type="int",
                  dest="start")
parser.add_option("-n", "--num-samples", type="int",
                  dest="samples")
parser.add_option("-N", "--num-tx-samples", type="int",
                  dest="txsamples")
parser.add_option("-p", "--plateau-index", type="int",
                  dest="plateau_index")

(options, args) = parser.parse_args()

tx1 = sp.fromfile(open("/tmp/tx1.dat"), dtype=sp.complex64)
tx2 = sp.fromfile(open("/tmp/tx2.dat"), dtype=sp.complex64)
rx1 = sp.fromfile(open("/tmp/rx1.dat"), dtype=sp.complex64)
rx2 = sp.fromfile(open("/tmp/rx2.dat"), dtype=sp.complex64)
f_sc1 = sp.fromfile(open("/tmp/f_sc_1.dat"), dtype=sp.float32)
f_sc2 = sp.fromfile(open("/tmp/f_sc_2.dat"), dtype=sp.float32)
tx_sig1  = sp.fromfile(open("/tmp/tx_sig1.dat"), dtype=sp.complex64)
tx_sig2  = sp.fromfile(open("/tmp/tx_sig2.dat"), dtype=sp.complex64)
rx_sig1  = sp.fromfile(open("/tmp/rx_sig1.dat"), dtype=sp.complex64)
rx_sig2  = sp.fromfile(open("/tmp/rx_sig2.dat"), dtype=sp.complex64)
f_s0 = []
f_s1 = []
for chan in range(num_streams):
  f_s0.append(sp.fromfile(open("/tmp/corr_" + str(chan + 1) + "_0.dat"), dtype=sp.float32))
  f_s1.append([])
  for ac_id in range(num_access_codes*num_streams):
    f_s1[chan].append(sp.fromfile(open("/tmp/corr_" + 
                                       str(chan + 1) + "_" +
                                       str(ac_id + 1) +
                                       ".dat"), dtype=sp.float32))

if(options.start != None):
  start = options.start
elif(options.plateau_index != None):
  start = options.plateau_index - 200
else:
  start = 0

if(options.samples != None):
  end = start + options.samples
else:
  end = start + len(tx1)

if(options.txsamples != None):
  txend = 0 + options.txsamples
else:
  txend = len(tx1)

print len(f_sc1), len(f_sc2)

tx, tx_axarr = plt.subplots(2, sharex=True)
title = "TX Signal"
tx_axarr[0].set_title(title)
tx_axarr[0].plot([x.real for x in tx1[:txend]], 'r', label="Ch1-Real")
tx_axarr[0].plot([x.imag for x in tx1[:txend]], 'b', label="Ch1-Imag")
tx_axarr[1].plot([x.real for x in tx2[:txend]], 'g', label="Ch2-Real")
tx_axarr[1].plot([x.imag for x in tx2[:txend]], 'y', label="Ch2-Imag")
tx_axarr[0].legend(loc=4)
tx_axarr[1].legend(loc=4)

rx, rx_axarr = plt.subplots(2, sharex=True)
title = "RX Signal"
rx_axarr[0].set_title(title)
rx_axarr[0].plot([x.real for x in rx1[start:end]], 'r', label="Ch1-Real")
rx_axarr[0].plot([x.imag for x in rx1[start:end]], 'b', label="Ch1-Imag")
rx_axarr[1].plot([x.real for x in rx2[start:end]], 'g', label="Ch2-Real")
rx_axarr[1].plot([x.imag for x in rx2[start:end]], 'y', label="Ch2-Imag")
rx_axarr[0].legend(loc=4)
rx_axarr[1].legend(loc=4)

corr, corr_axarr = plt.subplots(num_streams, sharex=True)
title = "Access Code correlations"
corr_axarr[0].set_title(title)
for chan in range(num_streams):
  corr_axarr[chan].plot(f_s0[chan], label= "S0 - Ch" + str(chan + 1))
  for ac_id in range(num_access_codes*num_streams):
    corr_axarr[chan].plot(f_s1[chan][ac_id],
                          label = "Ch" + str(chan + 1) + " - AC" + str(ac_id + 1))
#  corr_axarr[chan].legend(bbox_to_anchor=(1.1, 1.1))

title = "Sync Output"
plt.figure(4)
plt.plot(f_sc1[len(f_sc1) - 64*4:], 'b', label="Ch1")
plt.plot(f_sc2[len(f_sc2) - 64*4:], 'r', label="Ch2")

sig, sig_ax = plt.subplots(4, sharex=True)
sig_ax[0].plot([x.real for x in tx_sig1], label="Real")
sig_ax[0].plot([x.imag for x in tx_sig1], label="Imag")
sig_ax[1].plot([x.real for x in rx_sig1], label="Real")
sig_ax[1].plot([x.imag for x in rx_sig1], label="Imag")
sig_ax[2].plot([x.real for x in tx_sig2], label="Real")
sig_ax[2].plot([x.imag for x in tx_sig2], label="Imag")
sig_ax[3].plot([x.real for x in rx_sig2], label="Real")
sig_ax[3].plot([x.imag for x in rx_sig2], label="Imag")
sig_ax[0].legend(loc=0)
sig_ax[1].legend(loc=0)
plt.show()

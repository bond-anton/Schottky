from __future__ import division
from os.path import join, dirname
import numpy as np
import pandas as pd
import mpmath as mp
from matplotlib import pyplot as plt

V_p = 3
V_rb = -2

T_start = 40.0
T_stop = 80.0
T_step = 1.0

V_p = 3.0
V_rb = -2.0

t1 = np.array([1e-6, 5e-6, 2e-6, 1e-5, 2e-5, 5e-5, 1e-5])
a = 10
sampling_times = [np.array([t1_i, a * t1_i]) for t1_i in t1]
T_range = np.linspace(T_start, T_stop, abs((T_start - T_stop) / T_step) + 1, endpoint=True)

data_dir = join(dirname(__file__), '03-kinetics')
DLTS = []
for T in T_range:
    csv_file_name = '03_Au_nSi_BW_transient' + '_%02.2fVp_%02.2fVrb_%03.2fK' % (V_p, abs(V_rb), T) + '.csv'
    csv_file_name = join(data_dir, csv_file_name)

    print '\nReading in', csv_file_name
    df = pd.read_csv(csv_file_name)
    df.set_index('time', inplace=True)
    DLTS_T = []
    for sampling_time in sampling_times:
        df.drop_duplicates(take_last=True, inplace=True)
        try:
            df_resampled = df.reindex(np.union1d(df.index.values, sampling_time))
            df_resampled = df_resampled.apply(pd.Series.interpolate)
        except ValueError:
            print 'dont work'
            df_resampled = df
            pass
        #print df.loc[sampling_time].values[:,1]
        samples = df_resampled.loc[sampling_time].values[:, 1]
        result = samples[0] - samples[1]
        print 'T =', T, 'K, S =', result, samples
        DLTS_T.append(result)
    DLTS.append(DLTS_T)
DLTS = np.array(DLTS)
for i in range(len(t1)):
    plt.plot(T_range, DLTS[:, i], '-o', label='t1=%2.2gs' % t1[i])
#plt.legend()
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()


import os
import re

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pyteomics import mzxml, auxiliary

#D:\TapDev/Lund
#718417303 8137ame

mzxml_file = ['150210_1_01.mzXML','150210_2_01.mzXML','150210_3_01.mzXML','150210_4_01.mzXML','150210_5_01.mzXML','150210_6_01.mzXML','150210_7_01.mzXML','150210_8_01.mzXML','150210_9_01.mzXML','150210_11_01.mzXML','150210_13_01.mzXML','150210_14_01.mzXML','150210_15_01.mzXML','150210_16_01.mzXML','150210_17_01.mzXML','150210_18_01.mzXML','150210_19_01.mzXML','150210_20_01.mzXML']
a
spectra_adress = 'tablecreator/event_spectral_scan.txt'
spectra_index_data = load_data(spectra_adress)
def get_ms_spectra(peak_id):
    cleaned_index = spectra_index_data[peak_id].map(lambda id_string: id_string.replace('-','').split('}')[0])
    re_index = cleaned_index[cleaned_index != '']
    return [i for i in zip(re_index.index,re_index)]

id_spectra = get_ms_spectra(111)

file_name= '150210_1_01.mzXML'

def mz_main(scan_id, file_name):
    scan_id = str(scan_id)
    mzxml_file = mzxml.MzXML('mzXML/'+file_name)
    #star_bin = mzxml_file['15540']['offset']
    #end_bin = mzxml_file['15541']['offset']
    #name = mzxml_file.read()[star_bin:end_bin]
    raw_array_spectra = ''
    for item in mzxml_file.iterfind("scan[num='"+scan_id+"']"):
        raw_array_spectra = [item['m/z array'], item['intensity array']]
        break
    return raw_array_spectra 

def get_ms_spectra(db_aray):
    df_array = pd.DataFrame(db_aray).T
    df_array[0] = df_array[0].map(lambda x : round(x, 2))
    cleaned_df = df_array.groupby(0).mean()
    return cleaned_df

def get_peaks(raw_ms):
    non_zero = raw_ms[raw_ms[1]!=0]
    up_num = pd.DataFrame([0]+list(non_zero[1].iloc[:-1]), index = non_zero.index, columns = [1])
    down_num = pd.DataFrame(list(non_zero[1].iloc[1:]) + [0], index = non_zero.index, columns = [1])
    peaks = non_zero[((non_zero - up_num) > 0 ) & ((non_zero - down_num) > 0)].dropna()
    #return peaks
    sigma_j = np.sqrt(sum(peaks.apply(lambda x:x**2)))
    peaks = peaks.sort_values(1).tail(int(len(peaks)*0.1))
    return peaks/sigma_j

#ms_spectra_array = [mz_main(idsp, file_name) for idsp in range(8400,8405,1)]
ms_spectra_array = [mz_main(idsp[1], mzxml_file[idsp[0]]) for idsp in id_spectra]
ms_spectra_df = [get_ms_spectra(i) for i in ms_spectra_array]

#average intensity of spectra
avg_spectra = ms_spectra_df[0]
for each_spectrum in ms_spectra_df[1:]:
    avg_spectra = avg_spectra.add(each_spectrum, fill_value=0)
avg_spectra = avg_spectra/len(ms_spectra_df)
#plot
avg_spectra.plot()
plt.show()

'''
for i in range(5):
     ms_spectra.append(pd.read_csv(str(i)+' .csv',index_col=0))
 

plt.figure(1)
start_sub = closestNonPrime(len(ms_spectra_df)+1)
axs = [plt.subplot(start_sub[0],start_sub[1],i+1) for i in range(len(ms_spectra_df))]
for i in range(len(ms_spectra_df)):
    plt.sca(axs[i])
    axs[i] = ms_spectra_df[i].plot()
ms_spectra_df[3].plot()
plt.show()
'''
plt.figure(1)
axs = plt.subplot(111)
max_int = (ms_spectra_df[1].max(),ms_spectra_df[2].max())
Axes.set_ylim(bottom = -max_int, top=max_int)
axs.plot(ms_spectra_df[1].index, ms_spectra_df[1][1],'b',ms_spectra_df[2].index, -ms_spectra_df[2][1],'r')

ms_spectra_peaks = [get_peaks(i) for i in ms_spectra_df]


def matching(df1, df2, tol = 1000):
    if len(df1) > len(df2):
        trans = df1.copy()
        df1 = df2.copy()
        df2 = trans.copy()
    match_mz = []
    match_int = []
    match_tol = []
    for peak_m_z in df1.index:
        abs_tol = peak_m_z * tol
        can = df2.copy()
        can[0] = abs(df2.index - peak_m_z)*1e6/peak_m_z
        can_va = can.sort_values(0).iloc[0]
        match_mz.append(can_va.name)
        match_int.append(can_va[1])
        match_tol.append(can_va[0])
    res_data = df1.copy()
    res_data.columns = ['int']
    res_data['mz'] = df1.index
    res_data['match_mz'] = match_mz
    res_data['match_int'] = match_int
    res_data['match_tol'] = match_tol
    res_data.index = range(len(res_data))
    #return res_data
    res_data1 =  res_data[res_data.match_tol < tol]
    norm1 = np.sqrt(sum(res_data1['mz'].map(lambda x:x**2)))
    norm2 = np.sqrt(sum(res_data1['match_mz'].map(lambda x:x**2)))
    cos_score = res_data1['mz'].dot(res_data1['match_mz'])/(norm1*norm2)

    norm3 = np.sqrt(sum(res_data1.iloc[:,0].map(lambda x:x**2)))
    norm4 = np.sqrt(sum(res_data1.iloc[:,3].map(lambda x:x**2)))
    weight = res_data1.iloc[:,0].dot(res_data1.iloc[:,3])/(norm3*norm4)

    #aandb = len(res_data1)
    return res_data1,cos_score*weight#*np.sqrt(aandb/(len(df1)+len(df2)-aandb))


hm = [[matching(i,j)[1] for i in ms_spectra_peaks]for j in ms_spectra_peaks]
sns.heatmap(hm, cmap="YlGnBu")
sns.heatmap(hm, cmap="YlGnBu", vmin=0, vmax=1)
#plt.imshow(hm, interpolation='nearest')




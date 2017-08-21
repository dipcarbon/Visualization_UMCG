import re
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

class ExtractedIonChrome():
    def __init__(self, eic_path, text_path):
        self.eic_path = eic_path
        self.file_path = eic_path+'/FilesToMatch.txt'
        with open(self.file_path) as match_file:
            self.matched_index_dict = [line.strip().split() for line in match_file]
        self.refer_file_name = self.matched_index_dict[0][0].split('/')[-1].split('.')[0]
        #self.refer_file_name = '150210_9_01'
        self.contrl_group_list = [i[0].split('/')[-1].split('.')[0] for i in self.matched_index_dict if i[1] == '0']
        #self.contrl_group_list = ['150210_1_01.mesh','150210_2_01.mesh','150210_3_01.mesh','150210_4_01.mesh','150210_5_01.mesh','150210_6_01.mesh','150210_7_01.mesh','150210_8_01.mesh','150210_9_01.mesh']
        self.main_file_list = os.listdir(eic_path)
        self.tmap_list = tuple([i_list for i_list in self.main_file_list if i_list[-5:] == '.tmap'])
        self.eic_dict = {}

        with open(text_path) as file:
            text = file.read().split('\n\n')[:-1]
        for each_text in text:
            panel = each_text.split('\n')
            mass_index = panel.index('Mass')
            mass_list = [float(i) for i in panel[mass_index+1].split()]
            mass_len_list = len(mass_list)
            if mass_len_list%2:
                panel_key = round(mass_list[mass_len_list//2], 3)
            else:
                panel_key = round(sum(mass_list[(mass_len_list//2-1):mass_len_list//2])/2, 3)
            panel_name = panel[-1].split("\\")[-1]
            panel_rt = [float(i) for i in panel[mass_index+3].split()]
            panel_int = np.array([.0]*len(panel_rt))
            for panel_ints in panel[1:mass_index]:
                panel_int += np.array([float(i) for i in panel_ints.split()])/(mass_index-1)
            panel_package = {'int':panel_int, 'rt':panel_rt}
            if panel_key in self.eic_dict.keys():
                self.eic_dict[panel_key][panel_name] = panel_package
            else:
                self.eic_dict[panel_key] = {panel_name:panel_package}

        self.mass_list = tuple(self.eic_dict.keys())
        self.file_id_list = tuple([i.split('.')[0] for i in self.eic_dict[self.mass_list[0]].keys()])
        self.align_r_time_dict ={}
        for each_file in self.file_id_list:
            self.align_r_time_dict[each_file] = {'o_list':[], 't_list':[]}
            if each_file in self.refer_file_name:
                for tmap_each_list in self.tmap_list:
                    if len(tmap_each_list.split(each_file)) == 3:
                        each_name_file = tmap_each_list
            else:
                each_name_file = [i for i in self.tmap_list if each_file in i][0]

            with open(eic_path+'/'+each_name_file) as tmap_file:
                for line in tmap_file:
                    words = line.split('\t')
                    self.align_r_time_dict[each_file]['o_list'].append(words[0])
                    self.align_r_time_dict[each_file]['t_list'].append(words[1].strip())


    def rt_int_plot(self, mass_index_value, after_alignment=1, before_alignment=1):
        eic_box_dict = self.eic_dict[mass_index_value]
        red_patch = mpatches.Patch(color='red', label='Sample Group')
        blue_patch = mpatches.Patch(color='blue', label='Contrl Group')
        self.eic_dict[mass_index_value]
        if before_alignment and after_alignment:
            ax = plt.subplot(211)
            ax1 = plt.subplot(212)
        elif before_alignment:
            ax = plt.subplot(111)

        else:
            ax1 = plt.subplot(111)

        for text_item in eic_box_dict.items():
            ms_int = text_item[1]['int']
            ms_rt = text_item[1]['rt']
            align = self.align_r_time_dict[text_item[0].split('.')[0]]
            ms_align_rt = np.interp(ms_rt,align['o_list'],align['t_list'])

            if text_item[0] in self.contrl_group_list:
                color_code = 'b'
            else:
                color_code = 'r'
            if before_alignment:
                ax.plot(ms_rt, ms_int, color_code, linewidth=0.6,)
            if after_alignment:
                ax1.plot(ms_align_rt, ms_int, color_code, linewidth=0.6,)
        
        if before_alignment:
            ax.legend(handles=[red_patch, blue_patch])
            ax.set_xlabel('Orignal Retention Time')
            ax.set_ylabel('Intensity')
        if after_alignment:
            ax1.legend(handles=[red_patch, blue_patch])
            ax1.set_xlabel('Aligned Retention Time')
            ax1.set_ylabel('Intensity')
        plt.show()

'''
eic_path = "/Users/dippercheng/Desktop/groningen/MatlabEIVVisualisation"
text_path = "inputTIC_10MelindaSerumAscii.txt.txt"
spectra = ExtractedIonChrome(eic_path, text_path)
#mass list
mlist = spctra.mass_list
#ploting
spectra.rt_int_plot(mlist[0], after_alignment=1, before_alignment=1)
spectra.rt_int_plot(mlist[0], after_alignment=1, before_alignment=0)
spectra.rt_int_plot(mlist[0], after_alignment=0, before_alignment=1)
'''


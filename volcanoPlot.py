import os
import re
from collections import Counter

from Bio import SeqIO
from scipy import stats
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

protein_number_threshold = 2
marker_style = 'o'
my_point_size = 14
alpha_number = 0.6
contrl_group_id = list(range(9))
sample_group_id = list(range(9, 18))
colors_id = ('#5EA3ED','#FF9C00','#0AD670','#E54EC7','#6A4A3C','#00BFFF','#FFA07A','#FE4365')

fasta_adress = 'fasta/Uniprot_22072015_concatenated_target_decoy.fasta'
#fasta_adress = 'fasta/06-12043-database_concatenated_target_decoy.fasta'
ms_adress = 'metamatchout/merged_output_18_1D_50cm_mzRadius=0.010_TRadius=0.75_Fraction=0.75_mzRange=350-1800_rtRange=10-120_overlap=5.mpks'
protein_adress = 'tablecreator/fasta_sequence_id.txt'
peptide_adress = 'tablecreator/peptide_id.txt'

fasta_id = {}

def find_fasta_id(fasta_adress):
    global fasta_id
    FT_value = {'sp':True, 'tr':False}
    raw_fasta = SeqIO.parse(fasta_adress, "fasta")
    for record in raw_fasta:
        bargage = record.description.split('|')
        if bargage[0] in 'trsp':
            fasta_id[bargage[1]] = FT_value[bargage[0]]
    #return fasta_id


class AnnoteFinder(object):

    def __init__(self, xdata, ydata, annotes, ax=None, xtol=None, ytol=None):
        self.data = list(zip(xdata, ydata, annotes))
        if xtol is None:
            xtol = ((max(xdata) - min(xdata))/float(len(xdata)))/2
        if ytol is None:
            ytol = ((max(ydata) - min(ydata))/float(len(ydata)))/2
        self.xtol = xtol
        self.ytol = ytol
        if ax is None:
            self.ax = plt.gca()
        else:
            self.ax = ax
        self.drawnAnnotations = {}
        self.links = []
        self.drawnDetails = {}
        self.linksDetail = []
        self.tb_mode = ''


    def distance(self, x1, x2, y1, y2):
        return(np.sqrt((x1 - x2)**2 + (y1 - y2)**2))

    def __call__(self, event):
        tb = plt.get_current_fig_manager().toolbar
        print(event)
        if event.button == 1 and (tb.mode == self.tb_mode):
            print(1)
            clickX = event.xdata
            clickY = event.ydata
            if (self.ax is None) or (self.ax is event.inaxes):
                annotes = []
                #print(event.xdata, event.ydata)
                for x, y, a in self.data:
                    #print(x, y, a)
                    if ((clickX-self.xtol < x < clickX+self.xtol) and
                            (clickY-self.ytol < y < clickY+self.ytol)):
                        annotes.append(
                            (self.distance(x, clickX, y, clickY), x, y, a))
                if annotes:
                    annotes.sort()
                    distance, x, y, annote = annotes[0]
                    self.drawAnnote(event.inaxes, x, y, annote)
                    for l in self.links:
                        l.drawSpecificAnnote(annote)
        elif event.button == 1:
            self.tb_mode = tb.mode

        if event.button == 3:
            clickX = event.xdata
            clickY = event.ydata
            if (self.ax is None) or (self.ax is event.inaxes):
                annotes = []
                #print(event.xdata, event.ydata)
                for x, y, a in self.data:
                    #print(x, y, a)
                    if ((clickX-self.xtol < x < clickX+self.xtol) and
                            (clickY-self.ytol < y < clickY+self.ytol)):
                        annotes.append(
                            (self.distance(x, clickX, y, clickY), x, y, a))
                if annotes:
                    annotes.sort()
                    distance, x, y, annote = annotes[0]
                    self.drawDetail(event.inaxes, x, y, annote)
                    for l in self.linksDetail:
                        l.drawSpecificDetail(annote)

    def drawAnnote(self, ax, x, y, annote):
        if (x, y) in self.drawnAnnotations:
            markers = self.drawnAnnotations[(x, y)]
            for m in markers:
                m.set_visible(not m.get_visible())
            self.ax.figure.canvas.draw_idle()
        else:
            t = ax.text(x, y, "")
            m = ax.scatter([x], [y], s=my_point_size-12, marker=marker_style, c='r', zorder=100)
            self.drawnAnnotations[(x, y)] = (t, m)
            self.ax.figure.canvas.draw_idle()

    def drawSpecificAnnote(self, annote):
        annotesToDraw = [(x, y, a) for x, y, a in self.data if a == annote]
        for x, y, a in annotesToDraw:
            self.drawAnnote(self.ax, x, y, a)

    def drawDetail(self, ax, x, y, annote):
        if (x, y) in self.drawnDetails:
            markers = self.drawnDetails[(x, y)]
            for m in markers:
                m.set_visible(not m.get_visible())
            self.ax.figure.canvas.draw_idle()
        else:
            bbox_dict = dict(boxstyle='round', facecolor='#FF5554', alpha=0.95, edgecolor='white')
            t = ax.text(x+0.1, y+0.1, "  %s" % (annote),fontsize=8, bbox=bbox_dict, zorder=101)
            m = ax.scatter([x], [y], s=my_point_size+1, marker='', c='r', zorder=100)
            self.drawnDetails[(x, y)] = (t, m)
            self.ax.figure.canvas.draw_idle()

    def drawSpecificDetail(self, annote):
        annotesToDraw = [(x, y, a) for x, y, a in self.data if a == annote]
        for x, y, a in annotesToDraw:
            self.drawDetail(self.ax, x, y, a)

############################################################################


def load_data(file_id):
    if file_id == 'pep':
        load_adress = peptide_adress
    elif file_id == 'pro':
        load_adress = protein_adress
    else:
        load_adress = file_id
    line_file = []
    with open(load_adress) as f:
        for line in f:
            line_file.append(line)
        line_file_1 = line_file[1:]
    line_file_2 = [line_2.strip().split()[1:] for line_2 in line_file_1]
    ref_data = pd.DataFrame(line_file_2,).T
    return ref_data


def protein_string_clean(strings):
    strings_1 = re.sub('-','',strings)
    if not strings_1:
        return ''
    if '}' in strings_1:
        strings_1 = strings_1.split('}')[0]
    strings_group = strings_1.split(';')
    unq_strings_group = []
    for s in strings_group:
        if s and (s in fasta_id.keys()):
            if fasta_id[s]:
                unq_strings_group.append(s)
    unq_strings_group = list(set(unq_strings_group))
    if not unq_strings_group:
        return ''
    return ';'.join(unq_strings_group)


def peptide_string_clean(strings):
    strings_1 = re.sub('-','',strings)
    if not strings_1:
        return ''
    if '}' in strings_1:
        strings_1 = strings_1.split('}')[0]
    strings_group = strings_1.split(';')
    unq_strings_group = [s for s in strings_group if s]
    unq_strings_group = list(set(unq_strings_group))
    if not unq_strings_group:
        return ''
    alt_unq_strings_group = []
    for u in unq_strings_group:
        if '_' in u:
            alt_unq_strings_group.append(u.split('_')[0])
        else:
            alt_unq_strings_group.append(u)
    return ';'.join(alt_unq_strings_group)


def get_annotes(data_table, win_id=1):
    if win_id != 3:
        data_table = data_table[['pv', 'fc', 'pep', 'pro']]
        data_table.columns = ['P.V. ','F.C. ', 'PEP. ','PRO. ']
    else:
        data_table = data_table[['pv', 'fc', 'pro']]
        data_table.columns = ['P.V. ',' F.C. ',' PRO. ']
    annotes = [data_table.iloc[i].to_csv(sep=':')[:-1] for i in range(len(data_table))]
    return annotes


def data_parsing():
    find_fasta_id(fasta_adress)
    
    if not fasta_id:
        print('Loading Error')
        return None

    pep_ttap_data = load_data('pep').iloc[contrl_group_id].applymap(peptide_string_clean)
    pro_ttap_data = load_data('pro').iloc[contrl_group_id].applymap(protein_string_clean)
    with open(ms_adress,'r') as raw_meta_data:
        for raw_label_columns in raw_meta_data:
            break
    label_columns_pop = raw_label_columns.find('File')
    label_columns = raw_label_columns[label_columns_pop:].strip().split()
    meta_data = pd.read_csv(ms_adress, sep=' ', )[label_columns].T
    #caculate  p_value and fold_change
    p_value = pd.Series(stats.ttest_ind(meta_data.iloc[contrl_group_id], meta_data.iloc[sample_group_id])[1])
    fold_change = meta_data.iloc[sample_group_id].sum() / meta_data.iloc[contrl_group_id].sum()
    #log
    log10_p_value = p_value.map(lambda x: -np.log10(x))
    log2_fold_change = fold_change.map(lambda x: np.log2(x))
    #form, column0 means log10_p_value, column1 means log2_fold_change
    ms_table = pd.DataFrame([log10_p_value, log2_fold_change]).T
    ms_table.columns=['pv','fc']
    #find unique peptide
    uni_peptide_list = []
    for uni_proteins_id in pep_ttap_data:
        string_protein = ';'.join(pep_ttap_data[uni_proteins_id])
        protein_bag = [sp for sp in string_protein.split(';') if sp]
        cnt = Counter(protein_bag).most_common()
        if not cnt:
            uni_peptide_list.append('')
        elif cnt[0][1] < protein_number_threshold:
            uni_peptide_list.append('')
        else:
            uni_peptide_list.append(cnt[0][0])
    #get peak ids for unique peptide
    id_peptide_dict = {}
    uni_peptide_df = pd.DataFrame(uni_peptide_list)
    ms_table['pep'] = uni_peptide_df
    for each_uni_peptide_df in uni_peptide_df.groupby(0):
        if each_uni_peptide_df[0]:
            each_uni_string = list(each_uni_peptide_df[1].index)
            id_peptide_dict[each_uni_peptide_df[0]] = each_uni_string

    peptide_table = pd.DataFrame(list(id_peptide_dict.keys()),columns=['pep'])
    #
    def find_seq_lite(seque):
        peak_id = id_peptide_dict[seque][0]
        seq_ref_data = pro_ttap_data[peak_id]
        pep_ref_data  = pep_ttap_data[peak_id]
        pro_index = pep_ref_data [pep_ref_data == seque].index[0]
        candidate = seq_ref_data[pro_index]
        if ';' in candidate:
            return ''
        return candidate
    #
    peptide_table['pro'] = peptide_table['pep'].map(find_seq_lite)
    #add protein info to ms_plot
    ms_table['pro'] = ms_table['pep'].map(lambda x: ';'.join(peptide_table[peptide_table.pep==x]['pro']))
    ##get peak ids for each protein fasta id
    id_protein_dict = {}
    ##get peptide id for each protein fasta id
    id_protein_dict_pep = {}
    for each_pro_df in peptide_table.groupby('pro'):
        pep_of_pro_list = []
        for eatch_pep in each_pro_df[1].pep:
            pep_of_pro_list += id_peptide_dict[eatch_pep]
        id_protein_dict[each_pro_df[0]] = list(set(pep_of_pro_list))
        id_protein_dict_pep[each_pro_df[0]] = list(set(each_pro_df[1].pep))
    protein_tapp_list = [i for i in list(id_protein_dict.keys()) if i]
    protein_table = pd.DataFrame(protein_tapp_list, columns=['pro'])
    #parameter of valcono plot

    def fold_change_caculation(identified_item, p_or_p='pep'):
        if p_or_p == 'pep':
            id_dict = id_peptide_dict
        else:
            id_dict = id_protein_dict
        pro_data_slide = meta_data[id_dict[identified_item]]
        ttest_fold_fc = pro_data_slide.iloc[sample_group_id].sum().sum()/pro_data_slide.iloc[contrl_group_id].sum().sum()
        return np.log2(ttest_fold_fc)

    def ttest_caculation(identified_item, p_or_p='pep'):
        if p_or_p == 'pep':
            id_dict = id_peptide_dict
        else:
            id_dict = id_protein_dict
        pro_data_slide = meta_data[id_dict[identified_item]].T.mean()
        ttest_fold_pv = stats.ttest_ind(pro_data_slide[contrl_group_id], pro_data_slide[sample_group_id])[1]
        return - np.log10(ttest_fold_pv)

    peptide_table['pv'] = peptide_table['pep'].apply(ttest_caculation)
    peptide_table['fc'] = peptide_table['pep'].apply(fold_change_caculation)
    protein_table['pv'] = protein_table['pro'].apply(ttest_caculation, p_or_p='pro')
    protein_table['fc'] = protein_table['pro'].apply(fold_change_caculation, p_or_p='pro')
    #annoting
    ms_annotes = get_annotes(ms_table, win_id=1)
    peptide_annotes = get_annotes(peptide_table, win_id=2)
    protein_annotes = get_annotes(protein_table, win_id=3)

    return (ms_table, ms_annotes), (peptide_table, peptide_annotes), (protein_table, protein_annotes),\
    (id_peptide_dict, id_protein_dict_pep, id_protein_dict, ), meta_data


def ms_data_graph(data_fasta, data_annotes, if_unlabeled=1, if_labeled=1, win_id=1, selected =''):
    x = list(data_fasta['fc'])
    y = list(data_fasta['pv'])
    blue_data = data_fasta.loc[data_fasta['pro'] != '']
    grey_data = data_fasta.loc[data_fasta['pro'] == '']
    blue_label = 'Annotated'
    grey_label = 'Unannotated'
    if win_id == 2:
        blue_label = 'Unique'
        grey_label = 'Not_Unique'
    fig, ax = plt.subplots()
    #fig = Figure(figsize=(5, 4), dpi=100)
    #ax = f.add_subplot(111)
    if win_id != 3:
        if if_unlabeled:
            ax.scatter(grey_data['fc'], grey_data['pv'], c='#A8A7A2', s=my_point_size,zorder=4, alpha=alpha_number, label=grey_label)
        else:
            x = list(blue_data['fc'])
            y = list(blue_data['pv'])
        if if_labeled:
            ax.scatter(blue_data['fc'], blue_data['pv'], c='#383835', s=my_point_size,zorder=5, alpha=alpha_number, label=blue_label)
        else:
            x = list(grey_data['fc'])
            y = list(grey_data['pv'])
    else:
        ax.scatter(data_fasta['fc'], data_fasta['pv'], c='#383835', s=my_point_size,zorder=5, alpha=alpha_number, label='Protein')
    if selected:
        selection = selected.split(';')
        for sel_index in range(len(selection)):
            selection_data = data_fasta.loc[data_fasta['pro'] == selection[sel_index]]
            ax.scatter(selection_data['fc'], selection_data['pv'], c=colors_id[sel_index], s=my_point_size, zorder=50, alpha=alpha_number, label=selection[sel_index])
    af =  AnnoteFinder(x,y, data_annotes, ax=ax, xtol=0.075, ytol=0.075)
    plt.legend(bbox_to_anchor=(0.99, 0.99),loc=1, borderaxespad=0., fontsize='xx-small', edgecolor='white')
    fig.canvas.mpl_connect('button_press_event', af)
    plt.show()


def closestNonPrime(n):
    if n==0:
        return None, None
    def prime(num):
        if num == 1:
            return 1, 1,
        half = int(np.sqrt(num))
        prime_list = list(range(2, half + 1))
        prime_list.sort(reverse=True)
        for i in prime_list:
            if num % i == 0:
                if num/i/i > 3:
                     return int((half+1)), int((half+1))
                return int(i), int(num/i)
        return 0
    while 1:
        c = prime(n)
        if c:
            return c
        n += 1

def box_plot(protein_fasta_id, dicts, meta_data):
    id_peptide_dict, id_protein_dict_pep, id_protein_dict = dicts
    pep_list = id_protein_dict_pep[protein_fasta_id]
    start_sub = closestNonPrime(len(pep_list))
    plt.figure(protein_fasta_id)
    axs = [plt.subplot(start_sub[0],start_sub[1],i+1) for i in range(len(pep_list))]

    for i in range(len(pep_list)):
        peak_list = id_peptide_dict[pep_list[i]]
        contrl_group = meta_data[peak_list].iloc[contrl_group_id].mean()/100000
        sample_group = meta_data[peak_list].iloc[sample_group_id].mean()/100000
        plt.sca(axs[i])
        axs[i] = pd.DataFrame([contrl_group, sample_group]).T.boxplot(fontsize=10, labels=['',''])
        axs[i].set_title(pep_list[i], fontsize=8)
    plt.show()




'''
ms_tables, peptide_tables, protein_tables, dicts, meta_data = data_parsing()

selections = 'P69905;P11684;Q8IXQ6;Q92888'
ms_data_graph(ms_tables[0], ms_tables[1], win_id=1, selected=selections)
ms_data_graph(peptide_tables[0], peptide_tables[1], win_id=2, selected=selections)
ms_data_graph(protein_tables[0], protein_tables[1], win_id=3, selected=selections)

protein_fasta_id = 'P11684'
box_plot(protein_fasta_id, dicts, meta_data)
'''


#!/usr/bin/env python
# coding: utf-8

# In[4]:


import os
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import chisquare
from scipy.stats import multinomial
from scipy.stats import mannwhitneyu
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})


# In[5]:


def count_freq_detection(status_data,threshold,pop_groups):
    pop_dict = {}
    pop_freq_dict = {}
    prev_dict = {}
    pop_prev_dict = {}
    for pop in pop_groups:
        nested_array = status_data.filter(like=pop).values
        pop_array = np.ravel(nested_array)
        pop_array = pop_array[pop_array > threshold]
        pop_prevalence = [1 if max(sub_array) > threshold else 0 for sub_array in nested_array]
        pop_freq = sum(1 for sub_array in nested_array if max(sub_array) > threshold)
        if np.shape(status_data)[0] > 0:
            pop_freq_dict[pop]=pop_freq/np.shape(status_data)[0]
        else:
            pop_freq_dict[pop]=0
        pop_prev_dict[pop]=pop_freq
        prev_dict[pop] = pop_prevalence
        pop_dict[pop]=pop_array
    return(pop_dict,prev_dict,pop_prev_dict,pop_freq_dict)

def jitter_plot_detections(detection_dict,threshold,plot_title,size=(10,6), vertical = False, by = 'status'):
    df = pd.DataFrame(detection_dict)
    melted_df = df.reset_index().melt(id_vars='index', var_name='Condition', value_name='Value')
    exploded_df = melted_df.explode('Value')
    plt.figure(figsize=size)
    if not vertical:
        if by == 'status':
            plt.rcParams['figure.dpi'] = 300
            plt.rcParams['savefig.dpi'] = 300
            sns.stripplot(x='Condition', y='Value', data=exploded_df.reset_index(), hue='index', jitter=True,palette = pop_hues)
            plt.xlabel(r'\textbf{Cohort}')
            plt.ylabel(r'\textbf{Sequence Breadth}')
        elif by == 'pop':
            sns.stripplot(x='index', y='Value', data=exploded_df.reset_index(), hue='Condition', jitter=True)
            sns.set_palette("colorblind")
            plt.xlabel('Cohort')
            plt.ylabel(r'\textbf{Sequence Breadth}')
    else:
        sns.stripplot(x='Value', y='Condition', data=exploded_df.reset_index(), hue='index', jitter=True,palette = pop_hues)
        plt.xlabel('Sequence Breadth')
        plt.ylabel('Cohort')
    plt.title(f'{plot_title}, Detection cutoff {threshold}')
    plt.legend(title='Population',bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(f'jitterplot_{plot_title}.png', bbox_inches='tight')
    plt.close()    

def plot_detection_counts(prev_array_dict,threshold,plot_title,size = (10,6)):
    count_df = pd.DataFrame(prev_array_dict)
    melted_df = count_df.reset_index().melt(id_vars='index', var_name='Condition', value_name='Count')
    plt.figure(figsize=size)
    sns.barplot(x='Condition', y='Count', data=melted_df, hue='index', palette = pop_hues)
    plt.xlabel('Population')
    plt.ylabel('Detection Count')
    plt.title(f'{study}, Detection cutoff {threshold}')
    plt.legend(title='Condition')
    x_max = max(plt.xlim())
    y_center = max(melted_df.Count)
    plt.savefig(f'count_barplot_{plot_title}.png', bbox_inches='tight')
    plt.close()
    
def plot_detection_freq(freq_dict,threshold,plot_title,size = (10,6),flip=False):
    freq_df = pd.DataFrame(freq_dict)
    melted_df = freq_df.reset_index().melt(id_vars='index', var_name='Condition', value_name='Frequency')
    sns.set_palette("colorblind")
    plt.figure(figsize=size)
    if flip:
        sns.barplot(x='Condition', y='Frequency', data=melted_df, hue='index',palette = pop_hues)
    else:
        sns.barplot(x='index', y='Frequency', data=melted_df, hue='Condition')

    plt.xlabel('Population')
    plt.ylabel('Detection Frequency')
    plt.title(f'{study}, Detection cutoff {threshold}')
    plt.legend(title='Condition')
    x_center = np.mean(plt.xlim())
    y_center = max(melted_df.Frequency)
    plt.savefig(f'frequency_barplot_{plot_title}.png', bbox_inches='tight')
    plt.close()
    
def box_plot_detections(detection_dict,threshold,plot_title,size=(10,6)):
    #Plotting detection boxplots
    df = pd.DataFrame(detection_dict)
    melted_df = df.reset_index().melt(id_vars='index', var_name='Condition', value_name='Value')
    exploded_df = melted_df.explode('Value')
    exploded_df["Condition"]= exploded_df["Condition"].astype('str')
    plt.figure(figsize=size)
    sns.boxplot(x='Condition', y='Value', data=exploded_df, hue='index')
    sns.set_palette("colorblind")
    plt.xlabel('Condition')
    plt.ylabel('Value')
    plt.title(f'{study}, Detection cutoff {threshold}')
    plt.legend(title='Population',bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(f'boxplot_{plot_title}.png', bbox_inches='tight')
    plt.close()
    
def plot_cooccurence_matrix(prev_dict,threshold,plot_title,size=(10,6),count = False):
        #Creates cooccurence matrix from prevalance data
        prev_df = pd.DataFrame(prev_dict)
        co_occurrence_matrix = prev_df.T.dot(prev_df)
        total_count = co_occurrence_matrix.values.sum()
    
        # Normalize counts to get frequencies
        co_occurrence_matrix_freq = co_occurrence_matrix / co_occurrence_matrix.values.diagonal()
    
        plt.figure(figsize=size)
        if not count:
            sns.heatmap(co_occurrence_matrix_freq, annot=True, cmap="YlGnBu", fmt=".2f")
        else:
            sns.heatmap(co_occurrence_matrix, annot=True, cmap="YlGnBu", fmt=".2f")
        plt.xlabel("Species")
        plt.ylabel("Species")
        plt.title(f"Species Co-Occurrence Matrix {plot_title} Detection cutoff {threshold}")
        plt.savefig(f'Species Co-Occurrence Matrix {plot_title}.png', bbox_inches='tight')
        plt.close()
    
def mann_whitney_arrays(status_arrays):
    print('Mann-Whitney U Tests:')
    for i in range(0,len(status_arrays)-1):
        mw_status = list(status_arrays.keys())[i]
        for j in range(i+1,len(status_arrays)):
            status_two = list(status_arrays.keys())[j]
            for item in status_arrays[mw_status]:
                U1, p = mannwhitneyu(status_arrays[mw_status][item],status_arrays[status_two][item])
                print(f'{mw_status} vs {status_two}, {item}, p = {p}')
def just_srr(sample_name):
    return(sample_name.split('_')[0])


# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# In[21]:


import pandas as pd
import os
import numpy as np
import sqlite3
import csv


# In[22]:


#os.chdir('Downloads')


# In[ ]:





# In[23]:




def export_table_to_csv(db_file, table_name, csv_file):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # Execute SQL query to select all data from the table
    query = f"SELECT * FROM {table_name}"
    cursor.execute(query)

    # Fetch all rows from the table
    rows = cursor.fetchall()

    # Get column names
    column_names = [description[0] for description in cursor.description]

    # Write the data to a CSV file
    with open(csv_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        
        # Write the column headers
        csv_writer.writerow(column_names)
        
        # Write the data rows
        csv_writer.writerows(rows)

    # Close the database connection
    conn.close()


# In[4]:


def fix_split(defline):
    return(defline.split('_split_')[0])

def fix_contig(defline):
    return(defline.split('_genomic')[0])


# In[25]:


def detection_pop_csv():
#Export detection data from SQL to csv. Read into pandas, remove numbers from splits
    db_file = 'PROFILE.db'
    table_name = 'detection_splits'
    csv_file = 'detection_splits.csv'
    export_table_to_csv(db_file, table_name, csv_file)


# In[26]:


def main():
    dir_list = open('f_nucleatum_study_directories.txt','r').readlines()
    #Populations file should have isolate names and populations
    #Should be two-column csv
    pop_file = 'f_nuc_population_key.csv'
    pop_df = pd.read_csv(pop_file)
    dir_list


    for item in dir_list:
        os.chdir(item.strip())
        if os.path.exists('mean_detections_with_pops.csv'):
            os.chdir('..')
            continue
        detection_pop_csv()
        dataset = pd.read_csv('detection_splits.csv')
        dataset['item'] = dataset['item'].apply(fix_split)
        #Average detection of all splits in a given genome.
        #Attach population data to 
        mean_df = dataset.groupby(['layer', 'item'])['value'].mean().unstack()
        #print(mean_df.transpose())
        #print(mean_df.transpose().set_index)
        #pop_df = pd.read_csv(pop_file)
        #print(pop_df.columns)
        merged_df = pd.merge(pop_df,mean_df.transpose(), right_index = True ,left_on = 'layer', how = 'outer')
        #print(mean_df)
        merged_df.to_csv('mean_detections_with_pops.csv',index=False)
        os.chdir('..')


# In[8]:


if __name__ == '__main__':
    main()


# In[9]:


#mean_df.to_csv('mean_detections.txt',sep='\t')
#db_file = 'PROFILE.db'
#table_name = 'mean_coverage_splits'
#csv_file = 'mean_coverage_splits.csv'
#export_table_to_csv(db_file, table_name, csv_file)
#dataset = pd.read_csv('mean_coverage_splits.csv')
#dataset['item'] = dataset['item'].apply(fix_split)
#mean_df = dataset.groupby(['layer', 'item'])['value'].mean().unstack()
#mean_df.to_csv('mean_coverages.txt',sep='\t')


# In[ ]:





# In[12]:


#os.chdir('..')


# In[ ]:





# In[ ]:





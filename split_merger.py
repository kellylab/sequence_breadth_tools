import pandas as pd
import os
import numpy as np
import sqlite3
import csv


def export_table_to_csv(db_file, table_name, csv_file):
    #Find and extract table from SQLdb
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    query = f"SELECT * FROM {table_name}"
    cursor.execute(query)
    rows = cursor.fetchall()
    column_names = [description[0] for description in cursor.description]
    #Write to CSV
    with open(csv_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(column_names)
        csv_writer.writerows(rows)

    # Close the database connection
    conn.close()


def fix_split(defline):
    return(defline.split('_split_')[0])

def fix_contig(defline):
    return(defline.split('_genomic')[0])




def detection_pop_csv():
#Export detection data from SQL to CSV. Read into pandas, remove numbers from splits
    db_file = 'PROFILE.db'
    table_name = 'detection_splits'
    csv_file = 'detection_splits.csv'
    export_table_to_csv(db_file, table_name, csv_file)


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




if __name__ == '__main__':
    main()






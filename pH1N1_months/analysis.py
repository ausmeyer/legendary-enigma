# core library imports
import re, pandas as pd, numpy as np, random, os, shutil
from dateutil.parser import parse

# bio imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def separate_unpassaged(run_prefix):
    records = list(SeqIO.parse("pH1H1_human_northamerica_all_HA.fasta", "fasta"))
    new_records = []

    for record in records:
        if (bool(re.search('LUNG|P0|OR_|ORIGINAL|CLINICAL|DIRECT|TMK|RMK|RHMK|RII|PMK|R[1–9]|RX|S[1–9]|SX|SIAT|MDCK|C[1–9]|CX|C_[1–9]|M[1–9]|MX|X[1–9]|  ˆ X_$|ˆ S[1–9]_$|   ˆSX_$|SIAT2_SIAT1| SIAT3_SIAT1|not SIAT|SX|S[1–9]', record.description.split('|')[2].upper()))):
            record.seq = record.seq[record.seq.find('atg'):]
        	
            tmp = str(record.seq)
            codon_list = [tmp[i:i+3] for i in range(0, len(tmp), 3)]
            stop_list = [i for i, x in enumerate(codon_list) if x == 'tag' or x == 'tga' or x == 'taa']
            codon_list = codon_list[:stop_list[0]]
            
            record.seq = Seq(''.join(codon_list), generic_dna)
            new_records.append(record)
            
    SeqIO.write(new_records, run_prefix + '.fasta', 'fasta')

def order_records(seq_prefix):
	records = list(SeqIO.parse(seq_prefix + '.fasta', 'fasta'))
	
	full_dates = []
	month_dates = []
	ordered_records = []
	
	for record in records:
		record_date = record.description.split('|')[1].split('-')
		
		record_year = record_date[0].split(' ')[1]
		if len(record_date) > 1:
			record_month = record_date[1]
		if len(record_date) > 2:
			record_day = record_date[2]
		
		# skip records without months and keep data after april 2009
		if len(record_date) < 3:
			continue 
		elif int(record_year) == 2009 and int(record_month) < 4:
			continue
		else:
			ordered_records.append(record)
			
		full_dates.append(record_year + '-' + record_month + '-' + record_day)
		month_dates.append(record_year + '-' + record_month)
	
	sort_order = sorted(range(len(full_dates)), key=full_dates.__getitem__)
	full_dates = [full_dates[i] for i in sort_order]
	month_dates = [month_dates[i] for i in sort_order]
	ordered_records = [ordered_records[i] for i in sort_order]
	
	print('\nNumber of records in this set: ' + str(len(ordered_records)) + '\n')
	
	return(month_dates, ordered_records)

def generate_cut_points(month_dates):
	dates = pd.DataFrame({'dates':month_dates})
	month_counts = pd.crosstab(index = dates['dates'], columns = 'count')
	month_counts = pd.DataFrame({'dates':month_counts.index, 'count':[item for sublist in month_counts.values for item in sublist]})
	
	cumsums = month_counts['count'].cumsum()
	rollingsums = month_counts['count'].rolling(3).sum()
	cut_points = pd.DataFrame({                                                           \
        'dates':month_counts['dates'][2:].reset_index(drop=True),                         \
    	'high_cut_points':cumsums[2:].reset_index(drop=True),                             \
    	'low_cut_points':pd.concat([pd.Series([0]), cumsums[:-3]]).reset_index(drop=True),\
    	'verify_rolling_1':rollingsums[2:].reset_index(drop=True),                        \
    	'verify_rolling_2':cumsums[2:].reset_index(drop=True) - pd.concat([pd.Series([0]), cumsums[:-3]]).reset_index(drop=True)})
    
    # sanity check - make sure the cut points are the same by two different calculation strategies
	if list(map(int, cut_points['verify_rolling_1'])) != list(map(int, cut_points['verify_rolling_2'])):
		raise ValueError('The cut points are not at the correct indices')
    	
	return(cut_points)
    
def split_records(ordered_records, cut_points, sample_size, run_prefix):
	for index, row in cut_points.iterrows():
		sub_records = ordered_records[row['low_cut_points']:row['high_cut_points']]
		
		if len(sub_records) >= sample_size:
			sub_sample = random.sample(sub_records, sample_size)
		else:
			sub_sample = sub_records
		
		SeqIO.write(sub_sample, run_prefix + '/' + row['dates'] + '.fasta', 'fasta')

def main():
    run_prefix = 'unpassaged'
    
    #separate_unpassaged(run_prefix)
    (month_dates, ordered_records) = order_records(run_prefix)
    
    cut_points = generate_cut_points(month_dates)
    
    if os.path.isdir(run_prefix):
    	shutil.rmtree(run_prefix)
    
    os.mkdir(run_prefix)
    split_records(ordered_records, cut_points, sample_size = 25, run_prefix = run_prefix)
    
if __name__ == '__main__':
	main()
#!/usr/bin/env python

import subprocess
import re

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# CLASS: MONITORING
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
class Monitoring:
    def __init__(self,tot_raw):
        self.reads = [0] * 6
        self.reads[0]= tot_raw
        #{"raw","trimmed","phix","ncRNA","mapped","mapped_strict"}
    
    def count_seq_fq(self,fastq_file):
        return int(subprocess.check_output("wc -l %s"%(fastq_file), shell=True).split(' ')[0])/4
    
    def count_raw_reads(self):
        return self.reads[0]
        
    def count_processed_reads(self,trim_log):
        f = open(trim_log,'r')
        for line in f:
            if "Input Read Pairs: " in line: 
                reads = [int(x) for x in line.split(' ') if x.isdigit()]
                #0:raw 1:both_surv 2:fwd_surv 3:rev_surv 4:dropped
                break
        f.close()
        self.reads[1]= sum(reads[1:4])        
        return self.reads[1]        
        
    def count_non_Phix_reads(self, phix_ID_file):
        phiX_reads = int(subprocess.check_output("wc -l "+phix_ID_file, shell=True).split(' ')[0])
        self.reads[2]= self.reads[1]- phiX_reads   
        return self.reads[2]
    
    def count_non_ncRNA_reads(self,seq_paired_filtered, seq_single_filtered):
        pairs  = self.count_seq_fq(seq_paired_filtered)
        singles= self.count_seq_fq(seq_single_filtered)         
        self.reads[3]= pairs + singles  
        return self.reads[3]        
        
    def count_mapping_reads(self,mapping_log,before_filter):
        with open(mapping_log, 'r') as f:
            for line in f:
                if re.search('with itself and mate mapped',line):   
                    pairs_mapped = int(line.split(' ')[0])
                elif re.search('in total', line):
                    total_mapped = int(line.split(' ')[0])
            mapped_reads = pairs_mapped/2 +(total_mapped-pairs_mapped)
            f.close()
        if before_filter:
            self.reads[4]=mapped_reads
            return self.reads[4]
        else:
            self.reads[5]=mapped_reads
            return self.reads[5]

    def get_tot_percentage(self,count_reads):
        tot_percentage = round(float(count_reads)/self.reads[0]*100,2) 
        return (str(tot_percentage)+' %')    

    def get_percentage_prev(self, index):
        if index >0:        
            prev_percentage = round(float(self.reads[index]) /self.reads[index-1]*100,2)
            return str(prev_percentage)
        else: 
            raise Exception ("index must be >=1")
    
    def get_all_tot_percentage(self):
        tot_percentage = [round(float(x)/self.reads[0]*100,2) for x in self.reads]
        return tot_percentage
    
    def get_all_percentage_prev(self):
        prev_percentage =[0] * 6
        prev_percentage[0]=100.0
        for i in range(1,6):
            prev_percentage[i] = round(float(self.reads[i]) /self.reads[i-1]*100,2)
        return prev_percentage

    
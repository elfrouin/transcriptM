#!/usr/bin/env python


# lib
import os
from ruffus import *
import numpy
import ruffus.cmdline as cmdline
import subprocess
import tempfile
import csv
import re
import string
import collections
import shutil



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# LOCAL IMPORT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
from monitoring import Monitoring

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# CLASS PIPELINE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

class Pipeline :
    def __init__(self,args):
        self.args=args
        
        ### db
        # adapters
        if self.args.adapters=='truseq':
            self.adapters= os.path.join(self.args.db_path,'0-Adapters/TruSeq3-PE-2.fa')
        elif self.args.adapters=='nextra':
            self.adapters= os.path.join(self.args.db_path,'0-Adapters/NexteraPE-PE.fa')
        if not os.path.isfile(self.adapters):
            raise Exception ("The subdirectory 0-Adapters/ or the file %s does not exist in %s (db_path provided)"
                            %(os.path.basename(self.adapters),self.args.db_path))
            exit(1)
        # PhiX: reference genome 
        self.ref_genome_phiX= os.path.join(self.args.db_path,'2-PhiX/phiX.fa')
        if not os.path.isfile(self.ref_genome_phiX):
            raise Exception ("The subdirectory 2-PhiX/ or the file %s does not exist in %s (db_path provided)"
                            %(os.path.basename(self.ref_genome_phiX),self.args.db_path))
            exit(1)
        
        # gff files
        self.list_gff = list(numpy.sort(self.get_files(self.args.dir_bins ,'.gff')))
        if len(set([os.path.basename(x) for x in self.list_gff]))<len(self.list_gff):
            raise Exception ("--dir_bins args \nWarning: some gff files have the same name")
            exit(1)

        # prefix
        self.alias_pe = {}        
        for i in range(int(len(self.args.paired_end)/2)) :
            self.alias_pe[os.path.join(self.args.working_dir,'sample-'+str(i)+'_R1.fq.gz')] =self.args.paired_end[2*i]
            self.alias_pe[os.path.join(self.args.working_dir,'sample-'+str(i)+'_R2.fq.gz')] =self.args.paired_end[2*i+1]
        
        self.prefix_pe= {}
        for  i in range(int(len(self.args.paired_end)/2)) :
            value = self.longest_common_substring(os.path.basename(self.args.paired_end[2*i]),
                                             os.path.basename(self.args.paired_end[2*i+1]))
            if value.endswith(('.','_','-'),0):
                value=value[:-1]  
            elif value.endswith(('_R','-R','.R'),0):  
                value=value[:-2]                               
            self.prefix_pe['sample-'+str(i)]=value
        
        if len(set(self.prefix_pe.values()))< int(len(self.args.paired_end)/2):
            print [item for item, count in collections.Counter(self.prefix_pe.values()).items() if count > 1]
            raise Exception ("2 sets of paired-ends files have the same prefix. Rename one set. \n")
            exit(1)
         
        self.tot_pe= {}
        for i in range(int(len(self.args.paired_end)/2)) :
            count = int(subprocess.check_output("zcat %s | wc -l " %(self.args.paired_end[2*i]), shell=True).split(' ')[0])/4 
            self.tot_pe[self.prefix_pe['sample-'+str(i)]]=count
            print  ('\t').join([self.prefix_pe['sample-'+str(i)],'raw data','FastQC-check','raw reads',str(count),'100.00 %'])
        # log
        self.logger, self.logging_mutex = cmdline.setup_logging (__name__, args.log_file, args.verbose)
        

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # MISCELLANEOUS FONCTIONS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    def re_symlink (self,input_file, soft_link_name, logger, logging_mutex):
        """
        Helper function: relinks soft symbolic link if necessary
        """
        # Guard against soft linking to oneself: Disastrous consequences of deleting the original files
        if input_file == soft_link_name:
            logger.debug("Warning: No symbolic link made. You are using the original data directory as the working directory.")
            return
        # Soft link already exists: delete for relink?
        if os.path.lexists(soft_link_name):
        # do not delete or overwrite real (non-soft link) file
            if not os.path.islink(soft_link_name):
                raise Exception("%s exists and is not a link" % soft_link_name)
            try:
                os.unlink(soft_link_name)
            except:
                with logging_mutex:
                    logger.debug("Can't unlink %s" % (soft_link_name))
        with logging_mutex:
            logger.debug("os.symlink(%s, %s)" % (input_file, soft_link_name))
            #   symbolic link relative to original directory so that the entire path
            #       can be moved around with breaking everything
        os.symlink( os.path.relpath(os.path.abspath(input_file), os.path.abspath(os.path.dirname(soft_link_name))), soft_link_name)
    
        
    def get_files(self,directory,extension):
        """
        Helper function: list files with a specific extension in a directory and its subdirectories
        """
        list_files=[] 
        for root, dirs, files in os.walk(directory):
            for f in files:
                if f.endswith(extension):
                    list_files.append(os.path.join(root, f))
        return list_files
    
    def has_index(self,f,list_extension):
        """
        Helper function: check if a file f has index (if its directory also contains files ending with extensions given in a list)
        """
        index=True 
        for ext in list_extension:
            if not os.path.exists(f+ext):
                index=False
                break
        return index
    
    def longest_common_substring(self,S1, S2):
        """
        Helper function: find the longest common substring between 2 strings
        """
        M = [[0]*(1+len(S2)) for i in range(1+len(S1))]
        longest, x_longest = 0, 0
        for x in range(1,1+len(S1)):
            for y in range(1,1+len(S2)):
                if S1[x-1] == S2[y-1]:
                    M[x][y] = M[x-1][y-1] + 1
                    if M[x][y]>longest:
                        longest = M[x][y]
                        x_longest  = x
                else:
                    M[x][y] = 0
        return S1[x_longest-longest: x_longest]

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE STAGES FUNCTION
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #            
    def pipeline_stages(self):     
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE: STEP N_1
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Create symbolic link of inputs files in the working directory
        @mkdir(self.args.working_dir)
        @originate(self.alias_pe.keys(), self.logger, self.logging_mutex)
        def symlink_to_wd_metaT (soft_link_name, logger, logging_mutex):
            """
            Make soft link in working directory
            """
            
            input_file= self.alias_pe[soft_link_name]
            with logging_mutex:
                logger.info("Linking files %(input_file)s -> %(soft_link_name)s" % locals())
            self.re_symlink(input_file, soft_link_name, logger, logging_mutex)
            
            
        @mkdir(self.args.working_dir)        
        @transform(self.args.metaG_contigs, formatter(),
                        # move to working directory
                        os.path.join(self.args.working_dir,"{basename[0]}"+".fa"),
                        self.logger, self.logging_mutex)
        def symlink_to_wd_metaG (input_file, soft_link_name, logger, logging_mutex):
            """
            Make soft link in working directory
            """
            with logging_mutex:
                logger.info("Linking files %(input_file)s -> %(soft_link_name)s" % locals())
            self.re_symlink(input_file, soft_link_name, logger, logging_mutex)
        
        # check if bwa index are present 
        @active_if(self.has_index(self.args.metaG_contigs,['.amb','.bwt','.ann','.pac','.sa']))
        @mkdir(self.args.working_dir)        
        @transform([self.args.metaG_contigs+x for x in ['.amb','.bwt','.ann','.pac','.sa'] ], formatter(),
                        # move to working directory
                        os.path.join(self.args.working_dir,"{basename[0]}{ext[0]}"),
                        self.logger, self.logging_mutex)
        def symlink_to_wd_metaG_index (input_file, soft_link_name, logger, logging_mutex):
            """
            Make soft link in working directory
            """
            with logging_mutex:
                logger.info("Linking files %(input_file)s -> %(soft_link_name)s" % locals())
            self.re_symlink(input_file, soft_link_name, logger, logging_mutex)
            
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE: STEP N_2
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #    
    # First step in the QC process of the raw reads: trimming based on length and quality score     
        @collate(symlink_to_wd_metaT,
                 regex("R[12].fq.gz$"),
                 ["trimm_P1.fq.gz", "trimm_P2.fq.gz", "trimm_U1.fq.gz", "trimm_U2.fq.gz"],"trimmomatic.log",
                 self.logger, self.logging_mutex)
        def trimmomatic(input_files, output_file,log,logger, logging_mutex):
            """
            Trimmomatic. Trim and remove adapters of paired reads
            """   
            if len(input_files) != 2:
                raise Exception("One of read pairs %s missing" % (input_files,))
            cmd= "trimmomatic PE -threads %d -%s %s %s %s %s %s %s ILLUMINACLIP:%s:2:30:10 \
        LEADING:%d SLIDINGWINDOW:4:%d TRAILING:%d CROP:%d HEADCROP:%d MINLEN:%d 2> %s" %(self.args.threads,
                                                                                    self.args.phred,
                                                                                    input_files[0],
                                                                                    input_files[1],
                                                                                    output_file[0],
                                                                                    output_file[2],
                                                                                    output_file[1],
                                                                                    output_file[3],
                                                                                    self.adapters,
                                                                                    self.args.min_qc,
                                                                                    self.args.min_avg_qc,
                                                                                    self.args.min_qc,
                                                                                    self.args.crop,
                                                                                    self.args.headcrop,
                                                                                    self.args.min_len,
                                                                                    log)
            with logging_mutex:
                logger.info("Trim and remove adapters of paired reads of %(input_files)s" % locals())
                logger.debug("trimmomatic: cmdline\n"+ cmd)
            subprocess.check_call(cmd, shell=True)
           
            #  ~~~~ monitoring: count of reads  ~~~~ #  
            name_sample = self.prefix_pe[os.path.basename(input_files[0]).split('_R1.fq.gz')[0]]            
            stat= Monitoring(self.tot_pe[name_sample])
            ## processed reads
            processed_reads = stat.count_processed_reads(log)
            print ('\t').join([name_sample,'trimming','Trimmomatic','raw reads',str(processed_reads),stat.get_tot_percentage(processed_reads)])

                  
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE: STEP N_3
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #    
    # Second step in the QC process: remove phiX contamination
        @subdivide(trimmomatic,formatter(r"(.+)/(?P<BASE>.*)P1.fq.gz"),["{path[0]}/phiX.{BASE[0]}P1.bam",
                                                                        "{path[0]}/phiX.{BASE[0]}U1.bam",
                                                                        "{path[0]}/phiX.{BASE[0]}U2.bam"], self.logger, self.logging_mutex)
        def phiX_map (input_files, output_files, logger, logging_mutex):
            """
            BamM make. Map all reads against PhiX genome
            """
            cmd ="bamm make -d %s -c %s %s -s %s %s -o %s --threads %d -K --quiet" %(self.ref_genome_phiX,
                                                                                 input_files[0],
                                                                                 input_files[1],
                                                                                 input_files[2],
                                                                                 input_files[3],
                                                                                 self.args.working_dir,
                                                                                 self.args.threads)
            with logging_mutex:
                logger.info("Map reads [%s] against phiX genome"%(','.join(input_files)))
                logger.debug("phiX_map: cmdline\n"+ cmd)
            subprocess.check_call(cmd, shell=True)
            
        
        @transform(phiX_map, suffix(".bam"),".txt",self.logger, self.logging_mutex)
        def phiX_ID(input_files,output_files,logger, logging_mutex):
            """
            Samtools. Get the IDs of PhiX reads
            """
            cmd ="samtools view -F4 %s | awk {'print $1'} > %s "%(input_files,output_files)
            with logging_mutex:
                logger.info("Extract ID of phiX reads in %s" %(input_files))
                logger.debug("phiX_ID: cmdline\n"+ cmd)
            subprocess.check_call(cmd, shell=True)

        
        @collate(phiX_ID,formatter(r"phiX.(?P<BASE>.*)[UP][12].txt$"),'{path[0]}/{BASE[0]}phiX_ID.log','{BASE[0]}',self.logger, self.logging_mutex)
        def phiX_concat_ID(input_files, output_file,basename,logger, logging_mutex):
            """
            Concatenate all PhiX ID found previously
            """
            cmd ="cat %s %s %s | uniq > %s" %(input_files[0],
                                        input_files[1],
                                        input_files[2],
                                        output_file)
            with logging_mutex:
                logger.info("Concatenate all ID of phiX reads [%s]"%(','.join(input_files)))
                logger.debug("phiX_concat_ID: cmdline\n"+ cmd)
            subprocess.check_call(cmd, shell=True) 
           
           #  ~~~~ monitoring: count of reads  ~~~~ #                 
            name_sample = self.prefix_pe[os.path.basename(output_file).split('_trimm_phiX_ID.log')[0]]            
            stat= Monitoring(self.tot_pe[name_sample])
            ## non phiX reads
            trimm_file = os.path.join(self.args.working_dir,[f for f in os.listdir(self.args.working_dir) if re.search(r'%s.*trimmomatic.log'%(basename.split('_')[0]), f)][0])     
            processed_reads = stat.count_processed_reads(trimm_file)
            phiX_reads = int(subprocess.check_output("wc -l "+output_file, shell=True).split(' ')[0])
            non_phiX_reads = processed_reads - phiX_reads
            print ('\t').join([name_sample,'PhiX removal','bamM make','processed reads',str(non_phiX_reads),stat.get_tot_percentage(non_phiX_reads)])
 


        @subdivide(trimmomatic,regex(r"trimm_[UP][12].fq.gz"),["trimm_P1.fq.gz", "trimm_P2.fq.gz", "trimm_U1.fq.gz", "trimm_U2.fq.gz"])
        def QC_output(input_files,output_files):
            pass
        
        
        @transform(QC_output,suffix(".fq.gz"),add_inputs(phiX_concat_ID),"_phiX_ext.fq",self.logger, self.logging_mutex)
        def phiX_extract(input_files, output_files,logger, logging_mutex):
            """
            Remove PhiX reads
            """
            try:
                cmd ="fxtract -S -H -f %s -z -v %s > %s" %(input_files[1], input_files[0],output_files)
                with logging_mutex:
                    logger.info("Extract phiX reads in the file %s"%(input_files[0]))
                    logger.debug("phiX_extract: cmdline\n"+ cmd)
                subprocess.check_call(cmd, shell=True) 
            #flag -z if gzip input file
            except subprocess.CalledProcessError:
                cmd ="gzip  -cd %s > %s" %(input_files[0],output_files)
                with logging_mutex:
                    logger.info("No phiX reads in the file: %s"%(input_files[0]))
                    logger.debug("phiX_extract: cmdline\n"+ cmd)
                subprocess.check_call(cmd, shell=True) 

        
        
    
        
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE: STEP N_4
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #    
    # Third step in the QC process: remove rRNA, tRNA ...
        @subdivide(phiX_extract,formatter(), "{path[0]}/{basename[0]}_non_ncRNA.fq",
                   "{path[0]}/{basename[0]}_ncRNA.fq",self.logger, self.logging_mutex)
        def sortmerna (input_files, output_files, ncRNA_files,logger, logging_mutex):    
            """
            SortMeRNA. Remove non-coding RNA
            """
            cmd= "sortmerna --ref %s --reads %s --aligned %s --other %s --fastx -a %d --log" %(self.args.path_db_smr,
                                                                                               input_files,
                                                                                               ncRNA_files.split('.fq')[0],
                                                                                               output_files.split('.fq')[0],
                                                                                               self.args.threads)
            with logging_mutex:
                logger.info("Remove reads with SortMeRNA in %(input_files)s"%locals())
                logger.debug("sortmerna: cmdline\n"+ cmd)
            subprocess.check_call(cmd, shell=True)
    
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE: STEP N_5
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
    # mapping the reads to reference genome       
        #-move the singleton (generated with sortmerna) with the single ones
        @collate(sortmerna,regex(r"trimm_.*"),["concat_paired_R1.fq","concat_paired_R2.fq","concat_single.fq"],
                 "ID_single.txt", self.logger, self.logging_mutex)
        def concat_for_mapping(input_files, output_files,ID_single,logger, logging_mutex):
            """
            Prepare .fq files for the mapping stage
            """
            cmd_ID="comm -3 <(awk '"'{print substr($1,2) ;getline;getline;getline}'"' %s |sort) <(awk '"'{print substr($1,2) ;getline;getline;getline}'"' %s | sort) | awk  '{print $1 $2}' > %s"  %(input_files[0], input_files[1],ID_single)
            cmd_paired= "fxtract -H -S -f '%s' -v %s > %s;fxtract -H -S -f '%s' -v %s > %s"  %(ID_single,
                                                                                               input_files[0],
                                                                                               output_files[0],
                                                                                               ID_single,
                                                                                               input_files[1],
                                                                                               output_files[1])                                                          
            cmd_single="cat %s %s > %s; fxtract -H -S -f '%s' %s %s >> %s"  %(input_files[2],
                                                                              input_files[3],
                                                                              output_files[2],
                                                                              ID_single,
                                                                              input_files[0],
                                                                              input_files[1],
                                                                              output_files[2])
            with logging_mutex:
                logger.info("Find IDs of single reads generated with SortMeRNA in (%s,%s)" %(input_files[0], input_files[1]))
                logger.debug("concat_for_mapping: cmdline\n"+ cmd_ID)            
            subprocess.check_call(['bash','-c',cmd_ID])
            with logging_mutex:
                logger.info("Prepare paired reads files for the mapping (%s,%s)" %(input_files[0], input_files[1]))
                logger.debug("concat_for_mapping: cmdline\n"+ cmd_paired)  
            subprocess.check_call(cmd_paired, shell=True)
            with logging_mutex:
                logger.info("Prepare single reads files for the mapping (%s,%s)" %(input_files[2], input_files[3]))
                logger.debug("concat_for_mapping: cmdline\n"+ cmd_single)  
            subprocess.check_call(cmd_single, shell=True)
        
            #  ~~~~ monitoring: count of reads  ~~~~ #     
            name_sample = self.prefix_pe[os.path.basename(output_files[0]).split('_concat_paired_R1.fq')[0]]            
            stat= Monitoring(self.tot_pe[name_sample])
            ## non ncRNA reads
            non_ncRNA_reads = stat.count_seq_fq(output_files[0])+stat.count_seq_fq(output_files[2])
            print ('\t').join([name_sample,'remove ncRNA','SortMeRNA','filtered reads (1st)',
                               str(non_ncRNA_reads),stat.get_tot_percentage(non_ncRNA_reads)])


        # Map separately paired-end and singletons with 'BamM' and merge the results in one .bam file
        # WARNINGS
        #1. .bam files generated with 'BamM' only contain the mapped reads -> be carful with the interpretation of samtools flagstat
        #2. only one alignment per read is kept: the secondary and supplementary are removed
        @transform(concat_for_mapping,formatter(r"(.+)/(?P<BASE>.*)_concat_paired_R1.fq"),
                   add_inputs(symlink_to_wd_metaG),"{path[0]}/{BASE[0]}.bam",
                   ["{path[0]}/"+os.path.splitext(os.path.basename(self.args.metaG_contigs))[0]+".{basename[0]}.bam",
                    "{path[0]}/"+os.path.splitext(os.path.basename(self.args.metaG_contigs))[0]+".{basename[2]}.bam",
                    "{path[0]}/"+os.path.splitext(os.path.basename(self.args.metaG_contigs))[0]+"{BASE[0]}_merged.bam"],
                    "{path[0]}/{BASE[0]}_mapping.log",self.logger, self.logging_mutex)
        def map2ref (input_files, output_file, bams,flagstat,logger,logging_mutex):
            """
            BamM make. Map all metatranscriptomics reads against metagenomics contigs
            """
            if self.has_index(input_files[1],['.amb','.bwt','.ann','.pac','.sa']):
                # index already exists -> bamm Kept
                index_exists='K'
            else:
                # index doesn't exist yet -> bamm keep        
                index_exists='k'
            cmd ="bamm make -d %s -c %s %s -s %s -o %s --threads %d -%s --quiet" %(input_files[1],
                                                                            input_files[0][0],
                                                                            input_files[0][1],
                                                                            input_files[0][2],
                                                                            self.args.working_dir,
                                                                            self.args.threads,
                                                                            index_exists) 
            with logging_mutex:     
                logger.info("Map reads [%s] to the reference metagenome %s"%(",".join(input_files[0]),input_files[1]))
                logger.debug("map2ref: cmdline\n"+ cmd)  
            subprocess.check_call(cmd, shell=True)
            
            cmd2= "samtools merge -f %s %s %s ; samtools view -b -F2304 %s > %s " %(bams[2],bams[0],bams[1],bams[2],output_file)
            with logging_mutex:     
                logger.info("Concatenate %s and %s"%(bams[0],bams[1]))
                logger.debug("map2ref: cmdline\n"+ cmd2)  
            subprocess.check_call(cmd2, shell=True)
                
            cmd3= "samtools flagstat %s > %s " %(output_file,flagstat)
            with logging_mutex:     
                logger.info("Compute statistics of %(output_file)s (samtools flastat)"%locals())
                logger.debug("map2ref: cmdline\n"+ cmd3)  
            subprocess.check_call(cmd3, shell=True)
    
            #  ~~~~ monitoring: count of reads  ~~~~ #   
            name_sample = self.prefix_pe[os.path.basename(output_file).split('.bam')[0]]            
            stat= Monitoring(self.tot_pe[name_sample])
            ## reads filtered : mapped with high stringency
            mapped_reads = stat.count_mapping_reads(flagstat,True)
            print ('\t').join([name_sample,'alignment','BamM make','filtered reads (2nd)',
                               str(mapped_reads),stat.get_tot_percentage(mapped_reads)])
        
        
        @transform(map2ref,formatter('.bam'),"{path[0]}/{basename[0]}_filtered.bam", 
                   "{path[0]}/{basename[0]}_stringency_filter.log",self.logger, self.logging_mutex)
        def mapping_filter (input_file, output_file,flagstat,logger,logging_mutex):
            """
            BamM filter. Select reads which are mapped with high stringency
            """
            if self.args.no_mapping_filter :  
                pass 
            else :
                cmd= "bamm filter --bamfile %s --percentage_id %f --percentage_aln %f -o %s " %(input_file,
                                                                                          self.args.percentage_id,
                                                                                          self.args.percentage_aln,
                                                                                          self.args.working_dir)
                with logging_mutex:     
                    logger.info("Filter %(input_file)s" % locals())
                    logger.debug("mapping_filter: cmdline\n"+ cmd)  
                subprocess.check_call(cmd, shell=True)

                cmd3= "samtools flagstat %s > %s " %(output_file,flagstat)
                with logging_mutex:     
                    logger.info("Compute statistics of %(input_file)s (samtools flastat)" % locals())
                    logger.debug("mapping_filter: cmdline\n"+ cmd3)  
                subprocess.check_call(cmd3, shell=True)

                #  ~~~~ monitoring: count of reads  ~~~~ #   
                name_sample = self.prefix_pe[os.path.basename(output_file).split('_filtered.bam')[0]]            
                stat= Monitoring(self.tot_pe[name_sample])
                ## reads filtered : mapped with high stringency
                mapped_reads_f = stat.count_mapping_reads(flagstat,False)
                print ('\t').join([name_sample,'.bam filter','BamM filter','mapped reads',
                                   str(mapped_reads_f),stat.get_tot_percentage(mapped_reads_f)])

    
        
        
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE: STEP N_6 (fpkm)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
    # coverage if bins provided
        if self.args.no_mapping_filter :  
            bam_file = map2ref 
        else: 
            bam_file= mapping_filter
        @subdivide(bam_file,formatter(),'{path[0]}/*fpkm.csv','{path[0]}/coverage.csv' ,'{path[0]}/{basename[0]}_lib_size.log',self.args.dir_bins,self.logger, self.logging_mutex)
        def bam2fpkm(input_file, output_file,coverage_file,lib_size_log, dir_bins,logger, logging_mutex):
            """
            Dirseq (compute coverage values) +  coverage2fpkm
            """
    #                                                                       #
    ## add control! contigs in gff files must be present in metaG_contigs  ##
    #                   
            lib_size= int(subprocess.check_output("samtools view -c "+input_file, shell=True))
            f = open(lib_size_log, 'w')
            f.write(str(lib_size)) 
            f.close()                                                                #

            for i in range(len(self.list_gff)):
                cmd ="dirseq --bam %s --gff %s --ignore-directions -q>  %s " %(input_file,
                                                                           self.list_gff[i],
                                                                           coverage_file)
                with logging_mutex:     
                    logger.info("Calculte coverage from %s and %s"%(input_file,self.list_gff[i]))  
                    logger.debug("bam2fpkm: cmdline\n"+ cmd)                                       
                subprocess.check_call(cmd, shell=True)        
                
                if lib_size !=0:
                    cmd1= "sed 's/\t/|/g' %s | awk  -F '|' 'NR>=2 {$6= $6/%d*10e6}1' OFS='|' |  sed 's/|/\t/g' > %s ; rm %s " %(coverage_file,
                                                                                                                       lib_size,
                                                                                                                       input_file.split('.bam')[0]+'_'+ os.path.splitext(os.path.basename((self.list_gff[i])))[0]+'_fpkm.csv',
                                                                                                                       coverage_file )
                else:
                    cmd1= "cp %s %s; rm %s "%(coverage_file,input_file.split('.bam')[0]+'_'+ os.path.splitext(os.path.basename((self.list_gff[i])))[0]+'_fpkm.csv',coverage_file)
                    
                with logging_mutex:     
                    logger.info("Convert coverage to fpkm")
                    logger.debug("bam2fpkm: cmdline\n"+ cmd1)
                subprocess.check_call(cmd1, shell=True)                 
        
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE: STEP N_6 BIS (raw count)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
    # raw count if bins provided
        if self.args.no_mapping_filter :  
            bam_file = map2ref 
        else: 
            bam_file= mapping_filter
        @subdivide(bam_file,formatter(),'{path[0]}/*count.csv',self.list_gff,self.logger, self.logging_mutex)
        def bam2raw_count(input_file, output_file, list_gff,logger, logging_mutex):
            """
            Bedtools (count the number of mapped reads per gene)
            """                                                            #
            for i in range(len(self.list_gff)):
                gff_no_fasta= tempfile.NamedTemporaryFile(prefix='transcriptm', suffix='.gff')
                cmd0 = "sed '/^##FASTA$/,$d' %s > %s" %(self.list_gff[i], gff_no_fasta.name)
                subprocess.check_call(cmd0, shell=True)
                cmd ="bedtools intersect -c -a %s -b %s -bed >  %s " %( gff_no_fasta.name,
                                                                       input_file,
                                                                       input_file.split('.bam')[0]+'_'+os.path.splitext(os.path.basename((self.list_gff[i])))[0] +'_count.csv')
                with logging_mutex:     
                    logger.info("Calculte raw count from %s and %s "%(input_file,gff_no_fasta.name))  
                    logger.debug("bam2raw_count: cmdline\n"+ cmd)                                       
                subprocess.check_call(cmd, shell=True)        

 

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE: STEP N_7 (fpkm table)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
    # Concatenate all the FPKM results in a table
        @mkdir(self.args.output_dir)
        @merge(bam2fpkm,os.path.join(self.args.output_dir,os.path.basename(self.args.output_dir)+'_FPKM.csv'),self.logger, self.logging_mutex)
        def transcriptM_table (input_files, output_file, logger, logging_mutex): 
            """
            Create one table that contains RPKM values for each gene of each bin for the different samples
            """
            input_files=list(set(input_files))          
            fpkm_col= [list([]) for _ in xrange(int(len(self.args.paired_end)/2)+3)]       
            # headers of cols ->  0, n-1, n
            fpkm_col[0].append('bin_ID')
            fpkm_col[-2].append('gene location [contig:start:end]')
            fpkm_col[-1].append('annotation')      
        
            bins_path =[os.path.splitext((self.list_gff[i]))[0] for i in range(len(self.list_gff))]
            for b in bins_path :
                files_b= [f for f in input_files if re.search('_'+os.path.basename(b)+'_fpkm.csv', f)]  
                # first col: bins_name
                with open(files_b[0],'r') as csvfile:                    
                        reader = csv.reader(csvfile, delimiter='\t') 
                        next(reader)  # skip header 
                        for row in reader:
                            fpkm_col[0].append(b+'.gff')
                        csvfile.close() 
                # n-1 col: gene location
                with open(files_b[0],'r') as csvfile:                    
                        reader = csv.reader(csvfile, delimiter='\t') 
                        next(reader) # skip header 
                        for row in reader:
                            fpkm_col[-2].append(row[0]+':'+row[2]+':'+row[3])
                        csvfile.close() 
                # n col: annotation
                with open(files_b[0],'r') as csvfile:                    
                        reader = csv.reader(csvfile, delimiter='\t') 
                        next(reader) # skip header 
                        for row in reader:
                            fpkm_col[-1].append(row[6])
                        csvfile.close() 
                # remaining cols: FPKM
                for i in range(len(files_b)):
                    # create header of FPKM cols
                    if not fpkm_col[i+1]:
                        fpkm_col[i+1].append('FPKM_'+self.prefix_pe[os.path.basename(files_b[i]).split('_')[0]])
                    with open(files_b[i],'r') as csvfile:                    
                        reader = csv.reader(csvfile, delimiter='\t')                      
                        next(reader) # skip header  
                        for row in reader:
                            fpkm_col[i+1].append(row[5])
                        csvfile.close() 
         
            tab = numpy.array(fpkm_col)
            numpy.savetxt(output_file,numpy.transpose(tab),delimiter='\t', fmt="%s") 
            
            with logging_mutex:     
                logger.info("Create table that contains FPKM values for each gene of each bin given as input for the different samples: %s"%(','.join(self.prefix_pe.values())))    
            
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE: STEP N_7 BIS (raw count table)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
    # Concatenate all the raw count in a table
        @mkdir(self.args.output_dir)
        @merge(bam2raw_count,os.path.join(self.args.output_dir,os.path.basename(self.args.output_dir)+'_COUNT.csv'),self.logger, self.logging_mutex)
        def raw_count_table (input_files, output_file, logger, logging_mutex): 
            """
            Create one table that contains raw count values for each gene of each bin for the different samples
            """
            input_files=list(set(input_files))          
            count_col= [list([]) for _ in xrange(int(len(self.args.paired_end)/2)+3)]       
            # headers of cols ->  0, n-1, n
            count_col[0].append('bin_ID')
            count_col[-2].append('gene location [contig:start:end]')
            count_col[-1].append('annotation')      
        
            bins_path =[os.path.splitext((self.list_gff[i]))[0] for i in range(len(self.list_gff))]
            for b in bins_path :
                files_b= [f for f in input_files if re.search('_'+os.path.basename(b)+'_count.csv', f)]  
                # first col: bins_name
                with open(files_b[0],'r') as csvfile:                    
                        reader = csv.reader(csvfile, delimiter='\t') 
                        for row in reader:
                            count_col[0].append(b+'.gff')
                        csvfile.close() 
                # n-1 col: gene location
                with open(files_b[0],'r') as csvfile:                    
                        reader = csv.reader(csvfile, delimiter='\t') 
                        for row in reader:
                            count_col[-2].append(row[0]+':'+row[3]+':'+row[4])
                        csvfile.close() 
                # n col: annotation
                with open(files_b[0],'r') as csvfile:                    
                        reader = csv.reader(csvfile, delimiter='\t') 
                        for row in reader:
                            try: 
                                annotation = [x for x in row[8].split(";") if re.search('product',x)][0].split('=')[-1]
                            except IndexError :
                                annotation = 'unannotated'
                            count_col[-1].append(annotation)
                        csvfile.close() 
                # remaining cols: count
                for i in range(len(files_b)):
                    # create header of count cols
                    if not count_col[i+1]:
                        count_col[i+1].append('Count_'+self.prefix_pe[os.path.basename(files_b[i]).split('_')[0]])
                    with open(files_b[i],'r') as csvfile:                    
                        reader = csv.reader(csvfile, delimiter='\t')                      
                        for row in reader:
                            count_col[i+1].append(row[9])
                        csvfile.close() 
         
            tab = numpy.array(count_col)
            numpy.savetxt(output_file,numpy.transpose(tab),delimiter='\t', fmt="%s") 
            
            with logging_mutex:     
                logger.info("Create table that contains raw count values for each gene of each bin given as input for the different samples: %s"%(','.join(self.prefix_pe.values())))    
    
    
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE: TRACE FILE N_1
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Create the first output: a fastqc report of raw DATA
        subdir_1= os.path.join(self.args.output_dir,"FastQC_raw")
        # clean the dir (of previous run output)
        try:
            shutil.rmtree(subdir_1)
        except OSError:
            pass          
        @follows(bam2fpkm)
        @mkdir(subdir_1)
        @transform(symlink_to_wd_metaT, formatter(),
                        # move to output directory
                        os.path.join(subdir_1,"{basename[0]}"+"_fastqc.zip"),
                        self.logger, self.logging_mutex)
           
        def view_raw_data (input_file, soft_link_name, logger, logging_mutex):
            """
            Create a fastQC report in the ouptut directory
            """
            cmd ="fastqc %s -o %s --threads %d --quiet; rm %s/*.zip" %(input_file,
                                                                      subdir_1, 
                                                                      self.args.threads,
                                                                      subdir_1)
            with logging_mutex:
                logger.info("Create a fastqc report of raw %(input_file)s" % locals())
                logger.debug("view_raw_data: cmdline\n"+cmd)
            subprocess.check_call(cmd, shell=True)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE: TRACE FILE N_2
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Create the second output: a fastqc report of processed DATA
        subdir_2 = os.path.join(self.args.output_dir,"FastQC_processed")
        # clean the dir (of previous run output)
        try:
            shutil.rmtree(subdir_2)
        except OSError:
            pass          
        @follows(bam2fpkm)
        @mkdir(subdir_2)
        @transform(trimmomatic, formatter(),
                        # move to output directory
                        os.path.join(subdir_2,"{basename[0]}"+"_fastqc.zip"),
                        self.logger, self.logging_mutex)
           
        def view_processed_data (input_file, soft_link_name, logger, logging_mutex):
            """
            Create a fastQC report in the ouptut directory
            """
            cmd ="fastqc %s -o %s --threads %d --quiet; rm %s/*.zip" %(' '.join(input_file),
                                                                      subdir_2, 
                                                                      self.args.threads,
                                                                      subdir_2 )

            with logging_mutex:
                logger.info("Create a fastqc report of processed %(input_file)s" % locals())
                logger.debug("view_processed_data: cmdline\n"+cmd)
            subprocess.check_call(cmd, shell=True)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE:TRACE FILE N_3
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
    # get monitoring data
        subdir_3= os.path.join(self.args.output_dir,"log/")
        # clean the dir (of previous run output)
        try:
            shutil.rmtree(subdir_3)
        except OSError:
            pass        
        @follows(bam2fpkm)
        @mkdir(subdir_3)
        @transform(self.args.working_dir+'/*.log', formatter(".log"),  
                   os.path.join(subdir_3,"{basename[0]}"+".log"),
                   self.logger, self.logging_mutex)
        def save_log(input_files, output_files, logger, logging_mutex):
            """
            Save the log files, generated for different stages of the pipeline (in the temp directory)
            """
            cmd= "cp %s %s"   %(input_files, output_files)    
            with logging_mutex:
                logger.info("Save log files: %(input_files)s" % locals())
                logger.debug("save_log: cmdline\n"+cmd)                
            subprocess.check_call(cmd, shell=True)      
     
     
        subdir_4= os.path.join(self.args.output_dir,"reads_distribution") 
        @mkdir(subdir_4)              
        @collate(save_log,formatter(r"/log/(?P<BASE>.*)_((stringency_filter)|(lib_size)|(mapping)|(trimmomatic)|(trimm_((phiX_ID)|((U|P)(1|2)_phiX_ext_ncRNA)))).log$"),
                 subdir_4+"/{BASE[0]}_reads_stat",'{BASE[0]}')
        def logtable (input_files,output_file,basename):
            """
            Sums up the count of reads which are kept after each step 
            """
            name_sample = self.prefix_pe[basename]
            stat = Monitoring(self.tot_pe[name_sample])          
            # raw
            stat.reads[0]= stat.count_raw_reads()
            # processed reads
            trimm_file = [f for f in input_files if re.search(r'trimmomatic.log', f)][0]     
            stat.reads[1]= stat.count_processed_reads(trimm_file)         
            # non phiX reads
            phix_ID_file = [f for f in input_files if re.search(r'phiX_ID.log', f)][0]
            stat.reads[2]= stat.count_non_Phix_reads(phix_ID_file) 
            # non rRNA/tRNA/tmRNA reads
            list_fqfiles= self.get_files(self.args.working_dir,'.fq')
            pairs_filtered= [f for f in list_fqfiles if re.search(r'%s.+concat_paired_R1.fq'%(basename), f)][0]
            singles_filtered= [f for f in list_fqfiles if re.search(r'%s.+concat_single.fq'%(basename), f)][0]
            stat.reads[3]= stat.count_non_ncRNA_reads(pairs_filtered,singles_filtered)    
            # mapped reads
            mapping_log= [f for f in input_files if re.search(r'mapping.log', f)][0]  
            stat.reads[4]= stat.count_mapping_reads(mapping_log,True)      

            # reads mapped with a given stringency
            if self.args.no_mapping_filter :
                # save stat_table
                tab = numpy.array([[name_sample]*5,
                      ["raw data", "trimming","remove PhiX","remove ncRNA","alignment "],
                      ["FastQC-check", "Trimmomatic","bamM make","SortMeRNA","bamM make"],
                      ["raw reads", "raw reads","processed reads","filtered reads","filtered reads"],
                      map(str,stat.reads[:-1]) ,
                      map(str,stat.get_all_tot_percentage()[:-1]) ,
                      map(str,stat.get_all_percentage_prev()[:-1])])
                numpy.savetxt(output_file,numpy.transpose(tab),delimiter='\t', fmt="%s")                                  
            else:
                stringency_filter_log= [f for f in input_files if re.search(r'stringency_filter.log', f)][0]  
                stat.reads[5]=stat.count_mapping_reads(stringency_filter_log,False)                                   
                # save stat_table
                tab = numpy.array([[name_sample]*6,
                      ["raw data", "trimming","remove PhiX","remove ncRNA","alignment ",".bam filter"],
                      ["FastQC-check", "Trimmomatic","bamM make","SortMeRNA","bamM make","bamM filter"],
                      ["raw reads", "raw reads","processed reads","filtered reads","filtered reads","mapped reads"],
                      map(str,stat.reads) ,
                      map(str,stat.get_all_tot_percentage()) ,
                      map(str,stat.get_all_percentage_prev())])
                numpy.savetxt(output_file,numpy.transpose(tab),delimiter='\t', fmt="%s") 
       

        @merge(logtable,os.path.join(self.args.output_dir,'summary_reads'),self.logger, self.logging_mutex)
        def concatenate_logtables (input_files,output_file,logger,logging_mutex):
            """
            Concatenate the summuries of reads distribution from all samples
            """
            input_files = (' ').join(input_files)
            h= ["sample name","step name","tool used","input data","reads count","% total","% previous step"]  
            header= ('\t').join(h)            
            cmd= "cat %s > %s ;sed -i '1i%s' %s"   %(input_files, output_file,header,output_file)    
            with logging_mutex:
                logger.info("Concatenate summaries: %(input_files)s" % locals())
                logger.debug("concatenate_logtables: cmdline\n"+cmd)                
            subprocess.check_call(cmd, shell=True)  
            
            
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE:TRACE FILE N_4
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
    # save the processed metaT reads
        subdir_4 = os.path.join(self.args.output_dir,"processed_reads/")
        # clean the dir (of previous run output)
        try:
            shutil.rmtree(subdir_4)
        except OSError:
            pass  
        @mkdir(subdir_4)
        @transform(concat_for_mapping, formatter(),
                        # move to output directory
                        [os.path.join(subdir_4,"{basename[0]}{ext[0]}"),
                        os.path.join(subdir_4,"{basename[1]}{ext[1]}"),
                        os.path.join(subdir_4,"{basename[2]}{ext[2]}")],
                        self.logger, self.logging_mutex)
           
        def save_processed_reads (input_file, output_file, logger, logging_mutex):
            """
            Copy the processed reads in the ouptut directory
            """
            for i in range(3):
                cmd ="cp %s %s " %(input_file[i],output_file[i] )

                with logging_mutex:
                    logger.info("Copy the processed reads %s in the ouptut directory" %(input_file[i]))
                    logger.debug("save_processed_reads: cmdline\n"+cmd)
                subprocess.check_call(cmd, shell=True)
    


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # PIPELINE: RUN
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #         
        cmdline.run(self.args)
        
    
    
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # CLEAR FUNCTION: rename files ...
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #            
    def rename_files(self, directory):
    # rename all files in a directory with their 'real name'
        for f in os.listdir(directory):
            f_new=os.path.join(directory,string.replace(f,f.split('_')[0],self.prefix_pe[f.split('_')[0]]) )   
            f= os.path.join(directory ,f)          
            os.rename(f,f_new)        
        
    def clear(self):  
    # rename logfile
        log_dir = os.path.join(self.args.output_dir,"log/")
        self.rename_files(log_dir)          
    # rename fastqc report
        fastqc_raw = os.path.join(self.args.output_dir,"FastQC_raw/")
        self.rename_files(fastqc_raw) 
    # rename fastqc report
        fastqc_processed = os.path.join(self.args.output_dir,"FastQC_processed/")
        self.rename_files(fastqc_processed)
    # rename processed reads
        processed_dir =os.path.join(self.args.output_dir,"processed_reads/")
        self.rename_files(processed_dir)
        # clean dir
        reads_distrib_dir = os.path.join(self.args.output_dir,"reads_distribution/")
        try:
            shutil.rmtree(reads_distrib_dir)
        except OSError:
            pass    
        
        
        

#!/usr/bin/env python 

""" Implements the pipeline object for TranscriptM: a rapid-throughput pipeline of meta-transcriptomes. """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
debug = False
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Python Standard Library modules
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import collections   # high-performance container datatypes.
import csv           # Comma-Separated-Values file reading and writing.
import json          # JSON, an open-standard format that uses human-readable text to transmit data objects.
import os            # miscellaneous operating system interfaces.
import re            # regular expression operations.
import shutil        # high-level file operations.
import string        # common string operations.
import subprocess    # spawn new processes, connect to their input/output/error pipes.
import tempfile      # generates temporary files/dirs; works on all supported platforms.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# External, non-standard modules
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import extern        # version of Python's `subprocess`. See https://pypi.python.org/pypi/extern/0.1.0
import numpy         # array processing for numbers, strings, records, objects. See http://www.numpy.org/
import ruffus        # light-weight computational pipeline management. See http://www.ruffus.org.uk
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Locally written modules
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from monitoring import Monitoring  # implements class `Monitoring`.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exposed constants.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class const(object):
    """ exposes constants, such as valid stages of the pipeline. """
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """ Stage number values reflect stage order in the pipeline order,  
        and are expected to be positive integers or floats with <= 1 decimal place. """
    _valid_stages_dict = {
        "view_raw_reads":   1,
        "trim_raw_reads":   2,
        "phiX_removal":     3,
        "filter_rna":       4,
        "prep_for_mapping": 5.2,
        "map_to_reference": 5.5,
        "bam_stats":        6,
        "summary_tables":   7 }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @staticmethod
    def is_valid_stage(stage_name):
        """ boolean: whether the stage_name string represent a valid stage. """
        return stage_name in const._valid_stages_dict        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @staticmethod
    def valid_stages_str():
        """ returns a single string of valid stages sorted in pipeline order. """
        d = const._valid_stages_dict
        return ", ".join(sorted(d, key=d.__getitem__))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @staticmethod
    def stage_num(stage_name):
        """ returns the pipeline order number of the stage_name string. """
        return const._valid_stages_dict[stage_name]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The pipeline class
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class pipeline_object:
    """ a ruffus pipeline implementing the stages of TranscriptM. """
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # helpers: restarting and halting of stages.
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def require_run_stage(self, stage_num):
        ans = not ( \
            (self.restart_at and self.restart_at > round(stage_num, 1)) or \
            (self.halt_after and self.halt_after < round(stage_num, 1)))
        #print "*** testing stage {0:<3} : {1}".format(stage_num, "Run" if ans else "(no run)")
        return ans
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def require_stage_restart(self, stage_num):
        ans = self.restart_at and (self.restart_at == round(stage_num, 1))
        #print "*** testing stage {0:<3} : {1}".format(stage_num, "Restart from here" if ans else "(no restart)")
        return ans
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # helpers: restarting and halting of stages.
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, args):
        """ (automatically executed upon instantiation of an instance of this class). """
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def V_restart_halt_stages():
            """ initializes `self.halt_after` and `self.restart_at` stage number floats """
            self.halt_after = None
            self.restart_at = None
            if self.args.halt_after_stage:
                self.halt_after = float(const.stage_num(self.args.halt_after_stage))
            if self.args.restart_from_stage:
                self.restart_at = float(const.stage_num(self.args.restart_from_stage))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def V_metaG_contigs():
            """ if the metagenome contigs are needed for a required pipeline stage, checks that the file is specified
                and exists. 
            """
            if (self.restart_at != None and self.restart_at > const.stage_num("map_to_reference")) or \
                    (self.halt_after != None and self.halt_after < const.stage_num("map_to_reference")):
                pass
            else:
                if not self.args.metaG_contigs: 
                    raise Exception(
                        "error: argument --metaG_contigs is required when the pipeline "\
                        "includes stages that map to the reference genome/contigs.")
                if not os.path.isfile(self.args.metaG_contigs): 
                    raise Exception(
                        "File <{0}> specified for metaG_contigs does not exist.".format(self.args.metaG_contigs))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def V_dir_bins_list_gff():
            """ checks the dir_bins argment and instantiates and validates `self.list_gff`: the list of gff files. 
                if needed for a required pipeline stage.
            """
            if (self.restart_at != None and self.restart_at > const.stage_num("summary_tables")) or \
                    (self.halt_after != None and self.halt_after < const.stage_num("bam_stats")):
                pass
            else:
                if not self.args.dir_bins: 
                    raise Exception(
                        "error: argument --dir_bins is required when the pipeline "\
                        "includes stages that require it.")
                self.list_gff = list(numpy.sort(self.get_files(self.args.dir_bins, '.gff')))
                if len(set([os.path.basename(x) for x in self.list_gff])) < len(self.list_gff):
                    raise Exception ("--dir_bins args \nWarning: some gff files have the same name")
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def V_ref_genome_phiX():
            """ instantiates and validates `self.ref_genome_phiX`: the filepath of the reference genome file (PhiX)
                if needed for a required pipeline stage.
            """
            if (self.restart_at != None and self.restart_at > const.stage_num("phiX_removal")) or \
                    (self.halt_after != None and self.halt_after < const.stage_num("phiX_removal")):
                pass
            else:
                self.ref_genome_phiX = os.path.join(self.args.db_path, '2-PhiX/phiX.fa')
                if not os.path.isfile(self.ref_genome_phiX):
                    raise Exception(
                        "The subdirectory 2-PhiX/ or the file %s does not exist in %s (db_path provided)"\
                        %(os.path.basename(self.ref_genome_phiX), self.args.db_path))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def V_alias_pe():
            """ populates dict `self.alias_pe` of aliases symbolically linked filenames to the paired-end reads.
                key: working directory filenames (but these files will not yet be created).
                       format: `sample-n_R1.fq.gz`, `sample-n_R2.fq.gz`, 
                       where 'n' is a zero-based counter for each pair of filenames. 
                value: original filenames (of metatranscriptomic reads), as listed in argument `paired_end`.
                e.g. {'/srv/home/s4293029/tm_working/sample-0_R1.fq.gz': 
                        '/srv/home/s4293029/tm_data/20120800_P2M.1.fq.gz', ...}
            """
            self.alias_pe = {}  
            for i in range(int(len(self.args.paired_end)/2)): # Loop an index through the pairs in the `paired_end` list.
                self.alias_pe[os.path.join(self.args.working_dir, 'sample-'+str(i)+'_R1.fq.gz')] = \
                    self.args.paired_end[2*i]
                self.alias_pe[os.path.join(self.args.working_dir, 'sample-'+str(i)+'_R2.fq.gz')] = \
                    self.args.paired_end[2*i+1]
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def V_prefix_pe():
            """ populates dictionary `self.prefix_pe`:
                index: "sample-n", where 'n' is the zero-based counter for each *pair* of filenames passed in
                `paired_end`. 
                value: descriptor of the longest common substring (trimmed of some types of trailing characters) 
                for the pair of filenames.
                e.g. {'sample-0': '20120800_P2M', ...}
                Dependencies: self.args.paired_end.
            """
            self.prefix_pe = {} # Initialize the variable to populate.
            # Loop through the paired-end read filenames, in pairs.
            for  i in range(int(len(self.args.paired_end)/2)):
                # Find the longest common substring between the pair of filenames.
                value = self.longest_common_substring(
                    os.path.basename(self.args.paired_end[2*i]),
                    os.path.basename(self.args.paired_end[2*i+1]))
                # Trim the `value` to be rid of various trailing characters:
                if value.endswith(('.', '_', '-')):
                    value=value[:-1]  
                elif value.endswith(('_R', '-R', '.R')):
                    value=value[:-2]                               
                self.prefix_pe['sample-' + str(i)] = value
            # Verify that no pairs were parsed down to the same descriptor:
            if len(set(self.prefix_pe.values())) < int(len(self.args.paired_end)/2):
                print [item for item, count in collections.Counter(self.prefix_pe.values()).items() if count > 1]
                raise Exception ("Two sets of paired-ends files have the same prefix. Rename one set.")
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def V_SUBDIRS():
            """ sets up constant strings representing basename strings of subdirectories to be created 
                under the output directory. 
            """
            # names of subdirectories:
            self.SUBDIR_log                = os.path.join(self.args.output_dir, "log") # 6
            self.SUBDIR_FastQC_raw         = os.path.join(self.args.output_dir, "FastQC_raw") # 1
            self.SUBDIR_FastQC_processed   = os.path.join(self.args.output_dir, "FastQC_processed") # 2 
            self.SUBDIR_reads_distribution = os.path.join(self.args.output_dir, "reads_distribution") # 6
            self.SUBDIR_processed_reads    = os.path.join(self.args.output_dir, "processed_reads") # 5.2
            # list of subdirectories to have their contents processed at the end of the pipeline.
            self.SUBDIRS_for_content_renaming = [
                self.SUBDIR_log,
                self.SUBDIR_FastQC_raw,
                self.SUBDIR_FastQC_processed,
                self.SUBDIR_processed_reads
                ]
            # list of subdirectories to be deleted at the end of the pipeline.
            self.SUBDIRS_for_clearing = [
                self.SUBDIR_reads_distribution
                ]
            # Clean the subdirs (of previous run output) if the relevant stage is to be run.
            for (stage_name, subdir) in [
                    ("bam_stats", self.SUBDIR_log), 
                    ("view_raw_reads", self.SUBDIR_FastQC_raw),
                    ("trim_raw_reads", self.SUBDIR_FastQC_processed),
                    ("bam_stats", self.SUBDIR_reads_distribution),
                    ("prep_for_mapping", self.SUBDIR_processed_reads)]:
#                print (stage_name, subdir)
                if self.require_run_stage(const.stage_num(stage_name)):
                    try:
                        shutil.rmtree(subdir)
                        print "deleted:", subdir
                    except OSError:
                        pass          
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def setup_logging():
            """ sets up Ruffus' standard python logger, which can be synchronised across concurrent Ruffus tasks.
                Vars `self.logger` and `self.logging_mutex` are meant to be passed as parameters to each Ruffus job.  
                    `logger`: forwards logging calls across jobs.
                    `logging_mutex`: prevents different jobs which are logging simultaneously from being jumbled up.
            """   
            #import logging
            
            #self.logger, self.logging_mutex = ruffus.cmdline.setup_logging(__name__, args.log_file, args.verbose)
#            (self.logger, self.logging_mutex) = ruffus.cmdline.setup_logging(
#                    module_name = __name__, 
#                    log_file_name = args.log_file, # by default, Ruffus logs to STDERR.
#                    #verbose = 0)
#                    verbose = args.verbose) # By default, 1: Ruffus logs each task and each job as it is completed.
            #self.logger.setLevel = logging.INFO #logging.DEBUG
            
            
#            import ruffus.proxy_logger
#            (self.logger, self.logging_mutex) = ruffus.proxy_logger.make_shared_logger_and_proxy(
#                    ruffus.proxy_logger.setup_std_shared_logger,
#                    "my_logger",
#                    {"file_name" : args.log_file})
                    
            #self.logger, self.logging_mutex = ruffus.cmdline.setup_logging(__name__, args.log_file, int(self.args.verbose[0]))
            #print "self.args.verbose:", self.args.verbose, type(self.args.verbose)
            #print "self.args.verbose:", self.args.verbose, type(args.verbose)
            #self.logger, self.logging_mutex = ruffus.cmdline.setup_logging(__name__, args.log_file, self.args.verbose)

            #verbosity = 0 or int(self.args.verbose[0])
            #self.logger, self.logging_mutex = ruffus.cmdline.setup_logging(__name__, args.log_file, verbosity)

            #print type(self.logger)
            #self.logger.setLevel(logging.INFO) #logging.DEBUG            
            
            (self.logger, self.logging_mutex) = ruffus.cmdline.setup_logging(__name__, args.log_file, args.verbose)

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ''' function control: '''
        # Store arguments passed:
        self.args = args             
        # Validate and set up arguments:
        V_restart_halt_stages()
        V_metaG_contigs()
        V_dir_bins_list_gff()
        V_ref_genome_phiX()
        V_alias_pe()
        V_prefix_pe()
        # Handle output subdirectories:
        V_SUBDIRS()
        # logging:
        setup_logging()              
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # helpers
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def re_symlink(self, input_file, soft_link_name, logger, logging_mutex):
        """ (helper): relinks soft symbolic link if necessary.
            (This function from the ruffus website: http://www.ruffus.org.uk/faq.html?highlight=re_symlink).
        """
        # Guard against soft linking to oneself: Disastrous consequences of deleting the original files
        if input_file == soft_link_name:
            logger.debug("Warning: No symbolic link made. " + \
                "You are using the original data directory as the working directory.")
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
        # symbolic link relative to original dir so that the entire path can be moved around with breaking everything.
        os.symlink(
            os.path.relpath(os.path.abspath(input_file), os.path.abspath(os.path.dirname(soft_link_name))), 
            soft_link_name)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def get_files(self, directory, extension):
        """ (helper): returns a list of files with a specific extension in a directory and its subdirectories. """
        list_files = [] 
        for root, dirs, files in os.walk(directory):
            for f in files:
                if f.endswith(extension):
                    list_files.append(os.path.join(root, f))
        return list_files
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def has_index(self, f, list_extensions):
        """ (helper): check if a file f has index,
            (if its directory also contains files ending with extensions given in a list).
        """
        index = True 
        for ext in list_extensions:
            if not os.path.exists(f + ext):
                index = False
                break
        return index
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def longest_common_substring(self, S1, S2):
        """ (helper): returns the longest common substring between 2 strings. """
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
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def prt_progress(self, name_sample, p2, p3, p4, p5, p6):
        """ prints out the provided parameters as a line formated in standard column widths.
            This function probably belongs in module Monitoring, but is put here until that is sorted out.
        """
        print "{0:21} {1:16} {2:16} {3:20} {4:>12}  {5:>8}".format(name_sample, p2, p3, p4, p5, p6)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Pipeline Stages
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def build_pipeline_stages(self):
        """ builds the pipeline by adding ruffus tasks; and returns the built pipeline. """
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Implementation of individual pipeline stages
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 0: count_raw_reads
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """ This stage counts the total number of reads in each pair of paired-end-reads files. """
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def calc_tot_pe():
            """ populates dictionary `self.tot_pe` (which is used to monitor stats as the pipelin is run),
                    (index: common substring of pair of paired_end filenames; value: count of reads) 
                    e.g. {'20120800_P2M': 815272, ...},
                and prints that info to the console.
                Dependencies: `self.args.paired_end`, `self.prefix_pe`.
            """
            self.tot_pe = {}
            self.prt_progress("sample_name", "step_name", "tool_used", "input_data", "reads_count", " %_total")
            self.prt_progress("-----------", "---------", "---------", "----------", "-----------", "--------")
            for i in range(int(len(self.args.paired_end)/2)):
                count = int(subprocess.check_output(
                                "zcat %s | wc -l " % (self.args.paired_end[2*i]), 
                                shell=True)
                            .split(' ')[0])/4  #! is this method always accurate, e.g. if extra comment lines in file?
                self.tot_pe[self.prefix_pe['sample-'+str(i)]]=count
                self.prt_progress(\
                    self.prefix_pe['sample-'+str(i)], 'raw data', 'zcat | wc -l', 'raw reads', str(count), '100.00 %')
#                    self.prefix_pe['sample-'+str(i)], 'raw data', 'FastQC-check', 'raw reads', str(count), '100.00 %')
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 1: view_raw_reads
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """ This stage creates soft links (in the working directory) to the raw data files; 
            and outputs a FastQC report for each pair of paired-end-reads files. """        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def symlink_metaT(soft_link_name, logger, logging_mutex):
            """ Make soft link in working directory. """
            input_file = self.alias_pe[soft_link_name]
            with logging_mutex:
                logger.info("Linking files %(input_file)s -> %(soft_link_name)s"%locals())
            self.re_symlink(input_file, soft_link_name, logger, logging_mutex)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def view_raw_data(input_file, soft_link_name, logger, logging_mutex):
            """ generates a FastQC report for the raw data. """
            cmd = "fastqc %s -o %s --threads %d --quiet; rm %s/*.zip" %(
                input_file,
                self.SUBDIR_FastQC_raw, #! Consider passing as function parameter.
                self.args.threads,
                self.SUBDIR_FastQC_raw) #! Consider passing as function parameter.
            with logging_mutex:
                logger.info("Create a fastqc report of raw %(input_file)s"%locals())
                logger.debug("view_raw_data: cmdline\n" + cmd)
            extern.run(cmd)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 2: trim_raw_reads
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """ This stage does a trimmomatic removal of adapters of the paired reads; 
            and outputs a FastQC report for each pair of paired-end-reads files. """        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def trimmomatic(input_files, output_file, log, logger, logging_mutex):
            """ Trimmomatic. Trim and remove adapters of paired reads. """
            if len(input_files) != 2:
                raise Exception("One of read pairs %s missing" % (input_files,))  
            cmd = "trimmomatic PE "
            cmd += "-threads %d " % (self.args.threads)
            cmd += "-%s " % (self.args.phred)
            cmd += "%s " % (input_files[0])
            cmd += "%s " % (input_files[1])
            cmd += "%s " % (output_file[0])
            cmd += "%s " % (output_file[2])
            cmd += "%s " % (output_file[1])
            cmd += "%s " % (output_file[3])
            cmd += "ILLUMINACLIP:%s:2:30:10         " % (self.args.adaptersFile)
            cmd += "LEADING:%d " % (self.args.min_qc)
            cmd += "SLIDINGWINDOW:4:%d " % (self.args.min_avg_qc)
            cmd += "TRAILING:%d " % (self.args.min_qc)
            cmd += "CROP:%d " % (self.args.crop)
            cmd += "HEADCROP:%d " % (self.args.headcrop)
            cmd += "MINLEN:%d " % (self.args.min_len)
            cmd += "2> %s" % (log)
            with logging_mutex:
                logger.info("Trim and remove adapters of paired reads of %(input_files)s" % locals())
                logger.debug("trimmomatic: cmdline\n"+ cmd)
            extern.run(cmd)
            # monitoring: count of reads
            name_sample = self.prefix_pe[os.path.basename(input_files[0]).split('_R1.fq.gz')[0]]            
            stat = Monitoring(self.tot_pe[name_sample])
            ## processed reads
            processed_reads = stat.count_processed_reads(log)
            self.prt_progress(
                    name_sample, 'trimming', 'Trimmomatic', 'raw reads', 
                    str(processed_reads), stat.get_tot_percentage(processed_reads))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def view_processed_data(input_file, soft_link_name, logger, logging_mutex):
            """ Create a fastQC report in the output directory. """
            cmd = "fastqc %s -o %s --threads %d --quiet; rm %s/*.zip" %(
                    ' '.join(input_file),
                    self.SUBDIR_FastQC_processed,
                    self.args.threads,
                    self.SUBDIR_FastQC_processed)
            with logging_mutex:
                logger.info("Create a fastqc report of processed %(input_file)s" %locals())
                logger.debug("view_processed_data: cmdline\n"+cmd)  
            extern.run(cmd)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 3: phiX_removal
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """ This stage ... """        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def phiX_map(input_files, output_files, logger, logging_mutex):
            """ BamM make. Map all reads against PhiX genome. """
            cmd = "bamm make -d %s -c %s %s -s %s %s -o %s --threads %d -K --quiet" %( \
                    self.ref_genome_phiX,
                    input_files[0],
                    input_files[1],
                    input_files[2],
                    input_files[3],
                    self.args.working_dir,
                    self.args.threads)
            with logging_mutex:
                logger.info("Map reads [%s] against phiX genome"%(','.join(input_files)))
                logger.debug("phiX_map: cmdline\n"+ cmd)
            extern.run(cmd)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def phiX_ID(input_files, output_files, logger, logging_mutex):
            """ Samtools. Get the IDs of PhiX reads. """
            cmd = "samtools view -F4 %s | awk {'print $1'} > %s "%(input_files, output_files)
            with logging_mutex:
                logger.info("Extract ID of phiX reads in %s" %(input_files))
                logger.debug("phiX_ID: cmdline\n" + cmd)
            extern.run(cmd)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def phiX_concat_ID(input_files, output_file, basename, logger, logging_mutex):
            """ Concatenate all PhiX ID found previously. """
            cmd = "cat %s %s %s | uniq > %s" %(
                    input_files[0],
                    input_files[1],
                    input_files[2],
                    output_file)
            with logging_mutex:
                logger.info("Concatenate all ID of phiX reads [%s]"%(','.join(input_files)))
                logger.debug("phiX_concat_ID: cmdline\n" + cmd)
            extern.run(cmd) 
           #  monitoring: count of reads          
            name_sample = self.prefix_pe[os.path.basename(output_file).split('_trimm_phiX_ID.log')[0]]            
            stat= Monitoring(self.tot_pe[name_sample])
            ## non phiX reads
            trimm_file = os.path.join(self.args.working_dir,
                  [f for f in os.listdir(self.args.working_dir) \
                      if re.search(r'%s.*trimmomatic.log'%(basename.split('_')[0]), f)][0])     
            processed_reads = stat.count_processed_reads(trimm_file)
            phiX_reads = int(subprocess.check_output("wc -l "+output_file, shell=True).split(' ')[0])
            non_phiX_reads = processed_reads - phiX_reads
            self.prt_progress(
                name_sample, 'PhiX removal', 'bamM make', 'processed reads',
                str(non_phiX_reads), stat.get_tot_percentage(non_phiX_reads))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def QC_output(input_files, output_files):
            pass
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def phiX_extract(input_files, output_files, logger, logging_mutex):
            """ Remove PhiX reads """
            try:
                cmd = "fxtract -S -H -f %s -z -v %s > %s" %(input_files[1], input_files[0], output_files)
                with logging_mutex:
                    logger.info("Extract phiX reads in the file %s"%(input_files[0]))
                    logger.debug("phiX_extract: cmdline\n"+ cmd)
                extern.run(cmd) 
            #flag -z if gzip input file
            except subprocess.CalledProcessError:
                cmd ="gzip  -cd %s > %s" %(input_files[0], output_files)
                with logging_mutex:
                    logger.info("No phiX reads in the file: %s"%(input_files[0]))
                    logger.debug("phiX_extract: cmdline\n"+ cmd)
                extern.run(cmd) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 4: filter_rna
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """ This stage ... third step in the QC process: remove rRNA, tRNA ..."""        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def sortmerna(input_files, output_files, ncRNA_files, logger, logging_mutex):
            """ SortMeRNA. Remove non-coding RNA. """
            cmd = "sortmerna --ref %s --reads %s --aligned %s --other %s --fastx -a %d --log" %( \
                    self.args.path_db_smr,
                    input_files,
                    ncRNA_files.split('.fq')[0],
                    output_files.split('.fq')[0],
                    self.args.threads)
            with logging_mutex:
                logger.info("Remove reads with SortMeRNA in %(input_files)s"%locals())
                logger.debug("sortmerna: cmdline\n"+ cmd)
            extern.run(cmd)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 5.2: prep_for_mapping
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """ This stage ... 
            and saves those processed reads into the output directory. """        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def concat_for_mapping(input_files, output_files, ID_single, logger, logging_mutex):
            """ Prepare .fq files for the mapping stage. """
#            print "************** concat_for_mapping locals():", locals()
            cmd_ID = "comm -3 <(awk '"'{print substr($1,2) ;getline;getline;getline}'"' %s |sort) <(awk '"'{print substr($1,2) ;getline;getline;getline}'"' %s | sort) | awk  '{print $1 $2}' > %s" % \
                    (input_files[0], input_files[1], ID_single)            
#            cmd_ID = "comm -3 <(awk '"'{print substr($1,2) ;getline;getline;getline}'"' %s |sort) " + \
#                    "<(awk '"'{print substr($1,2) ;getline;getline;getline}'"' %s | sort) | awk  '{print $1 $2}' > %s"%(\
#                    input_files[0], 
#                    input_files[1], 
#                    ID_single)
            cmd_paired = "fxtract -H -S -f '%s' -v %s > %s;fxtract -H -S -f '%s' -v %s > %s"%(
                    ID_single,
                    input_files[0],
                    output_files[0],
                    ID_single,
                    input_files[1],
                    output_files[1])                                                          
            cmd_single = "cat %s %s > %s; fxtract -H -S -f '%s' %s %s >> %s" %(
                    input_files[2],
                    input_files[3],
                    output_files[2],
                    ID_single,
                    input_files[0],
                    input_files[1],
                    output_files[2])
            with logging_mutex:
                logger.info("Find IDs of single reads generated with SortMeRNA in (%s,%s)" %(
                    input_files[0], input_files[1]))
                logger.debug("concat_for_mapping: cmdline\n"+ cmd_ID)            
            #subprocess.check_call(['bash','-c',cmd_ID])
            extern.run(cmd_ID)
            with logging_mutex:
                logger.info("Prepare paired reads files for the mapping (%s,%s)" %(input_files[0], input_files[1]))
                logger.debug("concat_for_mapping: cmdline\n"+ cmd_paired)  
            subprocess.check_call(cmd_paired, shell=True)
            #extern.run(cmd_paired)
            with logging_mutex:
                logger.info("Prepare single reads files for the mapping (%s,%s)"%(input_files[2], input_files[3]))
                logger.debug("concat_for_mapping: cmdline\n"+ cmd_single)  
            subprocess.check_call(cmd_single, shell=True)
            #extern.run(cmd_single)
            #  monitoring: count of reads
            name_sample = self.prefix_pe[os.path.basename(output_files[0]).split('_concat_paired_R1.fq')[0]]            
            stat= Monitoring(self.tot_pe[name_sample])
            ## non ncRNA reads
            non_ncRNA_reads = stat.count_seq_fq(output_files[0])+stat.count_seq_fq(output_files[2])
            self.prt_progress(
                name_sample, 'remove ncRNA', 'SortMeRNA', 'filtered reads (1st)', 
                str(non_ncRNA_reads), stat.get_tot_percentage(non_ncRNA_reads))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def save_processed_reads(input_file, output_file, logger, logging_mutex):
            """ Copy the processed reads in the output directory. """
            for i in range(3):
                cmd = "cp %s %s " %(input_file[i], output_file[i])
                with logging_mutex:
                    logger.info("Copy the processed reads %s in the output directory"%(input_file[i]))
                    logger.debug("save_processed_reads: cmdline\n" + cmd)
                extern.run(cmd)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 5.5: map_to_reference
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """ This stage ... 
            ... """
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def symlink_metaG(input_file, soft_link_name, logger, logging_mutex):
            """ Make soft link in working directory. """
            with logging_mutex:
                logger.info("Linking files %(input_file)s -> %(soft_link_name)s"%locals())
            self.re_symlink(input_file, soft_link_name, logger, logging_mutex)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def symlink_metaG_index(input_file, soft_link_name, logger, logging_mutex):
            """ Make soft link in working directory. Intended to be used when bwa indexes are present. """
            with logging_mutex:
                logger.info("Linking files %(input_file)s -> %(soft_link_name)s"%locals())
            self.re_symlink(input_file, soft_link_name, logger, logging_mutex)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Map separately paired-end and singletons with 'BamM' and merge the results in one .bam file
        # WARNINGS
        #1. .bam files generated with 'BamM' only contain the mapped reads -> 
        #       be carful with the interpretation of samtools flagstat
        #2. only one alignment per read is kept: the secondary and supplementary are removed
        def map2ref(input_files, output_file, bams, flagstat, logger, logging_mutex):
            """ BamM make. Map all metatranscriptomics reads against metagenomics contigs. """
            if self.has_index(input_files[1], ['.amb','.bwt','.ann','.pac','.sa']):
                # index already exists -> bamm Kept
                index_exists = 'K'
            else:
                # index doesn't exist yet -> bamm keep        
                index_exists = 'k'
            cmd = "bamm make -d %s -c %s %s -s %s -o %s --threads %d -%s --quiet"%(
                    input_files[1],
                    input_files[0][0],
                    input_files[0][1],
                    input_files[0][2],
                    self.args.working_dir,
                    self.args.threads,
                    index_exists) 
            with logging_mutex:     
                logger.info("Map reads [%s] to the reference metagenome %s"%(",".join(input_files[0]), input_files[1]))
                logger.debug("map2ref: cmdline\n"+ cmd)  
            extern.run(cmd)
            cmd2 = "samtools merge -f %s %s %s ; samtools view -b -F2304 %s > %s " % \
                    (bams[2], bams[0], bams[1], bams[2], output_file)
            with logging_mutex:     
                logger.info("Concatenate %s and %s"%(bams[0],bams[1]))
                logger.debug("map2ref: cmdline\n"+ cmd2)  
            extern.run(cmd2)
            cmd3 = "samtools flagstat %s > %s " %(output_file,flagstat)
            with logging_mutex:     
                logger.info("Compute statistics of %(output_file)s (samtools flastat)"%locals())
                logger.debug("map2ref: cmdline\n"+ cmd3)  
            extern.run(cmd3)
            #  ~~~~ monitoring: count of reads  ~~~~ #   
            name_sample = self.prefix_pe[os.path.basename(output_file).split('.bam')[0]]            
            stat= Monitoring(self.tot_pe[name_sample])
            ## reads filtered : mapped with high stringency
            mapped_reads = stat.count_mapping_reads(flagstat,True)
            self.prt_progress(
                name_sample, 'alignment', 'BamM make', 'filtered reads (2nd)', 
                str(mapped_reads), stat.get_tot_percentage(mapped_reads))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def mapping_filter(input_file, output_file, flagstat, logger, logging_mutex):
            """ BamM filter. Select reads which are mapped with high stringency. """
            if self.args.no_mapping_filter:  
                pass 
            else:
                cmd = "bamm filter --bamfile %s --percentage_id %f --percentage_aln %f -o %s "%(
                    input_file,
                    self.args.percentage_id,
                    self.args.percentage_aln,
                    self.args.working_dir)
                with logging_mutex:     
                    logger.info("Filter %(input_file)s" % locals())
                    logger.debug("mapping_filter: cmdline\n"+ cmd)  
                extern.run(cmd)

                cmd3 = "samtools flagstat %s > %s " %(output_file, flagstat)
                with logging_mutex:     
                    logger.info("Compute statistics of %(input_file)s (samtools flastat)" % locals())
                    logger.debug("mapping_filter: cmdline\n"+ cmd3)  
                extern.run(cmd3)

                #  ~~~~ monitoring: count of reads  ~~~~ #   
                name_sample = self.prefix_pe[os.path.basename(output_file).split('_filtered.bam')[0]]            
                stat= Monitoring(self.tot_pe[name_sample])
                ## reads filtered : mapped with high stringency
                mapped_reads_f = stat.count_mapping_reads(flagstat,False)
                self.prt_progress(
                    name_sample, '.bam filter', 'BamM filter', 'mapped reads', 
                    str(mapped_reads_f), stat.get_tot_percentage(mapped_reads_f))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 6: bam_stats
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """ This stage ... 
            ... also has logging processes that probably belong elsewhere ..."""
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def bam2normalized_cov(input_file, output_file, coverage_file, dir_bins, logger, logging_mutex):
            """ Dirseq (compute coverage values) +  coverage2normalized_cov. """
            #                                                                       #
            ## add control! contigs in gff files must be present in metaG_contigs  ##
            #                   
            lib_size = int(subprocess.check_output("samtools view -c "+input_file, shell=True))
            #f = open(lib_size_log, 'w')
            #f.write(str(lib_size)) 
            #f.close()                                                                #
            #
            for i in range(len(self.list_gff)):
                cmd = "dirseq --bam %s --gff %s --ignore-directions -q>  %s "% \
                        (input_file, self.list_gff[i], coverage_file)
                with logging_mutex:     
                    logger.info("Calculate coverage from %s and %s"%(input_file,self.list_gff[i]))  
                    logger.debug("bam2normalized_cov: cmdline\n"+ cmd)                                       
                extern.run(cmd)        
                
                if lib_size != 0:
#                    cmd1 = "sed 's/\t/|/g' %s | awk  -F '|' 'NR>=2 {$6= $6/%d*10e6}1' OFS='|' |  " + \
#                            "sed 's/|/\t/g' > %s ; rm %s "%( \
                    cmd1 = "sed 's/\t/|/g' %s | awk  -F '|' 'NR>=2 {$6= $6/%d*10e6}1' OFS='|' |  sed 's/|/\t/g' > %s ; rm %s " %( \
                        coverage_file,
                        lib_size,
                        input_file.split('.bam')[0]+'_'+ 
                            os.path.splitext(os.path.basename((self.list_gff[i])))[0]+'_normalized_cov.csv',
                        coverage_file)
                else:
                    cmd1 = "cp %s %s; rm %s "%(
                        coverage_file,
                        input_file.split('.bam')[0]+'_'+ 
                            os.path.splitext(os.path.basename((self.list_gff[i])))[0]+'_normalized_cov.csv',
                        coverage_file)
                    
                with logging_mutex:     
                    logger.info("Convert coverage to normalized_cov")
                    logger.debug("bam2normalized_cov: cmdline\n" + cmd1)
                extern.run(cmd1)                
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def bam2raw_count(input_file, output_file, list_gff, logger, logging_mutex): 
            """ Bedtools (count the number of mapped reads per gene). """
            for i in range(len(self.list_gff)):
                gff_no_fasta= tempfile.NamedTemporaryFile(prefix='transcriptm', suffix='.gff')
                cmd0 = "sed '/^##FASTA$/,$d' %s > %s" %(self.list_gff[i], gff_no_fasta.name)
                extern.run(cmd0)
                cmd = "bedtools intersect -c -a %s -b %s -bed >  %s " % \
                    (gff_no_fasta.name,
                    input_file,
                    input_file.split('.bam')[0]+'_'+
                        os.path.splitext(os.path.basename((self.list_gff[i])))[0] +'_count.csv')
                with logging_mutex:     
                    logger.info("Calculate raw count from %s and %s "%(input_file, gff_no_fasta.name))  
                    logger.debug("bam2raw_count: cmdline\n"+ cmd)                                       
                extern.run(cmd)        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 7: summary_tables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """ This stage ... creates summary tables in the output directory ..."""
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def transcriptM_table(input_files, output_file, logger, logging_mutex):
            """ Create one table that contains RPKM values for each gene of each bin for the different samples. """
            # Concatenate all the normalized_cov results in a table
            input_files = list(set(input_files))      
            if len(input_files) == 0: 
                raise Exception(
                    "Incorrect input detected. Likely causes: "\
                    "\n\tOnly one sequence file sumitted\n\tAssembly file has been tampered with")
            #normalized_cov_col = [list([]) for _ in xrange(int(len(self.args.paired_end)/2)+3)]  
            normalized_cov_col = [list([]) for _ in range(int(len(self.prefix_pe))+3)]
            # headers of cols ->  0, n-1, n
            normalized_cov_col[0].append('bin_ID')
            normalized_cov_col[-2].append('gene location [contig:start:end]')
            normalized_cov_col[-1].append('annotation')      
        
            bins_path = [os.path.splitext((self.list_gff[i]))[0] for i in range(len(self.list_gff))]
            for b in bins_path :
                files_b= [f for f in input_files if re.search('_'+os.path.basename(b)+'_normalized_cov.csv', f)]  
                # first col: bins_name
                with open(files_b[0],'r') as csvfile:                
                        reader = csv.reader(csvfile, delimiter='\t') 
                        next(reader)  # skip header 
                        for row in reader:
                            normalized_cov_col[0].append(b+'.gff')
                        csvfile.close() 
                # n-1 col: gene location
                with open(files_b[0],'r') as csvfile:                    
                        reader = csv.reader(csvfile, delimiter='\t') 
                        next(reader) # skip header 
                        for row in reader:
                            normalized_cov_col[-2].append(row[0]+':'+row[2]+':'+row[3])
                        csvfile.close() 
                # n col: annotation
                with open(files_b[0],'r') as csvfile:                    
                        reader = csv.reader(csvfile, delimiter='\t') 
                        next(reader) # skip header 
                        for row in reader:
                            normalized_cov_col[-1].append(row[6])
                        csvfile.close() 
                # remaining cols: normalized_cov
                for i in range(len(files_b)):
                    # create header of normalized_cov cols
                    if not normalized_cov_col[i+1]:
                        normalized_cov_col[i+1].append( \
                            'normalized_cov_'+self.prefix_pe[os.path.basename(files_b[i]).split('_')[0]])
                    with open(files_b[i],'r') as csvfile:                    
                        reader = csv.reader(csvfile, delimiter='\t')                      
                        next(reader) # skip header  
                        for row in reader:
                            normalized_cov_col[i+1].append(row[5])
                        csvfile.close() 
         
            tab = numpy.array(normalized_cov_col)
            numpy.savetxt(output_file, numpy.transpose(tab), delimiter='\t', fmt="%s") 
            
            with logging_mutex:     
                logger.info("Create table that contains normalized_cov values for each gene of each bin given " + \
                    "as input for the different samples: %s" %(','.join(self.prefix_pe.values())))    
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def raw_count_table(input_files, output_file, logger, logging_mutex):
            """ Create one table that contains raw count values for each gene of each bin for the different samples.
            """
            # Concatenate all the raw count in a table
            input_files = list(set(input_files))          
            count_col = [list([]) for _ in xrange(int(len(self.args.paired_end)/2)+3)]       
            # headers of cols ->  0, n-1, n
            count_col[0].append('bin_ID')
            count_col[-2].append('gene location [contig:start:end]')
            count_col[-1].append('annotation')      
        
            bins_path =[os.path.splitext((self.list_gff[i]))[0] for i in range(len(self.list_gff))]
            for b in bins_path :
                files_b = [f for f in input_files if re.search('_'+os.path.basename(b)+'_count.csv', f)]  
                # first col: bins_name
                with open(files_b[0], 'r') as csvfile:                    
                        reader = csv.reader(csvfile, delimiter='\t') 
                        for row in reader:
                            count_col[0].append(b+'.gff')
                        csvfile.close() 
                # n-1 col: gene location
                with open(files_b[0], 'r') as csvfile:                    
                        reader = csv.reader(csvfile, delimiter='\t') 
                        for row in reader:
                            count_col[-2].append(row[0]+':'+row[3]+':'+row[4])
                        csvfile.close() 
                # n col: annotation
                with open(files_b[0], 'r') as csvfile:                    
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
            numpy.savetxt(output_file, numpy.transpose(tab), delimiter='\t', fmt="%s") 
            
            with logging_mutex:     
                logger.info("Create table that contains raw count values for each gene of each bin given as input " + \
                    "for the different samples: %s"%(','.join(self.prefix_pe.values())))    
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage Reporting (which probably doesn't belong in Stage 6).
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def save_log(input_files, output_files, logger, logging_mutex):
            """ Save the log files, generated for different stages of the pipeline (in the temp directory). """
            cmd = "cp %s %s"   %(input_files, output_files)    
            with logging_mutex:
                logger.info("Save log files: %(input_files)s" % locals())
                logger.debug("save_log: cmdline\n"+cmd)                
            extern.run(cmd)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def logtable(input_files, output_file, basename):
            """ Sums up the count of reads which are kept after each step. """
            name_sample = self.prefix_pe[basename]
            stat = Monitoring(self.tot_pe[name_sample])          
            # raw
            stat.reads[0]= stat.count_raw_reads()
            # processed reads
            trimm_file = [f for f in input_files if re.search(r'trimmomatic.log', f)][0]     
            stat.reads[1]= stat.count_processed_reads(trimm_file)         
            # non phiX reads
            phix_ID_file = [f for f in input_files if re.search(r'phiX_ID.log', f)][0]
            stat.reads[2] = stat.count_non_Phix_reads(phix_ID_file) 
            # non rRNA/tRNA/tmRNA reads
            list_fqfiles = self.get_files(self.args.working_dir,'.fq')
            pairs_filtered = [f for f in list_fqfiles if re.search(r'%s.+concat_paired_R1.fq'%(basename), f)][0]
            singles_filtered = [f for f in list_fqfiles if re.search(r'%s.+concat_single.fq'%(basename), f)][0]
            stat.reads[3] = stat.count_non_ncRNA_reads(pairs_filtered,singles_filtered)    
            # mapped reads
            mapping_log = [f for f in input_files if re.search(r'mapping.log', f)][0]  
            stat.reads[4] = stat.count_mapping_reads(mapping_log,True)      

            # reads mapped with a given stringency
            if self.args.no_mapping_filter :
                # save stat_table
                tab = numpy.array([[name_sample]*5,
                      ["raw data", "trimming","remove PhiX","remove ncRNA","alignment "],
                      ["FastQC-check", "Trimmomatic","bamM make","SortMeRNA","bamM make"],
                      ["raw reads", "raw reads","processed reads","filtered reads","filtered reads"],
                      map(str,stat.reads[:-1]),
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
                      map(str,stat.reads),
                      map(str,stat.get_all_tot_percentage()),
                      map(str,stat.get_all_percentage_prev())])
                numpy.savetxt(output_file, numpy.transpose(tab), delimiter='\t', fmt="%s") 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def concatenate_logtables(input_files, output_file, logger, logging_mutex):
            """ Concatenate the summaries of reads distribution from all samples. """
            input_files = (' ').join(input_files)
            h = ["sample name","step name","tool used","input data","reads count","% total","% previous step"]  
            header = ('\t').join(h)            
            cmd = "cat %s > %s ;sed -i '1i%s' %s" %(input_files, output_file, header, output_file)    
            with logging_mutex:
                logger.info("Concatenate summaries: %(input_files)s" % locals())
                logger.debug("concatenate_logtables: cmdline\n"+cmd)                
            extern.run(cmd)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # abeyance: saving information at the end of each stage, for halting and restarting of pipeline runs.
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        class abeyance(object):
            """ for disk storage and retrieval of output files from tasks in Ruffus pipeline. """
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            @staticmethod
            def store(input_filenames, output_filename, logger, logging_mutex):
                """ saves filenames (various structures of strings) to disk in JSON format; does a log entry. """
                with open(output_filename, 'w') as outfile:
                    json.dump(input_filenames, outfile)
                with logging_mutex:     
                    logger.info("Saved abeyance output filenames to " + output_filename)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            @staticmethod
            def retrieve(input_file, logger, logging_mutex):
                """ retrieves filenames (various structures of strings) from disk in JSON format; does a log entry. """
                try:
                    infile = open(input_file, 'r')
                except:
                    raise Exception(
                        "\nAttempt to open file <{0}> failed when attempting a restart at an intermittent stage. "\
                        "Ensure that the pipeline has already been run up to this point.".format(input_file))
                else:
                    with open(input_file, 'r') as infile:
                        output_filenames = json.load(infile)
                    with logging_mutex:     
                        logger.info("Retrieved abeyance output filenames from " + input_file)
                    return output_filenames
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            """ wrappers to provide unqualified names that are unique, as required by ruffus pipeline tasks. """
            # abeyance storage:
            @staticmethod
            def store_symlink_metaT(p1,p2,p3,p4):
                abeyance.store(p1,p2,p3,p4) # See for meaningful description.
            @staticmethod
            def store_trimmomatic(p1,p2,p3,p4):
                abeyance.store(p1,p2,p3,p4) # See for meaningful description.
            @staticmethod
            def store_phiX_extract(p1,p2,p3,p4):
                abeyance.store(p1,p2,p3,p4) # See for meaningful description.
            @staticmethod
            def store_sortmerna(p1,p2,p3,p4):
                abeyance.store(p1,p2,p3,p4) # See for meaningful description.
            @staticmethod
            def store_concat_for_mapping(p1,p2,p3,p4):
                abeyance.store(p1,p2,p3,p4) # See for meaningful description.
            @staticmethod
            def store_map2ref(p1,p2,p3,p4):
                abeyance.store(p1,p2,p3,p4) # See for meaningful description.
            @staticmethod
            def store_mapping_filter(p1,p2,p3,p4):
                abeyance.store(p1,p2,p3,p4) # See for meaningful description.
            @staticmethod
            def store_bam2normalized_cov(p1,p2,p3,p4):
                abeyance.store(p1,p2,p3,p4) # See for meaningful description.
            @staticmethod
            def store_bam2raw_count(p1,p2,p3,p4):
                abeyance.store(p1,p2,p3,p4) # See for meaningful description.
            @staticmethod
            def store_transcriptM_table(p1,p2,p3,p4):
                abeyance.store(p1,p2,p3,p4) # See for meaningful description.
            @staticmethod
            def store_raw_count_table(p1,p2,p3,p4):
                abeyance.store(p1,p2,p3,p4) # See for meaningful description.
            # abeyance retrieval:
            @staticmethod
            def retrieve_symlink_metaT(p1,p2,p3):
                return abeyance.retrieve(p1,p2,p3) # See for meaningful description.
            @staticmethod
            def retrieve_trimmomatic(p1,p2,p3):
                return abeyance.retrieve(p1,p2,p3) # See for meaningful description.
            @staticmethod
            def retrieve_phiX_extract(p1,p2,p3):
                return abeyance.retrieve(p1,p2,p3) # See for meaningful description.
            @staticmethod
            def retrieve_sortmerna(p1,p2,p3):
                return abeyance.retrieve(p1,p2,p3) # See for meaningful description.
            @staticmethod
            def retrieve_concat_for_mapping(p1,p2,p3):
                return abeyance.retrieve(p1,p2,p3) # See for meaningful description.
            @staticmethod
            def retrieve_map2ref(p1,p2,p3):
                return abeyance.retrieve(p1,p2,p3) # See for meaningful description.
            @staticmethod
            def retrieve_mapping_filter(p1,p2,p3):
                return abeyance.retrieve(p1,p2,p3) # See for meaningful description.
            @staticmethod
            def retrieve_bam2normalized_cov(p1,p2,p3):
                return abeyance.retrieve(p1,p2,p3) # See for meaningful description.
            @staticmethod
            def retrieve_bam2raw_count(p1,p2,p3):
                return abeyance.retrieve(p1,p2,p3) # See for meaningful description.
            @staticmethod
            def retrieve_transcriptM_table(p1,p2,p3):
                return abeyance.retrieve(p1,p2,p3) # See for meaningful description.
            @staticmethod
            def retrieve_raw_count_table(p1,p2,p3):
                return abeyance.retrieve(p1,p2,p3) # See for meaningful description.
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Build the pipeline.
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rpl = ruffus.Pipeline.pipelines["main"] # rpl: 'Ruffus PipeLine'
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 1
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.require_run_stage(1):
            rpl.originate(task_func = symlink_metaT,
                output = self.alias_pe.keys(), # soft-link filenames of metatranscriptomic paired-end reads.
                extras = [self.logger, self.logging_mutex]
                )\
            .follows(calc_tot_pe)\
            .mkdir(self.args.working_dir)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.transform(task_func = view_raw_data, # Create the first output: a fastqc report of raw DATA
                input = symlink_metaT, 
                filter = ruffus.formatter(),
                output = os.path.join(self.SUBDIR_FastQC_raw, "{basename[0]}"+"_fastqc.zip"),
                extras = [self.logger, self.logging_mutex]
                )\
            .mkdir(self.SUBDIR_FastQC_raw)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.merge(task_func = abeyance.store_symlink_metaT, 
                input = symlink_metaT, 
                output = os.path.join(self.args.working_dir, 'abeyance_symlink_metaT.output'),
                extras = [self.logger, self.logging_mutex]
                )\
            .follows(view_raw_data)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 2
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.require_stage_restart(2):
            decided_symlink_metaT = abeyance.retrieve_symlink_metaT(
                os.path.join(self.args.working_dir, 'abeyance_symlink_metaT.output'), self.logger, self.logging_mutex)
        else:
            decided_symlink_metaT = symlink_metaT
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.require_run_stage(2):
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.collate(task_func = trimmomatic,
                input = decided_symlink_metaT,
                filter = ruffus.regex("R[12].fq.gz$"),
                output = ["trimm_P1.fq.gz", "trimm_P2.fq.gz", "trimm_U1.fq.gz", "trimm_U2.fq.gz"],
                extras = [
                    "trimmomatic.log",
                    self.logger, self.logging_mutex]
                )\
            .follows(calc_tot_pe)\
            .mkdir(self.args.working_dir)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.transform(task_func = view_processed_data,
                input = trimmomatic, 
                filter = ruffus.formatter(),
                output = os.path.join(self.SUBDIR_FastQC_processed, "{basename[0]}"+"_fastqc.zip"),
                extras = [self.logger, self.logging_mutex]
                )\
            .mkdir(self.SUBDIR_FastQC_processed)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.merge(task_func = abeyance.store_trimmomatic,
                input = trimmomatic, 
                output = os.path.join(self.args.working_dir, 'abeyance_trimmomatic.output'),
                extras = [self.logger, self.logging_mutex]
                )\
            .follows(view_processed_data)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 3
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.require_stage_restart(3):
            decided_trimmomatic = abeyance.retrieve_trimmomatic(
                os.path.join(self.args.working_dir, 'abeyance_trimmomatic.output'), self.logger, self.logging_mutex)
        else:
            decided_trimmomatic = trimmomatic
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.require_run_stage(3):
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.subdivide(task_func = phiX_map,
                input = decided_trimmomatic,
                filter = ruffus.formatter(r"(.+)/(?P<BASE>.*)P1.fq.gz"),
                output = [
                    "{path[0]}/phiX.{BASE[0]}P1.bam",
                    "{path[0]}/phiX.{BASE[0]}U1.bam",
                    "{path[0]}/phiX.{BASE[0]}U2.bam"], 
                extras = [self.logger, self.logging_mutex]
                )
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.transform(task_func = phiX_ID,
                input = phiX_map,
                filter = ruffus.suffix(".bam"),
                output = ".txt",
                extras = [self.logger, self.logging_mutex]
                )
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.collate(task_func = phiX_concat_ID,
                input = phiX_ID, 
                filter = ruffus.formatter(r"phiX.(?P<BASE>.*)[UP][12].txt$"),
                output = '{path[0]}/{BASE[0]}phiX_ID.log',
                extras = [
                    '{BASE[0]}',
                    self.logger, self.logging_mutex]
                )\
            .follows(calc_tot_pe)\
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.subdivide(task_func = QC_output, 
                input = decided_trimmomatic,
                filter = ruffus.regex(r"trimm_[UP][12].fq.gz"),
                output = ["trimm_P1.fq.gz", "trimm_P2.fq.gz", "trimm_U1.fq.gz", "trimm_U2.fq.gz"]
                )
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.transform(task_func = phiX_extract, 
                input = QC_output,
                filter = ruffus.suffix(".fq.gz"), 
                add_inputs = ruffus.add_inputs(phiX_concat_ID), 
                output = "_phiX_ext.fq", 
                extras = [self.logger, self.logging_mutex]
                )
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.merge(task_func = abeyance.store_phiX_extract, 
                input = phiX_extract, 
                output = os.path.join(self.args.working_dir, 'abeyance_phiX_extract.output'),
                extras = [self.logger, self.logging_mutex]
                )
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 4
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.require_stage_restart(4):
            decided_phiX_extract = abeyance.retrieve_phiX_extract(
                os.path.join(self.args.working_dir, 'abeyance_phiX_extract.output'), self.logger, self.logging_mutex)
        else:
            decided_phiX_extract = phiX_extract
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.require_run_stage(4):
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.subdivide(task_func = sortmerna, 
                input = decided_phiX_extract,
                filter = ruffus.formatter(), 
                output = "{path[0]}/{basename[0]}_non_ncRNA.fq",
                extras = [
                    "{path[0]}/{basename[0]}_ncRNA.fq",
                    self.logger, self.logging_mutex]
                )
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.merge(task_func = abeyance.store_sortmerna, 
                input = sortmerna, 
                output = os.path.join(self.args.working_dir, 'abeyance_sortmerna.output'),
                extras = [self.logger, self.logging_mutex]
                )
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 5.2
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.require_stage_restart(5.2):
            decided_sortmerna = abeyance.retrieve_sortmerna(
                os.path.join(self.args.working_dir, 'abeyance_sortmerna.output'), self.logger, self.logging_mutex)
        else:
            decided_sortmerna = sortmerna
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.require_run_stage(5.2):
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.collate(task_func = concat_for_mapping, 
                input = decided_sortmerna,
                filter = ruffus.regex(r"trimm_.*"),
                output = ["concat_paired_R1.fq", "concat_paired_R2.fq", "concat_single.fq"],
                extras = [
                    "ID_single.txt", 
                    self.logger, self.logging_mutex]
                )\
            .follows(calc_tot_pe)\
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.transform(task_func = save_processed_reads, 
                input = concat_for_mapping, 
                filter = ruffus.formatter(),
                output = [
                    os.path.join(self.SUBDIR_processed_reads, "{basename[0]}{ext[0]}"),
                    os.path.join(self.SUBDIR_processed_reads, "{basename[1]}{ext[1]}"),
                    os.path.join(self.SUBDIR_processed_reads, "{basename[2]}{ext[2]}")],
                extras = [self.logger, self.logging_mutex]
                )\
            .mkdir(self.SUBDIR_processed_reads)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.merge(task_func = abeyance.store_concat_for_mapping, 
                input = concat_for_mapping, 
                output = os.path.join(self.args.working_dir, 'abeyance_concat_for_mapping.output'),
                extras = [self.logger, self.logging_mutex]
                )\
            .follows(save_processed_reads)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 5.5
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.require_stage_restart(5.5):
            decided_concat_for_mapping = abeyance.retrieve_concat_for_mapping(
                os.path.join(self.args.working_dir, 'abeyance_concat_for_mapping.output'), 
                self.logger, self.logging_mutex)
        else:
            decided_concat_for_mapping = concat_for_mapping
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.require_run_stage(5.5):
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.transform(task_func = symlink_metaG, 
                input = self.args.metaG_contigs, # filename of all contigs from the ref metagenome (in a fasta file).
                filter = ruffus.formatter(), 
                output = os.path.join(self.args.working_dir, "{basename[0]}"+".fa"), # put in working directory.
                extras = [self.logger, self.logging_mutex]
                )\
            .mkdir(self.args.working_dir)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.transform(task_func = symlink_metaG_index, 
                input = [self.args.metaG_contigs+x for x in ['.amb','.bwt','.ann','.pac','.sa']], 
                filter = ruffus.formatter(),
                output = os.path.join(self.args.working_dir,"{basename[0]}{ext[0]}"), # put in working directory.
                extras = [self.logger, self.logging_mutex]
                )\
            .mkdir(self.args.working_dir)\
            .active_if(self.has_index(self.args.metaG_contigs, ['.amb','.bwt','.ann','.pac','.sa']))
                # Intended to be used when bwa indexes are present.
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.transform(task_func = map2ref, 
                input = decided_concat_for_mapping, 
                filter = ruffus.formatter(r"(.+)/(?P<BASE>.*)_concat_paired_R1.fq"), 
                add_inputs = ruffus.add_inputs(symlink_metaG),
                output = "{path[0]}/{BASE[0]}.bam",
                extras = [
                ["{path[0]}/"+os.path.splitext(os.path.basename(self.args.metaG_contigs))[0]+".{basename[0]}.bam",
                 "{path[0]}/"+os.path.splitext(os.path.basename(self.args.metaG_contigs))[0]+".{basename[2]}.bam",
                 "{path[0]}/"+os.path.splitext(os.path.basename(self.args.metaG_contigs))[0]+"{BASE[0]}_merged.bam"],
                "{path[0]}/{BASE[0]}_mapping.log",
                self.logger, self.logging_mutex]
                )\
            .follows(symlink_metaG_index)\
            .follows(calc_tot_pe)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.merge(task_func = abeyance.store_map2ref, 
                input = map2ref, 
                output = os.path.join(self.args.working_dir, 'abeyance_map2ref.output'),
                extras = [self.logger, self.logging_mutex]
                )
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.transform(task_func = mapping_filter, 
                input = map2ref,
                filter = ruffus.formatter('.bam'),
                output = "{path[0]}/{basename[0]}_filtered.bam", 
                extras = [
                    "{path[0]}/{basename[0]}_stringency_filter.log", 
                    self.logger, self.logging_mutex]
                )\
            .follows(calc_tot_pe)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.merge(task_func = abeyance.store_mapping_filter, 
                input = mapping_filter, 
                output = os.path.join(self.args.working_dir, 'abeyance_mapping_filter.output'),
                extras = [self.logger, self.logging_mutex]
                )
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 6
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """ NB. flag `args.no_mapping_filter`:
            (a) determines whether mapping_filter() is performed (to transform the results of map2ref()).
                At the moment, mapping_filter() is still entered into, but the flag is checked again inside 
                the function so that nothing is performed.
            (b) determines whether the inputs for bam2raw_count() and bam2normalized_cov() are taken directly 
                from map2ref() or whether they are taken from the subsequent mapping_filter().
        """
        if self.require_stage_restart(6):
            if self.args.no_mapping_filter:  
                decided_bam_file = abeyance.retrieve_map2ref(
                    os.path.join(self.args.working_dir, 'abeyance_map2ref.output'), 
                    self.logger, self.logging_mutex)
            else: 
                decided_bam_file = abeyance.retrieve_mapping_filter(
                    os.path.join(self.args.working_dir, 'abeyance_mapping_filter.output'), 
                    self.logger, self.logging_mutex)
        else:
            if self.args.no_mapping_filter:  
                decided_bam_file = map2ref 
            else: 
                decided_bam_file = mapping_filter
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.require_run_stage(6):
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.subdivide(task_func = bam2normalized_cov, 
                input = decided_bam_file, 
                filter = ruffus.formatter(),
                output = '{path[0]}/*normalized_cov.csv',
                extras = [
                    '{path[0]}/coverage.csv',
                    self.args.dir_bins,
                    self.logger, self.logging_mutex]
                )
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.merge(task_func = abeyance.store_bam2normalized_cov, 
                input = bam2normalized_cov, 
                output = os.path.join(self.args.working_dir, 'abeyance_bam2normalized_cov.output'),
                extras = [self.logger, self.logging_mutex]
                )
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.subdivide(task_func = bam2raw_count, 
                input = decided_bam_file,
                filter = ruffus.formatter(),
                output = '{path[0]}/*count.csv',
                extras = [
                    self.list_gff,
                    self.logger, self.logging_mutex]
                )
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.merge(task_func = abeyance.store_bam2raw_count, 
                input = bam2raw_count, 
                output = os.path.join(self.args.working_dir, 'abeyance_bam2raw_count.output'),
                extras = [self.logger, self.logging_mutex]
                )
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #! Stage 6 reporting - NEED TO CHECK OUT WHERE THESE SHOULD GO.
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.transform(task_func = save_log, 
                input = self.args.working_dir+'/*.log', 
                filter = ruffus.formatter(".log"),  
                output = os.path.join(self.SUBDIR_log, "{basename[0]}"+".log"),
                extras = [self.logger, self.logging_mutex]
                )\
            .mkdir(self.SUBDIR_log)\
            .follows(bam2normalized_cov)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.collate(task_func = logtable, 
                input = save_log, 
                filter = 
                    ruffus.formatter(r"/log/(?P<BASE>.*)_((stringency_filter)|(mapping)|(trimmomatic)|" + \
                        "(trimm_((phiX_ID)|((U|P)(1|2)_phiX_ext_ncRNA)))).log$"),
                output = self.SUBDIR_reads_distribution + "/{BASE[0]}_reads_stat",
                extras = ['{BASE[0]}']
                )\
            .mkdir(self.SUBDIR_reads_distribution)\
            .follows(calc_tot_pe)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.merge(task_func = concatenate_logtables, 
                input = logtable, 
                output = os.path.join(self.args.output_dir, 'summary_reads'), 
                extras = [self.logger, self.logging_mutex]
                )
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Stage 7
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.require_stage_restart(7):
            decided_bam2normalized_cov = abeyance.retrieve_bam2normalized_cov(
                os.path.join(self.args.working_dir, 'abeyance_bam2normalized_cov.output'), 
                self.logger, self.logging_mutex)
            decided_bam2raw_count = abeyance.retrieve_bam2raw_count(
                os.path.join(self.args.working_dir, 'abeyance_bam2raw_count.output'), 
                self.logger, self.logging_mutex)
        else:
            decided_bam2normalized_cov = bam2normalized_cov
            decided_bam2raw_count = bam2raw_count
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.require_run_stage(7):
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.merge(task_func = transcriptM_table, # Concatenate all the normalized_cov results in a table
                input = decided_bam2normalized_cov, 
                output = os.path.join(self.args.output_dir, 
                                      os.path.basename(self.args.output_dir)+'_NORM_COVERAGE.csv'),
                extras = [self.logger, self.logging_mutex]
                )\
            .mkdir(self.args.output_dir)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.merge(task_func = abeyance.store_transcriptM_table, 
                input = transcriptM_table, 
                output = os.path.join(self.args.working_dir, 'abeyance_transcriptM_table.output'),
                extras = [self.logger, self.logging_mutex]
                )
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.merge(task_func = raw_count_table, # Concatenate all the raw count in a table
                input = decided_bam2raw_count,
                output = os.path.join(self.args.output_dir, os.path.basename(self.args.output_dir)+'_COUNT.csv'),
                extras = [self.logger, self.logging_mutex]
                )\
            .mkdir(self.args.output_dir)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            rpl.merge(task_func = abeyance.store_raw_count_table, 
                input = raw_count_table, 
                output = os.path.join(self.args.working_dir, 'abeyance_raw_count_table.output'),
                extras = [self.logger, self.logging_mutex]
                )
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Assign the ruffus pipeline to the object.
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.built_pipeline = rpl
        print "Pipeline built."
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # PIPELINE: RUN
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #for task in ruffus.pipeline_get_task_names():
        #    print task
        #ruffus.cmdline.run(self.args, target_tasks = [trimmomatic])
        #ruffus.cmdline.run(target_tasks = [trimmomatic])
        #rpl.run(target_tasks = [trimmomatic])
        #rpl.run(target_tasks = view_processed_data,  verbose = 1) # verbose = 0
        #rpl.run(target_tasks = [self.args.end_stages[0]], verbose = 0) # verbose = 0 defaults to 1
        #rpl.run(target_tasks = [symlink_metaG], verbose = 0) # verbose = 0 defaults to 1
        #multithread = 3
        #multiprocess = 5)
        
        #ruffus.cmdline.run(self.args, target_tasks = self.target_tasks)
#        if False and len(self.args.end_stages) == 0:
#            rpl.run(target_tasks = [self.args.end_stages[0]], verbose = 0) # verbose = 0 defaults to 1
#        else:
#            #rpl.run(target_tasks = [], verbose = 0) # verbose = 0 defaults to 1
#            rpl.run(target_tasks = self.args.end_stages, verbose = 0) # verbose = 0 defaults to 1
        #rpl.run(target_tasks = self.target_tasks, verbose = 0) # verbose = 0 defaults to 1

#        print "self.args.verbose:", self.args.verbose, type(self.args.verbose)
#        v = int(self.args.verbose[0])
#        print "v:", v, type(v)
#
#        v = 1 or int(self.args.verbose[0])
        #rpl.run(target_tasks = self.target_tasks, logger = self.logger, verbose = v) # verbose = 0 defaults to 1
#        rpl.run(target_tasks = self.target_tasks, logger = self.logger, verbose = v) # verbose = 0 defaults to 1
        #rpl.run(target_tasks = self.target_tasks, logger = self.logger, verbose = v) # verbose = 0 defaults to 1

        #verbosity = 0 or int(self.args.verbose[0])
#        self.built_pipeline.run(target_tasks = self.args.target_tasks, logger = self.logger, verbose = verbosity)
#        print "debug: about to run cmdline"
#        ruffus.cmdline.run(self.args)
        #self.built_pipeline.run(logger = self.logger, verbose = verbosity)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Run the pipeline
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def run_built_pipeline(self):
        assert self.built_pipeline != None, "The pipeline needs to be built before it can be run."
        print "Pipeline running..."
        ruffus.cmdline.run(self.args)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Post-pipeline activity
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def clear(self):
        """ tidies-up the program's output directory, (renames certain files, and deletes a subdirectory),
            intended to be performed after the pipeline stages have been run.
        """
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def rename_SUBDIRs_content():
            """ renames all files in output directories with their prefixes referring to their original name. """
            for subdir in self.SUBDIRS_for_content_renaming:
                if os.path.exists(subdir):
                    # Iterate through each file in that subdirectory, to rename the file.
                    for f in os.listdir(subdir): 
                        original_name = self.prefix_pe.get(f.split('_')[0])
                        # Allow for the case that the file has already been renamed in previous runs.
                        if original_name != None: 
                            # Replace the first part of the filenames, i.e. up to the first "_",
                            # by the first part of `self.prefix_pe`, i.e. up to the first "_". 
                            f_oldname = os.path.join(subdir, f)          
                            f_newname = os.path.join(subdir, string.replace(f, f.split('_')[0], original_name))
                            os.rename(f_oldname, f_newname)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        def clear_SUBDIRs():
            """ removes all nominated output directories. """
            for subdir in self.SUBDIRS_for_clearing:
                try:
                    shutil.rmtree(subdir)
                except OSError:
                    pass
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        rename_SUBDIRs_content()
        clear_SUBDIRs()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import operator

__author__ = 'Boris Zinshteyn' and 'Anthony Schuller'
"""
Originally intended for processing of DMS or other chemical probing data on rRNA
Now adapted for processing of cleaved rRNA by hydroxyl radical probing
Based on Alex Robertson's original RNA Bind n Seq pipeline, available on github
"""
import sys
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42  #leaves most text as actual text in PDFs, not outlines
import os
import argparse
import subprocess

import hrf_settings
import hrf_utils
import hrf_lib
import hrf_plotting

class hrf_seq_run:
    def __init__(self, settings, threads):
        self.threads = threads
        self.settings = settings
        self.remove_adaptor()
        self.trim_reads()
        self.generate_mapping_index()
        self.map_reads()
       # self.initialize_libs()
       # self.make_plots()
       # self.make_plots(exclude_constitutive=True)
       # self.make_tables()
       # self.make_tables(exclude_constitutive=True)
       # self.annotate_structures()
       # self.annotate_structures(exclude_constitutive=True)

    def remove_adaptor(self):
        if not self.settings.get_property('force_retrim'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.adaptorless_reads_exist():
                    break
            else:
                return

        if self.settings.get_property('trim_adaptor'):
            self.settings.write_to_log( 'trimming adaptors')
            hrf_utils.make_dir(self.rdir_path('adaptor_trimmed_reads'))
            map(lambda lib_setting: self.remove_adaptor_one_lib(lib_setting, self.threads),
                          self.settings.iter_lib_settings())
            self.settings.write_to_log( 'trimming adaptors done')

    def remove_adaptor_one_lib(self, lib_settings, threads):
        #TODO: implement paired end read trimming
        lib_settings.write_to_log('adaptor trimming')
        if self.settings.get_property('paired_end'):
            command_to_run = 'cutadapt --adapter %s --overlap 3 --discard-untrimmed --minimum-length %d %s --output %s 1>>%s 2>>%s' % (self.settings.get_property('adaptor_sequence'), self.settings.get_property('min_post_adaptor_length'),
                               lib_settings.get_fastq_file(), lib_settings.get_adaptor_trimmed_reads(), lib_settings.get_log(),
                               lib_settings.get_log())
        else:
            """
            -x specifies the 3' adaptor to trim from the forward read
            -y specifies the 5' adaptor to trim from the reverse read
            -m, --mode <str> trimming mode; 1) single-end -- head: 5' end; tail: 3' end; any: anywhere (tail)
                           2) paired-end -- pe: paired-end; mp: mate-pair; ap: amplicon (pe)


            -Q specifies the lowest acceptable mean read quality before trimming
            -l specifies the minimum post-trimming read length
            -L specifies the maximum post-trimming read length
            -o is the output prefix
            --threads specifies number of threads to use
            """
            command_to_run = 'skewer -x %s  -l %d -o %s --quiet --threads %s %s 1>>%s 2>>%s' % (
                    #option --compress did not function
                self.settings.get_property('read1_adaptor_sequence'),
                self.settings.get_property('min_post_adaptor_length'),
                lib_settings.get_adaptor_trimmed_reads(prefix_only=True)[0],
                threads,
                lib_settings.get_fastq_gz()[0],
                lib_settings.get_log(), lib_settings.get_log())

            subprocess.Popen(command_to_run, shell=True).wait()

            compression_command = "gzip %s-trimmed.fastq" % (
                lib_settings.get_adaptor_trimmed_reads(prefix_only=True)[0])
            subprocess.Popen(compression_command, shell=True).wait()
        lib_settings.write_to_log('adaptor trimming done')

    def trim_reads(self):
        #TODO: this step might be doable in skewer
        #TODO: this step needs to be paired-end compatible
        #TODO: it may be best to collapse identical mate pairs before removing the random barcode
        """
        Trim reads by given amount, removing potential random barcoding sequences from 5' end
        Trimming from 3' end can also help if mapping is problematic by reducing chance for indels to prevent mapping
        :return:
        """
        self.settings.write_to_log('trimming reads')
        if not self.settings.get_property('force_retrim'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.trimmed_reads_exist():
                    break
            else:
                return
        hrf_utils.make_dir(self.rdir_path('trimmed_reads'))
        hrf_utils.parmap(lambda lib_setting: self.trim_one_lib(lib_setting), self.settings.iter_lib_settings(),
                         nprocs = self.threads)
        self.settings.write_to_log('trimming reads complete')

    def trim_one_lib(self, lib_settings):
        lib_settings.write_to_log('trimming_reads')
        first_base_to_keep = self.settings.get_property('first_base_to_keep')  #the trimmer is 1-indexed. 1 means keep
                                                                               #  every base
        last_base_to_keep = self.settings.get_property('last_base_to_keep')  #Will keep entire 3' end if this is greater
                                                                             #than or equal to the read length

        subprocess.Popen('gunzip -c %s | fastx_trimmer -z -f %d -Q33 -l %d -o %s >>%s 2>>%s' % (lib_settings.get_adaptor_trimmed_reads()[0],
                                                                                      first_base_to_keep, last_base_to_keep,
                                                                                      lib_settings.get_trimmed_reads()[0],
                                                                                      lib_settings.get_log(),
                                                                                      lib_settings.get_log()), shell=True).wait()

        lib_settings.write_to_log('trimming_reads done')

    def generate_mapping_index(self):
        """
        builds a STAR index from the input fasta file
        """
        self.settings.write_to_log('building STAR index')
        if self.settings.get_property('force_index_rebuild') or not self.settings.rRNA_STAR_index_exists():
            hrf_utils.make_dir(self.rdir_path('star_indices'))
            hrf_utils.make_dir(self.settings.get_rRNA_STAR_index())


            subprocess.Popen('STAR --runThreadN %d --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --genomeSAindexNbases 4 1>>%s 2>>%s' %
                             (self.threads, self.settings.get_rRNA_STAR_index(), self.settings.get_rRNA_fasta(), self.settings.get_log(), self.settings.get_log()), shell=True).wait()
        self.settings.write_to_log('building STAR index complete')

    #adjust for paired ends

    def map_reads(self):
        """
        map all reads using STAR
        :return:
        """
        self.settings.write_to_log('mapping reads')
        if not self.settings.get_property('force_remapping'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.mapped_reads_exist():
                    break
            else:
                return
        hrf_utils.make_dir(self.rdir_path('mapped_reads'))
        map(lambda lib_setting: self.map_one_library(lib_setting, self.threads), self.settings.iter_lib_settings())
        self.settings.write_to_log( 'finished mapping reads')

    def map_one_library(self, lib_settings, threads):
        lib_settings.write_to_log('mapping_reads')
        command_to_run = 'STAR --limitBAMsortRAM 2000000000 --runThreadN %d --genomeDir %s --readFilesIn %s --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterMultimapNmax 1 --outWigType wiggle read1_5p --outFileNamePrefix %s --outReadsUnmapped FastX 1>>%s 2>>%s' %\
                         (threads, self.settings.get_rRNA_STAR_index(), lib_settings.get_trimmed_reads()[0],
                          lib_settings.get_mapped_reads_prefix(), lib_settings.get_log(), lib_settings.get_log())

        subprocess.Popen(command_to_run, shell=True).wait()
        #sort transcript-mapped bam file

        #subprocess.Popen('samtools sort %s %s.temp_sorted 1>>%s 2>>%s' % (lib_settings.get_transcript_mapped_reads(), lib_settings.get_transcript_mapped_reads(),
        #                                                                  lib_settings.get_log(), lib_settings.get_log()), shell=True).wait()
        #subprocess.Popen('mv %s.temp_sorted.bam %s' % (lib_settings.get_transcript_mapped_reads(),
        #                                                                  lib_settings.get_transcript_mapped_reads()), shell = True).wait()

        #subprocess.Popen('samtools index %s' % (lib_settings.get_transcript_mapped_reads()), shell = True).wait()
        subprocess.Popen('samtools index %s' % (lib_settings.get_mapped_reads()), shell=True).wait()
        lib_settings.write_to_log('mapping_reads done')


    def initialize_libs(self):
        self.settings.write_to_log('initializing libraries, counting reads')
        self.libs = []
        map(lambda lib_settings: self.initialize_lib(lib_settings), self.settings.iter_lib_settings())
        self.settings.write_to_log('initializing libraries, counting reads, done')

    def initialize_lib(self, lib_settings):
        lib = hrf_lib.ModLib(self, self.settings, lib_settings)
        self.libs.append(lib)

    def get_lib_from_name(self, normalizing_lib_name):
        for lib in self.libs:
            if lib.lib_settings.sample_name == normalizing_lib_name:
                return lib
        return None

    def get_normalizable_libs(self):
        normalizeable_libs = []
        for lib in self.libs:
            if lib.lib_settings.sample_name in self.settings.get_property('experimentals'):
                normalizeable_libs.append(lib)
        return normalizeable_libs

    def get_modified_libs(self):
        modified_libs = []
        for lib in self.libs:
            if (lib.lib_settings.sample_name in
                    self.settings.get_property('experimentals')) or (lib.lib_settings.sample_name in self.settings.get_property('with_mod_controls')):
                modified_libs.append(lib)
        return modified_libs

    def make_tables(self, exclude_constitutive=False):
        subfolders = ['raw', 'background_subtracted', 'control_subtracted', 'fold_change']
        for subfolder in subfolders:
            hrf_utils.make_dir(self.rdir_path('tables', subfolder))
            hrf_utils.make_dir(self.rdir_path('pickles', subfolder))
            hrf_utils.make_dir(self.rdir_path('tables', subfolder, 'exclude_constitutive'))
            hrf_utils.make_dir(self.rdir_path('pickles', subfolder, 'exclude_constitutive'))
        self.pickle_mutation_rates('mutation_rates.pkl', exclude_constitutive=exclude_constitutive)
        self.pickle_mutation_rates('back_subtracted_mutation_rates.pkl', subtract_background=True, exclude_constitutive=exclude_constitutive)
        self.pickle_mutation_rates('control_subtracted_mutation_rates.pkl', subtract_control=True, exclude_constitutive=exclude_constitutive)
        self.pickle_fold_changes('mutation_rate_fold_changes.pkl', exclude_constitutive=True)
        self.write_wigs('')
        self.write_wigs('back_subtract', subtract_background=True)
        self.write_wigs('control_subtract', subtract_control=True)
        self.write_mutation_rates_tsv('mutation_rates.tsv', exclude_constitutive=exclude_constitutive)
        self.write_mutation_rates_tsv('back_subtracted_mutation_rates.tsv', subtract_background=True, exclude_constitutive=exclude_constitutive)
        self.write_mutation_rates_tsv('control_subtracted_mutation_rates.tsv', subtract_control=True, exclude_constitutive=exclude_constitutive)
        self.write_combined_mutation_rates_tsv()
        self.write_combined_mutation_rates_tsv(exclude_constitutive=True)


    def write_mutation_rates_tsv(self, suffix, subtract_background=False, subtract_control=False, exclude_constitutive=False):
        if subtract_background or subtract_control:
            libs_to_write = self.get_normalizable_libs()
        else:
            libs_to_write = self.libs
        if subtract_background == False and subtract_control == False:
            prefix = 'raw'
        elif subtract_background == True and subtract_control == False:
            prefix = 'background_subtracted'
        elif subtract_background == False and subtract_control == True:
            prefix = 'control_subtracted'

        if exclude_constitutive:
            for lib in libs_to_write:
                lib.write_tsv_tables(os.path.join(self.rdir_path('tables', prefix, 'exclude_constitutive'),
                                                  lib.lib_settings.sample_name+'_'+suffix[:-4]+'_exclude_constitutive'+suffix[-4:]),
                                     subtract_background=subtract_background, subtract_control=subtract_control, exclude_constitutive=exclude_constitutive)
        else:
            for lib in libs_to_write:
                lib.write_tsv_tables(os.path.join(self.rdir_path('tables', prefix), lib.lib_settings.sample_name+'_'+suffix),
                                     subtract_background=subtract_background, subtract_control=subtract_control, exclude_constitutive=exclude_constitutive)

    def write_combined_mutation_rates_tsv(self, subtract_background=False, subtract_control=False, exclude_constitutive=False):
        if subtract_background and subtract_control:
            raise SyntaxError('Cannot subtract background and control simultaneously')

        if subtract_background or subtract_control:
            libs_to_write = list(self.get_normalizable_libs())
        else:
            libs_to_write = list(self.libs)

        if subtract_background == False and subtract_control == False:
            prefix = 'raw_'
        elif subtract_background == True and subtract_control == False:
            prefix = 'background_subtracted_'
        elif subtract_background == False and subtract_control == True:
            prefix = 'control_subtracted_'
        if exclude_constitutive:
            f = open(self.rdir_path('tables', prefix+'all_datasets_exclude_constitutive.tsv'), 'w')

        else:
            f = open(self.rdir_path('tables', prefix+'all_datasets.tsv'), 'w')
        f.write('rRNA\tposition\tnucleotide\t%s\n' % ('\t'.join([lib.lib_settings.sample_name for lib in libs_to_write])))
        for rRNA_name in sorted(self.settings.rRNA_seqs.keys()):
            for position in range(len(self.settings.rRNA_seqs[rRNA_name])):
                nuc_identity = self.settings.rRNA_seqs[rRNA_name][position]
                nuc_values = []
                for lib in libs_to_write:
                    nucleotide = lib.get_nucleotide(rRNA_name, position+1)
                    assert nucleotide.identity == nuc_identity
                    if not subtract_background and not subtract_control:
                        nuc_values.append(nucleotide.mutation_rate)
                    elif subtract_background:
                        nuc_values.append(nucleotide.get_back_sub_mutation_rate())
                    elif subtract_control:
                        nuc_values.append(nucleotide.get_control_sub_mutation_rate())
                assert len(nuc_values) == len(libs_to_write)
                if exclude_constitutive and nucleotide.exclude_constitutive:
                    f.write('%s\t%d\t%s\t%s\n' % (rRNA_name, position+1, nuc_identity, '\t'.join(['' for nuc_value in nuc_values])))
                else:
                    f.write('%s\t%d\t%s\t%s\n' % (rRNA_name, position+1, nuc_identity, '\t'.join([str(nuc_value) for nuc_value in nuc_values])))
        f.close()

    def write_wigs(self, suffix, subtract_background=False, subtract_control=False):
        hrf_utils.make_dir(self.rdir_path('wigs'))
        if subtract_background or subtract_control:
            libs_to_write = self.get_normalizable_libs()
        else:
            libs_to_write = self.libs
        #will also write a file to make batch import into mochiview easier
        f = open(os.path.join(self.rdir_path('wigs'), 'mochi_batch_'+suffix+'.txt'), 'w')
        f.write('SEQUENCE_SET\tFILE_NAME\tDATA_TYPE\tNAME\n')
        for lib in libs_to_write:
            f.write('<replace>\t%s\t<replace>\t%s\n' % (lib.lib_settings.sample_name+'_'+suffix+'.wig.gz', lib.lib_settings.sample_name+'_'+suffix))
            lib.write_mutation_rates_to_wig(os.path.join(self.rdir_path('wigs'), lib.lib_settings.sample_name+'_'+suffix),
                                      subtract_background=subtract_background, subtract_control=subtract_control)
        f.close()

    def pickle_mutation_rates(self, suffix, subtract_background=False, subtract_control=False, exclude_constitutive=False):
        if subtract_background:
            libs_to_pickle = self.get_modified_libs()
        elif subtract_control:
            libs_to_pickle = self.get_normalizable_libs()
        else:
            libs_to_pickle = self.libs

        if subtract_background == False and subtract_control == False:
            prefix = 'raw'
        elif subtract_background == True and subtract_control == False:
            prefix = 'background_subtracted'
        elif subtract_background == False and subtract_control == True:
            prefix = 'control_subtracted'

        if exclude_constitutive:
            for lib in libs_to_pickle:
                lib.pickle_mutation_rates(os.path.join(self.rdir_path('pickles', prefix, 'exclude_constitutive'),
                                                       lib.lib_settings.sample_name+'_'+suffix[:-4]+'_exclude_constitutive'+suffix[-4:]),
                                                       subtract_background=subtract_background, subtract_control=subtract_control,
                                                       exclude_constitutive=exclude_constitutive)
        else:
            for lib in libs_to_pickle:
                lib.pickle_mutation_rates(os.path.join(self.rdir_path('pickles', prefix), lib.lib_settings.sample_name+'_'+suffix),
                                          subtract_background=subtract_background, subtract_control=subtract_control, exclude_constitutive=exclude_constitutive)

    def pickle_fold_changes(self, suffix, exclude_constitutive=False):
        libs_to_pickle = self.get_normalizable_libs()

        if exclude_constitutive:
            for lib in libs_to_pickle:
                lib.pickle_mutation_fold_change(os.path.join(self.rdir_path('pickles', 'fold_change', 'exclude_constitutive'),
                                                       lib.lib_settings.sample_name + '_' + suffix[
                                                                                            :-4] + '_exclude_constitutive' + suffix[
                                                                                                                             -4:]), exclude_constitutive=exclude_constitutive)
        else:
            for lib in libs_to_pickle:
                lib.pickle_mutation_fold_change(
                    os.path.join(self.rdir_path('pickles', 'fold_change'), lib.lib_settings.sample_name + '_' + suffix), exclude_constitutive=exclude_constitutive)


    def make_plots(self, exclude_constitutive=False):
        if exclude_constitutive:
            hrf_utils.make_dir(self.rdir_path('plots', 'exclude_constitutive'))
            hrf_utils.make_dir(self.rdir_path('plots', 'exclude_constitutive', 'functional_groups'))
            hrf_utils.make_dir(self.rdir_path('plots', 'exclude_constitutive', 'interactive'))
            rdir = self.rdir_path('plots','exclude_constitutive')
            file_tag = '_exclude_constitutive'
        else:
            hrf_utils.make_dir(self.rdir_path('plots'))
            hrf_utils.make_dir(self.rdir_path('plots', 'interactive'))
            rdir = self.rdir_path('plots')
            file_tag = ''

        hrf_plotting.plot_mutated_nts_pie(self.libs, os.path.join(rdir, 'raw_mutation_fractions' + file_tag), exclude_constitutive=exclude_constitutive)
        hrf_plotting.plot_mutation_breakdown_pie(self.libs, os.path.join(rdir, 'raw_mutation_types' + file_tag), exclude_constitutive=exclude_constitutive)

        hrf_plotting.plot_mutated_nts_pie(self.libs,
                                          os.path.join(rdir, 'background_sub_mutation_fractions'+file_tag),
                                          subtract_background = True, exclude_constitutive=exclude_constitutive)
        hrf_plotting.plot_mutation_rate_cdfs(self.libs, os.path.join(rdir, 'mutation_rate_cdf' + file_tag),
                                             nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                             exclude_constitutive=exclude_constitutive)

        hrf_plotting.plot_mutation_rate_violins(self.libs, os.path.join(rdir, 'mutation_rate_cdf' + file_tag),
                                                nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                                exclude_constitutive=exclude_constitutive)

        hrf_plotting.ma_plots(self.get_normalizable_libs(), os.path.join(rdir, 'MA' + file_tag),
                              nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                              exclude_constitutive=exclude_constitutive)
        if self.settings.get_property('make_interactive_plots'):

                hrf_plotting.ma_plots_interactive(self.get_normalizable_libs(), os.path.join(rdir, 'interactive', 'MA' + file_tag),
                                                  nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                                  exclude_constitutive=False)


    def annotate_structures(self, exclude_constitutive=False):
        if exclude_constitutive:
            hrf_utils.make_dir(self.rdir_path('structures', 'protections_highlighted', 'exclude_constitutive'))
            hrf_utils.make_dir(self.rdir_path('structures', 'colored_by_change', 'exclude_constitutive'))
            file_tag = '_exclude_constitutive'
        else:
            hrf_utils.make_dir(self.rdir_path('structures', 'protections_highlighted'))
            hrf_utils.make_dir(self.rdir_path('structures', 'colored_by_change'))
            file_tag = ''
        if exclude_constitutive:
            hrf_plotting.highlight_structure(self.get_normalizable_libs(), self.rdir_path('structures', 'protections_highlighted', 'exclude_constitutive'),
                                             nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                             exclude_constitutive=exclude_constitutive)
        else:
            hrf_plotting.highlight_structure(self.get_normalizable_libs(), self.rdir_path('structures', 'protections_highlighted'),
                                             nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                             exclude_constitutive=exclude_constitutive)

    def collapse_identical_reads(self):
        """
        collapses all identical reads using FASTX toolkit
        :return:
        """
        self.settings.write_to_log('collapsing reads')
        if not self.settings.get_property('force_recollapse'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.collapsed_reads_exist():
                    break
            else:
                return
        hrf_utils.make_dir(self.rdir_path('collapsed_reads'))
        if self.settings.get_property('collapse_identical_reads'):
            hrf_utils.parmap(lambda lib_setting: self.collapse_one_fastq_file(lib_setting), self.settings.iter_lib_settings(),
                             nprocs = self.threads)
        else:
            hrf_utils.parmap(lambda lib_setting: self.fastq_to_fasta(lib_setting), self.settings.iter_lib_settings(),
                             nprocs = self.threads)
        self.settings.write_to_log('collapsing reads complete')

    def collapse_one_fastq_file(self, lib_settings):
        lib_settings.write_to_log('collapsing_reads')
        subprocess.Popen('gunzip -c %s | fastx_collapser -v -Q33 2>>%s | gzip > %s' % (lib_settings.get_fastq_file(),
                                                                                  lib_settings.get_log(),
                                                                                  lib_settings.get_collapsed_reads()
                                                                                  ), shell=True).wait()
        lib_settings.write_to_log('collapsing_reads_done')

    def fastq_to_fasta(self, lib_settings):
        lib_settings.write_to_log('fasta_conversion')
        subprocess.Popen('gunzip -c %s | fastq_to_fasta -v -Q33 2>>%s | gzip > %s' % (lib_settings.get_fastq_file(),
                                                                                  lib_settings.get_log(),
                                                                                  lib_settings.get_collapsed_reads()
                                                                                  ), shell=True).wait()
        lib_settings.write_to_log('fasta_conversion done')

    def get_barcode_match(self, barcode, barcodes):
        """
        takes a barcode and returns the one it matches (hamming <= 1)
        else
        empty string
        """
        if barcode in barcodes:
            return barcode
        for barcode_j in barcodes:
            if hrf_utils.hamming_N(barcode, barcode_j) <= self.settings.get_property('mismatches_allowed_in_barcode'):
                return barcode_j
        return ''

    def rdir_path(self, *args):
        return os.path.join(self.settings.get_rdir(), *args)

    def get_rdir_fhandle(self, *args):
        """
        returns a filehandle to the fname in the rdir
        """
        out_path = self.rdir_path(*args)
        out_dir = os.path.dirname(out_path)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        return hrf_utils.aopen(out_path, 'w')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("settings_file")
    parser.add_argument("--threads",
                        help="Max number of processes to use",
                        type = int, default = 3)
    args = parser.parse_args()

    return args

def main():
    """
    """
    args = parse_args()
    settings = hrf_settings.hrf_settings(args.settings_file)
    all_datasets = hrf_seq_run(settings, args.threads)


main()
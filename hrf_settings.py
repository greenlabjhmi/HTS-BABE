import os
import ConfigParser
import simplejson
import itertools
import shutil
import datetime

import hrf_utils

#TODO: make setting methods reflect new settings file
#TODO: remove superfluous settings methods
class hrf_settings:
    def __init__(self, settings_file):
        self.settings_file = settings_file
        self.process_settings(settings_file)
        self.rRNA_seqs = hrf_utils.convertFastaToDict(self.get_rRNA_fasta())

    def get_settings_file(self):
        return self.settings_file

    def get_property(self, property, default=None):
        try:
            if not property in self.settings and default != None:
                return default
            return self.settings[property]
        except:
            print self.settings
            raise  ValueError('cannot find %s' % property)

    def get_rdir(self):
        hrf_utils.make_dir(self.rdir)
        return self.rdir

    def iter_lib_settings(self):
        for i in range(len(self.sample_names)):
            yield hrf_lib_settings(self,
                                   self.sample_names[i],
                                   self.read1_fastq_gz_file_handles[i],
                                   self.read2_fastq_gz_file_handles[i])

    def process_settings(self, settings_file):
        """
        - reads the settings file and converts str to float, list, etc.
        - stores result in self.settings as a dict()
        """
        int_keys = [ 'first_base_to_keep', 'last_base_to_keep', 'min_post_adaptor_length']
        #float_keys = ['confidence_interval_cutoff', 'fold_change_cutoff']
        str_keys = ['experiment_name', 'read1_suffix', 'read2_suffix', 'read1_adaptor_sequence', 'read2_adaptor_sequence']
        boolean_keys = ['paired_end', 'collapse_identical_reads', 'force_read_resplit', 'force_remapping', 'force_recollapse',
                        'force_recount', 'force_index_rebuild', 'force_retrim', 'trim_adaptor',
                        'make_interactive_plots']
        list_str_keys = ['fastq_gz_prefixes', 'sample_names']
        extant_files = ['rrna_fasta', 'pymol_base_script', 'pymol_base_script_colorchange']
        config = ConfigParser.ConfigParser()
        config.read(settings_file)
        settings = {}
        for section in config.sections():
            for option in config.options(section):
                settings[option] = config.get(section, option)
                settings[section] = True
        for k in int_keys:
            settings[k] = int(settings[k])
        for k in str_keys:
            settings[k] = settings[k]
        #for k in float_keys:
        #    settings[k] = float(settings[k])
        for k in boolean_keys:
            if not settings[k].lower() in ['true', 'false']:
                raise ValueError(
                  'Boolean value %s must be "true" or "false"' % k)
            settings[k] = settings[k].lower() == 'true'
        #for k in list_float_keys:
        #    settings[k] = map(float, simplejson.loads(settings[k]))
        #for k in list_int_keys:
        #    settings[k] = map(int, simplejson.loads(settings[k]))
        for k in list_str_keys:
            settings[k] = simplejson.loads(settings[k])
        self.fqdir = settings['fastq_dir']
        self.sample_names = settings['sample_names']
        '''
        self.experimentals = settings['experimentals']
        self.no_mod_controls = settings['no_mod_controls']
        self.with_mod_controls = settings['with_mod_controls']
        self.exclude_constitutive = settings['exclude_constitutive']
        try:
            assert len(self.experimentals) == len(self.no_mod_controls)
            assert len(self.experimentals) == len(self.with_mod_controls)
        except:
            print 'error: experimentals, no_mod_controls, and with_mod_controls should all be the same length'
            print 'for mutation rate purposes, its ok to reuse a dataset here, it really doesnt matter'
        try:
            for sample_name in self.experimentals+self.no_mod_controls+self.with_mod_controls:
                assert sample_name in self.sample_names
        except:
            print sample_name, ' not in sample names, make sure you are using regular quotation marks'
        '''
        self.fastq_gz_read1_files = [fastq_gz_prefix + settings['read1_suffix'] for fastq_gz_prefix in
                                     settings['fastq_gz_prefixes']]
        self.fastq_gz_read2_files = [fastq_gz_prefix + settings['read2_suffix'] for fastq_gz_prefix in
                                     settings['fastq_gz_prefixes']]
        if settings['paired_end']:
            self.fastq_gz_files = self.fastq_gz_read1_files + self.fastq_gz_read2_files
        else:
            self.fastq_gz_files = self.fastq_gz_read1_files

        self.read1_fastq_gz_file_handles = [os.path.join(self.fqdir, fastq_gz_file)
                                            for fastq_gz_file in self.fastq_gz_read1_files]
        self.read2_fastq_gz_file_handles = [os.path.join(self.fqdir, fastq_gz_file)
                                            for fastq_gz_file in self.fastq_gz_read2_files]

        self.fastq_gz_file_handles = [os.path.join(self.fqdir, fastq_gz_file) for fastq_gz_file in self.fastq_gz_files]

        for file_handle in self.fastq_gz_file_handles:
            assert hrf_utils.file_exists(file_handle)
        for k in extant_files:
            assert hrf_utils.file_exists(settings[k])
        for file_handle in self.fastq_gz_file_handles:
            assert hrf_utils.file_exists(file_handle)
        self.settings = settings
        self.rdir = settings['results_dir']
        hrf_utils.make_dir(self.rdir)
        shutil.copy(settings_file, self.rdir)

    def get_rRNA_fasta(self):
        return self.get_property('rrna_fasta')

    def get_rRNA_STAR_index(self):
        index = os.path.join(
          self.get_rdir(),
          'star_indices',
          'rrna_index')
        return index

    def rRNA_STAR_index_exists(self):
        return hrf_utils.file_exists(self.get_rRNA_STAR_index())

    def get_log(self):
        log = os.path.join(
          self.get_rdir(),
          'log.txt')
        return log

    def write_to_log(self, text, add_time = True):
        f = open(self.get_log(), 'a')
        now = datetime.datetime.now()
        time = now.strftime("%Y-%m-%d %H:%M")
        if add_time:
            f.write('[%s] %s\n' % (time, text))
        else:
            f.write(text)
        f.close()

class hrf_lib_settings:
    def __init__(self, experiment_settings, sample_name, read1_fastq_gz_filehandle, read2_fastq_gz_filehandle):
        self.experiment_settings = experiment_settings
        self.sample_name = sample_name
        self.read1_fastq_gz_filehandle = read1_fastq_gz_filehandle
        self.read2_fastq_gz_filehandle = read2_fastq_gz_filehandle

    def get_property(self, property):
        return self.experiment_settings.get_property(property)

    def get_log(self):
        hrf_utils.make_dir(os.path.join(self.experiment_settings.get_rdir(), 'logs'))
        log = os.path.join(
          self.experiment_settings.get_rdir(),
          'logs',
          '%(sample_name)s.log' %
           {'sample_name': self.sample_name})
        return log

    def write_to_log(self, text, add_time = True):
        f = open(self.get_log(), 'a')
        now = datetime.datetime.now()
        time = now.strftime("%Y-%m-%d %H:%M")
        if add_time:
            f.write('[%s] %s\n' % (time, text))
        else:
            f.write(text)
        f.close()

    def get_fastq_gz(self):
        return self.read1_fastq_gz_filehandle, self.read2_fastq_gz_filehandle

    def get_collapsed_reads(self):
        collapsed_reads = os.path.join(
          self.experiment_settings.get_rdir(),
          'collapsed_reads',
          '%(sample_name)s.fasta.gz' %
           {'sample_name': self.sample_name})
        return collapsed_reads

    def get_adaptor_trimmed_reads(self, prefix_only=False):
            #default is false unless otherwise specified when calling fxn
        if prefix_only:
            adaptor_trimmed_reads_1 = os.path.join(
                self.experiment_settings.get_rdir(),
                'adaptor_trimmed_reads',
                 '%(sample_name)s_1.adaptor' %
                 {'sample_name': self.sample_name})
            adaptor_trimmed_reads_2 = os.path.join(
                self.experiment_settings.get_rdir(),
                'adaptor_trimmed_reads',
                '%(sample_name)s_2.adaptor' %
                {'sample_name': self.sample_name})
            return adaptor_trimmed_reads_1, adaptor_trimmed_reads_2
        else:
            adaptor_trimmed_reads_1 = os.path.join(
                self.experiment_settings.get_rdir(),
                'adaptor_trimmed_reads',
                '%(sample_name)s_1.adaptor-trimmed.fastq.gz' %
                {'sample_name': self.sample_name})
            adaptor_trimmed_reads_2 = os.path.join(
                self.experiment_settings.get_rdir(),
                'adaptor_trimmed_reads',
                '%(sample_name)s_2.adaptor-trimmed.fastq.gz' %
                {'sample_name': self.sample_name})
            return adaptor_trimmed_reads_1, adaptor_trimmed_reads_2

    def get_rRNA_mapping_stats(self):
        rRNA_mapping_stats = os.path.join(
          self.experiment_settings.get_rdir(),
          'mapping_stats',
          '%(sample_name)s.rRNA.txt' %
           {'sample_name': self.sample_name})
        return rRNA_mapping_stats

    def get_mapped_reads(self):
        mapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'mapped_reads', '%(sample_name)sAligned.sortedByCoord.out.bam' % {'sample_name': self.sample_name})
        return mapped_reads

    def get_mapped_reads_prefix(self):
        mapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'mapped_reads', '%(sample_name)s' % {'sample_name': self.sample_name})
        return mapped_reads


    def get_unmappable_reads(self):
        unmapped_reads = os.path.join(
          self.experiment_settings.get_rdir(),
          'unmapped_reads',
          '%(sample_name)s.unmappable.fasta.gz' %
           {'sample_name': self.sample_name})
        return unmapped_reads


    def get_trimmed_reads(self):
        trimmed_reads_1 = os.path.join(
          self.experiment_settings.get_rdir(),
           'trimmed_reads',
           '%(sample_name)s_1.trimmed.fastq.gz' %
            {'sample_name': self.sample_name})
        trimmed_reads_2 = os.path.join(
           self.experiment_settings.get_rdir(),
           'trimmed_reads',
           '%(sample_name)s_2.trimmed.fastq.gz' %
          {'sample_name': self.sample_name})
        return trimmed_reads_1, trimmed_reads_2

    def get_filtered_reads(self):
        trimmed_reads = os.path.join(
          self.experiment_settings.get_rdir(),
          'quality_filtered_reads',
          '%(sample_name)s.filtered.fastq.gz' %
           {'sample_name': self.sample_name})
        return trimmed_reads

    def get_counting_prefix(self):
        return os.path.join(
          self.experiment_settings.get_rdir(),
          'read_counts',
          '%(sample_name)s' %
           {'sample_name': self.sample_name})

    def get_read_5p_counts(self):
        return os.path.join(
          self.experiment_settings.get_rdir(),
          'read_counts',
          '%(sample_name)s.5p_ends.pkl' %
           {'sample_name': self.sample_name})

    def get_normalized_mutation_counts(self):
        return os.path.join(
          self.experiment_settings.get_rdir(),
          'normalized_mutation_counts',
          '%(sample_name)s.norm_mut.pkl' %
           {'sample_name': self.sample_name})

    def read_5p_counts_exists(self):
        return hrf_utils.file_exists(self.get_read_5p_counts())

    def get_positional_coverage(self):
        return os.path.join(
          self.experiment_settings.get_rdir(),
          'read_counts',
          '%(sample_name)s.coverage.pkl' %
           {'sample_name': self.sample_name})

    def positional_coverage_exists(self):
        return hrf_utils.file_exists(self.get_positional_coverage())

    def get_mutation_counts(self):
        return os.path.join(
          self.experiment_settings.get_rdir(),
          'read_counts',
          '%(sample_name)s.mutations.pkl' %
           {'sample_name': self.sample_name})

    def mutation_counts_exists(self):
        return hrf_utils.file_exists(self.get_mutation_counts())

    def counts_all_exist(self):
        return self.mutation_counts_exists() and self.positional_coverage_exists() and self.read_5p_counts_exists()

    def get_overall_contamination_summary(self):
        summary_file = os.path.join(
          self.experiment_settings.get_rdir(),
          'QC',
          '%(sample_name)s.contamination_summary.txt' %
           {'sample_name': self.sample_name})
        return summary_file

    def split_reads_exist(self):
        split_reads = self.get_split_reads()
        return hrf_utils.file_exists(split_reads)

    def collapsed_reads_exist(self):
        collapsed_reads = self.get_collapsed_reads()
        return hrf_utils.file_exists(collapsed_reads)

    def adaptorless_reads_exist(self):
        if self.experiment_settings.get_property('paired_end'):
            adaptorless_reads_1 = self.get_adaptor_trimmed_reads()[0]
            adaptorless_reads_2 = self.get_adaptor_trimmed_reads()[1]
            return hrf_utils.file_exists(adaptorless_reads_1) and hrf_utils.file_exists(adaptorless_reads_2)
        else:
            adaptorless_reads_1 = self.get_adaptor_trimmed_reads()[0]
            return hrf_utils.file_exists(adaptorless_reads_1)

    def primerless_reads_exist(self):
        primerless_reads = self.get_primer_trimmed_reads()
        return hrf_utils.file_exists(primerless_reads)

    def trimmed_reads_exist(self):
        if self.experiment_settings.get_property('paired_end'):
            trimmed_reads_1 = self.get_trimmed_reads()[0]
            trimmed_reads_2 = self.get_trimmed_reads()[1]
            return hrf_utils.file_exists(trimmed_reads_1) and hrf_utils.file_exists(trimmed_reads_2)
        else:
            trimmed_reads_1 = self.get_trimmed_reads()[0]
            return hrf_utils.file_exists(trimmed_reads_1)

    def filtered_reads_exist(self):
        filtered_reads = self.get_filtered_reads()
        return hrf_utils.file_exists(filtered_reads)

    def mapped_reads_exist(self):
        mapped_reads = self.get_mapped_reads()
        return hrf_utils.file_exists(mapped_reads)

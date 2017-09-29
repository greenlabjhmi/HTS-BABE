'''
When the GRCF sequences a pool on a full flowcell, they demultiplex each lane seperately.
This code will use the names to combine the files for the same sample from the 2 lanes

'''

import sys, os, subprocess

def files_matching_barcode(barcode, filenames):
    matches = []
    for filename in filenames:
        test_barcode = filename.split('_')[2]
        if test_barcode == barcode:
            matches.append(filename)
    return matches

fastq_folder, prefix, output_folder = sys.argv[1:]

filenames = os.listdir(fastq_folder)
fastq_gz_files = [filename for filename in filenames if filename.endswith('fastq.gz')]
print '# fastq.gz files:', len(fastq_gz_files)
barcodes = set()
for filename in fastq_gz_files:
    barcode = filename.split('_')[2]
    barcodes.add(barcode)
print '# barcodes:', len(barcodes)

for barcode in barcodes:
    matching_files = [os.path.join(fastq_folder, file) for file in files_matching_barcode(barcode, fastq_gz_files)]
    outfile = os.path.join(output_folder, '%s%s.fastq.gz' % (prefix, barcode))
    command =  'gzcat %s | gzip -c > %s' % (' '.join(matching_files), outfile)
    print command
    subprocess.Popen(command, shell=True).wait()

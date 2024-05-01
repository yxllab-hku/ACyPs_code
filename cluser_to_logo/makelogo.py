import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker as lm
import os, time, argparse, concurrent.futures
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO

try:
    import StringIO
except ImportError:
    from io import StringIO

start = time.perf_counter()


def makelogo(fasta_file,out_dir):
    mafft_cline = MafftCommandline(input=fasta_file)
    stdout, stderr = mafft_cline()
    prefix = os.path.basename(fasta_file)
    aligned_path = out_dir + f'/{prefix}.aln'
    with open(aligned_path, "w") as handle:
        handle.write(stdout)
    align = AlignIO.read(aligned_path, "fasta")
    seqs = []
    for seq_record in align:
        seqs.append(str(seq_record.seq))

    pre_aln_df = lm.alignment_to_matrix(sequences=seqs, to_type='counts', characters_to_ignore='.-X')
    pre_aln_logo = lm.Logo(pre_aln_df, font_name='sans',
                           fade_below=0.5, shade_below=0.5,
                           figsize=(25, 2),  # <= 50AA (15, 2)  <=100 (25,2)
                           stack_order='big_on_top',
                           color_scheme='NajafabadiEtAl2017', fade_probabilities=True,baseline_width=1.5,show_spines=False)
    plt.savefig(out_dir+ f'/{prefix}.png')
    plt.savefig(out_dir + f'/{prefix}.svg')
    plt.close()


if __name__ == '__main__':
    '''
    input aln sequences must be in one line.
    '''
    parser = argparse.ArgumentParser(description='Perform sequence alignment and then generate logo sequence .svg file to display sequence biases')
    parser.add_argument('-i', '--input_fasta', type=str, metavar='', required=False,
                        help='Folder contains to-be-processed protein or nucleotide sequence files, e.g., /data/protein/')
    parser.add_argument('-o', '--outdir', type=str, metavar='', required=True,
                        help='A  path to save all files, including .aln alignment file and .svg/.png seq log file. e.g., data/output')
    parser.add_argument('-n', '--num_threads', type=int, metavar='', required=False,default=1,
                        help='<-n 4> Number of threads, 1 thread  (default n=1)  will be used if not given')

    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    with concurrent.futures.ProcessPoolExecutor(args.num_threads) as executor:  # 24 threads
        for fasta_file in os.listdir(args.input_fasta):
            print(f'Processing {fasta_file}')
            executor.submit(makelogo,args.input_fasta+fasta_file, args.outdir)
            print(f'Finished {fasta_file}')

    finish = time.perf_counter()
    print(f'finished in {round(finish - start, 3)} seconds')
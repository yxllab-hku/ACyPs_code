import os
import time
import networkx
from Bio import SeqIO
from networkx.algorithms.components.connected import connected_components
import argparse

start = time.perf_counter()

# For parsing network generated by mmseq2, three functions bellow is a credit to
# Jochen Ritzel (https://stackoverflow.com/users/95612/jochen-ritzel)


def gen(cluster_list):
    with open(cluster_list, 'r') as f:
        cont = []
        for line in f:
            cont.append(line)
        s = ''.join(cont).replace('\n\t', '\t').split('\n')

    cluster = []
    for sub_cluster in s:
        if len(sub_cluster) != 0:
            cluster.append(sub_cluster.split('\t'))
    return cluster


def to_graph(list_of_connection):
    G = networkx.Graph()
    for part in list_of_connection:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also imlies a number of edges:
        G.add_edges_from(to_edges(part))
    return G


def to_edges(list_of_connection):
    """
        treat `list_of_connection` as a Graph and returns it's edges
        to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it = iter(list_of_connection)
    last = next(it)

    for current in it:
        yield last, current
        last = current
# G = to_graph(list_of_connection)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Parsing network file generated by mmseq2, input .tsv network file and .fasta file, output one-cluster-to-one-fasta_file. E.g., 0001_8253_protease.fasta referse to first cluster (0001) with 8253 sequences (nodes) sequence contains a member named protease')
    parser.add_argument('-i', '--input_tsv', type=str, metavar='', required=False,
                        help='Network file generated by mmseq2, in general, it is a .tsv file')
    parser.add_argument('-fasta', '--fasta_file', type=str, metavar='', required=False,
                        help='Original protein fasta file used for clustering, e.g., /data/test.fasta')
    parser.add_argument('-node', '--node_num', type=int, metavar='', required=False,
                        help='Minimal number (included) of nodes in cluster')
    parser.add_argument('-o', '--outdir', type=str, metavar='', required=True,
                        help='A  path to save all files. e.g., data/output')

    args = parser.parse_args()

    # pair = gen(r_filenameTSV)
    # G = to_graph(pair)
    # cluster_list = list(connected_components(G))
    # sorted_cluster = sorted(cluster_list,key=lambda a: len(a), reverse=True)
    pair = gen(args.input_tsv)
    G = to_graph(pair)
    cluster_list = list(connected_components(G))
    sorted_cluster = sorted(cluster_list,key=lambda a: len(a), reverse=True)
    id2seq = {}
    for seq_record in SeqIO.parse(args.fasta_file,'fasta'):
        id2seq[seq_record.id] = str(seq_record.seq)

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    n = 0
    minimal_node = args.node_num-1
    for i in sorted_cluster:
        if len(list(i)) > minimal_node:
            # print(len(i))
            n += 1
            with open(args.outdir + '/' + str(n).zfill(5) + '_' + str(len(i)) + '_' + list(i)[0] + '.fasta', 'w') as f_out:
                for each_seq in i:
                    f_out.write('>' + each_seq + '\n')
                    f_out.write(id2seq[each_seq] + '\n')
                f_out.close()

    finish = time.perf_counter()
    print(f'finished in {round(finish - start, 3)} seconds')
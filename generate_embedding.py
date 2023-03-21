import numpy as np
import random
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
import time
from sklearn.decomposition import PCA

from hsrl.graph import *
import hsrl.node2vec as node2vec
from hsrl.hier_samp import louvainModularityOptimization

def parse_args():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    parser.add_argument('--input',  
                        default='../datasets/BioGRID_network.txt',
                        help='input graph file')
    parser.add_argument('--output', 
                        default='../emds/BioGRID_network.txt',
                        help='output representation file')
    
    ## Embedding Size
    parser.add_argument('--representation-size',
                        default=256,
                        type=int,
                        help='number of latent dimensions to learn for each node')
    parser.add_argument('--method',
                        default='deepwalk',               
                        choices=['deepwalk', 'node2vec'],
                        help='the learning method')
    parser.add_argument('--hs_num',
                        default=3,
                        type=int,
                        help='the number of hierarchical sampling layers')
    
    ## Graph info
    parser.add_argument('--graph-format',
                        default='edgelist',
                        help='input graph format')
    parser.add_argument('--directed',
                        action='store_true',
                        help='treat graph as directed')
    parser.add_argument('--weighted',
                        action='store_true',
                        default=False,                        
                        help='treat graph as weighted')
    
    ## parameters for GraRep
    parser.add_argument('--Kstep',
                        default=4,
                        type=int,
                        help='use k-step transition probability matrix')
    
    ## parameters for deepwalk and note2vec
    parser.add_argument('--walk-length',
                        default=40,
                        type=int,
                        help='length of the random walk')
    parser.add_argument('--number-walks',
                        default=10,
                        type=int,
                        help='number of random walks to start at each node')
    parser.add_argument('--window-size',
                        default=5,
                        type=int,
                        help='window size of skipgram model')
    parser.add_argument('--workers',
                        default=10,
                        type=int,
                        help='number of parallel processes')
    parser.add_argument('--p',
                        default=1.0,
                        type=float)
    parser.add_argument('--q',
                        default=1.0,
                        type=float)
    
    args = parser.parse_args()

    return args

def buildGraph(args):
    g = Graph()
    g.read_edgelist(filename=args.input,
                    weighted=args.weighted,
                    directed=args.directed)

    return g

def buildModel(args, g, batch_size=128, hs=0):
    if args.method == 'deepwalk':
        model = node2vec.Node2vec(graph=g,
                                    path_length=args.walk_length,
                                    num_paths=args.number_walks,
                                    dim=args.representation_size,
                                    workers=args.workers,
                                    p=args.p,
                                    q=args.q,
                                    dw=True,
                                    window=args.window_size)
    elif args.method == 'node2vec':
        model = node2vec.Node2vec(graph=g,
                                    path_length=args.walk_length,
                                    num_paths=args.number_walks,
                                    dim=args.representation_size,
                                    workers=args.workers,
                                    p=args.p,
                                    q=args.q,
                                    window=args.window_size)
    
    return model

def merge_embeddings(embeddings1, embeddings2, graph1, graph2, comms, dim=None):
    embeddings = np.zeros((graph1.node_size, dim))

    for comm in comms.keys():
        for node in comms[comm]:
            embeddings[graph1.look_up_dict[node]] = np.concatenate((embeddings1[graph1.look_up_dict[node]], embeddings2[graph2.look_up_dict[comm]]))

    return embeddings

def get_embeddings(graph, embeddings):
    look_back = graph.look_back_list
    vectors = {}
    for i, embedding in enumerate(embeddings):
        vectors[look_back[i]] = embedding

    return vectors

def save_embeddings(args, graph, embeddings, filename):
    vectors = get_embeddings(graph, embeddings)
    fout = open(filename, 'w')
    for node, vec in vectors.items():
        fout.write('{} {}\n'.format(str(node), ' '.join([str(x) for x in vec])))

    fout.close()

    return vectors

def transPCA(emd_file, dim, pca_file,emd_size):
    protein_id = {}
    id_protein = {}
    emd = []
    with open(emd_file, 'r') as f:
        lines = f.readlines()
        for i,line in enumerate(lines):
            content = line.strip('\n').strip().split() 
            if len(content)-1 != emd_size:
                print("FFFFFFFFFFALSE", len(content), emd_size)
            emd_temp = []
            protein_name = content[0]
            new_id = len(protein_id)
            protein_id[protein_name] = new_id
            id_protein[new_id] = protein_name
            emd.append(content[1:])

    emd = np.array(emd)
    pca = PCA(n_components = dim)
    pca.fit(emd)
    print(pca.n_components_, emd.shape)
    emd = pca.transform(emd)
    print(emd.shape)
    aim_dim = {}
    for i, emd_temp in enumerate(emd):
        protein_name = id_protein[i]
        aim_dim[protein_name] = emd_temp
    
    with open(pca_file, 'w') as f:
        for protein_name in aim_dim.keys():
            emd_temp = aim_dim[protein_name]
            content = str(protein_name)
            for i in emd_temp:
                content += ' ' + str(i)
            f.writelines(content+'\n')

def main(args):
    g = buildGraph(args)
    print('hierarchical sampling...')
    t1 = time.time()    
    hier_graph = {}
    hier_graph[0] = g
    hier_emb = {}
    comms = {}
    
    k = args.hs_num
    model = buildModel(args, hier_graph[0])
    hier_emb[0] = model.embeddings
    for i in range(k):
        hier_graph[i+1], comms[i] = louvainModularityOptimization(hier_graph[i])       
        model = buildModel(args, hier_graph[i+1], hs=i+1)        
        hier_emb[i+1] = model.embeddings

    m = 0
    for j in range(k, 0, -1):      
        hier_emb[j-1] = merge_embeddings(hier_emb[j-1], hier_emb[j], hier_graph[j-1], hier_graph[j], comms[j-1], args.representation_size*(m+2))
        m += 1
        
    if args.output:
        vectors = get_embeddings(hier_graph[0], hier_emb[0])    
        all_emd_file = './emds/temp_embedding.txt'
        save_embeddings(args, hier_graph[0], hier_emb[0], all_emd_file)
        print("Final Embedding Size:{0}".format(args.representation_size//(k+1)))
        transPCA(all_emd_file, args.representation_size//(k+1) , args.output, emd_size)

    t2 = time.time()
    print('cost time: %s'%(t2-t1))

if __name__ == '__main__':
    random.seed(128)
    np.random.seed(128)
    main(parse_args())
    
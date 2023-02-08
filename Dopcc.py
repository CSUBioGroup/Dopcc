import time
import numpy as np
import math
from collections import defaultdict
import pandas as pd
from tqdm.autonotebook import tqdm,trange
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Core attachment based on Local and Global information.')
    parser.add_argument('--input', default='datasets/BioGRID_network.txt')
    parser.add_argument('--output', default='results/temp.txt')
    parser.add_argument('--emd', default='emds/BioGRID_network_1024_redu_256.txt')
    
    args = parser.parse_args()
    return args

def load_data(PPIfile):
    edge_num = 0
    relations = defaultdict(list)
    protein_id, id_protein = {}, {}
    
    with open(PPIfile, 'r') as f:
        content = f.readlines()
        for line in content:
            line = line.strip('\n').split('\t')
            protein_source = line[0]
            protein_destination = line[1]
            
            if protein_source not in protein_id.keys():
                newid = len(protein_id)
                protein_id[protein_source] = newid
                id_protein[newid] = protein_source
            if protein_destination not in protein_id.keys():
                newid = len(protein_id)
                protein_id[protein_destination] = newid
                id_protein[newid] = protein_destination
            
            protein_source_id = protein_id[protein_source]
            protein_destination_id = protein_id[protein_destination]
            if protein_source_id == protein_destination_id:
                continue
            if protein_destination_id not in relations[protein_source_id]:
                relations[protein_source_id].append(protein_destination_id)
                relations[protein_destination_id].append(protein_source_id)
                edge_num += 1
    return relations, protein_id, id_protein, edge_num

def cal_cosine(x, y):
    assert len(x) == len(y), "emd error"
    res = np.array( [[x[i]*y[i], x[i]*x[i], y[i]*y[i]] for i in range(len(x))] )
    cos = sum(res[:, 0]) / (np.sqrt(sum(res[:, 1])) * np.sqrt(sum(res[:, 2])))
    return cos        

def cal_similarity(v, u, relations):
    N_v = relations[v]
    N_u = relations[u]
    SN_v = set(N_v)
    SN_v.add(v)
    SN_u = set(N_u)
    SN_u.add(u)
    inter = SN_v & SN_u
    SS_vu = float(len(inter)) / math.sqrt(len(SN_v)*len(SN_u))
    return SS_vu

def cal_sumEcore(core_pc, relations, w):
    sum_Ecore = 0.0
    proteins = sorted(core_pc)
    proteins_set = set(proteins)
    for v in proteins:
        N_v = set(relations[v])
        inter = N_v & proteins_set
        inter = list(inter)
        inter = sorted(inter)
        for u in inter:
            if u > v:
                sum_Ecore += w[v][u]
    return sum_Ecore


def attachment_algorithm(core_pc, relations, w_nb, w_emd):
    attachment = set()
    overlapping_attachments = set()
    peripery_attachments = set()
    
    sum_Ecore_nb = cal_sumEcore(core_pc, relations, w_nb)
    w_avg_core_nb = 2.0*sum_Ecore_nb / len(core_pc)
    sum_Ecore_emd = cal_sumEcore(core_pc, relations, w_emd)
    w_avg_core_emd = 2.0*sum_Ecore_emd / len(core_pc)
    
    for u in core_pc:
        N_u = set(relations[u])
        attachment = attachment | N_u
    attachment = attachment - set(core_pc)
    
    attachment_nodes_dict_in_nb = {}
    attachment_nodes_dict_out_nb = {}
    attachment_nodes_dict_in_emd = {}
    attachment_nodes_dict_out_emd = {}
    
    candidate_attachment_proteins = []
    overlapping_node, peripery_node = [], []
    
    for p in attachment:
        count = 0
        N_p = relations[p]
        sum_p_core_in_nb = 0.0
        sum_p_core_out_nb = 0.0
        sum_p_core_in_emd = 0.0
        sum_p_core_out_emd = 0.0
        for t in N_p:
            if t in core_pc:
                sum_p_core_in_nb += w_nb[p][t]
                sum_p_core_in_emd += w_emd[p][t]
                count += 1
            else:
                sum_p_core_out_nb += w_nb[p][t]
                sum_p_core_out_emd += w_emd[p][t]
        if count>=2:
            candidate_attachment_proteins.append(p)
            attachment_nodes_dict_in_nb[p] = sum_p_core_in_nb
            attachment_nodes_dict_out_nb[p] = sum_p_core_out_nb
            attachment_nodes_dict_in_emd[p] = sum_p_core_in_emd
            attachment_nodes_dict_out_emd[p] = sum_p_core_out_emd
            
            if (sum_p_core_out_nb >= sum_p_core_in_nb) and (sum_p_core_out_emd >= sum_p_core_in_emd):
                overlapping_node.append(p)
            else:
                peripery_node.append(p)
    # save overlapping node
    if len(overlapping_node) > 0:
        for p in overlapping_node:
            if (attachment_nodes_dict_in_nb[p] >= 0.5*w_avg_core_nb) and (attachment_nodes_dict_in_emd[p] >= 0.5*w_avg_core_emd):
                overlapping_attachments.add(p)
    # save peripery node
    if len(peripery_node) > 0:
        sum_c_core_avg_nb = 0.0
        sum_c_core_avg_emd = 0.0
        for c in peripery_node:
            sum_c_core_avg_nb += attachment_nodes_dict_in_nb[c]
            sum_c_core_avg_emd += attachment_nodes_dict_in_emd[c]
        w_avg_cp_nb = sum_c_core_avg_nb / len(peripery_node)
        w_avg_cp_emd = sum_c_core_avg_emd / len(peripery_node)
        for p in peripery_node:
            if (attachment_nodes_dict_in_nb[p] >= w_avg_cp_nb) and (attachment_nodes_dict_in_emd[p] >= w_avg_cp_emd):
                peripery_attachments.add(p)
    attachment = peripery_attachments | overlapping_attachments
    attachment_pc = list(attachment)
    return attachment_pc

def Complex_algorithm(Core_PCs, relations, w_nb, w_emd):
    complexes = {}
    for v in tqdm(Core_PCs.keys(), desc='Generate Attachments'):
        core_pc = Core_PCs[v]
        attachment_pc = attachment_algorithm(core_pc, relations, w_nb, w_emd)
        complex_v = list(set(core_pc + attachment_pc))
        complex_v.sort()
        complexes[v] = complex_v
    return complexes

def cal_complex_similarity(cp1, cp2):
    inter = cp1 & cp2
    union = cp1 | cp2
    sm = len(inter) / len(union)
    return sm

def Complex_redundancy(complexes, relations):
    cp_redu = set()
    a = set(complexes.keys())
    b = set(complexes.keys())
    for u in tqdm(a, desc='Complex Redundancy'):
        if u not in cp_redu:
            for v in b:
                if v not in cp_redu:
                    if v > u:
                        set_u = set(complexes[u])
                        set_v = set(complexes[v])
                        sm = cal_complex_similarity(set_u, set_v)
                        if sm>=1:
                            cp_redu.add(v)
    complexes_new = {}
    for u in a:
        if u in cp_redu:
            continue
        else:
            complexes_new[u] = complexes[u]
    return complexes_new

def save_result(complexes, id_protein, result_file):
    with open(result_file, 'w') as f:
        rowid = 0
        for u in complexes.keys():
            cps = complexes[u]
            if len(cps) >= 3:
                rowid += 1
                content = str(rowid)
                for pid in cps:
                    content = content + ' ' + id_protein[pid]
                content +='\n'
                f.writelines(content)
    print('Save Results to file: ', result_file)

def cal_jcs(id_protein, relations, w):
    protein_num = len(id_protein)
    for u in trange(protein_num, desc='Calculating JCD'):
        N_u = relations[u]
        for v in range(u+1, protein_num):
            N_v = relations[v]
            if len(N_u)>1 or len(N_v)>1:
                inter = set(N_u) & set(N_v)
                union = set(N_u) | set(N_v)
                if len(inter) < 1:
                    w[u][v] = 0.0
                    w[v][u] = 0.0
                else:
                    w[u][v] = float(len(inter)) / len(union)
                    w[v][u] = float(len(inter)) / len(union)
            else:
                w[u][v] = 0.0
                w[v][u] = 0.0
    return w, relations

def construct_net_neighbor(id_protein, relations, w):
    w_jcs,relations = cal_jcs(id_protein, relations, w)
    w = w_jcs
    return w, relations

def construct_net_embedding(id_protein, protein_id, relations, w_emd, emd_file):
    emd = {}
    with open(emd_file, 'r') as f:
        content = f.readlines()
        for line in content:
            line = line.strip('\n').strip().split()
            emd_temp = []
            for i, val in enumerate(line):
                if i == 0:
                    protein_name = val
                else:
                    emd_temp.append(float(val))
            emd[protein_id[protein_name]] = emd_temp
        print("embedding size : %d"%(len(emd_temp)))
    
    protein_num = len(id_protein)
    for u in trange(protein_num, desc='Calculating Embedding similarity'):
        N_u = relations[u]
        for v in N_u:
            weight = cal_cosine(emd[u], emd[v])
            w_emd[u][v] = weight
    
    return w_emd, relations

def Core_algorithm(id_protein, relations, lamda):
    core_pcs = {}
    protein_num = len(id_protein)
    for v in trange(protein_num, desc='Generate Cores'):
        N_v = relations[v]
        Core_v = set()
        Core_v.add(v)
        for u in N_v:
            SS_vu = cal_similarity(v, u, relations)
            if SS_vu > lamda:
                Core_v.add(u)
        if len(Core_v) >=2 :
            core_pcs[v] = list(Core_v)
    return core_pcs

def main(args):
    PPI_file = args.input
    emd_file = args.emd
    result_file = args.output
    print(args)
    
    relations, protein_id, id_protein, edge_num = load_data(PPI_file)
    protein_num = len(protein_id)
    print('Protein_num: %d, Edge_num: %d'%(protein_num, edge_num))
    
    w_nb = np.zeros((protein_num, protein_num))
    w_nb, relations = construct_net_neighbor(id_protein, relations, w_nb)
    
    w_emd = np.zeros((protein_num, protein_num))
    w_emd, relations = construct_net_embedding(id_protein, protein_id, relations, w_emd, emd_file)
    
    lamda = 0.4
    Core_PCs = Core_algorithm(id_protein, relations, lamda)
    
    complexes = Complex_algorithm(Core_PCs, relations, w_nb, w_emd)
    
    complexes = Complex_redundancy(complexes, relations)
    
    save_result(complexes, id_protein, result_file = result_file)

    
    
if __name__ == '__main__':
    main(parse_args())
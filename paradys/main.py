import graph_tool.all as gt
import pandas as pd
import numpy as np
from scipy.stats import chi2
import multiprocessing
import os


# Import functions from other modules within the package
def dysregulation_condition(edge, patient_networks):
    '''
    Dysregulation condition

    Parameters
    ----------
    edge : string
    patient_networks : DataFrame

    Returns
    -------
    bool
        True if edge is present in the patient's dysregulatory network, False otherwise

    '''
    return edge in patient_networks['edge']

def mutation_condition(gene, patient_mutations):
    '''
    Mutation condition

    Parameters
    ----------
    gene : string
    patient_mutations : DataFrame

    Returns
    -------
    bool
        True if gene mutated for that given patient, False otherwise

    '''
    return gene in patient_mutations['gene']

def get_relations(kappa, patients, mutations, networks):
    ''''
    Returns and stores all connections between mutations and dysregulations

    Parameters
    ----------
    kappa : int
        Maximum path length between edge and mutated gene
    patients : list
        List of patient names.
    mutations : DataFrame
        DataFrame containing mutation data.
    networks : DataFrame
        DataFrame containing dysregulation network data.

    Returns
    -------
    data : DataFrame
        For every mutated gene for every patient, list of dysregulations connected to it in a max distance kappa

    '''
    data = []
    for p in patients:
        networks_p = networks[networks['patient'] == p]
        tfs = set(networks_p['tf'])  # Convert to set for faster membership testing

        mutations_p = mutations[mutations['patient'] == p]['gene']
        for mut in mutations_p:
            if mut in tfs:
                # kappa = 0
                targets_k = {mut}  # connected targets (use set for faster lookups)
                tfs_k = {mut}  # tfs to analyze (use set for faster lookups)

                for k in range(kappa):
                    if len(tfs_k) > 0:
                        new_tfs = set()
                        for tf in tfs_k:
                            if tf in tfs:
                                for g in networks_p[networks_p['tf'] == tf]['gene']:
                                    if g not in targets_k:
                                        targets_k.add(g)
                                        new_tfs.add(g)
                        tfs_k = new_tfs
                    else:
                        break

                data.append({
                    'connected targets': list(targets_k),
                    'mutations': mut,
                    'patients': p
                })

    data = pd.DataFrame(data)
    return data
    
def chi(edge, gene, patient, mutations, networks):
    '''
    Calculate chi2 statistic for a given gene and edge

    Parameters
    ----------
    edge : list
        List containing the edge (e.g., [source, target]).
    gene : string
        Gene to test.
    patient : string
        Patient name.
    mutations : DataFrame
        DataFrame containing mutation data.
    networks : DataFrame
        DataFrame containing dysregulation network data.

    Returns
    -------
    chi_val : float
        Chi-square value for that gene and edge.

    '''

    patient_mutations = mutations[mutations['patient'] == patient]
    patient_networks = networks[networks['patient'] == patient]

    # Create the observation matrix o using numpy vectorization
    o = np.zeros((2, 2), dtype=int)
    for _, row in patient_mutations.iterrows():
        dis = dysregulation_condition(edge, patient_networks)
        mut = mutation_condition(gene, row)
        o[dis][mut] += 1

    # Calculate the expected frequencies using the marginals
    row_sums = o.sum(axis=1)
    col_sums = o.sum(axis=0)
    total = o.sum()
    e = np.outer(row_sums, col_sums) / total

    # Perform the chi-square test using the chisquare function
    chi_val = chi2.sf(np.sum((o - e)**2 / e), 1)

    return chi_val

def create_statistics_p(p, d, scores, mutations, networks, patients, 
                        directed, tfs, edges_counts, output_dir):
    net = networks[networks['patient'] == p]
    data = []
    target_edges = set(net['edge'])  # Convert to set for faster lookups and avoid duplicates

    gene_sets_p = tfs[tfs['patients']==p]
    for mut in gene_sets_p['mutations']:
        scores_list = []
        dysregulations_list = []

        # Collect target genes for each TF-gene pair in gene_sets_p
        conn_targets_list = gene_sets_p[gene_sets_p['mutations'] == mut]['connected targets']
        print("Conn_targets:", conn_targets_list)
        mut_edges = net[net['gene'].isin(conn_targets_list)]
        print("Mut edges:", mut_edges)
        target_genes = mut_edges['gene']
        
        # Calculate chi-squared and p-value for each TF-gene-edge combination
        for target in target_genes:
            edge = [gene_sets_p.loc[gene_sets_p['connected targets'] == target, 'mutations'].values[0], target]
            target_edges.add(target)
            chi_g_e = chi(edge, target, p, mutations, networks)
            if chi_g_e > 3.84:
                p_value = chi2.sf(chi_g_e, 1)
                scores_list.append(p_value)
                dysregulations_list.append(edge)

        data.append({
            'driver': mut,
            'dysregulation': dysregulations_list,
            'p-value': scores_list
        })

    if len(data) > 0:
        data = pd.DataFrame(data)
        if scores:
            # Create line graph and apply PageRank algorithm to get drivers scores
            total_drivers = set(data['driver'])
            total_vertices = total_drivers.union(target_edges)
            total_vertices = list(total_vertices)
            if directed:
                h = gt.Graph(directed=True)
            else:
                h = gt.Graph(directed=False)
            h.add_vertex(len(total_vertices))

            vertex_weights = h.new_vertex_property('double')
            vertex_names = h.new_vertex_property('string')
            vertex_colors = h.new_vertex_property('double')
            edge_weights = h.new_edge_property('double')

            # Assign weights and names to vertices (genes and edges)
            for i, vertex in enumerate(total_vertices):
                vertex_weights[i] = edges_counts[vertex] / len(patients) if vertex in target_edges else 0
                vertex_names[i] = vertex
                vertex_colors[i] = 0.5 if vertex in target_edges else 1

            # Add edges between connected edges with weight 1
            for edge in target_edges:
                for edge2 in target_edges:
                    if edge[0] == edge2[1]:
                        e = h.add_edge(total_vertices.index(edge), total_vertices.index(edge2))
                        edge_weights[e] = 1

            # Add edges between genes and edges with weight -log(p-value)
            for driver in total_drivers:
                data_driver = data.loc[data['driver'] == driver]
                for edge in data_driver['dysregulation']:
                    e = h.add_edge(total_vertices.index(edge), total_vertices.index(driver))
                    edge_weights[e] = -np.log(data_driver.loc[data_driver['dysregulation'] == edge, 'p-value'].values[0])

            # Apply PageRank algorithm and store scores
            v_pr = h.new_vertex_property('float')
            pr = gt.centrality.pagerank(h, damping=d, pers=vertex_weights, weight=edge_weights, prop=v_pr)
            pr = pr.get_array()
            top_genes = [vertex_names[i] for i in np.argsort(pr)[:len(total_drivers)]]
            top_scores = [pr[i] for i in np.argsort(pr)[:len(total_drivers)]]
            df = pd.DataFrame({'Gene': top_genes, 'Score': top_scores})

            # Save score files to output directory.
            df.to_csv(output_dir+f'{p}_scores.csv', sep='\t')
        
        # Save per-patient driver-dysregulation file.
        data.to_csv(output_dir+f'{p}_drivers.csv', sep='\t')
    
    else:
        print('No driver genes found for patient ' + p)



def set_up_mutations(mutation_matrix):
    # Reorganize the mutation matrix for quicker access

    patients = mutation_matrix.index.tolist()
    genes = list(mutation_matrix.columns)[1:-1]


    return genes, patients


def set_up_networks(binarized_networks, patients_M):
    edges_counts = binarized_networks.sum().to_dict()
 
    patients_networks = set(binarized_networks.index.values)
    patients = list(patients_networks.intersection(patients_M))
    return edges_counts, patients

def process_patients(patients : list, kappa : int, d : float, scores : bool, 
                     mutations_path : str, networks_path : str, directed : bool,
                     output_dir : str):
    # Check if output directory exists.
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    # Read input data.
    mutations= pd.read_csv(mutations_path,sep=',')
    networks= pd.read_csv(networks_path,sep=',')
    
    # Compute kappa neighborhood of putative drivers.
    tfs = get_relations(kappa, patients, mutations, networks)
    edges_counts = networks['edge'].value_counts().to_dict()
    
    # Create argument list for concurrent task execution.
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-1)
    args_list = [(p, d, scores, mutations, networks, patients, directed, tfs, edges_counts, output_dir) for p in patients]

    # Start parallel analysis.
    pool.starmap(create_statistics_p, args_list)
    pool.close()
    pool.join()

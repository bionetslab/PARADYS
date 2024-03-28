import graph_tool.all as gt
import pandas as pd
import numpy as np
from scipy.stats import chisquare, chi2
import os

def build_dysregulation_graph(networks : dict, patient : str) -> tuple:
    """Build dysregulation network graph object from directed edge list.

    Args:
        networks (dict): Keys are patients, values are sets of edges in tuple format, i.e. ('V1', 'V2').
        patient (str): Patient ID to analyze.

    Returns:
        tuple: First entry is graph tools object, second entry is dict containing
        vertex i's label at index i, third entry is dictionary with keys as vertex labels and
        values as corresponding vertex ID.
    """
    g = gt.Graph(directed=True)
    patient_edges = networks[patient]
    vertex_prop = g.add_edge_list(patient_edges, hashed=True)
    # Turn vertex property into dict.
    vertex_to_label = {ind : vertex_prop[ind] for ind in range(g.num_vertices())}
    # Precompute and return inverse of vertex_to_label to save computation time later.
    label_to_vertex = {vertex_prop[ind] : ind for ind in range(g.num_vertices())}
        
    return g, vertex_to_label, label_to_vertex


def compute_neighbourhood(kappa : int, patient: str, mutations : dict, 
                          tfs : dict, dys_graph : gt.Graph, 
                          vertex_to_label : dict,
                          label_to_vertex : dict) -> dict:
    """Computes shortest-path neighborhood for all mutated TFs.

    Args:
        kappa (int): Size of maximal shortest-path distance.
        patient (str): ID of patient to analyze.
        mutations (dict): Dict storing per-patient mutations.
        tfs (dict): Dict storing per-patient existing transcription factors.
        dys_graph (gt.Graph): Graph tool graph representing patient's dysregulation network.
        vertex_to_label (gt.VertexPropertyMap): At index i, it stores vertex's i label.
        label_to_vertex (dict): Keys are vertex labels, index is vertex ID.
        
    Returns:
        dict: Keys are tuples of mutated TFs and its corresponding ID in the network graph and values 
        are tuples of gene labels and vertex ID in the network within kappa shortest-path neighborhood.
    """
    out_dict = dict()
        
    # Extract patient-specific mutations, tfs, and network.
    patient_tfs = tfs[patient]
    patient_muts = mutations[patient]
    
    # Extract putative drivers (i.e. TFs that are mutated).
    putative_drivers = patient_tfs.intersection(patient_muts)
    
    # If at least one putative driver exists, compute its kappa neighborhood.
    if len(putative_drivers) > 0:
        # Compute all putative driver neighborhoods.
        for driver in putative_drivers:
            driver_id = label_to_vertex[driver]
            
            # Compute shortest-path distances with maximal value of kappa.
            dist = gt.shortest_distance(dys_graph, source=dys_graph.vertex(driver_id), max_dist=kappa)
            np_dist = np.array(dist.get_array())
            
            # Extract all vertices with distances <= kappa.
            nbor_ids = np.where(np_dist <= kappa)[0]
            
            # Turn IDs into vertex string names (store both).
            nbors_set = [(vertex_to_label[x], x) for x in nbor_ids]
            nbors_set = set(nbors_set)
            
            # Append to per-patient dictionary.
            out_dict[(driver, driver_id)]=nbors_set
    
    return out_dict
            
        
def mutations_to_dict(mutations_df : pd.DataFrame) -> dict:
    """Transform mutations DataFrame into dictionary for faster lookup.

    Args:
        mutations_df (pd.DataFrame): Input mutation data (columns are 'patient' and 'gene').

    Returns:
        dict: Keys are patient IDs and values are sets of corresponding mutated genes.
    """
    all_patients = set(mutations_df['patient'])
    mutations_dict = {pat : set() for pat in all_patients}
    
    [mutations_dict[row[0]].add(row[1]) for row in zip(mutations_df['patient'], mutations_df['gene'])]
    
    return mutations_dict

def networks_to_dict(networks_df : pd.DataFrame) -> tuple:
    """Transform networks DataFrame edges and tfs into dictionaries for faster lookup.

    Args:
        networks_df (pd.DataFrame): Input dysregulation edge data (columns are 'patient', 'tf', 'gene').

    Returns:
        tuple: First entry is dict with keys as patient IDs and values are sets of edges tuples
        (in the format (TF, gene)); second entry is dict with keys as patients and values are
        sets of TFs.
    """  
    all_patients = set(networks_df['patient'])
    networks_dict = {pat : set() for pat in all_patients}
    tf_dict = {pat : set() for pat in all_patients}
    
    for row in zip(networks_df['patient'], networks_df['tf'], networks_df['gene']):
        patient = row[0]
        # Extract edge.
        source = row[1]
        target = row[2]
        networks_dict[patient].add((source, target))
        # Extract TF information.
        tf_dict[patient].add(source)
        
    return networks_dict, tf_dict


def compute_mutation_patient_dict(mutations : dict) -> dict:
    """Computes inverse dictionary of mutations dictionary.

    Args:
        mutations (dict): Keys are patients, values are sets of mutated genes.

    Returns:
        dict: Keys are mutated genes, values are sets of patients.
    """
    mut_patient_dict = dict()
    # Iterate over all patients' mutation sets and compute inverse dictionary.
    for pat in mutations.keys():
        patient_mutations = mutations[pat]
        for mut in patient_mutations:
            if not mut in mut_patient_dict.keys():
                mut_patient_dict[mut] = set()
            mut_patient_dict[mut].add(pat)
    
    return mut_patient_dict

def compute_edge_patient_dict(networks : dict) -> dict:
    """Computes inverse dictionary of networks dictionary.

    Args:
        networks (dict): Keys are patients, values are sets of edges.

    Returns:
        dict: Keys are edges, values are sets of patients.
    """
    edge_patient_dict = dict()
    for pat in networks.keys():
        patient_edges = networks[pat]
        for edge in patient_edges:
            if not edge in edge_patient_dict.keys():
                edge_patient_dict[edge] = set()
            edge_patient_dict[edge].add(pat)
    
    return edge_patient_dict

def compute_observation_matrix(edge : tuple, driver : str, mutation_patient_dict : dict,
                               edge_patient_dict : dict, all_patients : set) -> np.ndarray:
    """Compute frequencies of dysregulation and mutation in whole cohort.

    Args:
        edge (tuple): Dysregulation edge in tuple format (source, target).
        driver (str): Name of driver gene.
        mutation_patient_dict (dict): Keys are mutations, values are sets of patients.
        edge_patient_dict (dict): Keys are dysregulation edges, values are sets of patients.
        all_patients (set): Set of all patient IDs in cohort.

    Returns:
        np.ndarray: Observation matrix.
    """
    obs_matrix = np.zeros((2,2), dtype=int)
    
    # Extract sets of patients that have (not) mutation and (not) dysregulation.
    have_dys = edge_patient_dict[edge]
    not_have_dys = all_patients - edge_patient_dict[edge]
    
    have_mut = mutation_patient_dict[driver]
    not_have_mut = all_patients - mutation_patient_dict[driver]
    
    # Entry 0,0 corresponds to number of patients that have both mutation and dysregulation.
    obs_matrix[0,0] = len(have_dys.intersection(have_mut))
    
    # Entry 0,1 corresponds to number of patients that have dysregulation but not mutation.
    obs_matrix[0,1] = len(have_dys.intersection(not_have_mut))
    
    # Entry 1,0 corresponds to number of patients that don't have dysregulation but mutation.
    obs_matrix[1,0] = len(not_have_dys.intersection(have_mut))
    
    # Entry 1,1 corresponds to number of patients that don't have either.
    obs_matrix[1,1] = len(not_have_dys.intersection(not_have_mut))
    
    return obs_matrix


def chi_squared_test(edge : tuple, driver : str, mutation_patient_dict : dict,
                     edge_patient_dict : dict, all_patients : set) -> float:
    """Runs chi-squared test on given driver and dysregulation edge.

    Args:
        edge (tuple): Edge in tuple format (source, target).
        driver (str): Name of driver.
        mutation_patient_dict (dict): Keys are mutations, values are sets of patients.
        edge_patient_dict (dict): Keys are edges, values are sets of patients.
        all_patients (set): Set of all patient IDs in cohort.

    Returns:
        float: P value of chi-squared test.
    """
    # Build observation matrix.
    observation_mat = compute_observation_matrix(edge, driver, mutation_patient_dict, 
                                                 edge_patient_dict, all_patients)
    # Divide observation matrix into rows.
    obs_divided = list(observation_mat)
    
    # Calculate row and column sums.
    row_sums = np.array([np.sum(obs_divided, axis=1)])
    col_sums = np.array([np.sum(obs_divided, axis=0)])
    
    # Create total matrix sum.
    mat_sum = np.sum(obs_divided)
    
    # Build expectancy matrix from row and column sums.
    exp_matrix = np.dot(row_sums.T, col_sums) / mat_sum
    
    # Run chi-squared test.
    chisq, _ = chisquare(obs_divided, exp_matrix)
    chisq = np.sum(chisq)
    pval = 1 - chi2.cdf(chisq, 1)
    return pval    
    

def compute_per_patient_drivers(mutation_patient_dict : dict, edge_patient_dict : dict,
                                dys_graph : gt.Graph, putative_drivers : dict,
                                all_patients : set, vertex_to_label : dict) -> dict:
    """Compute per-patient drivers and dysregulations.

    Args:
        mutation_patient_dict (dict): Keys are mutated genes, values are sets of patients.
        edge_patient_dict (dict): Keys are edges, values are sets of patients.
        dys_graph (gt.Graph): Graph tool graph representing patients dysregulation network.
        putative_drivers (dict): Keys are mutated TFs, values are genes in kappa neighborhood.
        all_patients (set): Set of all patient IDs in cohort.
        vertex_to_label (dict): Dict with vertex i's label at position i.

    Returns:
        dict: Dict of drivers with keys 'driver', 'edge', 'pvalue'.
    """
    results = {key : [] for key in ['driver', 'dysregulation', 'pvalue']}
    
    # Process all putative drivers of given patient.
    for driver, nbors in putative_drivers.items():
        # Build corresponding dysregulation edges in neighborhood of driver.
        for vertex in nbors:
            # Second entry of tuple stores vertex ID in graph object.
            vertex_id = vertex[1]
            # Iterate over neighbors of this vertex.
            vertex_nbors = dys_graph.get_out_neighbors(vertex_id)
            for vertex_nbor in vertex_nbors:
                # First entry of tuple stores vertex name.
                source_node = vertex[0]
                target_node = vertex_to_label[vertex_nbor]
                edge = (source_node, target_node)
                # Test association between driver and dysregulation edge.
                pval = chi_squared_test(edge, driver[0], mutation_patient_dict, 
                                        edge_patient_dict, all_patients)
                # Add driver-dysregulation pair to results if significant.
                if pval < 0.05:
                    results['driver'].append(driver[0])
                    results['dysregulation'].append(edge)
                    results['pvalue'].append(pval)
    
    return results      

def compute_driver_scores(driver_results : pd.DataFrame, d : float, directed : bool) -> dict:
    """Compute personalized ranking of driver genes.

    Args:
        driver_results (pd.DataFrame): Personalized driver results (patient, driver, edge, pvalue).
        d (float): Damping factor for PageRank.
        directed (bool): Whether or not input networks are directed.

    Returns:
        dict: Keys are drivers, values are calculated impact scores.
    """
    drivers = set(driver_results['driver'])
    edges = set(driver_results['dysregulation'])
    nodes = list(drivers.union(edges))
    
    # Init graph object.
    graph = gt.Graph(directed=directed)
    
    # Add nodes to graph.
    vertices = graph.add_vertex(len(nodes))
    
    # Add labels to nodes.
    vertex_names = graph.new_vertex_property("string")
    vertex_weights = graph.new_vertex_property("double")
    edge_weights = graph.new_edge_property("double")
    
    # Set names and weights on vertices.
    for v, node in zip(vertices, nodes):
        vertex_names[v] = str(node)
        # TODO: update vertex weight of edges accordingly.
        vertex_weights[v] = len(driver_results[driver_results['driver']==node])/len(driver_results)
        
    # Add edges between connected dysregulation edges with weight 1.
    for source in edges:
        for target in edges:
            if source[0]==target[1]:
                e = graph.add_edge(nodes.index(source), nodes.index(target))
                edge_weights[e]=1
    
    # Add edges between genes and edges with weight -log(pval).
    for driver in drivers:
        # Extract corresponding dysregulation edges to driver.
        driver_edges = driver_results[driver_results['driver']==driver]
        for _, row in driver_edges.iterrows():
            edge = row['dysregulation']
            e = graph.add_edge(nodes.index(edge), nodes.index(driver))
            edge_weights[e]=-np.log(row['pvalue'])
    
    # Apply PageRank algorithm using vertex and edge weights.
    pagerank_prop = graph.new_vertex_property("float")
    pagerank = gt.pagerank(graph, damping=d, pers=vertex_weights,
                                      weight=edge_weights, prop=pagerank_prop)
    pagerank_res = pagerank.get_array()
    driver_names = {str(x) for x in drivers}
    sorted_drivers = [vertex_names[i] for i in np.argsort(pagerank_res) if vertex_names[i] in driver_names]
    sorted_scores = [pagerank_res[i] for i in np.argsort(pagerank_res) if vertex_names[i] in driver_names]
    
    out_dict = {'driver': sorted_drivers , 'score': sorted_scores}
    return out_dict 


def process_patients(patients : list, kappa : int, d : float, scores : bool, 
                     mutations_path : str, networks_path : str, directed : bool,
                     output_dir : str):
    """Main function for starting PARADYS analysis.

    Args:
        patients (list): List of patient ID to be analyzed.
        kappa (int): Size of shortest-path neighborhood for putative driver identification.
        d (float): Damping factor for PageRank algorithm.
        scores (bool): Whether to compute personalized driver impact scores.
        mutations_path (str): Path to input mutations file.
        networks_path (str): Path to input networks file.
        directed (bool): Whether dysregulation networks are directed.
        output_dir (str): Path to output directory for storing results.
    """
    
    # Check if output directory exists.
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    # Read input data.
    mutations_df = pd.read_csv(mutations_path,sep=',')
    networks_df = pd.read_csv(networks_path,sep=',')
    
    # Transform mutation and network data into dictionaries for faster lookup.
    mutations = mutations_to_dict(mutations_df)
    networks, tfs = networks_to_dict(networks_df)
    
    del mutations_df
    del networks_df
    
    # Init set of all patient IDs in cohort.
    all_patients = set(mutations.keys())
    assert all_patients==set(networks.keys()), "Mutation and network patient lists are not equal."
    
    # Precompute sets of patients that contain given mutation or dysregulation edge.
    # Saves computation time in later per-patient chi squared tests.
    mutation_patient_dict = compute_mutation_patient_dict(mutations)
    edge_patient_dict = compute_edge_patient_dict(networks)
    
    # Analyze each patient individually.
    for pat in patients:
        # Build dysregulation network graph object.
        dys_graph, vertex_to_label, label_to_vertex = build_dysregulation_graph(networks, pat)
        
        # Compute kappa neighborhood of putative drivers.
        putative_drivers = compute_neighbourhood(kappa, pat, mutations, tfs, 
                                                 dys_graph, vertex_to_label, label_to_vertex)
        # Compute significant driver-dysregulation pairs.
        driver_results = compute_per_patient_drivers(mutation_patient_dict, edge_patient_dict, 
                                    dys_graph, putative_drivers, all_patients, vertex_to_label)
        # Append drivers to final output dataframe.
        driver_df = pd.DataFrame(driver_results)
        driver_df.to_csv(output_dir+f'{pat}_drivers.csv', sep='\t', index=False)
        
        # Check if personalized ranking of drivers is desired.
        if scores:
            driver_scores = compute_driver_scores(driver_df, d, directed)
            scores_df = pd.DataFrame(driver_scores)
            # Save in personalized scores file.
            scores_df.to_csv(output_dir+f'{pat}_scores.csv', sep='\t', index=False)

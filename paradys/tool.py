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
        tuple: First entry is graph tools object, second entry is VertexPropertyMap containing
        vertex names to vertex ID translation.
    """
    g = gt.Graph(directed=True)
    patient_edges = networks[patient]
    vertex_names = g.add_edge_list(patient_edges, hashed=True)
    return g, vertex_names
    
def get_vertex_id(driver : str, vertex_labels : gt.VertexPropertyMap, max_id : int) -> int:
    """Finds corresponding ID of string vertex label in graph.

    Args:
        driver (str): String label of driver vertex.
        vertex_labels (gt.VertexPropertyMap): Array with vertex i's node label at index i.
        max_id (int): Length of vertex_labels, i.e. number of vertices.

    Returns:
        int: Vertex ID in graph object.
    """
    for vert in range(max_id):
        if vertex_labels[vert]==driver:
            return vert
    
    print("Warning: Index of given vertex name could not be found.")
    return -1

def compute_neighbourhood(kappa : int, patient: str, mutations : dict, 
                          tfs : dict, dys_graph : gt.Graph, 
                          vertex_translator : gt.VertexPropertyMap) -> dict:
    """Computes shortest-path neighborhood for all mutated TFs.

    Args:
        kappa (int): Size of maximal shortest-path distance.
        patient (str): ID of patient to analyze.
        mutations (dict): Dict storing per-patient mutations.
        tfs (dict): Dict storing per-patient existing transcription factors.
        dys_graph (gt.Graph): Graph tool graph representing patient's dysregulation network.
        vertex_translator (gt.VertexPropertyMap): At index i, it stores vertex's i label.
        
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
            driver_id = get_vertex_id(driver, vertex_translator, dys_graph.num_vertices())
            
            # Compute shortest-path distances with maximal value of kappa.
            dist = gt.shortest_distance(dys_graph, source=dys_graph.vertex(driver_id), max_dist=kappa)
            np_dist = np.array(dist.get_array())
            
            # Extract all vertices with distances <= kappa.
            nbor_ids = np.where(np_dist <= kappa)[0]
            
            # Turn IDs into vertex string names (store both).
            nbors_set = [(vertex_translator[x], x) for x in nbor_ids]
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
    mutations_dict = dict()
    all_patients = set(mutations_df['patient'])
    
    for patient in all_patients:
        mutations_patient = mutations_df[mutations_df['patient']==patient]
        genes = set(mutations_patient['gene'])
        mutations_dict[patient]=genes
    
    return mutations_dict

def networks_to_dict(networks_df : pd.DataFrame) -> dict:
    """Transform networks DataFrame into dictionary for faster lookup.

    Args:
        networks_df (pd.DataFrame): Input dysregulation edge data (columns are 'patient', 'tf', 'gene').

    Returns:
        dict: Keys are patient IDs and values are sets of edges tuples (in the format (TF, gene)).
    """  
    networks_dict = dict()
    all_patients = set(networks_df['patient'])
    
    for patient in all_patients:
        networks_patients = networks_df[networks_df['patient']==patient]
        edges = set()
        
        for _, row in networks_patients.iterrows():
            source = row['tf']
            target = row['gene']
            edges.add((source, target))
        
        networks_dict[patient] = edges
        
    return networks_dict

def extract_tfs(networks_df : pd.DataFrame) -> dict:
    """Extract set of transcription factors existing in per-patient network.

    Args:
        networks_df (pd.DataFrame):Input dysregulation edge data (columns are 'patient', 'tf', 'gene').

    Returns:
        dict: Keys are patients, values are sets of transcription factors.
    """
    all_patients = set(networks_df['patient'])
    tf_dict = dict()
    
    for patient in all_patients:
        networks_patient = networks_df[networks_df['patient']==patient]
        tfs = set(networks_patient['tf'])
        tf_dict[patient]=tfs
    
    return tf_dict   


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
                                all_patients : set, vertex_translator : gt.VertexPropertyMap) -> dict:
    """Compute per-patient drivers and dysregulations.

    Args:
        mutation_patient_dict (dict): Keys are mutated genes, values are sets of patients.
        edge_patient_dict (dict): Keys are edges, values are sets of patients.
        dys_graph (gt.Graph): Graph tool graph representing patients dysregulation network.
        putative_drivers (dict): Keys are mutated TFs, values are genes in kappa neighborhood.
        all_patients (set): Set of all patient IDs in cohort.
        vertex_translator (gt.VertexPropertyMap): Array with vertex i's label at index i.

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
                target_node = vertex_translator[vertex_nbor]
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
    networks = networks_to_dict(networks_df)
    tfs = extract_tfs(networks_df)
    
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
    out_data_frame = pd.DataFrame(columns=['patient', 'driver', 'dysregulation', 'pvalue'])
    for pat in patients:
        # Build dysregulation network graph object.
        dys_graph, vertex_translator = build_dysregulation_graph(networks, pat)
        
        # Compute kappa neighborhood of putative drivers.
        putative_drivers = compute_neighbourhood(kappa, pat, mutations, tfs, 
                                                 dys_graph, vertex_translator)
        # Compute significant driver-dysregulation pairs.
        driver_results = compute_per_patient_drivers(mutation_patient_dict, edge_patient_dict, 
                                    dys_graph, putative_drivers, all_patients, vertex_translator)
        # Append drivers to final output dataframe.
        driver_results['patient']=[pat]*(len(driver_results['driver']))
        driver_df = pd.DataFrame(driver_results)
        out_data_frame = pd.concat([out_data_frame, driver_df], ignore_index=True)
    
    # Write output dataframe for all patients.
    out_data_frame.to_csv(output_dir+'drivers.csv', sep='\t', index=False)
    
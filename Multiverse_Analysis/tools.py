import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import euclidean

# Custom Functions
def reshape_data(df, variable):
    """
    Reshape the data to wide format with each time point as a column (timeseries)
    
    Parameters:
    - df (dataframe): Input DataFrame containing the data to reshape
    - variable (str): The variable name to reshape into a wide format (timeseries)
    
    Returns:
    - reshaped_df: Reshaped dataframe, timeseries
    """
    reshaped_df = df.pivot(index='tempid', columns='time', values=variable)
    reshaped_df.columns = [f'time{col}' for col in reshaped_df.columns]
    return reshaped_df

def trans_mat(int_list):
    """
    Create a transition matrix from timeseries data
    
    Parameters:
    - int_list (list of integers): timseries data
    
    Returns:
    - transition_matrix: Transition matrix (array) where each cell (i, j) represents the 
                         probability of transitioning from state i+1 to state j+1
    """
    states = [1, 2, 3, 4, 5]
    n_states = 5
    transition_matrix = np.zeros((n_states, n_states), dtype=float)
    
    for i in range(1, len(int_list)):
        prev_state = int_list[i - 1] - 1
        next_state = int_list[i] - 1
        transition_matrix[prev_state, next_state] += 1
    
    for i in range(n_states):
        row_sum = transition_matrix[i, :].sum()
        if row_sum > 0:
            transition_matrix[i, :] /= row_sum

    return transition_matrix

def prop_trans_mat(int_list):
    """
    Create a proportional transition matrix from timeseries data
    
    Parameters:
    - int_list (list of integers): timeseries data
    
    Returns:
    - prop_transition_matrix: Proportional transition matrix where each cell (i, j) represents 
                              the proportion of transitions relative to the total number of transitions
    """
    states = [1, 2, 3, 4, 5]
    n_states = 5
    prop_transition_matrix = np.zeros((n_states, n_states), dtype=float)
    
    for i in range(1, len(int_list)):
        prev_state = int_list[i - 1] - 1
        next_state = int_list[i] - 1
        prop_transition_matrix[prev_state, next_state] += 1
    
    total_transitions = prop_transition_matrix.sum()
    if total_transitions > 0:
        prop_transition_matrix /= total_transitions

    return prop_transition_matrix

def compute_rsm_pearson(matrices):
    """
    Compute a similarity matrix using Pearson correlation on flattened matrices
    
    Parameters:
    - matrices (list of np.ndarray): List of transition matrices from each subject
    
    Returns:
    - similarity_matrix: Symmetric similarity matrix where each cell (i, j) represents the Pearson 
                         correlation between matrices i and j (i and j are subject numbers)
    """
    n_subjects = len(matrices)
    flattened_matrices = np.array([matrices[i].flatten() for i in range(n_subjects)])
    similarity_matrix = np.zeros((n_subjects, n_subjects))
    
    for i in range(n_subjects):
        for j in range(i, n_subjects):
            if np.all(flattened_matrices[i] == flattened_matrices[i][0]) or np.all(flattened_matrices[j] == flattened_matrices[j][0]):
                similarity_matrix[i, j] = np.nan
                similarity_matrix[j, i] = np.nan
            else:
                corr, _ = pearsonr(flattened_matrices[i], flattened_matrices[j])
                similarity_matrix[i, j] = corr
                similarity_matrix[j, i] = corr

    return similarity_matrix

def compute_rsm_spear(matrices):
    n_subjects = len(matrices)
    flattened_matrices = np.array([mat.flatten() for mat in matrices])
    similarity_matrix = np.zeros((n_subjects, n_subjects))
    
    for i in range(n_subjects):
        for j in range(i, n_subjects):
            try:
                if np.all(flattened_matrices[i] == flattened_matrices[i][0]) or np.all(flattened_matrices[j] == flattened_matrices[j][0]):
                    similarity_matrix[i, j] = np.nan
                    similarity_matrix[j, i] = np.nan
                else:
                    corr, _ = spearmanr(flattened_matrices[i], flattened_matrices[j])
                    similarity_matrix[i, j] = corr
                    similarity_matrix[j, i] = corr
            except Exception as e:
                print(f"Error at i={i}, j={j}: {e}")
                similarity_matrix[i, j] = np.nan
                similarity_matrix[j, i] = np.nan
    
    return similarity_matrix

def compute_rsm_euclid(matrices):
    """
    Compute a similarity matrix using Euclidean distance on flattened matrices
    
    Parameters:
    - matrices (list of np.ndarray): List of transition matrices from each subject
    
    Returns:
    - similarity_matrix: Symmetric similarity matrix where each cell (i, j) represents the Euclidean 
                         distance between matrices i and j (i and j are subject numbers)
    """
    n_subjects = len(matrices)
    flattened_matrices = np.array([matrices[i].flatten() for i in range(n_subjects)])
    similarity_matrix = np.zeros((n_subjects, n_subjects))
    
    for i in range(n_subjects):
        for j in range(i, n_subjects):
            if np.all(flattened_matrices[i] == flattened_matrices[i][0]) or np.all(flattened_matrices[j] == flattened_matrices[j][0]):
                similarity_matrix[i, j] = np.nan
                similarity_matrix[j, i] = np.nan
            else:
                corr = euclidean(flattened_matrices[i], flattened_matrices[j])
                similarity_matrix[i, j] = corr
                similarity_matrix[j, i] = corr

    return similarity_matrix

def calculate_entropy_per_subject(transition_matrices):
    """
    Calculate entropy for each subject based on their transition matrix (use with LOSO approaches)
    
    Parameters:
    - transition_matrices (list of np.ndarray): List of transition matrices for each subject
    
    Returns:
    - list of float:Entropy values for each subject
    """
    entropy_values = []
    for matrix in transition_matrices:
        flat_matrix = matrix.flatten()
        non_zero_probs = flat_matrix[flat_matrix > 0]
        entropy = -np.sum(non_zero_probs * np.log(non_zero_probs))
        entropy_values.append(entropy)
    
    return entropy_values

def calculate_entropy_differences(transition_matrices):
    """
    Compute absolute differences in entropy between each pair of subjects
    
    Parameters:
    - transition_matrices (list of np.ndarray): List of transition matrices for each subject
    
    Returns:
    - np.ndarray: Symmetric matrix of entropy differences
    """
    entropy_values = calculate_entropy_per_subject(transition_matrices)
    n_subjects = len(entropy_values)
    entropy_difference_matrix = np.zeros((n_subjects, n_subjects))
    
    for i in range(n_subjects):
        for j in range(i, n_subjects):
            diff = abs(entropy_values[i] - entropy_values[j])
            entropy_difference_matrix[i, j] = diff
            entropy_difference_matrix[j, i] = diff
    
    return entropy_difference_matrix

def indiv_ts_tmat(timeseries):
    """
    Create individual subject time series and transition matrices from a timeseries dataframe
    
    Parameters:
    - timeseries (dataframe): Input DataFrame with time series data indexed by subject
    
    Returns:
    - list of list: Individual time series for each subject
    - list of np.ndarray: Transition matrices for each subject
    """
    row_indexes = timeseries.index.tolist()
    timeseries_list = []
    transition_matrices = []
    for tempid in row_indexes:
        temp = timeseries.loc[tempid].dropna().astype(int).tolist()
        timeseries_list.append(temp)
        transmat = trans_mat(temp)
        transition_matrices.append(transmat)
    return timeseries_list, transition_matrices

def indiv_ts_pmat(timeseries):
    """
    Create individual time series and proportional transition matrices from a timeseries dataframe
    
    Parameters:
    - timeseries (dataframe): Input DataFrame with time series data indexed by subject
    
    Returns:
    - list of list: Individual time series for each subject
    - list of np.ndarray: Proportional transition matrices for each subject
    """
    row_indexes = timeseries.index.tolist()
    timeseries_list = []
    prop_tran_matrices = []
    for tempid in row_indexes:
        temp = timeseries.loc[tempid].dropna().astype(int).tolist()
        timeseries_list.append(temp)
        transmat = prop_trans_mat(temp)
        prop_tran_matrices.append(transmat)
    return timeseries_list, prop_tran_matrices

def loso_similarity_matrix_pearson(transition_matrix_list):
    """
    Compute Leave-One-Subject-Out (LOSO) similarity matrix using Pearson correlation
    
    Parameters:
    - transition_matrix_list: List of transition matrices for all subjects
    
    Returns:
    - results: List of Pearson correlations for each subject's matrix compared to the average of others
    """
    n_subjects = len(transition_matrix_list)
    transition_matrix_list = np.array(transition_matrix_list)
    results = []

    for i in range(n_subjects):
        avg_transition_matrix = np.mean(np.delete(transition_matrix_list, i, axis=0), axis=0)
        subject_matrix = transition_matrix_list[i].flatten()
        avg_matrix = avg_transition_matrix.flatten()
        correlation, _ = pearsonr(subject_matrix, avg_matrix)
        results.append(correlation)

    return results

def loso_similarity_matrix_spear(transition_matrix_list):
    """
    Compute Leave-One-Subject-Out (LOSO) similarity matrix using Spearman correlation
    
    Parameters:
    - transition_matrix_list: List of transition matrices for all subjects
    
    Returns:
    - results: List of Spearman correlations for each subject's matrix compared to the average of others
    """
    n_subjects = len(transition_matrix_list)
    transition_matrix_list = np.array(transition_matrix_list)
    results = []

    for i in range(n_subjects):
        avg_transition_matrix = np.mean(np.delete(transition_matrix_list, i, axis=0), axis=0)
        subject_matrix = transition_matrix_list[i].flatten()
        avg_matrix = avg_transition_matrix.flatten()
        correlation, _ = spearmanr(subject_matrix, avg_matrix)
        results.append(correlation)

    return results

def loso_similarity_matrix_euclid(transition_matrix_list):
    """
    Compute Leave-One-Subject-Out (LOSO) similarity matrix using Euclidean distance
    
    Parameters:
    - transition_matrix_list: List of transition matrices for all subjects
    
    Returns:
    - results: List of Euclidean distances for each subject's matrix compared to the average of others
    """
    n_subjects = len(transition_matrix_list)
    transition_matrix_list = np.array(transition_matrix_list)
    results = []

    for i in range(n_subjects):
        avg_transition_matrix = np.mean(np.delete(transition_matrix_list, i, axis=0), axis=0)
        subject_matrix = transition_matrix_list[i].flatten()
        avg_matrix = avg_transition_matrix.flatten()
        distance = euclidean(subject_matrix, avg_matrix)
        results.append(distance)

    return results

def reshape_data_to_daily_averages(df, variable, start_date_column='StartDate'):
    """
    Reshape the data to daily averages, with days as columns
    
    Parameters:
    - df (dataframe): Input dataframe
    - variable (str): The variable to calculate daily averages for.
    - start_date_column (str): Column name containing start dates. Default is 'StartDate'
    
    Returns:
    - reshaped_df: Reshaped dataframe with days as columns and average values for the variable
    """
    df[start_date_column] = pd.to_datetime(df[start_date_column])
    df['date'] = df[start_date_column].dt.date
    daily_avg = df.groupby(['tempid', 'date'])[variable].mean().reset_index()
    daily_avg['day'] = daily_avg.groupby('tempid')['date'].rank().astype(int)
    reshaped_df = daily_avg.pivot(index='tempid', columns='day', values=variable)
    reshaped_df.columns = [f'Day{col}' for col in reshaped_df.columns]
    
    return reshaped_df

def calculate_date_differences(start_days_list):
    """
    Calculate normalized date differences between all pairs of subjects
    
    Parameters:
    - start_days_list: List of start dates for each subject.
    
    Returns:
    - data_difference_matrix: Symmetric matrix of date differences between subjects
    """
    ordinal_days = [day.toordinal() for day in start_days_list]
    n_subjects = len(ordinal_days)
    date_difference_matrix = np.zeros((n_subjects, n_subjects))
    
    for i in range(n_subjects):
        for j in range(i, n_subjects):
            diff = abs(ordinal_days[i] - ordinal_days[j])
            date_difference_matrix[i, j] = diff
            date_difference_matrix[j, i] = diff

    return date_difference_matrix

def calculate_date_differences_loso(start_days_list):
    """
    Calculate normalized Leave-One-Subject-Out (LOSO) date differences
    
    Parameters:
    - start_days_list: List of start dates for each subject
    
    Returns:
    - loso_differences: LOSO differences for each subject
    """
    n_subjects = len(start_days_list)
    ordinal_days = [day.toordinal() for day in start_days_list]
    loso_differences = []

    for i in range(n_subjects):
        others = ordinal_days[:i] + ordinal_days[i + 1:]
        avg_others = np.mean(others)
        diff = abs(ordinal_days[i] - avg_others)
        loso_differences.append(diff)

    return loso_differences

def calculate_response_windows(response_data):
    """
    Calculate a Euclidean distance matrix based on response windows
    
    Parameters:
    - response_data: Input data containing response data and time information
    
    Returns:
    - distance_matrix: Symmetric matrix of Euclidean distances based on response patterns
    """
    response_data['TimeSlot'] = response_data['StartDate'].dt.strftime('%H:%M')
    time_slots = [f"{hour:02}:{minute:02}" for hour in range(10, 22) for minute in (0, 30)]
    subject_time_tally = response_data.groupby(['tempid', 'TimeSlot']).size().unstack(fill_value=0)
    subject_time_tally = subject_time_tally.reindex(columns=time_slots, fill_value=0)
    tally_array = subject_time_tally.values
    n_subjects = tally_array.shape[0]
    distance_matrix = np.zeros((n_subjects, n_subjects))
    
    for i in range(n_subjects):
        for j in range(i, n_subjects):
            dist = euclidean(tally_array[i], tally_array[j])
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist

    return distance_matrix

def calculate_response_windows_loso(response_data):
    """
    Calculate Leave-One-Subject-Out (LOSO) Euclidean distances for response windows
    
    Parameters:
    - response_data (dataframe): Input data containing response data and time information
    
    Returns:
    - loso_distances: Euclidean distances for each subject compared to the average of others
    """
    response_data['TimeSlot'] = response_data['StartDate'].dt.strftime('%H:%M')
    time_slots = [f"{hour:02}:{minute:02}" for hour in range(10, 22) for minute in (0, 30)]
    subject_time_tally = response_data.groupby(['tempid', 'TimeSlot']).size().unstack(fill_value=0)
    subject_time_tally = subject_time_tally.reindex(columns=time_slots, fill_value=0)
    tally_array = subject_time_tally.values
    n_subjects = tally_array.shape[0]
    loso_distances = []

    for i in range(n_subjects):
        others = np.delete(tally_array, i, axis=0)
        avg_others = np.mean(others, axis=0)
        dist_temp = euclidean(tally_array[i], avg_others)
        dist = dist_temp
        loso_distances.append(dist)

    return loso_distances

def loso_strat_euclidean(data, variables):
    """
    Calculate Leave-One-Subject-Out (LOSO) Euclidean distances for strategy variables
    
    Parameters:
    - data (dataframe): Input containing subject data.
    - variables (list of str): List of column names representing variables (strategies) to include in the calculation
    
    Returns:
    - results: Array of Euclidean distances for each subject compared to the average of others
    """
    results = []
    subject_ids = data['tempid'].values
    n_subjects = len(subject_ids)
    
    for i in range(n_subjects):
        current_subject = data.iloc[i]
        other_subjects = data[data['tempid'] != current_subject['tempid']][variables]
        mean_other = other_subjects.mean(axis=0)
        distance_temp = euclidean(current_subject[variables], mean_other)
        distance = distance_temp
        results.append(distance)
    
    return np.array(results)

def rsm_strat_euclidean(data, variables):
    """
    Compute a similarity matrix using Euclidean distances for strategy variables
    
    Parameters:
    - data (dataframe): Input containing subject data.
    - variables (list of str): List of column names representing variables (strategies) to include in the calculation.
    
    Returns:
    - np.ndarray: Symmetric matrix of Euclidean distances between subjects
    """
    subject_ids = data['tempid'].values
    n_subjects = len(subject_ids)
    rsm = np.zeros((n_subjects, n_subjects))
    
    for i in range(n_subjects):
        for j in range(n_subjects):
            distance = euclidean(data.loc[i, variables], data.loc[j, variables])
            rsm[i, j] = distance
    
    return rsm

def binarize_responses(
    sex, 
    gender, 
    eth, 
    race_1_black, 
    race_2_asian, 
    race_3_nathi_pacisl, 
    race_4_white, 
    race_5_amin_alnat, 
    race_6_other, 
    par_fin_help
):
    """
    Binarize demographic responses
    
    Parameters:
    - sex (int): 0 for Male, 1 for Female
    - gender (int): Integer representing gender identity (1-5)
    - eth (int): 1 for Yes, 2 for No, 3 for Prefer not to answer
    - race (multiple variables): 1 or 0 variables for each race asked about
    - par_fin_help (int): 0 for No, 1 for Yes
    
    Returns:
    - np.ndarray: Binarized response vector.
    """
    sex_bin = [0, 1] if sex == 1 else [1, 0]
    gender_bin = [1 if gender == i else 0 for i in range(1, 6)]
    hispanic_bin = [1 if eth == 1 else 0]
    race_bin = [race_1_black, race_2_asian, race_3_nathi_pacisl, race_4_white, race_5_amin_alnat, race_6_other]
    financial_bin = [1 if par_fin_help else 0]
    response_vector = (
        sex_bin +
        gender_bin +
        hispanic_bin +
        race_bin +
        financial_bin
    )

    return np.array(response_vector)

def process_row(row):
    """
    Process a single row of demographic data into a binarized response vector
    
    Parameters:
    - row: Row of demographic data
    
    Returns:
    - Binarized response vector for the given row
    """
    return binarize_responses(
        row['sex'],
        row['gender'],
        row['eth'],  
        row['race_1_black'],  
        row['race_2_asian'],
        row['race_3_nathi_pacisl'],
        row['race_4_white'],
        row['race_5_amin_alnat'],
        row['race_6_other'],
        row['par_fin_help']
    )


def loso_demo_euclidean(data, demographic_variables):
    """
    Calculate Leave-One-Subject-Out (LOSO) Euclidean distances for demographic variables
    
    Parameters:
    - data (dataframe): Input containing demographic data
    - demographic_variables (list of str): List of column names representing demographic variables
    
    Returns:
    - results: Array of Euclidean distances for each subject compared to the average of others
    """
    results = []
    subject_ids = data['tempid'].values
    n_subjects = len(subject_ids)

    for i in range(n_subjects):
        current_subject = data.iloc[i]
        current_vector = process_row(current_subject)

        other_subjects = data[data['tempid'] != current_subject['tempid']]
        other_vectors = other_subjects.apply(process_row, axis=1)
        mean_other = np.mean(np.stack(other_vectors), axis=0)

        distance_temp = euclidean(current_vector, mean_other)
        distance = distance_temp
        results.append(distance)

    return np.array(results)

def rsm_demo_euclidean(data, demographic_variables):
    """
    Compute a similarity matrix using Euclidean distances for demographic variables
    
    Parameters:
    - data  Input containing demographic data.
    - demographic_variables (list of str): List of column names representing demographic variables
    
    Returns:
    - results: Symmetric matrix of Euclidean distances between subjects based on demographic variables
    """
    subject_ids = data['tempid'].values
    n_subjects = len(subject_ids)
    rsm = np.zeros((n_subjects, n_subjects))

    binarized_vectors = data.apply(process_row, axis=1)

    for i in range(n_subjects):
        for j in range(n_subjects):
            distance = euclidean(binarized_vectors[i], binarized_vectors[j])
            rsm[i, j] = distance

    return rsm
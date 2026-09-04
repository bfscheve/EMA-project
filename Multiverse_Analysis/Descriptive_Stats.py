# Descriptive Statistics for our EMA data
import pandas as pd
import numpy as np

# CHANGE INPUT PATHS HERE
EMA_path = "Z:\Projects\EMA_Project\Data\EMA_UCLA_Data"
EMA_data_file = "ERT_EMA_DailySurvey.csv"
Data_path = f"{EMA_path}/{EMA_data_file}"
Intake_data_file = "ERT_EMA_Intake.csv"
Data_path2 = f"{EMA_path}/{Intake_data_file}"

# Import data
EMA_data_desc = pd.read_csv(Data_path)
EMA_data_desc.rename(columns={'Duration (in seconds)': 'duration'}, inplace=True)
EMA_data_desc['StartDate'] = pd.to_datetime(EMA_data_desc['StartDate'], format="%m/%d/%Y %H:%M")
Intake_data_desc = pd.read_csv(Data_path2)

# Count the number of times each tempid appears
tempid_counts = EMA_data_desc['tempid'].value_counts()
n_participants = tempid_counts.size
print(f"Total number of participants: {n_participants}")
print(f"Total number of EMA observations: {tempid_counts.sum()}")

# Count participants with more than 60 and more than 75 appearances
count_more_than_60 = (tempid_counts >= 60).sum()
count_more_than_75 = (tempid_counts >= 75).sum()

print(f"Number of participants with 60 or more EMAs: {count_more_than_60}")
print(f"Number of participants with 75 or more EMAs: {count_more_than_75}")
mean_emas = tempid_counts.mean()
print(f"Mean number of EMAs per participant: {mean_emas:.3f}")

# Calculate mean duration of EMAs
mean_duration = EMA_data_desc['duration'].mean()
print(f"Mean EMA duration (seconds): {mean_duration:.3f}")

# List of variables for which mean and standard deviation are to be calculated
variables = ['posAff', 'negAff', 'stress', 'selfER_1', 'selfER_2', 'selfER_3', 'selfER_4', 'selfER_5']

# Loop through each variable and calculate mean and standard deviation
for var in variables:
    mean_value = EMA_data_desc[var].mean()
    std_value = EMA_data_desc[var].std()
    print(f"{var}: Mean = {mean_value:.3f}, SD = {std_value:.3f}")

# Temporal range of dataset
print(f"Earliest survey StartDate: {EMA_data_desc['StartDate'].min()}")
print(f"Latest survey EndDate: {EMA_data_desc['EndDate'].max()}")

# Sex
n_female = np.sum(Intake_data_desc['sex'])
print(f"Number of female participants: {n_female}")
print(f"Percent female: {n_female/114:.3%}")

# Race/ethnicity
print("Race/Ethnicity breakdown:")
n_black = np.sum(Intake_data_desc['race_1_black'])
print(f"  Black: {n_black} ({n_black/114:.3%})")
n_asian = np.sum(Intake_data_desc['race_2_asian'])
print(f"  Asian: {n_asian} ({n_asian/114:.3%})")
n_nathi_pacisl = np.sum(Intake_data_desc['race_3_nathi_pacisl'])
print(f"  Native Hawaiian/Pacific Islander: {n_nathi_pacisl} ({n_nathi_pacisl/114:.3%})")
n_white = np.sum(Intake_data_desc['race_4_white'])
print(f"  White: {n_white} ({n_white/114:.3%})")
n_amin_alnat = np.sum(Intake_data_desc['race_5_amin_alnat'])
print(f"  American Indian/Alaska Native: {n_amin_alnat} ({n_amin_alnat/114:.3%})")
n_other = np.sum(Intake_data_desc['race_6_other'])
print(f"  Other: {n_other} ({n_other/114:.3%})")
n_eth = np.sum(Intake_data_desc['eth'])
print(f"  Hispanic/Latinx: {n_eth} ({n_eth/114:.3%})")
# Multiracial calculation
race_cols = [
    'race_1_black', 'race_2_asian', 'race_3_nathi_pacisl',
    'race_4_white', 'race_5_amin_alnat', 'race_6_other'
]
multiracial = Intake_data_desc[race_cols].sum(axis=1) > 1
n_multiracial = multiracial.sum()
print(f"  Multiracial: {n_multiracial} ({n_multiracial/114:.3%})")
declined = (Intake_data_desc[race_cols].sum(axis=1) == 0).sum()
print(f"  Declined to state race: {declined} ({declined/114:.3%})")

# Age
mean_age = np.mean(Intake_data_desc['age'])
std_age = np.std(Intake_data_desc['age'], ddof=1)
print(f"Mean age: {mean_age:.3f}")
print(f"Age standard deviation: {std_age:.3f}")

# EERQ
print("EERQ Subscale Scores:")
# Dictionary mapping subscales to their respective columns
erq_mapping = {
    'RP': ['ERQ_RP_01', 'ERQ_RP_02', 'ERQ_RP_03', 'ERQ_RP_04', 'ERQ_RP_05', 'ERQ_RP_06'],
    'DS': ['ERQ_DS_01', 'ERQ_DS_02', 'ERQ_DS_03', 'ERQ_DS_04', 'ERQ_DS_05'],
    'SP': ['ERQ_SP_01', 'ERQ_SP_02', 'ERQ_SP_03', 'ERQ_SP_04'],
    'SA': ['ERQ_SA_01', 'ERQ_SA_02', 'ERQ_SA_03', 'ERQ_SA_04'],
    'SS': ['ERQ_SS_01', 'ERQ_SS_02', 'ERQ_SS_03']
}

# Create new columns for each subscale sum and calculate mean and std
for subscale, columns in erq_mapping.items():
    # Sum for each subscale
    Intake_data_desc[subscale + '_sum'] = Intake_data_desc[columns].sum(axis=1)
    
    # Mean and standard deviation for each subscale
    mean_value = Intake_data_desc[subscale + '_sum'].mean()
    std_value = Intake_data_desc[subscale + '_sum'].std()
    
    # Print the results
    print(f"  {subscale}: Mean = {mean_value:.3f}, SD = {std_value:.3f}")
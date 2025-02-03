import psycopg2
import csv

# Connect to PostgreSQL database
connection = psycopg2.connect(
    dbname="chembl_alessandra",
    user="postgres",
    password=" Sn96(0O£TL)@",
    host="localhost",
    port="5432"
)
cursor = connection.cursor()

# Create table schema
create_table_query = '''
CREATE TABLE IF NOT EXISTS chemicals (
    SMILES TEXT,
    Mutagenicity_Ames_test_CONSENSUS_model_assessment TEXT,
    Mutagenicity_Ames_test_CONSENSUS_model_prediction TEXT,
    Mutagenicity_Ames_test_model_CAESAR_assessment TEXT,
    Mutagenicity_Ames_test_model_CAESAR_prediction TEXT,
    Mutagenicity_Ames_test_model_ISS_assessment TEXT,
    Mutagenicity_Ames_test_model_ISS_prediction TEXT,
    Mutagenicity_Ames_test_model_SarPy_IRFMN_assessment TEXT,
    Mutagenicity_Ames_test_model_SarPy_IRFMN_prediction TEXT,
    Mutagenicity_Ames_test_model_KNN_Read_Across_assessment TEXT,
    Mutagenicity_Ames_test_model_KNN_Read_Across_prediction TEXT,
    Mutagenicity_Ames_test_model_aromatic_amines_CONCERT_IRFMN_assessment TEXT,
    Mutagenicity_Ames_test_model_aromatic_amines_CONCERT_IRFMN_prediction TEXT,
    Developmental_Toxicity_model_CAESAR_assessment TEXT,
    Developmental_Toxicity_model_CAESAR_prediction TEXT,
    Developmental_Reproductive_Toxicity_library_PG_assessment TEXT,
    Developmental_Reproductive_Toxicity_library_PG_prediction TEXT,
    Carcinogenicity_model_CAESAR_assessment TEXT,
    Carcinogenicity_model_CAESAR_prediction TEXT,
    Carcinogenicity_model_ISS_assessment TEXT,
    Carcinogenicity_model_ISS_prediction TEXT,
    Carcinogenicity_model_IRFMN_ISSCAN_CGX_assessment TEXT,
    Carcinogenicity_model_IRFMN_ISSCAN_CGX_prediction TEXT,
    Carcinogenicity_model_IRFMN_Antares_assessment TEXT,
    Carcinogenicity_model_IRFMN_Antares_prediction TEXT,
    Carcinogenicity_oral_classification_model_IRFMN_assessment TEXT,
    Carcinogenicity_oral_classification_model_IRFMN_prediction TEXT,
    Carcinogenicity_oral_Slope_Factor_model_IRFMN_assessment TEXT,
    Carcinogenicity_oral_Slope_Factor_model_IRFMN_prediction TEXT,
    Carcinogenicity_inhalation_classification_model_IRFMN_assessment TEXT,
    Carcinogenicity_inhalation_classification_model_IRFMN_prediction TEXT,
    Carcinogenicity_inhalation_Slope_Factor_model_IRFMN_assessment TEXT,
    Carcinogenicity_inhalation_Slope_Factor_model_IRFMN_prediction TEXT,
    Carcinogenicity_in_male_rat_CORAL_assessment TEXT,
    Carcinogenicity_in_male_rat_CORAL_prediction TEXT,
    Carcinogenicity_in_female_Rat_CORAL_assessment TEXT,
    Carcinogenicity_in_female_Rat_CORAL_prediction TEXT,
    Acute_Toxicity_LD50_model_KNN_assessment TEXT,
    Acute_Toxicity_LD50_model_KNN_prediction TEXT,
    Hepatotoxicity_model_IRFMN_assessment TEXT,
    Hepatotoxicity_model_IRFMN_prediction TEXT,
    LogP_model_Meylan_Kowwin_assessment TEXT,
    LogP_model_Meylan_Kowwin_prediction TEXT,
    LogP_model_MLogP_assessment TEXT,
    LogP_model_MLogP_prediction TEXT,
    LogP_model_ALogP_assessment TEXT,
    LogP_model_ALogP_prediction TEXT
);
'''
cursor.execute(create_table_query)
connection.commit()

# Load CSV data into PostgreSQL
with open('data.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip header row
    for row in reader:
        cursor.execute(
            "INSERT INTO chemicals VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
            row
        )
connection.commit()

# Close the connection
cursor.close()
connection.close()
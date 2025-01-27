import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
import plotly.express as px

# Page configuration
st.set_page_config(page_title="Chemical Compound Explorer", layout="wide")
pd.set_option("styler.render.max_elements", 4842912)

# Data loading and caching
@st.cache_data
def load_data():
    df = pd.read_csv('planilha2.csv')
    return df

# Function to create molecule image
@st.cache_data
def create_molecule_image(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Draw.MolToImage(mol)
        return None
    except:
        return None

# Main application
def main():
    st.title('🧪 Chemical Compound Explorer (Leishmania)')
    
    # Load data
    try:
        df = load_data()
    except Exception as e:
        st.error(f"Error loading data: {str(e)}")
        return

    # Sidebar filters
    st.sidebar.header('📊 Filter Options')
    
    # Molecular weight filter
    mol_weight_range = st.sidebar.slider(
        'Molecular Weight',
        float(df['molecular_weight'].min()),
        float(df['molecular_weight'].max()),
        (float(df['molecular_weight'].min()), float(df['molecular_weight'].max()))
    )
    
    # LogP filter
    logp_range = st.sidebar.slider(
        'LogP',
        float(df['logp'].min()),
        float(df['logp'].max()),
        (float(df['logp'].min()), float(df['logp'].max()))
    )
    
    # Apply filters
    filtered_df = df[
        (df['molecular_weight'] >= mol_weight_range[0]) &
        (df['molecular_weight'] <= mol_weight_range[1]) &
        (df['logp'] >= logp_range[0]) &
        (df['logp'] <= logp_range[1])
    ]
    
    # Main content
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.write(f'📌 **Total Compounds:** {len(filtered_df)}')
        
        # Summary statistics
        st.subheader('📈 Summary Statistics')
        summary_stats = filtered_df[['molecular_weight', 'logp']].describe()
        st.dataframe(summary_stats.style.format("{:.2f}"))
    
    with col2:
        # Scatter plot
        st.subheader('🎯 Molecular Weight vs LogP')
        fig = px.scatter(
            filtered_df,
            x='molecular_weight',
            y='logp',
            hover_data=['compound_name'],
            title='Molecular Weight vs LogP Distribution'
        )
        st.plotly_chart(fig, use_container_width=True)
    
    # Data Analysis Section
    st.header('📊 Advanced Data Analysis')
    
    tab1, tab2, tab3 = st.tabs(["Distribution Plots", "Correlation Analysis", "Data Table"])
    
    with tab1:
        col3, col4 = st.columns(2)
        
        with col3:
            # Molecular Weight Distribution
            fig_mw = plt.figure(figsize=(10, 6))
            sns.histplot(data=filtered_df, x='molecular_weight', kde=True)
            plt.title('Molecular Weight Distribution')
            plt.xlabel('Molecular Weight')
            plt.ylabel('Count')
            st.pyplot(fig_mw)
        
        with col4:
            # LogP Distribution
            fig_logp = plt.figure(figsize=(10, 6))
            sns.histplot(data=filtered_df, x='logp', kde=True)
            plt.title('LogP Distribution')
            plt.xlabel('LogP')
            plt.ylabel('Count')
            st.pyplot(fig_logp)
    
    with tab2:
        # Correlation heatmap
        numeric_cols = filtered_df.select_dtypes(include=['float64', 'int64']).columns
        if len(numeric_cols) > 1:
            fig_corr = plt.figure(figsize=(10, 8))
            sns.heatmap(filtered_df[numeric_cols].corr(), annot=True, cmap='coolwarm', center=0)
            plt.title('Correlation Matrix')
            st.pyplot(fig_corr)
        else:
            st.info("Not enough numeric columns for correlation analysis")
    
    with tab3:
        # Interactive data table
        st.subheader('🔍 Compound Data')
        st.dataframe(
            filtered_df.style.highlight_max(axis=0, subset=['molecular_weight', 'logp']),
            height=400
        )
    
    # Chemical Structure Viewer (Optional)
    if st.checkbox('Show Chemical Structures'):
        st.subheader('⚗️ Chemical Structures')
        
        compound_names = filtered_df['compound_name'].unique()
        selected_compound = st.selectbox('Select a compound', compound_names)
        
        compound_row = filtered_df[filtered_df['compound_name'] == selected_compound].iloc[0]
        if 'smiles' in compound_row:
            mol_img = create_molecule_image(compound_row['smiles'])
            if mol_img:
                st.image(mol_img, caption=compound_row['compound_name'])
            else:
                st.warning(f"Could not generate structure for {compound_row['compound_name']}")

if __name__ == "__main__":
    main()
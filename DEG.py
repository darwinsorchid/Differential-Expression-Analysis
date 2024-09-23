# ---------------------------- Libraries -----------------------------------
import pandas as pd
import seaborn as sns
import numpy as np
import umap
import streamlit as st
from scipy.stats import f_oneway, dunnett
import statsmodels.api as sm
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.formula.api import ols
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA
from bokeh import plotting
from bokeh.models import HoverTool, ColumnDataSource, CategoricalColorMapper
from bokeh.palettes import cividis
import matplotlib.pyplot as plt
import sys

# Figure out how to make this work : sns.set_theme(style='white', context='notebook', rc={'figure.figsize':(14,10)})  

# ------------------------- Intro & File Upload -----------------------------
# Write header
st.markdown("# Differential Expression Analysis")

# Add info about project
st.markdown("An exploration of a differential expression analysis workflow.")

# Add file upload feature
file = st.file_uploader("Upload .tsv file")

# Check if file is uploaded
if file is None:
    st.markdown("Upload file to continue!")
    sys.exit(1)


# ---------------------- Load file and read data ----------------------------
genes_data = pd.read_csv(file, header =0, index_col=0, sep= '\t')

def create_boxplot(data):
    fig = plt.figure(figsize=(10,8))
    sns.boxplot(data = data, orient = 'h')
    plt.xticks(rotation=0, fontsize=8)
    st.pyplot(fig)

# Boxplot visualization before data quantile normalization 
st.markdown("## 1. Boxplot visualization before quantile normalization")

create_boxplot(genes_data)

# Add expander feature with info on figure
with st.expander("More info on figure"):
    st.write("This is a boxplot of..... [add info]")

# Display dataframe with streamlit
st.markdown("## Dataframe table")
st.dataframe(genes_data)


# ------------------------ Check For Missing Values ------------------------

st.markdown("## 2. Data cleaning: Check for missing values")

# Get pandas series of null values from original dataframe
missing_values = genes_data.isnull().sum()

# Convert series to pandas dataframe for better visualization
missing_values_df = pd.DataFrame({
    'Condition': missing_values.index,
    'Value': missing_values.values
})

# Display missing values dataframe with streamlit
st.table(missing_values_df) # table is static : no scrolling


# ------------------------- Check Data Distribution -------------------------

st.markdown("## 3. Data Distribution before quantile normalization")

def data_distribution_plots(data):
    # Create figure
    fig, axs = plt.subplots(nrows = 10, ncols = 7, figsize=(20, 15), sharey=True)

    # Flatten the axs array to make it easier to iterate
    axs = axs.flatten()

    # Add axes to figure with seaborn 
    for i, col in enumerate(data.columns):
        sns.kdeplot(data = data, x = col, ax = axs[i], color = "#009AEF")
        axs[i].set_title(col)

    # Remove empty subplots if any
    for j in range(i+1, len(axs)):
        fig.delaxes(axs[j])

    # Adjust the layout to prevent overlap
    plt.tight_layout()
    st.pyplot(fig)

# Display data distribution figure
data_distribution_plots(genes_data)


# ------------------------ Summarize numerical data --------------------------
st.markdown("## 4. Summarize Numerical Data")

# Create dataframe of min, max, mean, median, 1st quantile, 3rd quantile and standard deviation for every column/experimental condition
summary = genes_data.describe(include=[float, int]) # .describe() function excludes NA values by default
filtered_summary = summary.loc[['min','mean',  '25%', '50%', '75%', 'max', 'std']]

# Empty dictionary as placeholder
derived_mean = {}

# Calculate derived mean for every column of filtered_summary dataframe and save as dictionary
for col in filtered_summary.columns:
    d_mean = ((filtered_summary.loc['min',col] + filtered_summary.loc['max',col]) / 2).round(4)
    derived_mean[col] = d_mean

# Create pandas Series from derived_mean dictionary
new_row = pd.Series(derived_mean, name = 'derived_mean')
filtered_summary.loc['derived_mean'] = new_row

st.table(filtered_summary)


# ------------------------ Quantile Normalization -----------------------------
st.markdown("## 5. Data Normalization")

def quantile_normalization(df):
    ''' Function that goes through the necessary steps for data quantile normalization using pandas dataframe as input:
    1) Rank values of each column from smallest to largest
    2) Calculate average value of each rank --> MEAN
    3) Replace original data values with calculated average values
    ''' 
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()

    # Use Series, rank_mean, as a mapping for the ranks to get the normalized results:
    normalized_df = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()

    return normalized_df

norm_df = quantile_normalization(genes_data)
st.dataframe(norm_df) # might not need to display 

# Boxplot visualization after quantile normalization of data
st.markdown("## Boxplot after quantile normalization")
create_boxplot(norm_df)

# Data distribution after quantile normalization of data
st.markdown("## Data Distribution after quantile normalization")
data_distribution_plots(norm_df)


# ------------------------------- UMAP --------------------------------------
st.markdown("## 6. Dimension Reduction Algorithms")
st.markdown("### 6.1 Uniform Manifold Approximation and Projection (UMAP)")

# Keep only wild-type (WT) and transgenic (TG) conditions for simplicity
wt_tg_df = norm_df.iloc[ : , 1:23]

# Transpose dataframe for UMAP: genes as columns and samples as rows
wt_tg_df = wt_tg_df.transpose()
st.markdown("Transposed dataframe")
st.dataframe(wt_tg_df)

# Instantiate UMAP object
reducer = umap.UMAP(random_state = 42)

# UMAP algorithm learns a lower-dimensional representation of the input data
reducer.fit(wt_tg_df)

# Save embedding
embedding = reducer.transform(wt_tg_df)

# Verify that the result of calling transform is idenitical to accessing the embedding_ attribute
assert(np.all(embedding == reducer.embedding_))

# embedding.shape --> (22,2): Dataset of 22 samples and 2 feature columns [BY DEFAULT UMAP REDUCES DOWN TO 2D]

# Create list of names for samples' two categories: wild-type (wt) and transgenic(tg)
group = ['Wt'] * 10 + ['Tg'] * 12

# Convert categorical labels to numerical values (integers)
le = LabelEncoder()
group_encoded = le.fit_transform(group)

# ----------------------------- Interactive Plot with bokeh--------------------------------------

# Turn embedding to pandas DataFrame
#Each row in this DataFrame represents a single data point in the lower-dimensional space,
#with its corresponding coordinates in the 'x' and 'y' columns.
genes_df = pd.DataFrame(embedding, columns = ('x', 'y')) 

# Add column with sample names
genes_df['Samples'] = group
print(genes_df)

# Create a ColumnDataSource from genes dataframe
datasource_umap = ColumnDataSource(genes_df)

# Color mapping of sample names
color_mapping_umap = CategoricalColorMapper(factors = group, palette= cividis(3))

# Create interactive figure
p_umap = plotting.figure(title = "Interactive UMAP Visualization of WT and TG samples",
                    tools = ('pan, wheel_zoom, reset'),
                    width = 700, height = 600,
                    x_axis_label = 'UMAP1',
                    y_axis_label = 'UMAP2')

p_umap.title.text_font_size = "16px"

# Plot the data points as circles
p_umap.circle('x', 'y', size = 8, source = datasource_umap,
         color = dict(field='Samples', transform=color_mapping_umap),
         line_alpha=0.6, 
         fill_alpha=0.6)


# Add a HoverTool that shows the UMAP coordinates and sample name
p_umap.add_tools(HoverTool(tooltips="""
                      <div style="width: 100px;">
                        <div>
                            <div>
                                <span style="font-size: 12px; color: #233b42">Sample:</span>
                                <span style="font-size: 12px">@Samples</span>
                            </div>
                            <div>
                                <span style="font-size: 12px; color: #233b42">V1:</span>
                                <span style="font-size: 12px">@x</span>
                            </div>
                            <div>
                                <span style="font-size: 12px; color: #233b42">V2:</span>
                                <span style="font-size: 12px">@y</span>
                            </div>
                        </div>
                    """))

st.bokeh_chart(p_umap)


# ------------------------------- PCA --------------------------------------
st.markdown("### 6.2 Principal Component Analysis (PCA)")
pca_genes = PCA(n_components = 2)
principalComponents_genes = pca_genes.fit_transform(wt_tg_df)

#DataFrame that will have the principal component values of the samples
principal_gene_Df = pd.DataFrame(data = principalComponents_genes,
                                 columns = ('PC1', 'PC2'))
                                 
principal_gene_Df['Samples'] = group
st.dataframe(principal_gene_Df)

# Explained variance ratio (to see how much variance each component explains)
explained_variance = pca_genes.explained_variance_ratio_
st.markdown(f'Explained variability per principal component: {explained_variance}')

# Create a ColumnDataSource from genes dataframe
datasource_pca = ColumnDataSource(principal_gene_Df)

color_mapping_pca = CategoricalColorMapper(factors = group, palette= cividis(3))

# Create interactive figure
p_pca = plotting.figure(title = "Two first components of PCA",
                    tools = ('pan, wheel_zoom, reset'),
                    width = 700, height = 600,
                    x_axis_label = 'PC1',
                    y_axis_label = 'PC2')

p_pca.title.text_font_size = "16px"

# Plot the data points as circles
p_pca.circle('PC1', 'PC2', size = 8, source = datasource_pca,
         color = dict(field='Samples', transform=color_mapping_pca),
         line_alpha=0.6, 
         fill_alpha=0.6)


# Add a HoverTool that shows the UMAP coordinates and sample name
p_pca.add_tools(HoverTool(tooltips="""
                      <div style="width: 100px;">
                        <div>
                            <div>
                                <span style="font-size: 12px; color: #233b42">Sample:</span>
                                <span style="font-size: 12px">@Samples</span>
                            </div>
                            <div>
                                <span style="font-size: 12px; color: #233b42">PC1:</span>
                                <span style="font-size: 12px">@PC1</span>
                            </div>
                            <div>
                                <span style="font-size: 12px; color: #233b42">PC2:</span>
                                <span style="font-size: 12px">@PC2</span>
                            </div>
                        </div>
                    """))

st.bokeh_chart(p_pca)

# Add expander feature with info on figure
with st.expander("NOTE!"):
    st.write("PCA performed for the first two components explains 73% of variance,\nwhich is not an accurate representation.")

# --------------------------- STATISTICAL ANALYSIS -------------------------------
st.markdown("## 7. Statistical Analysis")
st.markdown("### 7.1 Group Treatments in Dataframe")

# Create matrix by excluding row names and column names
gene_matrix = np.matrix(norm_df)
st.dataframe(gene_matrix)

# Create groups
groups = (
    ["A_Wt"] * 10 +
    ["B_Tg"] * 13 +
    ["C_Proph_Ther_Rem"] * 3 +
    ["D_Ther_Rem"] * 10 +
    ["E_Ther_Hum"] * 10 +
    ["F_Ther_Enb"] * 10 +
    ["G_Ther_Cim"] * 10
)

# Convert the list into a pandas Series with categorical data type
groups_series = pd.Series(groups, dtype='category')

# Apply ANOVA on the first gene
# Create dataframe for gene1
gene1 = pd.DataFrame({"gene_expression": np.ravel(gene_matrix[1 , :]), "group":groups_series}) #np.ravel() to flatten multidimentional array

st.dataframe(gene1)

# Split the gene expression data by group
grouped_data = [gene1[gene1['group'] == g]['gene_expression'].values for g in gene1['group'].unique()]


# -------------------------------- ANOVA function on the first gene ------------------------------------
st.markdown("### 7.2 ANOVA")
f_stat, p_value = f_oneway(*grouped_data)

st.write("One-way ANOVA on gene1 results:")
st.write(f"F statistic: {f_stat} p-value: {p_value}")

model = ols('gene_expression ~ C(group)', data=gene1).fit()
aov_table = sm.stats.anova_lm(model, typ=2)
st.markdown("### Summary of ANOVA results")
st.dataframe(aov_table)

# Calculate mean expression value / group
group_mean_values = gene1.groupby('group')['gene_expression'].mean().reset_index()

# Rename columns
group_mean_values.columns = ['Group', 'Mean_Gene_Expression']

st.table(group_mean_values)

# -------------------------- Tukey's HSD post-hoc on the first (1st) gene ---------------------------
tukey = pairwise_tukeyhsd(endog=gene1['gene_expression'],
                          groups=gene1['group'],
                          alpha=0.05)
tukey

with st.expander("Interpretation"):
    st.write("Null hypothesis is not rejected:")
    st.write("No observed statistically significant difference between mean values across groups.")

# Convert tukey's summary to pandas dataframe to access specific metrics 
tukeys_df = pd.DataFrame(data = tukey.summary().data[1:], columns = tukey.summary().data[0])

# Summary for WT and TG conditions
st.write("Tukey's HSD summary for WT and TG conditions")
tukey_data = tukeys_df[(tukeys_df['group1'] == 'A_Wt')]
st.table(tukey_data.iloc[:, :4])


# ---------------------------- Perform Dunnett's post hoc test for gene1 ------------------------------
st.write("Dunnett's post hoc test for gene1")
samples = gene1.loc[gene1['group'] != 'A_Wt', 'gene_expression']
control = gene1.loc[gene1['group'] == 'A_Wt', 'gene_expression']

dunnett_test = dunnett(samples, control = control)

st.write(f"Dunnett's statistic: {dunnett_test.statistic}")
st.write(f"Dunnett's p-value: {dunnett_test.pvalue}")


# ------------------------------ Apply ANOVA on all genes ---------------------------------------------
st.write("One-way ANOVA on all genes of dataset:")

@st.cache_data
def process_gene(gene_data, groups_series):
    df = pd.DataFrame({"gene_expression": np.ravel(gene_data), "group": groups_series})
    model = ols('gene_expression ~ C(group)', data=df).fit()
    aov_table = sm.stats.anova_lm(model, typ=2)
    tukey = pairwise_tukeyhsd(endog=df['gene_expression'], groups=df['group'], alpha=0.05)
    tukeys_df = pd.DataFrame(data=tukey.summary().data[1:], columns=tukey.summary().data[0])
    tukey_data = tukeys_df[(tukeys_df['group1'] == 'A_Wt') & (tukeys_df['group2'] == 'B_Tg')]
    return aov_table, tukey_data

results_combined = pd.DataFrame()

#for i in range(len(gene_matrix[:, 1])):

for i, gene_data in enumerate(gene_matrix):
    aov_table, tukey_data = process_gene(gene_data, groups_series)

    # Create dataframe for each gene
    #df = pd.DataFrame({"gene_expression": np.ravel(gene_matrix[i , :]), "group":groups_series}) 

    # Apply ANOVA on each gene
    #model = ols('gene_expression ~ C(group)', data= df).fit()
    #aov_table = sm.stats.anova_lm(model, typ=2)

    # Apply Tukey's HSD post-hoc test on each gene
    #tukey = pairwise_tukeyhsd(endog=df['gene_expression'], groups=df['group'], alpha=0.05)

    # Turn tukey's summary to dataframe 
    #tukeys_df = pd.DataFrame(data = tukey.summary().data[1:], columns = tukey.summary().data[0])
    #tukey_data = tukeys_df[(tukeys_df['group1'] == 'A_Wt') & (tukeys_df['group1'] == 'B_Tg')]

    # Add gene index as a new column to both ANOVA and Tukey dataframes for tracking
    aov_table['gene'] = f'gene_{i + 1}'
    tukey_data['gene'] = f'gene_{i + 1}'

    # Concatenate ANOVA and Tukey's data for this gene
    combined_data = pd.concat([aov_table, tukey_data], axis=1, ignore_index=False)
    
    # Append the combined data to the main results dataframe
    results_combined = pd.concat([results_combined, combined_data], ignore_index=True)

print(results_combined)






    
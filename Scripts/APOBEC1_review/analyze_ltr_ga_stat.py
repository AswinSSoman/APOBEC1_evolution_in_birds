import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, pearsonr

# 1. Load data
df = pd.read_csv('all_ga_ltr_filtered.out', sep='\t')
df.columns = ['Group', 'Species', 'LTR_size', 'GA_edit_sites']

# 2. Calculate Editing Density (sites per Megabase)
df['Editing_Density'] = (df['GA_edit_sites'] / df['LTR_size']) * 1_000_000

# 3. Statistical Analysis
# Density vs Size (Spearman is best for non-linear relationships)
rho, p_spearman = spearmanr(df['LTR_size'], df['Editing_Density'])
# Raw Sites vs Size (for comparison)
r_pearson, p_pearson = pearsonr(df['LTR_size'], df['GA_edit_sites'])

# Save statistics to a file
with open('statistical_results.txt', 'w') as f:
    f.write("STATISTICAL ANALYSIS REPORT\n")
    f.write("===========================\n\n")
    f.write(f"Correlation between LTR Repertoire Size and Editing Density:\n")
    f.write(f"  - Spearman's Rho: {rho:.4f}\n")
    f.write(f"  - P-value: {p_spearman:.2e}\n\n")
    f.write(f"Correlation between LTR Repertoire Size and Raw Edit Sites:\n")
    f.write(f"  - Pearson's r: {r_pearson:.4f}\n")
    f.write(f"  - P-value: {p_pearson:.2e}\n\n")
    f.write("INTERPRETATION:\n")
    if p_spearman < 0.05 and rho > 0.5:
        f.write("Significant positive correlation found between size and density.\n")
        f.write("Conclusion: Editing is not a simple linear product of target number.\n")
        f.write("The probability of an average LTR sequence being edited increases\n")
        f.write("significantly as the total LTR content expands.\n")

# 4. Plotting
plt.figure(figsize=(7, 5))
sns.set_style("ticks")

# Distinct markers and color palette
markers = {"birds": "o", "mammals": "s", "invertebrate": "D", "non_mammalian_vertebrates": "^"}
palette = sns.color_palette("Set2")

plot = sns.scatterplot(
    data=df, 
    x='LTR_size', 
    y='Editing_Density', 
    hue='Group',
    style='Group',
    markers=markers,
    palette=palette,
    s=80,
    alpha=0.85,
    edgecolor='black',
    linewidth=0.6
)

# Formatting
plt.xscale('log')
plt.yscale('symlog', linthresh=1.0) 

# Adjusting limits to prevent clipping
plt.xlim(df['LTR_size'].min() * 0.5, df['LTR_size'].max() * 2)
plt.ylim(-0.5, df['Editing_Density'].max() * 2)

plt.xlabel('Total LTR Repertoire Size (bp)', fontsize=10)
plt.ylabel('GA Editing Density (sites/Mb)', fontsize=10)

# Compact legend inside the plot
plt.legend(title=None, loc='upper left', fontsize='8', frameon=True)

sns.despine()
plt.tight_layout()
plt.savefig('ltr_density_final.png', dpi=300, bbox_inches='tight')

print("Success: 'statistical_results.txt' and 'ltr_density_final.png' have been created.")

import pandas as pd


df = pd.read_csv('Orthogroup_Population_Stats.tsv', sep='\t')


populations = ['Cultivar', 'Xishuangbanna', 'Wild', 'Landrace']


summary = []

for pop in populations:
    mean_copy_col = f'{pop}_mean_copy'
    pav_col = f'{pop}_PAV'
    

    avg_copy = df[mean_copy_col].mean()
    

    avg_pav = df[pav_col].mean()
    

    dispensable_ratio = (df[pav_col] < 1).sum() / len(df)
    
    summary.append({
        'Population': pop,
        'Average_Copy_Number': avg_copy,
        'Average_PAV': avg_pav,
        'Dispensable_Gene_Ratio': dispensable_ratio
    })


summary_df = pd.DataFrame(summary)


print(summary_df)


summary_df.to_csv('Population_Diversity_Summary.tsv', sep='\t', index=False)

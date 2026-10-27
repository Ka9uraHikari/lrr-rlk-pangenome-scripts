import pandas as pd


region = pd.read_csv('/path/to/project/13_染色体定位/region.txt', 
                     sep='\t', header=None, names=['Sample', 'Population', 'Region'])


og = pd.read_csv('/path/to/project/09_OGG/RLKin39Genome/OrthoFinder/Results_Mar24/Orthogroups/Orthogroups.GeneCount.tsv',
                 sep='\t')


samples = region['Sample'].tolist()


og_samples = og[['Orthogroup'] + samples]


pop_map = dict(zip(region['Sample'], region['Population']))
pop_groups = region['Population'].unique()


results = []
for idx, row in og_samples.iterrows():
    og_id = row['Orthogroup']
    row_res = {'Orthogroup': og_id}
    for pop in pop_groups:
        pop_samples = [s for s in samples if pop_map[s]==pop]
        vals = row[pop_samples].values
        mean_copy = vals.mean()
        pav = (vals>0).sum() / len(vals)
        row_res[f'{pop}_mean_copy'] = mean_copy
        row_res[f'{pop}_PAV'] = pav
    results.append(row_res)

df_res = pd.DataFrame(results)
df_res.to_csv('Orthogroup_Population_Stats.tsv', sep='\t', index=False)
print("Done!")

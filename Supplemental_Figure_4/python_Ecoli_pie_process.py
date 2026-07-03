import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

path = 'goslim_counts_by_bracket_group.xlsx'
df = pd.read_excel(path)
df = df.rename(columns={c: c.strip() for c in df.columns})

col_goslim = next((c for c in df.columns if c.lower().replace('_','').replace('.','') in ['goslimcolumn','goslimcol','goslim']), None)
col_group  = next((c for c in df.columns if c.lower().replace('_','').replace('.','') in ['bracketgroup','group','bracket_group']), None)
col_term   = next((c for c in df.columns if c.lower().replace('_','').replace('.','') in ['goslimterm','term','goslim_term']), None)
col_count  = next((c for c in df.columns if c.lower().replace('_','').replace('.','') in ['count','counts','n']), None)

mask_p = df[col_goslim].astype(str).str.replace('_','.', regex=False).str.upper().str.contains('GOSLIM.P')
dfp = df[mask_p].copy()

agg_all = dfp.groupby([col_group, col_term], dropna=False)[col_count].sum().reset_index()
all_terms = sorted(agg_all[col_term].astype(str).unique())
cmap = plt.get_cmap('tab20', len(all_terms))
term_color = {term: cmap(i) for i, term in enumerate(all_terms)}

supergroups = {
    'Bulk 200 pg':                       ['200pg'],
    'SC intact stained':                 ['1x_intact_st'],
    'SC intact unstained':               ['1x_intact_us'],
    'SC spheroplast\nCephalexin stained': ['1x_protop_ceph_st'],
}

results = {}
for label, keywords in supergroups.items():
    mask = dfp[col_group].astype(str).str.lower().apply(lambda g: any(kw.lower() in g for kw in keywords))
    subset = dfp[mask]
    results[label] = subset.groupby(col_term, dropna=False)[col_count].sum().reset_index().sort_values(col_count, ascending=False)

fig, axes = plt.subplots(1, 4, figsize=(22, 7))

for idx, (label, sub) in enumerate(results.items()):
    ax = axes[idx]
    labels = sub[col_term].astype(str).tolist()
    sizes  = sub[col_count].astype(float).tolist()
    colors = [term_color[l] for l in labels]

    ax.pie(
        sizes, labels=None, colors=colors,
        autopct=lambda p: f"{p:.1f}%" if p >= 4 else '',
        startangle=90, counterclock=False,
        textprops={'fontsize': 12}
    )
    ax.axis('equal')
    ax.set_title(label, fontsize=13, fontweight='bold')

legend_handles = [mpatches.Patch(color=term_color[t], label=t) for t in all_terms]
fig.legend(
    handles=legend_handles, loc='lower center',
    ncol=4, fontsize=12, bbox_to_anchor=(0.5, -0.05),
    title='GO Biological Processes', title_fontsize=12, frameon=True
)

fig.suptitle('GO Biological Process', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig("process.png", bbox_inches='tight')

import pandas as pd
import numpy as np
import re

infile = 'final_Protein_Pivot_Report_with_GoSlim_unannotated_added.xlsx'
df = pd.read_excel(infile)

# Identify GoSlim columns
goslim_cols = [c for c in df.columns if isinstance(c, str) and 'goslim' in c.lower()]

# Identify bracket columns and group by base name after stripping leading [n]
bracket_cols = [c for c in df.columns if isinstance(c, str) and '[' in c]
base_map = {}
for c in bracket_cols:
    base = re.sub(r'^\s*\[\d+\]\s*', '', c).strip()
    base_map.setdefault(base, []).append(c)

# Non-empty / has a value: not NaN and not empty/whitespace

def has_value(s: pd.Series) -> pd.Series:
    notna = s.notna()
    out = pd.Series(False, index=s.index)
    if notna.any():
        ss = s[notna].astype(str)
        out.loc[notna] = ss.str.strip().ne('') & ss.str.lower().ne('nan')
    return out

results = []

for gcol in goslim_cols:
    goslim_has = has_value(df[gcol])
    goslim_val = df[gcol].where(goslim_has, np.nan).astype(object)
    goslim_val = goslim_val.apply(lambda x: x.strip() if isinstance(x, str) else x)

    for base, cols in sorted(base_map.items()):
        any_has = pd.Series(False, index=df.index)
        for c in cols:
            any_has = any_has | has_value(df[c])

        mask = goslim_has & any_has
        if int(mask.sum()) == 0:
            continue
        vc = goslim_val[mask].value_counts(dropna=True)
        for term, cnt in vc.items():
            results.append({
                'GoSlim_Column': gcol,
                'Bracket_Group': base,
                'GoSlim_Term': term,
                'Count': int(cnt)
            })

out_df = pd.DataFrame(results)

pivot_df = out_df.pivot_table(
    index=['GoSlim_Column', 'GoSlim_Term'],
    columns='Bracket_Group',
    values='Count',
    aggfunc='sum',
    fill_value=0
).reset_index()

outfile = 'goslim_counts_by_bracket_group.xlsx'
with pd.ExcelWriter(outfile, engine='openpyxl') as writer:
    out_df.sort_values(
        ['GoSlim_Column', 'Bracket_Group', 'Count'],
        ascending=[True, True, False]
    ).to_excel(writer, index=False, sheet_name='long')
    pivot_df.to_excel(writer, index=False, sheet_name='wide')

print('GoSlim columns:', goslim_cols)
print('Bracket groups:', len(base_map))
print('Rows in long output:', len(out_df))
print('Wrote', outfile)

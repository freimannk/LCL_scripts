import csv
import pandas as pd
import statsmodels.formula.api as smf

ge_norm_path = "GPL10558_SLE_healthy_gene_expression_wo_duplicate_probes_unique_subjects_INT_norm.tsv"
ge_df  = pd.read_csv(ge_norm_path, sep='\t')
ge_df.set_index("gene_id",inplace=True)
metadata = "USP18_lupus_study/GPL10558_SLE_healthy_patient_metadata.tsv"
metadata_df  = pd.read_csv(metadata, sep='\t')
metadata_df = metadata_df[["acc", "disease_state","age","gender","batch"]]
results_df = pd.DataFrame(columns=['gene', 'R-squared', 'Intercept', 'Coefficient SLE', 'P-value SLE','Coefficient gender', 'P-value gender','Coefficient age', 'P-value age','Coefficient batch', 'P-value batch'])
for gene_index, gene_row in ge_df.iterrows():
    merged_df = pd.merge(metadata_df, gene_row.rename('expression'), left_on="acc", right_index=True)
    formula = 'expression ~ disease_state + gender + age + batch'
    results = smf.ols(formula, data=merged_df).fit()
    results_df = results_df._append({
        'gene': gene_index,
        'R-squared': results.rsquared,
        'Intercept': results.params['Intercept'],
        'Coefficient SLE': results.params['disease_state[T.SLE]'],
        'P-value SLE': results.pvalues['disease_state[T.SLE]'],
        'Coefficient gender': results.params['gender[T.M]'],
        'P-value gender': results.pvalues['gender[T.M]'],
        'Coefficient age': results.params['age'],
        'P-value age': results.pvalues['age'],
        'Coefficient batch': results.params['batch'],
        'P-value batch': results.pvalues['batch']
    }, ignore_index=True)
results_df.to_csv(f"GPL10558_SLE_healthy_OLS_INT_normalized_w_batch.tsv", sep="\t", quoting=csv.QUOTE_NONE, index=False)




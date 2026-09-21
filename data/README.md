# Data

Save raw data files here.

## Versioned derived files

Most generated files live in `output/`, which is gitignored. Four are here
instead, and tracked:

| file | what it is |
|---|---|
| `GRM_Avena.rds`, `GRM_Pisum.rds` | genomic relationship matrices from `code/create_GRMs_T3.R` |
| `B4I_observations.rds`, `.csv.gz` | every observation from the selected trials, from `code/find_trials_with_B4I_accessions.R` |

They are derived, not raw, so strictly they belong in `output/`. They are here
because they are **inputs to everything downstream**, they take hours to
rebuild from T3, and versioning them means a fresh clone can run every model
and the whole simulation framework without a T3 login or a long download.

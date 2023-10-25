## Data
```
├─Data
│  │  facs_true_population.csv
│  │  mix_processed.csv
│  │  ref_13types.csv
│  │  tpm_mix_raw.csv
│  │
│  └─info
│          batch_info.csv
│          blood_biochemistry_values.csv
│          Mouse_stable2MouseMGI.csv
│          Sample_Summary.xlsx
```
-  ```facs_true_population.csv```: The ground truth cell type proportion matrix obtained by flow cytometry. (***P***)
-  ```mix_processed.csv```: Bulk gene expression for 11588 genes across 57 mouse liver injury samples. (***Y***)
-  ```ref_13types.csv```: Specific gene expression profiles for 13 cell types. (***X***)
-  ```tpm_mix_raw.csv```: TPM normalized gene expression processed at `/LiverDeconv/Sample_Codes/0_data_preprocessing.ipynb`.


- ```batch_info```: Batch information for sample acquisition.
- ```blood_biochemistry_values.csv```: ALT, AST and TBIL for each sample.
- ```Mouse_stable2MouseMGI.csv```: Transcripts ID <--> MGI symbols.
- ```Sample_Summary.xlsx```: Overview of the data acquired.

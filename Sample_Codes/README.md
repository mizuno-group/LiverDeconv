## Sample_Codes

-  ```0_data_preprocessing.ipynb```
   - How to down load raw file (fastq file).
   - How to obtain TPM normalized data.
   - Convert gene stable IDs to MGI symbols.
   - Performe various normlalization.
     - log2
     - cut off
     - batch normalization
     - quantile normalization
-  ```1_input_data.ipynb```
    - Mix data (Y)
    - LM13 (X)
    - True proportion (P)
    - Blood biochemistry values
-  ```2_simple_deconv_with_LM13.ipynb```
    - Detect differentially expressed genes (DEGs).
    - Conduct deconvolution based on Elastic Net.
-  ```3_reference_comb_optimization.ipynb```
    - Optimize the reference cell type combinations.
    - Evaluate the impact of the optimization.


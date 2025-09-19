# LiverDeconv
liver-specific deconvolution model for mouse  

<img src="https://github.com/mizuno-group/LiverDeconv/assets/92911852/ef517692-4c45-46f7-95dc-c2dd2c961928" width=800>

## Data
- RNA-seq data (fastq files) are available in GEO dataset. The accession number is [GSE237801](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237801).
- ```/Data``` directory
  -  **★```facs_true_population.csv```**: The ground truth cell type proportion matrix obtained by flow cytometry. (***P***)
  -  **★```mix_processed.csv```**: Bulk gene expression for 11588 genes across 57 mouse liver injury samples. (***Y***)
  -  ```ref_13types.csv```: Specific gene expression profiles for 13 cell types.
- ```SupplementaryFileS1.xlsx```: All cell name combinations including top 10 and bottom 10 used.
- **★ ```SupplementaryFileS2.xlsx```**: Bulk RNA-Seq-derived reference, a total of 503 DEGs across 13 cell types. (***X***)
  - The output of [Sample_Codes/2_simple_deconv_with_LM13.ipynb](https://github.com/mizuno-group/LiverDeconv/blob/main/Sample_Codes/2_simple_deconv_with_LM13.ipynb)
  - Used as the reference for benchmarking reference-based methods, including CIBERSORT
- ```SupplementaryFileS3.xlsx```: scRNA-Seq-derived reference a total of 184 DEGs across 10 cell types.
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
│ 
├─SupplementaryFileS1.xlsx
├─SupplementaryFileS2.xlsx
├─SupplementaryFileS3.xlsx
```

## Sample Code
- ```/Sample_Codes``` directory
  - ```0_data_preprocessing.ipynb```: guides the processing method of the acquired raw files and the storage location of the various data.
  - ```1_input_data.ipynb```: provides the shape of the data assumed for input.
  - ```2_simple_deconv_with_LM13.ipynb```: explains the procedure for performing a simple deconvolution method.
  - ```3_reference_comb_optimization.ipynb```: provides an example of reference optimization in the paper.

## Publication
- [peer-reviewed article](https://academic.oup.com/nargab/article/6/1/lqad111/7511730)
- [preprint](https://www.biorxiv.org/content/10.1101/2023.04.19.537436v3)  


## Citation
Please cite the following if you use any contents of this repository:  
  
```
@article{
    title={Investigation of the usefulness of liver-specific deconvolution method by establishing a liver benchmark dataset},
    author={Azuma I., Mizuno T., Morita K., Suzuki Y., & Kusuhara H.},
    journal={NAR Genomics and Bioinformatics},
    doi={10.1093/nargab/lqad111},
    URL={https://academic.oup.com/nargab/article/6/1/lqad111/7511730},
    year={2024}
}
```

## Authors
- [Iori Azuma](https://github.com/groovy-phazuma)  
    - main contributor  
- [Tadahaya Mizuno](https://github.com/tadahayamiz)  
    - correspondence  

## Contact
If you have any questions or comments, please feel free to create an issue on github here, or email us:  
- phazuma19980625[at]gmail.com  
- tadahaya[at]gmail.com  
    - lead contact  

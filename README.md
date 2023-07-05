# Grouping regulatory elements into enhancer, promoter, and ctcf sites
       
### Configuration:
- In ```config.yaml``` is indicated the *key*-*value* configuration with the respective documentation.

### Output folders:
- When running the pipeline, results will be automatically generated with all related outputs inside the folder named adequately inside the configuration file.

### Generting rule plots:
- ```snakemake --forceall --dag | dot -Tpdf > dagALL.pdf```
- ```snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf```

### Run snakemake:
- Run snakemake selecting number of cores (for parallelisation purpose) ```snakemake -s Snakefile --cores 1```



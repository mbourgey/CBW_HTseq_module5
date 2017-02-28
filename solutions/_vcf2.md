First, You can estimate the raw number of variant:

```
bcftools view  SVvariants/del.bcf| grep -v "^#" | wc -l

```

You should get 83. The merge will filter the low quality variants to retain only good candidate variants


Second you should notice that the there is no genotype in the resulting bcf.




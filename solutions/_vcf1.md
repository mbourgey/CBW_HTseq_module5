You can estimate the raw number of variant in each sample using the following command:

```
for i in SVvariants/*bcf ; do \ 
 echo $i ; \
 bcftools view $i | grep -v "^#" | wc -l  ;  \
done
```

You should get:

|NA12878|NA1289|NA12892|
|--|--|--|
|153|184|163|


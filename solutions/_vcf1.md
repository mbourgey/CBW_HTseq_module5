You can estimate the raw number of variant in each sample using the following command:

```bash
for bcf in SVvariants/*.bcf ; do
 basename $bcf .bcf
 bcftools view --no-header $i | grep -v "^#" | wc -l
done | paste - -
```

You should get:

|Sample|Count|
|--|--|
|NA12878|123|
|NA12891|155|
|NA12892|132|


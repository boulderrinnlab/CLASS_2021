#!/bin/bash
aws s3 sync /Shares/rinn_class/data/CLASS_2021/CLASS_2021/data/test_work/batch_1/results/bwa/mergedLibrary/bigwig  s3://bchm5631sp2021/bigwig/ --acl public-read
aws s3 sync /Shares/rinn_class/data/CLASS_2021/CLASS_2021/data/test_work/batch_2/results/bwa/mergedLibrary/bigwig  s3://bchm5631sp2021/bigwig/ --acl public-read
aws s3 sync /Shares/rinn_class/data/CLASS_2021/CLASS_2021/data/test_work/batch_3/results/bwa/mergedLibrary/bigwig  s3://bchm5631sp2021/bigwig/ --acl public-read
aws s3 sync /Shares/rinn_class/data/CLASS_2021/CLASS_2021/data/test_work/batch_4/results/bwa/mergedLibrary/bigwig  s3://bchm5631sp2021/bigwig/ --acl public-read
aws s3 sync /Shares/rinn_class/data/CLASS_2021/CLASS_2021/data/test_work/batch_5/results/bwa/mergedLibrary/bigwig  s3://bchm5631sp2021/bigwig/ --acl public-read
aws s3 sync /Shares/rinn_class/data/CLASS_2021/CLASS_2021/data/test_work/batch_6/results/bwa/mergedLibrary/bigwig  s3://bchm5631sp2021/bigwig/ --acl public-read

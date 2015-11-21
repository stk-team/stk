#!/bin/bash

nsamples=(16 64 128 256 512 1024 4096 10000)

for((ns=0;ns<8;ns++))
do
variance=$(./bin/variance-analysis-demo1-1 ./../results/ccvt-n${nsamples[ns]}-stratified/ccvt-sphere-stratified-p256-n${nsamples[ns]}-000001.txt ${nsamples[ns]} 60)
    echo ${nsamples[ns]} $variance >> variance-analysis-sphericalcap-theta60-ccvt.txt
done

#for((ns=0;ns<8;ns++))
#do
#    ./bin/variance-analysis-demo2-1 ./../results/ccvt-n${nsamples[ns]}-stratified/ccvt-sphere-stratified-p256-n${nsamples[ns]}-000001.txt ${nsamples[ns]} 4 0
#    echo ${nsamples[ns]} $variance >> variance-analysis-sphericalcap-ccvt.txt
#done

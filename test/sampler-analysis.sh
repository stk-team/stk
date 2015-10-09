#!/bin/sh

analyseSampler(){
	METHOD=$1
	N=$2
	M=$3
	FILENAME=./results/"$METHOD"T2-n$N-$M.mpts
	
	../bin/stk-sampler-"$METHOD"T2 -n $N -m $M -o "$FILENAME"
	../bin/stk-draw-pts -i "$FILENAME" -o ./results/"$METHOD"T2-pointset.png
	../bin/stk-fourier -i "$FILENAME" -O ./results/"$METHOD"T2-powerspec.png -R ./results/"$METHOD"T2-radialpowerspec.png -A ./results/"$METHOD"T2-anisotropy.png
}

analyseSampler whitenoise 4096 10
analyseSampler stratified 4096 10
analyseSampler fpo 4096 10
analyseSampler halton 4096 1
analyseSampler sobol 4096 1

#!/bin/sh -u

counter=1
while [ $counter -le $1 ]
do
	nohup Rscript -u single_gillespied.R "$counter" &
	counter=$(( counter+1 ))
done

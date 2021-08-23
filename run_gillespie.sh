#!/bin/sh

counter=1
While [ $counter -le $1]
do
	nohup Rscript random_walk.R "$counter" &
	((counter++))
done

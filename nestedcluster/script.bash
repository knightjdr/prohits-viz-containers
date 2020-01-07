#!/bin/bash

while getopts "m:p:" o; do
	case "${o}" in
		m)
      	matrix=${OPTARG}
      ;;
		p)
      	parameters=${OPTARG}
      ;;
    *)
      usage
    ;;
  esac
done

echo matrix:$matrix parameters:$parameters
echo running nestedcluster
/app/nestedcluster /files/$matrix /files/$parameters
echo generating images
Rscript /app/biclustering.R $matrix
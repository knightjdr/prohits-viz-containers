# prohits-viz-containers
 Containers for third-party tools used by ProHits-viz

## Nested cluster
> [source](https://sourceforge.net/projects/nestedcluster/)

### Build
```
docker build -t nestedcluster -f nestedcluster/Dockerfile .
```

### Usage
```
docker run -v $(pwd):/files/ nestedcluster -m matrix.txt -p parameters.txt
```

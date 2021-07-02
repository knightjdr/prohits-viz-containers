# Prohits-Viz containers
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

## RSVG
> [source](https://gitlab.gnome.org/GNOME/librsvg)

### Build
```
docker build -t rsvg -f rsvg/Dockerfile .
```

### Usage
```
docker run -v $(pwd):/files/ rsvg --format=png --output=./image.png --background-color=white --unlimited ./image.svg
```

## Utilities

### Python

#### Docker

1. Generate `requirements.txt` file
```
pipreqs utilities/python
```

2. Build
```
docker build -t pvutilitiespython -f utilities/python/Dockerfile .
```

#### Scripts

* saint summary statistics: `docker run -v $(pwd):/files/ pvutilitiespython /app/saint_stats.py -f 0.01 -s saint.txt`
* saint functional enrichment analysis: `docker run -v $(pwd):/files/ pvutilitiespython /app/saint_fea.py -f 0.01 -s saint.txt`

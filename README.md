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
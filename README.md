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

1. Generate `requirements.txt` file.
```
pipreqs utilities/python
```

2. Build
```
docker build -t pvutilitiespython -f utilities/python/Dockerfile .
```

#### Scripts

* crispr convert: `docker run -v $(pwd):/files/ --user $(id -u):$(id -g) pvutilitiespython /app/crispr_convert/main.py -f folder -t ranks`
* saint summary statistics: `docker run -v $(pwd):/files/ --user $(id -u):$(id -g) pvutilitiespython /app/saint_stats/main.py -f 0.01 -s saint.txt`
* saint functional enrichment analysis: `docker run -v $(pwd):/files/ --user $(id -u):$(id -g) pvutilitiespython /app/saint_fea/main.py -f 0.01 -s saint.txt`
* saint domain enrichment analysis: `docker run -v $(pwd):/files/ --user $(id -u):$(id -g) pvutilitiespython /app/saint_domain_enrich/main.py -b all -d domains.json -f 0.01 -g gene-db.json -i refseqp -s saint.txt`
* text biogrid network: `docker run -v $(pwd):/files/ --user $(id -u):$(id -g) pvutilitiespython /app/text_biogrid_network/main.py -k $access_key -f file.txt -g gene-db.json`
* text symbol fix: `docker run -v $(pwd):/files/ --user $(id -u):$(id -g) pvutilitiespython /app/text_symbol_fix/main.py -f file.txt -c "column1|column2"`

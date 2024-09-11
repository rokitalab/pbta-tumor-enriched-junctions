# neoepitope-identification

### To reproduce the code in this repository

1. Clone the repository:
```
git clone git@github.com:rokitalab/neoepitope-identification.git
```

2. Pull Docker container:
```
docker pull pgc-images.sbgenomics.com/d3b-bixu/splicing-neoepitopes:latest
```

3. Start the Docker container; from the `neoepitope-identification` folder, run:
```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/neoepitope-identification pgc-images.sbgenomics.com/d3b-bixu/splicing-neoepitopes:latest
```

4. Execute the shell within the docker image; from the `neoepitope-identification` folder, run: 
```
docker exec -ti <CONTAINER_NAME> bash
```
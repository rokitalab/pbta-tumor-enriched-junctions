# splicing-neoepitopes

### To reproduce the code in this repository

1. Clone the repository:
```
git clone git@github.com:d3b-center/splicing-neoepitopes.git
```

2. Pull Docker container:
```
docker pull pgc-images.sbgenomics.com/d3b-bixu/splicing-neoepitopes:latest
```

3. Start the Docker container; from the `splicing-neoepitopes` folder, run:
```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/splicing-neoepitopes pgc-images.sbgenomics.com/d3b-bixu/splicing-neoepitopes:latest
```

4. Execute the shell within the docker image; from the `splicing-neoepitopes` folder, run: 
```
docker exec -ti <CONTAINER_NAME> bash
```

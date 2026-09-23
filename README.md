# Identification of pediatric CNS tumor-enriched splicing events

### To reproduce the code in this repository

1. Clone the repository:
```
git clone git@github.com:rokitalab/tumor-enriched-splicing.git
```

2. Pull Docker container:
```
docker pull pgc-images.sbgenomics.com/rokita-lab/splicing-neoepitopes:latest
```

3. Start the Docker container

From the `tumor-enriched-splicing` folder, run:

```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/tumor-enriched-splicing pgc-images.sbgenomics.com/rokita-lab/splicing-neoepitopes:latest
```

Users can also run Rstudio in the project docker container from a web browser using the instructions below:

__Local Development in Rstudio__ (Max OS X and Linux users only)

```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/tumor-enriched-splicing pgc-images.sbgenomics.com/rokita-lab/splicing-neoepitopes:latest
```

Then, navigate to `localhost:8787` in your web browser. The username for login is `rstudio` and the password will be whatever password is set in the `docker run` command above (default: `pass`).

__Development using Amazon EC2, depending on your open ports__

```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 80:8787 -v $PWD:/home/rstudio/tumor-enriched-splicing pgc-images.sbgenomics.com/rokita-lab/splicing-neoepitopes:latest
```

Then, paste the instance IP address into your browser to start Rstudio. 

4. Execute the shell within the docker image; from the `tumor-enriched-splicing` folder, run: 
```
docker exec -ti <CONTAINER_NAME> bash
```

5. Run the `download-data.sh` shell script to obtain latest data files: 
```
bash download_data.sh
```

6. Navigate to an analysis module and run the shell script:
```
cd /home/rstudio/tumor-enriched-splicing/analyses/module_of_interest
```

### Below is the main directory structure listing the analyses and data files used in this repository

```
.
├── Dockerfile
├── LICENSE
├── README.md
├── analyses
│   ├── 00-pre-processing
│   ├── 01-tumor-specific-variants
│   ├── 02-translation
│   ├── find-gtex-aya
│   ├── process-ctrls-variants
│   ├── summarize-tumor-enriched-splicing
│   ├── tumor-specific-alt-splice-sites
│   └── tumor-specific-retained-introns
├── data
│   └── v8
├── doc
│   └── release-notes.md
├── download_data.sh
├── figures
│   └── theme_for_plots.R
├── scripts
│   └── download-original.sh
```


## Code Authors

Ryan Corbett ([@rjcorb](https://github.com/rjcorb)) and Ammar Naqvi ([@naqvia](https://github.com/naqvia))
**HapGUI**

HapGUI is a user-friendly GUI tool for performing haplotype, GWAS, and
haplo-pheno analysis designed to assist breeders, students, and
researchers. Follow the steps below to run HapGUI using Docker.

------------------------------------------------------------------------

## 📦 Prerequisites

Before running HapGUI, ensure that Docker is installed on your system.

1.  **Download Docker Desktop**  
    Visit the official Docker website to download Docker Desktop for
    your operating system:  
    👉 <https://www.docker.com/products/docker-desktop/>

2.  **Install Docker Desktop**  
    Follow the installation instructions provided in the video tutorial.

------------------------------------------------------------------------

## 🚀 Run HapGUI via Docker

Once Docker is installed, you can launch HapGUI by running the following
commands **once** in your terminal:

docker pull irrisah2012/hapgui:v1 
docker run -d -p 3232:3232 irrisah2012/hapgui:v1

## 🌐 Access the HapGUI Interface

Once the container is running:

1.  Open a web browser.
2.  **Paste the following link into the address bar**:  
    👉 <http://localhost:3232/>
3.  Press **Enter**.
4.  The HapGUI interface will open in a few seconds and is ready to use.

------------------------------------------------------------------------

## 📂 User Manual
 A user manual is available in the `tutorial/` (https://github.com/IRRI-South-Asia-Hub/HapGUI/tree/final/tutorial) folder of this repository for input file preparation. For running HapGUI please follow the Video Tutorial provided below.

------------------------------------------------------------------------

## 📺 Video Tutorial

A complete **video manual** is available to help you:

- Install Docker
- Run HapGUI in Docker
- Use the HapGUI and access results for respective analysis ( Phenotype analysis, GWAS, Et-WAS, Haplopheno analysis)

Please find the video manual (link) of this repository for easy offline reference.

------------------------------------------------------------------------
## 🛠️ Support

If you encounter any issues, please:

- Open a GitHub issue
  [here](https://github.com/IRRI-South-Asia-Hub/HapGUI/issues)



%%%%% README %%%%%%%%%%

-First, install the contents of bncmrk-dblts github to your folder of choice. Preferably, not your home folder, as the pipeline consumes lots of memory, but to rather /projects/ucar-lab/XX/bncmrk-dblts   (XX: whatever your username)

-Second, create a folder named Data. Within that folder, create folder per project name. e.g, CZI.PBMC, PBMC.8.HTO, Four.Cell.12.HTO

-Third, download HTO and RNA/umis data sets to their respective folders

-Fourth, run the pipeline within "organize_data_script.R" in bncmrk-dblts/R folder. You may need to modify the bncmrk-dblts homefolder location within script with respect to your configuration. Also, the new HTO data you upload may need to follow the logic within the script.

- "organize_data_script.R": will clean empty cells/barcodes not containing RNAs. Will find common genes among the given datasets, toss out rest of the genes. In our case you end up with ~18k genes. Sort the RNA matrix rows with respect to genes (each cell/barcode is a column at this point). Script saves these genes for later use, e.g as neural network features. It saves each rna/umi matrix in the following naming convention: "project_name" + "RNA.raw.matrix.sparse.Rds". e.g, Four.Cell.12.HTO.RNA.raw.matrix.sparse.Rds



-"organize_data_script.R": It will take HTOs, and by using either Generalized Mixture Models (GMM), or HTODemux, or MULTIseqDemux method to figure out the composition (annotations), and rather simply classifications. Save them to bncmrk-dblts/GroundTruths/"project_name" folder with naming conventions: "project_name" + "HTO.sample.annotations" + "classification method" + ".csv", and,  "project_name" + "HTO.sample.classifications" + "classification method" + ".csv". e.g, save "CZI.PBMC.HTO.sample.annotations.GMM.csv" files to bncmrk-dblts/GroundTruths/CZI.PBMC folder

-Fifth, modify "counts2software_script.R" by adding the projects names. e.g, if you have data for CZI.PBMC project, make sure you have this line in your script: 
"counts2software_pipe(pr_name="CZI.PBMC",foldersPresent=FALSE,numSamples = 10)". If you have more than one samples, such as in our CZI.PBMC data, make sure you add it. Because neither Seurat, nor other softwares running behind the scenes can't process large data (matrices greater than 40k column/row size). e.g, Even a simple operation such as matrix normalization crashes them.

Running this script will create DF file for DoubletFinder, SCR for Scrublet, DD for DoubletDecon. It creates input and output folders in each of those files named after softwares. Within each of these input and output files, it will create folders named after the project. For example, it will create the folder  "bncmrk-dblts/DD/input/CZI.PBMC". And, the script will create files files necessary to run the pipelines of those software in these "bncmrk-dblts"+"software name" + "input" + "project name" folders.  If you want to run them via batch jobs on cluster, you can find necessary PBS scripts in bncmrk-dblts/PBS folders

-Sixth, run preprocessing pipes through "project name_" + "preprocess_script.R" scripts. e.g, After some modification (add each of the project names to pr_names list), for the project CZI.PBMC, running "doubdec_preprocess_script.R" takes the input files from bncmrk-dblts/DD/input/CZI.PBMC, preprocess them, and generate the actual input files necessary to run them, save them in their respective input folders again. The actual workhorse are the function files named "software name_"+"preprocess_pipe.R". Script is running on them per project basis. If you want to run them via batch jobs on cluster, you can find necessary PBS scripts in bncmrk-dblts/PBS folders

-Seventh, run out pipes through "project name_" + "out_script.R" scripts. e.g, After some modification (add each of the project names to pr_names list), for the project CZI.PBMC, running "doubdec_out_script.R" takes the input files from bncmrk-dblts/DD/input/CZI.PBMC, outputs them, and generate output files, save them in their respective output folders. The actual workhorse are the function files named "software name_"+"out_pipe.R". Script is running on them per project basis.  If you want to run them via batch jobs on cluster, you can find necessary PBS scripts in bncmrk-dblts/PBS folders

The outpipe for scrublet is in bncmrk-dblts/Python. You can find the PBS script for it as well in PBS folder.

-There are various Jupyter Notebooks in bncmrk-dblts/R/Notebooks and bncmrk-dblts/Python/Notebooks folders. You can use them for analyzing ground truths and such. They are not that organized, mostly used for unit testing and experimentation.

-If you want to use Jupyter Notebooks with Python and R support, you need to create your own Anaconda container. Within your Anaconda container, you need to install needed R and Python dependencies by using conda and pipe. Obtain an interactive session by bashing/or typing the commands in interactive.sh. Then, type or bash the commands in jupyter.sh. Copy-paste the link given on terminal to your browser. You are using Jupyter notebooks now.


# CUT&RUN Snakemake Workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)

## Authors

* Robin Meyers (@robinmeyers)

## Usage

### Simple

#### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/robinmeyers/cutnrun-snakemake/releases).
If you intend to modify and further extend this workflow or want to work under version control, fork this repository as outlined in [Advanced](#advanced).

Clone this repositiory into a directory

```
$ git clone git@github.com:robinmeyers/cutnrun-snakemake.git my-project-name
$ cd my-project-name
```

Create and activate the conda environment

```
$ conda env create -f envs/conda.yaml -n cutnrun-snakemake
$ conda activate cutnrun-snakemake
```

Run the snakemake on test data

```
$ snakemake --directory .test
```

Examine the outputs of the workflow in the directory ```.test/outs/```


#### Step 2: Configure workflow

Configure the workflow by editing the files `config.yaml` and creating a samplesheet.

The samplesheet has the following format (other columns than these are ignored):
```txt
fastq,sample,control,condition
H3K27ac_D0_209_chr21,H3K27ac_D0_209,IgG_D0_209,H3K27ac_D0
H3K27ac_D2_209_chr21,H3K27ac_D2_209,IgG_D2_209,H3K27ac_D0
H3K4me3_D0_209_chr21,H3K4me3_D0_209,IgG_D0_209,H3K4me3_D0
H3K4me3_D2_209_chr21,H3K4me3_D2_209,IgG_D2_209,H3K4me3_D0
IgG_D0_209_chr21,IgG_D0_209,,IgG_D0
IgG_D2_209_chr21,IgG_D2_209,,IgG_D2
```
Control samples are included in the samplesheet with the "control" column blank.

Fastq files are identified by matching the prefix in the fastq column and the regex pattern: "(_S[0-9]+)?(_L[0-9]+)?_R1(_001)?.fastq.gz".

That is, they must contain R1/R2 for read1/read2 and end in the suffix .fastq.gz.


#### Step 3: Execute workflow

Ensure the conda environment is active:

```
$ conda activate cutnrun-snakemake
```

Test your configuration by performing a dry-run:

```
$ snakemake -np
```

Execute the workflow locally using `$N` cores:

```
$ snakemake -p --cores $N
```

The workflow can be executed on a cluster using snakemake cluster configuration. Install a [profile](https://github.com/Snakemake-Profiles) for your cluster's job submission system. Edit the defaults in the file `cluster.json` and run the workflow. For example:


```
$ snakemake -p --jobs $N --profile slurm --cluster-config cluster.json
```



### Updating the workflow

If you installed the workflow by cloning the github repo, you can pull latest updates to workflow with 

```$ git pull --rebase```

This will require you to first commit any changes you made to your configuration file before pulling new updates.

Then simply rerun the `snakemake` command.

### Advanced

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for the concrete project/run on your machine.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).
4. Modify the config, and any necessary sheets (and probably the workflow) as needed.
5. Commit any changes and push the project-branch to your fork on github.
6. Run the analysis.
7. Optional: Merge back any valuable and generalizable changes to the [upstream repo](https://github.com/robinmeyers/cutnrun-snakemake) via a [**pull request**](https://help.github.com/en/articles/creating-a-pull-request). This would be **greatly appreciated**.
8. Optional: Push results (plots/tables) to the remote branch on your fork.
9. Optional: Create a self-contained workflow archive for publication along with the paper (snakemake --archive).
10. Optional: Delete the local clone/workdir to free space.


## Testing

Tests cases are in the subfolder `.test`. They are automtically executed via continuous integration with Travis CI.

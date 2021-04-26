Here we are going to import the RNA-seq information from encode.
================================================================

Goal: to get a rational samplesheet for use in downstream analysis (ie., differential expression) and to rename the files in the format needed by nf-core/rnaseq. You need the read numbers to be a part of the filename to run with nextflow.
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

[Encode Query](%22https://www.encodeproject.org/report/?type=Experiment&status=released&assay_slims=Transcription&assay_slims=Transcription&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=HepG2&biosample_ontology.classification=cell+line&assay_title=total+RNA-seq&files.read_length=50&limit=all&advancedQuery=date_released:%5B2009-01-01%20TO%202021-12-31%5D%22)

``` bash
# Retrieve experiment info table
#wget -O samples.txt "https://www.encodeproject.org/report.tsv?type=Experiment&status=released&assay_slims=Transcription&assay_slims=Transcription&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=HepG2&biosample_ontology.classification=cell+line&assay_title=total+RNA-seq&files.read_length=50&limit=all&advancedQuery=date_released:[2009-01-01%20TO%202021-12-31]"

# Retrieve fastq file urls
#wget -O files.txt "https://www.encodeproject.org/batch_download/?type=Experiment&status=released&assay_slims=Transcription&assay_slims=Transcription&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=HepG2&biosample_ontology.classification=cell+line&assay_title=total+RNA-seq&files.read_length=50&limit=all&advancedQuery=date_released:[2009-01-01%20TO%202021-12-31]"


# You can either use the terminal and mkdir or create the directory in R.
#if(!dir.exists("/scratch/Shares/rinnclass/tardigrades/CLASS_2021/analysis/05_RNA-seq_expression/rnaseq/fastq")) dir.create("/scratch/Shares/rinnclass/tardigrades/CLASS_2021/analysis/05_RNA-seq_expression/rnaseq/fastq")
# Download the fastq files -- ~50 GB
#xargs -L 1 curl -O -J -L < files.txt
```

``` r
# We'll also rename this Accession column to clarify between experiment_accession and file_accession.
samples <- read.table("/scratch/Shares/rinnclass/tardigrades/CLASS_2021/analysis/05_RNA-seq_expression/samples.txt",
                      sep = "\t", skip = 1, header = T) %>%
  dplyr::rename(experiment_accession = Accession) 
```

We will need to clean this up as well as append new data to this (such as the md5 checksum for each file). You'll notice that there is no md5 column in this table, so that's why we need to supplement it with other data from ENCODE.
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

For that data we will use ENCODE's API to retrieve additional file information from their server.
-------------------------------------------------------------------------------------------------

``` r
library(httr)

# The base url we'll use for querying ENCODE is as follows:
# https://www.encodeproject.org/report.tsv?

# Let's look at an example request for this experiment accession: ENCSR541TIG
request_url <- "https://www.encodeproject.org/report.tsv?type=File&status=released&file_format=fastq&dataset=%2Fexperiments%2FENCSR541TIG%2F&field=accession&field=read_count&field=md5sum&field=controlled_by&field=paired_end&field=paired_with&field=replicate&field=target"

# Retrieving data for R session.
(resp <- GET(request_url))
```

    ## Response [https://www.encodeproject.org/report.tsv?type=File&status=released&file_format=fastq&dataset=%2Fexperiments%2FENCSR541TIG%2F&field=accession&field=read_count&field=md5sum&field=controlled_by&field=paired_end&field=paired_with&field=replicate&field=target]
    ##   Date: 2021-04-26 20:13
    ##   Status: 200
    ##   Content-Type: text/tsv; charset=UTF-8
    ##   Size: 900 B
    ## 2021-04-26 20:13:11.403666   https://www.encodeproject.org/report/?type=File&status=released&file_format=fastq&dataset=%2Fexperiments...
    ## Accession    Read count  MD5sum  Controlled by   Paired end identifier   Paired with Replicate   Target
    ## ENCFF969RUJ  20236753    21f1d2d06b53bd6ee456c67e629cc734        1   /files/ENCFF448RSL/ /replicates/00781b21-22c3-4be2-ac4f-54422fd5da0e/   
    ## ENCFF448RSL  20236753    9541757b0d6b243e2d4cccbd0cc463e6        2   /files/ENCFF969RUJ/ /replicates/00781b21-22c3-4be2-ac4f-54422fd5da0e/   
    ## ENCFF714MGS  24235894    9ec27b60e022922c0bf34b6b14b2bf0f        1   /files/ENCFF675FVZ/ /replicates/6ea50902-8c42-4940-aed5-e529a8dcfbd4/   
    ## ENCFF675FVZ  24235894    152e103a1911c625d206e11d88179884        2   /files/ENCFF714MGS/ /replicates/6ea50902-8c42-4940-aed5-e529a8dcfbd4/   

``` r
# This object contains the data we want, but also contains header information about the request, extract just the data (body) portion
body <- read_tsv(content(resp, "text"), skip = 1)
```

Generate a request url for our requirements from ENCODE.
========================================================

``` r
## This will generate a request URL in the format that ENCODE requires to retrieve each of the columns listed in the field default parameter (accession, read_count, md5sum, etc.)
contstruct_query <- function(experiment_accession,
                             base_url = "https://www.encodeproject.org/report.tsv?",
                             file_format = "fastq",
                             type = "File",
                             status = "released",
                             fields = c("accession", "read_count", "md5sum",
                                        "controlled_by", "paired_end",
                                        "paired_with", "replicate", "target")) {
  query <- paste(list(paste0("type=", type),
                      paste0("status=", status),
                      paste0("file_format=", file_format),
                      paste0("dataset=%2Fexperiments%2F", experiment_accession, "%2F"),
                      map_chr(fields, ~paste0("field=", .))) %>%
                   flatten(),
                 collapse = "&")
  url <- paste0(base_url, query)
  return(url)
}

# This function actually makes the request and returns the data only (without the response headers) in a data.frame format.
encode_file_info <- function(experiment_accession,
                             base_url = "https://www.encodeproject.org/report.tsv?",
                             file_format = "fastq",
                             type = "File",
                             status = "released",
                             fields = c("accession", "read_count", "md5sum",
                                        "controlled_by", "paired_end",
                                        "paired_with", "replicate", "target")) {
  path <- "report.tsv?"
  base_url <- modify_url("https://www.encodeproject.org/", path = path)
  url <- contstruct_query(experiment_accession,
                          base_url = base_url,
                          file_format,
                          type,
                          status,
                          fields)
  resp <- GET(url)
  if (http_error(resp)) {
    error_message <- content(resp, type = "text/html", encoding = "UTF-8") %>%
      xml_find_all("//p") %>%
      xml_text() %>%
      first()
    stop(
      sprintf(
        "ENCODE API request failed [%s]\n%s",
        status_code(resp),
        error_message
      ),
      call. = FALSE
    )
  }
  
  if (http_type(resp) != "text/tsv") {
    stop("API did not return text/tsv", call. = FALSE)
  }
  body <- read_tsv(content(resp, "text"), skip = 1) %>%
    clean_names()
  return(body)
}

dat <- encode_file_info("ENCSR541TIG")
```

Map the experiments with the ENCODE files.
==========================================

We can use map to call this encode\_file\_info function for each experiment accession and return each data.frame into a new column in our samples data.frame.
-------------------------------------------------------------------------------------------------------------------------------------------------------------

``` r
samples <- samples %>%
  mutate(file_info = map(experiment_accession, ~ encode_file_info(.x)))

# Note that in this new column we have a bunch of data.frames -- nested in our samples data.frame.

# We also need to tell it which column we want to unnest.
samples <- samples %>%
  unnest(cols = file_info)
```

Cleaning the data frame and organizing by accesion and rep number.
==================================================================

``` r
# Let's number our replicates and create a new column called sample id where we join the experiment accesion and the rep number
samples <- samples %>%
  group_by(experiment_accession) %>%
  mutate(rep_number = as.numeric(factor(replicate))) %>%
  # ?unite -- will combine two data columns into one. Be careful, by default it will remove the columns you're uniting.
  unite(sample_id, experiment_accession, rep_number, sep = "_rep")

# We're just selecting a subset of columns and 
samples <- samples %>%
   dplyr::select(sample_id, accession, Assay.title,
                Biosample.summary, md5sum, paired_end_identifier) %>%
  # Note that we're not removing the original columns
  unite(fastq_file, sample_id, paired_end_identifier, sep = "_read", remove = F)
```

Modifying the file names to match the nf-core RNA-seq pipeline.
===============================================================

``` r
samples <- samples %>%
  mutate(fq_extension = ".fastq.gz") %>%
  unite(fastq_file, fastq_file, fq_extension, sep = "", remove = F) %>%
  # This original file column will be used along with the new name column to rename the fastq files.
  unite(original_file, accession, fq_extension, sep = "")
```

Creating the bash commands with R.
==================================

``` r
# Rename the fastq files so that they contain the sample ID.
rename_script <- samples %>%
  dplyr::select(fastq_file, original_file) %>%
  mutate(command = "mv") %>%
  unite(command, command, original_file, fastq_file, sep = " ")
# The result of this is that each row is a bash command.

# We can write this out as a bash script with ?write_lines 
# We include a shebang header line so that the script is interpreted by bash.

write_lines(c("#!/bin/bash", rename_script$command), "/scratch/Shares/rinnclass/tardigrades/CLASS_2021/analysis/05_RNA-seq_expression/rnaseq/fastq/rename.sh")
# Here we use an R command to call bash and cd into that directory, then make the file executable, and then run it. We're not going to run it now because we have already run it before.
#system("cd /scratch/Shares/rinnclass/tardigrades/CLASS_2021/analysis/05_RNA-seq_expression/rnaseq/fastq; chmod u+x rename.sh; ./rename.sh")
```

Create a text file to run the md5sum check.
===========================================

``` r
# Let's create an md5.txt to run the checksums
md5 <- samples %>% 
  dplyr::select(md5sum, fastq_file) %>%
  unite(line, md5sum, fastq_file, sep = "  ")
write_lines(md5$line, "/scratch/Shares/rinnclass/tardigrades/CLASS_2021/analysis/05_RNA-seq_expression/rnaseq/fastq/md5.txt")
# We won't run this line now either.
#system("cd /scratch/Shares/rinnclass/tardigrades/CLASS_2021/analysis/05_RNA-seq_expression/rnaseq/fastq; md5sum -c md5.txt")
```

Format sample sheet that we will use downstream for further analysis of the read counts in R.
=============================================================================================

Add the header for the experiment sample.
-----------------------------------------

``` r
# Let's create the sample sheet that we will use later
# to do the RNA-seq analysis in R.
samplesheet <- samples %>%
  dplyr::rename(fastq = fastq_file,
                seq_type = Assay.title,
                sample_name = Biosample.summary) %>%
  # The minus sign will remove this column -- which we no longer need.
  dplyr::select(-original_file) 


# Now that we have it cleaned up, let's create one line for each replicate where the fastq read 1 and read 2 are in the same row.
# For this we will use the pivot wider function
# We need to tell the pivot_wider function which unique column combinations will specify each new row. 
suppressWarnings(samplesheet <- samplesheet %>%
  pivot_wider(id_cols = c("sample_id", "seq_type", "sample_name"),
              names_from = paired_end_identifier,
              values_from = c("fastq", "md5sum")))



# Harmonize the sample names with the design file.
samplesheet <- samplesheet %>%
  mutate(condition = gsub(" ", "_", sample_name) %>% tolower()) %>%
  separate(sample_id, into = c("experiment_accession", "replicate"), 
           remove = FALSE, sep = "_") %>%
  mutate(replicate = gsub("rep", "R", replicate)) %>%
  unite(sample_name, condition, replicate, sep = "_", remove = FALSE)

samplesheet$condition[samplesheet$condition == "hepg2"] <- "hepg2_total"
samplesheet <- samplesheet %>%
  mutate(cell_type = "hepg2",
         condition = gsub("hepg2_", "", condition)) %>%
  dplyr::select(sample_id, sample_name, replicate, condition,
                cell_type, seq_type, fastq_1, fastq_2, md5sum_1,
                md5sum_2)

write_csv(samplesheet, "/scratch/Shares/rinnclass/tardigrades/CLASS_2021/analysis/05_RNA-seq_expression/rnaseq/samplesheet.csv")
```

Create the design csv for the nf-core RNA-seq pipeline.
=======================================================

``` r
# This will have the following columns:
# group,replicate,fastq_1,fastq_2,strandedness
# The group will specify the experimental condition 
# We could use the experiment accesion or we could use the 
# sample_name.
# Let's create a group column using mutate -- and we'll clean up the names
# There's spaces in the names currently and we need to get rid of those.
design <- samplesheet %>%
  # We also need to add the proper file path
  mutate(group = gsub(" ", "_", sample_name) %>% tolower(),
         strandedness = "reverse",
         fastq_1 = paste0("fastq/", fastq_1),
         fastq_2 = paste0("fastq/", fastq_2)) %>%
  # We have the replicate number already, but it's just a part of the sample id
  # We can use the separate function (opposite of unite) to retrieve it.
  separate(sample_id, into = c("experiment_accession", "replicate"), sep = "_", remove = FALSE) %>%
  mutate(replicate = gsub("rep", "", replicate)) %>%
  # Now we gather just the columns we need
  dplyr::select(group, replicate, fastq_1, fastq_2, strandedness)
# Now we can write this out for input into nextflow
write_csv(design, "/scratch/Shares/rinnclass/tardigrades/CLASS_2021/analysis/05_RNA-seq_expression/rnaseq/design.csv")
```

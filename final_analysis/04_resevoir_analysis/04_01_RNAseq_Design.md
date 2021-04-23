Goal: To make the RNAseq design file required to run RNAseq to be analyzed in 04\_resevoir\_analysis.Rmd

\[Encode Link\] ("<https://www.encodeproject.org/report/?type=Experiment&status=released&assay_slims=Transcription&assay_slims=Transcription&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=HepG2&biosample_ontology.classification=cell+line&assay_title=total+RNA-seq&files.read_length=50&limit=all&advancedQuery=date_released>:\[2009-01-01%20TO%202021-12-31\]")

Retrieve RNAseq count data from encode:

``` bash
# You'll want to change directory to the rnaseq directory for this following bit.
# Retrieve experiment info table
# cd to rnaseq directory (/scratch/Shares/rinnclass/<your_dir>/rnaseq)
wget -O samples.txt "https://www.encodeproject.org/report.tsv?type=Experiment&status=released&assay_slims=Transcription&assay_slims=Transcription&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=HepG2&biosample_ontology.classification=cell+line&assay_title=total+RNA-seq&files.read_length=50&limit=all&advancedQuery=date_released:[2009-01-01%20TO%202021-12-31]"

# Retrieve fastq file urls
wget -O files.txt "https://www.encodeproject.org/batch_download/?type=Experiment&status=released&assay_slims=Transcription&assay_slims=Transcription&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=HepG2&biosample_ontology.classification=cell+line&assay_title=total+RNA-seq&files.read_length=50&limit=all&advancedQuery=date_released:[2009-01-01%20TO%202021-12-31]"

# Now cd to the rnaseq/fastq directory 
# If the fastq directory doesn't exist, let's create it.
# You can either use the terminal and mkdir or create the directory in R.
# if(!dir.exists("rnaseq/fastq")) dir.create("rnaseq/fastq")
# Download the fastq files -- ~50 GB
# xargs -L 1 curl -O -J -L < files.txt
```

    ## --2021-04-23 14:14:02--  https://www.encodeproject.org/report.tsv?type=Experiment&status=released&assay_slims=Transcription&assay_slims=Transcription&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=HepG2&biosample_ontology.classification=cell+line&assay_title=total+RNA-seq&files.read_length=50&limit=all&advancedQuery=date_released:[2009-01-01%20TO%202021-12-31]
    ## Resolving www.encodeproject.org (www.encodeproject.org)... 34.211.244.144
    ## Connecting to www.encodeproject.org (www.encodeproject.org)|34.211.244.144|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: unspecified [text/tsv]
    ## Saving to: ‘samples.txt’
    ## 
    ##      0K ......                                                 5.77M=0.001s
    ## 
    ## 2021-04-23 14:14:02 (5.77 MB/s) - ‘samples.txt’ saved [6478]
    ## 
    ## --2021-04-23 14:14:02--  https://www.encodeproject.org/batch_download/?type=Experiment&status=released&assay_slims=Transcription&assay_slims=Transcription&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=HepG2&biosample_ontology.classification=cell+line&assay_title=total+RNA-seq&files.read_length=50&limit=all&advancedQuery=date_released:[2009-01-01%20TO%202021-12-31]
    ## Resolving www.encodeproject.org (www.encodeproject.org)... 34.211.244.144
    ## Connecting to www.encodeproject.org (www.encodeproject.org)|34.211.244.144|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: unspecified [text/plain]
    ## Saving to: ‘files.txt’
    ## 
    ##      0K .                                                      3.17M=0.001s
    ## 
    ## 2021-04-23 14:14:02 (3.17 MB/s) - ‘files.txt’ saved [1996]

So now we have the samples (in this case, cellular fractionations from HepG2 cell line), and we want to read this into R to start making the design file.

``` r
# We'll also rename this Accession column to clarify between experiment_accession and file_accession.
samples <- read.table("/scratch/Shares/rinnclass/Group2/CLASS_2021/rnaseq/samples.txt",
                      sep = "\t", skip = 1, header = T) %>%
  dplyr::rename(experiment_accession = Accession) 
```

Now we are going to slowly build this samplesheet which will be used in the later RNAseq run. This requires the httr package.

``` r
# First we need a function  to query their api.
# A basic API query for a REST API is made in the form of an http request (just like when you visit a website)
# Instead of sending the information to your browser though, we'll just capture this information in text format.

# We'll use for this the R library httr
library(httr)
# The basic way to retrieve data from a rest API is with a GET request. 
(resp <- GET("http://httpbin.org/get"))
```

    ## Response [http://httpbin.org/get]
    ##   Date: 2021-04-23 20:14
    ##   Status: 200
    ##   Content-Type: application/json
    ##   Size: 366 B
    ## {
    ##   "args": {}, 
    ##   "headers": {
    ##     "Accept": "application/json, text/xml, application/xml, */*", 
    ##     "Accept-Encoding": "deflate, gzip", 
    ##     "Host": "httpbin.org", 
    ##     "User-Agent": "libcurl/7.29.0 r-curl/4.3 httr/1.4.2", 
    ##     "X-Amzn-Trace-Id": "Root=1-60832a8d-13cda40438e7bf0752e5e4eb"
    ##   }, 
    ##   "origin": "128.138.93.106", 
    ## ...

``` r
# This returns to us a response object which contains information about the request made as well as some accompanying data.
# In our case we want to query the ENCODE API using the experiment accession number and request information on the md5sum, number of reads, etc.
# The base url we'll use for querying ENCODE is as follows:
# https://www.encodeproject.org/report.tsv?
# You can see from this that we expect the data to be returned in a tsv -- tab separated values format.

# Let's look at an example request for this experiment accession: ENCSR541TIG
request_url <- "https://www.encodeproject.org/report.tsv?type=File&status=released&file_format=fastq&dataset=%2Fexperiments%2FENCSR541TIG%2F&field=accession&field=read_count&field=md5sum&field=controlled_by&field=paired_end&field=paired_with&field=replicate&field=target"
# Let's break this down into your parts.
# Just like a function, we have paramters and accomanying values separated by equal signs.
# https://www.encodeproject.org/report.tsv?type=File&status=released&file_format=fastq
# Here's where we request this particular accession
# dataset=%2Fexperiments%2FENCSR541TIG%2F

# And the field parameter is where we tell it which columns or which pieces of data we want to get.
# Here we're retrieving read_count, md5sum, controlled_by, paired_end, paired_with, replicate, and target
# field=read_count&field=md5sum&field=controlled_by&field=paired_end&field=paired_with&field=replicate&field=target


# Paste this in your browswer address bar and you'll get a preview of the data we're about to retrieve.
# https://www.encodeproject.org/report.tsv?type=File&status=released&file_format=fastq&dataset=%2Fexperiments%2FENCSR541TIG%2F&field=accession&field=read_count&field=md5sum&field=controlled_by&field=paired_end&field=paired_with&field=replicate&field=target

# Let's use the httr to retreive this same data into our r session
(resp <- GET(request_url))
```

    ## Response [https://www.encodeproject.org/report.tsv?type=File&status=released&file_format=fastq&dataset=%2Fexperiments%2FENCSR541TIG%2F&field=accession&field=read_count&field=md5sum&field=controlled_by&field=paired_end&field=paired_with&field=replicate&field=target]
    ##   Date: 2021-04-23 20:14
    ##   Status: 200
    ##   Content-Type: text/tsv; charset=UTF-8
    ##   Size: 900 B
    ## 2021-04-23 20:14:05.429030   https://www.encodeproject.org/report/?type=File&status=released&file_format=fastq&dataset=%2Fexperiments%2FEN...
    ## Accession    Read count  MD5sum  Controlled by   Paired end identifier   Paired with Replicate   Target
    ## ENCFF969RUJ  20236753    21f1d2d06b53bd6ee456c67e629cc734        1   /files/ENCFF448RSL/ /replicates/00781b21-22c3-4be2-ac4f-54422fd5da0e/   
    ## ENCFF714MGS  24235894    9ec27b60e022922c0bf34b6b14b2bf0f        1   /files/ENCFF675FVZ/ /replicates/6ea50902-8c42-4940-aed5-e529a8dcfbd4/   
    ## ENCFF448RSL  20236753    9541757b0d6b243e2d4cccbd0cc463e6        2   /files/ENCFF969RUJ/ /replicates/00781b21-22c3-4be2-ac4f-54422fd5da0e/   
    ## ENCFF675FVZ  24235894    152e103a1911c625d206e11d88179884        2   /files/ENCFF714MGS/ /replicates/6ea50902-8c42-4940-aed5-e529a8dcfbd4/   

``` r
# This object contains the data we want, but also contains header information about the request, so let's extract just the data (body) portion
body <- read_tsv(content(resp, "text"), skip = 1)
# We can see that this now gives us a data frame of just the requested information.
```

Now we are going to use some functions to retrieve the data from ENCODE required to make the samplesheet.

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
```

Now we will make a "matrix of matrices" that makes each cell of the sample sheet into its own embedded matrix that includes all of the information required to make the RNAseq design file.

``` r
samples <- samples %>%
  mutate(file_info = map(experiment_accession, ~ encode_file_info(.x)))

# Note that in this new column we have a bunch of data.frames -- nested in our samples data.frame.
# Let's retrieve one as an example
file_info <- samples$file_info[[1]]
```

To improve readability of the data frame that is being written, we can unnest some of the matrices from their nested cells.

``` r
# We also need to tell it which column we want to unnest.
samples <- samples %>%
  unnest(cols = file_info)
```

Now we have our completed sample sheet, we just need to clean this up and get it into the format that is required for nextflow.

``` r
samples <- read.table("/scratch/Shares/rinnclass/Group2/CLASS_2021/rnaseq/samples.txt",
                      sep = "\t", skip = 1, header = T) %>%
  dplyr::rename(experiment_accession = Accession) %>%
  mutate(file_info = map(experiment_accession, ~ encode_file_info(.x))) %>%
  unnest(file_info) %>% 
  group_by(experiment_accession) %>%
  mutate(rep_number = as.numeric(factor(replicate))) %>%
  unite(sample_id, experiment_accession, rep_number, sep = "_rep") %>%
  dplyr::select(sample_id, accession, Assay.title,
                Biosample.summary, md5sum, paired_end_identifier) %>%
  unite(fastq_file, sample_id, paired_end_identifier, sep = "_read", remove = F) %>%
  mutate(fq_extension = ".fastq.gz") %>%
  unite(fastq_file, fastq_file, fq_extension, sep = "", remove = F) %>%
  unite(original_file, accession, fq_extension, sep = "")
```

Now from this data frame, we can "write" the design file required for nf-seq (the nextflow pipeline for RNAseq).

``` r
# Rename the fastq files so that they contain the sample ID.
rename_script <- samples %>%
  dplyr::select(fastq_file, original_file) %>%
  mutate(command = "mv") %>%
  unite(command, command, original_file, fastq_file, sep = " ")
# The result of this is that each row is a bash command.

# We can write this out as a bash script with ?write_lines 
# We include a shebang header line so that the script is interpreted by bash.

write_lines(c("#!/bin/bash", rename_script$command), "/scratch/Shares/rinnclass/Group2/CLASS_2021/final_analysis/rnaseq/results/run.sh")
# Here we use an R command to call bash and cd into that directory, then make the file executable, and then run it.
# system("cd rnaseq/fastq; chmod u+x rename.sh; ./rename.sh")
```

Optionally, you can run the md5 checksums to make sure you have downloaded the correct data.

``` r
# Let's create an md5.txt to run the checksums
md5 <- samples %>% 
  dplyr::select(md5sum, fastq_file) %>%
  unite(line, md5sum, fastq_file, sep = "  ")
write_lines(md5$line, "/scratch/Shares/rinnclass/Group2/CLASS_2021/final_analysis/rnaseq/results/fastqc/md5.txt")
# Now let's run it.
# system("cd rnaseq/fastq; md5sum -c md5.txt")
```

Now we can write out everything we have done as a clean samplesheet.

``` r
# Let's create the sample sheet that we will use later
# to do the RNA-seq analysis in R.
samples <- samples %>%
  dplyr::rename(fastq = fastq_file,
                seq_type = Assay.title,
                sample_name = Biosample.summary) %>%
  # The minus sign will remove this column -- which we no longer need.
  dplyr::select(-original_file) 


# Now that we have it cleaned up, let's create one line for each replicate where the fastq read 1 and read 2 are in the same row.
# For this we will use the pivot wider function
# We need to tell the pivot_wider function which unique column combinations will specify each new row. 
samplesheet <- samples
samplesheet <- samplesheet %>%
  pivot_wider(id_cols = c("sample_id", "seq_type", "sample_name"),
              names_from = paired_end_identifier,
              values_from = c("fastq", "md5sum"))


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

write_csv(samplesheet, "/scratch/Shares/rinnclass/Group2/CLASS_2021/final_analysis/rnaseq/samplesheet.csv")
```

Finally, write out the final design file.

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
         fastq_1 = paste0("/scratch/Shares/rinnclass/data/rnaseq_fastq/", fastq_1),
         fastq_2 = paste0("/scratch/Shares/rinnclass/data/rnaseq_fastq/", fastq_2)) %>%
  # We have the replicate number already, but it's just a part of the sample id
  # We can use the separate function (opposite of unite) to retrieve it.
  separate(sample_id, into = c("experiment_accession", "replicate"), sep = "_", remove = FALSE) %>%
  mutate(replicate = gsub("rep", "", replicate)) %>%
  # Now we gather just the columns we need
  dplyr::select(group, replicate, fastq_1, fastq_2, strandedness)

design <- design %>%
  mutate(group = gsub("_r1", "", group),
         group = gsub("_r2", "", group))

design <- design %>% 
  arrange(group, replicate)

all(sapply(design$fastq_1, file.exists))
```

    ## [1] TRUE

``` r
all(sapply(design$fastq_2, file.exists))
```

    ## [1] TRUE

``` r
# Now we can write this out for input into nextflow
write_csv(design, "/scratch/Shares/rinnclass/Group2/CLASS_2021/final_analysis/rnaseq/design.csv")
```

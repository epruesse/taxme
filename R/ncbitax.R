library(data.table)

download.ncbitax <- function() {
  url <- "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
  fname <- "taxdump.tar.gz"
  tmpdir <- paste0(tempdir(), "/ncbitax/")
  dir.create(tmpdir)
  download.file(paste0(url, fname), paste0(tmpdir, fname))
  system(paste0("tar -C ", tmpdir, " -xzf ", paste0(tmpdir,fname)))
  tmpdir
}

read.ncbi <- function(dir, nrows=-1, use.cache=TRUE, clear.cache=FALSE) {
  if (length(list.files(pattern="ncbitax.Rds"))) {
    if (clear.cache) {
      unlink("ncbitax.Rds")
    } else if (use.cache) {
      return (readRDS("ncbitax.Rds"))
    }
  }
  if (missing(dir)) {
    dir <- download.ncbitax()
    unlink_dir <- TRUE
  } else {
    unlink_dir <- FALSE
  }
  warning("reading nodes...")
  nodes <- read.delim(paste0(dir, "/nodes.dmp"),
                      nrows=nrows,
                      sep="|",
                      strip.white=TRUE,
                      col.names = c('tax_id', 'parent_id','rank',1:11),
                      header=FALSE) %>%
    select(tax_id, parent_id, rank)

  warning("reading merged nodes...")
  merged <- read.delim(paste0(dir, "/merged.dmp"),
                       nrows=nrows,
                       sep="|",
                       strip.white=TRUE,
                       stringsAsFactors=FALSE,
                       col.names = c('tax_id', 'parent_id','rank',1),
                       header=FALSE) %>%
    select(tax_id, parent_id) %>%
    mutate(rank = "None")

  warning("reading names...")
  names <- read.delim(paste0(dir, "/names.dmp"),
                      nrows=nrows,
                      sep="|",
                      strip.white=TRUE,
                      quote="",
                      stringsAsFactors=FALSE,
                      col.names= c('tax_id', 'name', 'uname', 'name_cls',1),
                      header=FALSE) %>%
    filter(name_cls=="scientific name")  %>%
    select(tax_id, name)

  if (unlink_dir) {
    unlink(dir, recursive=TRUE)
  }

  nodes <- rbind(merged, nodes)
  nodes <- nodes[order(nodes$tax_id), ]
  ncbi <- merge(nodes, names, all.x=True) %>%
    mutate(rank = as.factor(rank))
  ncbi <- data.table(ncbi)
  setkey(ncbi, tax_id)

  if (use.cache) {
    saveRDS(ncbi, file="ncbitax.Rds")
  }
  assign(".ncbi", ncbi, envir = .GlobalEnv)
  ncbi
}


ncbi.lineage <- function(tax_id, ncbi=.ncbi) {
  lids = list(tax_id)
  setDT(lids, "tax_id")
  lin = list()
  lin <- append(lin, lids)
  while(max(lids) != 1) {
    parents <- ncbi[lids, parent_id]
    lids <- list(parents)
    setDT(lids, "tax_id")
    lin <- append(lin, lids)
  }
  data.table::transpose(lin)
}

ncbi.path <- function(t_tax_id, ncbi=.ncbi) {
  lapply(ncbi.lineage(t_tax_id),
         function(x) ncbi[list(x), name, rank])
}


norm_ranks = factor(c("species", "genus", "family", "order", "class", "phylum", "superkingdom"),
                    levels=levels(ncbi$rank))

ncbi.norm_path <- function(t_tax_id, ncbi=.ncbi) {
  lapply(ncbi.lineage(t_tax_id),
         function(x) ncbi[list(x)][list(norm_ranks), on="rank", name, rank])
}

ncbi.get_rank <- function(t_tax_id, t_rank="superkingdom", ncbi=.ncbi) {
  ranks <- ncbi[list("superkingdom"), on="rank", tax_id, name]
  taxranks <- unlist(sapply(
    ncbi.lineage(t_tax_id),
    function(x) { n = x %in% ranks$tax_id; if (length(n)) x[n] else 1 }
  ))
  ranks[list(taxranks),on="tax_id", name]
}

ncbi.group <- function(taxid, ranks=c(), nodes=c(), ncbi=.ncbi) {
  terminal = ncbi[1]
  if (length(ranks) > 0) {
    terminal <- rbind(ncbi[list(ranks), on="rank"], terminal)
  }
  if (length(nodes) > 0) {
    terminal <- rbind(ncbi[list(nodes), on="name"], terminal)
  }

  x<-ncbi.lineage(taxid)
  x<-lapply(x, match, terminal$tax_id)
  x<-sapply(x, detect, function(y) !is.na(y))
  terminal[x]$name
}

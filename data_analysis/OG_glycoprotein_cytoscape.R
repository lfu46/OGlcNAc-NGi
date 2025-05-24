#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "RCy3")
lapply(packages_names, require, character.only = TRUE)

#test
nodes <- data.frame(id=c("node 0","node 1","node 2","node 3"),
                    group=c("A","A","B","B"), # categorical strings
                    score=as.integer(c(20,10,15,5)), # integers
                    stringsAsFactors=FALSE)
edges <- data.frame(source=c("node 0","node 0","node 0","node 2"),
                    target=c("node 1","node 2","node 3","node 3"),
                    interaction=c("inhibits","interacts","activates","interacts"),  # optional
                    weight=c(5.1,3.0,5.2,9.9), # numeric
                    stringsAsFactors=FALSE)

#string
cytoscapePing()

string_cmd_1 <- paste0("string protein query query='", paste(OG_glycoprotein_Top_tb_HepG2$UniprotID, collapse = ","), "'")
commandsRun(string_cmd_1)

string_cmd_2 <- 'string retrieve enrichment background="genome" selectedNodesOnly=false'
commandsRun(string_cmd_2)

string_cmd_3 <- 'string show enrichment'
commandsRun(string_cmd_3)

string_cmd_4 <- 'string filter enrichment categories="InterPro Domains,STRING Clusters" removeOverlapping=true'
commandsRun(string_cmd_4)

closeSession(save.before.closing = FALSE)

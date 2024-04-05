library(pathview)

do_pathview <- function(
    path.id,
    genes.list) {
  genes.list <- as.vector(strsplit(genes.list, ",")[[1]])

  # convert.genes <- mapIds(org.Mm.eg.db, keys = genes.list,
  #      column = "ENTREZID", keytype = "SYMBOL")
  # names(convert.genes) <- NULL

  pv.out <- pathview(gene.data = genes.list, pathway.id = as.character(path.id), species = "mmu", same.layer = F, out.suffix = "pathview", gene.idtype = "ensembl")

  return(pv.out)
}

Pathview.ui <- tabItem(
  "Pathview",
  fluidRow(
    sidebarPanel(
      selectizeInput(
        inputId = "selectedGenes",
        label = "Please enter your genes names one by one into the tab below",
        choices = NULL,
        selected = NULL,
        multiple = TRUE,
        width = "100%",
        options = list(
          "plugins" = list("remove_button"),
          "create" = TRUE,
          "persist" = TRUE
        )
      )
    ),
    mainPanel(textOutput("pathview.genes"))
  ),
  fluidRow(
    sidebarPanel(
      selectizeInput(
        inputId = "selectedPath",
        label = "Please enter your pathway ID of interest",
        choices = NULL,
        selected = NULL,
        multiple = TRUE,
        width = "100%",
        options = list(
          "plugins" = list("remove_button"),
          "create" = TRUE,
          "persist" = TRUE
        )
      )
    ),
    mainPanel(textOutput("pathview.id"))
  ),
  fluidRow(box(actionButton("showPathview", "Show Pathview", width = 200), footer = "Click on the button to generate the pathview results for a given set of genes and an enrichment pathway.")),
  fluidRow(
    box(imageOutput("pathview", height = 1000, fill = TRUE), width = 12)
  ),
  fluidRow(
    box(imageOutput("pathview2", height = 1000, fill = TRUE), width = 12)
  )
)

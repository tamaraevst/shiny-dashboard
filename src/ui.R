## ui.R ##

# Ui part of the Shiny app

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Your Data", tabName = "yourdata", icon = icon("dashboard")),
    menuItem("PCA Analysis on Samples", tabName = "PCAsamples"),
    menuItem("PCA Analysis on Genes", tabName = "PCAplot"),
    menuItem("Enrichment Analysis", tabName = "Enrichment"),
    menuItem("Pathview Analysis", tabName = "Pathview")
  )
)

body <- dashboardBody(
  tabItems(
    # First tab content
    tabItem(
      tabName = "yourdata",
      # fluidRow(
      #    box(
      #       fileInput('file1', 'Choose CSV File',
      #       accept=c('text/csv',
      #                  'text/comma-separated-values,text/plain',
      #                  '.csv')),
      #       width=15
      #       )
      #         ),
      fluidRow(
        box(
          DT::dataTableOutput("contents"),
          style = "overflow-y: scroll;overflow-x: scroll;",
          width = 15
        )
      )
    ),
    PCA.samples.ui,
    PCA.ui,
    Enrichment.ui,
    Pathview.ui
  )
)

ui <- dashboardPage(
  dashboardHeader(title = "Dashboard for omics data"),
  sidebar,
  body,
  skin = "purple"
)

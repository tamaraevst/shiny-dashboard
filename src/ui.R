## ui.R ##

# Ui part of the Shiny app

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Your Data", tabName = "yourdata", icon = icon("dashboard")),
    menuItem("PCA Analysis", tabName = "PCAplot"),
    menuItem("Enrichment Analysis", tabName = "Enrichment")
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
          width = 15
        )
      )
    ),
    PCA.ui,
    Enrichment.ui
  )
)

ui <- dashboardPage(
  dashboardHeader(title = "Dashboard for omics data"),
  sidebar,
  body,
  skin = "purple"
)

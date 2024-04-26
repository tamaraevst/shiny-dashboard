## data.R ##

Data.server <- function(id, expression.matrix) {
  moduleServer(
    id,
    function(input, output, session) {
      output$contents <- DT::renderDataTable({
        DT::datatable(expression.matrix)
      })
    }
  )
}

Data.ui <- function(id) {
  tabItem(
    tabName = "yourdata",
    fluidRow(
      box(
        DT::dataTableOutput(NS(id, "contents")),
        style = "overflow-y: scroll;overflow-x: scroll;",
        width = 15
      )
    )
  )
}

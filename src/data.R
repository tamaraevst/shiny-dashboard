## data.R ##

Data.server <- function(id, expression.matrix) {
  moduleServer(
    id,
    function(input, output, session) {
      expression.matrix.filtered <- expression.matrix[, -c(1, 2)] %>%
        relocate(NAME)

      output$contents <- DT::renderDataTable({
        DT::datatable(expression.matrix.filtered)
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

library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(ggrepel)
library(mvtnorm)

options(shiny.maxRequestSize = 600 * 1024^2)

ui <- fluidPage(
  titlePanel("SexPeptID"),
  
  sidebarPanel(
    
    wellPanel(
      tags$h4("Select Species:"),
      selectizeInput(
        "species", 
        label = NULL, 
        choices = c(
          "Homo sapiens & hominins (Q99217/Q99218)"
        ), 
        selected = "Homo sapiens & hominins (Q99217/Q99218)",
        multiple = FALSE
      )
    ),
    
    wellPanel(
      tags$h4("Upload your data:"),
      fileInput("peptides", label = h5("peptides.txt"), accept = ".txt"),
      fileInput("msms", label = h5("msms.txt"), accept = ".txt")
    ),
    
    wellPanel(
      tags$h4("Graph"),
      downloadButton('downloadGraph', "Download the Graph")
    ),
    
    wellPanel(
      tags$h4("Table"),
      downloadButton('downloadTable', "Download the Table")
    )
  ),
  
  mainPanel(
    tabsetPanel(
      type = "tabs",
      
      tabPanel(
        "Graph and Table",
        plotOutput("plot1"),
        DTOutput("table1")
      ),
      
      tabPanel(
        "Peptide list",
        selectInput("raw_file_select", "Select Raw file:", choices = NULL),
        DTOutput("peptide_table")
      )
    )
  )
)

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    peptide_list = NULL,
    peptides_raw = NULL,
    checkbox_states = list()
  )
  
  rv_table <- reactiveVal(NULL)
  
  # ------------------- Load peptides.txt -------------------
  
  observeEvent(input$peptides, {
    
    req(input$peptides)
    species <- input$species
    
    data <- read.delim(input$peptides$datapath, sep = "\t", header = TRUE)
    
    if(species == "Homo sapiens & hominins (Q99217/Q99218)") {
      
      peptides_filter <- data[, c("Proteins","Sequence","Start.position","End.position")]
      
      AMELY <- peptides_filter %>%
        filter(grepl("Q99218", Proteins) & !grepl("Q99217", Proteins)) %>%
        mutate(SMI = 45) %>%
        filter(between(SMI, Start.position, End.position))
      
      AMELX <- peptides_filter %>%
        filter(grepl("Q99217", Proteins)) %>%
        mutate(SIR = 44, IR = 45) %>%
        filter(between(IR, Start.position, End.position)) %>%
        filter(between(SIR, Start.position, End.position))
      
      rv$peptides_raw <- list(AMELX = AMELX, AMELY = AMELY)
      
    } else {
      rv$peptides_raw <- list(AMELX = NULL, AMELY = NULL)
      showNotification("Currently only Homo sapiens & hominins is supported.", type = "warning")
    }
  })
  
  # ------------------- Load msms.txt and merge -------------------
  
  observeEvent(input$msms, {
    
    req(input$msms, rv$peptides_raw)
    
    msms <- read.delim(input$msms$datapath, sep = "\t", header = TRUE)
    
    msms_filter <- msms[, c("Sequence","Precursor.Intensity","Raw.file","Proteins","PEP")]
    matches_col <- msms[, c("Sequence","Proteins","Matches")]
    
    AMELX <- merge(rv$peptides_raw$AMELX, msms_filter, by = c("Sequence","Proteins"))
    AMELY <- merge(rv$peptides_raw$AMELY, msms_filter, by = c("Sequence","Proteins"))
    
    AMELX <- merge(AMELX, matches_col, by = c("Sequence","Proteins"))
    AMELY <- merge(AMELY, matches_col, by = c("Sequence","Proteins"), all.x = TRUE)
    
    if(nrow(AMELY) == 0) {
      showNotification("No AMELY peptides detected. Creating placeholder.", type = "warning")
      AMELY <- data.frame(
        Sequence = character(), Raw.file = character(), Proteins = character(),
        Start.position = numeric(), End.position = numeric(),
        Matches = character(), Precursor.Intensity = numeric(), PEP = numeric()
      )
    }
    
    # ------------------- Count b/y ions -------------------
    count_ions <- function(m, t) {
      if(is.na(m) || m == "") return(0)
      sum(grepl(paste0("^", t, "\\d+$"), unlist(strsplit(m, ";"))))
    }
    
    AMELX <- AMELX %>% rowwise() %>%
      mutate(`b ions` = paste0(count_ions(Matches,"b"), "/", nchar(Sequence)),
             `y ions` = paste0(count_ions(Matches,"y"), "/", nchar(Sequence))) %>% ungroup()
    
    AMELY <- AMELY %>% rowwise() %>%
      mutate(`b ions` = paste0(count_ions(Matches,"b"), "/", nchar(Sequence)),
             `y ions` = paste0(count_ions(Matches,"y"), "/", nchar(Sequence))) %>% ungroup()
    
    # ------------------- Aggregate -------------------
    AMELX <- AMELX %>% group_by(Sequence, Raw.file, Proteins, Start.position, End.position, PEP, `b ions`, `y ions`) %>%
      summarise(Precursor.Intensity = sum(Precursor.Intensity), .groups="drop")
    AMELY <- AMELY %>% group_by(Sequence, Raw.file, Proteins, Start.position, End.position, PEP, `b ions`, `y ions`) %>%
      summarise(Precursor.Intensity = sum(Precursor.Intensity), .groups="drop")
    
    peptide_list <- rbind(AMELX, AMELY)
    peptide_list$Use <- TRUE
    
    rv$peptide_list <- peptide_list
    rv$checkbox_states <- peptide_list %>% split(.$Raw.file) %>% 
      lapply(function(df) setNames(df$Use, seq_len(nrow(df))))
    
    updateSelectInput(session, "raw_file_select", choices = unique(peptide_list$Raw.file))
  })
  
  # ------------------- Peptide checkbox table -------------------
  peptide_proxy <- dataTableProxy("peptide_table")
  
  output$peptide_table <- renderDT(server = TRUE, {
    req(rv$peptide_list, input$raw_file_select)
    selected_file <- input$raw_file_select
    df <- rv$peptide_list %>% filter(Raw.file == selected_file)
    
    if(is.null(rv$checkbox_states[[selected_file]])){
      rv$checkbox_states[[selected_file]] <- setNames(df$Use, seq_len(nrow(df)))
    }
    
    df_display <- df
    df_display$Add <- sapply(seq_len(nrow(df_display)), function(i){
      id <- paste0("chk_", selected_file, "_", i)
      checked <- ifelse(isTRUE(rv$checkbox_states[[selected_file]][i]), "checked", "")
      sprintf('<input type="checkbox" id="%s" %s onclick="Shiny.setInputValue(\'%s\', this.checked, {priority: \'event\'})">', 
              id, checked, id)
    })
    
    # ✅ Utiliser la sélection classique par nom de vecteur
    df_display <- df_display[, c("Add", "Raw.file", "Sequence", "Proteins", "Start.position", "End.position", "b ions", "y ions")]
    
    datatable(
      df_display,
      escape = FALSE, 
      rownames = FALSE,
      options = list(scrollX=TRUE)
    )
  })
  
  # ------------------- Update checkbox states -------------------
  observe({
    req(rv$peptide_list, input$raw_file_select)
    selected_file <- input$raw_file_select
    df <- rv$peptide_list %>% filter(Raw.file == selected_file)
    n <- nrow(df)
    
    for(i in seq_len(n)){
      id <- paste0("chk_", selected_file, "_", i)
      if(!is.null(input[[id]])){
        rv$checkbox_states[[selected_file]][i] <- input[[id]]
      }
    }
  })
  
  # ------------------- Build summary table -------------------
  observe({
    req(rv$peptide_list)
    all_raw <- unique(rv$peptide_list$Raw.file)
    
    df_list <- lapply(all_raw, function(raw_file) {
      df <- rv$peptide_list %>% filter(Raw.file == raw_file)
      chk <- rv$checkbox_states[[raw_file]]
      df$Use <- chk
      df[df$Use == TRUE, ]
    })
    
    df <- do.call(rbind, df_list)
    if(nrow(df) == 0) { rv_table(data.frame()); return() }
    
    AMELX <- df[grepl("Q99217", df$Proteins), ]
    AMELY <- df[grepl("Q99218", df$Proteins), ]
    
    AMELX_int <- aggregate(Precursor.Intensity ~ Raw.file + Proteins, AMELX, sum)
    AMELY_int <- aggregate(Precursor.Intensity ~ Raw.file + Proteins, AMELY, sum)
    
    names(AMELX_int)[3] <- "Precursor.Intensity.AMELX"
    names(AMELY_int)[3] <- "Precursor.Intensity.AMELY"
    
    df_sum <- merge(AMELX_int, AMELY_int, by="Raw.file", all=TRUE)
    df_sum$`Protein X` <- df_sum$Proteins.x
    df_sum$`Protein Y` <- df_sum$Proteins.y
    df_sum$logAMELX <- log(df_sum$Precursor.Intensity.AMELX)
    df_sum$logAMELY <- ifelse(!is.na(df_sum$Precursor.Intensity.AMELY), log(df_sum$Precursor.Intensity.AMELY), 0)
    
    df_sum$`Biological sex` <- ifelse(is.na(df_sum$Precursor.Intensity.AMELY) | df_sum$Precursor.Intensity.AMELY==0,"Female","Male")
    male_idx <- which(df_sum$`Biological sex`=="Male")
    if(length(male_idx) > 1){
      mu <- mean(df_sum$logAMELY[male_idx])
      sigma <- sd(df_sum$logAMELY[male_idx])
      non_conclusive <- male_idx[dnorm(df_sum$logAMELY[male_idx], mean=mu, sd=sigma) < 0.05]
      df_sum$`Biological sex`[non_conclusive] <- "Non-conclusive"
    }
    if(!"Label" %in% names(df_sum)) df_sum$Label <- ""
    rv_table(df_sum)
  })
  
  # ------------------- Table1 output -------------------
  output$table1 <- renderDT({
    df <- rv_table()
    req(df)
    df$logAMELX <- sprintf("%.3f", df$logAMELX)
    df$logAMELY <- sprintf("%.3f", df$logAMELY)
    df <- df[, c("Raw.file","Label","Protein X","Protein Y","logAMELX","logAMELY","Biological sex")]
    
    datatable(df, rownames=FALSE, editable=list(target="cell", disable=list(columns=c(0,2,3,4,5,6))), options=list(scrollX=TRUE))
  })
  
  observeEvent(input$table1_cell_edit, {
    info <- input$table1_cell_edit
    df <- rv_table()
    req(df)
    if(info$col == 1) df$Label[info$row] <- info$value
    rv_table(df)
  })
  
  # ------------------- Plot -------------------
  output$plot1 <- renderPlot({
    df <- rv_table()
    req(df)
    ggplot(df, aes(logAMELX, logAMELY, colour=`Biological sex`)) +
      geom_point(size=4) +
      geom_text_repel(aes(label=Label, colour=`Biological sex`), size=4, max.overlaps=Inf, box.padding=0.5, point.padding=0.5, segment.size=0.7, min.segment.length=0, show.legend=FALSE) +
      geom_smooth(data=subset(df, `Biological sex`=="Male"), aes(logAMELX, logAMELY), method="lm", color="black", fill="gray", se=TRUE) +
      scale_color_manual(values=c("Female"="#fec000","Male"="#5b9cd4","Non-conclusive"="#ff3b3b")) +
      theme_minimal(base_size=16)
  })
  
  # ------------------- Downloads -------------------
  output$downloadGraph <- downloadHandler(
    filename = function() paste0("graph-DDA-SexID-", Sys.Date(), ".tiff"),
    content = function(file) ggsave(file, plot=last_plot(), width=12, height=8, device="tiff")
  )
  
  output$downloadTable <- downloadHandler(
    filename = function() paste0("table-DDA-SexID-", Sys.Date(), ".csv"),
    content = function(file) write.csv(rv_table(), file, row.names=FALSE)
  )
}

shinyApp(ui, server)

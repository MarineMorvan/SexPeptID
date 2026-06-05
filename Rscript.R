library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(ggrepel)


options(shiny.maxRequestSize = 800 * 1024^2)

ui <- fluidPage(
  
  titlePanel("SexPeptID"),
  
  sidebarPanel(
    
    wellPanel(
      tags$h4("Select Species:"),
      selectizeInput(
        "species",
        label = NULL,
        choices = c(
          "Homo sapiens & hominins (Q99217/Q99218)",
          "Bos taurus (P02817/Q99004)"
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
      downloadButton("downloadGraph", "Download the Graph")
    ),
    
    wellPanel(
      tags$h4("Table"),
      downloadButton("downloadTable", "Download the Table")
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
    checkbox_states = list(),
    protein_x = NULL,
    protein_y = NULL
  )
  
  rv_table <- reactiveVal(NULL)
  
  # ------------------- Load peptides.txt -------------------
  
  observeEvent(input$peptides, {
    
    req(input$peptides)
    
    species <- input$species
    
    data <- read.delim(input$peptides$datapath, sep = "\t", header = TRUE)
    
    required_cols <- c("Proteins","Sequence","Start.position","End.position")
    
    validate(
      need(all(required_cols %in% names(data)),
           paste("Missing required columns:",
                 paste(setdiff(required_cols, names(data)), collapse = ", ")))
    )
    
    peptides_filter <- data %>%
      dplyr::select(all_of(required_cols))
    
    extract_peptides <- function(data, protein_id, positions, exclude_id = NULL) {
      
      result <- data %>% filter(grepl(protein_id, Proteins))
      
      if (!is.null(exclude_id)) {
        result <- result %>% filter(!grepl(exclude_id, Proteins))
      }
      
      for (pos in positions) {
        result <- result %>% filter(Start.position <= pos & End.position >= pos)
      }
      
      result
    }

  # ------------------- Species logic -------------------

  extract_peptides <- function(data,
                             protein_id,
                             positions,
                             residues = NULL,
                             exclude_id = NULL) {
  
    result <- data %>%
      filter(grepl(protein_id, Proteins))
  
    if (!is.null(exclude_id)) {
      result <- result %>%
        filter(!grepl(exclude_id, Proteins))
    }
  
    for (pos in positions) {
      result <- result %>%
        filter(Start.position <= pos & End.position >= pos)
    }
  
    if (!is.null(residues)) {
      for (i in seq_along(positions)) {
        pos <- positions[i]
        aa  <- residues[i]
      
        result <- result %>%
          rowwise() %>%
          filter(
            substr(
              Sequence,
              pos - Start.position + 1,
              pos - Start.position + 1
            ) == aa
          ) %>%
          ungroup()
      }
    }
  
    result
  }

  if (species == "Homo sapiens & hominins (Q99217/Q99218)") {
  
    rv$protein_x <- "Q99217"
    rv$protein_y <- "Q99218"
    
    AMELX_44 <- extract_peptides(
      peptides_filter,
      protein_id = "Q99217",
      positions  = c(44, 45),
      residues   = c("S", "I")
    )
  
    AMELX_28 <- extract_peptides(
      peptides_filter,
      protein_id = "Q99217",
      positions  = c(28, 29),
      residues   = c("S", "I")
    )
  
    AMELX_58 <- extract_peptides(
      peptides_filter,
      protein_id = "Q99217",
      positions  = c(58, 59),
      residues   = c("S", "I")
    )
  
    AMELX <- bind_rows(
      AMELX_44,
      AMELX_28,
      AMELX_58
    ) %>%
      distinct()
  
    AMELY_45 <- extract_peptides(
      peptides_filter,
      protein_id = "Q99218",
      positions  = 45,
      residues   = "M",
      exclude_id = "Q99217"
    )
  
    AMELY_59 <- extract_peptides(
      peptides_filter,
      protein_id = "Q99218",
      positions  = 59,
      residues   = "M",
      exclude_id = "Q99217"
   )
  
    AMELY <- bind_rows(
      AMELY_45,
      AMELY_59
    ) %>%
      distinct()
  
    rv$peptides_raw <- list(
      AMELX = AMELX,
      AMELY = AMELY
    )
  
    showNotification(
      paste(
        nrow(AMELX), "AMELX peptides detected |",
        nrow(AMELY), "AMELY peptides detected"
      ),
      type = "message"
    )
  
  } else if (species == "Bos taurus (P02817/Q99004)") {
  
    rv$protein_x <- "P02817"
    rv$protein_y <- "Q99004"
  
    AMELX_44 <- extract_peptides(
      peptides_filter,
      protein_id = "P02817",
      positions  = 44,
      exclude_id = "Q99004"
    )
  
    AMELX_48 <- extract_peptides(
      peptides_filter,
      protein_id = "P02817",
      positions  = 48,
      exclude_id = "Q99004"
    )
  
    AMELX <- bind_rows(
      AMELX_44,
      AMELX_48
    ) %>%
      distinct()
  
    AMELY_44 <- extract_peptides(
      peptides_filter,
      protein_id = "Q99004",
      positions  = 44,
      exclude_id = "P02817"
    )
  
    AMELY_48 <- extract_peptides(
      peptides_filter,
      protein_id = "Q99004",
      positions  = 48,
      exclude_id = "P02817"
    )
  
    AMELY <- bind_rows(
      AMELY_44,
      AMELY_48
    ) %>%
      distinct()
  
    rv$peptides_raw <- list(
      AMELX = AMELX,
      AMELY = AMELY
    )
  
    showNotification(
      paste(
        nrow(AMELX), "AMELX peptides detected |",
        nrow(AMELY), "AMELY peptides detected"
      ),
      type = "message"
    )
  
  } else {
  
    rv$peptides_raw <- list(
      AMELX = NULL,
      AMELY = NULL
    )
  
    showNotification(
      "Currently only Homo sapiens & hominins and Bos taurus are supported.",
      type = "warning"
    )
  }
    })
  
  # ------------------- Load msms.txt -------------------
  
  observeEvent(input$msms, {
    
    req(input$msms, rv$peptides_raw)
    
    msms <- read.delim(input$msms$datapath, sep = "\t", header = TRUE)
    
    msms_filter <- msms[, c("Sequence","Precursor.Intensity","Raw.file","Proteins","PEP")]
    matches_col <- msms[, c("Sequence","Proteins","Matches")]
    
    AMELX <- merge(rv$peptides_raw$AMELX, msms_filter, by = c("Sequence","Proteins"))
    AMELY <- merge(rv$peptides_raw$AMELY, msms_filter, by = c("Sequence","Proteins"))
    
    AMELX <- merge(AMELX, matches_col, by = c("Sequence","Proteins"))
    AMELY <- merge(AMELY, matches_col, by = c("Sequence","Proteins"), all.x = TRUE)
    
    if (nrow(AMELY) == 0) {
      AMELY <- data.frame(
        Sequence = character(),
        Raw.file = character(),
        Proteins = character(),
        Start.position = numeric(),
        End.position = numeric(),
        Matches = character(),
        Precursor.Intensity = numeric(),
        PEP = numeric()
      )
    }

    # ------------------- Count b/y ions -------------------
    
    count_ions <- function(m,t){
      if(is.na(m)||m=="") return(0)
      sum(grepl(paste0("^",t,"\\d+$"), unlist(strsplit(m,";"))))
    }
    
    AMELX <- AMELX %>%
      rowwise() %>%
      mutate(
        `b ions` = paste0(count_ions(Matches,"b"),"/",nchar(Sequence)),
        `y ions` = paste0(count_ions(Matches,"y"),"/",nchar(Sequence))
      ) %>% ungroup()
    
    AMELY <- AMELY %>%
      rowwise() %>%
      mutate(
        `b ions` = paste0(count_ions(Matches,"b"),"/",nchar(Sequence)),
        `y ions` = paste0(count_ions(Matches,"y"),"/",nchar(Sequence))
      ) %>% ungroup()
    
    AMELX <- AMELX %>%
      group_by(Sequence,Raw.file,Proteins,Start.position,End.position,PEP,`b ions`,`y ions`) %>%
      summarise(Precursor.Intensity=sum(Precursor.Intensity),.groups="drop")
    
    AMELY <- AMELY %>%
      group_by(Sequence,Raw.file,Proteins,Start.position,End.position,PEP,`b ions`,`y ions`) %>%
      summarise(Precursor.Intensity=sum(Precursor.Intensity),.groups="drop")
    
    peptide_list <- rbind(AMELX, AMELY)

    # ------------------- Auto-uncheck PEP > 0.05 -------------------
    
    peptide_list$Use <- ifelse(peptide_list$PEP > 0.05, FALSE, TRUE)
    
    rv$peptide_list <- peptide_list
    
    rv$checkbox_states <- peptide_list %>%
      split(.$Raw.file) %>%
      lapply(function(df) setNames(df$Use, seq_len(nrow(df))))
    
    updateSelectInput(session, "raw_file_select",
                       choices = unique(peptide_list$Raw.file))
  })
  
  # ------------------- Print peptide counts in console -------------------

  observeEvent(rv$peptide_list, {

    req(rv$peptide_list)

    cat("\n--- Peptide counts per protein per Raw file ---\n")

    total_counts <- rv$peptide_list %>%
      group_by(Raw.file, Proteins) %>%
      summarise(total_peptides = n(), .groups = "drop")

    print(total_counts, n = Inf)

    cat("\n--- Peptide counts with PEP < 0.05 per protein per Raw file ---\n")

    retained_counts <- rv$peptide_list %>%
      filter(PEP < 0.05) %>%
      group_by(Raw.file, Proteins) %>%
      summarise(retained_peptides = n(), .groups = "drop")

    print(retained_counts, n = Inf)

  })

  # ------------------- Peptide checkbox table -------------------
  
  peptide_proxy <- dataTableProxy("peptide_table")
  
  output$peptide_table <- renderDT({
    
    req(rv$peptide_list, input$raw_file_select)
    
    selected_file <- input$raw_file_select
    
    df <- rv$peptide_list %>%
      filter(Raw.file == selected_file)
    
    df_display <- df
    
    df_display$Add <- sapply(seq_len(nrow(df_display)), function(i) {
      id <- paste0("chk_", selected_file, "_", i)
      checked <- ifelse(isTRUE(rv$checkbox_states[[selected_file]][i]), "checked", "")
      
      sprintf(
        '<input type="checkbox" id="%s" %s onclick="Shiny.setInputValue(\'%s\', this.checked, {priority: \'event\'})">',
        id, checked, id
      )
    })
    
    datatable(df_display[, c("Add","Raw.file","Sequence","Proteins",
                             "Start.position","End.position","PEP",
                             "b ions","y ions")],
              escape = FALSE, rownames = FALSE)
  })
  
  
  observe({
    
    req(rv$peptide_list, input$raw_file_select)
    
    selected_file <- input$raw_file_select
    df <- rv$peptide_list %>% filter(Raw.file == selected_file)
    
    for (i in seq_len(nrow(df))) {
      id <- paste0("chk_", selected_file, "_", i)
      
      if (!is.null(input[[id]])) {
        rv$checkbox_states[[selected_file]][i] <- input[[id]]
      }
    }
  })
  
  # ------------------- Build summary table -------------------

  summary_table <- reactive({

    req(rv$peptide_list)

  raw_status <- rv$peptide_list %>%
    group_by(Raw.file) %>%
    summarise(
      n_AMELX_raw = sum(grepl(rv$protein_x, Proteins)),
      n_AMELY_raw = sum(grepl(rv$protein_y, Proteins)),
      .groups = "drop"
    )

  df <- rv$peptide_list %>%
    filter(Use == TRUE)

  if (nrow(df) == 0) {

    return(
      raw_status %>%
        mutate(
          AMELX = 0,
          AMELY = 0,
          logAMELX = 0,
          logAMELY = 0,
          `P(male)` = 0,
          `Biological sex` = case_when(
            n_AMELY_raw == 0 & n_AMELX_raw > 0 ~ "Female",
            n_AMELY_raw > 0 ~ "Non-conclusive",
            TRUE ~ "Non-conclusive"
          ),
          Label = ""
        )
    )
  }

  AMELX <- df %>%
    filter(grepl(rv$protein_x, Proteins))

  AMELY <- df %>%
    filter(grepl(rv$protein_y, Proteins))

  status <- df %>%
    group_by(Raw.file) %>%
    summarise(

      AMELX_peptides = sum(grepl(rv$protein_x, Proteins)),
      AMELY_peptides = sum(grepl(rv$protein_y, Proteins)),

      AMELX_has_signal =
        any(is.finite(Precursor.Intensity[grepl(rv$protein_x, Proteins)])),

      AMELY_has_signal =
        any(is.finite(Precursor.Intensity[grepl(rv$protein_y, Proteins)])),

      .groups = "drop"
    )

    AMELX_int <- AMELX %>%
      group_by(Raw.file) %>%
      summarise(
        AMELX = sum(
          Precursor.Intensity[is.finite(Precursor.Intensity)],
          na.rm = TRUE
        ),
        .groups = "drop"
      )

    AMELY_int <- AMELY %>%
      group_by(Raw.file) %>%
      summarise(
        AMELY = sum(
          Precursor.Intensity[is.finite(Precursor.Intensity)],
          na.rm = TRUE
        ),
        .groups = "drop"
      )

    df_sum <- raw_status %>%
      full_join(AMELX_int, by = "Raw.file") %>%
      full_join(AMELY_int, by = "Raw.file") %>%
      full_join(status, by = "Raw.file")

    df_sum[is.na(df_sum)] <- 0

    eps <- 1

    df_sum <- df_sum %>%
      mutate(
        logAMELX = log(AMELX + eps),
        logAMELY = log(AMELY + eps)
      )

    k <- 0.000001

    df_sum <- df_sum %>%
      mutate(
        `P(male)` = 1 - exp(-k * AMELY)
      )

    df_sum <- df_sum %>%
      mutate(
        `Biological sex` = case_when(

          n_AMELY_raw == 0 & n_AMELX_raw > 0 ~ "Female",

          n_AMELY_raw > 0 & AMELY == 0 ~ "Non-conclusive",

          (AMELX_peptides > 0 & !AMELX_has_signal) |
            (AMELY_peptides > 0 & !AMELY_has_signal) ~ "Non-conclusive",

          `P(male)` > 0.9 ~ "Male",

          TRUE ~ "Non-conclusive"
        ),

        Label = ""
      )

    df_sum
  })

  observe({
    rv_table(summary_table())
  })
  
  # ------------------- Table1 output -------------------
  
  output$table1 <- renderDT({
    
    df <- rv_table()
    req(df)
    
    df$logAMELX <- sprintf("%.3f", df$logAMELX)
    df$logAMELY <- sprintf("%.3f", df$logAMELY)
    df$`P(male)` <- sprintf("%.3f", df$`P(male)`)
    
    datatable(df[, c("Raw.file","Label","logAMELX","logAMELY",
                     "P(male)","Biological sex")],
              rownames = FALSE,
              editable = list(
        target = "cell",
        disable = list(columns = c(0, 2, 3, 4, 5))
      ),
      options = list(scrollX = TRUE)
    )
  })

  observeEvent(input$table1_cell_edit, {
    
    info <- input$table1_cell_edit
    
    df <- rv_table()
    
    req(df)
    
    if (info$col == 1) {
      df$Label[info$row] <- info$value
    }
    
    rv_table(df)
  })

  
  # ------------------- Plot -------------------
  
  output$plot1 <- renderPlot({
    
    df <- rv_table()
    
    req(df)
    
    ggplot(
      df,
      aes(
        logAMELX,
        logAMELY,
        colour = `Biological sex`
      )
    ) +
      geom_point(size = 4) +
      geom_text_repel(
        aes(
          label = Label,
          colour = `Biological sex`
        ),
        size = 4,
        max.overlaps = Inf,
        box.padding = 0.5,
        point.padding = 0.5,
        segment.size = 0.7,
        min.segment.length = 0,
        show.legend = FALSE
      ) +
      geom_smooth(
        data = subset(df, `Biological sex` == "Male"),
        aes(logAMELX, logAMELY),
        method = "lm",
        color = "black",
        fill = "gray",
        se = TRUE
      ) +
      scale_color_manual(
        values = c(
          "Female" = "#fec000",
          "Male" = "#5b9cd4",
          "Non-conclusive" = "#ff3b3b"
        )
      ) +
      theme_minimal(base_size = 16)
  })

  # ------------------- Downloads -------------------
  
  output$downloadGraph <- downloadHandler(
    filename = function() {
      paste0("graph-SexPeptID-", Sys.Date(), ".tiff")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), width = 12, height = 8, device = "tiff")
    }
  )
  
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste0("table-SexPeptID-", Sys.Date(), ".csv")
    },
    content = function(file) {
      
      df <- rv_table()
      req(df)
      
      df <- df[, c(
        "Raw.file",
        "Label",
        "logAMELX",
        "logAMELY",
        "P(male)",
        "Biological sex"
      )]
      
      write.csv(df, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)

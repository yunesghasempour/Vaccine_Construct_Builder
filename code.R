library(shiny)
library(stringr)
library(gtools)
library(grid)
library(gridExtra)
library(shinyjs)

# Function to validate sequences
validate_sequence <- function(sequence, name) {
  sequence <- str_trim(sequence)
  if (sequence == "") return(TRUE)
  
  invalid_chars <- str_extract_all(sequence, "[^ACDEFGHIKLMNPQRSTVWY,]")[[1]]
  if (length(invalid_chars) > 0) {
    stop(sprintf(
      "Error in %s: Invalid characters found: %s\nOnly valid amino acids (ACDEFGHIKLMNPQRSTVWY) are allowed.",
      name,
      paste(unique(invalid_chars), collapse = ", ")
    ))
  }
  return(TRUE)
}

# Function to generate vaccine constructs permutations
generate_vaccine_permutations <- function(ctl_epitopes, htl_epitopes, bcell_epitopes,
                                          ctl_count, htl_count, bcell_count,
                                          n_adjuvants, c_adjuvants, 
                                          linker_ctl, linker_htl,
                                          linker_bcell, linker_adjuvant) {
  
  # Convert empty strings to empty vectors
  ctl_epitopes <- if(length(ctl_epitopes) == 1 && ctl_epitopes == "") character(0) else ctl_epitopes
  htl_epitopes <- if(length(htl_epitopes) == 1 && htl_epitopes == "") character(0) else htl_epitopes
  bcell_epitopes <- if(length(bcell_epitopes) == 1 && bcell_epitopes == "") character(0) else bcell_epitopes
  n_adjuvants <- if(length(n_adjuvants) == 1 && n_adjuvants == "") character(0) else n_adjuvants
  c_adjuvants <- if(length(c_adjuvants) == 1 && c_adjuvants == "") character(0) else c_adjuvants
  
  # Check requested counts against available epitopes
  if (length(ctl_epitopes) > 0 && length(ctl_epitopes) < ctl_count) {
    stop(paste("Not enough CTL epitopes. Requested:", ctl_count, "Available:", length(ctl_epitopes)))
  }
  if (length(htl_epitopes) > 0 && length(htl_epitopes) < htl_count) {
    stop(paste("Not enough HTL epitopes. Requested:", htl_count, "Available:", length(htl_epitopes)))
  }
  if (length(bcell_epitopes) > 0 && length(bcell_epitopes) < bcell_count) {
    stop(paste("Not enough B-cell epitopes. Requested:", bcell_count, "Available:", length(bcell_epitopes)))
  }
  
  # Generate permutations for each epitope type if count > 0
  ctl_perms <- if (length(ctl_epitopes) > 0 && ctl_count > 0) {
    as.data.frame(permutations(n=length(ctl_epitopes), r=ctl_count, v=ctl_epitopes))
  } else {
    data.frame(dummy = character(0))
  }
  
  htl_perms <- if (length(htl_epitopes) > 0 && htl_count > 0) {
    as.data.frame(permutations(n=length(htl_epitopes), r=htl_count, v=htl_epitopes))
  } else {
    data.frame(dummy = character(0))
  }
  
  bcell_perms <- if (length(bcell_epitopes) > 0 && bcell_count > 0) {
    as.data.frame(permutations(n=length(bcell_epitopes), r=bcell_count, v=bcell_epitopes))
  } else {
    data.frame(dummy = character(0))
  }
  
  # Helper function to join segments with a linker
  combine_with_linker <- function(segments, linker) {
    segments <- segments[segments != ""]
    if (length(segments) == 0) return("")
    paste(segments, collapse = linker)
  }
  
  # Determine adjuvant combinations
  n_adj_combinations <- if(length(n_adjuvants) > 0) n_adjuvants else ""
  c_adj_combinations <- if(length(c_adjuvants) > 0) c_adjuvants else ""
  
  # Generate adjuvant combinations
  if (length(n_adjuvants) == 0 && length(c_adjuvants) == 0) {
    adj_combinations <- list(list(n = "", c = ""))
  } else {
    adj_combinations <- list()
    for (n_adj in c("", n_adj_combinations)) {
      for (c_adj in c("", c_adj_combinations)) {
        if (!(n_adj == "" && c_adj == "")) {
          adj_combinations[[length(adj_combinations) + 1]] <- list(n = n_adj, c = c_adj)
        }
      }
    }
  }
  
  constructs <- vector("character")
  
  # Ensure at least one row for each type if empty
  if (nrow(ctl_perms) == 0) ctl_perms <- data.frame(dummy = "")
  if (nrow(htl_perms) == 0) htl_perms <- data.frame(dummy = "")
  if (nrow(bcell_perms) == 0) bcell_perms <- data.frame(dummy = "")
  
  # Build all possible constructs
  for (adj_combo in adj_combinations) {
    n_adj <- adj_combo$n
    c_adj <- adj_combo$c
    
    for (i in 1:nrow(ctl_perms)) {
      ctl_seg <- if ("dummy" %in% names(ctl_perms)) "" else 
        combine_with_linker(as.character(ctl_perms[i,]), linker_ctl)
      
      for (j in 1:nrow(htl_perms)) {
        htl_seg <- if ("dummy" %in% names(htl_perms)) "" else 
          combine_with_linker(as.character(htl_perms[j,]), linker_htl)
        
        for (k in 1:nrow(bcell_perms)) {
          bcell_seg <- if ("dummy" %in% names(bcell_perms)) "" else 
            combine_with_linker(as.character(bcell_perms[k,]), linker_bcell)
          
          # Combine all segments
          all_segments <- c(ctl_seg, htl_seg, bcell_seg)
          all_segments <- all_segments[all_segments != ""]
          
          if (length(all_segments) > 0 || n_adj != "" || c_adj != "") {
            construct <- paste(all_segments, collapse = linker_ctl)
            
            # Add adjuvants
            if (n_adj != "") {
              construct <- paste(n_adj, if(construct != "") linker_adjuvant else "", construct, sep="")
            }
            if (c_adj != "") {
              construct <- paste(construct, if(construct != "") linker_adjuvant else "", c_adj, sep="")
            }
            
            constructs <- c(constructs, construct)
          }
        }
      }
    }
  }
  
  return(unique(constructs))
}

# Sample data for testing
sample_data <- list(
  ctl = "SIINFEKL,GILGFVFTL,NLVPMVATV",
  htl = "FVNWLKDGV,LQTTIHDII,PKYVKQNTLKLAT",
  bcell = "QYIKANSKFIGITE,SLLMWITQC,EYLQAFTY"
)

# Function to validate epitope count inputs
validate_epitope_counts <- function(epitopes, count, name) {
  if (count > length(epitopes) && length(epitopes) > 0) {
    stop(sprintf(
      "Selected %s count (%d) exceeds available epitopes (%d).",
      name,
      count,
      length(epitopes)
    ))
  }
}

# UI definition (same as before)
ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$style(HTML("\n      .error-text { color: red; }\n      .help-text { color: gray; font-style: italic; }\n      .sample-btn { margin: 5px; }\n    "))
  ),
  
  titlePanel("Vaccine Construct Builder"),
  
  sidebarLayout(
    sidebarPanel(
      div(style="margin-bottom: 20px",
          actionButton("load_sample", "Load Sample Data", class="sample-btn"),
          actionButton("clear_all", "Clear All", class="sample-btn")
      ),
      
      textInput("ctl_epitopes", "CTL Epitopes (comma-separated):", ""),
      helpText("Example: SIINFEKL,GILGFVFTL"),
      numericInput("ctl_count", "Number of CTL Epitopes:", value = 1, min = 0),
      textOutput("ctl_error"),
      
      textInput("htl_epitopes", "HTL Epitopes (comma-separated):", ""),
      helpText("Example: FVNWLKDGV,LQTTIHDII"),
      numericInput("htl_count", "Number of HTL Epitopes:", value = 1, min = 0),
      textOutput("htl_error"),
      
      textInput("bcell_epitopes", "B-Cell Epitopes (comma-separated):", ""),
      helpText("Example: QYIKANSKFIGITE,SLLMWITQC"),
      numericInput("bcell_count", "Number of B-Cell Epitopes:", value = 1, min = 0),
      textOutput("bcell_error"),
      
      textInput("n_adjuvants", "N-Terminal Adjuvants (comma-separated):", ""),
      textInput("c_adjuvants", "C-Terminal Adjuvants (comma-separated):", ""),
      
      textInput("linker_ctl", "CTL Linker:", "KK"),
      textInput("linker_htl", "HTL Linker:", "GPGPG"),
      textInput("linker_bcell", "B-Cell Linker:", "GPGPG"),
      textInput("linker_adjuvant", "Adjuvant Linker:", "GPGPG"),
      
      div(style="margin-top: 20px",
          actionButton("generate", "Generate Constructs", class="btn-primary"),
          downloadButton("download_csv", "Download CSV"),
          downloadButton("download_pdf", "Download PDF")
      )
    ),
    
    mainPanel(
      div(
        h4("Results"),
        textOutput("construct_count"),
        tableOutput("constructs_table")
      )
    )
  )
)

# Server logic (same as before)
server <- function(input, output, session) {
  constructs_reactive <- reactiveVal(NULL)
  
  observe({
    ctl_epitopes <- str_trim(unlist(strsplit(input$ctl_epitopes, ",")))
    htl_epitopes <- str_trim(unlist(strsplit(input$htl_epitopes, ",")))
    bcell_epitopes <- str_trim(unlist(strsplit(input$bcell_epitopes, ",")))
    
    if (length(ctl_epitopes) > 0 && input$ctl_count > length(ctl_epitopes)) {
      output$ctl_error <- renderText({
        "Error: Selected CTL epitopes exceed available count"
      })
    } else {
      output$ctl_error <- renderText({ NULL })
    }
    
    if (length(htl_epitopes) > 0 && input$htl_count > length(htl_epitopes)) {
      output$htl_error <- renderText({
        "Error: Selected HTL epitopes exceed available count"
      })
    } else {
      output$htl_error <- renderText({ NULL })
    }
    
    if (length(bcell_epitopes) > 0 && input$bcell_count > length(bcell_epitopes)) {
      output$bcell_error <- renderText({
        "Error: Selected B-Cell epitopes exceed available count"
      })
    } else {
      output$bcell_error <- renderText({ NULL })
    }
  })
  
  observeEvent(input$load_sample, {
    updateTextInput(session, "ctl_epitopes", value = sample_data$ctl)
    updateTextInput(session, "htl_epitopes", value = sample_data$htl)
    updateTextInput(session, "bcell_epitopes", value = sample_data$bcell)
  })
  
  observeEvent(input$clear_all, {
    updateTextInput(session, "ctl_epitopes", value = "")
    updateTextInput(session, "htl_epitopes", value = "")
    updateTextInput(session, "bcell_epitopes", value = "")
    updateTextInput(session, "n_adjuvants", value = "")
    updateTextInput(session, "c_adjuvants", value = "")
    updateTextInput(session, "linker_ctl", value = "KK")
    updateTextInput(session, "linker_htl", value = "GPGPG")
    updateTextInput(session, "linker_bcell", value = "GPGPG")
    updateTextInput(session, "linker_adjuvant", value = "GPGPG")
  })
  
  observeEvent(input$generate, {
    tryCatch({
      ctl_epitopes <- str_trim(unlist(strsplit(input$ctl_epitopes, ",")))
      htl_epitopes <- str_trim(unlist(strsplit(input$htl_epitopes, ",")))
      bcell_epitopes <- str_trim(unlist(strsplit(input$bcell_epitopes, ",")))
      n_adjuvants <- str_trim(unlist(strsplit(input$n_adjuvants, ",")))
      c_adjuvants <- str_trim(unlist(strsplit(input$c_adjuvants, ",")))
      
      validate_sequence(input$ctl_epitopes, "CTL Epitopes")
      validate_sequence(input$htl_epitopes, "HTL Epitopes")
      validate_sequence(input$bcell_epitopes, "B-Cell Epitopes")
      
      constructs <- generate_vaccine_permutations(
        ctl_epitopes, htl_epitopes, bcell_epitopes,
        input$ctl_count, input$htl_count, input$bcell_count,
        n_adjuvants, c_adjuvants,
        input$linker_ctl, input$linker_htl,
        input$linker_bcell, input$linker_adjuvant
      )
      
      constructs_reactive(constructs)
      output$construct_count <- renderText({
        paste("Number of constructs generated:", length(constructs))
      })
      output$constructs_table <- renderTable({
        data.frame(Constructs = constructs)
      })
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        e$message,
        easyClose = TRUE
      ))
    })
  })
  
  output$download_csv <- downloadHandler(
    filename = function() { "constructs.csv" },
    content = function(file) {
      write.csv(data.frame(Constructs = constructs_reactive()), file, row.names = FALSE)
    }
  )
  
  output$download_pdf <- downloadHandler(
    filename = function() { "constructs.pdf" },
    content = function(file) {
      pdf(file)
      grid.newpage()
      grid.table(data.frame(Constructs = constructs_reactive()))
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)

# To run the app:
# shiny::runApp("final code1.R")

library(shiny)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringi)
library(RColorBrewer)
library(parallel)
library(plotly)

options(shiny.maxRequestSize = 5 * 1024^3) 

# Reverse complement DNA sequence function
rc <- function (z) {
  rc1 <- function(zz) {
    s <- strsplit(zz, split = "")[[1]]
    s <- rev(s)
    dchars <- strsplit("ACGTMRWSYKVHDBNI", split = "")[[1]]
    comps <- strsplit("TGCAKYWSRMBDHVNI", split = "")[[1]]
    s <- s[s %in% dchars]
    s <- dchars[match(s, comps)]
    s <- paste0(s, collapse = "")
    return(s)
  }
  z <- toupper(z)
  tmpnames <- names(z)
  res <- unname(sapply(z, rc1))
  if (!is.null(attr(z, "quality"))) {
    strev <- function(x) sapply(lapply(lapply(unname(x), charToRaw), rev), rawToChar)
    attr(res, "quality") <- unname(sapply(attr(z, "quality"), strev))
  }
  names(res) <- tmpnames
  return(res)
}

# Function to find and plot matches
find_and_plot_matches <- function(df, query_df) {
  withProgress(message = 'Generating plot...', value = 0, {
    
    numberOfCores <- ceiling(detectCores() * 0.75)
    
    # Add row index as ID if not present
    if (!is.null(df) && !"id" %in% names(df)) {
      df$id <- seq_len(nrow(df))
    }
    
    # Calculate maximum oligo length
    appended_len <- nchar(df$appended_oligopaint)
    max_length <- max(appended_len, na.rm = TRUE)
    
    # Function to process each query
    process_query <- function(query_seq) {
      locs <- stringi::stri_locate_first_fixed(df$appended_oligopaint, query_seq) / max_length
      locs_temp <- as_tibble(locs) %>%
        mutate(id = row_number()) %>%
        filter(!is.na(start) & !is.na(end)) %>%
        arrange(start, end, id) %>%
        mutate(group_id = cumsum(c(TRUE, abs(diff(id)) != 1))) %>%
        group_by(start, end, group_id) %>%
        reframe(
          y_start = min(id),
          y_end = max(id)
        ) %>%
        mutate("query_seq" = query_seq) %>%
        select(start, end, y_start, y_end, query_seq)
      locs_temp
    }
    
    # Use mclapply to process queries in parallel
    out_coords_list <- mclapply(query_df$query_seq, process_query, mc.cores = numberOfCores)
    
    remove_jumps <- function(df) {
      query_seq_temp <- unique(df$query_seq)
      df %>%
        arrange(y_start) %>%
        mutate(distance = y_end - y_start) %>%
        filter(distance < 40) %>%
        summarize(start = min(start), end = max(end), y_start = min(y_start), y_end = max(y_end)) %>%
        mutate(query_seq = query_seq_temp)
    }
    
    for(i in seq_along(out_coords_list)) {
      if(nrow(out_coords_list[[i]]) > 1) {
        out_coords_list[[i]] <- remove_jumps(df = out_coords_list[[i]])
      } else {
        next
      }
    }
    
    incProgress(0.5, detail = "Combining results")
    # Combine results into a single tibble
    out_coords <- bind_rows(out_coords_list)
    out_coords <- left_join(out_coords, query_df, by = "query_seq")
    
    incProgress(0.5, detail = "Preparing plot")
    # Check for matches
    if (!nrow(out_coords)) {
      message("No matches found at all.")
      return(invisible(NULL))
    }
    
    # Determine colors for plotting
    seq_used <- unique(out_coords$query_seq)
    num_used <- length(seq_used)
    
    colors_list <- if (num_used > 8) {
      colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(num_used)
    } else {
      RColorBrewer::brewer.pal(max(3, num_used), "Dark2")[1:num_used]
    }
    
    # Assign colors to sequences
    color_mapping <- setNames(colors_list, seq_used)
    out_coords$color <- color_mapping[out_coords$query_seq]
    out_coords$start <- out_coords$start * max_length
    out_coords$end <- out_coords$end * max_length
    
    # Prepare label coordinates
    label_coords_seq <- out_coords %>%
      group_by(query_seq, query_name) %>%
      summarise(
        x = mean(start + end) / 2,
        y = max(y_end) - ((y_end - y_start) * 0.25),
        .groups = 'drop'
      )
    
    label_coords_id <- out_coords %>%
      group_by(query_seq, query_name) %>%
      summarise(
        x = mean(start + end) / 2,
        y = max(y_end) - ((y_end - y_start) * 0.10),
        .groups = 'drop'
      )
    
    # Create the plot
    p <- ggplot(data = out_coords, aes(x = (start + end) / 2, y = (y_start + y_end) / 2)) +
      geom_tile(aes(width = end - start, height = y_end - y_start, fill = query_seq), color = "white") +
      scale_fill_manual(values = color_mapping) +
      geom_text(aes(label = (y_end - y_start) + 1), size = 2, color = "black") +
      geom_text(data = label_coords_seq, aes(x = x, y = y, label = query_seq), size = 1, vjust = 0) +
      geom_text(data = label_coords_id, aes(x = x, y = y, label = query_name), size = 2, vjust = 0) +
      labs(x = "Oligopaint length", y = "Oligopaint number in the library") +
      theme_minimal() +
      scale_x_continuous(limits = c(0, max_length)) +
      theme(legend.position = "none")
    
    # Convert ggplot to plotly
    ggplotly(p)
  })
}

# UI
ui <- fluidPage(
  titlePanel("Oligopaints Insight"),
  fluidRow(
    column(4,
           fileInput("oligopaintFile", "Upload appended_oligopaint File", accept = c(".csv", ".tsv", ".txt"))
    ),
    column(4,
           fileInput("queryFile", "Upload Query File", accept = c(".csv", ".tsv", ".txt"))
    ),
    column(4,
           actionButton("plotButton", "Generate Plot", style = "margin-top: 25px;")
    )
  ),
  fluidRow(
    column(12,
           # Use plotlyOutput for an interactive plot
           plotlyOutput("sequencePlot", height = "600px"),
           # Add a loading spinner
           conditionalPanel(
             condition = "$('html').hasClass('shiny-busy')",
             tags$div("Loading...", id = "loadmessage", style = "position: fixed; top: 50%; left: 50%; margin-top: -50px; margin-left: -50px;")
           )
    )
  )
)

# Server
server <- function(input, output) {
  
  observeEvent(input$plotButton, {
    req(input$oligopaintFile)
    req(input$queryFile)
    
    # Progress bar setup
    withProgress(message = 'Processing...', value = 0, {
      
      # Read input files based on the extension
      read_file <- function(file) {
        ext <- tools::file_ext(file$name)
        switch(ext,
               csv = read_csv(file$datapath),
               tsv = read_tsv(file$datapath),
               txt = read_delim(file$datapath, delim = "\t"),
               stop("Invalid file type"))
      }
      
      incProgress(1/3, detail = "Reading library file")
      df <- read_file(input$oligopaintFile)
      
      incProgress(2/3, detail = "Reading query file")
      query_df <- read_file(input$queryFile)
      
      incProgress(3/3, detail = "Generating plot")
      output$sequencePlot <- renderPlotly({
        validate(
          need("appended_oligopaint" %in% colnames(df), "The df does not contain the column 'appended_oligopaint'"),
          need("query_seq" %in% colnames(query_df), "The query_df does not contain the column 'query_seq'"),
          need("query_name" %in% colnames(query_df), "The query_df does not contain the column 'query_name'")
        )
        
        find_and_plot_matches(df, query_df)
        
      })
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
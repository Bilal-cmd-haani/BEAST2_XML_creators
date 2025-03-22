

#C:/04072024_Firstyear_Analyses
#C:/Users/s2340217/OneDrive - University of Edinburgh/1st year/PhD work/20062024_Finalwork

# Clear the workspace
rm(list = ls())

# Load necessary libraries
library(xml2)
library(Biostrings)

# Specify paths
template_file <- "C:/multiplexmlfiles/S1a_19992014_sim_1.xml"  # Path to your original XML template file
#clean_template_file <- "C:/multiplexmlfiles/clean_example.xml"  # Path to save the cleaned XML template file
clean_template_file <- "C:/multiplexmlfiles/clean_example.xml"  # Path to the cleaned XML template file
fasta_folder <- "C:/multiplexmlfiles"  # Path to the folder containing FASTA files
output_folder <- "C:/multiplexmlfiles"  # Path to the output folder where new XML files will be saved

################################################################################

# Function to parse and clean the XML template
clean_xml_template <- function(xml_file, output_file) {
  cat("Parsing and cleaning XML template...\n")
  
  # Read the XML file
  xml_doc <- tryCatch({
    read_xml(xml_file)
  }, error = function(e) {
    stop("Error reading XML file: ", e$message)
  })
  
  # Remove existing sequences
  sequence_nodes <- xml_find_all(xml_doc, ".//data[@spec='Alignment']/sequence")
  if (length(sequence_nodes) > 0) {
    xml_remove(sequence_nodes)
    cat("Removed existing sequences.\n")
  } else {
    cat("No existing sequences to remove.\n")
  }
  
  # Clear date attributes
  date_nodes <- xml_find_all(xml_doc, ".//trait[@traitname='date']")
  if (length(date_nodes) > 0) {
    for (date_node in date_nodes) {
      xml_set_attr(date_node, "value", "")
    }
    cat("Cleared existing dates.\n")
  } else {
    cat("No date nodes found to clear.\n")
  }
  
  # Clear type (host) attributes
  type_nodes <- xml_find_all(xml_doc, ".//typeTrait[@traitname='type']")
  if (length(type_nodes) > 0) {
    for (type_node in type_nodes) {
      xml_set_attr(type_node, "value", "")
    }
    cat("Cleared existing type nodes (hosts).\n")
  } else {
    cat("No type nodes (hosts) found to clear.\n")
  }
  
  # Replace all occurrences of 'S1a_19992014_sim_1' with a placeholder
  xml_text <- as.character(xml_doc)
  if (grepl("S1a_19992014_sim_1", xml_text)) {
    xml_text <- gsub("S1a_19992014_sim_1", "PLACEHOLDER_TEMPLATE_NAME", xml_text)
    xml_doc <- read_xml(xml_text)
    cat("Replaced placeholder with 'PLACEHOLDER_TEMPLATE_NAME'.\n")
  } else {
    cat("No placeholder 'S1a_19992014_sim_1' found in XML.\n")
  }
  
  # Save the cleaned XML template
  tryCatch({
    write_xml(xml_doc, output_file)
    cat("Saved cleaned XML template to:", output_file, "\n")
  }, error = function(e) {
    stop("Error saving XML file: ", e$message)
  })
}

#################################################################################

# Function to parse the cleaned XML template
parse_xml_template <- function(xml_file) {
  cat("Parsing XML template...\n")
  tryCatch({
    xml_doc <- read_xml(xml_file)
    return(xml_doc)
  }, error = function(e) {
    stop("Error parsing XML template: ", e)
  })
}

# Function to read FASTA sequences and clean white space
read_fasta <- function(fasta_file) {
  cat("Reading FASTA file: ", fasta_file, "\n")
  fasta_content <- readLines(fasta_file)
  sequences <- list()
  current_sequence <- ""
  current_id <- NULL
  
  for (line in fasta_content) {
    if (startsWith(line, ">")) {
      if (!is.null(current_id)) {
        sequences[[current_id]] <- gsub("\\s", "", current_sequence) # Remove whitespace
      }
      current_id <- sub("^>", "", line)
      current_sequence <- ""
    } else {
      current_sequence <- paste0(current_sequence, line)
    }
  }
  
  if (!is.null(current_id)) {
    sequences[[current_id]] <- gsub("\\s", "", current_sequence) # Remove whitespace
  }
  
  cleaned_sequences <- DNAStringSet(unlist(sequences))
  names(cleaned_sequences) <- names(sequences)
  return(cleaned_sequences)
}

# Function to parse sequence information from the name
parse_sequence_info <- function(sequence_name) {
  parts <- strsplit(sequence_name, "_")[[1]]
  if (length(parts) < 3) {
    stop("Sequence name does not contain enough parts to extract ID, host, and date: ", sequence_name)
  }
  list(id = parts[1], host = parts[2], date = parts[3])
}

# Function to replace sequence, date, and type (host) information in the XML template
replace_info_in_xml <- function(xml_doc, cleaned_sequences, template_name) {
  # Find sections in the XML template
  sequence_section <- xml_find_first(xml_doc, ".//data[@spec='Alignment']")
  date_section <- xml_find_first(xml_doc, ".//trait[@traitname='date']")
  type_section <- xml_find_first(xml_doc, ".//typeTrait[@traitname='type']")
  
  if (is.null(sequence_section) || is.null(date_section) || is.null(type_section)) {
    stop("One or more required sections are missing in the XML template.")
  }
  
  # Prepare new data
  date_values <- c()
  type_values <- c()
  sequences <- list()
  
  for (i in seq_along(cleaned_sequences)) {
    full_seq_name <- names(cleaned_sequences)[i]
    seq_value <- as.character(cleaned_sequences[[i]])
    seq_info <- parse_sequence_info(full_seq_name)
    
    cat("Parsed sequence info: ", seq_info$id, ", ", seq_info$host, ", ", seq_info$date, "\n")
    
    seq_id <- seq_info$id
    seq_host <- seq_info$host
    seq_date <- seq_info$date
    
    # Collect date and type (host) values
    date_values <- c(date_values, sprintf("%s=%s", full_seq_name, seq_date))
    type_values <- c(type_values, sprintf("%s=%s", full_seq_name, seq_host))
    sequences <- c(sequences, list(list(id = full_seq_name, value = seq_value, taxon = full_seq_name)))
  }
  
  # Set date and type (host) attributes
  if (!is.null(date_section)) {
    xml_set_attr(date_section, "value", paste(date_values, collapse = ","))
    cat("Updated dates.\n")
  }
  
  if (!is.null(type_section)) {
    xml_set_attr(type_section, "value", paste(type_values, collapse = ","))
    cat("Updated types (hosts).\n")
  }
  
  # Add new sequences
  if (!is.null(sequence_section)) {
    for (seq in sequences) {
      cat("Adding sequence for ID: ", seq$id, "\n")
      sequence_node <- xml_add_child(sequence_section, "sequence")
      xml_set_attrs(sequence_node, c(id = seq$id, spec = "Sequence", taxon = seq$taxon, totalcount = "4", value = seq$value))
    }
    cat("Updated sequences.\n")
  }
  
  # Replace the placeholder with the template name
  xml_text <- as.character(xml_doc)
  xml_text <- gsub("PLACEHOLDER_TEMPLATE_NAME", template_name, xml_text)
  xml_doc <- read_xml(xml_text)
  
  return(xml_doc)
}

# Function to create a new XML file for each FASTA file
create_xml_for_fasta <- function(template_file, fasta_file, output_file) {
  cat("Processing FASTA file: ", fasta_file, "\n")
  cleaned_sequences <- read_fasta(fasta_file)
  
  # Parse the cleaned XML template afresh for each FASTA file
  xml_template <- parse_xml_template(template_file)
  
  template_name <- tools::file_path_sans_ext(basename(fasta_file))
  
  cat("Replacing sequence, date, and type (host) information in the XML template...\n")
  xml_doc <- replace_info_in_xml(xml_template, cleaned_sequences, template_name)
  
  cat("Saving new XML to file: ", output_file, "\n")
  write_xml(xml_doc, output_file)
  
  cat("Finished processing: ", fasta_file, "\n")
}

# Main function to process all FASTA files in a directory
generate_beast_xmls <- function(template_file, fasta_folder, output_folder) {
  cat("Generating BEAST XMLs...\n")
  
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
  }
  
  fasta_files <- list.files(fasta_folder, pattern = "\\.fasta$", full.names = TRUE)
  
  for (fasta_file in fasta_files) {
    output_file <- file.path(output_folder, paste0(tools::file_path_sans_ext(basename(fasta_file)), ".xml"))
    create_xml_for_fasta(template_file, fasta_file, output_file)
  }
  
  cat("All files processed.\n")
}

################################################################################
# Clean the XML template
clean_xml_template(template_file, clean_template_file)

# Run the main function
generate_beast_xmls(clean_template_file, fasta_folder, output_folder)

################################################################################
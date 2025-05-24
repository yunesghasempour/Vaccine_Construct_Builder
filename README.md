# Vaccine_Construct_Builder
A Shiny application for generating and downloading all possible vaccine construct sequences by combining CTL, HTL, B-cell epitopes and optional N- and C-terminal adjuvants with customizable linkers.
## Vaccine Construct Builder Shiny App

---

### Features

* **Sequence validation**: Ensures only valid amino acid characters are used (ACDEFGHIKLMNPQRSTVWY).
* **Permutation generation**: Creates all unique combinations (permutations) of specified epitope counts.
* **Adjuvant support**: Optionally prepend or append N-terminal and/or C-terminal adjuvants.
* **Custom linkers**: Define peptide linkers for CTL, HTL, B-cell segments and adjuvants.
* **Downloadable outputs**: Export results as CSV or PDF.
* **Sample data loader**: Quickly populate fields with example epitopes.
* **Error handling**: Clear on-screen messages when counts or sequences are invalid.

---

### Installation

1. Install R (>= 4.0)
2. Install required packages:

   ```r
   install.packages(c("shiny", "stringr", "gtools", "grid", "gridExtra", "shinyjs"))
   ```
3. Save the provided `final_code1.R` script to your working directory.
4. In R console, run:

   ```r
   shiny::runApp("final_code1.R")
   ```

---

### Usage

1. Start the app by running `final_code1.R`.
2. **Enter epitopes** (comma-separated) for:

   * CTL
   * HTL
   * B-cell
3. **Set counts** for each epitope type.
4. **(Optional)** Provide N-terminal and C-terminal adjuvant sequences.
5. **Customize linkers** between segments.
6. Click **Generate Constructs**.
7. View the number of constructs and table of sequences.
8. Download results via **Download CSV** or **Download PDF**.

---

### File Structure

This is a single-script Shiny app. All UI, server logic, and helper functions are contained within `final_code1.R`.

---

### License

This project is licensed under the MIT License.

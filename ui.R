# Copyright (C) 2023 Wenjie SUN <sunwjie@gmail.com>
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.

library(shiny)
library(plotly)

meta_data_var = c("lum_score", "bas_score", "PC_1", "PC_2", "PC_3", "PC_4", "PC_5", "time")

# Define UI for application that draws a histogram
shinyUI(
    # fluidPage(
    navbarPage(
        "Bonjour Embryonic Mammary gland",
        tabPanel("About", 
            HTML(r"(
                <p>This online visualization tool provides an interactive 3D trajectory and pseudotime plot to explore the Mammary Epithelial Single Cell expression data from E13.5 to P0 generated in the publication “Cell fate specification underlies positional cues during branching morphogenesis of the embryonic mammary epithelium” by Carabana et al (Institut Curie, Paris). Please check the paper for details regarding Materials and Methods.</p>
                <h2>3D trajectory plot</h2>
                <p>The 3D trajectory plot is derived from Figure 1G, which shows the differentiation path of embryonic MECs towards a luminal or basal fate. 
                </br>
                </br>
                The second principal component in PCA was attributed to the developmental stage (y-axis) and plotted against the basal and luminal scores used in this publication (x-axis). In the default view, colours represent each cell cluster. Single-cell gene expression levels can be visualized in the 3D trajectory by searching the official gene symbol in the “Search" box. The obtained plot is colour-coded by the normalized gene expression level.</p>
                <h2>Pseudotime gene expression plot</h2>
                <p>The pseudotime gene expression plot is derived from Figure 2. 
                </br>
                </br>
                The pseudotime is represented in the x-axis, with the luminal differentiation trajectory towards the right side and basal differentiation trajectory to the left side. Normalized gene expression levels are represented in the y-axis. Colours represent each cell cluster. You can search for a gene of interest by adding the official gene name in the “Search” box.</p>
                <h2>Note</h2>
                <p>No cookies have been used on this website. For questions, issues or comments, please contact: silvia.fre@curie.fr or wenjie.sun@curie.fr</p>
                )")
        ),
        ## choose gene for visualization
        tabPanel("3D Trajectory",
            sidebarLayout(
                sidebarPanel(
                    textInput("col", "Gene (or keep it empty to display cluster type)", value="", placeholder = "Input a gene symbol"),
                    submitButton("Submit")
                    ),
                mainPanel(
                    plotlyOutput("trajectory_3d", height="800px")
                )
            )
            ),
        tabPanel("Pseudotime",
            sidebarLayout(
                sidebarPanel(
                    textInput("gene", "Gene", value="Krt8", placeholder = "Input a gene symbol"),
                    submitButton("Submit")
                    ),
                mainPanel(
                    plotOutput("pseudotime")
                )
            )
        )
    )
)


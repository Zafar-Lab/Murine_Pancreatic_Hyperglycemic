library(readxl)
library(forestplot)
library(dplyr)
library(ggplot2)
file_name<- "Raw_associations_CCR2_2.xlsx"
sheet_names <- excel_sheets(file_name)
sheet_name <- sheet_names[15]
data <- read_excel("Raw_associations_CCR2_2.xlsx", sheet = sheet_name)

# Create a new variable to determine the fill
data$point_fill <- ifelse(data$pValue < 0.05, "filled", "empty")
tiff("rs113507038.tiff", units="in", width = 5, height = 5, res=300)
ggplot(data=data, aes(y=phenotype, x=beta,
                             xmin=lower, 
                             xmax=upper, shape = point_fill, fill = point_fill)) +
  geom_point(size = 4, color = "black", stroke = 1.5) + 
  geom_errorbarh(height=.1)+
  scale_shape_manual(values = c("filled" = 16, "empty" = 1)) +  # 16 is filled, 1 is empty
  scale_fill_manual(values = c("filled" = "blue", "empty" = NA)) +  # Blue for filled, NA for empty
  theme_minimal()+
    labs(title = sheet_name)+
  geom_vline(xintercept=0, color='green', linetype='dashed',alpha=0.8)

dev.off()



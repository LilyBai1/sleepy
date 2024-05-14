# most used functions:
# simple : runs mixed anova but also returns eta squared, which is standard when reporting stats
# get_top8_components : principal component analysis, returns top 8 principal components
# get_top8_cum_var : creates cum variance plot for pca
# calculate_and_save_correlation : calculates and saves correlation

install.packages("ggplot2")
library(ggplot2)

 if (!requireNamespace('BiocManager', quietly = TRUE))
        install.packages('BiocManager')

 BiocManager::install('PCAtools')
install.packages("PCAtools", force = TRUE)
library(PCAtools)


install.packages("car")
install.packages("psych")

library(car)
library(psych)

#*** Start Creating Matrix ***#
# builds the matrices and data frames used by subsequent stat functions
# set directory
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)

# read region of interest (ROI) values from csv
ROI_df <- read.csv(paste(datadir, '/PET_SEGMENTATION.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)
# make sure that in PET_SEGMENTATION.csv that the brain regions are rows, and the subjects are columns. first 8 columns are pre, the next 8 are post SD

# create matrices with pre sleep dep (pre) and sleep dep (SD) values for the ROI values
pre_ROI_matrix <- data.matrix(ROI_df[,2:9])
SD_ROI_matrix <- data.matrix(ROI_df[,10:17])

# indicate the ROI values using this vector
ROI_names <- ROI_df[,1]

# indicate the ROI on the matrices themselves
rownames(pre_ROI_matrix) <- ROI_names
rownames(SD_ROI_matrix) <- ROI_names

#*** End Creating Matrix ***#



#*** old PCA function, skip the functions plot_pca_result & get_eigenvalues, just run get_top8_components ***#

# generally based on: https://chat.openai.com/share/71e454aa-b5f5-4366-80e6-b9fd6a48fbdf
# do principal component analysis on 1 genomic matrix "data"
plot_pca_result <- function(data){
	#Perform PCA
  # Normalize the data: subtract mean and divide by standard deviation
  normalized_data <- scale(data)
	# the test subjects are the variables, which have to be in rows
	# the genes are the samples, which have to be in columns
	pca_result <- pca(t(normalized_data), scale. = TRUE)

	#Create a data frame for plotting
	pca_df <- data.frame(Sample = row.names(pca_result$rotated), 
                     PC1 = pca_result$rotated[,1], 
                     PC2 = pca_result$rotated[,2])

	#Plot the first two principal components
	ggplot(data = pca_df, aes(x = PC1, y = PC2)) +
  		geom_point() +
  		theme_classic() +
  		xlab("PC1") +
  		ylab("PC2") +
  		ggtitle("PCA Plot of Genomic Data")
}

get_eigenvalues <- function(data){
	#Perform PCA
	# the test subjects are the variables, which have to be in rows
	# the genes are the samples, which have to be in columns
	# removes bottom 10% of variation
	pca_result <- pca(t(data), scale = TRUE)
	
	# Get variances
	print('pca_result$variance')
	print(pca_result$variance)
	print('eigen(cov(data))')
	print(eigen(cov(data)))
	
}



# calculate the top 8 principal component analysis factors for matrix "data" 
get_top8_components <- function(data){

  # Normalize the data: subtract mean and divide by standard deviation
  normalized_data <- scale(data)
	#Perform PCA
	pca_result <- prcomp(t(normalized_data), scale. = TRUE)
	
	# Get the top 8 principal components
	# Does this do what it is supposed to?
	top8_pcs <- pca_result$x[,1:8]

	# View the top 8 principal components
	print(top8_pcs)
  write.csv(top8_pcs, "top8_ROI_PC.csv")
  
	# Get the loadings for the top 8 principal components
	loadings <- pca_result$rotation[,1:8]

	# View the loadings
	print(loadings)
  write.csv(loadings, "loadings.csv")
	
	# Get the eigenvalues from the PCA result
	eigenvalues <- pca_result$sdev^2
	# Create a scree plot
	plot(1:length(eigenvalues), eigenvalues, type = "b", xlab = "Principal Component", ylab = "Eigenvalue", main = "Scree Plot")
    
    # Calculate the proportion of variance explained by each principal component
	variance_proportion <- eigenvalues / sum(eigenvalues)
	
	# Create a scree plot with variance proportions
	plot(1:length(variance_proportion), variance_proportion, type = "b", xlab = "Principal Component", ylab = "Proportion of Variance", main = "Scree Plot: Variance Proportions")
}
#be sure you run "Start Creating Matrix" first to build the appropriate matrices
combined_ROI_matrix <- cbind(pre_ROI_matrix, SD_ROI_matrix)
get_top8_components(combined_ROI_matrix)



#11/27/2023 cumulative variance proportion plot
# you should first run the function: get_top8_components
get_top8_cum_var <- function(data){

  # Normalize the data: subtract mean and divide by standard deviation
  normalized_data <- scale(data)
	# Perform PCA
	pca_result <- prcomp(t(normalized_data), scale. = TRUE)
	
	# Get the eigenvalues from the PCA result
	eigenvalues <- pca_result$sdev^2
    
	# Calculate the proportion of variance explained by each principal component
	variance_proportion <- eigenvalues / sum(eigenvalues)
	
	# Calculate the cumulative variance
	cumulative_variance <- cumsum(variance_proportion)
	
	# Create a dataframe for plotting
	cumulative_variance_df <- data.frame(PC = 1:length(cumulative_variance), CumulativeVariance = cumulative_variance)
	
	# Create a ggplot2 scree plot with cumulative variance
	ggplot(cumulative_variance_df, aes(x = PC, y = CumulativeVariance)) +
		geom_line() +
		geom_point() +
		xlab("Principal Component") +
		ylab("Cumulative Proportion of Variance Explained") +
		ggtitle("Cumulative Variance Plot")
}

get_top8_cum_var(combined_ROI_matrix)



#***create difference matrix for ROI, followed by PCA on difference matrix*******#
ROI_difference_matrix <- pre_ROI_matrix - SD_ROI_matrix
colnames(ROI_difference_matrix) <- c("MoEm00F_1", "McDe00F_2", "AnAr00F_3", "MoKr00F_4", "CrBr00M_5", "CrCh00M_6", "MaJa00M_7", "WiAn00M_8")
print(ROI_difference_matrix)
diff_pca_result1 <- get_top8_components(ROI_difference_matrix)

#*** End PRINCIPAL COMPONENT ANALYSIS for each matrix ***#




#*** PCA-based unsupervised feature extraction (works as of 02/02/2024)***#

PCA_based_ufe <- function(data, k = 3){
  # Transpose the data frame
  t_data <- t(data)

  # Convert the transposed object to a data frame, ensuring strings are not converted to factors
  t_data_df <- as.data.frame(t_data, stringsAsFactors = FALSE)
  
  # Set the first row as column names
  colnames(t_data_df) <- as.character(unlist(t_data_df[1,]))
  
  # Remove the first row from the data frame
  t_data_df <- t_data_df[-1,]
  
  # Convert all columns to numeric
  t_data_df[] <- lapply(t_data_df, function(x) as.numeric(as.character(x)))
  
  # Proceed with normalization and PCA as before, now using t_data_df
  normalized_data <- scale(t_data_df)
  pca_result <- prcomp(normalized_data, scale. = TRUE)
  
  # Get the top 8 principal components
  top8_pcs <- pca_result$x[,1:8]
  
  # Get the loadings for the top 8 principal components
  loadings <- pca_result$rotation[,1:8]
  
  # Assign row names to loadings for plotting
  rownames(loadings) <- colnames(t_data_df)
  
  # Perform k-means clustering on the loadings for PC1 and PC2
  cluster_data <- loadings[,c(1, 2)]  # select PC1 and PC2 loadings
  clusters <- kmeans(cluster_data, centers = k)
  
  # Create a plot of loadings for PC1 and PC4 with cluster information
  plot(cluster_data, col = clusters$cluster, pch = 19, xlab = "Loadings for PC1", ylab = "Loadings for PC2", main = "Loadings Plot: PC1 vs PC2 with Clustering")
  text(cluster_data[,1], cluster_data[,2], labels = rownames(loadings), cex = 0.7, pos = 4)
  
  # Optional: Save plots and data
  print(top8_pcs)
  print(loadings)
}

datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)
ROI_df <- read.csv(paste(datadir, '/PET_SEGMENTATION.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)
# make sure that in PET_SEGMENTATION.csv that the brain regions are rows, and the subjects are columns. first 8 columns are pre, the next 8 are post SD
PCA_based_ufe(ROI_df, k = 3)



# Function that filters the regions of interest you want to look at from the , presumably from one of the PC loadings for PCA
# and then combines it with existing data frame that feeds into creating correltation matrix 02/01/2024

library(dplyr)
library(tidyr)

process_ROI_data <- function(data_directory, PC1_ROI_filename, ROI_fixed_filename, exhibit_D_filename, exhibit_F_filename) {
  # Set the working directory
  setwd(data_directory)
  
  # Read data
  PC_loading_ROI <- read.csv(file.path(data_directory, PC1_ROI_filename), header=TRUE, stringsAsFactors=FALSE)
  ROI_df <- read.csv(file.path(data_directory, ROI_fixed_filename), header=TRUE, stringsAsFactors=FALSE)
  
  # Filter ROI data
  filtered_ROI <- ROI_df %>%
    dplyr::filter(BrainRegion %in% PC_loading_ROI$brain_region)
  write.csv(filtered_ROI, "filtered_ROI.csv", row.names=FALSE)
  
  # Read the ROI data and replace "Pre" and "Post" with 1 and 2 respectively
  df_ROI <- read.csv('filtered_ROI.csv', header=TRUE, stringsAsFactors=FALSE) %>%
    mutate(Condition = ifelse(Condition == 'Pre', 1, 2))
  
  df_ROI_wide <- df_ROI %>%
    pivot_wider(
      names_from = BrainRegion, 
      values_from = MetVal
    ) %>%
    arrange(Patient_ID, Condition)
  
  # Read the exhibit_D data
  df_exhibit_D <- read.csv(file.path(data_directory, exhibit_D_filename), header=TRUE, stringsAsFactors=FALSE)
  
  # Combine the data
  combined_data <- merge(df_exhibit_D, df_ROI_wide, by.x = c("subject", "Condition"), by.y = c("Patient_ID", "Condition"), all.x = TRUE) %>%
    select(-Gender.y) %>% # Assuming there's a Gender.y column to exclude
    rename (Gender = Gender.x)

  # Output the combined data
  print(combined_data)
  write.csv(combined_data, file.path(data_directory, exhibit_F_filename), row.names=FALSE)
}

# Example usage, arguments: (data directory, list of things you want to find from list A, List A, data frame you want to add these regions to, output file name)
process_ROI_data("C:/Users/lilyb/Desktop/New folder", "Pc_123_ROI_loading_regions.csv", "ROI_fixed_fixed.csv", "exhibit_D.csv", "exhibit_I_3.csv")
#exhibit_D = original data with just demographics and PC1-8 of both genes and ROI
#exhibit_F = exhibit_D + PC1 ROI top 16 loadings
#exhibit_G = exhibit_F + PC2 ROI top 18 loadings
#exhibit_H = exhibit_G + PC3 ROI top 16 loadings
#exhibit_I = exhibit_H + PC4 ROI top loadings





#***************correlation matrix with P values*********************

#10_22_2023

install.packages("Hmisc")
library("Hmisc")

datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)
exhibit_D <- read.csv(paste(datadir, '/exhibit_D.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)

correlation_matrix_D_pval <- rcorr(as.matrix(exhibit_D))
print(correlation_matrix_D_pval)

correlation_matrix_D_pval_df <-do.call(data.frame, correlation_matrix_D_pval)

write.csv(correlation_matrix_D_pval_df, "correlation_matrix_D2_pval.csv", row.names=FALSE, quote=FALSE)



#************END CORRELATION MATRIX WITH P VALUES************************




#function to create correlation matrix -01/31/2024
library("Hmisc")
calculate_and_save_correlation <- function(df, output_filename) {
  # Calculate the correlation matrix with p-values
  correlation_matrix_pval <- rcorr(as.matrix(df))
  print(correlation_matrix_pval)
  
  # Convert the correlation matrix to a data frame
  correlation_matrix_pval_df <- do.call(data.frame, correlation_matrix_pval)

  # Round numeric columns
  numeric_columns <- sapply(correlation_matrix_pval_df, is.numeric)
  correlation_matrix_pval_df[numeric_columns] <- lapply(correlation_matrix_pval_df[numeric_columns], function(x) round(x, 5))

  
  # Write the data frame to a CSV file
  write.csv(correlation_matrix_pval_df, output_filename, row.names=FALSE, quote=FALSE)
}
# Usage example:
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)

# Read your data into a data frame
exhibit_D <- read.csv(paste(datadir, '/exhibit_I_3.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)

# Call the function with your data frame and desired output file name
calculate_and_save_correlation(exhibit_D, "correlation_matrix_I3_pval.csv")

# this is for corr matrix for roi pc1, pc2, and pc4 filtered regions
exhibit_H <- read.csv(paste(datadir, '/exhibit_F_2.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)
calculate_and_save_correlation(exhibit_H, "correlation_matrix_F2_pval.csv")
#finished


# making a difference matrix using roi_fixed_fixed, suitable for correlation matrix with pvt, sss, and gene pca ************ 11/16/2023
library(dplyr)

# set directory
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)

# read region of interest (ROI) values from csv
ROI_df <- read.csv(paste(datadir, '/ROI_fixed_fixed.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)

# Convert factors to characters to avoid issues in grouping
ROI_df$Patient_ID <- as.character(ROI_df$Patient_ID)
ROI_df$BrainRegion <- as.character(ROI_df$BrainRegion)

# Calculate the change in metabolism for each subject and brain region
metabolism_changes <- ROI_df %>%
  group_by(BrainRegion, Patient_ID) %>%
  summarise(Change_in_Metabolism = MetVal[Condition == "Post"] - MetVal[Condition == "Pre"]) %>%
  ungroup() # Ungroup for further manipulation if necessary
 # arrange(BrainRegion) #IT'S BROKEN... DO NOT UNCOMMENT THIS

# Create a new CSV file with the changes
write.csv(metabolism_changes, "metabolism_changes.csv", row.names = FALSE)

# add a new column with each subjects' change in PVT
library(tidyr)
# Use pivot_wider to spread BrainRegion into separate columns for each Patient_ID
df_wide <- metabolism_changes %>%
  pivot_wider(names_from = BrainRegion, values_from = Change_in_Metabolism)

# rename the first column to 'subject' IMPORTANT BECAUSE THE FIRST ROW NEEDS TO BE THE SAME IN ORDER TO MERGE IT WITH EXHIBIT D CHANGE
df_wide <- df_wide %>%
  rename(subject = names(df_wide)[1])

# The resulting data frame 'df_wide' will have each Patient_ID as a row and each BrainRegion as a column
write.csv(df_wide, "metabolism_changes_t.csv", row.names = FALSE)

# Read the CSV files into data frames
exhibit_D_for_merge <- read.csv("exhibit_D_change.csv")
df_wide_selected <- read.csv("metabolism_changes_t.csv")

# perform the join operation without the first column of df_wide
ROI_change_with_exhibitD_change_df <- left_join(exhibit_D_for_merge, df_wide_selected, by = "subject")

# remove columns 2, 3, and 5 because these two columns are gender, age, and condition, which do not change.
ROI_change_with_exhibitD_change_df <- ROI_change_with_exhibitD_change_df %>% select(-2, -3, -5)

write.csv(ROI_change_with_exhibitD_change_df, "Combined_exhibit_D_ROI_change.csv", row.names = FALSE)

# correlation matrix for this new df. includes p values for correlations
install.packages("Hmisc")
library("Hmisc")

datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)
exhibit_E <- read.csv(paste(datadir, '/Combined_exhibit_D_ROI_change.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)

correlation_matrix_E_pval <- rcorr(as.matrix(exhibit_E))
print(correlation_matrix_E_pval)

correlation_matrix_E_pval_df <-do.call(data.frame, correlation_matrix_E_pval)

write.csv(correlation_matrix_E_pval_df, "correlation_matrix_E_pval.csv", row.names=FALSE, quote=FALSE)

#**************finish difference correlation************************





#******************* difference correlation, but gender segregated********************** 11/27/2023
library(dplyr)

datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)

# read met diff data, transposed
df_wide_selected <- read.csv("metabolism_changes_t.csv", header = TRUE)

# Subset the data for females (subjects 1 to 4)
females <- df_wide_selected[1:4, ]

# Subset the data for males (subjects 5 to 8)
males <- df_wide_selected[5:8, ]

# If you want to save these to new CSV files
write.csv(females, 'females_metabolism_changes_t.csv', row.names = FALSE)
write.csv(males, 'males_metabolism_changes_t.csv', row.names = FALSE)

# separate exhibit D change by gender , exhibit_D_change was generated manually
exhibit_D_change <- read.csv(paste(datadir, '/exhibit_D_change.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)

exD_male_change <- exhibit_D_change[5:8, ]
exD_female_change <- exhibit_D_change[1:4, ]

write.csv(exD_female_change, 'exD_female_change.csv', row.names = FALSE)
write.csv(exD_male_change, 'exD_male_change.csv', row.names = FALSE)

exD_male_for_merge <- read.csv("exD_male_change.csv")
exD_female_for_merge <- read.csv("exD_female_change.csv")

# perform the join operation without the first column of df_wide
ROI_change_with_exhibitD_change_df_F <- left_join(exD_female_for_merge, females, by = "subject")
ROI_change_with_exhibitD_change_df_M <- left_join(exD_male_for_merge, males, by = "subject")

# remove columns 2, 3, and 5 because these two columns are gender, age, and condition, which do not change.
ROI_change_with_exhibitD_change_F_df <- ROI_change_with_exhibitD_change_df_F %>% select(-2, -3, -5)
ROI_change_with_exhibitD_change_M_df <- ROI_change_with_exhibitD_change_df_M %>% select(-2, -3, -5)
write.csv(ROI_change_with_exhibitD_change_F_df, "ROI_change_with_exhibitD_change_F_df.csv", row.names = FALSE)
write.csv(ROI_change_with_exhibitD_change_M_df, "ROI_change_with_exhibitD_change_M_df.csv", row.names = FALSE)

# correlation matrix for this new df. includes p values for correlations
install.packages("Hmisc")
library("Hmisc")

gender_change_corr <- function(data){
  correlation_diff_matrix_Gender_pval <- rcorr(as.matrix(data))
  print(correlation_diff_matrix_Gender_pval)

  correlation_diff_matrix_Gender_pval_df <-do.call(data.frame, correlation_diff_matrix_Gender_pval)

  # Create a new filename by appending "corr" to the data_name
  output_filename <- paste0(data_name, "_corr.csv")

  # Write the data frame to a CSV file with the new filename
  write.csv(correlation_diff_matrix_Gender_pval_df, output_filename, row.names = FALSE, quote = FALSE)

}

gender_change_corr(ROI_change_with_exhibitD_change_F_df)
gender_change_corr(ROI_change_with_exhibitD_change_M_df)
#END





#******************start gender segregated correlation matrix*******************
# set directory
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)

# read region of interest (ROI) values from csv
#ROI_df <- read.csv(paste(datadir, '/PET_SEGMENTATION.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)
# make sure that in PET_SEGMENTATION.csv that the brain regions are rows, and the subjects are columns. first 8 columns are pre, the next 8 are post SD

# create matrices for both genders
#female_subjects <- ROI_df[, c(1, 2:5, 10:13)]
#male_subjects <- ROI_df[, -c(2:5, 10:13)]

#write.csv(female_subjects, 'female_subjects.csv', row.names = FALSE)
#write.csv(male_subjects, 'male_subjects.csv', row.names = FALSE)

# ahaha just realized that I should use the roi_fixed_fixed matrix instead of the raw pet segmentation matrix

ROI_df2 <- read.csv(paste(datadir, '/ROI_fixed_fixed.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)
# Create a logical vector indicating rows for female subjects
is_female <- ROI_df2$Gender == 'F'

# Subset the data frame for female subjects
female_subjects2 <- ROI_df2[is_female, ]

# Subset the data frame for male subjects
male_subjects2 <- ROI_df2[!is_female, ]

# write to CSV
write.csv(female_subjects2, 'female_subjects2.csv', row.names = FALSE)
write.csv(male_subjects2, 'male_subjects2.csv', row.names = FALSE)


library(tidyr)
library(dplyr)
# Use pivot_wider to spread BrainRegion into separate columns for each Patient_ID, and then rename patient_ID to subject
df_wide_female <- female_subjects2 %>%
  pivot_wider(names_from = BrainRegion, values_from = MetVal)
df_wide_female <- df_wide_female %>%
  rename(subject = names(df_wide_female)[1])

df_wide_male <- male_subjects2 %>%
  pivot_wider(names_from = BrainRegion, values_from = MetVal)
  df_wide_male <- df_wide_male %>%
  rename(subject = names(df_wide_male)[1])


# write to CSV
write.csv(df_wide_female, 'female_subjects2_t.csv', row.names = FALSE)
write.csv(df_wide_male, 'male_subjects2_t.csv', row.names = FALSE)

# separate exhibit D by gender
exhibit_D <- read.csv(paste(datadir, '/exhibit_D.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)

# Separate the subjects by gender
exD_male_subjects2 <- exhibit_D[exhibit_D$Gender == 1, ]
exD_female_subjects2 <- exhibit_D[exhibit_D$Gender == 2, ]

# write to CSV
write.csv(exD_female_subjects2, 'exD_female.csv', row.names = FALSE)
write.csv(exD_male_subjects2, 'exD_male.csv', row.names = FALSE)

# Read the CSV files into data frames
exD_male_for_merge <- read.csv("exD_male.csv")
exD_female_for_merge <- read.csv("exD_female.csv")

ROI_male_for_merge <- read.csv("male_subjects2_t.csv")
ROI_female_for_merge <- read.csv("female_subjects2_t.csv")
# Select the columns to keep from the ROI data frames (excluding columns 2 and 3)
ROI_male_for_merge_selected <- select(ROI_male_for_merge, -c(2, 3))
ROI_female_for_merge_selected <- select(ROI_female_for_merge, -c(2, 3))

# perform the join operation, make sure subjects column match names
ROI_exhibitD_F_df <- left_join(exD_female_for_merge, ROI_female_for_merge_selected, by = "subject")
ROI_exhibitD_M_df <- left_join(exD_male_for_merge, ROI_male_for_merge_selected, by = "subject")
# write to CSV
write.csv(ROI_exhibitD_F_df, 'ROI_exhibitD_F_df.csv', row.names = FALSE)
write.csv(ROI_exhibitD_M_df, 'ROI_exhibitD_M_df.csv', row.names = FALSE)

# correlation matrix with P values for this new df. includes p values for correlations
install.packages("Hmisc")
library("Hmisc")

datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)
exhibit_D_F <- read.csv(paste(datadir, '/ROI_exhibitD_F_df.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)
exhibit_D_M <- read.csv(paste(datadir, '/ROI_exhibitD_M_df.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)

correlation_matrix_D_F_pval <- rcorr(as.matrix(exhibit_D_F))
correlation_matrix_D_M_pval <- rcorr(as.matrix(exhibit_D_M))

corr_matrix_D_F_pval_df <-do.call(data.frame, correlation_matrix_D_F_pval)
corr_matrix_D_M_pval_df <-do.call(data.frame, correlation_matrix_D_M_pval)

write.csv(corr_matrix_D_F_pval_df, "corr_matrix_D_F_pval.csv", row.names=FALSE, quote=FALSE)
write.csv(corr_matrix_D_M_pval_df, "corr_matrix_D_M_pval.csv", row.names=FALSE, quote=FALSE)

#END (gender segregated correlation matrices)


#***********************sorting significant results from correlation matrix*******************************
# no need to load any additional libraries if working in base R
# Function to list significant correlations
list_significant_correlations <- function(data) {
  # Initialize an empty vector to store results
  significant_correlations <- character()
  
  # Loop through each row of the data frame
  for (i in 1:nrow(data)) {
    # Extract the brain area or variable name
    brain_area <- as.character(data[i, 1])
    # Replace periods with spaces in the brain area name
    brain_area_cleaned <- gsub("\\.", " ", brain_area)
    
    # Loop through columns, skipping the first one (brain area names)
    for (j in 2:ncol(data)) {
      column_name <- colnames(data)[j]
      
      # Check if this is an r value column and has a corresponding p-value column
      if (grepl("^r\\.", column_name) && (paste0('P', substr(column_name, 2, nchar(column_name))) %in% colnames(data))) {
        r_value <- data[i, column_name]
        p_column_name <- paste0('P', substr(column_name, 2, nchar(column_name)))
        p_value <- data[i, p_column_name]
        
        # Check if the correlation is significant (p < 0.05)
        if (!is.na(p_value) && p_value < 0.05) {
          # Create a string describing the significant correlation
          correlation_str <- sprintf("significant correlations were found between %s and %s (r(16) = %.2f, p = %.2f)",
                                      brain_area_cleaned, gsub("^r\\.", "", column_name), r_value, p_value)
          # Append the string to the results vector
          significant_correlations <- c(significant_correlations, correlation_str)
        }
      }
    }
  }
  
  # Return the significant correlations
  if (length(significant_correlations) == 0) {
    return("No significant correlations were found.")
  } else {
    return(significant_correlations)
  }
}
  
#you can probably just plug in any of the correlation matrices generated to run this, but I recommend abridging the corr mat to a smaller size
#basically lists out the significant correlations in the format r(16) = X, p = Y, useful for writing manuscript.
#usage example:
  #read in the data frame
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)
sig_df <- read.csv(paste(datadir, '/correlation_matrix_I4_pval.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)
  #call function
significant_results <- list_significant_correlations(sig_df)
print(significant_results)
#  END sorting out and printing significant correlations



# function to transform data from long to wide
library(tidyverse)
wide_to_long <- function(df, cols_to_long, names_to = "name", values_to = "value") {
  long_df <- df %>% 
    pivot_longer(cols = starts_with("PC"),
                 names_to = "Principal_Component",
                 values_to = "Value")
  return(long_df)
}
#read in the data frame
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)
genepca <- read.csv(paste(datadir, '/genePCA.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)
  #call function
dataframetemp <- wide_to_long(genepca)
write.csv(dataframetemp, "genePCA_long.csv", row.names=FALSE, quote=TRUE)

#finished converting pca data to long



#*** start hierarchical clustering***#
hier_clustering <- function(data) {
	
	#compute dissimilarity matrix
	dis_matrix <- dist(data, method = "euclidean")
	hierarchical_results <- hclust(dis_matrix, method = "complete") #method = "ward.D"
	plot(hierarchical_results, main = "Dendogram of Hierarchical Clustering")

	#cutting dendogram
	clusters <- cutree(hierarchical_results, k = 3) # k = 3 cluster, can change this to whatever number of clusters I need
}

hier_clustering(ROI_difference_matrix) #FOR EXAMPLE

#*** END hierarchical clustering***#





#***************START MIXED REPEATED MEASURES ANOVA WITH STOREYQ FDR CORRECTION*****************************

#with P-values  mixed ANOVA with ROI brain region values
#ok this one actually works, ignore other versions. Not sure why this works. good luck 09/28/2023
library(lme4)
library(car)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")

BiocManager::install("qvalue")
library(qvalue)

if (!require("lmerTest")) {
    install.packages("lmerTest")
}
library(lmerTest)

# set directory
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)

# Read in the CSV data
df <- read.csv('ROI_fixed_fixed.csv', header = TRUE, stringsAsFactors = FALSE) 

get_lmer_value <- function(brain_region){
  
  # Subset data for the current brain region
  sub_data <- df[df$BrainRegion == brain_region, ]

  # Set up model formula
  model_formula <- MetVal ~ Gender * Condition + (1|Patient_ID)

  # Needed for Type III
  type3 <- list(Patient_ID = contr.sum, Gender = contr.sum, Condition = contr.sum)
  
  # Fit the mixed model
  mixed_model <- lmer(model_formula, data = sub_data, contrasts=type3)
  mixed_anova <- Anova(mixed_model, test="F", type=3)

  # Get f and Pr(>F), 1=F value, 4=p value of R. unlist creates vector instead of list, c indicates vector
  anova_values <- unlist(mixed_anova[2:4,c(1,4)])
  
  # Compute some means (this part is based on your original code, but may need adjustment)
  met_values <- sub_data$MetVal
  important_data <- data.frame(c(anova_values, mean(met_values[1:8]), mean(met_values[9:16]), mean(met_values[1:8 * 2]), mean(met_values[1:8 * 2 - 1]), mean(met_values[1:4 * 2 ]), mean(met_values[1:4 * 2 - 1]), mean(met_values[5:8 * 2 ]), mean(met_values[5:8 * 2 - 1])))
  
  # Return as a row
  return(t(important_data))
}

# Apply the function across unique brain regions
brain_regions <- unique(df$BrainRegion)
lmer_list <- lapply(brain_regions, get_lmer_value)
lmer_df <- do.call(rbind, lmer_list)

# Assign column names
colnames(lmer_df) <- c('Gender F', 'Condition F', 'Gender:Condition F', 'Gender P', 'Condition P', 'Gender:Condition P', 'Mean Female', 'Mean Male', 'Mean Pre-SD', 'Mean Post-SD', 'Mean Female Pre-SD', 'Mean Female Post-SD', 'Mean Male Pre-SD', 'Mean Male Post-SD')

#save results
write.csv(lmer_df, "unga_bunga_get_lmer_value_ROI_fixed_fixed.csv", row.names=FALSE, quote=FALSE) #has issues with row names



# set directory
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)

# Read in the CSV data
df <- read.csv('ROI_fixed_fixed.csv', header = TRUE, stringsAsFactors = FALSE) 

# Make sure to install these packages if you haven't already
install.packages("afex")
install.packages("lsr")

library(afex)
library(lsr)

conduct_mixed_anova <- function(data) {
  # Assuming your data frame is structured with the following columns:
  # Patient_ID, Gender, Condition, BrainRegion, MetVal

  # Convert factors to factor type if they are not already
  data$Patient_ID <- as.factor(data$Patient_ID)
  data$Gender <- as.factor(data$Gender)
  data$Condition <- as.factor(data$Condition)
  data$BrainRegion <- as.factor(data$BrainRegion)
  
  # Conduct the mixed ANOVA
  result <- aov_ez("Patient_ID", "MetVal", data, between = c("Gender", "BrainRegion"), within = c("Condition"), type = 3)
  
  # Calculate partial eta squared
  eta_sq <- etaSquared(result)
  
  return(list(ANOVA = result, Partial_Eta_Squared = eta_sq))
}

your_data_frame <- read.csv("ROI_fixed_fixed.csv")
results <- conduct_mixed_anova(your_data_frame)
print(results$ANOVA)
print(results$Partial_Eta_Squared)

# Setting column names for the 9 ANOVA result columns
colnames(lmer_df) <- c(
  'F_Gender', 'p_Gender', 'EtaSq_Gender',
  'F_Condition', 'p_Condition', 'EtaSq_Condition',
  'F_GenderCondition', 'p_GenderCondition', 'EtaSq_GenderCondition'
)

# Add brain region names
brainregion_lmer_df <- data.frame('BrainRegion' = brain_regions, lmer_df)

#rename rows to brain region
brainregion_lmer_df_renamed <- data.frame(brainregion_lmer_df)
rownames(brainregion_lmer_df_renamed) <- brainregion_lmer_df$BrainRegion

#save results
write.csv(brainregion_lmer_df_renamed, "unga_bunga_fixed_cond_.csv", row.names=FALSE, quote=FALSE)

#mental health deteriorates


# Storey Q Correction for p-values
# Read in the CSV data
df_2_for_FDR <- read.csv('unga_bunga_fixed_cond_part_squares.csv', header = TRUE, stringsAsFactors = FALSE) 

#qvalue_truncp function was removed from the qvalue package pepesmoge, so I define it here before calling it later
#@export
qvalue_truncp <- function(p, fdr.level = NULL, pfdr = FALSE, lfdr.out = TRUE, pi0 = NULL, ...) {
  # Argument checks
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  } else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 1)) {
    stop("'fdr.level' must be in (0, 1].")
  }
  p <- p / max(p)
  # Calculate pi0 estimate
  if (is.null(pi0)) {
    pi0s <- pi0est(p, ...)
  } else {
    if (pi0 > 0 && pi0 <= 1)  {
      pi0s = list()
      pi0s$pi0 = pi0
    } else {
      stop("pi0 is not (0,1]")
    }
  }
  
  # Calculate q-value estimates
  m <- length(p)
  i <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  if (pfdr) {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m / (i * (1 - (1 - p[o]) ^ m))))[ro]
  } else {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m /i ))[ro]
  }
  qvals_out[rm_na] <- qvals
  # Calculate local FDR estimates
  if (lfdr.out) {
    lfdr <- lfdr(p = p, pi0 = pi0s$pi0, ...)
    lfdr_out[rm_na] <- lfdr
  } else {
    lfdr_out <- NULL
  }
  
  # Return results
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, fdr.level = fdr.level,
                   significant = (qvals <= fdr.level),
                   pi0.lambda = pi0s$pi0.lambda, lambda = pi0s$lambda,
                   pi0.smooth = pi0s$pi0.smooth)
  } else {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, pi0.lambda = pi0s$pi0.lambda,
                   lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth)
  }
  class(retval) <- "qvalue"
  return(retval)
}

# Define the function to compute Storey Q correction the previously defined function
qvalue_corrected <- function(p_vals) {
  q_obj <- qvalue_truncp(p_vals)
  return(q_obj$qvalues)
}

# Apply Storey Q correction on p-values
df_2_for_FDR$Gender_Q <- qvalue_corrected(df_2_for_FDR$`Gender.P`)
df_2_for_FDR$Condition_Q <- qvalue_corrected(df_2_for_FDR$`Condition.P`)
df_2_for_FDR$GenderCondition_Q <- qvalue_corrected(df_2_for_FDR$`Gender.Condition.P`)

# Write the modified dataframe to a new CSV file
output_path <- "unga_bunga_storey_Q_corrected_fixed_cond_part_sq.csv" # You can change this to the desired output path
write.csv(df_2_for_FDR, output_path, row.names = FALSE)

cat("Q-value correction applied and results saved to:", output_path)





#Ok so the code below runs mixed anova again, but also gives partial eta squared, which is a stat that is commonly presented in scientific papers in anova tables. 12/06/2023
#the function simple should do the mixed anova for brain regions and ROI, the simple1 function is for genes. 

library(lme4)
library(car)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")

BiocManager::install("qvalue")
library(qvalue)

if (!require("lmerTest")) {
    install.packages("lmerTest")
}
library(lmerTest)

simple <- function(brain_region){
  
  # Subset data for the current brain region
  sub_data <- df[df$BrainRegion == brain_region, ]

  if(nrow(sub_data) == 0){
    return(NA)  # Returns NA if no data for the brain region
  }

  # Set up model formula
  model_formula <- MetVal ~ Gender * Condition + (1|Patient_ID)

  # Define contrasts for type III
  type3_contrasts <- list(Patient_ID = contr.sum, Gender = contr.sum, Condition = contr.sum)
  
  # Fit the mixed model with explicit contrasts
  mixed_model <- lmer(model_formula, data = sub_data, contrasts = type3_contrasts)
  mixed_anova <- Anova(mixed_model, test="F", type=3)

  # Convert the ANOVA output to a dataframe
  anova_values <- data.frame(mixed_anova)

  # Check column names
  print(colnames(anova_values)) 

  residual_df_colname <- 'Df.res' # Replace with the correct column name if needed

  # Calculate eta squared
  anova_values$eta_sq <- with(anova_values, ifelse(F > 0, (Df * F) / (Df * F + get(residual_df_colname)), 0))

  # Replace NA with 0 in eta squared
  anova_values$eta_sq[is.na(anova_values$eta_sq)] <- 0

  # Add brain region to the data frame
  anova_values$BrainRegion <- brain_region
  
  # Return as a row
  return(t(anova_values))
}

# set directory
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)

# Read in the CSV data
df <- read.csv('ROI_fixed_fixed.csv', header = TRUE, stringsAsFactors = FALSE) 
# Apply the function across unique brain regions
brain_regions <- unique(df$BrainRegion)
lmer_list <- lapply(brain_regions, simple)
lmer_df <- do.call(rbind, lmer_list)

#save results
write.csv(lmer_df, "mixed_anova_eta_squares.csv", row.names = TRUE)


#*******************************END MIXED REPEATED MEASURES ANOVA WITH STOREYQ FDR CORRECTION********************************


# The things below this are mostly experimental functions. the accuracy of the tests are questionable for the sleep dep data set but
# feel free to use on whatever you may need it on. 



 
#**************Simple 2 group T-test with affyRNAdeg slopes for condition****************** 10/05/2023

Group1 <- c(1.353166, 1.794559, 1.623671, 2.080824, 2.25435, 1.691016, 1.828324, 1.019396)
Group2 <- c(2.732662, 2.194767, 2.027774, 2.43369, 2.582972, 2.263606, 1.323386, 1.843323)

Group1_mean <- mean(Group1)
Group2_mean <- mean(Group2)
T.test_result <- t.test(Group1, Group2, paired = TRUE)

Group1_sd <- sd(Group1)
Group2_sd <- sd(Group2)

cat("Mean Pre Sleep Deprivation AffydegRNA slopes: ", Group1_mean, "\n")
cat("Mean Post Sleep Deprivation AffydegRNA slopes: ", Group2_mean, "\n")
cat("SD Pre Sleep Deprivation AffydegRNA slopes: ", Group1_sd, "\n")
cat("SD Post Sleep Deprivation AffydegRNA slopes: ", Group2_sd, "\n")

print(T.test_result)

#**************Simple 2 group T-test with affyRNAdeg slopes for gender****************** 10/15/2023

female <- c(1.353166, 1.794559, 1.623671, 2.080824, 2.732662, 2.194767, 2.027774, 2.43369)
male <- c(2.25435, 1.691016, 1.828324, 1.0193962, 2.582972, 2.263606, 1.323386, 1.843323)

female_mean <- mean(female)
male_mean <- mean(male)
T.test_result_gender <- t.test(female, male, paired = FALSE)

female_sd <- sd(female)
male_sd <- sd(male)

cat("Mean Pre Sleep Deprivation AffydegRNA slopes: ", female_mean, "\n")
cat("Mean Post Sleep Deprivation AffydegRNA slopes: ", male_mean, "\n")
cat("SD Pre Sleep Deprivation AffydegRNA slopes: ", female_sd, "\n")
cat("SD Post Sleep Deprivation AffydegRNA slopes: ", male_sd, "\n")

print(T.test_result_gender)
# end





#******************grouping subjects based off their sensitivity to sleep deprivation, perhaps k-means clustering

library(tidyverse)

#create df from PVT data
sensitivity_df <- data.frame(
  subject = 1:8,
  change_in_PVT = c(85.1100134, 107.2399866, 281.3100049, 170.5200061, 50.8999939, 15.2599946, 26.3700104, 17.4700164)
)

num_clusters <- 2
set.seed(13)


clusters <- kmeans(sensitivity_df$change_in_PVT, centers = num_clusters)

print(clusters)

#add cluster assignments back to original data frame
sensitivity_df$Cluster <- as.factor(clusters$cluster)

print(sensitivity_df)

#charactaristics of cluster
cluster_summary <- sensitivity_df %>%
group_by(Cluster) %>%
summarise(
  Count = n(),
  Avg_change_in_PVT = mean(change_in_PVT)
)
print(cluster_summary)

#plot clusters
library(ggplot2)
ggplot(sensitivity_df, aes(x = change_in_PVT, y = subject))
# end



#elbow method for determining number of clusters

set.seed(123)
# Check the number of subjects you have
num_subjects <- nrow(sensitivity_df)
wss <- sapply(1:(num_subjects - 1), function(k) {
  kmeans(sensitivity_df$change_in_PVT, centers = k, nstart = 50)$tot.withinss
})

# Plot the results
plot(1:(num_subjects - 1), wss, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K", 
     ylab = "Total within-clusters sum of squares")





#**************************clustering with principal components****************************
#create df from PCA ROI data, specifically pc6
sensitivity_PC6_roi_df <- data.frame(
  subject = 1:8,
  change_in_PC6_roi = c(4.365581187, 0.58652279, 0.4042733, 2.894853216, 1.646589819, 3.606740265, 0.190583845, 4.598098631)
)

num_clusters <- 2
set.seed(13)


clusters <- kmeans(sensitivity_PC6_roi_df$change_in_PC6_roi, centers = num_clusters)

print(clusters)

#add cluster assignments back to original data frame
sensitivity_PC6_roi_df$Cluster <- as.factor(clusters$cluster)

print(sensitivity_PC6_roi_df)

#charactaristics of cluster
cluster_summary <- sensitivity_PC6_roi_df %>%
group_by(Cluster) %>%
summarise(
  Count = n(),
  Avg_change_in_PC6_roi = mean(change_in_PC6_roi)
)
print(cluster_summary)

#plot clusters
library(ggplot2)
ggplot(sensitivity_PC6_roi_df, aes(x = change_in_PC6_roi, y = subject))

#elbow method for determining number of clusters

set.seed(123)
# Check the number of subjects you have
num_subjects <- nrow(sensitivity_PC6_roi_df)
wss <- sapply(1:(num_subjects - 1), function(k) {
  kmeans(sensitivity_PC6_roi_df$change_in_PC6_roi, centers = k, nstart = 50)$tot.withinss
})

# Plot the results
plot(1:(num_subjects - 1), wss, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K", 
     ylab = "Total within-clusters sum of squares")


#***********************************AUTOMATION OF KiMEANS CLUSTERING  11/15/2023 ******************************************


#create difference matrix for all the principal components using exhibit D
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)
exhibit_D <- read.csv(paste(datadir, '/exhibit_D.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)

# create data frames with pre sleep dep (pre) and sleep dep (SD) values for the ROI values
pre_PC_matrix <- exhibit_D[1:8, ]
SD_PC_matrix <- exhibit_D[9:16, ]

# Calculate the change for each variable
# Assuming that the first column is 'subject' and should be excluded from the subtraction
changes <- SD_PC_matrix[, -1] - pre_PC_matrix[, -1]

# Create a new data frame to store the changes
# The subject column is retained from the pre sleep deprivation data frame
changes_df <- data.frame(subject = pre_PC_matrix[, 'subject'], changes)

# View the changes data frame
print(changes_df)

# Write the modified dataframe to a new CSV file
output_path <- "exhibit_D_change.csv" # You can change this to the desired output path
write.csv(changes_df, output_path, row.names = FALSE)

cat("Change in exhibit D factors by condition saved to:", output_path)


# clustering for each column seperately
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)
exhibit_D_change1 <- read.csv(paste(datadir, '/exhibit_D_change1.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)
k <- 2
# Loop through each column you want to cluster
auto_k_means_clust <- function(data){
for (i in names(data)) {
  # Assuming the columns you want to cluster are numeric
  if (is.numeric(data[[i]])) {
    print(paste("Clustering on", i))
    clusters <- kmeans(data[[i]], centers = k)
    print(clusters)
  }
}
}

auto_k_means_clust(exhibit_D_change1)

# Load necessary library
library(readr)

datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)
exhibit_D_change <- read.csv(paste(datadir, '/exhibit_D_change.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)

# Remove non-numeric columns if necessary (e.g., subject identifier)
# Assuming the first column 'subject' is an identifier and should be excluded
#data_for_clustering <- exhibit_D_change[, -c(1:6)]
data_for_clustering <- exhibit_D_change[,7:22]

# Determine the number of clusters (k)
# You might need to use methods like the elbow method to determine a good number for k
k <- 2

# Perform k-means clustering
# Note: It's important to scale your data because k-means is distance-based and features on larger scales can unduly influence the algorithm
data_scaled <- scale(data_for_clustering)  # Scaling the data
clusters <- kmeans(data_scaled, centers = k)

# The result is a list containing cluster assignments, cluster centers, and more.
print(clusters)

# Add the cluster assignments to your data frame
exhibit_D_change$cluster <- clusters$cluster

# Examine the first few rows of the data frame to see the cluster assignments
head(exhibit_D_change)

output_path <- "exhibit_D_change_clustering.csv" # You can change this to the desired output path
write.csv(exhibit_D_change, output_path, row.names = FALSE)

cat("Change in exhibit D factors by condition saved to:", output_path)

#check for clustering tendency of my data
# Install factoextra package if it's not already installed
if (!require(factoextra)) install.packages("factoextra")
library(factoextra)

datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)
exhibit_D_change <- read.csv(paste(datadir, '/exhibit_D_change.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)

columns_to_keep <- c(6:22)
data_for_clustering <- exhibit_D_change[, columns_to_keep]

# Calculate the Hopkins statistic
set.seed(42)  # Set a random seed for reproducibility
hopkins_stat <- get_clust_tendency(data_for_clustering, n = 7)

# Print the Hopkins statistic
print(hopkins_stat$hopkins_stat)

#**************************************DBSCAN (Density-Based Spatial Clustering of Applications with Noise) 11/15/2023
if (!require(dbscan)) install.packages("dbscan")
library(dbscan)

datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)
exhibit_D_change <- read.csv(paste(datadir, '/exhibit_D_change.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)

# Remove non-numeric columns if necessary (e.g., subject identifier)
# Assuming the first column 'subject' is an identifier and should be excluded
#data_for_clustering <- exhibit_D_change[, -c(1:6)]
data_for_clustering <- exhibit_D_change[,6:22]
data_scaled <- scale(data_for_clustering)

# Choose k - typically minPts - 1
k <- 5  # Adjust based on your choice of minPts

# Compute the k-distance
k_distance <- kNNdist(data_scaled, k = k)

# Sort the distances in descending order
k_distance_sorted <- sort(k_distance, decreasing = TRUE)

# Plot the k-distance graph
plot(k_distance_sorted, main = "k-Distance Graph", xlab = "Points sorted by distance", ylab = "k-distance")

# Run DBSCAN
# eps and minPts are crucial parameters and need to be chosen carefully, use the elbow of the k-distance graph to determine eps.
# These are example values and should be tuned for your specific dataset
set.seed(42)  # For reproducibility
dbscan_result <- dbscan(data_scaled, eps = 4.5 , minPts = 2)

# View the cluster assignments
print(dbscan_result)

# Add cluster assignments to your data
data_for_clustering$cluster <- dbscan_result$cluster

# Examine the modified data frame
head(data_for_clustering)
# if a cluster asignment is 0. then it means that the algorithm has labeled that point as noise. 
output_path <- "exhibit_D_change_DBSCAN_clustering.csv" # You can change this to the desired output path
write.csv(data_for_clustering, output_path, row.names = FALSE)

cat("DBSCAN clustering saved to:", output_path)





######## START MULTIPLE LINEAR REGRESSION, LASSO TECHNIQUE **********************************

install.packages("glmnet")
library(glmnet)

#create df from PVT data, the response variable
response_vector <- c(85.1100134, 107.2399866, 281.3100049, -170.5200061, 50.8999939, 15.2599946, 26.3700104, 17.4700164)

#create df for predictors, which will be the change in PC value from pre and post, for both gene and ROI
predictors <- read.csv('predictors.csv', header = TRUE, stringsAsFactors = FALSE) 

#make sure to code categorical variables as dummy variables

#fittig Lasso model
set.seed(123) # Set seed for reproducibility

# cv.glmnet() function performs cross-validation to choose the lambda parameter that results in the best model fit
cvfit <- cv.glmnet(as.matrix(predictors), response_vector, alpha = 1)
# The `alpha = 1` argument specifies Lasso regression. alpha = 0 argument specifies ridge. anything in between represents an Elastic Net

plot(cvfit)

#determine best lambda
best_lambda <- cvfit$lambda.min

#run model with best lambda
lasso_model <- glmnet(as.matrix(predictors), response_vector, alpha = 1, lambda = best_lambda)

#see coefficients of model
coef(lasso_model)


#make predictions with model, needs new observations
#predictions <- predict(lasso_model, newx = new_predictors_matrix)
#mse <- mean((predictions - actual_response)^2)





#replacing probe IDs from Ram's GE PCA loadings into gene symbols

if (!requireNamespace("biomaRt", quietly = TRUE))
    install.packages("biomaRt")
library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
probe_ids <- unique(c(df$PC1_Rownames, df$PC2_Rownames, ...))  # Include all columns with probe IDs

gene_info <- getBM(attributes = c('affy_hg_u133_plus_2', 'hgnc_symbol'),
                   filters = 'affy_hg_u133_plus_2',
                   values = probe_ids,
                   mart = mart)


gene_symbols <- setNames(gene_info$hgnc_symbol, gene_info$affy_hg_u133_plus_2)
df$PC1_Rownames <- gene_symbols[df$PC1_Rownames]
df$PC2_Rownames <- gene_symbols[df$PC2_Rownames]
# Repeat for all PC columns







if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("hgu133a2.db")



if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")



replace_probe_with_symbol <- function(df, num_pc, gene_symbols) {
  # Loop over the range of PCs
  for (i in 1:num_pc) {
    # Create the column name for the current PC
    pc_colname <- paste0("PC", i, "_Rownames")
    
    # Check if the column exists in the dataframe
    if (!pc_colname %in% names(df)) {
      warning(paste("Column", pc_colname, "does not exist in dataframe. Skipping."))
      next
    }
    
    # Replace the probe IDs with gene symbols for the current PC
    # It uses the gene_symbols named vector for the replacement
    df[[pc_colname]] <- sapply(df[[pc_colname]], function(probe_id) {
      if (!is.na(probe_id) && probe_id %in% names(gene_symbols)) {
        return(gene_symbols[probe_id])
      } else {
        return(NA)  # Return NA if the probe_id is not found in the gene_symbols
      }
    })
  }
  return(df)
}

datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)
actual_pca_loadings <- read.csv(paste(datadir, '/actual_pca_loadings.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)


replace_probe_with_symbol(df, 8, actual_pca_loadings)

# Assume gene_symbols is a named vector where names are probe IDs and values are gene symbols
# Example usage:
# df <- replace_probe_with_symbol(df, 8, gene_symbols)











datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)
GE_PC1 <- read.csv(paste(datadir, '/GE_PC1.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)

#download dataset and check contents
BiocManager::install("GEOquery")
require(GEOquery)
BiocManager::install("Biobase")
require(Biobase)
gset <- getGEO("GSE12056", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

dim(exprs(gset))

rownames(exprs(gset))[1:50]


#annotate these first 50 and create a 'lookup' table of annotation that can be used to rename your Affy IDs to gene names
#require("biomaRt")
#mart <- useMart("ENSEMBL_MART_ENSEMBL")
#mart <- useDataset("hsapiens_gene_ensembl", mart)
#annotLookup <- getBM(
  #mart=mart,
  #attributes=c(
    #"affy_hg_u133_plus_2",
    #"ensembl_gene_id",
    #"gene_biotype",
    #"external_gene_name"),
  #filter = "affy_hg_u133_plus_2",
  #values = rownames(exprs(gset))[1:50], uniqueRows=TRUE)

#head(annotLookup, 20)

#indicesLookup <- match(rownames(gset), annotLookup$affy_hg_u133_plus_2)
#head(annotLookup[indicesLookup, "external_gene_name"])

#double check the order is correct
#dftmp <- data.frame(rownames(gset), annotLookup[indicesLookup, c("affy_hg_u133_plus_2", "external_gene_name")])
#head(dftmp, 20)

#table(dftmp[,1] == dftmp[,2])

#replace row names
#rownames(gset) <- paste(annotLookup[indicesLookup, "external_gene_name"], c(1:length(indicesLookup)), sep="_")





if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("hgu133a2.db")



if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")

library(biomaRt)

# Assuming you've already connected to the BioMart database and selected the dataset
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)

# Set the directory where your "GE_PC1.csv" is located
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)

# Load your "GE_PC1.csv" file
GE_PC1 <- read.csv(paste(datadir, '/GE_PC1.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)

# Retrieve annotations for all probe IDs in GE_PC1
annotLookup <- getBM(
  attributes = c(
    "affy_hg_u133_plus_2",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "affy_hg_u133_plus_2",
  values = GE_PC1$PC1_Rownames, # Use the probe IDs from your file
  mart = mart,
  uniqueRows = TRUE)

# Match the probe IDs from your GE_PC1 with the annotations
indicesLookup <- match(GE_PC1$PC1_Rownames, annotLookup$affy_hg_u133_plus_2)

# Check for any probes not found in the annotation data
missing_probes <- is.na(indicesLookup)
if (any(missing_probes)) {
  warning("Some probes are not found in the annotation data and will be left as is.")
}

# Replace the probe IDs in 'PC1_Rownames' with the corresponding gene symbols
# For probes not found in 'annotLookup', keep the original probe ID
GE_PC1$Gene_Symbols <- annotLookup$external_gene_name[indicesLookup] # Create a new column for gene symbols
GE_PC1$Gene_Symbols[missing_probes] <- GE_PC1$PC1_Rownames[missing_probes] # Fill in missing values with original probe IDs

# Now GE_PC1 has a new column 'Gene_Symbols' with gene symbols where matches were found, and the original probe IDs otherwise
# The order of rows in GE_PC1 is preserved



annotate_probes <- function(num_files, directory) {
  # Load required libraries
  library(biomaRt)
  
  # Connect to the BioMart database
  mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  
  # Iterate through the files
  for (i in 1:num_files) {
    # Construct the file name
    file_name <- file.path(directory, paste("GE_PC", i, ".csv", sep = ""))
    
    # Read the CSV file
    GE_data <- read.csv(file_name, header = TRUE, stringsAsFactors = FALSE)
    
    # Retrieve annotations for all probe IDs in GE_data
    annotLookup <- getBM(
      attributes = c(
        "affy_hg_u133_plus_2",
        "ensembl_gene_id",
        "gene_biotype",
        "external_gene_name"),
      filter = "affy_hg_u133_plus_2",
      values = GE_data$PC1_Rownames, # Assumes the column name is the same across files
      mart = mart,
      uniqueRows = TRUE)
    
    # Match the probe IDs with the annotations
    indicesLookup <- match(GE_data$PC1_Rownames, annotLookup$affy_hg_u133_plus_2)
    
    # Check for any probes not found in the annotation data
    missing_probes <- is.na(indicesLookup)
    if (any(missing_probes)) {
      warning(paste("Some probes are not found in the annotation data for file", i, "and will be left as is."))
    }
    
    # Replace the probe IDs with the corresponding gene symbols
    GE_data$Gene_Symbols <- annotLookup$external_gene_name[indicesLookup]
    GE_data$Gene_Symbols[missing_probes] <- GE_data$PC1_Rownames[missing_probes]
    
    # Save the annotated data back to a new CSV file
    write.csv(GE_data, file.path(directory, paste("Annotated_GE_PC", i, ".csv", sep = "")), row.names = FALSE)
  }
}

# Usage:
# Specify the number of files and the directory where the files are located
annotate_probes(num_files = 10, directory = 'C:/Users/lilyb/Desktop/New folder')








#the same ANOVA with eta squared but for gene expression data
library(lme4)
library(car)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")

BiocManager::install("qvalue")
library(qvalue)

if (!require("lmerTest")) {
    install.packages("lmerTest")
}
library(lmerTest)

simple1 <- function(gene){
  
  # Subset data for the current brain region
  sub_data <- df[df$Gene == gene, ]

  if(nrow(sub_data) == 0){
    return(NA)  # Returns NA if no data for the brain region
  }

  # Set up model formula
  model_formula <- Expression ~ Gender * Condition + (1|Patient_ID)
  
  # Define contrasts for type III
  type3_contrasts <- list(Patient_ID = contr.sum, Gender = contr.sum, Condition = contr.sum)
  
  # Fit the mixed model with explicit contrasts
  mixed_model <- lmer(model_formula, data = sub_data, contrasts = type3_contrasts)
  mixed_anova <- Anova(mixed_model, test="F", type=3)

  # Convert the ANOVA output to a dataframe
  anova_values <- data.frame(mixed_anova)

  # Check column names
  print(colnames(anova_values)) 

  residual_df_colname <- 'Df.res' # Replace with the correct column name if needed

  # Calculate eta squared
  anova_values$eta_sq <- with(anova_values, ifelse(F > 0, (Df * F) / (Df * F + get(residual_df_colname)), 0))

  # Replace NA with 0 in eta squared
  anova_values$eta_sq[is.na(anova_values$eta_sq)] <- 0

  # Add brain region to the data frame
  anova_values$Gene <- gene
  
  # Return as a row
  return(t(anova_values))
}

#CONDITION GENES ANOVA with eta squared
# set directory
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)

# Read in the CSV data
df <- read.csv('condition_genes.csv', header = TRUE, stringsAsFactors = FALSE) 

#condition genes
brain_regions <- unique(df$Gene)
lmer_list <- lapply(brain_regions, simple1)
lmer_df <- do.call(rbind, lmer_list)

#save results
write.csv(lmer_df, "gene_anova_cond.csv", row.names = TRUE)


#GENDER GENES ANOVA with eta squared
# set directory
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)

# Read in the CSV data
df <- read.csv('Gender_genes.csv', header = TRUE, stringsAsFactors = FALSE) 

# gender genes
brain_regions <- unique(df$Gene)
lmer_list <- lapply(brain_regions, simple1)
lmer_df <- do.call(rbind, lmer_list)

#save results
write.csv(lmer_df, "gene_anova_gender.csv", row.names = TRUE)

# Read in the CSV data
df <- read.csv('Interaction_genes.csv', header = TRUE, stringsAsFactors = FALSE) 


brain_regions <- unique(df$Gene)
lmer_list <- lapply(brain_regions, simple1)
lmer_df <- do.call(rbind, lmer_list)

#save results
write.csv(lmer_df, "gene_anova_interaction.csv", row.names = TRUE)







library(lme4)
library(car)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")

BiocManager::install("qvalue")
library(qvalue)

if (!require("lmerTest")) {
    install.packages("lmerTest")
}
library(lmerTest)

simple2 <- function(Principal_Component, df){
  
  # Subset data
  sub_data <- df[df$Principal_Component == Principal_Component, ]
  sub_data$subject <- factor(sub_data$subject)
  sub_data$Gender <- factor(sub_data$Gender)
  sub_data$Condition <- factor(sub_data$Condition)

  if(nrow(sub_data) == 0){
    return(NA)  # Returns NA if no data for the brain region
  }

  # Set up model formula
  model_formula <- as.formula(paste("Value ~ Gender * Condition + (1|subject)"))
  
  # Define contrasts for type III
  type3_contrasts <- list(Gender = contr.sum, Condition = contr.sum)

  # Fit the mixed model with explicit contrasts
  mixed_model <- lmer(model_formula, data = sub_data, contrasts = type3_contrasts)
  mixed_anova <- Anova(mixed_model, test="F", type=3)

  # Convert the ANOVA output to a dataframe
  anova_values <- data.frame(mixed_anova)

  # Check column names
  print(colnames(anova_values)) 

  residual_df_colname <- 'Df.res' # Replace with the correct column name if needed

  # Calculate eta squared
  anova_values$eta_sq <- with(anova_values, ifelse(F > 0, (Df * F) / (Df * F + get(residual_df_colname)), 0))

  # Replace NA with 0 in eta squared
  anova_values$eta_sq[is.na(anova_values$eta_sq)] <- 0

  # Add brain region to the data frame
  anova_values$Principal_Component <- Principal_Component
  
  # Return as a row
  return(t(anova_values))
}


# set directory
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)

# Read in the CSV data
df <- read.csv('genePCA_long.csv', header = TRUE, stringsAsFactors = TRUE) 
# Create a list of unique principal components
Principal_components <- unique(df$Principal_Component)

# Perform analysis for each principal component
lmer_results <- lapply(Principal_components, simple2, df = df)

# Combine all results into a single data frame
lmer_df <- do.call(rbind, lmer_results)

# Save the results to a CSV file
write.csv(lmer_df, "genePCA_anova.csv", row.names = TRUE)





#plotting TUBB2A and NEBL against Sex and condition
install.packages("ggplot2")
library(ggplot2)
# set directory
datadir <- file.path('C:/Users/lilyb/Desktop/New folder')
setwd(datadir)

#fill out data frame
data_TUBB2A <- read.csv('TUBB2A.csv')
data_TUBB2A$Gender <- as.factor(data_TUBB2A$Gender)

ggplot(data_TUBB2A, aes(x = Condition, y = GeneExpression, color = Gender)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x * Gender, se = FALSE) +  # Interaction term added here
    theme_minimal() +
    labs(title = "Gene Expression by Condition and Gender",
         x = "Condition",
         y = "Gene Expression Level TUBB2A")

data_NEBL <- read.csv('NEBL.csv')
data_NEBL$Gender <- as.factor(data_NEBL$Gender)

ggplot(data_NEBL, aes(x = Condition, y = GeneExpression, color = Gender)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x * Gender, se = FALSE) +  # Interaction term added here
    theme_minimal() +
    labs(title = "Gene Expression by Condition and Gender",
         x = "Condition",
         y = "Gene Expression Level NEBL")
        
        
# installs required libraries
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")

# do principal component analysis
library(PCAtools)

###
# Load DE_subset_40000 and applied_subset_40000
# before running this.
###

DE_subset <- DE_subset_40000
applied_subset <- applied_subset_40000
generow_names <- generow_names_40000


# create matrices with pre and SD values for the gene expression (GE) values
pre_GE_matrix <- data.matrix(applied_subset[,1:8 * 2])
SD_GE_matrix <- data.matrix(applied_subset[,1 + 0:7 * 2])
all_gene_matrix <- data.matrix(applied_subset)

# create matrices for female and male increase in gene expression from before to after sleep deprivation (SD - pre)
female_gene_change_matrix <- data.matrix(data.frame('P01_Increase' = DE_subset[,1] - DE_subset[,2], 'P02_Increase' = DE_subset[,3] - DE_subset[,4], 'P03_Increase' = DE_subset[,5] - DE_subset[,6], 'P04_Increase' = DE_subset[,7] - DE_subset[,8]))
male_gene_change_matrix <- data.matrix(data.frame('P05_Increase' = DE_subset[,9] - DE_subset[,10], 'P06_Increase' = DE_subset[,11] - DE_subset[,12], 'P07_Increase' = DE_subset[,13] - DE_subset[,14], 'P08_Increase' = DE_subset[,15] - DE_subset[,16]))
all_gene_change_matrix <- data.matrix(data.frame(female_gene_change_matrix, male_gene_change_matrix))

# rename the matrix of which you want a pca
# into "data"
data <- all_gene_matrix

# Perform PCA
# the test subjects are the variables, which have to be in rows
# the genes are the samples, which have to be in columns
pca_result <- pca(t(data), center = TRUE, scale = TRUE)
pca_result_loadings <- pca_result$loadings

# Start Using Residuals #

# patient_order <- c(01Post.CEL, 01Pre.CEL, 02Post.CEL, 02Pre.CEL, 03Post.CEL, 03Pre.CEL, 04Post.CEL, 04Pre.CEL, 05Post.CEL, 05Pre.CEL, 06Post.CEL, 06Pre.CEL, 07Post.CEL, 07Pre.CEL, 08Post.CEL, 08Pre.CEL)
condition <- factor(rep(c('Post', 'Pre'), 8)) #Pre or Post Treatment
gender <- factor(c(rep('F', 8), rep('M', 8))) #gender
affyRNAdeg_slope <- c(2.732662, 1.353166, 2.194767, 1.794559, 2.027774, 1.623671, 2.43369, 2.080824, 2.582972, 2.25435, 2.263606, 1.691016, 1.323386, 1.828324, 1.843323, 1.019396)

# gets the value of the residuals for 1 row
get_residual_value <- function(generow){
	# Set up model formula
	model_formula <- genes ~ degradation
	lm_data <- data.frame(degradation=affyRNAdeg_slope,genes=generow)

	# Fit the mixed model
	mixed_model <- lm(model_formula, data = lm_data)
	summary_resid <- summary(mixed_model)$residuals
	return(summary_resid)
}

# applied_subset_fp <- file.path('C:/Users/Lily/Desktop/gene_anova/applied_subset_40000.csv')
# Set directory
datadir <- file.path('C://Users/lilyb/Desktop/New folder/eh')
setwd(datadir)
applied_subset <- read.csv('221728_x_at.csv', row.names=1, header = TRUE, stringsAsFactors = FALSE)
# applied_subset$ProbeID <-NULL

# gets the value of the residuals for all the rows
genes_degradation_residuals <- t(apply(applied_subset,1,get_residual_value))

# Perform PCA
# the test subjects are the variables, which have to be in rows
# the genes are the samples, which have to be in columns
affyRNAdeg_pca_result <- pca(t(genes_degradation_residuals), center = TRUE, scale = TRUE)
affyRNAdeg_pca_result_rotated <- data.frame(affyRNAdeg_pca_result$rotated, generow_names)

# get all the values

pc1_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,1], decreasing=TRUE),],50)
pc2_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,2], decreasing=TRUE),],50)
pc3_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,3], decreasing=TRUE),],50)
pc4_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,4], decreasing=TRUE),],50)
pc5_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,5], decreasing=TRUE),],50)
pc6_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,6], decreasing=TRUE),],50)
pc7_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,7], decreasing=TRUE),],50)
pc8_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,8], decreasing=TRUE),],50)
pc9_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,9], decreasing=TRUE),],50)
pc10_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,10], decreasing=TRUE),],50)
pc11_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,11], decreasing=TRUE),],50)
pc12_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,12], decreasing=TRUE),],50)
pc13_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,13], decreasing=TRUE),],50)
pc14_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,14], decreasing=TRUE),],50)
pc15_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,15], decreasing=TRUE),],50)
pc16_values <- head(affyRNAdeg_pca_result_rotated[order(affyRNAdeg_pca_result_rotated[,16], decreasing=TRUE),],50)


actual_pca_loadings <- data.frame(
"PC1_Rownames" = pc1_values$Genename, "PC1_Values"=pc1_values[,1],
"PC2_Rownames" = pc2_values$Genename, "PC2_Values"=pc2_values[,2],
"PC3_Rownames" = pc3_values$Genename, "PC3_Values"=pc3_values[,3],
"PC4_Rownames" = pc4_values$Genename, "PC4_Values"=pc4_values[,4],
"PC5_Rownames" = pc5_values$Genename, "PC5_Values"=pc5_values[,5],
"PC6_Rownames" = pc6_values$Genename, "PC6_Values"=pc6_values[,6],
"PC7_Rownames" = pc7_values$Genename, "PC7_Values"=pc7_values[,7],
"PC8_Rownames" = pc8_values$Genename, "PC8_Values"=pc8_values[,8],
"PC9_Rownames" = pc9_values$Genename, "PC9_Values"=pc9_values[,9],
"PC10_Rownames" = pc10_values$Genename, "PC10_Values"=pc10_values[,10],
"PC11_Rownames" = pc11_values$Genename, "PC11_Values"=pc11_values[,11],
"PC12_Rownames" = pc12_values$Genename, "PC12_Values"=pc12_values[,12],
"PC13_Rownames" = pc13_values$Genename, "PC13_Values"=pc13_values[,13],
"PC14_Rownames" = pc14_values$Genename, "PC14_Values"=pc14_values[,14],
"PC15_Rownames" = pc15_values$Genename, "PC15_Values"=pc15_values[,15],
"PC16_Rownames" = pc16_values$Genename, "PC16_Values"=pc16_values[,16]
)

# End Using Residuals #

# Cumulative Variance Plot #

get_top8_cum_var <- function(data){
	# Perform PCA
	pca_result <- prcomp(t(data), scale. = TRUE)
	
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

get_top8_components <- function(data){
	#Perform PCA
	pca_result <- prcomp(t(data), scale. = TRUE)
	
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


# Start Making ANOVA #

# These stuff are important.
library(lme4)
library(car)

# Reminder:
patient <- factor(rep(c(1,2,3,4,5,6,7,8), each=2)) # patient ID
condition <- factor(rep(c('Post', 'Pre'), 8)) #Pre or Post Treatment
gender <- factor(c(rep('F', 8), rep('M', 8))) #gender

# finds the values indicated using through lmer and mean
# anova_values contains the six values: 'Gender F', 'Condition F', 'Gender:Condition F', 'Gender Pr(>F)', 'Condition Pr(>F)', 'Gender:Condition Pr(>F)'
# 14 desired values: 'Gender F', 'Condition F', 'Gender:Condition F', 'Gender Pr(>F)', 'Condition Pr(>F)', 'Gender:Condition Pr(>F)', 'Mean Male', 'Mean Female', 'Mean Pre-SD', 'Mean Post-SD', 'Mean Male Pre-SD', 'Mean Male Post-SD', 'Mean Female Pre-SD', 'Mean Female Post-SD'
get_lmer_value <- function(generow){
	# Set up model formula
	model_formula <- Generow ~ Gender * Condition + (1|PatientID)
	lmer_data <- data.frame(PatientID=patient,Condition=condition,Gender=gender,Generow=generow)
	# Needed for Type III: https://rdoodles.rbind.io/2020/10/type-3-anova-in-r-an-easy-way-to-publish-wrong-tables/
	type3 <- list(patient=contr.sum, gender=contr.sum, condition=contr.sum, generow=contr.sum)
	# Fit the mixed model
	mixed_model <- lmer(model_formula, data = lmer_data, contrasts=type3)
	mixed_anova <- Anova(mixed_model, test="F", type=3)
	# Get f and Pr(>F), 1=F value, 4=p value of R. unlist creates vector instead of list, c indicates vector
	anova_values <- unlist(mixed_anova[2:4,c(1,4)])
	# Put these values and the other values in a single argument
	important_data <- data.frame(c(anova_values, mean(generow[1:8]), mean(generow[9:16]), mean(generow[1:8 * 2]), mean(generow[1:8 * 2 - 1]), mean(generow[1:4 * 2]), mean(generow[1:4 * 2 - 1]), mean(generow[5:8 * 2]), mean(generow[5:8 * 2 - 1])))
	# have to return this as a row, so transpose it
	return(t(important_data))
}

#this one attempts to correct for unequal variances by applying log transformation
get_lmer_value2 <- function(generow){
    # Apply logarithmic transformation to Generow
    # Avoiding log(0) issues by adding a small constant (1e-6)
    generow <- log(generow + 1e-6)
    
    # Set up model formula
    model_formula <- generow ~ Gender * Condition + (1|PatientID)
    
    # Prepare the data frame with the transformed Generow
    lmer_data <- data.frame(PatientID=patient, Condition=condition, Gender=gender, Generow=generow)
    
    # Needed for Type III contrasts
    type3 <- list(patient=contr.sum, gender=contr.sum, condition=contr.sum, generow=contr.sum)
    
    # Fit the mixed model
    mixed_model <- lmer(model_formula, data = lmer_data, contrasts=type3)
    mixed_anova <- Anova(mixed_model, test="F", type=3)
    
    # Extract F value and p-value from the ANOVA table
    anova_values <- unlist(mixed_anova[2:4,c(1,4)])
    
    # Calculating means of original Generow data for specific subgroups
    important_data <- data.frame(c(anova_values, 
                                   mean(generow[1:8]), mean(generow[9:16]), 
                                   mean(generow[1:8 * 2]), mean(generow[1:8 * 2 - 1]), 
                                   mean(generow[1:4 * 2]), mean(generow[1:4 * 2 - 1]), 
                                   mean(generow[5:8 * 2]), mean(generow[5:8 * 2 - 1])))
    
    # Return as a row, so transpose it
    return(t(important_data))
}


# have to transpose the matrix created by this (long time, ignore warnings like "boundary (singular) fit: see ?isSingular")
residuals_lmer_df <- t(apply(genes_degradation_residuals,1,get_lmer_value))

# assign column names based on the value represented
colnames(residuals_lmer_df) <- c('Gender F', 'Condition F', 'Gender:Condition F', 'Gender Pr(>F)', 'Condition Pr(>F)', 'Gender:Condition Pr(>F)', 'Mean Female', 'Mean Male', 'Mean Pre-SD', 'Mean Post-SD', 'Mean Female Pre-SD', 'Mean Female Post-SD', 'Mean Male Pre-SD', 'Mean Male Post-SD')

# find bonferroni values
bonferroni_values <- data.frame("Gender Pr(>F)" = p.adjust(residuals_lmer_df[,4]), "Condition Pr(>F)" = p.adjust(residuals_lmer_df[,5]), "Gender:Condition Pr(>F)" = p.adjust(residuals_lmer_df[,6]))

# End Making ANOVA #

# Start Using StoreyQ #

# Use code from https://github.com/StoreyLab/qvalue
# This stuff installs it. You should only run it if library(qvalue) doesn't work
#install.packages("devtools")
#library("devtools")
#install_github("jdstorey/qvalue")

library(qvalue)

# Finds Storey q-values of the p-values
residuals_lmer_df_qvalues <- qvalue(residuals_lmer_df[,4:6])$qvalues

# End Using StoreyQ #
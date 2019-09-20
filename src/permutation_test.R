library(dplyr)
library(caret)
library(pracma)
library(kernlab)
library(lattice)

# wd is working directory
# dw is download directory

arcsinh <- function(x) log(x+sqrt(x*x+1))

permutation_test <- function(wd, dw, samples, tries = 10, cores = 28){
  
  deps <- c('dplyr', 'caret', 'pracma', 'kernlab', 'lattice')

	setwd(wd)

	z <- read.csv("arcsinh_edata.csv", header=TRUE, row.names = 1)
	c <- arcsinh(dplyr::select(read.csv("full_model.csv", header = TRUE, row.names = 1), lifespan))
	df <- read.csv("cluster_table_classic_arcsinh.csv", header=TRUE, row.names = 1)
	scores <- read.csv("resultados_arcsinh_DESeq2.csv", header=TRUE, row.names = 1)
	colnames(z) <- df$gene

	colors_list <- unique(scores$cluster_color)

	permutation_table <- data.frame(cluster_color = character(), signature = character(), r2_svm = numeric(), sigma = numeric(), C = numeric())

	test_rows <- which(rownames(z) %in% samples)

	z_test <- z[test_rows, ]
	c_test <- c[test_rows, ]
	c_test <- as.data.frame(arcsinh(c_test))
	names(c_test) <- c("lifespan")
	rownames(c_test) <- rownames(z_test)
		
	z_train <- z[-test_rows, ]
	c_train <- c[-test_rows, ]
	c_train <- as.data.frame(arcsinh(c_train))
	names(c_train) <- c("lifespan")
	rownames(c_train) <- rownames(z_train)

	
	library(dplyr)
	library(caret)
	library(pracma)
	library(kernlab)
	library(lattice)
	
	
	for(counter in 1:tries){
		
		print("iteration number: ", counter)
		
		list_ranks <- data.frame(cluster_color = character(), signature = character(), r2_svm = numeric(), sigma = numeric(), C = numeric())

		for(pig in colors_list){

			cluster <- dplyr::filter(df, clusterColor == pig)$gene
			z_norm <- dplyr::select(z_train, as.character(cluster))

			target <- c_train[,1]

			n_genes <- length(unlist(strsplit(as.character(scores[scores$cluster_color == pig, ]$signature),  ",")))

			sig <- sample(colnames(z_norm), n_genes)

			TrainCtrl1 <- trainControl(method = "repeatedcv", number = 30, repeats=5,verbose = FALSE)
			SVMgrid <- expand.grid(sigma = linspace(.03, 0.9, n = 15), C = linspace(1, 5, n = 15))
			modelSvmRRB <- train(z_norm[, sig], as.numeric(c_train$lifespan), method="svmRadial", trControl=TrainCtrl1, tuneGrid = SVMgrid, verbose=FALSE)

			TrainCtrl1 <- trainControl(method = "repeatedcv", number = 30,repeats=5,verbose = FALSE)
			SVMgrid <- expand.grid(sigma = modelSvmRRB$bestTune$sigma, C = modelSvmRRB$bestTune$C)
			modelSvmRRB <- train(z_norm[, sig], as.numeric(c_train$lifespan), method="svmRadial", trControl=TrainCtrl1, tuneGrid = SVMgrid, verbose=FALSE)

			z_norm_test <- dplyr::select(z_test, as.character(cluster))

			PredictedTest <- predict(modelSvmRRB, z_norm_test[, sig])
			real <- c_test

			list_ranks <- rbind(list_ranks, data.frame(cluster_color = pig, signature = paste(sig, collapse = ","), r2_svm = cor(real$lifespan,PredictedTest), sigma = modelSvmRRB$bestTune$sigma, C = modelSvmRRB$bestTune$C))
		}    
			permutation_table <- rbind(permutation_table, list_ranks)
			setwd(dw)
			write.csv(permutation_table, paste("permutation_table_X",counter,".csv", sep=""))		
	}

}

samples = c("counts_mauratus_liver", "counts_mauratus_kidney",
         "counts_fdamarensis_brain", "counts_fdamarensis_kidney",
         "counts_fdamarensis_liver", "counts_uamericanus_liver", "counts_uamericanus_kidney",
         "counts_uamericanus_brain", "counts_hsapiens_liver",
         "counts_hsapiens_kidney", "counts_hsapiens_heart",
         "counts_hsapiens_lung")
data = "/home/rstudio/DESeq2_preprocessing"
output = "/data/results/cross-species"
tries = 10
permutation_test(data, output, samples, tries)
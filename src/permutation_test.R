library(dplyr)
library(caret)
library(pracma)
library(kernlab)
library(lattice)


permutation_test <- function(wd,dw, iter){ 

	setwd(wd)

	z <- read.csv("ORQ_edata.csv", header=TRUE, row.names = 1)
	c <- log(dplyr::select(read.csv("full_model.csv", header = TRUE, row.names = 1), lifespan))
	df <- read.csv("cluster_table_classic_ORQ.csv", header=TRUE, row.names = 1)
	scores <- read.csv("resultados_ORQ_DESeq2_log.csv", header=TRUE, row.names = 1)
	colnames(z) <- df$gene

	colors_list <- unique(scores$cluster_color)

	permutation_table <- data.frame(cluster_color = character(), signature = character(), r2_svm = numeric(), sigma = numeric(), C = numeric())

	test_rows <- which(rownames(z) %in% c("counts_mauratus_liver", "counts_mauratus_kidney", 
										  "counts_fdamarensis_brain", "counts_fdamarensis_kidney", 
										  "counts_fdamarensis_liver", "counts_uamericanus_liver", "counts_uamericanus_kidney", 
										  "counts_uamericanus_brain", "counts_hsapiens_liver", 
										  "counts_hsapiens_kidney", "counts_hsapiens_heart", 
										  "counts_hsapiens_lung"))

	z_test <- z[test_rows, ]
	c_test <- c[test_rows, ]
	c_test <- as.data.frame(log(c_test))
	names(c_test) <- c("lifespan")
	rownames(c_test) <- rownames(z_test)
		
	z_train <- z[-test_rows, ]
	c_train <- c[-test_rows, ]
	c_train <- as.data.frame(log(c_train))
	names(c_train) <- c("lifespan")
	rownames(c_train) <- rownames(z_train)


	for(counter in 1:iter){
		
		list_ranks <- data.frame(cluster_color = character(), signature = character(), r2_svm = numeric(), sigma = numeric(), C = numeric())

		for(pig in colors_list){

			cluster <- dplyr::filter(df, clusterColor == pig)$gene
			z_norm <- dplyr::select(z_train, as.character(cluster))

			target <- c_train[,1]

			n_genes <- length(unlist(strsplit(as.character(scores[scores$cluster_color == pig, ]$signature),  ",")))

			sig <- sample(colnames(z_norm), n_genes)

			TrainCtrl1 <- trainControl(method = "repeatedcv", number = 30, repeats=5,verbose = FALSE)
			SVMgrid <- expand.grid(sigma = linspace(.03, 0.9, n = 20), C = linspace(1, 5, n = 30))
			modelSvmRRB <- train(z_norm[, sig], as.numeric(c_train$lifespan), method="svmRadial", trControl=TrainCtrl1, tuneGrid = SVMgrid, verbose=FALSE)

			TrainCtrl1 <- trainControl(method = "repeatedcv", number = 30,repeats=5,verbose = FALSE)
			SVMgrid <- expand.grid(sigma = modelSvmRRB$bestTune$sigma, C = modelSvmRRB$bestTune$C)
			modelSvmRRB <- train(z_norm[, sig], as.numeric(c_train$lifespan), method="svmRadial", trControl=TrainCtrl1, tuneGrid = SVMgrid, verbose=FALSE)

			z_norm_test <- dplyr::select(z_test, as.character(cluster))

			PredictedTest <- exp(predict(modelSvmRRB, z_norm_test[, sig]))
			real <- exp(c_test)

			list_ranks <- rbind(list_ranks, data.frame(cluster_color = pig, signature = paste(sig, collapse = ","), r2_svm = cor(real$lifespan,PredictedTest), sigma = modelSvmRRB$bestTune$sigma, C = modelSvmRRB$bestTune$C))
		}    
			permutation_table <- rbind(permutation_table, list_ranks)
		
	}

	setwd(dw)

	write.csv(permutation_table, "permutation_table_X10.csv")

}

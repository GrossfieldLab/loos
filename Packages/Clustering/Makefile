IEIGEN=-I/anaconda3/envs/loos/include
CC=clang++

cluster-kgs: cluster-kgs.o Clustering.o HAC.o AverageLinkage.o ClusteringUtils.o
	$(CC) -O3 -o cluster-kgs $(IEIGEN) cluster-kgs.o Clustering.o HAC.o AverageLinkage.o ClusteringUtils.o

cluster-kgs.o: cluster-kgs.cpp Clustering.hpp ClusteringUtils.hpp 
	$(CC) -c -O3 $(IEIGEN) cluster-kgs.cpp

Clustering.o: ClusteringUtils.hpp HAC.hpp AverageLinkage.hpp KGS.hpp
	$(CC) -c -O3 Clustering.hpp

ClusteringUtils.o: Clustering.hpp ClusteringUtils.hpp
	$(CC) -c -O3 $(IEIGEN) ClusteringUtils.cpp

HAC.o: Clustering.hpp HAC.hpp HAC.cpp
	$(CC) -c -O3 $(IEIGEN) HAC.cpp

AverageLinkage.o: Clustering.hpp HAC.hpp AverageLinkage.hpp AverageLinkage.cpp
	$(CC) -c -O3 $(IEIGEN) AverageLinkage.cpp

KGS.o: Clustering.hpp HAC.hpp Averagelinkage.hpp 
	$(CC) -c $(IEIGEN) -O3 KGS.cpp


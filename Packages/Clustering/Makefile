IEIGEN=-I/home/louis/anaconda3/envs/loos/include

cluster-kgs: cluster-kgs.o Clustering.o HAC.o AverageLinkage.o ClusteringUtils.o
	g++ -o cluster-kgs $(IEIGEN) -O3 cluster-kgs.o Clustering.o HAC.o AverageLinkage.o ClusteringUtils.o

cluster-kgs.o: cluster-kgs.cpp Clustering.hpp ClusteringUtils.hpp 
	g++ -c -O3 $(IEIGEN) cluster-kgs.cpp

Clustering.o: ClusteringUtils.hpp HAC.hpp AverageLinkage.hpp KGS.hpp
	g++ -c -O3 Clustering.hpp

ClusteringUtils.o: Clustering.hpp ClusteringUtils.hpp
	g++ -c -O3 $(IEIGEN) ClusteringUtils.cpp

HAC.o: Clustering.hpp HAC.hpp HAC.cpp
	g++ -c -O3 $(IEIGEN) HAC.cpp

AverageLinkage.o: Clustering.hpp HAC.hpp AverageLinkage.hpp AverageLinkage.cpp
	g++ -c -O3 $(IEIGEN) AverageLinkage.cpp

KGS.o: Clustering.hpp HAC.hpp Averagelinkage.hpp 
	g++ -c $(IEIGEN) -O3 KGS.cpp

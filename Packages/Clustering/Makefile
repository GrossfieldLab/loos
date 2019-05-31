IEIGEN=-I/home/louis/anaconda3/envs/loos/include

cluster-kgs: cluster-kgs.o HAC.o AverageLinkage.o ClusteringUtils.o
	g++ -v -o cluster-kgs $(IEIGEN) -O3 cluster-kgs.o HAC.o AverageLinkage.o ClusteringUtils.o

cluster-kgs.o: cluster-kgs.cpp Clustering.hpp ClusteringUtils.hpp 
	g++ -c -O3 $(IEIGEN) cluster-kgs.cpp

#Clustering.o: ClusteringUtils.o HAC.o AverageLinkage.o KGS.o Clustering.hpp
#	g++ -c -O3 Clustering.hpp

ClusteringUtils.o: ClusteringUtils.hpp ClusteringUtils.cpp
	g++ -c -O3 $(IEIGEN) ClusteringUtils.cpp

HAC.o: HAC.hpp HAC.cpp
	g++ -c -O3 $(IEIGEN) HAC.cpp

AverageLinkage.o: HAC.o AverageLinkage.cpp AverageLinkage.hpp
	g++ -c -O3 $(IEIGEN) AverageLinkage.cpp

KGS.o: HAC.o AverageLinkage.o KGS.hpp KGS.cpp
	g++ -c $(IEIGEN) -O3 KGS.cpp

#get network
cut -f1,2 $1 > net.txt
#calculate network properties for each node separately
R --no-save --args net.txt < pRoperties.R 1> Prop.out 2>Prop.err
rm -f p1 p2 Prop.err Prop.out

#find topological modules for cluster bigger than 100 nodes
R --no-save --args net.txt fgc $2 < Modularity.R 1> Module.out 2>Module.err

clu=${1%.txt}-ClusterMemberships.txt
nei=${1%.txt}-Neighbors.txt
pro=${1%.txt}-NetworkProperties.txt

mv ClusterMemberships.txt $clu
sed 's/^/"/' NetworkProperties.txt | sed 's/\t/"\t/' > $pro
sed 's/ /\t/' Neighbors.txt | sed 's/^/"/' | sed 's/\t/"\t/' > $nei

#calculate Participation
#Two files are needed: 
#  1) network neihgborhood  
#  2) table with following columns: 
#	a) number (same order as in niegborhood file)
#	b) name
#	c) membership
#	d) in degree
#	e) out degree

#get list of clusters
sed '1d' $clu | cut -f1 | sort -n -u > clusters.txt
#iterate clusters
while read F; do
  #divide data by cluster
  grep -w -e "^$F" $clu | cut -f 2- > p1
  cut -f2 p1 | sed 's/^/"/' | sed 's/$/"/' > nodes.txt
  grep -F -f nodes.txt $nei | sed 's/"//g' > neighborhood.txt
  grep -F -f nodes.txt $pro | cut -f 2- > p2
  paste p1 p2 > properties.txt
  perl NodeParticipation.pl -n neighborhood.txt -p properties.txt
  mv participation.txt ${1%.txt}-$F-participation.txt
  mv nodes.txt $F-nodes.txt
  mv neighborhood.txt $F-neighborhood.txt
  mv properties.txt $F-properties.txt
  rm -f p1 p2 nodes.txt neighborhood.txt properties.txt
done< clusters.txt

egrep "R5|R6|R7" *participation.txt | cut -f2 > all_hubs.txt

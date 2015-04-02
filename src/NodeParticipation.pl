#!/usr/bin/perl -w
# # Copyright:: Copyright 2014 Babraham Institute
# # # Authors::   Wiktor Jurkowski <wiktor.jurkowski@uni.lu>
use strict;
use warnings;
use Getopt::Std;
my %opt;
getopts('hn:p:',\%opt);

#Two files are needed: 1) neihgborhood from R script 2) table with following columns: a) number (same order as in niegborhood file, b) name c) membership d) in degree e) out degree
die &usage() if (is_here_any($opt{n},$opt{p}));
&usage() if $opt{h};

my @sasiedzi=open_file($opt{n});
my @members=open_file($opt{p});
open(OUT, "> participation.txt") or die "Can not open an output file: $!";#output file

my (%memberships,%degin,%degout,%gsymbol,%other,%modules,%within,%party,%zety,%noderole);

foreach my $lin(@members){
  my @tab=split(/\s+/,$lin);
  my $id=shift @tab;
  $gsymbol{$id}=shift @tab;
  $memberships{$id}= shift @tab;	
  $degin{$id}= shift @tab;
  $degout{$id}= shift @tab;
  $other{$id}=join("\t",@tab);
  $modules{$memberships{$id}}++;
}

#ANALYZE LOCAL PROPERTIES OF EACH NODE
#count within node's neighberhood
foreach my $i(@sasiedzi){
  my @nb=split("\t",$i);
  my $id=$nb[1];#node ID
  my $hood=$nb[2];#use neighbors IDs
  my @tab=split(/\s+/,$hood);
  my(%sasmemb,%counts,$sumr,$sumr1,$sumr2);
  my $degree=scalar(@tab);
  my $din=$degin{$id};
  my $dout=$degout{$id};
  if($din == 0){$sumr1=1};
  if($dout == 0){$sumr2=1};
  #for eaach neighboor get membership of neighbor taken from the list of memberships of all nodes
  foreach my $k(@tab){
    my @memb=split("-",$memberships{$k});#split multiple memberships
    #counts number of neighbors of same type/module
    foreach my $j(@memb){$counts{$j}++;}
    $degree=$degree+scalar(@memb);#increase total node degree by number of multiple instances of neighbors membership
  }
  #calculate ratios for given node  
  my @tmemb=split("-",$memberships{$id});#split target node membership
  foreach my $tm(@tmemb){$within{$id}{$tm}=0;}
  foreach my $key (keys %counts){#for each modules present in the neighberhood
    #if ith node has mebership of given module
    if(exists $within{$id}{$key}){$within{$id}{$key}=$counts{$key};}#how many links it has within its own module (number of links to nodes of that module)
    my $r=($counts{$key}/$degree)**2;
    $sumr+=$r;
    if($din > 0){
      my $r1=($counts{$key}/$din)**2;
      $sumr1+=$r1;
    }
    if($dout > 0){
      my $r2=($counts{$key}/$dout)**2;
      $sumr2+=$r2;
    }
  }
  my $participation=1-$sumr;
  my $participation_in=1-$sumr1;
  my $participation_out=1-$sumr2;
  $party{$id}=$participation;
}

#calculates within-module degree z-score for each node
my (%suma,%sumak,%avwithin);
#this sections checks all modules and sum within module degree k (or square of k) over all nodes in that module
#each membership instance of given node contributes to overall within module degree
foreach my $j(keys %memberships){
  my @memb=split("-",$memberships{$j});#split node membership
  foreach my $k(@memb){
    if(exists $modules{$k}){#membership of node is of ith module
      $suma{$k}+=$within{$j}{$k};# sum all within module degree
      my $kw=($within{$j}{$k})**2;
      $sumak{$k}+=$kw;
    }
  }
}

#for multiple membership sum all within degree sums ans squares of sums as well as within module degree
#this gives averages of iver all modules in which given node is involved
foreach my $j(keys %memberships){#iterates all nodes
  my $sum=0;
  my $sumk=0;
  my $mbn=0;
  my $sw=0;
  my @memb=split("-",$memberships{$j});#split target node membership
  foreach my $k(@memb){#iterate instances of membership
    $sum=+$suma{$k};
    $sumk=+$sumak{$k};
    $mbn=+$modules{$k};
    $sw=+$within{$j}{$k};
  }
  $avwithin{$j}=$sw/scalar(@memb);
  if($sum == 0){$zety{$j}=0;}
  else{
    my $ave=$sum/$mbn;
    my $avek=$sumk/$mbn;
    my $kave=($ave)**2;
    my $bigs=sqrt($avek-$kave);
#print "cece\t$j\t$gsymbol{$j}\t$ave\t$avek\t$bigs\n";
    my $z=0;
    $z=($avwithin{$j}-$ave)/$bigs if $bigs gt 0;
    $zety{$j}=$z;
  }
}


#CLASSIFY NODES
#non-hubs	z<2.5
  #R1 ultra peripheral		P=<0.05
  #R2 peripheral		0.05<P=<0.62
  #R3 satellite connectors	0.62<P=<0.80
  #R4 kinless nodes		P>0.8
#hubs		z=>2.5
  #R5 provinical hubs		P=<0.30
  #R6 connector hubs	 	0.30<P=<0.75
  #R7 global hubs		P>0.75
printf OUT "nodeID\tGene\tParticipation\tDegree_z-score\tRole\tBetweenness\tEccentricity\tCloseness\n";
foreach my $j(keys %zety){
  my $z=$zety{$j};
  my $P=$party{$j};
  my $sy=$gsymbol{$j};
  if($z < 2.5){
    if($P <= 0.05){$noderole{$sy}="R1";}
    elsif($P > 0.05 and $P <= 0.62){$noderole{$sy}="R2";}
    elsif($P > 0.62 and $P <= 0.80){$noderole{$sy}="R3";}
    elsif($P > 0.80){$noderole{$sy}="R4";}
  }
  else{
    if($P <= 0.30){$noderole{$sy}="R5";}
    elsif($P > 0.30 and $P <= 0.75){$noderole{$sy}="R6";}
    elsif($P > 0.75){$noderole{$sy}="R7";}
  }
  printf OUT "%d\t%s\t%6.2f\t%6.2f\t%s\t%s\n",$j,$sy,$P,$z,$noderole{$sy},$other{$j};
}


#=============================
#FUNCTIONS
sub open_file{
  my ($file_name)=@_;
  open(INP1, "< $file_name") or die "Can not open an input file: $!";
  my @file1=<INP1>;
  close (INP1);
  chomp @file1;
  return @file1;
}

sub usage(){
print STDERR << "EOF";
Usage: NodeParticipation.pl -n [neighbours] -p [properties]
 -h : help message
EOF
exit;
}

sub is_here_any { ( grep $_, @_ ) < 1 }



#!/usr/bin/perl

#! MIFView is part of IsoMIF (See Below). It allows the superimposition of two proteins
#! based on MIF similarities identified with IsoMIF.

#! IsoMIF is a program to identify molecular interaction field similarities between proteins
#! Copyright (C) 2015 - Matthieu Chartier (under the supervision or Rafael Najmanovich)

#! This program is free software: you can redistribute it and/or modify
#! it under the terms of the GNU General Public License as published by
#! the Free Software Foundation, either version 3 of the License, or
#! (at your option) any later version.

#! This program is distributed in the hope that it will be useful,
#! but WITHOUT ANY WARRANTY; without even the implied warranty of
#! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#! GNU General Public License for more details.

#! You should have received a copy of the GNU General Public License
#! along with this program.  If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;

#by Matthieu Chartier
#Description
#This program generates the mifView pml files

my @probesInt=();
my @ca=();
my @pseudo=();
my @protgrid=();
my @probes=();
my @probesLab=();
my @grid=();
my @pbColors=("aquamarine","brightorange","blue","red","limegreen","lightmagenta");
my @gridLab=("200","150","100","050");
my @gridColors=("br2","palegreen","hydrogen","cyan");
my $mifViewFolder="./";
my $mifName="";
my $mifFolder="";
my $mifFile="";
my $nbpb=0;
my $grids=0;
my $gride=3;
my $ff="";
my $ot="";

#Read command line
for(my $i=0; $i<=$#ARGV; $i++){
  if($ARGV[$i] eq "-m"){ $mifFile=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-o"){ $mifViewFolder=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-t"){ $ot=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-h"){
    print "##################\nWelcome to pipeIsoMifView\n##################\n";
    print "-m         <path to mif file>\n";
    print "-o         <mifView output directory>\n";
    print "-t         <output tag name>\n";
    print "-h         <print help menu>\n";
    exit;
  }
}

if($mifFile eq ""){
  print "Error: Missing mif File\n";
  exit;
}

# $mifViewFolder=&get_dirs("/Users/matthieuchartier/hive/","mifView") if($mifViewFolder eq "");

if($mifFile=~/\/([a-z0-9_-]+)\.mif$/i){
  $mifName=$1;
  $mifFolder=$mifFile;
  $mifFolder=~s/\/$mifName\.mif$//;
}elsif($mifFile=~/^([a-z0-9_-]+)\.mif$/i){
  $mifName=$1;
  $mifFolder="./";
}else{
  print "\nError. Mif filename must contain a-z 0-9 _ or - characters.\n";
  exit;
}

my $thinness=0;
open IN, "<".$mifFile or die "cant open mif file $mifFile";
while(my $line=<IN>){
  if($line=~/^#protein_grid_distance\s+([\.0-9-]+)\s+to\s+([\.0-9-]+)/){
    $thinness=abs($2-$1);
    next;
  }
  if($line=~/^#ATOM/){
    if($line=~/#ATOM\s+([a-z0-9]+)\s+([0-9]+)\s+([a-z0-9]+)\s+([0-9]+)\s+([a-z0-9]{1})\s+([\.0-9-]+)\s+([\.0-9-]+)\s+([\.0-9-]+)\s+([0-9]{1})\s+([0-9]{1})$/i){
      if($10==1 && $3 eq "CA"){
        push @ca, "$1 $2 $3 $4 $5 $6 $7 $8";
      }
    }
    next;
  }elsif($line=~/^#PSEUDO/){
    if($line=~/#PSEUDO\s+([a-z]+)\s+([\.0-9-]+)\s+([\.0-9-]+)\s+([\.0-9-]+)$/i){
      push @pseudo, "$1 $2 $3 $4";
    }
    next;
  }elsif($line=~/^#probe\[([0-9]+)\]\s+([0-9a-z\.-]+)$/i){
    $probesLab[$1]=$2;
    $nbpb++;
    next;
  }elsif($line=~/^#zip ([0-9]+)$/i){
    $grids=$gride=$1 if($1!=-1);
    next;
  }elsif($line=~/^#ff ([0-9a-z]+)$/i){
    $ff=$1;
    next;
  }elsif($line=~/^#PG\s+([\.0-9-]+)\s+([\.0-9-]+)\s+([\.0-9-]+)$/){
      push @protgrid, "$1 $2 $3";
    next;
  }elsif($line=~/^#/){
    next;
  }elsif($line=~/^$/){
    next;
  }
  $line=~s/^\s+//g;
  $line=~s/\s+$//g;
  my @info=split(/\s+/,$line);

  #Store vrtx potential interaction
  # print $info[0]." ".$info[1]." ".$info[2]."\n";
  for(my $i=3; $i<$#info-4; $i+=3){
    my $pbid=(($i/3)-1);
    my $g0=$#info-4;
    my $g1=$#info-3;
    my $g2=$#info-2;
    my $g3=$#info-1;
    # next if($info[$i+1]>(-50));
    push @{$probes[$pbid]}, ($info[0],$info[1],$info[2],$info[$g0],$info[$g1],$info[$g2],$info[$g3],$info[$i+1],$info[$i+2]) if($info[$i]==1);
    # print $pbid." ".$probesLab[$pbid]." ".$info[$i+1]." ".$info[$g0]." ".$info[$g1]." ".$info[$g2]." ".$info[$g3]."\n";
  }

  #Store vrtx grid presence
  for(my $i=$#info-4; $i<$#info; $i++){
    my $gid=$i-($#info-4);
    my $gp=$info[$i];
    my $bu=$info[$#info];
    push @{$grid[$gid]}, ($info[0],$info[1],$info[2],$bu) if($gp==1);
    # print $gid." ".$gp." ".$bu."\n";
  }
}
close IN;

# print $ff."\n";

open NPML, ">".$mifViewFolder."/".$mifName.$ot.".pml";

#print protein
print NPML "feedback disable,all,output\n";
open IN, "<".$mifFolder."/".$mifName."_cpy.pdb";
print NPML "cmd.read_pdbstr(\"\"\"";
while(my $line=<IN>){
  chomp($line);
  print NPML $line."\\\n";
}
print NPML "TER \\\n\"\"\",\"".$mifName."\")\nset connect_mode,1\n";
close IN;

#Print pseudocenters
my $it=0;
if(@pseudo){
  print NPML "cmd.read_pdbstr(\"\"\"";
  foreach my $p (@pseudo){
    my @s=split(/\s+/,$p);
    printf NPML "HETATM%5d  N   %3s A0000    %8.3f%8.3f%8.3f  0.00 10.00           N\\\n",$it,uc($s[0]),$s[1],$s[2],$s[3];
    $it++;
  }
  print NPML "TER \\\n\"\"\",\"pseudocenters\")\nhide nonbonded\nset connect_mode,1\n";  
}

#print mif points for each probe and each grid resolution
for(my $i=0; $i<$nbpb; $i++){ #Loop each probe
  if($#{$probes[$i]}){
    for(my $g=$grids; $g<=$gride; $g++){ #Loop each grid resolution
      $probesInt[$i][$g][0]=$it;
      print NPML "cmd.read_pdbstr(\"\"\"";
      for(my $j=0; $j<@{$probes[$i]}; $j+=9){ #For each node
        if($probes[$i][$j+3+$g]==1){ #If its in this grid resolution
          printf NPML "HETATM%5d  N   %3s A0000    %8.3f%8.3f%8.3f  0.00 10.00           N\\\n",$it,$i,$probes[$i][$j],$probes[$i][$j+1],$probes[$i][$j+2];
          $it++ unless($it==99999);
        }
      }
      print NPML "TER \\\n\"\"\",\"".$probesLab[$i]."_".$gridLab[$g]."\")\n";
      $probesInt[$i][$g][1]=$it-1;
    }
  }
}
# print "grids $grids gride $gride\n";
#print grid points
for(my $i=$grids; $i<=$gride; $i++){
  my $it=0;
  if(exists $grid[$i] && scalar @{$grid[$i]} > 0){
    print NPML "cmd.read_pdbstr(\"\"\"";
    for(my $j=0; $j<@{$grid[$i]}; $j+=4){       
      printf NPML "HETATM%5d  N   %3s A0000    %8.3f%8.3f%8.3f  0.00%6.2f           N\\\n",$it,$gridLab[$i],$grid[$i][$j],$grid[$i][$j+1],$grid[$i][$j+2],$grid[$i][$j+3];
      $it++ unless($it==99999);
    }
    print NPML "TER \\\n\"\"\",\"".$gridLab[$i]."\")\n";
  }
}
print NPML "\n";
#print protein grid points
if(scalar @protgrid){
  print NPML "cmd.read_pdbstr(\"\"\"";
  for(my $i=0; $i<@protgrid; $i++){
    my @s=split(/\s+/,$protgrid[$i]);
    printf NPML "HETATM%5d  N   %3s A0000    %8.3f%8.3f%8.3f  0.00 10.00           N\\\n",0,"PGD",$s[0],$s[1],$s[2];
  }
  print NPML "TER \\\n\"\"\",\"".$mifName."_grid\")\n";
}
print NPML "\n";
print NPML "feedback enable,all,output\norient\nshow cartoon, ".$mifName."\nremove (resn HOH)\nshow sticks, HET & ".$mifName."\ncolor white,".$mifName."_grid\nshow nonbonded,".$mifName."_grid\n";
print NPML "\n";
for(my $i=0; $i<$nbpb; $i++){
  if(@{$probes[$i]}){
    for(my $g=$grids; $g<=$gride; $g++){
      if($probesInt[$i][$g][0]!=$probesInt[$i][$g][1]){
        print NPML "show spheres, ".$probesLab[$i]."_".$gridLab[$g]."\nset sphere_scale,0.2,".$probesLab[$i]."_".$gridLab[$g]."\nrebuild\n";
        print NPML "color $pbColors[$i],".$probesLab[$i]."_".$gridLab[$g]."\n";
        print NPML "hide nonbonded,".$probesLab[$i]."_".$gridLab[$g]."\n\n";
      }
    }
  }
}
print NPML "\n";
for(my $i=0; $i<3; $i++){
  if(exists $grid[$i] && scalar @{$grid[$i]} > 0){
    if($grid[$i][0]!=$grid[$i][1]){
      print NPML "color $gridColors[$i],".$gridLab[$i]."\nshow nonbonded,".$gridLab[$i]."\n";
    }
  }
}

# Print Ca atoms
# if(@ca){
#   print NPML "create ca, id ";
#   for(my $i=0; $i<@ca; $i++){
#     my @s=split(/\s+/,$ca[$i]);
#     print NPML "$s[3]";
#     print NPML "+" unless($i==$#ca);
#   }
#   print NPML " & ".$mifName."\n";
# }

# my @pbColors=("aquamarine","brightorange","blue","red","limegreen","lightmagenta");
# print NPML "color deepteal, ca\nset sphere_scale, 0.3, ca\nshow spheres, ca\n";
print NPML "\nset sphere_scale, 0.3, pseudocenter\ncolor aquamarine, resn HYD & pseudocenters\ncolor brightorange, resn ARM & pseudocenters\n";
print NPML "color blue, resn DON & pseudocenters\ncolor red, resn ACC & pseudocenters\ncolor limegreen, resn DOA & pseudocenters\nshow spheres, pseudocenters\n";

# print NPML "set sphere_scale, 0.3, 100\nshow spheres, 100\n";
# for(my $i=0; $i<15; $i++){
#   if($i<7){
#     my $c=((1-((7-$i)/7)));
#     print NPML "set_color blue".$i.", [".$c.",".$c.",1]\n";
#     print NPML "color blue".$i.", 100 & b=".$i."\n";  
#   }else{
#     my $c=((1-(($i-7)/7)));
#     print NPML "set_color red".$i.", [1,".$c.",".$c."]\n";
#     print NPML "color red".$i.", 100 & b=".$i."\n";
#   }
# }

# print PML "delete ".$mifName."\nhide lines\nshow sticks, HET & ".$mifName."_cpy\n";
# print PML "save ".$mifViewFolder."/".$mifName.".pse\n";
close NPML;

########################################
#   SUBS
########################################

sub get_dirs{
  my $rootDir=$_[0];
  my $base=$_[1];
  my $listDir=$rootDir.$base."/";
  my $dirString;

  #Get the mifs directory
  opendir my $dh, $listDir or die "$0: opendir: $!";
  my @dirs=();
  my $count=0;
  print "\n\n";
  while (defined(my $name = readdir $dh)) {
    next unless -d $listDir."$name";
    next if($name eq "..");
    print "$count $name\n";
    push @dirs, $name;
    $count++;
  }
  print "\n".$listDir;
  print "\nEnter the number of the desired dir:";
  my $answer=<STDIN>;
  if ($dirs[$answer] eq ".") {
    $dirString=$listDir;
    } else {
      $dirString=$listDir.$dirs[$answer]."/";
    }

    print "Path is: $dirString";
    return($dirString);
  }

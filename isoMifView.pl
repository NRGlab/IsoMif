#!/usr/bin/perl

#! IsoMIFView is part of IsoMIF (See Below). It allows the superimposition of two proteins
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

my @probesLab=("HYD","ARM","DON","ACC","POS","NEG");
my $matchIn="";
my $outDir="./";
my $prefix="";
my $cg=0;
my $res=1;
my $tcg=0;
my @ca=();
my @pseudo=();
my @rot=();
my @cen=();
my @va=();
my @vb=();
my @mifV1=();
my @mifV1int=();
my @mifV2=();
my @mifV2int=();
my $sm="taninorm";
my $cid="";
my $detori="";

#Read command line
for(my $i=0; $i<=$#ARGV; $i++){
  if($ARGV[$i] eq "-m"){ $matchIn=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-o"){ $outDir=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-p"){ $prefix=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-p1"){ $p1Path=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-p2"){ $p2Path=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-m1"){ $m1Path=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-m2"){ $m2Path=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-g"){ $cg=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-s"){ $sm=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-c"){ $cid=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-h"){
    print "##################\nWelcome to pipeIsoMifView\n##################\n";
    print "-m         <path to isoMif file>\n";
    print "-o         <isoMifView output directory>\n";
    print "-p         <prefix for output files>\n";
    print "-p1        <protein 1 path>\n";
    print "-p2        <protein 2 path>\n";
    print "-m1        <mif 1 path>\n";
    print "-m2        <mif 2 path>\n";
    print "-g         <coarse grain step>\n";
    print "-s         <similarity measure to find best clique>\n";
    print "-c         <clique ID to use for superimposition>\n";
    print "-h         <print help menu>\n";
    exit;
  }
}

$res=$cg;
$res=1 if($cg<0);

if($outDir eq ""){
  $outDir=&get_dirs("/Users/matthieuchartier/hive/","matchView");
}

my @probeNames=();
my @pbColors=("aquamarine","brightorange","blue","red","limegreen","lightmagenta");

my %best=();
my $cc=0;
if($sm ne ""){
  open IN, "<".$matchIn or die "Cant open match file";
  while($line=<IN>){
    if($line=~/^REMARK CLIQUE CG/){
      my @s=split(/\s+/,$line);
      next unless ($s[3] eq $cg);

        # fprintf(fpout,"REMARK CLIQUE CG %d NODES %d NODESM %d NODESMW %6.4f NORMNODES %6.4f NORMNODESRMSD %6.4f TANI %5.4f TANIM %5.4f TANIMW %5.4f TANINORM %5.4f NRG %.3f SS1 %d SS2 %d SS1M %d SS2M %d LIGRMSD %5.3f\n");
        $data{$k1}{$k2}{"nodes"}[$c]=$s[5];
        $data{$k1}{$k2}{"nodesm"}[$c]=$s[7];
        $data{$k1}{$k2}{"nodesmw"}[$c]=$s[9];
        $data{$k1}{$k2}{"normnodes"}[$c]=$s[11];
        $data{$k1}{$k2}{"normnodesrmsd"}[$c]=$s[13];
        $data{$k1}{$k2}{"tani"}[$c]=$s[15];
        $data{$k1}{$k2}{"tanim"}[$c]=$s[17];
        $data{$k1}{$k2}{"tanimw"}[$c]=$s[19];
        $data{$k1}{$k2}{"taninorm"}[$c]=$s[21];
        $data{$k1}{$k2}{"nrg"}[$c]=$s[23];
        $data{$k1}{$k2}{"ss1"}[$c]=$s[25];
        $data{$k1}{$k2}{"ss2"}[$c]=$s[27];
        $data{$k1}{$k2}{"ss1m"}[$c]=$s[29];
        $data{$k1}{$k2}{"ss2m"}[$c]=$s[31];
        $data{$k1}{$k2}{"ligrmsd"}[$c]=$s[33];
      if($cc==0){
        $best{"nodes"}[0]=$cc;
        $best{"nodes"}[1]=$s[5];
        $best{"nodesm"}[0]=$cc;
        $best{"nodesm"}[1]=$s[7];
        $best{"normnodes"}[0]=$cc;
        $best{"normnodes"}[1]=$s[11];
        $best{"normnodesrmsd"}[0]=$cc;
        $best{"normnodesrmsd"}[1]=$s[13];
        $best{"tani"}[0]=$cc;
        $best{"tani"}[1]=$s[15];
        $best{"tanim"}[0]=$cc;
        $best{"tanim"}[1]=$s[17];
        $best{"taninorm"}[0]=$cc;
        $best{"taninorm"}[1]=$s[21];
        $best{"nrg"}[0]=$cc;
        $best{"nrg"}[1]=$s[23];
        $best{"ligrmsd"}[0]=$cc;
        $best{"ligrmsd"}[1]=$s[33];
      }else{
        if($s[5]>$best{"nodes"}[1]){
          $best{"nodes"}[1]=$s[5];
          $best{"nodes"}[0]=$cc;
        }
        if($s[7]>$best{"nodesm"}[1]){
          $best{"nodesm"}[1]=$s[7];
          $best{"nodesm"}[0]=$cc;
        }
        if($s[11]>$best{"normnodes"}[1]){
          $best{"normnodes"}[1]=$s[11];
          $best{"normnodes"}[0]=$cc;
        }
        if($s[13]>$best{"normnodesrmsd"}[1]){
          $best{"normnodesrmsd"}[1]=$s[13];
          $best{"normnodesrmsd"}[0]=$cc;
        }
        if($s[15]>$best{"tani"}[1]){
          $best{"tani"}[1]=$s[15];
          $best{"tani"}[0]=$cc;
        }
        if($s[17]>$best{"tanim"}[1]){
          $best{"tanim"}[1]=$s[17];
          $best{"tanim"}[0]=$cc;
        }
        if($s[21]>$best{"taninorm"}[1]){
          $best{"taninorm"}[1]=$s[21];
          $best{"taninorm"}[0]=$cc;
        }
        if($s[23]>$best{"nrg"}[1]){
          $best{"nrg"}[1]=$s[23];
          $best{"nrg"}[0]=$cc;
        }
        if($s[33]<$best{"ligrmsd"}[1]){
          $best{"ligrmsd"}[1]=$s[33];
          $best{"ligrmsd"}[0]=$cc;
        }
      }
      $cc++;
    }
  }
  close IN;
  print "Best ".$sm." ".$best{$sm}[0]." ".$best{$sm}[1]."\n";
}

$cc=0;
my $flagStore=0;
#Retrieve the nodes and other info from match file
open IN, "<".$matchIn or die "Cant open match file";
while($line=<IN>){
  
  chomp($_);
  if($line=~/^REMARK\s+mif_file_1:\s+([-_\.\/a-z0-9]+)/i){
    $mifFilePath1=$1;
    $m1Path=$mifFilePath1 unless($m1Path ne "");
    if($mifFilePath1=~/\/([-_\.a-z0-9]+)$/i){
      $mif1=$1;
    }elsif($mifFilePath1=~/^([-_\.a-z0-9]+)$/i){
      $mif1=$1;
    }
    unless($p1Path ne ""){
      $p1Path=$mifFilePath1;
      $p1Path=~s/\.mif/_cpy\.pdb/;
    }
  }

  if($line=~/^REMARK\s+mif_file_2:\s+([-_\.\/a-z0-9]+)/i){
    $mifFilePath2=$1;
    $m2Path=$mifFilePath2 unless($m2Path ne "");
    if($mifFilePath2=~/\/([-_\.a-z0-9]+)$/i){
      $mif2=$1;
    }elsif($mifFilePath2=~/^([-_\.a-z0-9]+)$/i){
      $mif2=$1;
    }
    unless($p2Path ne ""){
      $p2Path=$mifFilePath2;
      $p2Path=~s/\.mif/_cpy\.pdb/;
    }
  }

  if($line=~/^REMARK CLIQUE CG ([0-9-]+)/){
    $tcg=$1;
    last if($flagStore==1);
    $flagStore=1 if(($sm ne "" && $best{$sm}[0]==$cc && $cid eq "") || ($cid ne "" && $cc==$cid));
    $cc++ if($tcg==$cg);
  }

  if($line=~/^REMARK ROTMAT\s+([-\.\/0-9\s+]+)$/i && $tcg==$cg){
    my @rd=split(/\s+/,$1);
    $rot[0][0]=$rd[0]; $rot[0][1]=$rd[1]; $rot[0][2]=$rd[2];
    $rot[1][0]=$rd[3]; $rot[1][1]=$rd[4]; $rot[1][2]=$rd[5];
    $rot[2][0]=$rd[6]; $rot[2][1]=$rd[7]; $rot[2][2]=$rd[8];
  }

  if($line=~/^REMARK DETORI\s+([-0-9]+)$/i && $tcg==$cg){
    $detori=$1;
  }

  if($line=~/^REMARK CENTRES\s+([-\.\/0-9\s+]+)$/i && $tcg==$cg){
    my @rd=split(/\s+/,$1);
    $cen[0][0]=$rd[0]; $cen[0][1]=$rd[1]; $cen[0][2]=$rd[2];
    $cen[1][0]=$rd[3]; $cen[1][1]=$rd[4]; $cen[1][2]=$rd[5];
  }

  if($line!~/^REMARK/ && $tcg==$cg && $tcg==-1){
    my @l=split(/\s+/,$line);
    push @ca, "$l[0];$l[1];$l[4];$l[8];$l[9];$l[12]";
  }elsif($line=~/^([a-z]{3})\s+([0-9\.-]+)\s+([0-9\.-]+)\s+([0-9\.-]+)\s+([a-z]{3})\s+([0-9\.-]+)\s+([0-9\.-]+)\s+([0-9\.-]+)/ && $cg==-3){
    push @pseudo, "$1;$2;$3;$4;$5;$6;$7;$8";
  }elsif($line=~/^A\s+([0-9\.-]+)\s+([0-9\.-]+)\s+([0-9\.-]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)/ && $cg==-2){
    $line=~s/^A\s+//;
    my @s=split(/\s+/,$line);
    for(my $p=0; $p<6; $p++){
      push @{$va[$p]}, ($s[0],$s[1],$s[2]) if($s[$p+3]==1);
    }
  }elsif($line=~/^B\s+([0-9\.-]+)\s+([0-9\.-]+)\s+([0-9\.-]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)/ && $cg==-2){
    $line=~s/^B\s+//;
    my @s=split(/\s+/,$line);
    for(my $p=0; $p<6; $p++){
      push @{$vb[$p]}, ($s[0],$s[1],$s[2]) if($s[$p+3]==1);
    }
  }elsif($line!~/^REMARK/ && $tcg==$cg){
    next unless($flagStore==1);
    my @l=split(/\s+/,$line);
    $pb=$l[0];
    $data[$tcg][$pb].="$l[1];$l[2];$l[3];$l[4];$l[5];$l[6]\n";
  }
}
close IN;

for(my $i=0; $i<3; $i++){
  for(my $j=0; $j<3; $j++){
    print $rot[$i][$j]." ";
  }
}
print "\n";
print "$cen[0][0] $cen[0][1] $cen[0][2]\n";
print "$cen[1][0] $cen[1][1] $cen[1][2]\n";
print "detori: ".$detori."\n";

&storeMif($m1Path,\@mifV1) if(-e $m1Path);
&storeMif($m2Path,\@mifV2) if(-e $m2Path);

sub storeMif{
  open IN, "<".$_[0] or die "cant open mif file $mifFile";
  # print $_[0]."\n";
  while(my $line=<IN>){
    if($line=~/^#probe\[([0-9]+)\]\s+([0-9a-z\.-]+)$/i){
      $probeNames[$1]=$2;
      next;
    }
    next if($line=~/^#ATOM/);
    next if($line=~/^#/);
    next if($line=~/^$/);
    $line=~s/^\s+//g;
    $line=~s/\s+$//g;
    my @info=split(/\s+/,$line);
    #Store vrtx potential interaction
    for(my $i=3; $i<$#info-4; $i+=3){
      my $pbid=($i/3)-1;
      my $g0=$#info-4;
      my $g1=$#info-3;
      my $g2=$#info-2;
      my $g3=$#info-1;
      # print "$pbid $info[0] $info[1] $info[2]\n" if($info[$i]==1);
      push @{${$_[1]}[($i/3)-1]}, ($info[0],$info[1],$info[2],$info[$g0],$info[$g1],$info[$g2],$info[$g3]) if($info[$i]==1);
    }
    # #Store vrtx grid presence
    # for(my $i=9; $i<13; $i++){
    #     push @{$grid[$i-9]}, ($info[0],$info[1],$info[2]) if($info[$i]==1);
    # }
  }
  close IN;
}

$mif1=~s/\.mif//;
$mif2=~s/\.mif//;

my $tag=$mif1."_".$mif2;
$tag=$prefix."_".$tag if($prefix ne "");

my $mif1str=&printMif(1,\@mifV1,\@mifV1int);
my $mif2str=&printMif(2,\@mifV2,\@mifV2int);

sub printMif{
  # open OUT, ">".$outDir.$tag."_mif".$_[0].".pdb" or die "Cant open mif1 out file";
  my $it=0;
  # print "\n\nprinting mif $_[0]\n";
  my $pdbstr="cmd.read_pdbstr(\"\"\"";
  for(my $i=0; $i<@{$_[1]}; $i++){ #Loop each probe
    if(scalar @{$_[1]}){
      ${$_[2]}[$i][$res][0]=$it;
      for(my $j=0; $j<@{${$_[1]}[$i]}; $j+=7){ #For each node
        if(${$_[1]}[$i][$j+3+$res]==1){ #If its in the g grid resolution
          my $coor=();
          my $ncoor=();
          $coor[0]=$ncoor[0]=${$_[1]}[$i][$j];
          $coor[1]=$ncoor[1]=${$_[1]}[$i][$j+1];
          $coor[2]=$ncoor[2]=${$_[1]}[$i][$j+2];
          if($_[0]==1){
            for(my $i=0; $i<3; $i++){
              $ncoor[$i]=$cen[1][$i];
              for(my $j=0; $j<3; $j++){
                $ncoor[$i]+=($coor[$j]-$cen[0][$j])*$rot[$i][$j];
              }
            } 
          }
          $pdbstr.=sprintf("HETATM%5d  N   %3s A0000    %8.3f%8.3f%8.3f  0.00 10.00           N\\\n",$it,$i,$ncoor[0],$ncoor[1],$ncoor[2]);
          # printf OUT "HETATM%5d  N   %3s A0000    %8.3f%8.3f%8.3f  0.00 10.00           N\n",$it,$probeNames[$i],$ncoor[0],$ncoor[1],$ncoor[2];          
          $it++;
        }
      }
      ${$_[2]}[$i][$res][1]=$it-1;
      # print "probe $i ${$_[2]}[$i][$res][0] ${$_[2]}[$i][$res][1]\n";
    }
  }
  $pdbstr.="TER \\\n\"\"\",\"".$tag."_mif".$_[0]."\")\n";
  # close OUT;
  return($pdbstr);
}

#Create protein file 1
my $p1str="cmd.read_pdbstr(\"\"\"";
# open PDB1OUT, ">".$outDir.$tag."_1.pdb" or die "Cant open ".$outDir.$tag."_1.pdb";
open PDB1IN, "<".$p1Path or die "Cant open ".$p1Path;
while (my $line = <PDB1IN>) {
  if($line=~/^ATOM/ or $line=~/^HETATM/){
    my @coor=();
    my @ncoor=();
    my $b4=substr($line,0,30);
    my $after=substr($line,54);
    chomp($after);
    $coor[0]=substr($line,30,8);
    $coor[1]=substr($line,38,8);
    $coor[2]=substr($line,46,8);
    $coor[0]=~s/\s+//g;
    $coor[1]=~s/\s+//g;
    $coor[2]=~s/\s+//g;
    # print $line;
    for(my $i=0; $i<3; $i++){
      $ncoor[$i]=$cen[1][$i];
      # printf("\n%10.5f",$ncoor[$i]);
      for(my $j=0; $j<3; $j++){
        $ncoor[$i]+=($coor[$j]-$cen[0][$j])*$rot[$i][$j];
      }
      # printf("\n%10.5f\n",$ncoor[$i]);
    }
    # printf PDB1OUT $b4."%8.3f%8.3f%8.3f".$after."\n",$ncoor[0],$ncoor[1],$ncoor[2];
    $p1str.=sprintf($b4."%8.3f%8.3f%8.3f".$after."\\\n",$ncoor[0],$ncoor[1],$ncoor[2]);
    # printf $b4."%8.3f%8.3f%8.3f".$after."\n",$ncoor[0],$ncoor[1],$ncoor[2];
  }else{
    # print PDB1OUT $line;
  }
}
$p1str.="TER \\\n\"\"\",\"".$mif1."\")\n";;
close PDB1IN;
# close PDB1OUT;

open IN, "<".$p2Path;
my $p2str="cmd.read_pdbstr(\"\"\"";
while(my $line=<IN>){
  chomp($line);
  $p2str.=$line."\\\n";
}
close IN;
$p2str.="TER \\\n\"\"\",\"".$mif2."\")\n";

open NPML, ">".$outDir.$tag.".pml";
print NPML $p1str.$p2str."show cartoon\nhide lines\nset connect_mode,1\n".$mif1str.$mif2str;

# system("cp ".$p2Path." ".$outDir.$tag."_2.pdb");
# open PML3, ">".$outDir.$tag.".pml" or die "Cant open ".$tag.".pml";

# print PML3 "load ".$outDir.$tag."_1.pdb, $mif1\nload ".$outDir.$tag."_2.pdb, $mif2\nremove (resn HOH)\nshow cartoon\nhide lines\n";
if($cg==-1){
  foreach my $nod (@ca){
    my @s=split(/;/,$nod);
    # print PML3 "show lines, resi $s[1] & chain $s[2] & ".$mif1."\nshow lines, resi $s[4] & chain $s[5] & ".$mif2."\n";
    print NPML "show lines, resi $s[1] & chain $s[2] & ".$mif1."\nshow lines, resi $s[4] & chain $s[5] & ".$mif2."\n";
  }
}elsif($cg==-3){
  my $ps1="";
  my $ps2="";
  foreach my $p (@pseudo){
    my @s=split(/;/,$p);
    my $id=0;
    my @coor=();
    my @ncoor=();
    $coor[0]=$s[1];
    $coor[1]=$s[2];
    $coor[2]=$s[3];
    for(my $i=0; $i<3; $i++){
      $ncoor[$i]=$cen[1][$i];
      for(my $k=0; $k<3; $k++){
        $ncoor[$i]+=($coor[$k]-$cen[0][$k])*$rot[$i][$k];
      }
    }
    $ps1.=sprintf("HETATM%5d  CA  %3s A        %8.3f%8.3f%8.3f  0.00 10.00           C  \\\n",$id,uc($s[0]),$ncoor[0],$ncoor[1],$ncoor[2]);
    $ps2.=sprintf("HETATM%5d  CA  %3s A        %8.3f%8.3f%8.3f  0.00 10.00           C  \\\n",$id,uc($s[4]),$s[5],$s[6],$s[7]);
    $id++;
  }
  print NPML "cmd.read_pdbstr(\"\"\"".$ps1."TER \\\n\"\"\",\"".$tag."_1_pseudo\")\n";
  print NPML "cmd.read_pdbstr(\"\"\"".$ps2."TER \\\n\"\"\",\"".$tag."_2_pseudo\")\n";
}elsif($cg==-2){
  # print PML3 "set connect_mode,1\nload ".$outDir.$tag."_1_nodes.pdb\nset connect_mode,1\nload ".$outDir.$tag."_2_nodes.pdb\n";
  
  my $id=0;
  my $strnodes="";
  # open NODES1, ">".$outDir.$tag."_1_nodes.pdb" or die "Cant open ".$outDir.$tag."_1_nodes.pdb";
  print NPML "set connect_mode,1\ncmd.read_pdbstr(\"\"\"";
  for(my $p=0; $p<6; $p++){
    if(scalar @{$va[$p]}>0){
      my $start=$id;
      for(my $i=0; $i<@{$va[$p]}; $i+=3){
        # printf NODES1 "HETATM%5d  CA  NRG A        %8.3f%8.3f%8.3f  0.00 10.00           C  \n",$id,$va[$p][$i],$va[$p][$i+1],$va[$p][$i+2];
        printf NPML "HETATM%5d  CA  NRG A        %8.3f%8.3f%8.3f  0.00 10.00           C  \\\n",$id,$va[$p][$i],$va[$p][$i+1],$va[$p][$i+2];
        $id++;
        # my $sd=1000;
        # # print "\n$p - $va[$p][$i] $va[$p][$i+1] $va[$p][$i+2]";
        # for(my $j=0; $j<@{$vb[$p]}; $j+=3){
        #   my $dist=dist3D($va[$p][$i],$va[$p][$i+1],$va[$p][$i+2],$vb[$p][$j],$vb[$p][$j+1],$vb[$p][$j+2]);
        #   $sd=$dist if($dist<$sd);
        #   if($dist<=2.01){
        #     # print " -> $vb[$p][$j] $vb[$p][$j+1] $vb[$p][$j+2] = $dist";
        #     last;  
        #   }
        # }
        # print " (sd: $sd)";
      }
      my $stop=$id-1;
      # print PML3 "create ".$probeNames[$p]."_".$mif1.", id $start-$stop & ".$tag."_1_nodes\nset sphere_scale,0.25,".$probeNames[$p]."_".$mif1."\nshow spheres, ".$probeNames[$p]."_".$mif1."\nrebuild\ncolor ".$pbColors[$p].", ".$probeNames[$p]."_".$mif1."\n";
      $strnodes.="create ".$probeNames[$p]."_".$mif1.", id $start-$stop & ".$tag."_1_nodes\nset sphere_scale,0.25,".$probeNames[$p]."_".$mif1."\nshow spheres, ".$probeNames[$p]."_".$mif1."\nrebuild\ncolor ".$pbColors[$p].", ".$probeNames[$p]."_".$mif1."\n";
    }
  }
  # close NODES1;
  print NPML "TER \\\n\"\"\",\"".$tag."_1_nodes\")\n";

  $id=0;
  print NPML "set connect_mode,1\ncmd.read_pdbstr(\"\"\"";
  # open NODES2, ">".$outDir.$tag."_2_nodes.pdb" or die "Cant open ".$outDir.$tag."_2_nodes.pdb";
  for(my $p=0; $p<6; $p++){
    if(scalar @{$vb[$p]}>0){
      my $start=$id;
      for(my $i=0; $i<@{$vb[$p]}; $i+=3){
        # printf NODES2 "HETATM%5d  CA  NRG A        %8.3f%8.3f%8.3f  0.00 10.00           C  \n",$id,$vb[$p][$i],$vb[$p][$i+1],$vb[$p][$i+2];
        printf NPML "HETATM%5d  CA  NRG A        %8.3f%8.3f%8.3f  0.00 10.00           C  \\\n",$id,$vb[$p][$i],$vb[$p][$i+1],$vb[$p][$i+2];
        
        $id++;
      }
      my $stop=$id-1;
      # print PML3 "create ".$probeNames[$p]."_".$mif2.", id $start-$stop & ".$tag."_2_nodes\nset sphere_scale,0.15,".$probeNames[$p]."_".$mif2."\nshow spheres, ".$probeNames[$p]."_".$mif2."\nrebuild\ncolor ".$pbColors[$p].", ".$probeNames[$p]."_".$mif2."\n";
      $strnodes.="create ".$probeNames[$p]."_".$mif2.", id $start-$stop & ".$tag."_2_nodes\nset sphere_scale,0.15,".$probeNames[$p]."_".$mif2."\nshow spheres, ".$probeNames[$p]."_".$mif2."\nrebuild\ncolor ".$pbColors[$p].", ".$probeNames[$p]."_".$mif2."\n";
    }
  }
  # close NODES2;
  print NPML "TER \\\n\"\"\",\"".$tag."_2_nodes\")\n";
  print NPML $strnodes;

}else{

  # print PML3 "set connect_mode,1\nload ".$outDir.$tag."_1_nodes.pdb\nset connect_mode,1\nload ".$outDir.$tag."_2_nodes.pdb\n";

  my $ids=0;
  # open NODES1, ">".$outDir.$tag."_1_nodes.pdb" or die "Cant open ".$outDir.$tag."_1_nodes.pdb";
  # open NODES2, ">".$outDir.$tag."_2_nodes.pdb" or die "Cant open ".$outDir.$tag."_2_nodes.pdb";
  my $str1="cmd.read_pdbstr(\"\"\"";
  my $str2="cmd.read_pdbstr(\"\"\"";
  my $strSel="";

  for(my$j=0; $j<@{$data[$cg]}; $j++){ #For each probe
    my @nodes=split(/\n/,$data[$cg][$j]);
    if(@nodes){
      my $start=$ids;
      foreach $node (@nodes){
        my @info=split(/;/,$node);
        my @coor=();
        my @ncoor=();
        $coor[0]=$info[0];
        $coor[1]=$info[1];
        $coor[2]=$info[2];
        # printf("%8.3f %8.3f %8.3f\n",$coor[0],$coor[1],$coor[2]);
        for(my $i=0; $i<3; $i++){
          $ncoor[$i]=$cen[1][$i];
          for(my $k=0; $k<3; $k++){
            $ncoor[$i]+=($coor[$k]-$cen[0][$k])*$rot[$i][$k];
          }
        }
        # printf("%8.3f %8.3f %8.3f\n",$ncoor[0],$ncoor[1],$ncoor[2]);
        # printf NODES1 "HETATM%5d  CA  NRG A        %8.3f%8.3f%8.3f  0.00 10.00           C  \n",$ids,$ncoor[0],$ncoor[1],$ncoor[2];
        # printf NODES2 "HETATM%5d  CA  NRG A        %8.3f%8.3f%8.3f  0.00 10.00           C  \n",$ids,$info[3],$info[4],$info[5];
        $str1.=sprintf("HETATM%5d  CA  NRG A        %8.3f%8.3f%8.3f  0.00 10.00           C  \\\n",$ids,$ncoor[0],$ncoor[1],$ncoor[2]);
        $str2.=sprintf("HETATM%5d  CA  NRG A        %8.3f%8.3f%8.3f  0.00 10.00           C  \\\n",$ids,$info[3],$info[4],$info[5]);

        $ids++;
      }
      my $stop=$ids-1;
      # print PML3 "create ".$probeNames[$j]."_".$mif1.", id $start-$stop & ".$tag."_1_nodes\nset sphere_scale,0.25,".$probeNames[$j]."_".$mif1."\nshow spheres, ".$probeNames[$j]."_".$mif1."\nrebuild\ncolor ".$pbColors[$j].", ".$probeNames[$j]."_".$mif1."\n";
      # print PML3 "create ".$probeNames[$j]."_".$mif2.", id $start-$stop & ".$tag."_2_nodes\nset sphere_scale,0.15,".$probeNames[$j]."_".$mif2."\nshow spheres, ".$probeNames[$j]."_".$mif2."\nrebuild\ncolor ".$pbColors[$j].", ".$probeNames[$j]."_".$mif2."\n";
      $strSel.="create ".$probeNames[$j]."_".$mif1.", id $start-$stop & ".$tag."_1_nodes\nset sphere_scale,0.25,".$probeNames[$j]."_".$mif1."\nshow spheres, ".$probeNames[$j]."_".$mif1."\nrebuild\ncolor ".$pbColors[$j].", ".$probeNames[$j]."_".$mif1."\n";
      $strSel.="create ".$probeNames[$j]."_".$mif2.", id $start-$stop & ".$tag."_2_nodes\nset sphere_scale,0.15,".$probeNames[$j]."_".$mif2."\nshow spheres, ".$probeNames[$j]."_".$mif2."\nrebuild\ncolor ".$pbColors[$j].", ".$probeNames[$j]."_".$mif2."\n";
    }
  }
  $str1.="TER \\\n\"\"\",\"".$tag."_1_nodes\")\n";
  $str2.="TER \\\n\"\"\",\"".$tag."_2_nodes\")\n";
  print NPML $str1.$str2.$strSel;
}

my $mstr1=&printMifPml(\@mifV1,\@mifV1int,1,$mif1,0.25);
my $mstr2=&printMifPml(\@mifV2,\@mifV2int,2,$mif2,0.15);

print NPML "set connect_mode,1\n".$mstr1.$mstr2;

sub printMifPml{
  my $str="";
  # print PML3 "set connect_mode,1\nload ".$outDir.$tag."_mif".$_[2].".pdb\n";
  for(my $i=0; $i<@{$_[0]}; $i++){
    if(@{${$_[0]}[$i]}){
        if(${$_[1]}[$i][$res][0]!=${$_[1]}[$i][$res][1]){
          $str.="create mif_".$_[3]."_".$probeNames[$i].", id ${$_[1]}[$i][$res][0]-${$_[1]}[$i][$res][1] & ".$tag."_mif".$_[2]."\n";
          $str.="show spheres, mif_".$_[3]."_".$probeNames[$i]."\nset sphere_scale,".$_[4].",mif_".$_[3]."_".$probeNames[$i]."\nset sphere_transparency,0.6,mif_".$_[3]."_".$probeNames[$i]."\nrebuild\n";
          $str.="color $pbColors[$i],mif_".$_[3]."_".$probeNames[$i]."\n";
          $str.="hide nonbonded,mif_".$_[3]."_".$probeNames[$i]."\n";
        }
    }
  }
  $str.="delete ".$tag."_mif".$_[2]."\n";
  # print PML3 "delete ".$tag."_mif".$_[2]."\n";
  return($str);
}
print NPML "remove hydrogens\nshow sticks, HET\ndelete ".$tag."_1_nodes\ndelete ".$tag."_2_nodes\n";
# print NPML "hide nonbonded\n";

my @pseudoLab=("HYD","ARM","DON","ACC","DOA");

print NPML "set sphere_scale, 0.25, ".$tag."_1_pseudo\ncolor aquamarine, resn HYD & ".$tag."_1_pseudo\ncolor brightorange, resn ARM & ".$tag."_1_pseudo\n";
print NPML "color blue, resn DON & ".$tag."_1_pseudo\ncolor red, resn ACC & ".$tag."_1_pseudo\ncolor limegreen, resn DOA & ".$tag."_1_pseudo\nshow spheres, ".$tag."_1_pseudo\n";
print NPML "set sphere_scale, 0.15, ".$tag."_2_pseudo\ncolor aquamarine, resn HYD & ".$tag."_2_pseudo\ncolor brightorange, resn ARM & ".$tag."_2_pseudo\n";
print NPML "color blue, resn DON & ".$tag."_2_pseudo\ncolor red, resn ACC & ".$tag."_2_pseudo\ncolor limegreen, resn DOA & ".$tag."_2_pseudo\nshow spheres, ".$tag."_2_pseudo\n";

# print PML3 "remove hydrogens\nhide nonbonded\nshow sticks, HET\ndelete ".$tag."_1_nodes\ndelete ".$tag."_2_nodes\nsave ".$outDir.$tag.".pse";
# close NODES1;
# close NODES2;
# close PML3;
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
  $answer=<STDIN>;
  if ($dirs[$answer] eq ".") {
    $dirString=$listDir;
  } else {
    $dirString=$listDir.$dirs[$answer]."/";
  }
  
  print "Path is: $dirString";
  return($dirString);
}

sub dist3D{
  my $dist=sqrt( ($_[0]-$_[3]) * ($_[0]-$_[3]) + ($_[1]-$_[4]) * ($_[1]-$_[4]) + ($_[2]-$_[5]) * ($_[2]-$_[5]) );
  return($dist);
}

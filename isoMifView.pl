#!/usr/bin/perl

#by Matthieu Chartier
#Description
#This program takes a match file and creates a .pml file to visualize results

my @probesLab=("HYD","ARM","DON","ACC","POS","NEG");
my $matchIn="";
my $outDir="./";
my $prefix="";
my $cg=0;
my $res=1;
my $tcg=0;
my @ca=();
my @rot=();
my @cen=();
my @va=();
my @vb=();
my @mifV1=();
my @mifV1int=();
my @mifV2=();
my @mifV2int=();

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
    print "-h         <print help menu>\n";
    exit;
  }
}

$res=$cg;
$res=1 if($cg<0);

if($outDir eq ""){
  $outDir=&get_dirs("/Users/matthieuchartier/hive/","matchView");
}

my @probeNames=("HYD","ARM","DON","ACC","POS","NEG");
my @pbColors=("aquamarine","brightorange","blue","red","limegreen","lightmagenta");

#Retrieve the nodes and other info from match file
open IN, "<".$matchIn or die "Cant open match file";
while($line=<IN>){
  
  chomp($_);
  if($line=~/^REMARK\s+mif_file_1:\s+([-_\.\/a-z0-9]+)/i){
    $mifFilePath1=$1;
    $m1Path=$mifFilePath1 unless($m1Path ne "");
    if($mifFilePath1=~/\/([_\.a-z0-9]+)$/i){
      $mif1=$1;
    }elsif($mifFilePath1=~/^([_\.a-z0-9]+)$/i){
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
    if($mifFilePath2=~/\/([_\.a-z0-9]+)$/i){
      $mif2=$1;
    }elsif($mifFilePath2=~/^([_\.a-z0-9]+)$/i){
      $mif2=$1;
    }
    unless($p2Path ne ""){
      $p2Path=$mifFilePath2;
      $p2Path=~s/\.mif/_cpy\.pdb/;
    }
  }

  if($line=~/^REMARK CLIQUE CG ([0-9-]+) NODES ([0-9]+) TANI ([0-9\.]+) SS1 ([0-9]+) SS2 ([0-9]+)$/){
    $tcg=$1;
  }

  if($line=~/^REMARK ROTMAT\s+([-\.\/0-9\s+]+)$/i && $tcg==$cg){
    my @rd=split(/\s+/,$1);
    $rot[0][0]=$rd[0]; $rot[0][1]=$rd[1]; $rot[0][2]=$rd[2];
    $rot[1][0]=$rd[3]; $rot[1][1]=$rd[4]; $rot[1][2]=$rd[5];
    $rot[2][0]=$rd[6]; $rot[2][1]=$rd[7]; $rot[2][2]=$rd[8];
  }

  if($line=~/^REMARK CENTRES\s+([-\.\/0-9\s+]+)$/i && $tcg==$cg){
    my @rd=split(/\s+/,$1);
    $cen[0][0]=$rd[0]; $cen[0][1]=$rd[1]; $cen[0][2]=$rd[2];
    $cen[1][0]=$rd[3]; $cen[1][1]=$rd[4]; $cen[1][2]=$rd[5];
  }

  if($line!~/^REMARK/ && $tcg==$cg && $tcg==-1){
    my @l=split(/\s+/,$line);
    push @ca, "$l[0];$l[1];$l[4];$l[8];$l[9];$l[12]";
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
    my @l=split(/\s+/,$line);
    $pb=$l[0];
    $data[$tcg][$pb].="$l[1];$l[2];$l[3];$l[4];$l[5];$l[6]\n";
  }
}
close IN;

&storeMif($m1Path,\@mifV1) if(-e $m1Path);
&storeMif($m2Path,\@mifV2) if(-e $m2Path);

sub storeMif{
  open IN, "<".$_[0] or die "cant open mif file $mifFile";
  while(my $line=<IN>){
    next if($line=~/^#ATOM/);
    next if($line=~/^#/);
    next if($line=~/^$/);
    $line=~s/^\s+//g;
    $line=~s/\s+$//g;
    my @info=split(/\s+/,$line);
    #Store vrtx potential interaction
    for(my $i=3; $i<9; $i++){
      push @{${$_[1]}[$i-3]}, ($info[0],$info[1],$info[2],$info[9],$info[10],$info[11],$info[12]) if($info[$i]==1);
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
  my $pdbstr="cmd.read_pdbstr(\"\"\"";
  for(my $i=0; $i<6; $i++){ #Loop each probe
    # print $i." $#{${$_[1]}[$i]}\n";
    if($#{${$_[1]}[$i]}){
      ${$_[2]}[$i][$res][0]=$it;
      for(my $j=0; $j<@{${$_[1]}[$i]}; $j+=7){ #For each node
        if(${$_[1]}[$i][$j+3+$res]==1){ #If its in this grid resolution
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
          $pdbstr.=sprintf("HETATM%5d  N   %3s A0000    %8.3f%8.3f%8.3f  0.00 10.00           N\\\n",$it,$probesLab[$i],$ncoor[0],$ncoor[1],$ncoor[2]);
          # printf OUT "HETATM%5d  N   %3s A0000    %8.3f%8.3f%8.3f  0.00 10.00           N\n",$it,$probesLab[$i],$ncoor[0],$ncoor[1],$ncoor[2];          
          $it++;
        }
      }
      ${$_[2]}[$i][$res][1]=$it-1;
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
print NPML $p1str.$p2str."remove (resn HOH)\nshow cartoon\nhide lines\nset connect_mode,1\n".$mif1str.$mif2str;

# system("cp ".$p2Path." ".$outDir.$tag."_2.pdb");
# open PML3, ">".$outDir.$tag.".pml" or die "Cant open ".$tag.".pml";

# print PML3 "load ".$outDir.$tag."_1.pdb, $mif1\nload ".$outDir.$tag."_2.pdb, $mif2\nremove (resn HOH)\nshow cartoon\nhide lines\n";
if($cg==-1){
  foreach my $nod (@ca){
    my @s=split(/;/,$nod);
    # print PML3 "show lines, resi $s[1] & chain $s[2] & ".$mif1."\nshow lines, resi $s[4] & chain $s[5] & ".$mif2."\n";
    print NPML "show lines, resi $s[1] & chain $s[2] & ".$mif1."\nshow lines, resi $s[4] & chain $s[5] & ".$mif2."\n";
  }
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
  for(my $i=0; $i<6; $i++){
    if(@{${$_[0]}[$i]}){
        if(${$_[1]}[$i][$res][0]!=${$_[1]}[$i][$res][1]){
          $str.="create mif_".$_[3]."_".$probesLab[$i].", id ${$_[1]}[$i][$res][0]-${$_[1]}[$i][$res][1] & ".$tag."_mif".$_[2]."\n";
          $str.="show spheres, mif_".$_[3]."_".$probesLab[$i]."\nset sphere_scale,".$_[4].",mif_".$_[3]."_".$probesLab[$i]."\nset sphere_transparency,0.6,mif_".$_[3]."_".$probesLab[$i]."\nrebuild\n";
          $str.="color $pbColors[$i],mif_".$_[3]."_".$probesLab[$i]."\nhide nonbonded,mif_".$_[3]."_".$probesLab[$i]."\n";
          # print PML3 "create mif_".$_[3]."_".$probesLab[$i].", id ${$_[1]}[$i][$res][0]-${$_[1]}[$i][$res][1] & ".$tag."_mif".$_[2]."\n";
          # print PML3 "show spheres, mif_".$_[3]."_".$probesLab[$i]."\nset sphere_scale,".$_[4].",mif_".$_[3]."_".$probesLab[$i]."\nset sphere_transparency,0.6,mif_".$_[3]."_".$probesLab[$i]."\nrebuild\n";
          # print PML3 "color $pbColors[$i],mif_".$_[3]."_".$probesLab[$i]."\nhide nonbonded,mif_".$_[3]."_".$probesLab[$i]."\n";
        }
    }
  }
  $str.="delete ".$tag."_mif".$_[2]."\n";
  # print PML3 "delete ".$tag."_mif".$_[2]."\n";
  return($str);
}
print NPML "remove hydrogens\nhide nonbonded\nshow sticks, HET\ndelete ".$tag."_1_nodes\ndelete ".$tag."_2_nodes";
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
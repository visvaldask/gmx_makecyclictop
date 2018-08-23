#!/usr/bin/perl -w
#
# by Visvaldas Kairys, Vilnius University, Institute of Biotechnology
# for questions/comments please send me an email visvaldas.kairys (at) bti.vu.lt
# see instructions below, or in the README.md file

use strict;

if(@ARGV!=1){
	print "The FIRST step a top file has to be created using pdb2gmx, with the pdb augmented with two extra residues:\n";
	print "1. the last ORIGINAL residue copypasted before the first ORIGINAL residue, with residue numbering ADJUSTED; and\n";
	print "2. the first ORIGINAL residue copypasted after the last ORIGINAL residue, with residue numbering ADJUSTED.\n";
	print "For example, if you want cyclic peptide FAIR, create a pdb file RFAIRF.\n";
	print "The atom numbering could be left as is (ie in a mess), pdb2gmx will sort it out.\n";
	print "Run the resulting pdb file through something like \"pdb2gmx -f myfile.pdb -o myfile_gmx.pdb\", but what you choose for the termini is not that important.\n";
	print "Just make sure terminal handling doesn't create extra residues, like small capping residues .\n";
	print "\n";
	print "This script deals with the top file created (will work also with itp), i.e. it is the SECOND step.\n";
	print "==================================================\n";
	die "Usage: $0 topfile\n";
}

#just one connecting bond closing the cycle is possible to handle with this script

#the last residue pasted in the beginning and the first in the end

#need to use the coords from the GROMACS generated pdb file, not the original


(my $root = $ARGV[0]) =~ s/(.+)\.[^.]+$/$1/;
my ($ext) = $ARGV[0] =~ /(\.[^.]+)$/;



my $modfile="$root" . "_cyc" . "$ext";
print "extension $ext\n";
print "old file $ARGV[0], newly created top file: $modfile\n";
############
my $connectlastdel='NA';  #number of connection atom in the last residue to be deleted

open(MODF,">$modfile") or die "Error while opening modfile $!\n";
open(TOPF,"<$ARGV[0]") or die "Error while opening $ARGV[0] $!\n";
my %at2res=();
my $numatomsfirstres=0;
my $orignumatoms=0;
#also get the first and last residues, which will be the ones to remove
my $firstres="NA";
my $lastres="NA";
my $curres="NA";
my %atomsbyres=();
while(<TOPF>){
	#print $_;
	#if(/^\[ atoms \]/){
	if(/^\[ atoms \]/ .. /^\s*$/){  #until blank line
		#print $_;
		unless(/^\[ atoms \]/ or /^\s*$/ or /^;/){
			#print $_;
			chomp; my @tmp=split;
			$curres=$tmp[2];
			$at2res{$tmp[0]}=$tmp[2];
			if(exists $atomsbyres{$tmp[2]}){
				$atomsbyres{$tmp[2]}+=1;
			}else{
				$atomsbyres{$tmp[2]}=1;
			}
			if($firstres eq "NA"){
				$firstres=$tmp[2];  #capture number of the first residue
			}
			#print "number $tmp[0] residue $tmp[2] atombyres count $atomsbyres{$tmp[2]} firstres $firstres\n";
			$orignumatoms++;
			#if($tmp[2] ==0){
			#	$numatomsfirstres=$numatomsfirstres+1;
			#}
			#search for connectivity atom numbers
			#if($tmp[2] ==0 and ($tmp[4] eq "$cycleatomname2")){
			#	print "number $tmp[0] residue $tmp[2] atomname $tmp[4] will be used to make a cycle\n"; #not used
			#}
			#if($tmp[2] ==5 and ($tmp[4] eq "$cycleatomname1")){
			#	print "number $tmp[0] residue $tmp[2] atomname $tmp[4] will be used to make a cycle\n";
			#	#$connectlastdel=$tmp[0];  #number of tail connection atom in the last residue to be deleted
			#}

		}
	}
	if(/^\[ bonds \]/){
		last;  #read atoms, need to go back for printing
	}
}
close(TOPF);
$lastres=$curres;
print "original number of atoms  $orignumatoms\n";
print "first res $firstres and last res $lastres will be removed, and a cycle will be made\n";
$numatomsfirstres=$atomsbyres{$firstres};
print "first res has $numatomsfirstres atoms!\n";
#need to add the number of atoms of the residues to be retained
my %tobedeleted=();
my $retained=0;
foreach my $key (keys %atomsbyres){
	unless ($key eq $firstres or $key eq $lastres){
		$retained+=$atomsbyres{$key};
		$tobedeleted{$key}=0;
	}else{
		$tobedeleted{$key}=1;
	}
}
print "retained atoms in the cyclic complex: $retained\n";
$connectlastdel=$orignumatoms-$atomsbyres{$lastres}+1; #the first atom number of the last residue in the old file
#print "num $num\n";
print "connectlastdel $connectlastdel\n"; 
my %newnumat=(); #renumbered atoms  stored here
foreach my $i (1 .. $orignumatoms){
	$newnumat{$i}=$i-$numatomsfirstres;
}
my $j=1;
foreach my $i ($connectlastdel .. $orignumatoms){
	$newnumat{$i}=$j;
	$j++;
}
#foreach my $i (1 .. $orignumatoms){
#	print "oldnumber $i new $newnumat{$i}\n";
#}

#now rewind
#seek(TOPF,0,0); #rewind
#second read instead of rewind
open(TOPF,"<$ARGV[0]") or die "Error while opening $ARGV[0] $!\n";
my $qtot=0;
while(<TOPF>){
	#print $_;
	if(/^\[ atoms \]/ .. /^\s*$/){  #until blank line
		unless(/^\[ atoms \]/ or /^\s*$/ or /^;/){
			#print $_;
			chomp; my @tmp=split;
			my $newat1=$newnumat{$tmp[0]};
			#print "number $tmp[0] newnumber $newat1 residue $tmp[2]\n";
			#unless($newat1 <= 0){
			unless($tobedeleted{$tmp[2]}){  #residues to be deleted
				#print "to be printed! number $tmp[0] newnumber $newat1 residue $tmp[2]\n";
				$qtot+=$tmp[6];
				#printf ("%6d %10s %6d %6s %6s %6d %10s %10s   ; qtot %10.3f\n",$newat1,$tmp[1],$tmp[2],$tmp[3],$tmp[4],$newat1,$tmp[6],$tmp[7],$qtot);
				printf MODF ("%6d %10s %6d %6s %6s %6d %10s %10s   ; qtot %10.3f\n",$newat1,$tmp[1],$tmp[2],$tmp[3],$tmp[4],$newat1,$tmp[6],$tmp[7],$qtot);
			}else{
				#print "atom skipped\n";
			}
		}else{
			print MODF $_;
		}
	}elsif(/^\[ bonds \]/ .. /^\s*$/){  #until blank line
		unless(/^\[ bonds \]/ or /^\s*$/ or /^;/){
			chomp; my @tmp=split;
			my $res1=$at2res{$tmp[0]};
			my $res2=$at2res{$tmp[1]};
			my $newat1=$newnumat{$tmp[0]};
			my $newat2=$newnumat{$tmp[1]};
			#print "bond: atom $tmp[0] to $tmp[1], residues $res1 to $res2\n";
			if($tobedeleted{$res1}+$tobedeleted{$res2} == 2){
				#print "To be deleted!\n";
			}
			if($tobedeleted{$res1}+$tobedeleted{$res2} == 1 ){
				#print $_;
				unless($newat1 <= 0 or $newat2 <= 0){
					#print "retained,making a cycle!\n";
					#print "new CYCLE line:\n";
					#printf ("%5d %5d %5d\n",$newat1,$newat2,$tmp[2]);
					printf MODF ("%5d %5d %5d ;cycle\n",$newat1,$newat2,$tmp[2]);
				}
			}
			if($tobedeleted{$res1}+$tobedeleted{$res2} == 0 ){
				#print $_;
				#print "retained,just need to renumber\n";
				#print "new line:\n";
				#printf ("%5d %5d %5d\n",$newat1,$newat2,$tmp[2]);
				printf MODF ("%5d %5d %5d\n",$newat1,$newat2,$tmp[2]);
			}
		}else{
			print MODF $_;
		}
	}elsif(/^\[ pairs \]/ .. /^\s*$/){  #until blank line
		unless(/^\[ pairs \]/ or /^\s*$/ or /^;/){
			chomp; my @tmp=split;
			my $res1=$at2res{$tmp[0]};
			my $res2=$at2res{$tmp[1]};
			my $newat1=$newnumat{$tmp[0]};
			my $newat2=$newnumat{$tmp[1]};
			#print "pairs: atom $tmp[0] to $tmp[1], residues $res1 to $res2\n";
			if($tobedeleted{$res1}+$tobedeleted{$res2} == 2){
				#print "To be deleted!\n";
			}
			if($tobedeleted{$res1}+$tobedeleted{$res2} == 1 ){
				#print $_;
				unless($newat1 <= 0 or $newat2 <= 0){
					#print "retained,making a cycle!\n";
					#print "new CYCLE line:\n";
					#printf ("%5d %5d %5d\n",$newat1,$newat2,$tmp[2]);
					printf MODF ("%5d %5d %5d ;cycle\n",$newat1,$newat2,$tmp[2]);
				}
			}
			if($tobedeleted{$res1}+$tobedeleted{$res2} == 0 ){
				#print $_;
				#print "retained,just need to renumber\n";
				#print "new line:\n";
				#printf ("%5d %5d %5d\n",$newat1,$newat2,$tmp[2]);
				printf MODF ("%5d %5d %5d\n",$newat1,$newat2,$tmp[2]);
			}
		}else{
			print MODF $_;
		}
	}elsif(/^\[ angles \]/ .. /^\s*$/){  #until blank line
		unless(/^\[ angles \]/ or /^\s*$/ or /^;/){
			chomp; my @tmp=split;
			my $res1=$at2res{$tmp[0]};
			my $res2=$at2res{$tmp[1]};
			my $res3=$at2res{$tmp[2]};
			my $newat1=$newnumat{$tmp[0]};
			my $newat2=$newnumat{$tmp[1]};
			my $newat3=$newnumat{$tmp[2]};
			#print "angles: atom $tmp[0] $tmp[1] $tmp[2], residues $res1 $res2 $res3\n";
			if($tobedeleted{$res1}+$tobedeleted{$res2}+$tobedeleted{$res3} == 3){
				#print "To be deleted!\n";
			}
			if($tobedeleted{$res1}+$tobedeleted{$res2}+$tobedeleted{$res3} == 1 or $tobedeleted{$res1}+$tobedeleted{$res2}+$tobedeleted{$res3} == 2 ){
				#print $_;
				unless($newat1 <= 0 or $newat2 <= 0 or $newat3 <= 0){
					#print "retained,making a cycle!\n";
					#print "new CYCLE line:\n";
					#printf ("%5d %5d %5d %5d\n",$newat1,$newat2,$newat3,$tmp[3]);
					printf MODF ("%5d %5d %5d %5d ;cycle\n",$newat1,$newat2,$newat3,$tmp[3]);
				}else{
					#print "cycle, but deleted!\n";
				}
			}
			if($tobedeleted{$res1}+$tobedeleted{$res2}+$tobedeleted{$res3} == 0 ){
				#print $_;
				#print "retained,just need to renumber\n";
				#print "new line:\n";
				#printf ("%5d %5d %5d %5d\n",$newat1,$newat2,$newat3,$tmp[3]);
				printf MODF ("%5d %5d %5d %5d\n",$newat1,$newat2,$newat3,$tmp[3]);
			}
		}else{
			print MODF $_;
		}
	}elsif(/^\[ dihedrals \]/ .. /^\s*$/){  #until blank line
		unless(/^\[ dihedrals \]/ or /^\s*$/ or /^;/){
			chomp; my @tmp=split;
			my $stuff="UNDEF";
			my $res1=$at2res{$tmp[0]};
			my $res2=$at2res{$tmp[1]};
			my $res3=$at2res{$tmp[2]};
			my $res4=$at2res{$tmp[3]};
			my $newat1=$newnumat{$tmp[0]};
			my $newat2=$newnumat{$tmp[1]};
			my $newat3=$newnumat{$tmp[2]};
			my $newat4=$newnumat{$tmp[3]};
			if (exists $tmp[5]){ #check params like "improper_O_C_X_Y"
				$stuff=$tmp[5];
			}else{
				$stuff="   ";  #make empty param if absent
			}
			#print "dihedrals: atom $tmp[0] $tmp[1] $tmp[2] $tmp[3], residues $res1 $res2 $res3 $res4\n";
			if($tobedeleted{$res1}+$tobedeleted{$res2}+$tobedeleted{$res3}+$tobedeleted{$res4} == 4){
				#print "To be deleted!\n";
			}
			if($tobedeleted{$res1}+$tobedeleted{$res2}+$tobedeleted{$res3}+$tobedeleted{$res4} != 4 and $tobedeleted{$res1}+$tobedeleted{$res2}+$tobedeleted{$res3}+$tobedeleted{$res4} != 0 ){
				#print $_;
				unless($newat1 <= 0 or $newat2 <= 0 or $newat3 <= 0 or $newat4 <= 0){
					#print "retained,making a cycle!\n";
					#print "new CYCLE line:\n";
					#printf ("%5d %5d %5d %5d %5d\n",$newat1,$newat2,$newat3,$newat4,$tmp[4]);
					#printf ("%5d %5d %5d %5d %5d %s\n",$newat1,$newat2,$newat3,$newat4,$tmp[4],$stuff);
					#printf MODF ("%5d %5d %5d %5d %5d ;cycle\n",$newat1,$newat2,$newat3,$newat4,$tmp[4]);
					printf MODF ("%5d %5d %5d %5d %5d %s ;cycle\n",$newat1,$newat2,$newat3,$newat4,$tmp[4],$stuff);
				}else{
					#print "cycle, but deleted!\n";
				}
			}
			if($tobedeleted{$res1}+$tobedeleted{$res2}+$tobedeleted{$res3}+$tobedeleted{$res4} == 0 ){
				#print $_;
				#print "retained,just need to renumber\n";
				#print "new line:\n";
				#printf ("%5d %5d %5d %5d %5d\n",$newat1,$newat2,$newat3,$newat4,$tmp[4]);
				#printf ("%5d %5d %5d %5d %5d %s\n",$newat1,$newat2,$newat3,$newat4,$tmp[4],$stuff);
				#printf MODF ("%5d %5d %5d %5d %5d\n",$newat1,$newat2,$newat3,$newat4,$tmp[4]);
				printf MODF ("%5d %5d %5d %5d %5d %s\n",$newat1,$newat2,$newat3,$newat4,$tmp[4],$stuff);
			}
		}else{
			print MODF $_;
		}
	}else{
			print MODF $_;
	}


}


close(TOPF);
close(MODF);

print "the modified topology file $modfile created.\n";
print "please use the pdb file created by pdb2gmx as described in the beginning, WITH THE EXTRA RESIDUES REMOVED!!!\n";

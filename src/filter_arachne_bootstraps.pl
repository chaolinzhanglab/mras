#!/usr/bin/perl -w

use strict;

use File::Basename;
use Carp;
use Data::Dumper;
use Getopt::Long;

use Common;

my $prog = basename ($0);

my $miThreshold = 0;
my $regulonSizeThreshold = 0;

GetOptions ("mi:f"=>\$miThreshold,
		"regulon-size:i"=>\$regulonSizeThreshold);

if (@ARGV != 2)
{
	print "$prog [options] <in.dir> <out.dir>\n";
	print " --mi           [float]: mutual information cutoff\n";
	print " --regulon-size [int]: set max regulon size cutoff\n"; 
	exit (1);
}


my ($inDir, $outDir) = @ARGV;

if (-f $outDir || -d $outDir)
{
	Carp::croak "$outDir already exists\n";
}
else
{
	system("mkdir $outDir") 
}

my $dh;

opendir ($dh, $inDir) || Carp::croak "cannot open dir $inDir to read\n";

my $i = 0;
my $iter = 0;
while (my $f = readdir($dh))
{
	next if $f eq '.' || $f eq '..';
	next unless $f=~/^bootstrapNetwork/;

	print "$i: $f ...\n";
	$i++;

	my ($fin, $fout);
	
	my %regulators;
	open ($fin, "<$inDir/$f") || Carp::croak "cannot open file $inDir/$f to read\n";
	
	my $header = <$fin>;
	chomp $header;

	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;

		my ($r, $t, $mi, $p) = split("\t", $line);
		
		if (not exists $regulators{$r})
		{
			$regulators{$r}{'iter'} = $iter;
			$iter++;
		}
		
		if ($miThreshold > 0)
		{
			next if $mi < $miThreshold;
		}

		push @{$regulators{$r}{'targets'}}, {t=>$t, mi=>$mi};
	}
	close ($fin);

	open ($fout, ">$outDir/$f") || Carp::croak "cannot open file $outDir/$f to write\n";
	
	print $fout $header, "\n";
	my @regulatorIds = sort {$regulators{$a}{'iter'} <=> $regulators{$b}{'iter'}} keys %regulators;
	foreach my $r (@regulatorIds)
	{
		my @targets = sort {$b->{'mi'} <=> $a->{'mi'}} @{$regulators{$r}{'targets'}};
	
		#Carp::croak Dumper (\@targets), "\n";

		my $regulonSize = @targets;

		#Carp::croak "size = $regulonSize\n";

		if ($regulonSizeThreshold > 0 && $regulonSize > $regulonSizeThreshold)
		{
			$regulonSize = $regulonSizeThreshold;
		}
		#Carp::croak "size = $regulonSize\n";

		for (my $k = 0; $k < $regulonSize; $k++)
		{
			my $t = $targets[$k];
			print $fout join ("\t", $r, $t->{'t'}, $t->{'mi'}), "\n";
		}
	}

	close ($fout);
}


closedir ($dh);

#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Basename;

my %options;
getopts("i:s:m:", \%options);
my $infile = $options{i} or &usage;
my $sampleId = $options{s} or &usage;
my $minAltCov = $options{m} || 0;
sub usage {die "USAGE: " . basename($0) . " [-i mpileup infile] [-s id] [-m minimum alternative coverage (default, 0)]\n";} 

print "##fileformat=VCFv4.1\n";
print "##source=$0\n";
print "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n";
print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";
print "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";

print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sampleId\n";

# Open pileup
open (IN, $infile) or die "ERROR: could not open $infile.\n";
while (<IN>) {
	chomp;

	# Ignore zero coverage (why is this even reported???)
	next if /^\S+\t\d+\t\w\t0/;

	# Extract info from pileup line.
	/^(\S+)\t(\d+)\t(\w)\t(\d+)\t(\S+)\t(\S+)$/ or die "ERROR: regex failed with line '$_'.";
	my $chr = $1;
	my $pos = $2;
	my $ref = uc($3);
	my $cov = $4;
	
	# get all alleles at current position.
	my @chars = getAlleles($5);
	
	# Now count each type of base encountered.
	my %data;
	$data{A} = 0;
	$data{C} = 0;
	$data{G} = 0;
	$data{T} = 0;
	$data{N} = 0;
	foreach my $char (@chars) {
		$char = $ref if $char =~ /[\.\,]/;
		$data{uc($char)}++;       
                }

	# calculate true coverage.
	my $trueCov = $data{A} + $data{C} + $data{G} + $data{T} + $data{N};

	# determine coverage of ref allele.
	my $refCov = $data{uc($ref)};
	
	# determine alt allele and coverage of alt allele. (define alt allele as next highest freq)
	
	# First assume no alt:
	my $alt = $ref;
	my $altCov = 0;
	
	# next determine alt and its frequency.
	# Ignore N
	foreach my $base (qw(A C G T)) {
		next if $base eq $ref;	
		if ($data{$base} > $altCov) {
			$altCov = $data{$base};
			$alt = $base;
			}
		}

	# only print variant if exceeds threshold. 
	next if $altCov < $minAltCov;

	# Don't print if no coverage.
	next if $trueCov == 0;

	my $altFreq = $altCov / $trueCov;
	my $refFreq = $refCov / $trueCov;

	# Determine genotype (crude determination)
	my $GT = "0/0";
	
	if ($altFreq > 50) {
		$GT = "1/1";
		}
	elsif ($altFreq > 0) {
		$GT = "0/1";
		}
	

	printf(
		"%s\t%s\t%s\t%s\t%s\t%s\t%s\tDP=%d;AF=%.3f\tGT:AD:DP\t%s:%d,%d:%d\n",
		$chr,	# CHR
		$pos,	# POS
		".", 	# ID
		$ref,	# REF
		$alt,	# ALT
		$cov,	# QUAL
		"PASS",	# FILTER
		# INFO
		$trueCov, # DP
		$altFreq, # AF
		# FORMAT
		$GT,
		$refCov,
		$altCov,
		$trueCov
		);

	}
close (IN) or die "ERROR: could not close $infile.\n";

###############################################################################

sub getAlleles_DEFUNCT {

        my $alleleString = shift;


        # Remove non-base characters:
        # Start of read symbol plus ascii representation of mapping quality.
        $alleleString =~ s/\^.//g;
        # End of read symbol.
        $alleleString =~ s/\$//g;

	while ($alleleString =~ /(\+\d+)/g) {
		warn "- $1\n";
		}

        return split("", $alleleString);

        } # End of method.


###############################################################################

sub getAlleles {
	
	my $string = shift;
	my @old = split("", $string);


	my @new;

	my $chr;
	while ($chr = shift(@old)) {
		
		# Ignore read starts (with associated quality score)
		if ($chr eq "^") {
			shift(@old);
			}

		# Ignore read ends.
		elsif ($chr eq "\$") {
			# Do nothing	
			}
	
		# Dealing with insert
		#-------------------------------------------------------------
		elsif ($chr eq "+") {
			process("insertion", \@old);
			# Do nothing with insert for now
			}
		#------------------------------------------------------------	

		# Dealing with deletion
		#------------------------------------------------------------
		elsif ($chr eq "-") {
			process("deletion", \@old);
			# Do nothing with deletion for now
			}
		#------------------------------------------------------------
		

		# Store valid bases
		else {
			push(@new, $chr);
			}


		}

	return @new;

	} # Ends of method.

###############################################################################

sub process {
	
	# See indelFreqs.pl for full use of this method. In current situation
	# only used to remove indel characters from string.

        my $type = shift; # UNUSED
        my $array = shift;

        # Get size of indel (and remove number from array)
        my $size = "";
        while ($array->[0] =~ /\d/) {
                $size .= shift(@$array);
                }

        # Get indel (and remove from array)
        my $indel = "";
        for (my $j = 0; $j < $size; $j++) {
                $indel .= shift(@$array);
                }

        # Store indel (upper case only)
       	# SEE indelFreqs.pl tool for example of what can be done with this info.
	# FOR CURRENT TOOL, INFO NOT REQUIRED. 

        } # End of method.

###############################################################################


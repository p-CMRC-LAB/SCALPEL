#!/usr/bin/perl -w

######################################
# code written by mireya plass       #
# contact: mplass at idibell dot cat #
######################################

use strict;
use Getopt::Long;



my $fasta_file;
my $strand = 1;
my $help;

&GetOptions('f=s'   => \$fasta_file,  ## input annotation file
	    's:1'   => \$strand,      ## strand
	    'help'  => \$help,        ## help
            );

## define options
if (defined($help) || ! $fasta_file){
    print STDERR "
    USAGE: ./find_polyAs.pl -a fasta_file
      -f:     Multifasta file to read. COMPULSORY

      -s:     integer defining whether only sense (0) or both strands (1) should be scanned
              Default: both strands

      -h:     shows this help\n\n";
    
    exit (1);
}

my %seq = %{read_fasta ($fasta_file)};

foreach my $chr (keys %seq){
    my $seq = $seq{$chr};
    my $splus=$seq;
    $splus=~tr/ATGCN/10000/;
    my @plus  = @{find_matches ($splus,0)};
    &print_uniq (\@plus, $chr,"+");
  
    if ($strand == 1){
	my $sminus=$seq;
	$sminus=~tr/ATGCN/01000/;
	my @minus  = @{find_matches ($sminus,1)};
	&print_uniq (\@minus, $chr,"-");
    }
}




########################################################################################
### FUNCTIONS ##########################################################################
########################################################################################

sub read_fasta{
    my $file =shift;
    my $seq ="";    
    my $id; 
    my %hash;
    my $i = 0;
    
    open (ARX, "<$file") || die ("cannot open $file for reading");

    while (<ARX>){
	chomp;
	if ($_=~m/^>/){
	    if ($i >0){
		$hash{$id} = $seq; 
		$seq= ""; 
	    }
	    my @p = split;
	    $id=$p[0];
	    $id =~s/>//;
	    $i=1; 
	}
	elsif ($_=~m/\w+/){
	    $seq .=$_;    
	}    
    }
    close (ARX);
    $hash{$id}= $seq;
    
    return \%hash;
}

########################################
sub find_matches{
    my @seq = split (//,$_[0]);
    my $strand = $_[1];
    my $len = 9;
    my $seqlen= scalar (@seq) -$len;
    
    my @matches;
    my $this_sum=0;

    for (my $i = 0; $i < $seqlen; $i++){
	my $save=0;
	if ($i == 0){
	    $this_sum =sum(@seq[$i..$i+$len]);
	}
	else{
	    $this_sum+= $seq[$i+$len]-$seq[$i-1];
	}
	if ($this_sum >= 6){
	    if (($seq[$i] == 1 && $strand == 0 ) || ($seq[$i+$len] == 1 && $strand == 1)){ ## possible match
		$save= check_value(@seq[$i..$i+$len]);
		if ($save == 1){
		    push (@matches, $i);
		}
	    }
	}
    }
    return \@matches;
}

########################################
sub sum{
    my @array = @_;
    my $total = 0;

    foreach my $a (@array){
	$total +=$a;
    }
    return $total;
}

########################################
sub print_uniq{
    my @m = @{$_[0]};
    my $chr = $_[1];
    my $strand = $_[2];
    my $len = 9;

    my @final;

    for (my $i = 0; $i < @m; $i++){
	my @p = ($m[$i],$m[$i]+10);
	if ($i == 0){
	    push (@final, \@p);
	}
	elsif($p[0] <= $final[-1][1]){
	    $final[-1][1]  = $p[1];
	}
	else{
	    push (@final, \@p);
	}
    }

    foreach my $f (@final){
	print join ("\t", $chr, @{$f},".",".",$strand), "\n";
    }
}

########################################
sub check_value{
    my @v = @_;
    my $all = 0;
    my $max_all= 0;

    for (my $j = 0; $j <@v; $j++){
	if ($v[$j] == 1){
	    $all++;
	    if ($all > $max_all){
		$max_all = $all;
	    }
	}
	else{
	    $all=0;
	}
    }
    if ($max_all == 6){
	return 1;
    }
    else{
	return 0;
    }
}


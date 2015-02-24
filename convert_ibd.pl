#!/usr/local/bin/perl
#
# Cared for by Filipe G. Vieira <>
#
# Copyright Filipe G. Vieira
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

    convert_ibd.pl v0.0.1

=head1 SYNOPSIS

    perl file_utils.pl [-h] -ind_file /path/to/ind_file -pos_file /path/to/pos_file [-ibd_pos_file /path/to/ibd_pos_file] [-ibd_coord_file /path/to/ibd_coord_file]

    OPTIONS:
       -help            This help screen
       -ind_file        File with individual information (ID on first column)
       -pos_file        TSV file with position information (CHR,POS)
       -ibd_pos_file    File with regions as IBD state per position (one indiv per line, one 0/1 per site)
       -ibd_coord_file  BED file with the coordinates of the IBD regions (CHR,START,END,ID)

=head1 DESCRIPTION

    This script will convert between the two IBD formats, 
    that is, between the one based on region coordinates 
    and the one based on 0/1 per position.

=head1 AUTHOR

    Filipe G. Vieira - fgarrettvieira _at_ gmail _dot_ com

=head1 CONTRIBUTORS

    Additional contributors names and emails here

=cut

# Let the code begin...

use strict;
use warnings;
use Getopt::Long;
$| = 1;

my ($ind_file, $pos_file, $ibd_pos_file, $ibd_coord_file, $ids);
my ($s, $ind_id, $chr, $start_pos, $end_pos, @buf, @inds, @sites, %ibd);

#Get the command-line options
&GetOptions('h|help'            => sub { exec('perldoc',$0); exit; },
	    'i|ind_file=s'      => \$ind_file,
            'p|pos_file=s'      => \$pos_file,
            'ibd_pos_file:s'    => \$ibd_pos_file,
	    'ibd_coord_file:s'  => \$ibd_coord_file,
    );



if($ibd_pos_file && $ibd_coord_file) {
    print(STDERR "ERROR: both IBD_POS and IBD_COORD files provided. What do you want to do?");
    exit(-1)
}



# Read POS file
open(FILE, $pos_file) or die("ERROR: cannot read POSITIONS file: $pos_file\n");
while(<FILE>) {
    chomp();
    my ($chr, $pos) = split(/[\t ]/);
    push(@sites, {'chr' => $chr, 'pos' => $pos});
}
my $n_sites = $#sites+1;
close(FILE);



# Read IND file
open(FILE, $ind_file) or die("ERROR: cannot read INDIVIDUALS file: $ind_file\n");
while(<FILE>) {
    my ($ind) = split(/[\t ]/);
    chomp($ind);
    push(@inds, $ind);
}
close(FILE);



if($ibd_pos_file) {
    my $curr_ind = -1;
    # IBD_POS file provided, will convert to REG...
    open(FILE, $ibd_pos_file) or die("ERROR: cannot read IBD_POSITIONS file: $ibd_pos_file\n");
    while(<FILE>) {
	chomp();
	# Skip lines starting with "//" (IBD file first line)
	next if(m/^\/\//);
	# Update current indiviudal
	$curr_ind++;
	# Skip individual if there is no valid ID
	next unless ($inds[$curr_ind]);
	# Stop reading file when no more indiv
	last if($curr_ind > $#inds);
	$s = 0;

	while( (($s = index($_, "1", $s)) != -1) ) {
	    $ind_id = $inds[$curr_ind];
	    $chr = $sites[$s]{'chr'};
	    $start_pos = $sites[$s]{'pos'}-1;
	    
	    for(; $s<=$#sites; $s++) {
		if($s == $#sites || $sites[$s+1]{'chr'} != $chr || substr($_, $s+1, 1) == 0) {
		    $end_pos = $sites[$s]{'pos'};
		    print($chr."\t".$start_pos."\t".$end_pos."\t".$ind_id."\t".($end_pos-$start_pos)."\n");
		    $s++;
		    last;
		}
	    }
	}
    }
    close(FILE);
} elsif($ibd_coord_file) {
    # IBD_COORD file provided, will convert to POS...
    $ibd{$_} = "0"x$n_sites for (@inds);

    open(FILE, $ibd_coord_file) or die("ERROR: cannot read IBD_COORDINATES file: $ibd_coord_file\n");
    while(<FILE>) {
	chomp();
	($chr, $start_pos, $end_pos, $ind_id) = split(/[\t ]/);

	next unless($ind_id && $chr && $start_pos && $end_pos);
	next unless( grep {$_ eq $ind_id} keys(%ibd) );
	# Increment start_pos by one since BED files are 0-based
	$start_pos++;

	for($s=0; $s<=$#sites; $s++) {
	    if($sites[$s]->{'chr'} == $chr && $sites[$s]->{'pos'} >= $start_pos && $sites[$s]->{'pos'} <= $end_pos) {
		substr($ibd{$ind_id}, $s, 1, "1");
	    }
	}
    }
    close(FILE);

    # Print to STDOUT
    print($ibd{$_}."\n") for (@inds);

} else {
    print(STDERR "ERROR: no IBD_POS or IBD_COORDINATES files provided. What do you want to do?");
    exit(-1)
}

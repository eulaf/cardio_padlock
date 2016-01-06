#!/usr/bin/env perl

=head1 NAME

add_blast.pl -- add number of blast matches to primers.txt

=head1 SYNOPSIS

B<add_blast.pl> <all_primers.txt and blastn outfmt 7 file(s)>

Run blast locally using output format 7.

blastn -db ucsc.hg19.fasta -evalue 0.05 -max_target_seqs 15 -outfmt 7 -task blastn -query <sample.all_primers.fa> -out <sample.all_primers.blastn7.txt>

=cut

#############################################################################

=head1 REQUIRES

Cwd,
FileHandle,
File::Basename,
File::Spec,
Getopt::Long,
Pod::Usage,
#Common,

=cut

use Cwd qw/abs_path/;
use FileHandle;
use File::Basename;
use File::Spec;
use FindBin;
use Getopt::Long;
require Pod::Usage;
use strict;
use warnings; no warnings 'once';

#############################################################################

### CONSTANTS

my $PROGRAM = $0;
$PROGRAM =~ s;^.*/;;;

my $VERSION = '$Revision: 1.2 $';
$VERSION =~ s/^\$Revision: //;
$VERSION =~ s/\$$//;
$VERSION = 0 unless $VERSION;
$VERSION = sprintf "%3.1f", $VERSION/10;

#############################################################################

=head1 DESCRIPTION

Add number of blast matches to all_primers.txt files. The all_primers.txt files
are created by call_primers.pl.  Blast should be run using a command such as:

blastn -db ucsc.hg19.fasta -evalue 0.01 -max_target_seqs 25 -outfmt 7 -task blastn -query <sample.all_primers.fa> -out <sample.all_primers.blastn7.txt>


=cut

#############################################################################

=head1 OPTIONS

=over 4

=item B<--outdir> <dir>

Directory to save output files.  Default is same location as input
primer file.

=item B<--dbug>

Print debugging messages.

=back

=cut


sub setup {
    Pod::Usage::pod2usage(-verbose => 1) unless @ARGV;

    my $man = 0;
    my $help = 0;

    Getopt::Long::Configure('pass_through');
    Getopt::Long::GetOptions('help|?' => \$help, man => \$man);
    Pod::Usage::pod2usage(-verbose => 1,
                          -exitstatus => 0) if $help;
    Pod::Usage::pod2usage(-verbose => 2,
                          -exitstatus => 0) if $man;

    my $p = {
        OUTDIR => 0,
        FORCE => 0,
        DBUG => 0,
    };

    Getopt::Long::Configure('no_pass_through');
    Getopt::Long::GetOptions(
        'version' => sub { print "$PROGRAM, version $VERSION\n\n";
            Pod::Usage::pod2usage(-verbose => 0,
                                  -exitstatus => 0) },
        'outdir=s' => \$$p{OUTDIR},
        'force' => \$$p{FORCE},
        'dbug' => \$$p{DBUG},
    ) || die "\n";
    @ARGV || die "Need all_primers.txt and blast file(s)\n";
    my @blastfiles = grep(/blast/i, @ARGV);
    @blastfiles || die "No blast files\n";
    my @primerfiles = grep(/primers.txt/i, @ARGV);
    @primerfiles || die "No primer files\n";
    my $infiles = pair_inputfiles(\@blastfiles, \@primerfiles);
    ($infiles, $p);
}

sub pair_inputfiles {
    my $blastfiles = shift;
    my $primerfiles = shift;

    my %infile;
    foreach my $input (['primers', $primerfiles],
                      ['blast', $blastfiles]) {
        my $tag = $$input[0];
        my $filelist = $$input[1];
        foreach my $fname (@$filelist) {
            my $base = basename($fname);
            $base =~ s/\.all_primers.*$//;
            $infile{$base}{$tag} = $fname;
        }
    }
    my $num = 0;
    my $tot = keys %infile;
    foreach my $k (keys %infile) {
        my @have = keys %{ $infile{$k} };
        if (scalar(@have)==2) { $num++; } 
        else {
            warn "$k only have @have\n";
            delete $infile{$k}
        }
    }
    warn scalar(keys %infile)." paired input files\n";
    return \%infile;
}

################################################################################

### MAIN PROGRAM

my ($infiles, $p) = setup();

my %dbsnp;
foreach my $base (sort keys %$infiles) {
    my $pfile = $$infiles{$base}{primers};
    my $bfile = $$infiles{$base}{blast};
    my $blastnums = parse_blast($bfile);
    add_blastnums($pfile, $blastnums, $p);
}

##############################################################################

sub add_blastnums {
    my $pfile = shift;
    my $blastnums = shift;

    my $outfile = $pfile;
    $outfile =~ s/\.txt$//;
    $outfile .= ".blast_counts.txt";
    if ($$p{OUTDIR}) {
        my $base = basename($outfile);
        $outfile = File::Spec->catfile($$p{OUTDIR}, $base);
    }
    warn "Writing $outfile\n";
    my $ofh = new FileHandle ">$outfile" or die ">$outfile: $!";
    my $fh = new FileHandle "$pfile" or die ">$pfile: $!";
    my $head = $fh->getline;
    my @fields = split(/\t/, $head);
    my @outfields;
    foreach my $f (@fields) {
        push(@outfields, $f);
        if ($f =~ /num_snps/) {
            push(@outfields, 'num_blast_hits_left');
            push(@outfields, 'num_blast_hits_right');
        }
    }
    print $ofh join("\t", @outfields);
    foreach my $l ($fh->getlines) {
        my @v = split(/\t/, $l);
        my %d = map { $_=>shift @v } @fields;
        my $pname = join("_", $d{target_name}, $d{segment},
                     $d{primer_num}.$d{template_strand});
        foreach my $side ('_left', '_right') {
            $d{'num_blast_hits'.$side} = $$blastnums{$pname.$side}{hits};
        }
        my @row = map { defined($d{$_}) ? $d{$_} : '' } @outfields;
        print $ofh join("\t", @row);
    }
    $fh->close;
    $ofh->close;
}

sub parse_blast {
    my $bfile = shift;

    my $fh = new FileHandle $bfile or die "$bfile: $!";
    my %bcounts;
    my $pname = '';
    foreach my $l (<$fh>) {
        if ($l =~ /^# Query: (\S+)/) {
            $pname = $1;
        } elsif ($l =~ /^# (\d+) hits found/) {
            $bcounts{$pname}{hits} = $1 if $pname;
            $pname = '';
        }
    }
    $fh->close;
    return \%bcounts;
}

#-----------------------------------------------------------------------------

sub rev_complement {
    my $seq = shift;

    $seq =~ tr/ATCGatcg/TAGCtagc/;
    $seq = reverse($seq);
    return $seq;
}

sub enumdate {
    my $long = shift;
    my ($se,$mi,$hr,$day,$mo,$yr) = localtime(time);
    my $enumdate = sprintf "%04d%02d%02d", $yr+1900, $mo+1, $day;
    unless ($long) { $enumdate =~ s/^..//; }
    ($enumdate);
}


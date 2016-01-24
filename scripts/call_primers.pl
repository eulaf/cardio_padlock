#!/usr/bin/env perl

=head1 NAME

call_primers.pl -- call primers using primer3

=head1 SYNOPSIS

B<call_primers.pl> <fa file(s)>

B<call_primers.pl> *.gene.fa *.sequence.fa

=cut

# 1/6/2016 Do not need to check for RE sites. Remove optimum product size.
#          Change range to 300-425 instead of 350-425 bp
#          Change flank to 10 instead of 20 bp
#############################################################################

=head1 REQUIRES

Cwd,
FileHandle,
File::Basename,
File::Spec,
Getopt::Long,
Pod::Usage

Also requires the primer3 program found at http://primer3.sourceforge.net.

=cut

use Cwd qw/abs_path/;
use FileHandle;
use File::Basename;
use File::Spec;
use Getopt::Long;
require Pod::Usage;
use strict;
use warnings; no warnings 'once';

BEGIN {
    my ($base, $path, $suff) = fileparse(abs_path($0));
    my $up = File::Spec->updir();
    my $lib = File::Spec->catfile(abs_path($path), $up, 'lib'); 
    my $lib2 = File::Spec->catfile(abs_path($path), 'lib'); 
    unshift(@INC, 'lib', $lib, $lib2);
}

#############################################################################

### CONSTANTS

my $PROGRAM = $0;
$PROGRAM =~ s;^.*/;;;
my $VERSION = 1.2;

#############################################################################

=head1 DESCRIPTION

Design primers where left and right primers are both on same strand
for every exon in given genes.  Uses primer3 (http://primer3.sourceforge.net)
to design primers.  By default, this script expects to find the primer3_config
directory in /opt. 

If the size of the exon plus primers is less than the desired amplicon size
(350-400 bp), the program will attempt to find primers that amplify the
entire target region.  If no primers are found, the target region will be
divided into 2 segments and the program will try to find two pairs of primers
that cover the target region.

If the target sequence is larger than the desired amplicon size, the
sequence will be divided into a minimal number of segments and primers found 
to amplify each segment.

=head2 Input files

Input sequences should be the *.gene.fa files created by get_exons_fa.pl, 
*.sequence.fa files created by get_region_fa.pl or fasta files with the same 
format.  These fasta files contain sequences in lowercase and uppercase.  
The sequences in uppercase are the target sequences to be amplified.  
The lowercase sequences are flanking sequences to the target sequence.

CHANGE 2015/12: To work with sequences created by RepeatMasker (which creates
sequences in all uppercase), the input fasta files should alternate between
flanking sequence and target sequence, starting and ending with flanking 
sequence.

The *.snps.txt files created by the get_exons_fa.pl program are also used.
They should be located in the same directory as the *.gene.fa and 
*.sequence.fa files.

=head2 Output files

A directory called "p3files" will be created to store primer3 input files
(*.p3in.txt) and output files (*.p3out.txt).  If you rerun the program
with the same label as a previous run, then the program will parse existing
p3out.txt files rather than running primer3 again.  

By default, output files will be labelled with today's date YYMMDD.
For each target sequence, the program finds up to 50 primer pairs.

=item YYMMDD.primers_summary.txt

This tab-delimited file reports number of segments used to cover each exon
and number of segments not covered by primers for all genes.

=item YYMMDD.primers_found.txt

This tab-delimited file reports number of primers found for each segment
for each exon for all genes.

=item YYMMDD.geneID.primers_found.txt

This tab-delimited file reports the number of primers found for each target
segment for all exons and genes.

=item YYMMDD.geneID.results.txt

Similar to YYMMDD.results.txt but specific for this geneID.  This 
tab-delimited file reports the number of primers found for each target
segment for all exons in geneID.

=item YYMMDD.geneID.all_primers.txt

This tab-delimited file lists all primer pairs and their sequences, amplicon
sequence, location in human genome, number of snps found in the primers, snp 
names, product size, primer3 penalty score, GC content, and more.

=item YYMMDD.geneID.all_primers.fa

This fasta file contains all primer sequences for geneID.

=item YYMMDD.geneID.all_amplicons.fa

This fasta file contains all amplicon sequences for the geneID.  

=cut

#############################################################################

=head1 OPTIONS

=over 4

=item B<--od>

=item B<--outdir>

Output directory for files.

=item B<--label> <str>

Label to add to output primer3 files.  Default label
is the date YYMMDD.

=item B<--force>

Overwrite existing primer3 output files.

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

    my %P3_PARAMS = (
        PRIMER_PICK_LEFT_PRIMER=>0,
        PRIMER_PICK_INTERNAL_OLIGO=>0,
        PRIMER_PICK_RIGHT_PRIMER=>0,
        PRIMER_MAX_NS_ACCEPTED=>0,
        PRIMER_EXPLAIN_FLAG=>1,
        PRIMER_MIN_TM=>50,
        PRIMER_OPT_TM=>65,
        PRIMER_MAX_TM=>87,
    );

    my $p = {
        P3EXE => 'primer3_core',
        P3CONFIG_PATH => '/opt/primer3_config/',
        PRODUCT_MIN_SIZE => 300,
        PRODUCT_OPT_SIZE => 375,
        PRODUCT_MAX_SIZE => 425,
        PRIMER_MIN_SIZE=>18,
        PRIMER_OPT_SIZE=>23,
        PRIMER_MAX_SIZE=>28,
        PRIMER_NUM_RETURN => 30,
        PRIMER_NUM_RETURN_right => 30,
        PRIMER_NUM_RETURN_left => 15,
        OUTDIR => undef,
        OUTLABEL => enumdate(),
        TARGET_FLANK => 10, # num bases flanking target 
                            # region to include in amplicon
        P3FILES_SUBDIR => 'p3files',
        MIN_MAF => 0.10,
        FORCE => 0,
        P3PARAMS => \%P3_PARAMS,
        DBUG => 0,
    };

    Getopt::Long::Configure('no_pass_through');
    Getopt::Long::GetOptions(
        'version' => sub { print "$PROGRAM, version $VERSION\n\n";
            Pod::Usage::pod2usage(-verbose => 0,
                                  -exitstatus => 0) },
        'outdir|od=s' => \$$p{OUTDIR},
        'label=s' => \$$p{OUTLABEL},
        'force' => \$$p{FORCE},
        'dbug' => \$$p{DBUG},
    ) || die "\n";
    $$p{FLANKSIZE} = $$p{PRODUCT_MAX_SIZE} - 2*$$p{TARGET_FLANK} -
        $$p{PRIMER_MIN_SIZE};
#    $P3_PARAMS{PRIMER_PRODUCT_OPT_SIZE} = $$p{PRODUCT_OPT_SIZE};
    $P3_PARAMS{PRIMER_PRODUCT_SIZE_RANGE} = 
        sprintf "%d-%d", $$p{PRODUCT_MIN_SIZE}, $$p{PRODUCT_MAX_SIZE};
    $P3_PARAMS{PRIMER_THERMODYNAMIC_PARAMETERS_PATH} = $$p{P3CONFIG_PATH};
    $P3_PARAMS{PRIMER_MIN_SIZE} = $$p{PRIMER_MIN_SIZE};
    $P3_PARAMS{PRIMER_OPT_SIZE} = $$p{PRIMER_OPT_SIZE};
    $P3_PARAMS{PRIMER_MAX_SIZE} = $$p{PRIMER_MAX_SIZE};
    if ($$p{OUTDIR}) {
        $$p{P3FILES_SUBDIR} = File::Spec->catfile($$p{OUTDIR},
                $$p{P3FILES_SUBDIR});
    }
    unless (-d $$p{P3FILES_SUBDIR}) {
        mkdir ($$p{P3FILES_SUBDIR}) || die "$$p{P3FILES_SUBDIR}: $!";
    }
    @ARGV || die "Need FASTA file(s)\n";
    my @fafiles;
    foreach my $f (@ARGV) {
        if ($f =~ /\.fa/i) {
            push(@fafiles, $f);
        } else {
            my $fh = new FileHandle $f or die "$f: $!";
            my @l = grep(/\w/, $fh->getlines);
            chomp @l;
            foreach my $l (@l) {
                my @f = <$l.gene.fa>;
#                print "$l\t".join(", ", @f)."\n";
                push(@fafiles, @f);
            }
            $fh->close;
        }
    }
    warn scalar(@fafiles). " input files\n";
    warn "Label for output files: $$p{OUTLABEL}\n";
    (\@fafiles, $p);
}
################################################################################

### MAIN PROGRAM

my ($fafiles, $p) = setup();

foreach my $fafile (@$fafiles) {
    my ($x, $filelabel) = split(/\./, basename($fafile, '.gene.fa',
                '.gene.fa.masked', '.sequence.fa')); 
    unless ($filelabel) { $filelabel = $x }
    my ($targets, $seqhash) = parse_fasta($fafile);
    parse_snpfile($fafile, $seqhash, $p);
    $targets = classify_exons($targets, $seqhash, $p);
    my $primers = pick_primers($targets, $seqhash, $filelabel, $p);
    create_exon_results($filelabel, $primers, $seqhash, $p);
    print_primers($filelabel, $primers, $seqhash, $p);
}
warn $$p{RES_FIELDS}."\n";
warn join("\n", @{ $$p{RESULTS} })."\n";

# final results file with number of primers found for each segment
my $outfile = "primers_found.txt";
if ($$p{OUTLABEL}) { $outfile = $$p{OUTLABEL}.".$outfile"; }
if ($$p{OUTDIR}) { $outfile = File::Spec->catfile($$p{OUTDIR}, $outfile) }
warn "\nWriting $outfile\n";
my $ofh = new FileHandle ">$outfile" or die ">$outfile: $!";
print $ofh $$p{RES_FIELDS}."\n";
print $ofh join("\n", @{ $$p{RESULTS} })."\n";
$ofh->close;

$outfile = "primers_summary.txt";
if ($$p{OUTLABEL}) { $outfile = $$p{OUTLABEL}.".$outfile"; }
if ($$p{OUTDIR}) { $outfile = File::Spec->catfile($$p{OUTDIR}, $outfile) }
warn "\nWriting $outfile\n";
$ofh = new FileHandle ">$outfile" or die ">$outfile: $!";
print $ofh join("\t", qw/gene target_name target_length_with_flank num_segments
            num_segments_missing_primers/)."\n";
my $tot_targets = 0;
my $tot_length = 0;
my $tot_segments = 0;
my $tot_missing = 0;
foreach my $row (@{ $$p{SUMMARY} }) {
    print $ofh join("\t", @$row)."\n";
    $tot_targets++;
    $tot_length += $$row[2];
    $tot_segments += $$row[3];
    $tot_missing += $$row[4];
}
print $ofh join("\t", "Totals", $tot_targets, $tot_length, $tot_segments,
        $tot_missing)."\n";
$ofh->close;


print "\nDone\n";

#-----------------------------------------------------------------------------

sub parse_fasta {
    my $fafile = shift;
    my $seqhash = shift || {};

    print STDERR "\nParsing $fafile: ";
    $/ = "\n>";
    my $is_masked = $fafile =~ /.masked$/;
    my $fh = new FileHandle $fafile or die "$fafile: $!";
    my @recs = $fh->getlines;
    $fh->close;
    $/="\n";

    print STDERR scalar(@recs)." sequences\n";
    my @targets;
    my $prevseqid = undef;
    my $flag = 1;
    foreach my $rec (@recs) {
        $rec =~ s/>//g;
        my @lines = split(/\n/, $rec);
        my $head = shift @lines;
        my $seq = join('', @lines);
        $seq =~ s/\s//g;
        my ($seqid, $loc, $desc) = split(/\s/, $head, 3);
        my @name = split(/:/, $seqid);
        my $name = pop @name;
        my ($chrom, $start, $end, $strand) = split(/:/, $loc);
        if ($is_masked) {
            # repeatmasker sequences are all uppercase
            # assume alternating non-coding/coding.
            $flag = !$flag;
        } else {
            # sequences to create primers for are in uppercase
            # while introns and flanking seqs are in lowercase
            $flag = $seq =~ /[ACTG]+/ ? 1 : 0;
        }
        push(@targets, $seqid) if $flag;
        $$seqhash{$seqid} = { need_primers=>$flag, prevseq=>$prevseqid,
            desc=>$desc, seq=>$seq, seqlen=>length($seq), name=>$name, 
            chrom=>$chrom, start=>$start, end=>$end, strand=>$strand,
        };
        if ($prevseqid) { $$seqhash{$prevseqid}{nextseq} = $seqid; }
        $prevseqid = $seqid;
    }
    (\@targets, $seqhash);
}

# Add any snp info to seqhash
sub parse_snpfile {
    my $fafile = shift;
    my $seqhash = shift;
    my $p = shift;

    my $snpfile = $fafile; $snpfile =~ s/.masked$//;
    $snpfile =~ s/(gene_?.*|targets|sequence).fa$/snps.txt/;
    return {} unless -f $snpfile;
    print STDERR "\nParsing $snpfile\n";
    my $fh = new FileHandle $snpfile or die "$snpfile: $!";
    my @lines = $fh->getlines;
    $fh->close;
    chomp @lines;
    my @fields = map { s/ .*//; $_ } split(/\t/, shift @lines);

    my %snps;
    my $num = 0;
    foreach my $l (@lines) {
        my @v = split(/\t/, $l);
        my %d = map { $_=>shift @v } @fields;
        if ($d{MAF} && $d{MAF} >=$$p{MIN_MAF}) {
            my $mafv = $d{MAF} ? $d{MAF} : 'n/a';
            $d{desc} = "$d{snp} $d{chrom}:$d{pos} (MAF: $mafv)";
            if ($$p{DBUG}) { print STDERR $d{desc}."\n"; }
            $num++;
            $snps{$d{chrom}}{$d{pos}} = \%d;
        }
    }
    warn "  Got $num snps\n";
    my $tot = 0;
    foreach my $seqid (keys %$seqhash) {
        my $snphash = $snps{$$seqhash{$seqid}{chrom}};
        next unless $snphash;
        my $start = $$seqhash{$seqid}{start};
        my $end = $$seqhash{$seqid}{end};
        my $seqsnps = snps_in_interval($snphash, $start, $end);
        my $numsnps = scalar(keys %$seqsnps);
        $tot += $numsnps;
        warn "  $seqid: $numsnps snps\n";
        $$seqhash{$seqid}{snps} = $seqsnps;
    }
    warn "  Total $tot snps placed\n";
    (\%snps);
}

sub snps_in_interval {
    my $snphash = shift;
    my $start = shift;
    my $end = shift;

    my @pos = grep ($_>=$start && $_<=$end, keys %$snphash);
    my %seqsnps = map { $_=>$$snphash{$_} } @pos;
    (\%seqsnps);
}

#-----------------------------------------------------------------------------

sub pick_primers {
    my $targets = shift;
    my $seqhash = shift;
    my $filelabel = shift; #for summary
    my $p = shift;

    warn "\nPicking primers\n";
    my $limit = scalar(@$targets);
    my $num = 0;
    my %ppairs;
    my $minpairs = 50; 
    foreach my $seqid (@$targets) {
        last if ++$num > $limit;
        warn "\n($num/$limit)==== SeqID $seqid\n";
        my ($targetlen, @seqs) = get_sequence_template($seqid, $seqhash, $p);
        my $minsize = $targetlen + 2*$$p{PRIMER_OPT_SIZE};
        my $numpieces = int($targetlen/$$p{PRODUCT_MAX_SIZE}+1);
        warn "  Targetlen $targetlen\t(minsize $minsize)\n";
        my %iterres;
        for(my $iter=0; $iter<50; $iter++) {
            my $segsize = $targetlen/$numpieces;
            last if $segsize < 100 && $iter>3;# stop if sizes get too small
            warn "  Finding primers for $numpieces pieces (iter=$iter)\n";
            my $strand = $$seqhash{$seqid}{strand};
            my $templates = $numpieces==1 ? [\@seqs] :
                split_sequence_template($p, $numpieces, $strand, $seqid, @seqs);
            my $numtem = scalar(@$templates);
            for(my $i=0; $i<@$templates; $i++) {
                my $let = letterLabel($i, $numtem);
                my $segment = "$numpieces$let";
                my $label = $numpieces==1 ? '' : "_$numpieces$let";
                my $seqs = $$templates[$i];
                my $primers = pick_primers_for_template($p, 'F', $seqid, 
                                          $label, $seqhash, @$seqs);
                ($primers, my $bestpenalty) = sort_and_clean_primerlist(
                                                   $primers, $minpairs);
                if (@$primers < $minpairs || $bestpenalty>3) {
                    # try to find more primers on reverse strand
                    my $primersR = pick_primers_for_template($p, 'R', $seqid,
                                              "_$segment", $seqhash, @$seqs);
                    ($primers) = sort_and_clean_primerlist([@$primers, 
                                                   @$primersR], $minpairs);
                }
                $iterres{$numpieces}{seg}{$segment} = $primers;
            }
            # determine num failures
            my $noprimers = 0;
            foreach my $seg (sort keys %{ $iterres{$numpieces}{seg} }) {
                my $numprimers = scalar(@{ $iterres{$numpieces}{seg}{$seg} });
                warn "$seqid\tNum pieces $numpieces\t$seg ($segsize bp)".
                     "\tPrimers found: $numprimers\n";
                unless ($numprimers) { $noprimers++ }
            }
            my $rate = $iterres{$numpieces}{fail} = $noprimers/$numpieces;
            warn "$seqid\tNum pieces $numpieces\tFailure rate: $rate\n";
            if ($noprimers==0) { # success - found primers for all segments
                $ppairs{$seqid} = $iterres{$numpieces}{seg};
                last; # stop iteration
            }
            $numpieces++; #increase number of pieces each iteration
        }
        unless ($ppairs{$seqid}) { # if none completely successful,  
            # assign results with lowest failures
            my @best =  sort { $iterres{$a}{fail}<=>$iterres{$b}{fail} }
                keys %iterres;
            warn "Best: $best[0] pieces\tRates: ".join(", ", @best)."\n";
            $ppairs{$seqid} = $iterres{$best[0]}{seg};
        }
        # report summary of results
        my $noprimers = 0;
        foreach my $segment (sort keys %{ $ppairs{$seqid} }) {
            my $numprimers = scalar(@{ $ppairs{$seqid}{$segment} });
            warn "Primers found: $seqid\t$segment\t$numprimers\n";
            unless ($numprimers) { $noprimers++ }
        }
        my $numtargets = scalar(keys %{$ppairs{$seqid}});
        warn "SeqID $seqid\tLength: $targetlen\tNum targets: $numtargets\t".
            "Missing primers: $noprimers\n";
        push(@{ $$p{SUMMARY} }, [ $filelabel, $seqid, $targetlen, 
                    $numtargets, $noprimers]);
    }
    return \%ppairs;
}

sub letterLabel {
    my $num = shift;
    my $totnum = shift;
        
    my $numlet = 0;
    while ($totnum>1) {
        $numlet++;
        $totnum /= 26;
    }
    my $let = '';
    for(my $i=0; $i<$numlet; $i++) {
        my $l = $num % 26;
        $let = chr(ord('a')+$l).$let;
        $num = int(($num-$l)/26);
    }
    return $let;        
}

sub get_sequence_template {
    my $seqid = shift;
    my $seqhash = shift;
    my $p = shift;

    warn "== Getting sequence template for $seqid\n";
    my $exonseq = $$seqhash{$seqid}{seq};
    my $exonlen = $$seqhash{$seqid}{seqlen};
    my $flanksize = $$p{FLANKSIZE};
    warn "  Exonlen $exonlen flank $flanksize\n";
    my %seqs = (prevseq=>'', nextseq=>'');
    foreach my $side ('prevseq', 'nextseq') {
        my $thisid = $seqid;
        while (length($seqs{$side})<$flanksize) {
            my $id = $$seqhash{$thisid}{$side};
            last unless $id;
            my $seq = $$seqhash{$id}{seq};
            warn "  $side: ".length($seq)." bp\t$id\n";
            if ($side eq 'prevseq') {
                $seqs{$side} = $seq . $seqs{$side};
            } else {
                $seqs{$side} .= $seq;
            }
            $thisid = $id;
        }
    }
    my $lfseq = substr($seqs{prevseq}, -$flanksize, $flanksize);
    my $rtseq = substr($seqs{nextseq}, 0, $flanksize);
    my $strand = $$seqhash{$seqid}{strand};
    my $lfchrpos = $strand > 0 ? 
        $$seqhash{$seqid}{start} - length($lfseq) :
        $$seqhash{$seqid}{end} + length($lfseq);
    my $lfflank = substr($lfseq, -$$p{TARGET_FLANK}, $$p{TARGET_FLANK});
    substr($lfseq, -$$p{TARGET_FLANK}, $$p{TARGET_FLANK}) = '';
    my $rtflank = substr($rtseq, 0, $$p{TARGET_FLANK});
    substr($rtseq, 0, $$p{TARGET_FLANK}) = '';
    my $targetseq = $lfflank.$exonseq.$rtflank;
    my $targetlen = length($targetseq);
    printf STDERR "  lfseq %d\ttarget %d\trtseq %d\tlfchrpos %d\n",
          length($lfseq), $targetlen, length($rtseq), $lfchrpos;
    if ($$p{DBUG}) {
        warn "  Left seq starting at $lfchrpos (strand $strand)\n";
        warn "$lfseq\n";
        warn "  targetseq\n$targetseq\n";
        warn "  Right seq\n$rtseq\n";
    }
    die "No lf seq for $seqid" unless $lfseq;
    die "No target seq for $seqid" unless $targetseq;
    die "No rt seq for $seqid" unless $rtseq;
    ($targetlen, $lfseq, $targetseq, $rtseq, $lfchrpos);
}

sub split_sequence_template {
    my $p = shift;
    my $numpieces = shift;
    my $strand = shift;
    my $seqid = shift;
    my $flankL = shift;
    my $targetseq = shift;
    my $flankR = shift;
    my $lfchrpos = shift;

    warn "== Splitting sequence template into $numpieces pieces\n";
    my $targetlen = length($targetseq);
    my $piecelen = int($targetlen/$numpieces+0.6);
    warn "  Each piece len: $piecelen\n";
    my $seq = $flankL.$targetseq.$flankR;
    my @templates;
    my $offset = length($flankL);
    my $flanksize = $$p{FLANKSIZE};
    my $piecepos = $lfchrpos;
    warn "  TargetLen $targetlen pieces $numpieces piecelen $piecelen".
        " lfchrpos $lfchrpos\n" if $$p{DBUG};
    for(my $i=0; $i<$numpieces; $i++) {
        warn "$i) targetlen $targetlen pieces $numpieces piecelen $piecelen\n";
        # if last piece, then use remainder of target
        my $lenP = ($i+1==$numpieces) ? $targetlen : $piecelen;
        my $offsetL = $offset-$flanksize; $offsetL = 0 if $offsetL<0;
        my $lenL = $offset-$offsetL;
        my $offsetR = $offset+$lenP;
        warn "$i) lenP $lenP offset $offset offsetR $offsetR\n" if $$p{DBUG};
        my $pieceL = substr($seq, $offsetL, $lenL);
        my $piece = substr($seq, $offset, $lenP);
        my $pieceR = substr($seq, $offsetR, $flanksize);
        $piecepos = $strand>0 ? $lfchrpos+$offsetL : $lfchrpos-$offsetL;
        my @template = ($pieceL, $piece, $pieceR, $piecepos);
        push(@templates, \@template);
        if ($$p{DBUG}) {
            warn "$i) L $offsetL, $lenL\n";
            warn "$i) T $offset, $piecelen\n";
            warn "$i) R $offsetR, $flanksize\n";
            warn "$i) remaining $targetlen lfchrpos $piecepos\n";
        }
        $offset += $piecelen;
        $targetlen -= $piecelen; #length of target remaining
    }
    if ($$p{DBUG}) {
        warn ">fullseq $lfchrpos\n$seq\n";
        for(my $i=0; $i<$numpieces; $i++) {
            my $t = $templates[$i];
            my ($lf, $ts, $rt, $pos) = @$t;
            warn ">piece$i $pos\n$lf\n$ts\n$rt\n";
        }
    }
    (\@templates);
}

#-----------------------------------------------------------------------------

sub pick_primers_for_template {
    my $p = shift;
    my $FR = shift;
    my $seqid = shift;
    my $segment = shift;
    my $seqhash = shift;
    my @seqs = @_;

    my $segid = $seqid.$segment;
    my $pcrprimers = pick_pcr_primer_pair($p, $FR, $segid, @seqs);
    my $primerpairs = convert_primers_pcr2padlock($pcrprimers);
    my $primerlist = primer_pair_info($p, $FR, $primerpairs, $segid,
                                      $$seqhash{$seqid}, @seqs);
    ($primerlist);
}

sub pick_pcr_primer_pair {
    my $p = shift;
    my $FR = shift;
    my $segid = shift;
    my $flankL = shift;
    my $targetseq = shift;
    my $flankR = shift;

    warn "== Picking primers on $FR for $segid\n";
    if ($$p{OUTLABEL}) {
        $$p{FILELABEL} = "$$p{OUTLABEL}.$segid.pcr-$FR";
        $$p{FILELABEL} =~ s/[\/:]/_/g;
        my $outfile = "$$p{FILELABEL}.p3out.txt";
        if ($$p{P3FILES_SUBDIR}) {
            $outfile = File::Spec->catfile($$p{P3FILES_SUBDIR}, $outfile);
        }
        if (-f $outfile && !$$p{FORCE}) {
            warn "  Parsing old output file $outfile\n";
            my $primers = parse_primer3_output($outfile, $p);
            return filter_primer_pairs($primers, $FR);
        }
    }
    # write primer3 parameters to input file
    my $infile = $$p{FILELABEL} ? "$$p{FILELABEL}.p3in.txt" :
                               "PRIMER_PARAMS.tmp";
    if ($$p{P3FILES_SUBDIR}) {
        $infile = File::Spec->catfile($$p{P3FILES_SUBDIR}, $infile);
    }
    if ($$p{DBUG}) { warn "  Writing param file $infile\n"; }
    my $pphi = new FileHandle ">$infile" or die ">$infile: $!";
    $pphi->autoflush(1);
    # compile parameters for primer3
    my $seq = $flankL.$targetseq.$flankR;
    my $target_start = length($flankL);
    my $target_len = length($targetseq);
    my %p3in = %{ $$p{P3PARAMS} };
    $p3in{SEQUENCE_ID} = $segid.$FR;
#    $p3in{SEQUENCE_TEMPLATE} = replace_REsite_with_N($seq);
    $p3in{SEQUENCE_TEMPLATE} = $seq;
    $p3in{SEQUENCE_TARGET} = "$target_start,$target_len";
    $p3in{PRIMER_PICK_LEFT_PRIMER} = 1,
    $p3in{PRIMER_PICK_RIGHT_PRIMER} = 1,
    $p3in{PRIMER_NUM_RETURN} = $$p{PRIMER_NUM_RETURN};
    # Note that for LPP primers, the right primer is on the same
    # strand as the left primer, which is the opposite of traditional
    # primer pairs as picked by primer3
    print_p3input_params($pphi, \%p3in);
    $pphi->close();
    my $primers = run_primer3($infile, $p);
    return filter_primer_pairs($primers, $FR);
}

# Do not need to check for RE sites 1/6/2016
sub replace_REsite_with_N {
    my $seq = shift;

    $seq =~ s/GAGTC|GACTC/NNNNN/gi; # MlyI
    $seq =~ s/GGTCTC|GAGACC/NNNNNN/gi; #BsaI
    return $seq
}

sub hasREsite {
    my $ps = shift;
    my $lfrt = shift;
#    $ps = $lfrt eq 'left' ?  "CTCG" . RC_seqstr($ps) . "GCGG" :
#        "TCTCA" . $ps . "GGCGC";
    if ($ps =~/(GAGTC)/i) {return "MlyI $1";} #MlyI
    elsif ($ps =~/(GACTC)/i) {return "MlyI $1";}  #MlyI
    elsif ($ps =~/(GGTCTC)/i) {return "BsaI $1";} #BsaI
    elsif ($ps =~/(GAGACC)/i) {return "BsaI $1";} #BsaI
    else {return '';}       
}

sub filter_primer_pairs {
    my $pprimers = shift;
    my $FR = shift;

    my @pprimers;
    warn "== Filtering primer pairs $FR\n";
    foreach my $k (sort keys %$pprimers) {
        my $numprimers = @{ $$pprimers{$k} };
        warn "  Total of $numprimers primers: $k\n" if $$p{DBUG};
        for(my $i=0; $i<@{ $$pprimers{$k} }; $i++) {
            my @reSite;
            my @pairseqs;
            foreach my $side ('left', 'right') {
                my $phash = $$pprimers{$k}[$i]{$side};
                my $pseq = $$phash{sequence};
                push(@pairseqs, $pseq);
#                my $reSite = hasREsite($pseq);
#                if ($reSite) { push(@reSite, "$side $pseq $reSite") }
            }
            if (@reSite) {
                my $num = $i+1;
                warn "$num) ".join(";\t", @reSite)."\n" if $$p{DBUG}
            } else {
                push(@pprimers, $$pprimers{$k}[$i]);
            }
        }
    }
    my $numgood = @pprimers;
    warn "  Total of $numgood good primer pairs\n";
    (\@pprimers);
}

sub convert_primers_pcr2padlock {
    my $pcrprimers = shift;

    foreach my $pcr (@$pcrprimers) {
        my $rseq = rev_complement($$pcr{right}{sequence});
        $$pcr{right}{sequence} = $rseq;
    }
    return $pcrprimers;
}

sub primer_pair_info {
    my $p = shift;
    my $FR = shift;
    my $ppairs = shift;
    my $segid = shift;
    my $seqinfo = shift;
    my $flankL = shift;
    my $targetseq = shift;
    my $flankR = shift;
    my $leftchrpos = shift;

    my $seqtemplate = $flankL.$targetseq.$flankR;
    my $plusstrand = $$seqinfo{strand} > 0 ? 1 : 0;
    my $start = $$seqinfo{start};
    my $end = $$seqinfo{end};
    my $snphash = $$seqinfo{snps} || {};
#    warn "$segid: plusstrand? $plusstrand start $start end $end\n";
    my $num = 0;
    my @pairs;
    foreach my $pp (@$ppairs) {
        my $lseq = $$pp{left}{sequence};
        my $lstart = $$pp{left}{start}; # zero-based indexing
        my $llen = $$pp{left}{len};
        my $lchrstart = $$pp{left}{chr_start} = $plusstrand ? 
            $leftchrpos + $lstart : $leftchrpos - $lstart;
        my $lchrend = $plusstrand ? $lchrstart+$llen-1 : $lchrstart-$llen+1;

        my $rseq = $$pp{right}{sequence};
        my $rstart = $$pp{right}{start}; # zero-based indexing
        my $rlen = $$pp{right}{len};
        my $rchrend = $$pp{right}{chr_end} = $plusstrand ? 
            $leftchrpos + $rstart : $leftchrpos - $rstart;
        my $rchrstart = $plusstrand ? $rchrend-$rlen+1 : $rchrend+$rlen-1;
        my $amplicon_seq = get_amplicon($seqtemplate, $lstart, $rstart);
        my $amplicon_gc = gc_perc($amplicon_seq);
        my $amplicon_len = length($amplicon_seq);
        if ($amplicon_len != $$pp{pair}{product_size}) {
            die "Problem with $segid amplicon length ".
                "$amplicon_len != $$pp{pair}{product_size}".
                "amplicon seq\n$amplicon_seq\n";
        }
        my ($perc_cov, $bp_cov) = exon_coverage($lchrstart, $rchrstart, 
                         $start, $end, $$p{TARGET_FLANK});
        my $lsnps = snps_in_interval($snphash, $lchrstart, $lchrend);
        my @lsnps = map { $$lsnps{$_}{desc} } sort {$a<=>$b} keys %$lsnps;
        my $rsnps = snps_in_interval($snphash, $rchrstart, $rchrend);
        my @rsnps = map { $$rsnps{$_}{desc} } sort {$a<=>$b} keys %$rsnps;
        my $totsnps = @lsnps + @rsnps;
        $num++;
        my ($ll, $rl) = ('left','right');
        if ($FR eq 'R') { 
            # On reverse strand, left primer is actually right primer 
            # and vice versa.  Amplicon will be in reverse comp direction
            ($ll, $rl) = ('right','left');
            $amplicon_seq = rev_complement($amplicon_seq);
        }
        my %pp = (
            "${ll}_primer_seq"=>$lseq, "${ll}_primer_len"=>$$pp{left}{len},
            "${ll}_primer_start"=>$lchrstart, "${ll}_primer_end"=>$lchrend,
            "${ll}_primer_gc"=>$$pp{left}{gc_percent}, 
            "${ll}_primer_tm"=>$$pp{left}{tm},
            "${ll}_snps"=>join(", ", @lsnps),
            "${rl}_primer_seq"=>$rseq, "${rl}_primer_len"=>$$pp{right}{len},
            "${rl}_primer_start"=>$rchrstart, "${rl}_primer_end"=>$rchrend,
            "${rl}_primer_gc"=>$$pp{right}{gc_percent}, 
            "${rl}_primer_tm"=>$$pp{right}{tm},
            "${rl}_snps"=>join(", ", @rsnps),
            amplicon_seq=>$amplicon_seq,
            amplicon_gc=>$amplicon_gc,
            amplicon_len=>$amplicon_len,
            product_size=>$$pp{pair}{product_size},
            penalty=>$$pp{pair}{penalty},
            template_strand=>$FR, 
            perc_cov=>$perc_cov, bases_cov=>$bp_cov, 
            num_snps=>$totsnps,
        );
        push(@pairs, \%pp);
        if ($$p{DBUG}) {
            warn "leftchrpos $leftchrpos\n";
            warn "seqtemplate:\n$seqtemplate\n";
            warn "  $num)\n".hash_debug(\%pp)."\n";
        }
    }
    (\@pairs);
}

sub get_amplicon {
    my $seqtemplate = shift;
    my $start = shift;
    my $end = shift;

    my $len = $end - $start + 1;
    my $amplicon = substr($seqtemplate, $start, $len);
    ($amplicon);
}

sub exon_coverage { # target plus flanking covered by amplicon
    my $amp_s = shift;
    my $amp_e = shift;
    my $exon_s = shift;
    my $exon_e = shift;
    my $flank = shift;

    my @amp = sort {$a<=>$b} ($amp_s, $amp_e);
    my @ex = sort {$a<=>$b} ($exon_s, $exon_e);
    $ex[0] -= $flank;
    $ex[1] += $flank;
    my $s = $amp[0]<$ex[0] ? $ex[0] : $amp[0];
    my $e = $amp[1]>$ex[1] ? $ex[1] : $amp[1];
    my $totsize = $ex[1]-$ex[0]+1;
    my $bp_cov = $e-$s+1;
    my $perc = sprintf "%4.1f", $bp_cov*100/$totsize;
    ($perc, $bp_cov);
}

#-----------------------------------------------------------------------------

sub sort_and_clean_primerlist {
    my $ppairs = shift;
    my $numpairs = shift;
    # sort by pair penalty; best pair first
    my @sortedpairs = sort { $$a{penalty} <=> $$b{penalty} } @$ppairs;
    my %seen;
    my @goodpairs;
    # remove any pairs that use the same primers
    foreach my $pp (@sortedpairs) {
        my $lseq = $$pp{left_primer_seq};
        my $rseq = $$pp{right_primer_seq};
        my $k = join(":", $lseq, $rseq);
        next if $seen{$k};
        $seen{$k} = $pp;
        push(@goodpairs, $pp);
        last if scalar(@goodpairs)==$numpairs;
    }
    my $best_penalty = @goodpairs ? $goodpairs[0]{penalty} : 99;
    (\@goodpairs, $best_penalty);
}

sub create_exon_results {
    my $filelabel = shift;
    my $primers = shift;
    my $seqhash = shift;
    my $p = shift;

    my $outfile = "$filelabel.primers_found.txt";
    if ($$p{OUTLABEL}) { $outfile = $$p{OUTLABEL}.".$outfile"; }
    if ($$p{OUTDIR}) { $outfile = File::Spec->catfile($$p{OUTDIR}, $outfile) }
    warn "==Writing $outfile\n";
    my $ofh = new FileHandle ">$outfile" or die ">$outfile: $!";
    my @fields = qw/seqlen strand chrom start end/;
    my $fieldstr = join("\t", qw/segment num_primers target_length strand chrom
            target_start target_end/);
    print $ofh $fieldstr."\ttarget_name\n";
    $$p{RES_FIELDS} = join("\t", "gene", "target_name", $fieldstr);
    my @seqids = sort { $$seqhash{$a}{start}<=>$$seqhash{$b}{start} }
                 keys %$primers;
    if ($$seqhash{$seqids[0]}{strand}<0) { @seqids = reverse @seqids }
    foreach my $seqID (@seqids) {
        my @row = map { $$seqhash{$seqID}{$_} } @fields;
        foreach my $segment (sort keys %{ $$primers{$seqID} }) {
            my $numprimers = @{ $$primers{$seqID}{$segment} };
            my $rowstr = join("\t", $segment, $numprimers, @row);
            print $ofh $rowstr."\t$seqID\n";
            push(@{ $$p{RESULTS} }, "$filelabel\t$seqID\t".$rowstr);
        }
    }
    $ofh->close;
}

sub print_primers {
    my $filelabel = shift;
    my $primers = shift;
    my $seqhash = shift;
    my $p = shift;

    my $label = "$filelabel.";
    if ($$p{OUTLABEL}) { $label = "$$p{OUTLABEL}.$label"; }
    if ($$p{OUTDIR}) { $label = File::Spec->catfile($$p{OUTDIR}, $label); }
    my $outfile = $label."all_primers.txt";
    my $fafile = $label."all_primers.fa";
    my $ampfile = $label."all_amplicons.fa";
    warn "Writing $outfile, $fafile, and $ampfile\n";
    my $ofh = new FileHandle ">$outfile" or die ">$outfile: $!";
    my $fafh = new FileHandle ">$fafile" or die ">$fafile: $!";
    my $ampfh = new FileHandle ">$ampfile" or die ">$ampfile: $!";

    my @fields = qw/template_strand penalty product_size perc_cov num_snps
        left_primer_seq left_adaptor_seq left_primer_len left_primer_start 
        left_primer_end left_primer_gc left_primer_tm left_snps
        right_primer_seq right_adaptor_seq right_primer_len right_primer_start 
        right_primer_end right_primer_gc right_primer_tm right_snps
        amplicon_seq amplicon_gc amplicon_len
        /;
    my $fstr = join("\t", @fields);
    $fstr =~ s/num_snps/num_snps (MAF>=$$p{MIN_MAF})/;
    print $ofh "target_name\tsegment\tprimer_num\t$fstr\t".
            join("\t", qw/chrom target_start target_end exon_strand 
            target_length/)."\n";
    foreach my $seqID (sort {$$seqhash{$a}{start}<=>$$seqhash{$b}{start}} 
            keys %$primers) {
        my $chrom = $$seqhash{$seqID}{chrom};
        my $strand = $$seqhash{$seqID}{strand};
        my $start = $$seqhash{$seqID}{start}-$$p{TARGET_FLANK};
        my $end = $$seqhash{$seqID}{end}+$$p{TARGET_FLANK};
        my $seqlen = $$seqhash{$seqID}{seqlen};
        foreach my $segment (sort keys %{ $$primers{$seqID} }) {
            my $id = $seqID; 
            $id .= "_$segment" unless $segment eq '1a';
            my $num = 0;
            foreach my $pp (@{ $$primers{$seqID}{$segment} }) {
                                         # remove restriction site 2015/12/15
                $$pp{left_adaptor_seq} = #"CATCGTGAGTCACTCG". 
                    rev_complement($$pp{left_primer_seq})."gcggccgcCTATAGTGT";
                                         # remove restriction site 2015/12/15
                $$pp{right_adaptor_seq} = #"GTACGAGGTCTCA".
                    $$pp{right_primer_seq}."GGCGCGCCTCCCTTTAG";
                my @row = map { $$pp{$_} } @fields;
                $num++;
                print $ofh join("\t", $seqID, $segment, $num, @row, $chrom,
                        $start, $end, $strand, $seqlen)."\n";
                my $FR = $$pp{template_strand};
                my $commentL = chr_coord_str($$seqhash{$seqID}{chrom}, $FR,
                    $$pp{left_primer_start}, $$pp{left_primer_end}, $strand);
                my $commentR = chr_coord_str($$seqhash{$seqID}{chrom}, $FR,
                    $$pp{right_primer_start},$$pp{right_primer_end}, $strand);
                my $commentAmp;
                if ($FR eq 'R') {
                    $commentAmp = chr_coord_str($$seqhash{$seqID}{chrom}, $FR,
                        $$pp{left_primer_end}, $$pp{right_primer_start}, $strand);
                } else {
                    $commentAmp = chr_coord_str($$seqhash{$seqID}{chrom}, $FR,
                        $$pp{left_primer_start}, $$pp{right_primer_end}, $strand);
                }
                print $fafh ">${id}_$num${FR}_left\t$commentL\n".
                    "$$pp{left_primer_seq}\n";
                print $fafh ">${id}_$num${FR}_right\t$commentR\n".
                    "$$pp{right_primer_seq}\n";
                print $ampfh ">${id}_$num${FR}\t$commentAmp\n".
                    "$$pp{amplicon_seq}\n";
            }
        }
    }
    $ofh->close();
    $fafh->close();
    $ampfh->close();
}

sub chr_coord_str {
    my $chrom = shift;
    my $FR = shift;
    my $start = shift;
    my $end = shift;
    my $strand = shift;

    $strand *= -1 if $FR eq 'R';
    ($start, $end) = sort {$a<=>$b} ($start, $end);
    if ($start > $end) {
        my $tmp = $start;
        $start = $end;
        $end = $tmp;
        $strand = -1*$strand;
    }
    return sprintf "%s:%d:%d:%d", $chrom, $start, $end, $strand;
}

#-----------------------------------------------------------------------------

sub classify_exons {
    my $exonids = shift;
    my $seqhash = shift;

    my @targets;
    my $summary = 
        join("\t", qw/num exon_size lf_flank_size rt_flank_size name/)."\n";
    my $num = 0;
    my $gene = '?';
    foreach my $exonid (@$exonids) {
        $gene = $exonid; $gene =~ s/:.*//;
#        warn "exonid $exonid\n";
        $num++;
        my $seqlen = $$seqhash{$exonid}{seqlen};
        my $previd = $$seqhash{$exonid}{prevseq};
        my $nextid = $$seqhash{$exonid}{nextseq};
        my ($prevseqlen, $nextseqlen) = ('N/A', 'N/A');
        if ($previd) { $prevseqlen = $$seqhash{$previd}{seqlen} }
        if ($nextid) { $nextseqlen = $$seqhash{$nextid}{seqlen}; }
        my @row = ($seqlen, $prevseqlen, $nextseqlen);
        $summary .= "$num)\t".join("\t", @row, $exonid)."\n";
        next unless $$seqhash{$exonid}{need_primers};

        if ($$seqhash{$nextid}{need_primers}) {
            my $newid = merge_seqs($exonid, $nextid, $seqhash);
            push(@targets, $newid);
        } else {
            push(@targets, $exonid);
        }
    }
    my $diff = scalar(@$exonids)-scalar(@targets);
    if ($diff) {
        warn "\nChanges $diff\n";
        return classify_exons(\@targets, $seqhash);
    }
    warn "\nTotal sequence targets: ".scalar(@targets)."\t$gene\n";
    warn $summary;
    return \@targets;
}


sub merge_seqs {
    my $seqid1 = shift;
    my $seqid2 = shift;
    my $seqhash = shift;

    warn "Merging $seqid1 $seqid2\n";
    $$seqhash{$seqid1}{need_primers} = 0;
    $$seqhash{$seqid2}{need_primers} = 0;
    my ($newstart) = sort {$a<=>$b} 
        ($$seqhash{$seqid1}{start}, $$seqhash{$seqid2}{start});
    my ($newend) = sort {$b<=>$a} 
        ($$seqhash{$seqid1}{end}, $$seqhash{$seqid2}{end});
    my $newseqlen = $newend-$newstart+1;
    my $newseq = $$seqhash{$seqid1}{seq}.$$seqhash{$seqid2}{seq};
    my $strand = $$seqhash{$seqid1}{strand};
    if (length($newseq) != $newseqlen) {
        warn "Problem merging $seqid1, $seqid2 ".
            "($newstart, $newend, $strand)\n";
        die "expect seqlen $newseqlen but is ".length($newseq)." bp\n$newseq";
    }
    my @name1 = split(/(ENSE\d+)_/, $$seqhash{$seqid1}{name});
    my @name2 = split(/(ENSE\d+)_/, $$seqhash{$seqid2}{name});
    my @type1 = split(/-/, pop(@name1));
    my @type2 = split(/-/, pop(@name2));
    my %types;
    foreach my $t (@type1, @type2) { $types{$t} = 1; }
    my @ens1 = grep(/ENSE/, (map { split(/-/, $_) } @name1));
    my @ens2 = grep(/ENSE/, (map { split(/-/, $_) } @name2));
    my %ens;
    foreach my $e (@ens1, @ens2) { $ens{$e} = 1; }
    my $newname = join("-", sort keys %ens)."_".join("-", sort keys %types);
    my $newseqid = join(":", $$seqhash{$seqid1}{chrom}, $newstart, $newend,
            $strand, $newname);
    $$seqhash{$newseqid} = { need_primers=>1, 
        prevseq=>$$seqhash{$seqid1}{prevseq},
        nextseq=>$$seqhash{$seqid2}{nextseq},
        desc=>$$seqhash{$seqid1}{desc}, seq=>$newseq, seqlen=>$newseqlen, 
        name=>$newname, chrom=>$$seqhash{$seqid1}{chrom}, 
        start=>$newstart, end=>$newend, strand=>$strand,
    };
    ($newseqid);
}


#-----------------------------------------------------------------------------

sub print_p3input_params {
    my $ofh = shift;
    my $p3params = shift;

    # build Boulder IO string:
    my @fields = reverse sort keys %$p3params;
    my @primer_params = map { my $f = uc $_; $f."=".$$p3params{$_} } @fields;
    my $primer_params = join("\n", @primer_params)."\n";
    print $ofh "$primer_params=\n";

}

sub run_primer3 {
    my $infile = shift;
    my $p = shift;

    warn "== Running primer3\n";
    # run primer3 and save any errors to error file
    my $primer_prog = "$$p{P3EXE} -strict_tags";
    my $errfile = $$p{FILELABEL} ? "$$p{FILELABEL}.err" : "PRIMER3.err";
    if ($$p{OUTDIR}) { $errfile = File::Spec->catfile($$p{OUTDIR}, $errfile) }
    unlink($errfile) if -f $errfile;
    my $cmd = "$primer_prog < $infile 2> $errfile";
    warn "  Running $cmd\n";
    my $ppho = new FileHandle "$cmd |";
    local($/) = undef;
    my $primer_output = <$ppho>;
    $ppho->close();

    if (-s $errfile) { # have errors
        die "  ERROR output from primer3 in $errfile\n";
    } else { unlink($errfile); } # remove empty file

    my $primers;
    if ($$p{FILELABEL}) { # save copy of primer3 output to file
        my $outfile = "$$p{FILELABEL}.p3out.txt";
        if ($$p{P3FILES_SUBDIR}) {
            $outfile = File::Spec->catfile($$p{P3FILES_SUBDIR}, $outfile);
        }
        if ($$p{DBUG}) { warn "Writing $outfile\n"; }
        my $ofh = new FileHandle ">$outfile" or die ">$outfile: $!";
        print $ofh $primer_output;
        $ofh->close;
        $primers = parse_primer3_output($outfile, $p);
    } else {
        $primers = parse_primer3_output($primer_output, $p);
    }
    ($primers);
}


sub parse_primer3_output {
    my $out = shift;
    my $p = shift;
    my $primerhash = shift || {};

    my @recs;
    my $isfile = $out !~ /\n/ && -f $out;
    if ($isfile) { # if input is a file
        $/ = "\n=\n";
        my $ofh = new FileHandle $out or die "$out: $!";
        @recs = $ofh->getlines;
        $ofh->close;
        $/ = "\n";
    } else {
        @recs = split(/\n=\n/, $out);
    }
    if ($recs[$#recs] !~ /=$/) {
        warn "Reading $out\n" if $isfile;
        die "Truncated output from primer3\n$recs[$#recs]\n";
    }

    warn "  Got ".scalar(@recs)." primer record(s)\n";
    foreach my $rec (@recs) {
        my @lines = split(/\n/, $rec);
        my ($seqid) = grep(s/SEQUENCE_ID=//, @lines);
        my ($rtprimer) = grep(s/SEQUENCE_PRIMER_REVCOMP=//, @lines);
        my ($lfprimer) = grep(s/SEQUENCE_PRIMER=//, @lines);
        my ($target) = grep(s/SEQUENCE_TEMPLATE=//, @lines);
        my @primers;
        foreach my $l (@lines) {
            next unless $l =~ s/^PRIMER_(LEFT|RIGHT|PAIR)_(\d+)_?(.*)=//;
            my $type = lc $1;
            my $primer_num = $2;
            my $field = lc $3;
            if ($field) {
                $primers[$primer_num]{$type}{$field} = $l unless
                    $primers[$primer_num]{$type}{$field};
            } else { # start and length of primer
                my ($start, $len) = split(/,/, $l);
                next if $primers[$primer_num]{$type}{len};
                $primers[$primer_num]{$type}{start} = $start;
                $primers[$primer_num]{$type}{len} = $len;
                $primers[$primer_num]{$type}{target} = $target;
            }
        }
        my $key = $seqid;
        $key .= ":rt-$rtprimer" if $rtprimer;
        $key .= ":lf-$lfprimer" if $lfprimer;
        $$primerhash{$key} = \@primers;
        warn "  ".scalar(@primers)." primers for $key\n" if $$p{DBUG};
    }
    ($primerhash);
}

sub gc_perc {
    my $seq = shift;

    my $totlen = length($seq);
    my @gc = $seq =~ /[gcGC]/g;
    my $numgc = scalar(@gc);
    my $gcperc = $numgc*100/$totlen;
    ($gcperc); 
}

sub rev_complement {
    my $seq = shift;

    $seq =~ tr/ATCGatcg/TAGCtagc/;
    $seq = reverse($seq);
    return $seq;
}

##############################################################################

sub hash_debug {
    my $h = shift;
    my $maxlevel = shift || 1;
    my $level = shift || 0;
    my $indent = shift || '  ';

    my @k = sort keys %$h;
    my $str = '';
    foreach my $k (@k) {
        my $v = $$h{$k};
        $str .= sprintf "$indent$k\t%s\n", defined($v) ? $v : '';
        if ($v && $v =~ /HASH/ && $level < $maxlevel) {
            $str .= hash_debug($v, $maxlevel, $level+1, $indent.'  ' );
        }
    }
    ($str);
}

sub enumdate {
    my $long = shift;
    my ($se,$mi,$hr,$day,$mo,$yr) = localtime(time);
    my $enumdate = sprintf "%04d%02d%02d", $yr+1900, $mo+1, $day;
    unless ($long) { $enumdate =~ s/^..//; }
    ($enumdate);
}


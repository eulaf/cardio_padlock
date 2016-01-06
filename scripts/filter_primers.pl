#!/usr/bin/env perl

=head1 NAME

filter_primers.pl -- pick best primer pairs from list

=head1 SYNOPSIS

B<filter_primers.pl> <all_primers.blast_counts.txt file(s)>

B<filter_primers.pl> *.all_primers.blast_counts.txt

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
use Getopt::Long;
require Pod::Usage;
use strict;
use warnings; no warnings 'once';

BEGIN {
    my ($base, $path, $suff) = fileparse(abs_path($0));
    my $up = File::Spec->updir();
    my $lib = File::Spec->catfile($path, $up, 'lib'); 
    my $lib2 = File::Spec->catfile($path, 'lib'); 
    unshift(@INC, 'lib', $lib, $lib2);
}

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

Input is the all_primers.txt files from call_primers.pl.

This script will pick the best primer pair for each template.  The best
is determined as the primer pair with the fewest SNPs and lowest primer3
penalty score.

=head2 Output files

=item YYMMDD.geneID.filtered_primers.txt

This tab-delimited text file describes each primer pair chosen, one pair
per line.  The file includes the primer sequences, amplicon sequence,
location in human genome, number of SNPs in primers, SNP names, and more.

=item YYMMDD.geneID.filtered_amplicons.fa

This fasta file contains the amplicons for selected primer pairs.

=item YYMMDD.geneID.filtered_primers.fa

This fasta file contains the primer sequences for selected primer pairs.

=item YYMMDD.filtered_primers.txt, YYMMDD.filtered_amplicons.fa,
    YYMMDD.filtered_primers.fa

Compilation of all geneID files into one file.

=cut

#############################################################################

=head1 OPTIONS

=over 4

=item B<--od> <dir>

=item B<--outdir> <dir>

Directory for output files.

=item B<--label> <str>

Label for output files.

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
        BLAST0 => 3, # num blast match default when 0
        LABEL => undef,
        FORCE => 0,
        OUTDIR => undef,
        DBUG => 0,
    };

    Getopt::Long::Configure('no_pass_through');
    Getopt::Long::GetOptions(
        'version' => sub { print "$PROGRAM, version $VERSION\n\n";
            Pod::Usage::pod2usage(-verbose => 0,
                                  -exitstatus => 0) },
        'label=s' => \$$p{LABEL},
        'force' => \$$p{FORCE},
        'blast0=i' => \$$p{BLAST0},
        'od|outdir=s' => \$$p{OUTDIR},
        'dbug' => \$$p{DBUG},
    ) || die "\n";
    @ARGV || die "Need all_primers.txt file(s)\n";
    unless ($$p{LABEL}) {
        $$p{LABEL} = basename($ARGV[0]);
        $$p{LABEL} =~ s/\..*//;
    }
    (\@ARGV, $p);
}
################################################################################

### MAIN PROGRAM

my ($infiles, $p) = setup();

my @best = ();
my %notcov = ();
my ($fields, $head);
foreach my $infile (@$infiles) {
    (my $primers, $fields, $head) = parse_all_primers($infile, $p);
    my $pset = get_best_primer_pairs($primers, $p);
    my ($best, $notcov) = check_primer_set($pset, $p);
    push(@best, @$best);
    %notcov = (%$notcov, %notcov);
    print_best_primers($best, $fields, $head, $infile, $notcov, $p);
}
if ($#$infiles) {
    print_best_primers(\@best, $fields, $head, '', \%notcov, $p);
}

#-----------------------------------------------------------------------------

sub create_pkey {
    my $d = shift;
    my $k = sprintf "%5s %11d %s %s", $$d{chrom}, $$d{target_start}, 
                                   $$d{target_name}, $$d{exon_strand};
    return $k;
}

sub pkey_components {
    my $k = shift;
    $k =~ s/^\s+|\s+$//;
    my ($chrom, $start, $name, $strand) = split(/\s+/, $k, 4);
    return ($chrom, $start, $name, $strand);
}

sub parse_all_primers {
    my $file = shift;
    my $p = shift;

    print STDERR "\nParsing $file\n";
    my $fh = new FileHandle $file or die "$file: $!";
    my @lines = $fh->getlines;
    $fh->close;
    chomp @lines;
    my $head = shift @lines;
    my @fields = map { s/ .*//; $_ } split(/\t/, $head);

    my %primers;
    my $num = 0;
    foreach my $l (@lines) {
        my @v = split(/\t/, $l);
        my %d = map { $_=>shift @v } @fields;
        $num++;
        my $k = create_pkey(\%d);
        my $seg = $d{segment};
        my $FR = $d{template_strand};
        my $numleft = $d{num_blast_hits_left} || $$p{BLAST0};
        my $numright = $d{num_blast_hits_right} || $$p{BLAST0};
        my $blasthits = $numleft * $numright;
        # hash enables sorting by fields: num snp, blast hits, penalty
        push(@{ $primers{$k}{$seg}{$d{num_snps}}{$blasthits}{$d{penalty}}{$FR} }, \%d);
    }
    my $numtargets = keys %primers;
    warn "  Got $num primers for $numtargets targets\n";
    (\%primers, \@fields, $head);
}


sub get_best_primer_pairs {
    my $primers = shift;
    my $p = shift;

    my %best;
    warn "Getting best primers by num_snps, blast hits, penalty, FR\n";
    foreach my $k (sort keys %$primers) {
        foreach my $segment (sort keys %{$$primers{$k}}) {
            my @numsnps = sort {$a<=>$b} keys %{$$primers{$k}{$segment}};
            my $leastsnps = shift @numsnps;
            my $h = $$primers{$k}{$segment}{$leastsnps};
            my @blasthits = sort {$a<=>$b} keys %$h;
            my $leastblast = shift @blasthits;
            my $hh = $$primers{$k}{$segment}{$leastsnps}{$leastblast};
            my @penalties = sort {$a<=>$b} keys %$hh;
            my $leastpen = $penalties[0];
            my @strand = sort keys %{$$hh{$leastpen}}; #prefer F to R
            my @primers = sort {$$b{perc_cov}<=>$$a{perc_cov}} 
                          @{ $$hh{$leastpen}{$strand[0]} };
            warn "$k\t$segment\tnum_snp=$leastsnps\tblast_hits=$leastblast\t".
                "penalty=$leastpen\t".
                "strand $strand[0]\t".scalar(@primers)." primers\n";
            $best{$k}{$segment} = $primers[0];
        }
    }
    (\%best);
}

sub amplified_region_coords {
# this is the area covered by the amplicon excluding primers
    my $pp = shift;
    my @coords = sort {$a<=>$b} map { $$pp{$_} } qw/left_primer_start 
        left_primer_end right_primer_start right_primer_end/;
    shift @coords; pop @coords; # get inner coordinates
    my $s = shift @coords;
    my $e = pop @coords;
    ($s+1, $e-1);
}

sub check_primer_set {
    my $pset = shift;
    my $p = shift;

    warn "Checking primers for overlaps\n";
    my @best = ();
    my %intervals_uncovered;
    my $covered_to = 0;
    TARGET: foreach my $k (sort keys %$pset) {
        my @covered = ();
        my %coords;
        my ($s, $e);
        # mark where each primer covers the target
        foreach my $segment (sort keys %{$$pset{$k}}) {
            my $pp = $$pset{$k}{$segment};
            unless (@covered) {
                $s = $$pp{target_start};
                $e = $$pp{target_end};
                # check if part of target already covered by previous 
                # primer set
                if ($covered_to >= $e) {
                    warn "$k already covered!\n";
                    next TARGET;
                } elsif ($covered_to >= $s) {
                    $s = $covered_to+1;
                }
            }
            my ($s1, $e1) = amplified_region_coords($pp);
            $coords{$segment} = [$s1, $e1];
            for(my $i=$s1; $i<=$e1; $i++) {
                next if $i<$s or $i>$e;
                push(@{ $covered[$i-$s] }, $segment);
            }
        }
        # determine if primer covers any of the target uniquely
        my %uniq;
        my @notcovered;
        for(my $loc=$s; $loc<=$e; $loc++) {
            my $i = $loc-$s;
            my $num = $covered[$i] ? scalar(@{$covered[$i]}) : 0;
            if ($num==1) {
                if ($$p{DBUG}) {
                    my $segment = $covered[$i][0];
                    my ($s1, $e1) = @{ $coords{$segment} };
                    warn "$k: $segment is unique at $loc ($s1, $e1)\n";
                }
                $uniq{$covered[$i][0]}++;
            } elsif ($num==0) {
                push(@notcovered, $loc);
                warn "$k base $loc not covered by primers\n";
            }
        }
        my @chosen = sort keys %uniq;
        foreach my $segment (sort keys %{$$pset{$k}}) {
            next if grep(/^$segment$/, @chosen);
            my $patt = "^".join('$|^', @chosen).'$';
            my $pp = $$pset{$k}{$segment};
            my ($s1, $e1) = @{ $coords{$segment} };
            my $uniq = 0;
            for(my $i=$s1; $i<=$e1; $i++) {
                next if $i<$s or $i>$e;
                # num other chosen primers covering this base
                my $num = grep(/$patt/, @{ $covered[$i-$s] });
                $uniq++ if $num==0;
            }
            if ($uniq) { 
                push(@chosen, $segment);
            } else {
                warn "$k Not needed: $segment ($s1, $e1)\n";
            }
        }
        foreach my $segment (sort @chosen) {
            my $pp = $$pset{$k}{$segment};
            my ($s1, $e1) = amplified_region_coords($pp);
            $covered_to = $e1 if $covered_to < $e1;
            push(@best, $pp);
        }
        if (@notcovered) {
            my @intervals_uncovered = ();
            my $tot_uncovered = 0;
            my $s = shift @notcovered;
            my $e = $s;
            foreach my $pos (@notcovered) {
                if ($pos==$e+1) {
                    $e = $pos;
                } else {
                    if ($e>$s) {
                        push(@intervals_uncovered, [$s, $e]);
                        $tot_uncovered += $e - $s + 1;
                    } elsif ($e==$s) {
                        push(@intervals_uncovered, [$s, $s]);
                        $tot_uncovered += 1;
                    }
                    $s = $pos;
                    $e = $pos;
                }
            }
            if ($e>$s) {
                push(@intervals_uncovered, [$s, $e]);
                $tot_uncovered += $e - $s + 1;
            } elsif ($e==$s) {
                push(@intervals_uncovered, [$s, $s]);
                $tot_uncovered += 1;
            }
            for(my $i=0; $i<@intervals_uncovered; $i++) {
                my $loc = $intervals_uncovered[$i];
                my $tot = $$loc[1]-$$loc[0]+1;
                warn "$k bases not covered by primers: @$loc ($tot bp)\n";
                my ($chrom, $start, $target, $strand) = pkey_components($k);
                unless ($intervals_uncovered{$k}) { 
                    $intervals_uncovered{$k} = '' 
                }
                my $num = scalar(@intervals_uncovered)>1 ? sprintf ".%d",$i+1 : '';
                warn "Target $target Num $num INTS ".scalar(@intervals_uncovered)."\n";
                $intervals_uncovered{$k} .= join("\t", $target.$num, $chrom, 
                                    $$loc[0], $$loc[1], $strand, $tot)."\n";
            }
        }
    }
    (\@best, \%intervals_uncovered);
}

sub rev_complement {
    my $seq = shift;

    $seq =~ tr/ATCGatcg/TAGCtagc/;
    $seq = reverse($seq);
    return $seq;
}

##############################################################################

sub print_best_primers {
    my $best = shift;
    my $fields = shift;
    my $head = shift;
    my $infile = shift;
    my $notcov = shift;
    my $p = shift;

    my $label;
    if ($infile) {
        $label = basename($infile); $label =~ s/\.all_primers.*//;
        if ($$p{LABEL}) { $label =~ s/^\d+\.//; $label = "$$p{LABEL}.$label"; }
    } elsif ($$p{LABEL}) { $label = $$p{LABEL} 
    } else { $label = enumdate() }
    if ($$p{OUTDIR}) { $label = File::Spec->catfile($$p{OUTDIR}, $label) }
    my $outfile .= "$label.filtered_primers.txt"; 
    my $ampfile .= "$label.filtered_amplicons.fa";
    my $pfile .= "$label.filtered_primers.fa";
    warn "\nWriting $outfile, $ampfile and $pfile\n";
    my $ofh = new FileHandle ">$outfile" or die ">$outfile: $!";
    my $afh = new FileHandle ">$ampfile" or die ">$ampfile: $!";
    my $pfh = new FileHandle ">$pfile" or die ">$pfile: $!";
    print $ofh "$head\n";
# add FR info to names
    foreach my $pp (@$best) {
        my @row = map { $$pp{$_} } @$fields;
        print $ofh join("\t", @row)."\n";
        my $chrom = $$pp{chrom};
        my $strand = $$pp{exon_strand};
        my $FR = $$pp{template_strand};
        if ($FR eq 'R') { $strand *= -1 }
        my $seg = $$pp{segment} eq '1a' ? '' : "_$$pp{segment}";
        my $name = "$$pp{target_name}${seg}_${FR}";
        my @amploc = $FR eq 'R' ?
            sort {$a<=>$b} ($$pp{left_primer_end}, $$pp{right_primer_start}) :
            sort {$a<=>$b} ($$pp{left_primer_start}, $$pp{right_primer_end});
        my $ampcomment = join(":", $chrom, @amploc, $strand);
        print $afh ">$name\t$ampcomment\n".
            "$$pp{amplicon_seq}\n";
        my @lloc = sort {$a<=>$b} ($$pp{left_primer_start}, 
                $$pp{left_primer_end});
        my @rloc = sort {$a<=>$b} ($$pp{right_primer_start}, 
                $$pp{right_primer_end});
        my $lfcomment = join(":", $chrom, @lloc, $strand);
        my $rtcomment = join(":", $chrom, @rloc, $strand);
        print $pfh ">${name}_left\t$lfcomment\n$$pp{left_primer_seq}\n";
        print $pfh ">${name}_right\t$rtcomment\n$$pp{right_primer_seq}\n";
    }
    $ofh->close;
    $afh->close;
    $pfh->close;
    if ($notcov && %$notcov) {
        my $misfile .= "$label.filtered_missing.txt"; 
        warn "Have uncovered regions. Writing $misfile\n";
        my $mfh = new FileHandle ">$misfile" or die ">$misfile: $!";
        print $mfh "target_name\tchrom\tstart\tend\tstrand".
                   "\tnum_bases_not_covered\n";
        foreach my $k (sort keys %$notcov) {
            print $mfh  $$notcov{$k};
        }
        $mfh->close;
    }
}

sub enumdate {
    my $long = shift;
    my ($se,$mi,$hr,$day,$mo,$yr) = localtime(time);
    my $enumdate = sprintf "%04d%02d%02d", $yr+1900, $mo+1, $day;
    unless ($long) { $enumdate =~ s/^..//; }
    ($enumdate);
}


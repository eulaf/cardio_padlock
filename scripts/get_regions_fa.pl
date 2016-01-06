#!/usr/bin/env perl

=head1 NAME

get_regions_fa.pl -- get sequences from Ensembl for specified regions of the
human genome

=head1 SYNOPSIS

B<get_regions_fa.pl> <location file>

The input file should be a tab-delimited file with columns: name,
chromosome, start, end, strand.

Script fetches the sequence and flanking sequence from Ensembl so must be 
run on a computer with an internet connection

=cut

#############################################################################

=head1 REQUIRES

Cwd,
FileHandle,
File::Basename,
File::Spec,
Getopt::Long,
Pod::Usage,
Bio::EnsEMBL::Registry,

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
    my $lib = abs_path(File::Spec->catfile($path, $up, 'lib')); 
    my $lib2 = File::Spec->catfile($path, 'lib'); 
    # path to ensembl api modules
    foreach my $l ($lib, $lib2, 'lib') {
        next unless -d $l;
        my $bioperl = File::Spec->catfile($l, 'bioperl-live'); 
        my $ensAPI = File::Spec->catfile($l, 'ensembl', 'modules'); 
        my $ensAPIvar = File::Spec->catfile($l, 'ensembl-variation', 'modules'); 
        unshift(@INC, $bioperl, $ensAPI, $ensAPIvar);
    }
}

# This script uses the Ensembl perl API
# http://www.ensembl.org/info/docs/api/index.html
use Bio::EnsEMBL::ApiVersion; # get version number
use Bio::EnsEMBL::Registry;

#############################################################################

### CONSTANTS

my $PROGRAM = $0;
$PROGRAM =~ s;^.*/;;;

my $VERSION = '$Revision: 1 $';
$VERSION =~ s/^\$Revision: //;
$VERSION =~ s/\$$//;
$VERSION = 0 unless $VERSION;
$VERSION = sprintf "%3.1f", $VERSION/10;

#############################################################################

=head1 DESCRIPTION

Ensembl API found at
http://www.ensembl.org/info/docs/api/api_installation.html.

This script uses the Ensembl API to access the Ensembl database to retrieve
regions of human genome sequences.  Regions are specified by chromosome,
start, end, and strand.  Flanking sequence will also be retrieved.

=head2 Input

The input file should be a tab-delimited file specifying the location
of the sequence to retrieve.

=head2 Output files

=item geneID.exon_summary.txt

This tab-delimited text file lists all exons in all transcripts for the 
given geneID.  Each unique gene/transcript/exon is listed, one per line,
along with biotype, start, end, size, chromosome, and strand information
from Ensembl.  The 'class' column lists whether the exon is considered 
unique, short or a duplicate by the merging algorithm.  The last 4 columns 
in the file describe the merged exon that covers the gene/transcript/exon.

=item geneID.gene.fa

This fasta file contains the entire gene sequence for geneID,
including 500 bp of flanking sequence and all introns.  The sequences
alternate between lowercase letters (for non-target sequences
such as flanking and intron sequences) and uppercase letters (for
target sequences).  

=item geneID.exons.fa

This fasta file contains all target sequences for geneID.  These
sequences are a subset of the sequences in geneID.gene.fa (only the
uppercase ones).

=cut

#############################################################################

=head1 OPTIONS

=over 4

=item B<--od>

=item B<--outdir>

Output directory for files.

=item B<--flank> <int>

Length of flanking regions.  (default: 500)

=item B<--sep>

Save each sequence in a separate file.

=item B<--label>

Add this label to output files.

=item B<--maf> <float>

Get SNPs with a MAF >= this value. (default: 0)

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
        RSSUMMARY => [qw/snp chrom pos MAF allele consequence_types/],
        SNPS => 1, # get snps in region of gene
        LABEL => undef,
        OUTDIR => undef,
        FLANK => 500, # bp flanking sequence to include
        MAF => 0, # minimum MAF
        ONEPER => 0,
        DBUG => 0,
    };

    Getopt::Long::Configure('no_pass_through');
    Getopt::Long::GetOptions(
        'version' => sub { print "$PROGRAM, version $VERSION\n\n";
            Pod::Usage::pod2usage(-verbose => 0,
                                  -exitstatus => 0) },
        'label=s' => \$$p{LABEL},
        'outdir|od=s' => \$$p{OUTDIR},
        'sep' => \$$p{ONEPER},
        'maf=f' => \$$p{MAF},
        'flank=i' => \$$p{FLANK},
        'dbug' => \$$p{DBUG},
    ) || die "\n";
    @ARGV || die "Need location file\n";
    ($ARGV[0], $p);
}
################################################################################

### MAIN PROGRAM

my ($infile, $p) = setup();
my $locs = parse_locations($infile, $p);
my $ofh;
my $sfh;
unless ($$p{ONEPER}) {
    my $outfile = basename($infile); $outfile =~ s/\..*$//;
    $ofh = get_filehandle($outfile.".sequence.fa", $p);
    $sfh = get_filehandle($outfile.".snps.txt", $p);
    warn "Writing $outfile.sequence.fa, $outfile.snps.txt\n";
    print $sfh join("\t", @{ $$p{RSSUMMARY} });
}

printf STDERR "\nAPI version %s (%s)\n", software_version(),
    'http://www.ensembl.org/info/docs/api/api_installation.html';
warn "Loading registry\n";
my $registry = registry();
warn "Getting slice adaptor (human, core, slice)\n";
my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
warn "  Getting vf adaptor\n";
my $vf_adaptor = $registry->get_adaptor('human', 'variation',
                                        'variationfeature');
my $numlocs = scalar(@$locs);
warn "\nProcessing $numlocs regions\n";
my $num = 0;
my $numfound = 0;
my @notfound = ();
# process regions
foreach my $loc (@$locs) {
    my $name = $$loc{name};
    my $chrom = $$loc{chrom};
    $chrom =~ s/chr//;
    my $strand = $$loc{strand};
    if ($strand =~ /^(\+|1|F)$/) { $strand = 1 }
    elsif ($strand =~ /^(-|-1|R)$/) { $strand = -1 }
    my ($start, $end) = sort {$a<=>$b} ($$loc{start}, $$loc{end});
    my $s = $start - $$p{FLANK}; if ($s<1) { $s=1 }
    my $e = $end + $$p{FLANK};
    my $region = "$chrom:$start:$end:$strand";
    my $seqreg = "$chrom:$s:$e:$strand";
    printf STDERR "\nRegion %2d) $name $region, with flank $seqreg\n", ++$num;
    my $slice = $slice_adaptor->fetch_by_region('chromosome', $chrom, $s, $e);
    if ($slice) { # slice coordinates start at 1 (offset by $s+1)
        if ($$p{ONEPER}) {
            my $outfile = $name; $outfile =~ s/[^-\w_]/_/g;
            $ofh = get_filehandle($outfile.".sequence.fa", $p);
            $sfh = get_filehandle($outfile.".snps.txt", $p);
            print $sfh join("\t", @{ $$p{RSSUMMARY} });
        }
        my $flankL = lc $slice->subseq(1, $start-$s);
        my $target = uc $slice->subseq($start-$s+1, $end-$s+1);
        my $flankR = lc $slice->subseq($end-$s+2, $e-$s+1);
        if ($strand>0) {
            printf $ofh ">${name}:flank_left\tchr$chrom:%d:%d:%d\n$flankL\n", 
                $s, $start-1, $strand;
            printf $ofh ">${name}\tchr$chrom:%d:%d:%d\n$target\n", $start, 
                $end, $strand;
            printf $ofh ">${name}:flank_right\tchr$chrom:%d:%d:%d\n$flankR\n", 
                $end+1, $e, $strand;
        } else { #reverse
            printf $ofh ">${name}:flank_left\tchr$chrom:%d:%d:%d\n%s\n", 
                $end+1, $e, $strand, rev_complement($flankR);
            printf $ofh ">${name}\tchr$chrom:%d:%d:%d\n%s\n", $start, $end, 
                $strand, rev_complement($target);
            printf $ofh ">${name}:flank_right\tchr$chrom:%d:%d:%d\n%s\n", $s, 
                $start-1, $strand, rev_complement($flankL);
        }
        # get snps
        my $vf = $vf_adaptor->fetch_all_by_Slice($slice);
        my $numsnps = @$vf;
        my %snps;
        warn "    Finding SNPs at position: $s - $e (flank $$p{FLANK})\n";
        warn "    Got ".scalar(@$vf)." snps\n" if $$p{DBUG};
        foreach my $snp (@$vf) {
            my $snp_s = $snp->start + $s - 1;
            my $snpinfo = { 
                snp=>$snp->variation_name(), 
                chrom=>'chr'.$snp->seq_region_name, 
                allele =>$snp->allele_string,
                consequence_types => join(",",@{$snp->consequence_type()}),
                MAF=>$snp->minor_allele_frequency, 
                pos=>$snp_s }; 
            my $maf = $$snpinfo{MAF} || 0;
            next if $maf < $$p{MAF};
            if ($snps{$snp_s}) {
                $$snpinfo{snp} = $snps{$snp_s}{snp}.";".$$snpinfo{snp};
                my $allele = sprintf "%s;%s", $snps{$snp_s}{allele} || '',
                                              $$snpinfo{allele}||'';
                my $maf = sprintf "%s;%s", map {defined($_) ? $_ : ''} 
                    ($snps{$snp_s}{MAF}, $$snpinfo{MAF});
                if ($$p{DBUG}) {
                    warn "    Multiple snps at location $snp_s: $$snpinfo{snp}\n";
                    warn "      allele $allele MAF $maf\n";
                }
            }
            $snps{$snp_s} = $snpinfo;
        }
        foreach my $chrpos (sort {$a<=>$b} keys %snps) {
            my $seqpos = $chrpos - $s + 1;
            my @row = map { defined($snps{$chrpos}{$_}) ?
                $snps{$chrpos}{$_} : '' } @{ $$p{RSSUMMARY} };
            print $sfh join("\t", @row)."\n";
        }
        warn "  Total of  ".scalar(keys %snps)." snps\n";
    } else {
        push(@notfound, $name);
        warn "    Data not found for $name\n";
        next;
    }
    if ($$p{ONEPER}) {
        $ofh->close if $ofh;
        $sfh->close if $sfh;
    }
}

$ofh->close if $ofh;
$sfh->close if $sfh;
print "\nDone\n";

##############################################################################

sub parse_locations {
    my $file = shift;
    my $p = shift;

    warn "Reading location file $file\n";
    my $fh = new FileHandle $file or die "$file: $!";
    my @lines = grep(/\w/, map { split(/[\r\n]+/, $_) } $fh->getlines);
    chomp @lines;
    $fh->close;
#    warn "num lines ".scalar(@lines)."\n";
    my @fields = qw/name chrom start end strand/;
    my %fields_left = map { $_=>1 } @fields;
    if ($lines[0] =~ /start/i) {
        my $head = shift @lines;
        my @f = split(/\t/, $head);
        for(my $i=0; $i<@f; $i++) {
            if ($f[$i] =~ /^chr/i && $fields_left{chrom}) {
                $f[$i] = 'chrom';
                $fields_left{chrom} = 0;
            } elsif ($f[$i] =~ /start/i && $fields_left{start}) {
                $f[$i] = 'start';
                $fields_left{start} = 0;
            } elsif ($f[$i] =~ /end/i && $fields_left{end}) {
                $f[$i] = 'end';
                $fields_left{end} = 0;
            } elsif ($f[$i] =~ /strand|sense/i && $fields_left{strand}) {
                $f[$i] = 'strand';
                $fields_left{strand} = 0;
            } elsif (($f[$i] =~ /name/i or $i==0) && $fields_left{name}) {
                $f[$i] = 'name';
                $fields_left{name} = 0;
            }
        }
        @fields = @f;
    }
    my @notfound = grep($fields_left{$_}, @fields);
    if (@notfound) {
        die "Could not find fields (".join(", ", @notfound).") in $file\n";
    }
    my @data;
    if ($$p{DBUG}) { warn join("\t", @fields)."\n" }
    foreach my $l (@lines) {
        my @vals = split(/\s*\t\s*/, $l);
        my %d = map { $_=>shift @vals } @fields;
        push(@data, \%d);
        if ($$p{DBUG}) {
            my @row = map { $d{$_} } @fields;
            warn "Region: ".join("\t", @row)."\n";
        }
    }
    my $num = scalar(@data);
    warn "  Got $num regions\n";
    (\@data);
}


sub registry {
    my $registry = 'Bio::EnsEMBL::Registry';
    $registry->load_registry_from_db(
        -host => 'ensembldb.ensembl.org',
        -user => 'anonymous');
    return $registry;
}

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

sub rev_complement {
    my $seq = shift;
    my @pieces = split(/\n/, $seq, 3);
    my $seqidx = ($pieces[0] =~ /^>/) ? 1 : 0;
    $pieces[$seqidx] =~ tr/ATCGatcg/TAGCtagc/;
    $pieces[$seqidx] = reverse($pieces[$seqidx]);
    return join("\n", @pieces);
}

sub get_filehandle {
    my $filename = shift;
    my $p = shift;

    my $outfile = $$p{LABEL} ? "$$p{LABEL}.$filename" : $filename;
    if ($$p{OUTDIR}) { $outfile = File::Spec->catfile($$p{OUTDIR}, $outfile) }
    warn "  Writing $outfile\n";
    my $ofh = new FileHandle ">$outfile" or die ">$outfile: $!";
    return $ofh;
}


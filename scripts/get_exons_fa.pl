#!/usr/bin/env perl

=head1 NAME

get_exons_fa.pl -- get gene/exon sequences from Ensembl

=head1 SYNOPSIS

B<get_exons_fa.pl> <gene(s) or gene file>

Input is list of genes or a text file with genes listed one per line.

Script fetches gene information from Ensembl and must be run on a 
computer with an internet connection

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
human gene sequences.  The script merges all coding and processed transcripts
into a set of coordinates that encompasses all transcripts.

The sequences include exons and introns and 500 bp (default setting) of
flanking sequence.  Each gene sequence will have the intron sequences,
exon sequences, and flanking sequences as separate lines in the file.  It can
be assumed that the first line after the header is the 5' flanking sequence,
second line is the first exon, third line is first intron, etc.  The
intron and flanking sequences are lowercase while utr/exon sequences are
in uppercase.

=head2 Input

Input can be a list of genes or Ensembl geneIDs on the command line or
in a text file with one gene/Ensembl geneID per line.

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

=item B<--label>

Add this label to output files.

=item B<--canonical>

Only use canonical transcript.

=item B<--coding>

Only use protein coding transcripts.

=item B<--maf> <float>

Get SNPs with a MAF >= this value. (default: 0)

=item B<--rmutrs>

Remove UTRs from exons.

=item B<--dbug>

Print debugging messages.

=back

=cut

#=item B<--snpfile> <file>

#Limit SNPs to those in this file.

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
        MAKE_SEQS => 1, #create fasta files
        SNPS => 1, # get snps in region of gene
        EXONSUMMARY => 1, # create exon summary file
        RSSUMMARY => 1, # create snp file
        LABEL => undef,
        OUTDIR => undef,
        SNPFILE => undef, # restrict snps to those in this file
        FLANK => 500, # bp flanking sequence to include
        RM_UTRS => 0, # do not include utrs with exon seqs
        CODING => 0, # only use protein coding transcripts
        CANONICAL => 0, # only use canonical transcript
        MAF => 0, # minimum MAF
        DBUG => 0,
    };

    Getopt::Long::Configure('no_pass_through');
    Getopt::Long::GetOptions(
        'version' => sub { print "$PROGRAM, version $VERSION\n\n";
            Pod::Usage::pod2usage(-verbose => 0,
                                  -exitstatus => 0) },
        'canonical' => \$$p{CANONICAL},
        'coding' => \$$p{CODING},
        'rmutrs' => \$$p{RM_UTRS},
        'exonsummary!' => \$$p{EXONSUMMARY},
        'rssummary!' => \$$p{RSSUMMARY},
        'label=s' => \$$p{LABEL},
        'outdir|od=s' => \$$p{OUTDIR},
        'snpfile=s' => \$$p{SNPFILE},
        'flank=i' => \$$p{FLANK},
        'maf=f' => \$$p{MAF},
        'dbug' => \$$p{DBUG},
    ) || die "\n";
    @ARGV || die "Need argument\n";
    my @genes;
    if (-f $ARGV[0]) { 
        warn "\nReading genes from file $ARGV[0]\n";
        # one gene per line, first column of file
        my $fh = new FileHandle $ARGV[0] or die "$ARGV[0]: $!";
        my @lines = grep(/\w/, map { split(/[\r\n]+/, $_) } $fh->getlines);
        $fh->close;
        chomp @lines;
        my %seen; #remove any duplicates
        foreach my $l (@lines) {
            my ($gene) = split(/\s/, $l, 2);
            next unless $gene;
            next if $gene =~ /^\s*(gene|name)\s*$/i; 
            push(@genes, $gene) unless $seen{$gene};
            $seen{$gene} = 1;
        }
    } else { @genes = @ARGV }
    if ($$p{EXONSUMMARY}) {
        $$p{EXONSUMMARY} = [qw/gene transcript exon exon_biotype 
            exonStart exonEnd exon_size chromosome strand class 
            target_start target_end ENSEMBL_name nickname/];
    }
    if ($$p{RSSUMMARY}) {
        $$p{RSSUMMARY} = [qw/snp chrom pos MAF allele consequence_types/];
    }
    warn "Num genes: ".scalar(@genes)."\n";
#    exit if $$p{DBUG};
    (\@genes, $p);
}
################################################################################

### MAIN PROGRAM

my ($genes, $p) = setup();
my $snphash = $$p{SNPFILE} ? read_snp_file($$p{SNPFILE}) : undef;

printf STDERR "\nAPI version %s (%s)\n", software_version(),
    'http://www.ensembl.org/info/docs/api/api_installation.html';
warn "Loading registry\n";
my $registry = registry();
warn "Getting gene adaptor (human, core, gene)\n";
my $gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );

my $numgenes = scalar(@$genes);
warn "\nProcessing $numgenes genes\n";
my $num = 0;
my $numfound = 0;
my @notfound = ();
# process genes: find exon coordinates and get sequences
foreach my $genename (@$genes) {
    printf STDERR "\nGene %2d) $genename\n", ++$num;
    my ($geneslice, $data) = get_gene_from_ensembl($genename, $gene_adaptor, 
            $registry, $snphash, $p);
    if ($geneslice) {
        $numfound++;
    } else {
        push(@notfound, $genename);
        warn "    Data not found for $genename\n";
        next;
    }
    create_sequences($genename, $geneslice, $data, $p);
}
if ($snphash) {
    my $numused = grep($$snphash{$_}{found}, keys %$snphash);
    my $tot = scalar(keys %$snphash);
    warn "Total snps used: $numused\n";
    warn "Total snps: $tot\n";
}
warn "\nFound $numfound/$numgenes genes\n";
if (@notfound) {
    my $numnotfound = scalar(@notfound);
    warn "Genes not found ($numnotfound):\n  ".join("\n  ", @notfound)."\n";
}

print "\nDone\n";

##############################################################################

sub read_snp_file {
    my $file = shift;

    warn "Reading snp file $file\n";
    my $fh = new FileHandle $file or die "$file: $!";
    my @lines = grep(/\srs\d+/, map { split(/[\r\n]+/, $_) } $fh->getlines);
    $fh->close;
#    warn "num lines ".scalar(@lines)."\n";
    my %snphash = map { $_ =~ s/.*?\srs/rs/; my ($s) = split(/\s+/, $_,2); 
                     $s =~ s/_[\w]$//; ($s => {line=>$_}) } @lines;
    warn "Parsed ".scalar(keys %snphash)." snps from file\n";
#    warn join("\n", sort keys %snphash)."\n";
    (\%snphash);
}

sub registry {
    my $registry = 'Bio::EnsEMBL::Registry';
    $registry->load_registry_from_db(
        -host => 'ensembldb.ensembl.org',
        -user => 'anonymous');
    return $registry;
}

sub get_gene_from_ensembl {
    my $gene_input = shift;
    my $gene_adaptor = shift;
    my $registry = shift;
    my $snps = shift;
    my $p = shift;

    # Query Ensembl to find gene by genename
    my ($geneid, $geneENS, $genename);
    if ($gene_input =~ /^ENSG\d+/) {
        $geneENS = $gene_adaptor->fetch_by_stable_id($gene_input);
    } else {
        $geneENS = $gene_adaptor->fetch_by_display_label(lc $gene_input);
    }
    if ($geneENS) { 
        my $gsumm = feature2string($geneENS);
        $geneid = $geneENS->stable_id();
        $genename = $geneENS->display_xref->display_id();
        warn "  Found $gsumm\n";
    } else { 
        warn "  Gene $gene_input - Not found.\n";
        return;
    }
    
    # Get exon locs of all transcripts
    my $transcripts = $geneENS->get_all_Transcripts();
    my $numtrans = @$transcripts;
    warn "  $numtrans transcripts\t$gene_input\n";
    next unless $numtrans;
    my %transcripts;
    while ( my $transcript = shift @{$transcripts} ) {
        my $stable_id  = $transcript->stable_id();
        print STDERR "    $genename transcript\t$stable_id; ";
        warn hash_debug($transcript) if $$p{DBUG};
        if ($$p{CANONICAL}) { 
            warn "    CANONICAL\n" if $$p{DBUG};
            next unless $transcript->is_canonical; 
        }
        # Get exon locations minus UTRs
        my $type = $transcript->biotype;# eq 'protein_coding' ? 'exon' : 'transcript';
        if ($$p{CODING}) { next unless $type eq 'protein_coding'; }
        my $data = get_exon_locs($genename, $geneid, $transcript, $type, $p);
        $transcripts{$stable_id} = $data;
    }
    return unless %transcripts;
    my $summfh = *STDERR;
    if ($$p{EXONSUMMARY}) { # write exon summary to file
        $summfh = get_filehandle($gene_input.'.exon_summary.txt', $p);
        print $summfh join("\t", @{ $$p{EXONSUMMARY} })."\n";
    }
    my $merged_locs = merge_locs($summfh, \%transcripts);
    $merged_locs = merge_locs($summfh, \%transcripts, 'utrs', $merged_locs);
    $summfh->close();

    # get SNPs for this gene if requested
    $$merged_locs{snps} = get_snps($geneENS, $registry, $merged_locs, 
                                      $snps, $p) if $$p{SNPS};

    ($geneENS->slice, $merged_locs);
}

sub exondebug {
    my $exonlist = shift;
    my $type = shift || '';

    my $num = 0;
    my $str = '';
    foreach my $e (@$exonlist) {#dbug
        $str .= sprintf "%d) %7d-%7d\t%s\t$type\n", ++$num, 
            $$e{start}, $$e{end}, $$e{name};
    }
    ($str);
}

sub merge_types {
    my %seen;
    my @types;
    foreach my $type (@_) {
        my @t = split(/-/, $type);
        foreach my $t (@t) { 
            push(@types, $t) unless $seen{$t};
            $seen{$t} = 1 
        }
    }
    return join("-", @types);
}

sub merge_locs {
    my $ofh = shift;
    my $transcripts = shift;
    my $type = shift || 'exons';
    my $merged = shift || {};

    my @superlist; #list encompassing all variants 
    my @tsid  = sort keys %$transcripts;
    my $numtrans = scalar(@tsid);
    warn "  Merging $type in $numtrans transcripts\n";
    my $tnum = 0;
    foreach my $tsid (@tsid) {
        $tnum++;
        my @exons = $$transcripts{$tsid}{plusstrand} ? 
            @{ $$transcripts{$tsid}{$type} } : 
            reverse @{ $$transcripts{$tsid}{$type} };
        warn "  --Start    Transcript $tnum/$numtrans) $tsid\n".exondebug(\@exons, $type) if $$p{DBUG};
        if (@superlist) {
            foreach my $ex (@exons) { #see where this exon belongs
                print STDERR "EX$$ex{start}-$$ex{end}:\t" if $$p{DBUG};
                my $found = 0;
                for(my $i=0; $i<@superlist; $i++) {# in the superexon
                    my $se = $superlist[$i];
                    if ($$ex{end} < $$se{start}) {#add
                        splice(@superlist, $i, 0, {%$ex});
        warn "    Adding     $$ex{name} ($$ex{end}<$$se{start}) $i\n".exondebug(\@superlist, $type) if $$p{DBUG};
                        $found=1;
                    } elsif ($$ex{end} < $$se{end}) {
                        if ($$ex{start} < $$se{start}) { #merge
        warn "    Merging $$se{name} with $$ex{name} ($$ex{start}<$$se{start})\n".exondebug(\@superlist, $type) if $$p{DBUG};
                            $$se{start} = $$ex{start};
                            $$se{size} = $$se{end}-$$se{start}+1;
                            $$se{name} = $$ex{name}.'-'.$$se{name};
                            $$se{type} = merge_types($$ex{type}, $$se{type});
                        } else { 
                            $$se{type} = merge_types($$ex{type}, $$se{type});
                        }
                        $found=1;
                    } elsif ($$ex{end} == $$se{end}) {
                        if ($$ex{start} < $$se{start}) { #replace
                            my $type = merge_types($$ex{type}, $$se{type});
        warn "    Replacing $superlist[$i]{name} with $$ex{name} ($$ex{start}<$$se{start}) $i\n".exondebug(\@superlist, $type) if $$p{DBUG};
                            $superlist[$i] = {%$ex};
                            $superlist[$i]{type} = $type;
                        } else { 
                            if ($$ex{start} == $$se{start}) {
                                $$se{type} = merge_types($$ex{type}, $$se{type});
                            } elsif ($$ex{start} > $$se{start}) { # short
                            }
                        }
                        $found=1;
                    } else { # $$ex{end} > $$se{end}
                        if ($$ex{start} <= $$se{start}) { #replace
        warn "    Replace $superlist[$i]{name} with $$ex{name} (end $$ex{end}>$$se{end}; start $$ex{start}<=$$se{start})\n".exondebug(\@superlist, $type) if $$p{DBUG};
                            my $type = merge_types($$ex{type}, $$se{type});
                            $superlist[$i] = {%$ex};
                            $superlist[$i]{type} = $type;
                            $found=1;
                        } elsif ($$ex{start} < $$se{end}) { #merge
        warn "    Merge $superlist[$i]{name}and $$ex{name} ($$ex{start}<$$se{end})\n".exondebug(\@superlist, $type) if $$p{DBUG};
                            $$se{end} = $$ex{end};
                            $$se{size} = $$se{end}-$$se{start}+1;
                            $$se{name} = $$se{name}.'-'.$$ex{name};
                            $$se{type} = merge_types($$ex{type}, $$se{type});
                            $found=1;
                        }
                    }
                    last if $found;
                }
                unless ($found) {
        warn "    Add to end $$ex{name}\n".exondebug(\@superlist, $type) if $$p{DBUG};
                    push(@superlist, $ex);
                }
            }
        } else {
            @superlist = map { {%$_} } @exons;
        }
        warn "  Done       $tsid\n".exondebug(\@superlist, $type).
            "  =end\n" if $$p{DBUG};
    }
    # check for units that should be merged together
    @superlist = sort { $$a{start} <=> $$b{start} } @superlist;
    my $totexons = @superlist;
    my $label = $type; $label =~ s/s$//;
    warn "  Cleaning merged $label list\n";
    print STDERR exondebug(\@superlist, $type)."\n";
    for(my $loop=0; $loop<49; $loop++) {
        printf STDERR "    LOOP %d tot$type=$totexons\n", $loop+1;
        for(my $i=$#superlist; $i>0; $i--) {
            my $j = $i-1;
            my $se1 = $superlist[$j];
            my $se2 = $superlist[$i];
            if ($$se1{end}+1>=$$se2{start}) {
                printf STDERR "      $label %d/%d %s (%d,%d)\t%s (%d,%d)\n",
                    $j+1, $i+1, $$se1{name}, $$se1{start}, $$se1{end},
                    $$se2{name}, $$se2{start}, $$se2{end} if $$p{DBUG};
                if ($$se1{end}>= $$se2{end}) { #replace
                    printf STDERR "    --delete $label%d\n", $i+1;
                } else { #merge
                    $$se1{end} = $$se2{end};
                    $$se1{name} = merge_types($$se1{name}, $$se2{name});
                    $$se1{type} = merge_types($$se1{type}, $$se2{type});
                    warn "      --merge $$se1{start}, $$se1{end} : $$se1{name}\n";
                }
                splice(@superlist, $i, 1);
            }
        }
        my $exontot = scalar(@superlist);
        printf STDERR "    Done loop %d num$type=%d\n", $loop+1, $exontot;
        print STDERR exondebug(\@superlist, $type)."\n";
        last if $totexons == $exontot;
        $totexons = $exontot;
    }
    unless (%$merged) {
        %$merged = map { $_=>$$transcripts{$tsid[0]}{$_} } qw/chrom strand
            plusstrand genename geneid/;
    }
    warn "  Strand $$merged{strand}\n" if $$p{DBUG};
    $$merged{$type} = $$merged{plusstrand} ? 
        [ sort { $$a{start} <=> $$b{start} } @superlist ] :
        [ sort { $$b{start} <=> $$a{start} } @superlist ];
    if ($type eq 'exons') {
        $$merged{start} = $superlist[0]{start};
        $$merged{end} = $superlist[$#superlist]{end};
        warn "  ".scalar(@{ $$merged{$type} })." merged exons\n";
    }
    add_exon_numbering(\@superlist, $$merged{plusstrand}, $$merged{genename});
    if ($$p{DBUG}) { exondebug($$merged{$type}, $type); }
    warn hash_debug($merged) if $$p{DBUG};
    exon_summary($p, $ofh, $transcripts, $merged) if $type eq 'exons';
    ($merged);
}

# $$transcripts{tsid} = {
#    chrom=>$chr, transcript=>$transcript, cds_size=>length($cds),
#    start=>$cds_locs[0], end=>$cds_locs[$#cds_locs],
#    stable_id=>$transcript->stable_id(), strand=>$sstr, 
#    genename=>$genename, geneid=>$geneid, plusstrand=>$strand,
#    exons=>\@exons, utrs=>[], utr5size=>$size5, utr3size=>$size3 };
sub exon_summary {
    my $p = shift;
    my $ofh = shift;
    my $transcripts = shift;
    my $merged = shift;

    my $geneid = $$merged{geneid} || '';
    my $strand = $$merged{strand};
    my $chrom = $$merged{chrom};
    my $superexons = $$merged{exons};
    my @fields = qw/gene_id transcript_id exon_id exon_start exon_end 
        exon_size chromosome strand class superexon_start superexon_end 
        superexon_id/;

    my @tsid  = sort keys %$transcripts;
    warn "  Creating exon summary from ".scalar(@tsid)." transcripts\n";
    foreach my $tsid (@tsid) {
        my @row = ($$transcripts{$tsid}{geneid}, $tsid);
        my @exons = $$transcripts{$tsid}{plusstrand} ? 
            @{ $$transcripts{$tsid}{exons} } : 
            reverse @{ $$transcripts{$tsid}{exons} };
        foreach my $exon (@exons) {
            my @info = map { $$exon{$_} } qw/name type_num start end size/;
            my @m = grep($$_{start}<=$$exon{start} &&
                         $$_{end}>=$$exon{end}, @$superexons);
            my $class = '';
            if (@m>1) {
                warn "$tsid exon\n".exondebug([$exon], 'has too many matches')."\n";
                die "$tsid matches\n".exondebug(\@m, "match to $$exon{name}")."\n";
            } elsif (!@m) { 
                warn "merged\n".exondebug($superexons, 'no match')."\n";
                die "$tsid exon \n".exondebug([$exon], 'no match');
            } else { # one match
                my $m = $m[0];
                if ($$exon{start}==$$m{start} && $$exon{end}==$$m{end}) {
                    # matching coordinates 
                    $class = $$exon{name} eq $$m{name} ? 'unique' :
                        'duplicate';
                } elsif ($$m{name} =~ /\b$$exon{name}\b/) {
                    $class = 'merged';
                } else {
                    $class = 'short';
                }
                warn "$tsid exon \n".exondebug([$exon], $class) if $$p{DBUG};
                warn "match\n".exondebug(\@m, "match $class")."\n" if
                    $$p{DBUG};
                push(@info, $chrom, $strand, $class, $$m{start}, $$m{end},
                        $$m{name}, $$m{nickname});
            }
            print $ofh join("\t", @row, @info)."\n";
        }
    }
}


sub get_exon_locs {
    my $genename = shift;
    my $geneid = shift;
    my $transcript = shift;
    my $type = shift;
    my $p = shift;

    # find UTR sizes so we can exclude them
    my $utr5 = $transcript->five_prime_utr();
    my $utr3 = $transcript->three_prime_utr();
    my $size5 = $utr5 ? length($utr5->seq) : 0;
    my $size3 = $utr3 ? length($utr3->seq) : 0;

    # determine all exon and intron sizes
    my $strand = $transcript->strand() < 0 ? 0 : 1;
    my $slice = $transcript->slice;
    my $chr = $slice->seq_region_name();
    my $sstr = $strand ? '+' : '-';
    my $totsize = 0;
    my @exons;
    my $numexons = my @exonlist = @{ $transcript->get_all_Exons() };
    my @cds_locs;
    foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
        if ($$p{DBUG}) {
            warn "** Exon \n";
            warn hash_debug($exon);
            warn "** done\n";
        }
        my $name = $exon->stable_id();
        my $start = $exon->start();
        my $end = $exon->end();
        my $size = $end-$start+1;
        $totsize += $size;
        push(@cds_locs, $start, $end);
        my %exoninfo = ( start=>$start, end=>$end, size=>$size, 
                         name=>$name, type=>$type );
        push(@exons, \%exoninfo);
        printf STDERR "\t    Exon%3d   ($sstr)\t$start-$end\t%8d bp\n",
               scalar(@exons), $size if $$p{DBUG};
    }
    @cds_locs = sort {$a<=>$b} @cds_locs;
    my $cdna = $transcript->spliced_seq();
    my $cds = $transcript->translateable_seq();
    my %info = (chrom=>$chr, transcript=>$transcript, cds_size=>length($cds),
                start=>$cds_locs[0], end=>$cds_locs[$#cds_locs],
                stable_id=>$transcript->stable_id(),
                strand=>$sstr, plusstrand=>$strand,
                genename=>$genename, geneid=>$geneid, 
                exons=>\@exons, utrs=>[],
                utr5size=>$size5, utr3size=>$size3);
    trim_utrs(\%info, $totsize, $transcript, $strand, $p);
    add_exon_numbering(\@exons, $strand);
    print STDERR "\tNum exons $numexons\tCDS size $info{cds_size}\n";
    if ($$p{DBUG}) {
        printf STDERR "\t$genename cDNA: %d bp\n", length($cdna);
        printf STDERR "\t$genename CDS: %d bp\n", $info{cds_size};
        printf STDERR "\t$genename 5' UTR: $size5 bp\n";
        printf STDERR "\t$genename 3' UTR: $size3 bp\n";
    }
    (\%info);
}

# add numbering to exon, utr or processed_transcript
sub add_exon_numbering {
    my $exons = shift;
    my $plusstrand = shift;
    my $genename = shift;

    my %count;
    my @exons = $plusstrand ? @$exons : reverse @$exons;
    foreach my $ex (@exons) { 
        my @t = split(/-/, $$ex{type});
        if (grep(/protein/, @t)) {
            $count{protein_coding}++;
            $$ex{num} = $count{protein_coding}; 
        } elsif (my @u = grep(/utr\d+/, @t)) {
            $u[0] =~ /(utr\d+)/;
            $count{$1}++;
            $$ex{num} = $count{$1}; 
        } else {
            $count{$t[0]}++;
            $$ex{num} = $count{$t[0]}; 
        }
        $$ex{type_num} = "$$ex{type}_$$ex{num}"; 
    }
    # add nicknames
    if ($genename) {
        foreach my $ex (@exons) {
            if ($$ex{type} =~ /protein_coding/) {
                $$ex{nickname} = "$genename:exon_$$ex{num}"
            } elsif ($$ex{type} =~ /(utr\d+)/) {
                my $maxnum = $count{$1};
                $$ex{nickname} = "$genename:".$1;
                if ($maxnum>1) { $$ex{nickname}.= "_$$ex{num}" }
            } else {
                my @t = split(/-/, $$ex{type});
                my $maxnum = $count{$t[0]};
                $$ex{nickname} = "$genename:$$ex{type}";
                if ($maxnum>1) { $$ex{nickname}.= "_$$ex{num}" }
            }
        }
    }
}

sub trim_utrs {
    my $info = shift;
    my $totexonsize = shift;
    my $transcript = shift;
    my $strand = shift;
    my $p = shift;

    my $size5 = $$info{utr5size};
    my $size3 = $$info{utr3size};
    my $totsize = 0;
    my $cds_end = $totexonsize-$size3;
    my @exons;
    my @utrs;
    foreach my $exon (@{ $$info{exons} }) {
        $totsize += $$exon{size};
        if ($totsize<=$size5) {#this exon is in 5'UTR
            $$exon{type} = 'utr5';
            $$p{RM_UTRS} ? push(@utrs, $exon) : push(@exons, $exon);
        } elsif ($totsize-$$exon{size}<$size5) {#end of 5'UTR & start of CDS
            my $exonlen = $totsize-$size5;#need to split exon
            if ($$p{RM_UTRS}) {
                (my $utr, $exon) = split_segment($exon, $exonlen, $strand);
                $$utr{type} = 'utr5';
                push(@utrs, $utr);
                push(@exons, $exon);
            } else {
                $$exon{type} = 'utr5-'.$$exon{type};
                push(@exons, $exon);
            }
        } elsif ($totsize>$cds_end) {#in 3'UTR
            if ($totsize-$$exon{size}<$cds_end) {#end of CDS
                if ($$p{RM_UTRS}) {
                    my $utrlen = $totsize-$cds_end;
                    ($exon, my $utr) = split_segment($exon, $utrlen, $strand);
                    $$utr{type} = 'utr3';
                    push(@exons, $exon);
                    push(@utrs, $utr);
                } else {
                    $$exon{type} .= '-utr3';
                    push(@exons, $exon);
                }
            } else {
                $$exon{type} = 'utr3';
                $$p{RM_UTRS} ? push(@utrs, $exon) : push(@exons, $exon);
            }
        } else {#in CDS
            push(@exons, $exon);
        }
    }
    $$info{numexons} = scalar(@exons);
    $$info{exons} = \@exons;
    $$info{utrs} = \@utrs;
    my $i=0;
    foreach my $utr (@utrs) {
        if (($strand && $$utr{start}>$exons[0]{start}) ||
            (!$strand && $$utr{start}<$exons[0]{start})) {
            for(;$i<@exons; $i++) {
                my $exon = $exons[$i];
                warn "\t  Exon  \t$$exon{start}-$$exon{end}\n" if $$p{DBUG};
            }
        }
        warn "\t  UTR   \t$$utr{start}-$$utr{end}\n" if $$p{DBUG};
    }
}

sub get_snps {
    my $geneENS = shift;
    my $registry = shift;
    my $data = shift;
    my $snphash = shift;
    my $p = shift;

    warn "  Getting vf adaptor\n";
    my $vf_adaptor = $registry->get_adaptor('human', 'variation',
                                        'variationfeature');
    my $s = $$data{start}-$$p{FLANK};
    my $e = $$data{end}+$$p{FLANK};
    my $slice = $geneENS->slice->sub_Slice($s, $e);

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
    warn "  Total of  ".scalar(keys %snps)." snps\n";
    \%snps;
}

# retrieve snps in region of this gene
sub create_sequences {
    my $gene = shift;
    my $slice = shift;
    my $data = shift;
    my $p = shift;
    my $flank = $$p{FLANK};
    my $limit = 60;

    unless ($slice) {
        warn "COULD not create sequence for $gene!\n";
        return;
    }
    my $rsfh = *STDERR;
    if ($$p{RSSUMMARY}) { # write rs snps to file
        $rsfh = get_filehandle($gene.'.snps.txt', $p);
        print $rsfh join("\t", @{ $$p{RSSUMMARY} });
        print $rsfh " (MAF>=$$p{MAF})" if $$p{MAF};
        print $rsfh "\n";
    }
    my @seq;
    my @nseq;
    my @transcripts;
    my $intron_s = 0;
    my $snps = $$data{snps} || {};
    # If gene on reverse strand, we must reverse complement it
    my $chrom = $$data{chrom};
    my $geneinfo = "$$data{genename} $$data{geneid}";
    my $strand = $$data{plusstrand} ? 1 : -1;
    my @exons = $$data{plusstrand} ? 
        @{ $$data{exons} } : reverse @{ $$data{exons} };
    my $numcoding = grep($$_{type} =~ /protein_coding/, @exons);
    my $exonnum = $$data{plusstrand} ? 1 : $numcoding;
    my $intronnum = $exonnum;
    my ($start, $end); # start and end coordinates of seq relative to Ensembl
    my @exon_locs = ();
    for(my $i=0; $i<@exons; $i++) {
        my $ex = $exons[$i];
        if ($$p{DBUG}) {
            warn "create_sequences:exon $i\n";
            warn hash_debug($ex)."\n";
        }
        if ($i==0) {
            my ($s, $e) = ($$ex{start}-$flank, $$ex{start}-1);
            my $label = $$data{plusstrand} ? 'flank_left' : 'flank_right';
            my $seqcomment = "chr$chrom:$s:$e:$strand\t$geneinfo";
            my $seqhead = ">$gene:$label\t$seqcomment\n";
            $start = $s;
            printf STDERR "  Adding $label:\t%d-%d (%d bp)\n", 
                $s, $e, $flank if $$p{DBUG};
            my $seq = lc $slice->subseq($s, $e);
            my $err = check_seq($seq, $s, $e, $gene, $label);
            push(@seq, $seqhead.$seq.$err);
            my $nseq = n_snps($snps, $seq, $s, $e, $rsfh, $p); 
            push(@nseq, $seqhead.$nseq.$err);
        }
        if (($intron_s>0) && ($$ex{start}-$intron_s>0)) {
            my $size = $$ex{start}-$intron_s;
            my ($s, $e) = ($intron_s, $$ex{start}-1);
            my $label = "intron_$intronnum";
            my $seqcomment = "chr$chrom:$s:$e:$strand\t$geneinfo";
            my $seqhead = ">$gene:$label\t$seqcomment\n";
            printf STDERR "  Adding intron:\t%d-%d (%d bp)\n",
                $s, $e, $size if $$p{DBUG};
            if ($s>=$e) {
                warn "exstart $$ex{start}, intronstart $intron_s\n";
            }
            my $seq = lc $slice->subseq($s, $e);
            my $err = check_seq($seq, $s, $e, $gene, $intronnum);
            push(@seq, $seqhead.$seq.$err);
            my $nseq = n_snps($snps, $seq, $s, $e, $rsfh, $p);
            push(@nseq, $seqhead.$nseq.$err);
            $intronnum = $$data{plusstrand} ? $intronnum+1 : $intronnum-1;
        }
        my $type = $$ex{type};
        printf STDERR "  Adding exon%02d:\t%d-%d (%d bp)\n", $exonnum,
            $$ex{start}, $$ex{end}, $$ex{size} if $$p{DBUG};
        my $label = "$$ex{name}_${type}";
        if ($type =~ /protein_coding/) {
            $label .= "_$exonnum";
            $exonnum = $$data{plusstrand} ? $exonnum+1 : $exonnum-1;
        }
        my $seqcomment = "chr$chrom:$$ex{start}:$$ex{end}:$strand\t$geneinfo";
        my $seqname = $$ex{nickname} ? $$ex{nickname} : "$gene:$label";
        my $seqhead = ">$seqname\t$seqcomment\n";
        my $seq = uc $slice->subseq($$ex{start}, $$ex{end});
        my $err = check_seq($seq, $$ex{start}, $$ex{end}, $gene, $label);
        push(@seq, $seqhead.$seq.$err);
        push(@transcripts, $seqhead.$seq.$err);
        my $nseq = n_snps($snps, $seq, $$ex{start}, $$ex{end}, $rsfh, $p);
        push(@nseq, $seqhead.$nseq.$err);
        $intron_s = $$ex{end}+1;
    }
    $end = $intron_s+$flank-1;
    my $label = $$data{plusstrand} ? 'flank_right' : 'flank_left';
    my $seqcomment = "chr$chrom:$intron_s:$end:$strand\t$geneinfo";
    my $seqhead = ">$gene:$label\t$seqcomment\n";
    printf STDERR "  Adding $label:\t%d-%d (%d bp)\n", $intron_s, 
        $end, $flank if $$p{DBUG};
    my $seq = lc $slice->subseq($intron_s, $end);
    my $nseq = n_snps($snps, $seq, $intron_s, $end, $rsfh, $p);
    $rsfh->close();
    push(@seq, $seqhead.$seq);
    push(@nseq, $seqhead.$nseq);
    if (!$$data{plusstrand}) { # need to reverse complement
        @seq = map { rev_complement($_) } reverse @seq;
        @nseq = map { rev_complement($_) } reverse @nseq;
        @transcripts = map { rev_complement($_) } reverse @transcripts;
    }
    if ($$p{MAKE_SEQS}) {
        my @createlist = (["$gene.gene.fa", \@seq], ["$gene.exons.fa", \@transcripts]);
#        if ($$p{SNPS}) { push(@createlist, ["$gene.nseq.fa", \@nseq]) }
        foreach my $type (@createlist) {
            my ($fafile, $seqdat) = @$type;
            my $ofh = get_filehandle($fafile, $p);
            print $ofh join("\n", @$seqdat)."\n";
            $ofh->close;
        }
    }
}



# Make sure sequence is of expected length
# In past, sometimes queries to ensembl have silently failed
# creating incomplete sequences
sub check_seq {
    my $seq = shift;
    my $s = shift;
    my $e = shift;
    my $gene = shift;
    my $type = shift;

    my $len = $e-$s+1;
    my $seqlen = length($seq);
    my $err = '';
    if ($len != $seqlen) {
        $err = "\nError printing seq for $gene $type ($s, $e)\n";
        $err .= "Length of seq is $seqlen bp but expect $len bp\n";
        warn "$err";
    }
    ($err);
}

sub n_snps {
    my $snps = shift;
    my $seq = shift;
    my $chr_s = shift;
    my $chr_e = shift;
    my $rsfh = shift;
    my $p = shift;

    my @snppos = grep($_>=$chr_s-1 && $_<$chr_e, keys %$snps);
    my $n = $seq =~ /[ACTG]/ ? 'N' : 'n';
    foreach my $chrpos (sort {$a<=>$b} @snppos) {
        my $seqpos = $chrpos - $chr_s + 1;
#        my $base = substr($seq, $seqpos, 70);
        if ($$p{RSSUMMARY}) {
            my @row = map { defined($$snps{$chrpos}{$_}) ?
                $$snps{$chrpos}{$_} : '' } @{ $$p{RSSUMMARY} };
            print $rsfh join("\t", @row)."\n";
        }
        substr($seq, $seqpos, 1) = $n;
    }
    ($seq);
}

# Split a segment into two segments where the second segment
# has the given length.
# Use for removing UTRs.
sub split_segment {
    my $segment = shift;
    my $len = shift;
    my $strand = shift;

    my $segBlen = $strand ? $len : $$segment{size}-$len;
    my %segA = %$segment;
    my %segB = %$segment;
    $segA{end} = $segA{end}-$segBlen;
    $segB{start} = $segB{end}-$segBlen+1;
    $segA{size} = $segA{size}-$segBlen;
    $segB{size} = $segBlen;
    $strand ? (\%segA, \%segB) : (\%segB, \%segA);
}


sub feature2string { # for debugging
    my $feature = shift;

    my $stable_id  = $feature->stable_id();
    my $name       = $feature->display_xref->display_id() || '';
    my $desc       = $feature->description() || '';
    my $seq_region = $feature->slice->seq_region_name();
    my $start      = $feature->start();
    my $end        = $feature->end();
    my $strand     = $feature->strand();
    my $size = $end - $start + 1;
    if (length($desc)>60) { $desc = substr($desc, 0, 60) . '...' }
#    my $str = hash_debug($feature);
#    return $str;
    return ("ID $stable_id;\tName: $name\n".
        sprintf "    Desc: $desc\n    Seq_region: %s:%d-%d\tStrand(%+d) %d bp",
        $seq_region, $start, $end, $strand, $size);
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


#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $vcf_a_file;
my $vcf_b_file;
my $output_f1_vcf_file;
my $help;

GetOptions(
    "vcf_a=s" => \$vcf_a_file,
    "vcf_b=s" => \$vcf_b_file,
    "out=s"   => \$output_f1_vcf_file,
    "help|h"  => \$help,
) or die usage();

usage() if $help || !$vcf_a_file || !$vcf_b_file || !$output_f1_vcf_file;

sub usage {
    print STDERR <<EOF;
Usage: $0 --vcf_a <species_A.vcf> --vcf_b <species_B.vcf> --out <f1_output.vcf>

Description:
  Creates an F1 generation VCF from two parent VCF files.
  Assumes 5 individuals in species A VCF and 5 individuals in species B VCF.
  F1_1 = A1 x B1, F1_2 = A2 x B2, ..., F1_5 = A5 x B5.
  The output F1 VCF will only contain the GT field in FORMAT.

Required arguments:
  --vcf_a  FILE   VCF file for species A (5 individuals)
  --vcf_b  FILE   VCF file for species B (5 individuals)
  --out    FILE   Output VCF file for the 5 F1 individuals

Optional arguments:
  --help, -h      Show this help message and exit.
EOF
    exit 1;
}

my %variants_data; # Stores GT data and variant info
my @A_sample_names;
my @B_sample_names;
my %header_lines; # To collect ## lines

print STDERR "Reading Species A VCF: $vcf_a_file\n";
open(my $VCF_A, "<", $vcf_a_file) or die "Cannot open $vcf_a_file: $!";
while (my $line = <$VCF_A>) {
    chomp $line;
    if ($line =~ /^##/) {
        $header_lines{$line} = 1 unless $line =~ /^##FORMAT=<ID=GT/; # Collect header, exclude GT format for now
        next;
    }
    if ($line =~ /^#CHROM/) {
        my @fields = split /\t/, $line;
        @A_sample_names = @fields[9 .. $#fields];
        if (scalar(@A_sample_names) != 5) {
            die "Error: Species A VCF ($vcf_a_file) does not contain exactly 5 samples in the header. Found " . scalar(@A_sample_names) . ".\n";
        }
        next;
    }

    my @fields = split /\t/, $line;
    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info_str, $format_str) = @fields[0..8];
    my @genotype_fields_A = @fields[9 .. $#fields];

    if (scalar(@genotype_fields_A) != 5 && scalar(@A_sample_names) == 5) {
         print STDERR "Warning: Data line in $vcf_a_file at $chrom:$pos has " . scalar(@genotype_fields_A) . " genotype columns, expected 5. Skipping line.\n";
         next;
    } elsif (scalar(@genotype_fields_A) != scalar(@A_sample_names)) {
         print STDERR "Warning: Data line in $vcf_a_file at $chrom:$pos has inconsistent number of genotype columns. Skipping line.\n";
         next;
    }


    my $variant_key = "$chrom\t$pos\t$ref\t$alt";
    $variants_data{$variant_key}{INFO_FIELDS} = {
        CHROM  => $chrom, POS => $pos, ID => $id, REF => $ref, ALT => $alt,
        QUAL   => $qual, FILTER => $filter, INFO => $info_str
    };

    my @format_keys = split /:/, $format_str;
    my $gt_idx_A = -1;
    for my $i (0 .. $#format_keys) {
        if ($format_keys[$i] eq 'GT') {
            $gt_idx_A = $i;
            last;
        }
    }
    die "Error: GT field not found in FORMAT string ($format_str) in $vcf_a_file at $chrom:$pos\n" if $gt_idx_A == -1;

    for my $i (0 .. $#A_sample_names) {
        my $gt_full_A = $genotype_fields_A[$i];
        my $gt_val_A = (split /:/, $gt_full_A)[$gt_idx_A] // './.'; # Default to missing if GT not fully specified
        $variants_data{$variant_key}{A_PARENTS}[$i] = $gt_val_A;
    }
}
close $VCF_A;

print STDERR "Reading Species B VCF: $vcf_b_file\n";
open(my $VCF_B, "<", $vcf_b_file) or die "Cannot open $vcf_b_file: $!";
while (my $line = <$VCF_B>) {
    chomp $line;
    if ($line =~ /^##/) {
        $header_lines{$line} = 1 unless $line =~ /^##FORMAT=<ID=GT/;
        next;
    }
    if ($line =~ /^#CHROM/) {
        my @fields = split /\t/, $line;
        @B_sample_names = @fields[9 .. $#fields];
        if (scalar(@B_sample_names) != 5) {
            die "Error: Species B VCF ($vcf_b_file) does not contain exactly 5 samples in the header. Found " . scalar(@B_sample_names) . ".\n";
        }
        next;
    }

    my @fields = split /\t/, $line;
    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info_str, $format_str) = @fields[0..8];
    my @genotype_fields_B = @fields[9 .. $#fields];

    if (scalar(@genotype_fields_B) != 5 && scalar(@B_sample_names) == 5) {
         print STDERR "Warning: Data line in $vcf_b_file at $chrom:$pos has " . scalar(@genotype_fields_B) . " genotype columns, expected 5. Skipping line.\n";
         next;
    } elsif (scalar(@genotype_fields_B) != scalar(@B_sample_names)) {
         print STDERR "Warning: Data line in $vcf_b_file at $chrom:$pos has inconsistent number of genotype columns. Skipping line.\n";
         next;
    }

    my $variant_key = "$chrom\t$pos\t$ref\t$alt";
    unless (exists $variants_data{$variant_key}{INFO_FIELDS}) {
        $variants_data{$variant_key}{INFO_FIELDS} = {
            CHROM  => $chrom, POS => $pos, ID => $id, REF => $ref, ALT => $alt,
            QUAL   => $qual, FILTER => $filter, INFO => $info_str
        };
    }

    my @format_keys = split /:/, $format_str;
    my $gt_idx_B = -1;
    for my $i (0 .. $#format_keys) {
        if ($format_keys[$i] eq 'GT') {
            $gt_idx_B = $i;
            last;
        }
    }
    die "Error: GT field not found in FORMAT string ($format_str) in $vcf_b_file at $chrom:$pos\n" if $gt_idx_B == -1;
    
    for my $i (0 .. $#B_sample_names) {
        my $gt_full_B = $genotype_fields_B[$i];
        my $gt_val_B = (split /:/, $gt_full_B)[$gt_idx_B] // './.';
        $variants_data{$variant_key}{B_PARENTS}[$i] = $gt_val_B;
    }
}
close $VCF_B;

print STDERR "Generating F1 VCF: $output_f1_vcf_file\n";
open(my $OUT_F1, ">", $output_f1_vcf_file) or die "Cannot open $output_f1_vcf_file for writing: $!";

print $OUT_F1 "##fileformat=VCFv4.2\n";
# }
print $OUT_F1 "##source=$0\n"; # Indicate script source
print $OUT_F1 "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";

my @f1_sample_names = map { "F1_Individual_$_" } (1..5);
print $OUT_F1 join("\t", "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", @f1_sample_names) . "\n";

my @sorted_variant_keys = sort {
    my ($chrom_a, $pos_a) = ($a =~ /^([^\t]+)\t(\d+)/);
    my ($chrom_b, $pos_b) = ($b =~ /^([^\t]+)\t(\d+)/);
    return ($chrom_a cmp $chrom_b) || ($pos_a <=> $pos_b);
} keys %variants_data;

for my $key (@sorted_variant_keys) {
    my $variant_site_data = $variants_data{$key};
    my $info_fields = $variant_site_data->{INFO_FIELDS};

    my @f1_gts_for_this_line;

    for my $i (0 .. 4) { # For F1_1 to F1_5 (corresponds to A_parent[i] x B_parent[i])
        my $gt_a_parent = ($variant_site_data->{A_PARENTS} && defined $variant_site_data->{A_PARENTS}[$i])
                        ? $variant_site_data->{A_PARENTS}[$i] : './.';
        my $gt_b_parent = ($variant_site_data->{B_PARENTS} && defined $variant_site_data->{B_PARENTS}[$i])
                        ? $variant_site_data->{B_PARENTS}[$i] : './.';
        
        my $allele_a = (split /[\/\|]/, $gt_a_parent)[0];
        my $allele_b = (split /[\/\|]/, $gt_b_parent)[0];
        
        $allele_a = '.' if !defined $allele_a || $allele_a eq '';
        $allele_b = '.' if !defined $allele_b || $allele_b eq '';

        my $f1_gt;
        if ($allele_a eq '.' || $allele_b eq '.') {
            $f1_gt = './.';
        } else {
            my @parent_alleles_for_f1;
            if ($allele_a =~ /^\d+$/ && $allele_b =~ /^\d+$/) {
                @parent_alleles_for_f1 = sort {$a <=> $b} ($allele_a, $allele_b);
            } else { 
                print STDERR "Warning: Non-numeric alleles ('$allele_a', '$allele_b') found at $info_fields->{CHROM}:$info_fields->{POS}. Using string sort.\n"
                    unless ($allele_a =~ /^\d+$/ && $allele_b =~ /^\d+$/); # Warn only if truly not digits
                @parent_alleles_for_f1 = sort {$a cmp $b} ($allele_a, $allele_b);
            }
            $f1_gt = join '/', @parent_alleles_for_f1;
        }
        push @f1_gts_for_this_line, $f1_gt;
    }

    print $OUT_F1 join("\t",
        $info_fields->{CHROM}, $info_fields->{POS}, $info_fields->{ID},
        $info_fields->{REF}, $info_fields->{ALT}, $info_fields->{QUAL},
        $info_fields->{FILTER}, $info_fields->{INFO}, "GT",
        @f1_gts_for_this_line
    ) . "\n";
}

close $OUT_F1;
print STDERR "F1 VCF generation complete: $output_f1_vcf_file\n";

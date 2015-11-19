#!/bin/env perl

use FindBin;               
use lib "$FindBin::Bin/..";

use strict;

use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . "/../GFFLib";
use lib $FindBin::Bin . "/../Orthologia";


use GD::SVG;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Getopt::Std;
use GFFLib::GFFFile;
use Orthologia::Ort;

my $DEFAULT_WIDTH  = 1000;
my $DEFAULT_HEIGHT = 600;
my $DEFAULT_COLOR  = "white";
my $MARGIN         = 20;

my $CHROM_HEIGHT = 20;
my $CHROM_SPACER = 200;

# Position of gene caption relative to top left corner of 
# rectangle representing the gene
my $TEXT_H_POS_RELATIVE_GENE = 0;
my $TEXT_V_POS_RELATIVE_GENE = -15;


# Space between scaffold/chrom in terms of base pairs
# L. loa
my $CHROM_WIDTH_SPACER_IN_BP = 48000;

# A. darlingi
# my $CHROM_WIDTH_SPACER_IN_BP = 12000;

my $debug_projections = 1;
my $debug = 0;

STDOUT->autoflush(1);

my @pallete = (
	{ red => 200, green => 80,  blue => 80 },
	{ red => 80,  green => 200, blue => 80 },
	{ red => 80,  green => 80,  blue => 200 },
	{ red => 200, green => 200, blue => 80 },
	{ red => 80,  green => 200, blue => 200 },
	{ red => 200, green => 80,  blue => 200 },
);



# A. darlingi
#my $gene_color_pallete_index = 1;

# L. loa
my $gene_color_pallete_index = 2;

my $usage =
"$0 <orts file> <list chrom/scaffolds> <list gff files> <list fasta files> <svg out>\n";

die $usage if ( scalar(@ARGV) != 5 );

my $orts_file  = $ARGV[0];
my $list_chrom = $ARGV[1];
my $list_gff   = $ARGV[2];
my $list_fasta = $ARGV[3];
my $outputFile = $ARGV[4];

my $width = $DEFAULT_WIDTH;

my $overallMaxLength = 0;

my @fasta_seq;
my @fasta_seqnames;



######################################
# Reads Chrom/Scaffolds list

#FORMAT
# <org name on cluster file> <scaffold/chrom>(start:end)orientation ...

# "(start:end)orientation" are optional. No need them if the whole chrom/scaffold will be rendered or if it will be rendered
# in its current orientation
# A start or end equals to -1 indicates to the program to render the chrom/scaffold untils its beginning (start=-1)
# or untils its end (end=-1)  

# Loa_loa_V2 7000000145608817(1:200)- 7000000145609071(-1:400) 7000000145608793
# C_elegans_WS224 7000000183869869-

# * A dash (minus) suffix indicates sequences that should be render as reverse complement

# Open list of chrom
my @chroms;
my @species;

my $cont_species = 0;

# Hash used to get rid second copy of the same contig
# Only the first instance of a contig will be rendered
my @chrom_registered;

print STDERR "Reading list of chrom/scaffolds...\n";
open LIST_CHROM, "$list_chrom" or die "Unable to open $list_chrom\n";
while (<LIST_CHROM>) {
	my $line = $_;
	chomp($line);

	my @temp_chroms = split /\s/, $line;

	# Add species name to @species vector
	push( @species, shift @temp_chroms );

	foreach my $curr_chrom (@temp_chroms) {

		my $orientation = "+";

		if ( $curr_chrom =~ /-$/ ) {
			$orientation = "-";
			$curr_chrom =~ s/-$//;
		}
		elsif ( $curr_chrom =~ /\+$/ ) {
			$orientation = "+";
			$curr_chrom =~ s/\+$//;
		}

		my $start;
		my $end;

		if ( ( $start, $end ) = ( $curr_chrom =~ /\((\d+):(\d+)\)/ ) ) {
			$curr_chrom =~ s/\([\W\w]+//;
		}else{
			$start = -1;
			$end   = -1;
		}
		
		my $start_str = $start;
		my $end_str = $end;
		
		$start_str = 'undefined' if $start_str == -1;
		$end_str = 'undefined' if $end_str == -1; 

			print STDERR
"\tChrom.: $curr_chrom  Orient.: $orientation  start: $start_str  end: $end_str\n";

		#getc();

		# Add chrom/contig to $chroms[$cont_species] vector only if
		# this is the first instance of this chrom/contig
		push(
			@{ $chroms[$cont_species] },
			{
				chrom       => $curr_chrom,
				orientation => $orientation,
				start       => $start,
				end         => $end
			}
		) if not defined $chrom_registered[$cont_species]{$curr_chrom};

		# Register the chrom by placing its index on $chrom in the hash below
		my $index = scalar(@{$chroms[$cont_species]}) - 1;
		$chrom_registered[$cont_species]{$curr_chrom} = $index;
		print "Registering chrom \'$curr_chrom\' of species $cont_species as index: $index\n" if $debug;
		getc() if $debug;
	}

	$cont_species++;
}
close(LIST_CHROM);

#map {print $_} @{$chroms[0]};
#getc();
#foreach my $curr_seq_name ( @{$chroms[ 0 ] } ){
#	print $curr_seq_name . " ";
#}
#getc();

######################################
# Reads GFFs

$cont_species = 0;
my @genesHash;
open LIST_GFF, "$list_gff" or die "Unable to open $list_gff\n";
while (<LIST_GFF>) {
	my $gff_file = $_;
	chomp($gff_file);

	print STDERR "Reading GFF from species $cont_species ...\n";
	my $gffFile = GFFFile::new($gff_file);
	$gffFile->read();
	$genesHash[$cont_species] = $gffFile->get_genes_hash();

	$cont_species++;
}
close(LIST_GFF);

######################################
# Reads FASTAs

$cont_species = 0;
my @fasta_seq;
my @fasta_seq_names;
my $overallMaxLength = 0;
my @start_pixel;
my @sum_chrom_length_arr;
my $pixel_bp_ratio;

open LIST_FASTA, "$list_fasta" or die "Unable to open $list_fasta\n";
while (<LIST_FASTA>) {
	my $fasta_file = $_;
	chomp($fasta_file);

	my $sum_chrom_length = 0;

	print STDERR "Reading FASTA from species $cont_species ...\n";
	my $inSeqIO = Bio::SeqIO->new( -file => $fasta_file, '-format' => 'Fasta' );
	while ( my $inSeq = $inSeqIO->next_seq() ) {

		my $inSeqID = $inSeq->id();

		# Next if this is not a chosen chrom
		my $found = 0;
		map { $found = 1 if $inSeqID eq $_->{chrom} }
		  @{ $chroms[$cont_species] };
		next if $found == 0;
		

		my $curr_chrom = $inSeqID;
		
		print STDERR "\tParsing chrom/scaffold $curr_chrom\n";

		# Mark that this scaffold was found
		$chroms[$cont_species]
			  [ $chrom_registered[$cont_species]{$curr_chrom} ]->{found} =  1;

		if ( !defined( $inSeq->seq() ) || length( $inSeq->seq() ) == 0 ) {
			print STDERR $inSeq->id()
			  . " have an empty string as sequence. Sequence discarded\n";
		}
		else {
			my $seq_length = $inSeq->length();
			my $seq_name   = $inSeq->id();

			my $start =
			  $chroms[$cont_species]
			  [ $chrom_registered[$cont_species]{$curr_chrom} ]->{start};
			my $end =
			  $chroms[$cont_species]
			  [ $chrom_registered[$cont_species]{$curr_chrom} ]->{end};

			# If START was not defined in the text file
			if ( $start == -1 ) {
				$start = 1;
				$chroms[$cont_species]
				  [ $chrom_registered[$cont_species]{$curr_chrom} ]->{start} =
				  1;
			}

			# If END was not defined in the text file
			if ( $end == -1 ) {
				$end = $seq_length;
				$chroms[$cont_species]
				  [ $chrom_registered[$cont_species]{$curr_chrom} ]->{end} =
				  $seq_length;
			}
			my $drawing_length = abs( $end - $start ) + 1;

			push( @{ $fasta_seq_names[$cont_species] }, $seq_name );
			$fasta_seq[$cont_species]{$seq_name}{len} = $seq_length;

			$fasta_seq[$cont_species]{$seq_name}{start} = $start;
			$fasta_seq[$cont_species]{$seq_name}{end}   = $end;

			$fasta_seq[$cont_species]{$seq_name}{drawing_length} =
			  $drawing_length;
			  
			print "Registering \'$seq_name\' start: $start end: $end\n" if $debug;
			getc() if $debug;

			$sum_chrom_length += $drawing_length + $CHROM_WIDTH_SPACER_IN_BP;
		}
	}
	
	# Checking if all scaffolds/chroms were found
	my $error_message = '';
	foreach my $checking_chrom ( keys %{$chrom_registered[$cont_species]} ){ 
		$error_message .= "Error: Unable to find chrom/scaffold $checking_chrom on file $fasta_file\n"
		  if ( not defined( $chroms[$cont_species]
				  [ $chrom_registered[$cont_species]{$checking_chrom} ]->{found} ) );
	}	
	die $error_message if ( $error_message ne '' ); 

	$sum_chrom_length -= $CHROM_WIDTH_SPACER_IN_BP;
	$sum_chrom_length_arr[ $cont_species ] = $sum_chrom_length;

	die "Error: No chrom/contigs selected for species $cont_species!!!\n"
	  if $sum_chrom_length == 0;


	$overallMaxLength = $sum_chrom_length
	  if ( $overallMaxLength < $sum_chrom_length );

	$inSeqIO->close();
	$cont_species++;
}
close(LIST_FASTA);

$pixel_bp_ratio = ( $DEFAULT_WIDTH - 2 * $MARGIN ) / $overallMaxLength;

for ( my $cont_species = 0 ; $cont_species < scalar(@chroms) ; $cont_species++ ){
	$start_pixel[ $cont_species ] = int( ($overallMaxLength - $sum_chrom_length_arr[ $cont_species ]) / 2 * $pixel_bp_ratio );
} 

#################
# Reading Orts
my $orts = Ort::new( $orts_file, "gene", "RBH" );
$orts->read();

# create a new image
my $im = new GD::SVG::Image( $DEFAULT_WIDTH, $DEFAULT_HEIGHT );

# allocate some colors
my $white = $im->colorAllocate( 255, 255, 255 );
my $gray  = $im->colorAllocate( 200, 200, 200 );
my $black = $im->colorAllocate( 0,   0,   0 );
my $red   = $im->colorAllocate( 255, 0,   0 );
my $blue  = $im->colorAllocate( 0,   0,   255 );

# make the background transparent and interlaced
#$im->transparent($white);
#$im->interlaced('true');

my $chrom_cont = 0;
my %chrom_height_offset;
my %chrom_color;

# Ref seq
print STDERR "Rendering chroms ...\n";
my @species_height_offset;

for ( my $cont_species = 0 ; $cont_species < scalar(@chroms) ; $cont_species++ )
{
	my $width_offset = $MARGIN + $start_pixel[ $cont_species ];

	foreach my $curr_chrom ( @{ $chroms[$cont_species] } ) {
		my $curr_seq_name = $curr_chrom->{chrom};
		print STDERR
		  "\tRendering species $cont_species seq. $curr_seq_name ...\n";
		my $curr_seq_length =
		  $fasta_seq[$cont_species]{$curr_seq_name}{drawing_length};

		my $height_offset =
		  $MARGIN + ( $cont_species * ( $CHROM_SPACER + $CHROM_HEIGHT ) );
		$species_height_offset[$cont_species] = $height_offset;

		$fasta_seq[$cont_species]{$curr_seq_name}{x_offset} = $width_offset;
		$fasta_seq[$cont_species]{$curr_seq_name}{y_offset} = $height_offset;

		$im->rectangle(
			$width_offset,
			$height_offset,
			int( $curr_seq_length * $pixel_bp_ratio ) +
			  $width_offset,
			$height_offset + $CHROM_HEIGHT,
			$black
		);

		$width_offset +=
		  int( ( $curr_seq_length + $CHROM_WIDTH_SPACER_IN_BP ) *
			  $pixel_bp_ratio );

	}
}

######################
# Rendering genes
print STDERR "Rendering genes...\n";
my %gene_figure;
my @genes_per_genome;

my $genes_color = $im->colorAllocateAlpha(
	$pallete[$gene_color_pallete_index]{red},
	$pallete[$gene_color_pallete_index]{green},
	$pallete[$gene_color_pallete_index]{blue}, 0
);

for ( my $cont_species = 0 ; $cont_species < scalar(@chroms) ; $cont_species++ )
{
	my $gffGenesRef = $genesHash[$cont_species];

	foreach my $currGene ( keys %{$gffGenesRef} ) {

		my $contig_chrom = $gffGenesRef->{$currGene}->get_chrom();
		my $id           = $gffGenesRef->{$currGene}->get_id();
		my $alias        = $gffGenesRef->{$currGene}->get_attribute("alias");
		my $name         = $gffGenesRef->{$currGene}->get_name();
		my $gene_start   = $gffGenesRef->{$currGene}->get_start();
		my $gene_end     = $gffGenesRef->{$currGene}->get_end();
		my $strand       = $gffGenesRef->{$currGene}->get_strand();

		my $found = 0;
		my $orientation;
		my $chrom_start;
		my $chrom_end;

		map {
			if ( $contig_chrom eq $_->{chrom} )
			{
				$found = 1;

				$orientation = $_->{orientation};
				$chrom_start = $_->{start};
				$chrom_end   = $_->{end};

			}
		} @{ $chroms[$cont_species] };

		# Next if this is not a chosen chrom
		next if $found == 0;

# Invert orientation of chromosome if that's was indicated in the chrom list input file
		my $start = $gene_start;
		my $end   = $gene_end;

		# Next if gene not part of the region being depicted
		next if $chrom_start > $start;
		next if $chrom_end < $end;

		my $curr_chrom_start = $fasta_seq[$cont_species]{$contig_chrom}{start};

		my $curr_chrom_end = $fasta_seq[$cont_species]{$contig_chrom}{end};

		print STDERR
"Species $cont_species $contig_chrom $orientation chrom.start: $curr_chrom_start chrom.end: $curr_chrom_end  gene start: $start gene end: $end\n" if $debug;
		getc() if $debug;


		if ( $orientation eq "+" ) {
			$start = $start - $curr_chrom_start + 1;
			$end   = $end - $curr_chrom_start + 1;

		}
		elsif ( $orientation eq "-" ) {
			my $temp = $start;
			$start = $curr_chrom_end - $end + 1;
			$end   = $curr_chrom_end - $temp + 1;
		}

		print STDERR
"Species $cont_species $contig_chrom $orientation chrom.start: $curr_chrom_start chrom.end: $curr_chrom_end  gene start: $start gene end: $end\n" if $debug;
		getc() if $debug;

		my $width_offset  = $fasta_seq[$cont_species]{$contig_chrom}{x_offset};
		my $height_offset = $fasta_seq[$cont_species]{$contig_chrom}{y_offset};

		my @rectangle;
		if (   ( $strand eq "+" && $orientation eq "+" )
			|| ( $strand eq "-" && $orientation eq "-" ) )
		{
			@rectangle = (
				$width_offset + int( $start * $pixel_bp_ratio ),
				$height_offset + 1,
				$width_offset + int( $end * $pixel_bp_ratio ),
				$height_offset + int( $CHROM_HEIGHT / 2 ) - 1
			);

		}
		elsif (( $strand eq "-" && $orientation eq "+" )
			|| ( $strand eq "+" && $orientation eq "-" ) )
		{
			@rectangle = (
				$width_offset + int( $start * $pixel_bp_ratio ),
				$height_offset + int( $CHROM_HEIGHT / 2 ) + 1,
				$width_offset + int( $end * $pixel_bp_ratio ),
				$height_offset + $CHROM_HEIGHT - 1
			);
		}
		else {
			die "unrecognized strand $strand\n";
		}

		@{ $gene_figure{$id}{rectangle} } = (
			$rectangle[0], $height_offset, $rectangle[2],
			$height_offset + $CHROM_HEIGHT
		);
		$gene_figure{$id}{species} = $cont_species;
		push( @{ $genes_per_genome[$cont_species] }, $id );

		my $poly = new GD::SVG::Polygon;
		
		
		
		$poly->addPt( $rectangle[0], $rectangle[1] );
		$poly->addPt( $rectangle[2], $rectangle[1] );

		$poly->addPt( $rectangle[2], $rectangle[3] );
		$poly->addPt( $rectangle[0], $rectangle[3] );

		print STDERR
"Species $cont_species: adding gene $id Name: $name Chrom/Scaffold: $contig_chrom Start: $gene_start End: $gene_end\n";
		print STDERR
"Coords. LEFT: $rectangle[0]  RIGHT: $rectangle[2]  TOP: $rectangle[1] BOTTOM: $rectangle[3] \n\n" if $debug;

		$im->filledPolygon( $poly, $genes_color );
		
		#if( $name =~ "KPC" ){
			$im->stringUp(gdSmallFont, $rectangle[0] + $TEXT_H_POS_RELATIVE_GENE,
						$height_offset + $TEXT_V_POS_RELATIVE_GENE , $id ,$black);
		#}
		

#		if( $name =~ "KPC" ){
#			$im->string(gdSmallFont, $rectangle[0] + $TEXT_H_POS_RELATIVE_GENE,
#						$height_offset + $TEXT_V_POS_RELATIVE_GENE , "KPC" ,$black);
#			$im->filledPolygon( $poly, $red );
#		}else{
#			$im->filledPolygon( $poly, $black );
#		}
			

	}
}

###########################
# Drawing ortholog line/projections

# Projection or line color
my $projection_fill_color;
my $projection_outline;

my %homologous_block_color;
my $genome_header_height = 15;

print STDERR "\n\nDrawing orhtology projections ...\n";

# Looping through genome numbers, user input. The index indicates the source of projections
for (
	my $curr_index = 0 ;
	$curr_index < ( scalar(@chroms) - 1 ) ;
	$curr_index++
  )
{

	my $curr_genome = $curr_index;

	# Loop through genes of each source genome
	foreach my $curr_gene ( @{ $genes_per_genome[$curr_genome] } ) {

		my $gene_name = $curr_gene;

		print STDERR "Gene name: $gene_name\n" if $debug_projections == 1;

		my @src_rectangle = @{ $gene_figure{$gene_name}{rectangle} };

		print STDERR "Src rectangle: $src_rectangle[ 0 ], $src_rectangle[ 1 ]\n"
		  if $debug_projections == 1;

# Draw projection to ortholog in genome just below, if there are no orhtologs in that genome
# try the next genome until it finds at least one ortholog

		my $dest_genome = $curr_index + 1;

		print STDERR "Dst genome: $dest_genome\n" if $debug_projections == 1;
		getc() if $debug_projections == 1;

		my $org_name = $species[$dest_genome];
		my @ort_genes = $orts->get_orts( $gene_name, $org_name );

		if( $debug_projections == 1 ){
			print STDERR "Number of member from >>>$org_name<<< in the ort cluster:"
		  	. scalar(@ort_genes) . "\n";
		  	$orts->print_cluster( $gene_name );
		  	getc();		  	
		}
		  

		my $draw_at_least_one = 0;

		foreach my $curr_ort_gene (@ort_genes) {

			print STDERR "Ort name: $curr_ort_gene\n"
			  if $debug_projections == 1;
			next if $curr_ort_gene eq "0";

			next if not defined( $gene_figure{$curr_ort_gene}{rectangle} );
			my @dst_rectangle = @{ $gene_figure{$curr_ort_gene}{rectangle} };

			print STDERR
			  "Dst rectangle: $dst_rectangle[ 0 ], $dst_rectangle[ 1 ]\n"
			  if $debug_projections == 1;

# Parameters are R,G,B and alpha for transparency.
# R,G,B goes from 0 to 255 (complete saturated)
# The alpha value may range from 0 (opaque) to 127 (transparent).
# The alphaBlending function changes the way this alpha channel affects the resulting image.
			$projection_fill_color =
			  $im->colorAllocateAlpha( 125, 125, 125, 110 );
			$projection_outline = $im->colorAllocateAlpha( 80, 80, 80, 80 );

			my $poly = new GD::SVG::Polygon;
			$poly->addPt( $src_rectangle[0], $src_rectangle[3] + 1 );
			$poly->addPt( $dst_rectangle[0], $dst_rectangle[1] - 1 );

			$poly->addPt( $dst_rectangle[2], $dst_rectangle[1] - 1 );
			$poly->addPt( $src_rectangle[2], $src_rectangle[3] + 1 );

			$im->filledPolygon( $poly, $projection_fill_color );
			$im->polygon( $poly, $projection_outline );

			getc() if $debug_projections == 1;

		}
	}
}

open OUTPUT, ">" . $outputFile;
print OUTPUT $im->svg;
close OUTPUT;

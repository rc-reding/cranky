include { MERGE_BARCODES } from './modules/barcode_processing.nf'

include {ASSEMBLE_ONT; MAP_REFERENCE; VARIANT_CALL_CLAIR3;
		 VARIANT_CALL_SAMTOOLS; VARIANT_CALL_LONGSHOT} from './modules/assembly.nf'

include {CALCULATE_COVERAGE_DEPTH; QC_ONT; FILTER_READS;
		 FIND_MLST; FIND_wgMLST} from './modules/assembly_utils.nf'

// Cannot call a process more than once, rename to re-use
include {FIND_MLST as FIND_MLST_REF} from './modules/assembly_utils.nf'
include {QC_ONT as FILTERED_QC_ONT} from './modules/assembly_utils.nf'

include {GENERATE_CONSENSUS; PREALIGN_GENOMES; DISTANCE_MATRIX;
		 DISTANCE_MATRIX as DISTANCE_MATRIX_wo_CTRL;
		 ALIGN_FROM_VCF; GENERATE_PHYLOGENY_GUBBINS; GENERATE_PHYLOGENY;
		 CORRECT_RECOMBINATION} from './modules/phylogeny.nf'

include {PLOT_PHYLOGENY; PLOT_QC; PLOT_COVERAGE} from './modules/plot_figures.nf'


// MERGE BARCODES IF MULTIPLE RUNS PER BARCODE

workflow preprocess_reads {
        data = Channel.from("$params.input")
        MERGE_BARCODES(data)
}

// MAIN WORKFLOW

workflow reference_mapping {
    // Input
    take:
		reads
		ref_genome

    // Function
    main:
		mapped_asmbl = MAP_REFERENCE(reads, ref_genome,
								     "$params.output/assembly")
		barcode_depth = CALCULATE_COVERAGE_DEPTH(mapped_asmbl.bam,
												 "$params.output/depth")
		PLOT_COVERAGE("$params.output/depth",
					  barcode_depth.cov.collect(),
				      "$params.output/depth")
    // Output
    emit:
    	asmbl = mapped_asmbl.bam
		coverage = barcode_depth.cov
		depth = barcode_depth.tsv
		depth_asmbl = barcode_depth.depth
}


workflow denovo_assembly{
	// Input
	take:
		reads

	// Main
	main:
		denovo = ASSEMBLE_ONT(reads, "$params.output/assembly/de_novo")

		asmbl = denovo.assembly.map { it -> it[1] }  // Get rid of 'barcode'

		// Compute depth using de novo assembly as reference
		mapped = reference_mapping(reads, asmbl)
									
	// Output
	emit:
		asmbl = denovo.assembly
		coverage = mapped.coverage
		depth = mapped.depth
		depth_asmbl = mapped.depth_asmbl
}


workflow quality_control{
    // Channels
    reads = Channel.fromPath("$params.input/reads_merged/*")
            .map { it -> tuple( it.baseName.replace('.fastq',''), it ) }

    // Function
    main:
		filtered = FILTER_READS(reads, "$params.output/reads_merged/filtered")
		control = QC_ONT(reads, "$params.output/qc")
		control_filtered = FILTERED_QC_ONT(filtered.reads,
										   "$params.output/qc/filtered")
		PLOT_QC("$params.output/qc",
				control_filtered.qc.collect(),
				"$params.output/qc")
    emit:
    	filtered.reads
}


workflow variant_caller {
	// Input
	take:
		asmbl
		ref_genome
		depth
		depth_asmbl

	// Main
	main:
        variants = VARIANT_CALL_CLAIR3(asmbl.join(depth_asmbl),
									   ref_genome, "$params.output/vcf")
        //variants = VARIANT_CALL_LONGSHOT(asmbl, ref_genome,
		//							     "$params.output/vcf")

		consensus = GENERATE_CONSENSUS(variants.vcf.join(depth).join(depth_asmbl),
								       asmbl.map{it -> it[1]}, ref_genome,
								       "$params.output/assembly/consensus")

	// Output
	emit:
		variants = variants.vcf
		consensus = consensus.fa
}


workflow MLST {
	// Input
	take:
		consensus
		ref_genome

	// Main
	main:
		if ( "$params.whole_genome" == true) {
			MLST = FIND_MLST(consensus, "$params.output/mlst/wg")
			// Whole-Genome MLST (wgMLST) - NOTE: Requires a reference genome in the DB
			wgMLST = FIND_wgMLST(consensus.collect(), "$params.output/assembly/de_novo",
								 ref_genome, "$params.output/mlst/wg/wgmlst")
			
		} else {
		    FIND_MLST_REF(tuple("reference", ref_genome), "$params.output/mlst")
			MLST = FIND_MLST(consensus, "$params.output/mlst")
			// Whole-Genome MLST (wgMLST) - NOTE: Requires a reference genome in the DB
			//wgMLST = FIND_wgMLST(consensus.collect(), "$params.output/assembly/consensus",
			//					 ref_genome, "$params.output/mlst/wgmlst")
		}

		// wgMLST.msa.collect()
		// wgMLST.msa_wo_controls.collect()
		MLST.mlst.collect()
							
	// Output
	emit:
		classic = MLST.mlst
		//wgMLST_msa = wgMLST.msa
		//wgMLST_msa_wo_controls = wgMLST.msa_wo_controls
}


workflow phylogeny {
	// Input
	take:
		consensus
		variants
		coverage
		ref_genome

	// Main
	main:
		MLST = MLST(consensus, ref_genome)
		// Assign outputs
		MLST_classic = MLST.classic
		//wgMLST_msa = MLST.wgMLST_msa  // Used to compute SNP distance
		//wgMLST_msa_wo_controls = MLST.wgMLST_msa_wo_controls  // Used to compute SNP distance


		// Gubbins allows SH-test, ClonalFrameML does not.
		// Length correction similar otherwise, with Gubbins being
		// a bit more conservative than ClonalFrameML
		if ( "$params.whole_genome" == true ) {
			dist = DISTANCE_MATRIX(wgMLST_msa.collect(), "$params.output/phylogeny/wg")
							
			dist_wo_controls = DISTANCE_MATRIX_wo_CTRL(wgMLST_msa_wo_controls.collect(),
													   "$params.output/phylogeny/wg/no_controls")
							
			// Phylogeny
			phy = GENERATE_PHYLOGENY_GUBBINS(wgMLST_msa.collect(),
											 "$params.output/phylogeny/wg")
			// Maximum-Likelihood (ML) method with RaxML
			//phy_ML = GENERATE_PHYLOGENY(wgMLST_msa.collect(), "$params.output/phylogeny/wg")
			//phy = CORRECT_RECOMBINATION(phy_ML.tree, wgMLST_msa.collect(),
			//						    "$params.output/phylogeny/wg")
			
			// SNP distance required for plotting, wait for completion
			dist.snp.collect()
			dist_wo_controls.snp.collect()

			// Plot
			PLOT_PHYLOGENY("$params.output/phylogeny/wg", phy.tree,
					       "$params.output/mlst/wg", "$params.output/phylogeny/wg")
		} else {
		  	rehead_variants = PREALIGN_GENOMES(coverage,
											   variants.collect(),
											   "$params.output/vcf",
										 	   "$params.output/vcf/rehead")
						
			// Alignment must wait for all reads to be pre-processed
			barcodes = rehead_variants.preprocessed.map{ it -> it[1] }.collect()
		
			msa = ALIGN_FROM_VCF(barcodes, ref_genome,
							     "$params.output/vcf/rehead",
							     "$params.output/alignments")

			dist = DISTANCE_MATRIX(msa.alignment.collect(), "$params.output/phylogeny")
							
			dist_wo_controls = DISTANCE_MATRIX_wo_CTRL(msa.alignment_wo_controls.collect(),
													   "$params.output/phylogeny/no_controls")
							
			// Phylogeny
			phy = GENERATE_PHYLOGENY_GUBBINS(msa.alignment.collect(),
											 "$params.output/phylogeny")
			// Maximum-Likelihood (ML) method with RaxML
			//phy_ML = GENERATE_PHYLOGENY(msa.alignment.collect(), "$params.output/phylogeny")
			//phy = CORRECT_RECOMBINATION(phy_ML.tree, msa.alignment.collect(),
			//						    "$params.output/phylogeny")
			
			// SNP distance required for plotting, wait for completion
			dist.snp.collect()
			dist_wo_controls.snp.collect()

			// Plot
			PLOT_PHYLOGENY("$params.output/phylogeny", phy.tree,
					       "$params.output/mlst", "$params.output/phylogeny")
		}


}


workflow outbreaker {
	// Channel
	    ref_genome = file("$params.reference")

	// Main
	main:
		filtered_reads = quality_control()

		if ( "$params.whole_genome" == true ) {
			denovo_assembly(filtered_reads)

			// Assign outputs - WAAAAY TOO HACKY... WTF NEXTFLOW?!
			consensus_seq = denovo_assembly.out.asmbl
			coverage = denovo_assembly.out.coverage
			depth = denovo_assembly.out.coverage

			// TODO: ue minimap2 / something to sort contig based on references
			variants = ''
		} else {
			reference_mapping(filtered_reads, ref_genome)

			// Assign outputs - WAAAAY TOO HACKY... WTF NEXTFLOW?!
			asmbl = reference_mapping.out.asmbl
			coverage = reference_mapping.out.coverage
			depth = reference_mapping.out.depth
			depth_asmbl = reference_mapping.out.depth_asmbl

			variant_caller(asmbl, ref_genome, depth, depth_asmbl)
			// Assign outputs
			variants  = variant_caller.out.variants
			consensus_seq  = variant_caller.out.consensus
		
		}
		
		phylogeny(consensus_seq, variants, coverage, ref_genome)
}

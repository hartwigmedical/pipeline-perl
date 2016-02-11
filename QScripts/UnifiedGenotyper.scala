package org.broadinstitute.gatk.queue.qscripts

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodsCalculationModel
import org.broadinstitute.gatk.tools.walkers.genotyper.OutputMode

class VariantCaller extends QScript {
    // Create an alias 'qscript' to be able to access variables in the VariantCaller.
    // 'qscript' is now the same as 'VariantCaller.this'
    qscript =>

    // Required arguments. All initialized to empty values.
    @Input(doc="The reference file for the bam files.", shortName="R", required=true)
    var referenceFile: File = _

    @Input(doc="One or more bam files.", shortName="I")
    var bamFiles: List[File] = Nil

    @Input(doc="Output core filename.", shortName="O", required=true)
    var out: File = _

    @Argument(doc="Maxmem.", shortName="mem", required=true)
    var maxMem: Int = _

    @Argument(doc="Number of cpu threads per data thread", shortName="nct", required=true)
    var numCPUThreads: Int = _

    @Argument(doc="Number of scatters", shortName="nsc", required=true)
    var numScatters: Int = _

    @Argument(doc="Genotype likelihoods calculation model to employ (SNP, INDEL or BOTH)", shortName="glm", required=true)
    var glm: String = _

    @Argument(doc="Make reference calls", shortName="rc", required=false)
    var refcalls: String = _

    @Argument(doc="Minimum phred-scaled confidence to call variants", shortName="stand_call_conf", required=true)
    var standCallConf: Int = _ //30 //default: best-practices value

    @Argument(doc="Minimum phred-scaled confidence to emit variants", shortName="stand_emit_conf", required=true)
    var standEmitConf: Int = _ //10 //default: best-practices value

    // The following arguments are all optional.
    @Input(doc="An optional file with known SNP sites.", shortName="D", required=false)
    var dbsnpFile: File = _

    @Input(doc="An optional file with targets intervals.", shortName="L", required=false)
    var targetFile: File = _

    @Argument(doc="Amount of padding (in bp) to add to each interval", shortName="ip", required=false)
    var intervalPadding: Int = 0

    def script() {
	val unifiedGenotyper = new UnifiedGenotyper

	// All required input
	unifiedGenotyper.input_file = bamFiles
	unifiedGenotyper.reference_sequence = referenceFile

	unifiedGenotyper.scatterCount = numScatters
	unifiedGenotyper.memoryLimit = maxMem
	unifiedGenotyper.num_cpu_threads_per_data_thread = numCPUThreads

	unifiedGenotyper.stand_emit_conf = standEmitConf
	unifiedGenotyper.stand_call_conf = standCallConf

	unifiedGenotyper.out = qscript.out + ".raw_variants.vcf"

	if(refcalls == "yes"){
	    unifiedGenotyper.output_mode = OutputMode.EMIT_ALL_CONFIDENT_SITES
	}
	


	//SNP INDEL or BOTH
	if (glm == "SNP") {
	    unifiedGenotyper.genotype_likelihoods_model = GenotypeLikelihoodsCalculationModel.Model.SNP
	} else if (glm == "INDEL") {
	    unifiedGenotyper.genotype_likelihoods_model = GenotypeLikelihoodsCalculationModel.Model.INDEL
	} else if (glm == "BOTH") {
	    unifiedGenotyper.genotype_likelihoods_model = GenotypeLikelihoodsCalculationModel.Model.BOTH
	}

	// Optional input
	if (dbsnpFile != null) {
	    unifiedGenotyper.D = dbsnpFile
	}
	if (targetFile != null) {
	    unifiedGenotyper.L :+= targetFile
	    unifiedGenotyper.ip = intervalPadding
	}

	//add function to queue
	add(unifiedGenotyper)
    }
}
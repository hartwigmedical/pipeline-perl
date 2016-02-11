package org.broadinstitute.gatk.queue.qscripts

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._

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

    @Argument(doc="Ploidy (number of chromosomes) per sample", shortName="ploidy", required=false)
    var samplePloidy: Int = 2

    def script() {
	val haplotypeCaller = new HaplotypeCaller

	// All required input
	haplotypeCaller.input_file = bamFiles
	haplotypeCaller.reference_sequence = referenceFile
	haplotypeCaller.out = qscript.out + ".raw_variants.vcf"

	haplotypeCaller.scatterCount = numScatters
	haplotypeCaller.memoryLimit = maxMem
	haplotypeCaller.num_cpu_threads_per_data_thread = numCPUThreads

	haplotypeCaller.stand_emit_conf = standEmitConf
	haplotypeCaller.stand_call_conf = standCallConf

	// Optional input
	if (dbsnpFile != null) {
	    haplotypeCaller.D = dbsnpFile
	}
	if (targetFile != null) {
	    haplotypeCaller.L :+= targetFile
	    haplotypeCaller.ip = intervalPadding
	}
	
	haplotypeCaller.sample_ploidy = samplePloidy
	
	//add function to queue
	add(haplotypeCaller)
    }
}
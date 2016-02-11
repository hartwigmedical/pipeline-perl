package org.broadinstitute.gatk.queue.qscripts

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._

import org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType

class VariantCaller extends QScript {
    // Create an alias 'qscript' to be able to access variables in the VariantCaller.
    // 'qscript' is now the same as 'VariantCaller.this'
    qscript =>

    // Required arguments. All initialized to empty values.
    @Input(doc="The reference file for the bam files.", shortName="R", required=true)
    var referenceFile: File = _

    @Input(doc="One or more input gVCF files", shortName="V")
    var gvcfFiles: List[File] = Nil

    @Input(doc="Output core filename.", shortName="O", required=true)
    var out: File = _

    @Argument(doc="Maxmem.", shortName="mem", required=true)
    var maxMem: Int = _

    @Argument(doc="Number of cpu threads per data thread", shortName="nct", required=true)
    var numCPUThreads: Int = _

    @Argument(doc="Number of scatters", shortName="nsc", required=true)
    var numScatters: Int = _

    def script() {
	val genotypeGVCFs = new GenotypeGVCFs
	
	genotypeGVCFs.V = gvcfFiles
	genotypeGVCFs.reference_sequence = referenceFile
	genotypeGVCFs.scatterCount = numScatters
	genotypeGVCFs.num_threads = numCPUThreads //for now use numCPUThreads, maybe change to new numDataThreads variable
	
	genotypeGVCFs.out = qscript.out + ".raw_variants.vcf"
	
	add(genotypeGVCFs)
    }
}

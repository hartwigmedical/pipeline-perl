package org.broadinstitute.gatk.queue.qscripts

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.extensions.cancer.MuTect

class SomaticCaller extends QScript {
    // Create an alias 'qscript' to be able to access variables in the SomaticCaller.
    // 'qscript' is now the same as 'SomaticCaller.this'
    qscript =>

    // Required arguments. All initialized to empty values.
    @Input(doc="The reference file for the bam files.", shortName="R", required=true)
    var referenceFile: File = _

    @Input(doc="input tumor BAM file - or list of BAM files", shortName="tb", required=true)
    var tumor_bam: File = _

    @Input(doc="input normal BAM file - or list of BAM files", shortName="nb", required=true)
    var normal_bam: File = _
    
    @Input(doc="Output core filename.", shortName="O", required=true)
    var out: File = _

    @Argument(doc="Maxmem.", shortName="mem", required=true)
    var maxMem: Int = _

    @Argument(doc="Number of scatters", shortName="nsc", required=true)
    var numScatters: Int = _

    @Input(doc="An optional file with known SNP sites.", shortName="D", required=false)
    var dbsnpFile: File = _
    
    @Input(doc="COSMIC sites to use (must be in VCF format)", shortName="C", required=false)
    var cosmicFile: File = _

    // The following arguments are all optional.
    @Input(doc="An optional file with targets intervals.", shortName="L", required=false)
    var targetFile: File = _

    @Argument(doc="Amount of padding (in bp) to add to each interval", shortName="ip", required=false)
    var intervalPadding: Int = 0

    def script() {
	val mutect = new MuTect

	// All required input
	mutect.input_file :+= new TaggedFile(tumor_bam, "tumor")
	mutect.input_file :+= new TaggedFile(normal_bam, "normal")
	mutect.reference_sequence = referenceFile
	mutect.out = qscript.out + ".call_stats.txt"
	mutect.vcf = qscript.out + ".vcf"

	mutect.scatterCount = numScatters
	mutect.memoryLimit = maxMem

	mutect.dbsnp :+= dbsnpFile
	mutect.cosmic :+= cosmicFile
	
	// Optional input
	if (targetFile != null) {
	    mutect.L :+= targetFile
	    mutect.ip = intervalPadding
	}
	
	//add function to queue
	add(mutect)
    }
}
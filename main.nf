#!/usr/bin/env/ nextflow

nextflow.enable.dsl = 2
def helpMessage() {
    log.info"""
    Nanopore Sequencing De Novo Assembly Pipeline

    Takes an input of raw reads produced by Oxford Nanopore sequencing, trims
    the reads using Porechop, and performs De Novo Assembly using the assembler
    Flye.

    USAGE: nextflow run main.nf [OPTIONS] [--inputFASTQs|--inputFASTAs] INPUT_DIR --output OUTPUT_DIR --model MEDAKA_MODEL
    OPTIONS:

    --inputFASTQs INPUT_DIR - [Required (Cannot be supplied with --inputFASTAs)] A directory containing FASTQ files produced by ONT Sequencing. The pipeline will assemble these files and then classify them.

    --inputFASTAs INPUT_DIR - [Required (Cannot be supplied with --inputFASTQs)] A directory containing already assembled FASTA files. In this case, the piepline will skip assembly and simply classify these sequences.

    --output OUTPUT_DIR - [Required] A directory to place output files (If not existing, pipeline will create)
    
    --model MEDAKA_MODEL - [Required if --inputFASTQs was supplied] Medaka requires a 'model' which corresponds to the flow-cell type/basecaller parameters used to correct for errors typical to that technology. Use 'medaka tools list_models' to find this.

    OPTIONAL:
        --trimONTAdapters - Enables ONT Adapter/Barcode trimming using Porechop [Default = off] (Not compatible with --inputFASTAs option)

        --minReadLen INT - If supplied, the pipeline will perform length filtering using Chopper excluding reads less than this size and sets the minimum length for Miniasm assembly [Default = 200 bp] (Disabled when --inputFASTAs option supplied)

        --maxReadLen INT - If supplied, the pipeline will perform legnth filtering using Chopper excluding reads greater than this size [Default = off] (Not compatible with --inputFASTAs option)

        --minReadsToAssemble INT - By default, the pipeline will assemble to the HA/NA segments with the most reads aligned to them. Supplying this parameter changes this to assemble reads to any segment with more than the supplied number of reads [Default = OFF] (Not compatible with --inputFASTAs option)

        --minCov INT - The minimum coverage below which a position will be masked [Default = 10] (Disabled when --inputFASTAs option supplied).
"""
}

// If the help parameter is supplied, link display the help message
// and quit the pipeline
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

def createSummaryHeader (inType, trim, minLenFilt, maxLenFilt) {
    
    if (inType == 'FASTQ') {

        FinalHeader = "Sample,Raw Reads,Average Raw Read Length,"

        if (trim != false) {
            FinalHeader = FinalHeader + "Read Post Porechop Trimming,Average Trimmed Read Length,"
        }

        if (minLenFilt != 0 || maxLenFilt != 0) {
            FinalHeader = FinalHeader + "Reads Post Length Filter,Average Filtered Read Length,"
        }

        FinalHeader = FinalHeader + "References Assembled To,Reads Aligned,Average Read Depth,SNPs,Indels,Masked Sites,Coverage,Typing"

    }
    else if (inType == "FASTA") {
        FinalHeader = "Sample,Number of Sequences,Sequence Coverage,Typing"
    }

    return FinalHeader
}

// Defines input parameters. Setting parameters
// to false by default allows us to check whether
// these parameters have been set by the user (useful
// for required parameters or binary parameters)
params.inputFASTQs = false
params.inputFASTAs = false
params.output = false
params.threads = 1
params.model = false
params.trimONTAdapters = false
params.minReadsToAssemble = false
params.minReadLen = 200
params.maxReadLen = 0
params.minCov = 10

include { Setup } from "./modules.nf"
include { QC_Report } from "./modules.nf"
include { QC_Report_Filtered } from "./modules.nf"
include { Collect_Raw_Read_Stats } from "./modules.nf"
include { Collect_FASTA_Stats } from "./modules.nf"
include { Porechop_Trimming } from "./modules.nf"
include { Length_Filtering } from "./modules.nf"
include { Select_References_For_Assembly } from "./modules.nf"
include { Minimap2_Alignment } from "./modules.nf"
include { Medaka_Alignment_Polish } from "./modules.nf"
include { Call_Variants } from "./modules.nf"
include { Generate_Consensus } from "./modules.nf"
include { Abricate_Typing } from "./modules.nf"
include { Write_Summary } from "./modules.nf"


// Checks the input parameters
inDir = ''
inputFiles_ch = ""
inputType = ""
if ((params.inputFASTAs == false) && (params.inputFASTQs == false)) {
    // If neither the --inputFASTQs option or --inputFASTAs option was supplied,
    // notify the user and exit
    println "ERROR: Either the --inputFASTQs or --inputFASTAs option needs to be supplied. Please revise the command arguments."
    exit(1)
}
else if ((params.inputFASTAs != false) && (params.inputFASTQs != false)) {
    // If both the --inputFASTQs option and --inputFASTAs option are supplied,
    // this is not allowed. Notify the user and exit.
    println "ERROR: Both the --inputFASTQs and --inputFASTAs options were supplied. Please only supply one option."
    exit(1)
}
else if (params.inputFASTQs != false) {
    
    if (!(file(params.inputFASTQs).isDirectory())) {
        // If the input directory is not set, notify the user and exit.
        println "ERROR: ${params.inputFASTQs} is not an existing directory."
        exit(1)
    }
    else {
    
        // If the parameter is set, convert the value provided to a file type
        // to get the absolute path, and then convert back to a string to be
        // used in the pipeline.
        inDir = file(params.inputFASTQs).toString()

        inputType = "FASTQ"

        println "Input Directory: ${inDir}\n"

        // Create a channel for hte input files.
        inputFiles_ch = Channel
        .fromPath("${inDir}/*.fastq*")
        // The .fromPath() function spits out each file, and we can use
        // the build-in nextflow functionality to get the name of the file
        // and place it into a tuple along with its corresponding file. That way
        // we can use the file name in the pipeline to name any output files. The
        // getSimpleName() function is used as opposed to getName(), as the
        // getName() function only removes the last .extension and would
        // cause issues if the fastqs were gzipped.
        .map { it -> [it.getSimpleName(), it]}
    }

}
else if (params.inputFASTAs != false) {

    if (!(file(params.inputFASTAs).isDirectory())) {
        // If the input directory is not set, notify the user and exit.
        println "ERROR: ${params.inputFASTAs} is not an existing directory."
        exit(1)
    }
    else {

        // If the parameter is set, convert the value provided to a file type
        // to get the absolute path, and then convert back to a string to be
        // used in the pipeline.
        inDir = file(params.inputFASTAs).toString()

        inputType = "FASTA"

        println "Input Directory: ${inDir}\n"

        // Create a channel for hte input files.
        inputFiles_ch = Channel
        .fromPath("${inDir}/*.fasta*")
        // The .fromPath() function spits out each file, and we can use
        // the build-in nextflow functionality to get the name of the file
        // and place it into a tuple along with its corresponding file. That way
        // we can use the file name in the pipeline to name any output files. The
        // getSimpleName() function is used as opposed to getName(), as the
        // getName() function only removes the last .extension and would
        // cause issues if the fastqs were gzipped.
        .map { it -> [it.getSimpleName(), it]}
    }
}

// Checks the output parameter.
outDir = ''
if (params.output == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No output directory provided. Pipeline requires an output directory."
    exit(1)
}
else {
    // If the parameter is set, ensure that the directory provided ends
    // in a trailing slash (to keep things consistent throughout) the
    // pipeline code.
    outDir = file(params.output).toString()
    println "Output Directory: ${outDir}\n"
}

// Checks the model parameter.
model = ''
if ((params.model == false) && (params.inputFASTQs != false)){
    // If the parameter is not set, notify the user and exit.
    //
    // This parameter requires a special input, so an example value is provided to the user
    // along with instructions to view the possible models using a medaka command.
    println "ERROR: no ONT model provided. Medaka requires a model in the format:"
    println ""
    println "{pore}_{device}_{caller variant}_{caller version}"
    println ""
    println "To see existing models enter: medaka tools list_models"
    exit(1)
}
else if ((params.model != false) && (params.inputFASTQs != false)){
    // If the parameter is provided, set the model variable equal to the
    // provided value.
    model = params.model
}
else if ((params.model == false) && (params.inputFASTAs != false)) {
    model = "DISABLED"
}
else if ((params.model != false) && (params.inputFASTAs != false)) {
    model = "DISABLED"
    println "NOTE: The --inputFASTAs option was supplied. Thus, the --model option will be ignored.\n"
}


// Handles the --trimONTAdapters parameter. Updates
// a variable to be written to the analysis parameters file
// if trimming was enabled.
trimmingEnabled = "DISABLED"
if ((params.trimONTAdapters != false) && (params.inputFASTQs != false)){
    trimmingEnabled = "ENABLED"
}
else if ((params.trimONTAdapters != false) && (params.inputFASTAs != false)) {
    println "NOTE: The --inputFASTAs option was supplied. Thus, the --trimONTAdapters option will be ignored.\n"
}


// Handles the --minReadLen parameter.
minReadLenVal = "DISABLED"
minReadLenParam = ""
if ((params.minReadLen != 0) && (params.inputFASTQs != false)) {
    // Updates a variable to be written to the analysis parameters file
    // if a minimum read length threshold was supplied.
    minReadLenVal = params.minReadLen + " bp"

    // As well if a value was supplied, it creates a variable
    // containing teh string to be passed to chopper. By doing so,
    // this allows us to handle the case when a minimum length
    // was supplied but not a maximum length (as the maximum length
    // parameter string will be blank and not affect the command).
    minReadLenParam = "-l ${params.minReadLen}"
}
else if ((params.minReadLen != 0) && (params.inputFASTAs != false)) {
    println "NOTE: The --inputFASTAs option was supplied. Thus, the minimum read length filter will be ignored.\n"
}


// Handles the --maxReadLen parameter.
maxReadLenVal = "DISABLED"
maxReadLenParam = ""
if ((params.maxReadLen != 0) && (params.inputFASTQs != false)) {

    // Updates a variable to be written to the analysis parameters file
    // if a maximum read length threshold was supplied.
    maxReadLenVal = params.maxReadLen + " bp"

    // As well if a value was supplied, it creates a variable
    // containing teh string to be passed to chopper. By doing so,
    // this allows us to handle the case when a minimum length
    // was supplied but not a maximum length (as the maximum length
    // parameter string will be blank and not affect the command).
    maxReadLenParam = "--maxlength ${params.maxReadLen}"
}
else if ((params.maxReadLen != 0) && (params.inputFASTAs != false)) {
    println "NOTE: The --inputFASTAs option was supplied. Thus, the maximum read length filter will be ignored.\n"
}

minReadsToAssembleVal = "Assemble to Segments with Highest Number of Reads"
minReadsToAssembleParam = ""
if ((params.minReadsToAssemble != false) && (params.inputFASTQs != false)) {
    minReadsToAssembleVal = "Assebmle any segments with more than ${params.minReadsToAssemble} reads aligned"
    minReadsToAssembleParam = "--minReads ${params.minReadsToAssemble}"
}
else if ((params.minReadsToAssemble != false) && (params.inputFASTAs != false)) {
    println "NOTE: The --inputFASTAs option was supplied. Thus, the --minReadsToAssemble option will be ignored.\n"
    minReadsToAssembleVal = "DISABLED"
}
else if ((params.minReadsToAssemble == false) && (params.inputFASTAs != false)) {
    minReadsToAssembleVal = "DISABLED"
}

// Generates the summary header based on the parameters supplied. Allows for the header
// to be dynamic depending on the statistics collected.
summaryHeader = createSummaryHeader(inputType, params.trimONTAdapters, params.minReadLen, params.maxReadLen)

workflow {

    Setup( inputType, trimmingEnabled, minReadLenVal, maxReadLenVal, model, minReadsToAssembleVal, summaryHeader, outDir )

    if (params.inputFASTQs != false) {

        QC_Report( inputFiles_ch, outDir )

        Collect_Raw_Read_Stats( inputFiles_ch )

        if (params.trimONTAdapters != false && (params.minReadLen != 0 || params.maxReadLen != 0)) {
            Porechop_Trimming( Collect_Raw_Read_Stats.out[0], outDir, Collect_Raw_Read_Stats.out[1] )

            Length_Filtering( Porechop_Trimming.out[0], minReadLenParam, maxReadLenParam, outDir, Porechop_Trimming.out[2] )

            QC_Report_Filtered( Length_Filtering.out[0], outDir )

            // Performs De Novo Assembly of the filtered reads using Miniasm
            Select_References_For_Assembly( Length_Filtering.out[0], minReadsToAssembleParam, baseDir, outDir, Length_Filtering.out[1] )

        }
        else if (params.trimONTAdapters != false && params.minReadLen == 0 && params.maxReadLen == 0) {

            Porechop_Trimming( Collect_Raw_Read_Stats.out[0], outDir, Collect_Raw_Read_Stats.out[1] )

            QC_Report_Filtered( Porechop_Trimming.out[0], outDir )

            // Performs De Novo Assembly of the trimmed reads using Miniasm
            Select_References_For_Assembly( Porechop_Trimming.out[0], minReadsToAssembleParam, baseDir, outDir, Porechop_Trimming.out[2] )
        }
        else if (params.trimONTAdapters == false && (params.minReadLen != 0 || params.maxReadLen != 0)) {
            Length_Filtering( Collect_Raw_Read_Stats.out[0], minReadLenParam, maxReadLenParam, outDir, Collect_Raw_Read_Stats.out[1] )

            QC_Report_Filtered( Length_Filtering.out[0], outDir )

            // Performs De Novo Assembly of the filtered reads using Flye
            Select_References_For_Assembly( Length_Filtering.out[0], minReadsToAssembleParam, baseDir, outDir, Length_Filtering.out[1] )
        }
        else {
            // Performs De Novo Assembly of the reads using Flye
            Select_References_For_Assembly( Collect_Raw_Read_Stats.out[0], minReadsToAssembleParam, baseDir, outDir, Collect_Raw_Read_Stats.out[1] )
        }

        Minimap2_Alignment( Select_References_For_Assembly.out[0], outDir, Select_References_For_Assembly.out[2] )

        Medaka_Alignment_Polish( Minimap2_Alignment.out[0], model, outDir, Minimap2_Alignment.out[1] )

        Call_Variants( Medaka_Alignment_Polish.out[0], params.minCov, outDir, Medaka_Alignment_Polish.out[1] )

        Generate_Consensus( Call_Variants.out[0], baseDir, outDir, params.minCov, Call_Variants.out[1])

        // Corrects the assembled contigs using Medaka
        //Medaka_Correct( Miniasm_Assembly.out[0], outDir, params.threads, model, Miniasm_Assembly.out[1] )

        Abricate_Typing( Generate_Consensus.out[0], baseDir, outDir, Generate_Consensus.out[2])
    }
    else if (params.inputFASTAs) {

        Collect_FASTA_Stats(inputFiles_ch, baseDir)

        Abricate_Typing( Collect_FASTA_Stats.out[0], baseDir, outDir, Collect_FASTA_Stats.out[1] )
    }

    Write_Summary( Abricate_Typing.out[1], outDir )

}
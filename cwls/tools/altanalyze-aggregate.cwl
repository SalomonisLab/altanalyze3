cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: DockerRequirement
  dockerPull: altanalyze:latest
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_prefix = function(ext) {
        if (inputs.output_prefix) {
          return inputs.output_prefix+ext;
        }
        return "results"+ext;
    };


inputs:

  junction_bed_file:
    type:
      - "null"
      - File
      - type: array
        items: File
    inputBinding:
      prefix: "--juncounts"
    doc: |
      Path the junction counts files. Coordinates are treated as 0-based.
      If provided with --intcounts, the number and the order should
      correspond to --intcounts.

  intron_bed_file:
    type:
      - "null"
      - File
      - type: array
        items: File
    inputBinding:
      prefix: "--intcounts"
    doc: |
      Path the intron counts files. Coordinates are treated as 0-based.
      If provided with --juncounts, the number and the order should
      correspond to --juncounts.

  aliases:
    type:
      - "null"
      - string
      - type: array
        items: string
    inputBinding:
      prefix: "--aliases"
    doc: |
      Column names to be used for the loaded counts. The number of provided aliases
      should be equal to the number of --intcounts and/or --juncounts files.
      Default: rootname of --intcounts and/or --intcounts files.

  reference_file:
    type: File?
    inputBinding:
      prefix: "--ref"
    doc: |
      Path to the gene model reference file. Coordinates are treated as 1-based.
      Required if --juncounts parameter was provided.

  chrom_list:
    type:
      - "null"
      - string
      - type: array
        items: string
    inputBinding:
      prefix: "--chr"
    doc: |
      Select chromosomes to process.
      Default: only main chromosomes

  export_to_bed:
    type: boolean?
    inputBinding:
      prefix: "--bed"
    doc: |
      Export annotated coordinates as BED file.
      Default: False

  log_level:
    type:
    - "null"
    - type: enum
      symbols:
      - "fatal"
      - "error"
      - "warning"
      - "info"
      - "debug"
    inputBinding:
      prefix: "--loglevel"
    doc: |
      Logging level.
      Default: info

  threads:
    type: int?
    inputBinding:
      prefix: "--threads"
    doc: |
      Number of threads to run in parallel where applicable.
      Default: 1

  processes:
    type: int?
    inputBinding:
      prefix: "--cpus"
    doc: |
      Number of processes to run in parallel where applicable.
      Default: 1

  output_prefix:
    type: string?
    default: ""
    inputBinding:
      prefix: "--output"
      valueFrom: $(get_output_prefix("_agg"))
    doc: |
      Output prefix.
      Default: results


outputs:

  aggregated_bed_file:
    type: File?
    outputBinding:
      glob: $(get_output_prefix("_agg.bed.gz"))
    secondaryFiles:
    - .tbi

  aggregated_h5ad_file:
    type: File
    outputBinding:
      glob: $(get_output_prefix("_agg.h5ad"))

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["altanalyze3", "aggregate"]

stdout: $(get_output_prefix("_agg_stdout.log"))
stderr: $(get_output_prefix("_agg_stderr.log"))


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "AltAnalyze3 Aggregate"
s:name: "AltAnalyze3 Aggregate"
s:alternateName: "AltAnalyze3 Aggregate"

s:downloadUrl: https://raw.githubusercontent.com/SalomonisLab/altanalyze3/master/cwls/tools/altanalyze-aggregate.cwl
s:codeRepository: https://github.com/SalomonisLab/altanalyze3
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  AltAnalyze3 Aggregate


s:about: |
  usage: altanalyze3 aggregate [-h] [--juncounts [JUNCOUNTS [JUNCOUNTS ...]]] [--intcounts [INTCOUNTS [INTCOUNTS ...]]]
                               [--aliases [ALIASES [ALIASES ...]]] [--ref REF] [--chr [CHR [CHR ...]]] [--bed]
                               [--loglevel {fatal,error,warning,info,debug}] [--threads THREADS] [--cpus CPUS] [--tmp TMP]
                               [--output OUTPUT]

  optional arguments:
    -h, --help            show this help message and exit
    --juncounts [JUNCOUNTS [JUNCOUNTS ...]]
                          Path the junction counts files. Coordinates are treated as 0-based. If provided with --intcounts,
                          the number and the order should correspond to --intcounts.
    --intcounts [INTCOUNTS [INTCOUNTS ...]]
                          Path the intron counts files. Coordinates are treated as 0-based. If provided with --juncounts,
                          the number and the order should correspond to --juncounts.
    --aliases [ALIASES [ALIASES ...]]
                          Column names to be used for the loaded counts. The number of provided aliases should be equal to
                          the number of --intcounts and/or --juncounts files.
                          Default: rootname of --intcounts and/or --intcounts files.
    --ref REF             Path to the gene model reference file. Coordinates are treated as 1-based. Required if --juncounts
                          parameter was provided.
    --chr [CHR [CHR ...]]
                          Select chromosomes to process. Default: only main chromosomes
    --bed                 Export annotated coordinates as BED file. Default: False
    --loglevel {fatal,error,warning,info,debug}
                          Logging level. Default: info
    --threads THREADS     Number of threads to run in parallel where applicable. Default: 1
    --cpus CPUS           Number of processes to run in parallel where applicable. Default: 1
    --tmp TMP             Temporary files location. Default: tmp
    --output OUTPUT       Output prefix. Default: results
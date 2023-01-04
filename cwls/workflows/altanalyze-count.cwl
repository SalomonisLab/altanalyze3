cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement


inputs:

  alignment_file:
    type: File[]
    secondaryFiles:
    - .bai
    doc: |
      Path to the coordinate-sorted indexed BAM file

  reference_file:
    type: File
    doc: |
      Path to the gene model reference file. Coordinates are treated as 1-based.

  overlap_bp:
    type: int?
    doc: |
      5' and 3' overlap that read should have over a splice-site to be counted.
      Default: 10

  strandness:
    type:
    - "null"
    - type: enum
      symbols:
      - "auto"
      - "forward"
      - "reverse"
      - "unstranded"
    doc: |
      Strand specificty of the RNA library.'unstranded' - reads from the
      left-most end of the fragment (in transcript coordinates) map to the
      transcript strand, and the right-most end maps to the opposite strand.
      'forward' - same as 'unstranded' except we enforce the rule that the
      left-most end of the fragment (in transcript coordinates) is the first
      sequenced (or only sequenced for single-end reads). Equivalently, it is
      assumed that only the strand generated during second strand synthesis
      is sequenced. Used for Ligation and Standard SOLiD. 'reverse' - same as
      'unstranded' except we enforce the rule that the right-most end of the
      fragment (in transcript coordinates) is the first sequenced (or only
      sequenced for single-end reads). Equivalently, it is assumed that only
      the strand generated during first strand synthesis is sequenced. Used
      for dUTP, NSR, and NNSR. Default: first 'auto' (try to detect strand
      from the XS tag of the read), then downgrade to 'unstranded'.
      Default: auto

  chrom_list:
    type:
      - "null"
      - string
      - type: array
        items: string
    doc: |
      Select chromosomes to process.
      Default: only main chromosomes

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
    doc: |
      Logging level.
      Default: info

  threads:
    type: int?
    doc: |
      Number of threads to run in parallel where applicable.
      Default: 1

  processes:
    type: int?
    doc: |
      Number of processes to run in parallel where applicable.
      Default: 1


outputs:

  intron_bed_file:
    type: File[]
    outputSource: intcount/intron_bed_file
    doc: "Intron counts BED file"

  junction_bed_file:
    type: File[]
    outputSource: juncount/junction_bed_file
    doc: "Junction counts BED file"

  aggregated_bed_file:
    type: File
    outputSource: aggregate/aggregated_bed_file
    doc: "Aggergated counts BED file"

  aggregated_h5ad_file:
    type: File
    outputSource: aggregate/aggregated_h5ad_file
    doc: "Aggergated counts H5AD file"

  intcount_stdout_log:
    type: File[]
    outputSource: intcount/stdout_log
    doc: "Stdout log from intcount step"

  intcount_stderr_log:
    type: File[]
    outputSource: intcount/stderr_log
    doc: "Stderr log from intcount step"

  juncount_stdout_log:
    type: File[]
    outputSource: juncount/stdout_log
    doc: "Stdout log from juncount step"

  juncount_stderr_log:
    type: File[]
    outputSource: juncount/stderr_log
    doc: "Stderr log from juncount step"

  aggregate_stdout_log:
    type: File
    outputSource: aggregate/stdout_log
    doc: "Stdout log from aggregate step"

  aggregate_stderr_log:
    type: File
    outputSource: aggregate/stderr_log
    doc: "Stderr log from aggregate step"


steps:

  intcount:
    run: ../tools/altanalyze-intcount.cwl
    in:
      alignment_file: alignment_file
      reference_file: reference_file
      overlap_bp: overlap_bp
      strandness: strandness
      chrom_list: chrom_list
      log_level: log_level
      threads: threads
      processes: processes
    scatter:
    - alignment_file
    out:
    - intron_bed_file
    - stdout_log
    - stderr_log
    
  juncount:
    run: ../tools/altanalyze-juncount.cwl
    in:
      alignment_file: alignment_file
      chrom_list: chrom_list
      log_level: log_level
      threads: threads
      processes: processes
    scatter:
    - alignment_file
    out:
    - junction_bed_file
    - stdout_log
    - stderr_log

  aggregate:
    run: ../tools/altanalyze-aggregate.cwl
    in:
      junction_bed_file: juncount/junction_bed_file
      intron_bed_file: intcount/intron_bed_file
      reference_file: reference_file
      chrom_list: chrom_list
      export_to_bed:
        default: true
      log_level: log_level
      threads: threads
      processes: processes
    out:
    - aggregated_bed_file
    - aggregated_h5ad_file
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "AltAnalyze3 Count"
s:name: "AltAnalyze3 Count"
s:alternateName: "AltAnalyze3 Count"

s:downloadUrl: https://raw.githubusercontent.com/SalomonisLab/altanalyze3/master/cwls/workflows/altanalyze-count.cwl
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
  AltAnalyze3 Count
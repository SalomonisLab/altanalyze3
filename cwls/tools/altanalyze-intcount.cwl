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
        var root = inputs.alignment_file.basename.split('.').slice(0,-1).join('.');
        return (root == "")?inputs.alignment_file.basename+ext:root+ext;
    };
- class: InitialWorkDirRequirement
  listing: |
    ${
      return [
        {
          "entry": inputs.alignment_file,
          "entryname": inputs.alignment_file.basename,
          "writable": true
        }
      ]
    }



inputs:

  alignment_file:
    type: File
    secondaryFiles:
    - .bai
    inputBinding:
      prefix: "--bam"
    doc: |
      Path to the coordinate-sorted indexed BAM file

  reference_file:
    type: File
    inputBinding:
      prefix: "--ref"
    doc: |
      Path to the gene model reference file. Coordinates are treated as 1-based.

  overlap_bp:
    type: int?
    inputBinding:
      prefix: "--span"
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
    inputBinding:
      prefix: "--strandness"
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
    inputBinding:
      prefix: "--chr"
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
      valueFrom: $(get_output_prefix("_int"))
    doc: |
      Output prefix.
      Default: results


outputs:

  intron_bed_file:
    type: File
    outputBinding:
      glob: "*.bed"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["altanalyze3", "intcount"]

stdout: $(get_output_prefix("_int_stdout.log"))
stderr: $(get_output_prefix("_int_stderr.log"))


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "AltAnalyze3 Intron Count"
s:name: "AltAnalyze3 Intron Count"
s:alternateName: "AltAnalyze3 Intron Count"

s:downloadUrl: https://raw.githubusercontent.com/SalomonisLab/altanalyze3/master/cwls/tools/altanalyze-intcount.cwl
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
  AltAnalyze3 Intron Count


s:about: |
  usage: altanalyze3 intcount [-h] --bam BAM --ref REF [--span SPAN]
                              [--strandness {auto,forward,reverse,unstranded}]
                              [--chr [CHR [CHR ...]]] [--savereads]
                              [--loglevel {fatal,error,warning,info,debug}] [--threads THREADS]
                              [--cpus CPUS] [--tmp TMP] [--output OUTPUT]

  optional arguments:
    -h, --help            show this help message and exit
    --bam BAM             Path to the coordinate-sorted indexed BAM file
    --ref REF             Path to the gene model reference file. Coordinates are treated as 1-based.
    --span SPAN           5' and 3' overlap that read should have over a splice-site to be
                          counted
    --strandness {auto,forward,reverse,unstranded}
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
                          from the XS tag of the read), then downgrade to 'unstranded'
    --chr [CHR [CHR ...]]
                          Select chromosomes to process. Default: only main chromosomes
    --savereads           Export processed reads into the BAM file. Default: False
    --loglevel {fatal,error,warning,info,debug}
                          Logging level. Default: info
    --threads THREADS     Number of threads to run in parallel where applicable. Default: 1
    --cpus CPUS           Number of processes to run in parallel where applicable. Default: 1
    --tmp TMP             Temporary files location. Default: tmp
    --output OUTPUT       Output prefix. Default: results
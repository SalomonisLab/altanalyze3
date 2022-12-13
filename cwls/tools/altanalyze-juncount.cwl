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
      valueFrom: $(get_output_prefix("_jun"))
    doc: |
      Output prefix.
      Default: results


outputs:

  junction_bed_file:
    type: File
    outputBinding:
      glob: "*.bed"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["altanalyze3", "juncount"]

stdout: $(get_output_prefix("_jun_stdout.log"))
stderr: $(get_output_prefix("_jun_stderr.log"))


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "AltAnalyze3 Junction Count"
s:name: "AltAnalyze3 Junction Count"
s:alternateName: "AltAnalyze3 Junction Count"

s:downloadUrl: https://raw.githubusercontent.com/SalomonisLab/altanalyze3/master/cwls/tools/altanalyze-juncount.cwl
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
  AltAnalyze3 Junction Count


s:about: |
  usage: altanalyze3 juncount [-h] --bam BAM [--chr [CHR [CHR ...]]] [--savereads]
                              [--loglevel {fatal,error,warning,info,debug}]
                              [--threads THREADS] [--cpus CPUS] [--tmp TMP]
                              [--output OUTPUT]

  optional arguments:
    -h, --help            show this help message and exit
    --bam BAM             Path to the coordinate-sorted indexed BAM file
    --chr [CHR [CHR ...]]
                          Select chromosomes to process. Default: only main chromosomes
    --savereads           Export processed reads into the BAM file. Default: False
    --loglevel {fatal,error,warning,info,debug}
                          Logging level. Default: info
    --threads THREADS     Number of threads to run in parallel where applicable. Default: 1
    --cpus CPUS           Number of processes to run in parallel where applicable. Default: 1
    --tmp TMP             Temporary files location. Default: tmp
    --output OUTPUT       Output prefix. Default: results
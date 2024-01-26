Design
------

Objects
=======

cwl_dependency_aliases
~~~~~~~~~~~~~~~~~~~~~~

Alias for importing subworkflows (a mapping table between step name in subworkflow and actual name (e.g. CWL ID or Nextflow's Process name).


dag.input_sources
~~~~~~~~~~~~~~~~~

Structure:

    <step_name>: {
        "in_emit_dict": {
            <cwl_id>: <output_element>
            ...
        }
        "in_sink_dict": {
            <step_name>: {


"in_emit_dict" describes list of outputs from predecessors which become inputs for the step. 
"in_sink_dict" describes list of inputs to the step, along with their source (i.e. predecessor steps). 
This is what we mean by caller's & callee's perspective in the paper. 

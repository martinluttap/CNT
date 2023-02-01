from pathlib import Path
from string import Template
from typing import List

import os
import sys

# import json5

from cwl_utils.parser import load_document_by_uri, save

class CWLFile:
    def __init__(self, cwl_path) -> None:
        self.cwl_path = cwl_path
        self.cwl_object = None
        self.cwl_class = None
        self.cwl_inputs = None # List
        self.cwl_outputs = None
        self.cwl_steps = None
        self.cwl_id = None
        self.container_uri = None
        self.cwl_name = None
        self.cwl_dependencies = None
        self.cwl_dependency_aliases = None
        self.cwl_basecommand = None

        # Initialization methods
        self.get_cwl_object(cwl_path)
        self.get_id()
        self.get_container_uri()
        self.get_name()
        self.get_class()
        self.get_inputs()
        self.get_outputs()
        self.get_steps()
        self.get_dependencies()
        self.get_dependency_aliases()
        self.get_basecommand()
        # self.get_id()

    def get_cwl_object(self, cwl_path): 
        # File Input - This is the only thing you will need to adjust or take in as an input to your function:
        # cwl_file = Path("./gdc_dnaseq.bamfastq_align.workflow.cwl")  # or a plain string works as well
        cwl_file = Path(cwl_path)

        # Import CWL Object
        cwl_obj = load_document_by_uri(cwl_file)

        # Export CWL Object into a built-in typed object
        saved_obj = save(cwl_obj)

        self.cwl_object = saved_obj

    def get_id(self):
        self.cwl_id = self.cwl_object['id']

    def get_container_uri(self):
        container_uri = None
        if "requirements" in self.cwl_object:
            for req in self.cwl_object['requirements']:
                if req['class'] == "DockerRequirement":
                    container_uri = req['dockerPull']
        
        self.container_uri = container_uri

    def get_name(self):
        self.cwl_name = self.cwl_object['id'].split('#')[-1]

    def get_class(self):
        self.cwl_class = self.cwl_object['class']

    def get_inputs(self):
        self.cwl_inputs = self.cwl_object['inputs']

    def get_outputs(self):
        self.cwl_outputs = self.cwl_object['outputs']

    def get_steps(self):
        if 'steps' in self.cwl_object:
            self.cwl_steps = self.cwl_object['steps']
        else:
            assert(self.cwl_object['class'] != "Workflow")

    def get_dependencies(self):
        if self.cwl_steps:
            self.cwl_dependencies = [e['run'] for e in self.cwl_steps]

    def get_dependency_aliases(self):
        if self.cwl_steps:
            dep_aliases = {}
            for e in self.cwl_steps:
                dep_aliases[e['run']] = e['id']

            self.cwl_dependency_aliases = dep_aliases

    def get_basecommand(self):
        if "baseCommand" in self.cwl_object:
            self.cwl_basecommand = self.cwl_object["baseCommand"]

    def print_cwl_object(self):
        print(self.cwl_object)

    def print_class(self):
        print(self.cwl_object['class'])

    def print_inputs(self):
        print(self.cwl_object['inputs'])

        for e in self.cwl_inputs:
            print(e)

    def print_id(self):
        print(self.cwl_object['id'])

    def print_outputs(self):
        outputs = self.construct_outputs()
        print(outputs)

    def print_steps(self):
        print(self.cwl_steps)

    def print_dependencies(self):
        for e in self.cwl_dependencies:
            print(e)

    def print_body(self):
        if self.cwl_class == "Workflow":
            self.print_template_workflow()
        else:
            self.print_template_process()

    def print_template_process(self):
        file_dir = os.path.dirname(os.path.realpath(__file__))
        in_template = open(file_dir + "/templates/template_process.txt", "r")
        template = Template(in_template.read())

        pr_input_list = [
            e for e in self.cwl_inputs
        ]

        pr_output_list = [
            e for e in self.cwl_outputs
        ]

        
        script = self.construct_script()
        # pr_output = self.construct_outputs()
        # pr_script_list = ' '.join(self.cwl_basecommand) if self.cw

        pr_input = self.format_input_list(pr_input_list)
        pr_output = self.format_output_list(pr_output_list)
        pr_script = """
        \"\"\"
        {}
        \"\"\"
        """.format(script)
        
        subs = {
            "process_id": str.upper(self.cwl_name),
            "container": "",  # By default we let container empty
            "input": pr_input,
            "output": pr_output,
            "script": pr_script
        }
        if self.container_uri:
            subs['container'] = "container \"{}\"\n".format(self.container_uri)
        
        result = template.substitute(subs) 

        print(result)

    def construct_script(self):
        ret_script = ' '.join(self.cwl_basecommand)
        # If inputs does not have position binding, 
        # we assign the same priority using a big value.
        for e in self.cwl_inputs:
            if "inputBinding" in e:
                if "position" in e:
                    pass
            else:
                e['inputBinding'] = {
                    'position': 9999,
                }
        # Here we construct arguments for basecommand with the right order
        # sorted_inputs = sorted(self.cwl_inputs, key=lambda x: x['inputBinding']['position'])
        sorted_inputs = self.cwl_inputs
        
        arguments = [' ']
        for e in sorted_inputs:
            assert('inputBinding' in e)
            if 'prefix' in e['inputBinding']:
                arguments.append(e['inputBinding']['prefix'])
            arguments.append("${{{}}}".format(e['id']))
        
        ret_script += ' '.join(arguments)

        return ret_script

    def construct_outputs(self):
        ret_outputs = []
        # for e in self.cwl_outputs:
        #     output = {
        #         'id': e['id'],
        #         'glob': e['outputBinding']['glob'] if 'glob' in e['outputBinding'] else None,
        #         'type': 


        return ret_outputs

    def print_template_workflow(self):

        in_template = open("{}/templates/template_workflow.txt".format(os.getcwd()), "r")
        template = Template(in_template.read())
        
        wf_take_list = [
            e['id'] for e in self.cwl_inputs
        ]
        wf_take_bindings = []

        wf_main_list = [
            "{}()".format(e['id'].upper()) for e in self.cwl_steps
        ]
        wf_main_bindings = []

        wf_emit_list = [
            e['id'] for e in self.cwl_outputs
        ]
        wf_emit_bindings = []

        wf_take = self.format_list(wf_take_list, wf_take_bindings)
        wf_main = self.format_list(wf_main_list, wf_main_bindings)
        wf_emit = self.format_list(wf_emit_list, wf_emit_bindings)

        subs = {
            "workflow_id": str.upper(self.cwl_name),
            "take": wf_take,
            "main": wf_main,
            "emit": wf_emit
        }

        result = template.substitute(subs) 

        print(result)

    def print_template_header(self):
        for e in self.cwl_dependencies:
            lib_name = e.strip().split('/')[-1].strip('.cwl')
            lib_alias = self.cwl_dependency_aliases[e]
            lib_source = e.replace(".cwl", ".nf")
            print("INCLUDE {{ {} as {} }} from {}".format(
                lib_name.upper(),
                lib_alias.upper(),
                lib_source
            ))

    def format_input_list(self, inp_list: List) -> str:
        ret_str = ""
        for e in inp_list:
            ret_str += self.get_inp_type(e) + " {}\n\t".format(e['id']) 
        
        return ret_str

    def get_inp_type(self, inp):
        inp_type = None
        if self.cwl_class == 'Workflow':
            inp_type = ''
        else: # CommandLineTool or ExpressionTool
            inp_type = "val"

        return inp_type


    def get_out_type(self, out):
        out_type = None
        if self.cwl_class == 'Workflow':
            out_type = ''
        else: # CommandLineTool or ExpressionTool
            if out['type'] == "File" or "*" in out['outputBinding']['glob']:
                out_type = "path"
            else:
                out_type = "val"

        return out_type

    def format_output_list(self, out_list: List) -> str:        
        ret_str = ""
        for e in out_list:
            ret_str += self.get_out_type(e) + " \"{}\"\t, emit: {}\n\t".format(e['outputBinding']['glob'], e['id']) 
        
        return ret_str


    def format_list(self, inp_list: List, bindings: List) -> str:
        inp_type = '' if self.cwl_class == 'Workflow' else 'val '
        
        ret_str = inp_type + inp_list[0] + "({})".format(bindings[0]) if len(bindings) else ""

        if len(inp_list) > 1:
            assert len(inp_list) == len(bindings)
            for e, b in zip(inp_list[1:], bindings[1:]):
                ret_str += "\n\t\t" + inp_type + e + "({})".format(b)

        return ret_str

    def print_script(self):
        # self.print_header()
        self.print_body()

def convert_to_nextflow(cwl_path):
    # File Input - This is the only thing you will need to adjust or take in as an input to your function:
    # cwl_file = Path("./gdc_dnaseq.bamfastq_align.workflow.cwl")  # or a plain string works as well
    cwl_file = Path(cwl_path)

    # Import CWL Object
    cwl_obj = load_document_by_uri(cwl_file)

    # Export CWL Object into a built-in typed object
    saved_obj = save(cwl_obj)

    print(str_obj)

if __name__ == "__main__":
    # print(sys.argv[1])
    cwl_path = sys.argv[1]

    cwl_file = CWLFile(cwl_path)
    if cwl_file.cwl_class == "Workflow":
        cwl_file.print_steps()
        cwl_file.print_dependencies()
        cwl_file.print_template_header()
        cwl_file.print_body()
    else:
        cwl_file.print_body()

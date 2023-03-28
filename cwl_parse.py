from argparse import ArgumentParser
from itertools import count
from pathlib import Path
from string import Template
from typing import List

import ast 
import json
import logging
import os
import re
import subprocess
import sys
import yaml

# Local libraries
from libs.DAGDrawer import DagDrawer

# Third-party libraries
from cwl_utils.parser import load_document_by_uri, save
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

# Logging
logger = logging.getLogger(__name__)
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.DEBUG)

# Argument parsing
parser = ArgumentParser(description="Automated CWL <-> X Translator")
parser.add_argument("cwl_source", help="Source CWL file")
parser.add_argument(
    "--import_yaml",
    action="store_true",
    help="Read YAML import_yaml",
)
parser.add_argument(
    "--inspect_docker",
    action="store_true",
    help="Optimize translation by allowing inspection of docker images. \
        WARNING: This option will pull the docker images and does not clean them up.",
)
parser.add_argument(
    "--unit_test",
    action="store_true",
    help="Append a Workflow(...) call for stand-alone execution.",
)
parser.add_argument(
    "--with_tool_config",
    action="store_true",
    help="Generate <cwl_source>.nextflow.config for docker usage.",
)
parser.add_argument(
    "--show_dag",
    action="store_true",
    help="Print the DAG of workflow.",
)
parser.add_argument(
    "--with_wf_config",
    action="store_true",
    help="Generate <cwl_source>.nextflow.config with all tool configs associated with this cwl.",
)
parser.set_defaults(inspect_docker=False)
parser.set_defaults(unit_test=False)
parser.set_defaults(with_tool_config=False)
args = parser.parse_args()
class CWLFile:
    def __init__(self, cwl_path) -> None:
        self.cwl_path = cwl_path
        self.cwl_object = None
        self.cwl_class = None
        self.cwl_inputs = None  # List
        self.cwl_outputs = None
        self.cwl_steps = None
        self.cwl_id = None
        self.container_uri = None
        self.cwl_name = None
        self.cwl_dependencies = None
        self.cwl_dependency_aliases = None
        self.cwl_basecommand = None
        self.cwl_arguments = None
        self.cwl_stdout = None
        self.param_dict = {}
        self.step_dep = {}
        self.dag = None
        self.step_files = {}

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
        self.get_arguments()
        self.get_stdout()
        if self.cwl_class == "Workflow":
            self.dag = nx.DiGraph()
            self.get_dag()
        # self.get_id()

    def get_cwl_object(self, cwl_path):
        # File Input - This is the only thing you will need to adjust or take in
        # as an input to your function: cwl_file =
        # Path("./gdc_dnaseq.bamfastq_align.workflow.cwl")  # or a plain string
        # works as well
        cwl_file = Path(cwl_path)

        # Import CWL Object
        cwl_obj = load_document_by_uri(cwl_file)

        # Export CWL Object into a built-in typed object
        saved_obj = save(cwl_obj)

        self.cwl_object = saved_obj

    def get_id(self):
        self.cwl_id = self.cwl_object["id"].replace('-', '_')

    def get_container_uri(self):
        container_uri = None
        if "requirements" in self.cwl_object:
            for req in self.cwl_object["requirements"]:
                if req["class"] == "DockerRequirement":
                    container_uri = req["dockerPull"]

        self.container_uri = container_uri

        # At this moment, latest alpine moment does not have /bin/bash,
        # meanwhile Nextflow hardcoded shell entrypoint as /bin/bash.
        # This means we can't run files which use 'alpine:latest'.
        # For now we change it to 'ubuntu:latest'.
        if self.container_uri and "alpine" in self.container_uri:
            self.container_uri = "ubuntu:latest" 

    def get_name(self):
        self.cwl_name = self.cwl_object["id"].replace('-', '_').split("#")[-1]

    def get_class(self):
        self.cwl_class = self.cwl_object["class"]

    def get_inputs(self):
        self.cwl_inputs = self.cwl_object["inputs"]

    def get_outputs(self):
        self.cwl_outputs = self.cwl_object["outputs"]

    def get_steps(self):
        if "steps" in self.cwl_object:
            self.cwl_steps = self.cwl_object["steps"]
        else:
            assert self.cwl_object["class"] != "Workflow"


    def get_dag(self):
        cwl_object_short_id = self.cwl_object['id'].split('#')[-1]
        # print(json.dumps(self.cwl_object['steps'], indent=4))
        for step in self.cwl_object['steps']:
            step_inp_dep = set()
            in_emit_dict = {}
            in_sink_dict = {}
            for in_elem in step['in']:
                assert(type(in_elem) is dict)
                if 'source' in in_elem:
                    prev_steps = []
                    if type(in_elem['source']) is list:
                        for e in in_elem['source']:
                            prev_step = e.split('/')[-2]
                            elem_name = e.split('/')[-1]
                            # If intermediate input comes from a step, 
                            # we take the step name. Else if it
                            # is an element of "inputs", we take 
                            # the element's name.
                            if prev_step != cwl_object_short_id:
                                prev_steps.append(prev_step)
                            else:
                                prev_steps.append(elem_name)
                            in_emit_dict[prev_step] = elem_name
                            in_sink_dict[in_elem['id']] = {
                                'step_source': prev_step,
                                'emit_source': elem_name
                            }
                    else:
                        prev_step = in_elem['source'].split('/')[-2]
                        elem_name = in_elem['source'].split('/')[-1]
                        if prev_step != cwl_object_short_id:
                            prev_steps.append(prev_step)
                        else:
                            prev_steps.append(elem_name)
                        in_emit_dict[prev_step] = elem_name
                        in_sink_dict[in_elem['id']] = {
                            'step_source': prev_step,
                            'emit_source': elem_name
                        }
                    # Add everything to dependency dict
                    for prev_step in prev_steps:
                        step_inp_dep.add(prev_step)

            self.step_dep[step['id']] = {
                'in': list(step_inp_dep),
                'in_emit_dict': in_emit_dict,
                'in_sink_dict': in_sink_dict,
                'run': step['run']
            }

            self.dag.add_node(step['id'])
            self.dag.add_nodes_from(step_inp_dep)
            self.dag.add_edges_from([(dep, step['id']) for dep in step_inp_dep])
            # Construct input source
            # input_sources = { dep_name: (lambda dep_name:  for input_name, input_source in }
            nx.set_node_attributes(self.dag, { step['id']: { 'in_emit_dict': self.step_dep[step['id']]['in_emit_dict'],
                                                                     'in_sink_dict':  self.step_dep[step['id']]['in_sink_dict']
                                                                     }}, 
                                   name='input_sources')

        # Also add outputs to NetworkX DAG
        output_nodes = []
        connections = []
        for out in self.cwl_object['outputs']:
            output_nodes.append(out['id'])
            prev_step = out['outputSource'].split('/')[-2]
            connections.append((prev_step, out['id']))
        
        self.dag.add_nodes_from(output_nodes)
        self.dag.add_edges_from(connections)

        # Here we show the order of executions
        all_nodes = {}
        for key, val in self.step_dep.items():
            if key not in all_nodes:
                all_nodes[key] = set()
            for node in val['in']:
                if node not in all_nodes:
                    all_nodes[node] = set()
                all_nodes[key].add(node)

    def show_dag(self):
        cwl_filename = self.cwl_object['id'].split('#')[0].split('/')[-1]
        # Also draw DAG to a PNG
        dag_drawer = DagDrawer(title=cwl_filename[:-len(".cwl")],
                               VERTICAL_SPACE=14,
                               HORIZONTAL_SPACE=20)
        dag_drawer.draw(self.dag)

    def get_dependencies(self):
        if self.cwl_steps:
            self.cwl_dependencies = [e["run"] for e in self.cwl_steps]

    def get_dependency_aliases(self):
        if self.cwl_steps:
            dep_aliases = {}
            for e in self.cwl_steps:
                if type(e["run"]) is str:
                    dep_aliases[e["run"]] = e["id"].replace('-', '_')

            self.cwl_dependency_aliases = dep_aliases

    def get_basecommand(self):
        if "baseCommand" in self.cwl_object:
            self.cwl_basecommand = self.cwl_object["baseCommand"]

    def get_arguments(self):
        if "arguments" in self.cwl_object:
            self.cwl_arguments = self.cwl_object["arguments"]
    
    def get_stdout(self):
        if "stdout" in self.cwl_object:
            self.cwl_stdout = self.cwl_object["stdout"]
    
    def print_cwl_object(self):
        print(self.cwl_object)

    def print_class(self):
        print(self.cwl_object["class"])

    def print_inputs(self):
        print(self.cwl_object["inputs"])

        for e in self.cwl_inputs:
            print(e)

    def print_id(self):
        print(self.cwl_object["id"].replace('-', '_'))

    def print_outputs(self):
        outputs = self.construct_outputs()
        print(outputs)

    def print_steps(self):
        for e in self.cwl_steps:
            print(e)

    def print_dependencies(self):
        for e in self.cwl_dependencies:
            print(e)

    def print_body(self):
        if self.cwl_class == "Workflow":
            alias_name = self.cwl_path.split('/')[-1][:-len(".cwl")]
            self.print_template_workflow(alias_name=str.upper(alias_name))
        # else:
        #     self.print_template_process()

    def print_template_config(self):
        file_dir = os.path.dirname(os.path.realpath(__file__))
        in_template = open(file_dir + "/templates/template_config.txt", "r")
        template = Template(in_template.read())

        subs = {
            "process_id": str.upper(self.cwl_name),
            "container_uri": self.container_uri,
        }

        results = template.substitute(subs)

        outfile = open("{}.nextflow.config".format(
            os.path.dirname(os.path.abspath(__file__)) + '/' + str.upper(self.cwl_name)),
            'w'
        )
        outfile.write(results)

    def print_template_process(self):
        file_dir = os.path.dirname(os.path.realpath(__file__))
        in_template = open(file_dir + "/templates/template_process.txt", "r")
        template = Template(in_template.read())

        pr_input_list = [e for e in self.cwl_inputs]

        pr_output_list = [e for e in self.cwl_outputs]

        opt_declaration = self.construct_process_opt_declaration()
        script = self.construct_script()

        pr_input = self.format_input_list(pr_input_list)
        pr_output = self.format_output_list(pr_output_list)
        pr_script = """
        {}

        \"\"\"
        {}
        \"\"\"
        """.format(
            opt_declaration, script
        )

        subs = {
            "process_id": str.upper(self.cwl_name),
            "container": "",  # By default we let container empty
            "input": pr_input,
            "output": pr_output,
            "script": pr_script,
        }
        if self.container_uri:
            subs["container"] = 'container "{}"\n'.format(self.container_uri)

        result = template.substitute(subs)

        print(result)

    def construct_process_opt_declaration(self):
        # We assume the value of all optional inputs at this point are "". The
        # input YAML might give 'null' as value for optional inputs, so this
        # function expect the null --> "" conversion was done elsewhere before.
        #
        # We handle optional inputs by checking whether their value is "". If it
        # isn't, we declare a new variable with concatenated prefix like shown
        # in this pattern:
        # https://nextflow-io.github.io/patterns/optional-input/ , e.g.: def
        # opt_adapters = opt_adapters != "" ? "--adapters ${opt_adapters}" : ""
        declaration_lines = []
        count_opt_input_file = 0 # Assume ordering of keys does not change
        for e in self.cwl_inputs:
            # assert('inputBinding' in e)
            prefix = ""
            separator = " "
            if "inputBinding" in e:
                if "prefix" in e["inputBinding"]:
                    prefix = e["inputBinding"]["prefix"]
                if "separator" in e["inputBinding"]:
                    is_separator = e['inputBinding']['separator']
                    separator = " " if is_separator else ""
            if self.check_is_optional(e):
                opt_name = "opt_" + e["id"].replace('-', '_')
                js_method_call = ".join(' ')" if self.check_is_array_of_items(e) else ""
                opt_declaration_val = "{}{}${{{}{}}}".format(prefix, separator, opt_name, js_method_call)

                # If input_type is path, we check its name against NO_FILE else
                # if it's val we simply check whether it's empty. Apparently
                # Nextflow complains if >1 optional input use file("NO_FILE"),
                # so we number them.
                val_to_check = None
                val_against = None
                if self.check_is_file(e) or self.check_is_dir(e):
                    count_opt_input_file += 1
                    val_to_check = opt_name + ".name"
                    val_against = '"NO_FILE{}"'.format(count_opt_input_file)
                else:
                    val_to_check = opt_name
                    val_against = '""'
                line = 'def {} = {} != {} ? "{}" : ""'.format(
                    opt_name, val_to_check, val_against, opt_declaration_val
                )
                declaration_lines.append(line)

        ret_declaration = "\n\t\t".join(declaration_lines)

        return ret_declaration

    def construct_workflow_all_params(self):
        all_yaml_path = open("/home/martinlp/nextflow/wf-translation-suite/wfe-log-parser/all.yaml", 'r')

        all_params_yaml = yaml.safe_load(all_yaml_path)
        print(all_params_yaml)

    def construct_workflow_opt_declaration(self):
        # For each optional input, we check whether it is a file or value. If it
        # is a file, we assign it with "file('NO_FILE<#>')", where <#> E {0,N; n
        # = number of optional input files}. This is due to Nextflow not
        # allowing us to use "file('NO_FILE')" for more than one variable. If it
        # is a val, we simply assign ""

        declaration_lines = []
        count_opt_input_file = 0
        for e in self.cwl_inputs:
            if self.check_is_optional(e):
                opt_name = "opt_" + e["id"].replace('-', '_')
                # If input_type is path, we check its name against NO_FILE else
                # if it's val we simply check whether it's empty. Apparently
                # Nextflow complains if >1 optional input use file("NO_FILE"),
                # so we number them.
                param_name = None
                default_val = None
                param_name = "params['{}']".format(e["id"].replace('-', '_'))
                if self.check_is_file(e) or self.check_is_dir(e):
                    count_opt_input_file += 1
                    default_val = 'file("NO_FILE{}")'.format(count_opt_input_file)
                else:
                    default_val = '""'
                line = '{} = {} ? {}{} : {}'.format(
                    opt_name, param_name, param_name,  (
                        ".path" if (self.check_is_file(e) or \
                                    self.check_is_dir(e)) 
                                else ""
                    ), 
                    default_val
                )
                declaration_lines.append(line)

        ret_declaration = "\n\t".join(declaration_lines)

        return ret_declaration

    def construct_workflow_arg(self, input_elem):
        arg_name = None
        if self.check_is_optional(input_elem):
            arg_name = "opt_" + input_elem['id']
        else:
            arg_name = "params['{}']".format(input_elem['id'])  + (
                            ".path" if (self.check_is_file(input_elem)  or \
                                        self.check_is_dir(input_elem))
                                    else ""
                        )
        return arg_name

    def handle_cadet_suffix(self, input_elem, suffix):
        ret = ""
        assert(suffix in input_elem['secondaryFiles'])
        if '^' in suffix:
            ret += ".take(params['{}'].path.lastIndexOf('.'))".format(input_elem['id'])

        return ret
        
    def construct_workflow_process_args(self):
        # If input is optional, we pass opt_<input_name> as argument, which is
        # defined by construct_workflow_opt_declaration. Else we pass
        # params['{}'].format(<input_name>). 
        arg_lines = []
        for e in self.cwl_inputs:
            arg = self.construct_workflow_arg(e)
            arg_lines.append(arg)
            # If input does not request secondary files,
            # we simply a it to the list of arguments.
            # Else we construct new arguments which takes 
            # the input's value plus the secondary files suffix.
            # Note that these new arguments might not be included 
            # in script construction, so we're adding the logic 
            # to include them directly here.
            if 'secondaryFiles' in e:
                for suffix in e['secondaryFiles']:
                    secondary_arg = self.construct_workflow_arg(e) + "{} + \"{}\"".format(
                            self.handle_cadet_suffix(e, suffix),
                            suffix.lstrip('^')
                        )
                    arg_lines.append(secondary_arg)
                
        ret_args = ',\n\t'.join(arg_lines) 
        return ret_args


    def construct_script(self):
        ret_script = " ".join(self.cwl_basecommand) \
                    if type(self.cwl_basecommand) is list \
                    else self.cwl_basecommand

        # First we take only inputs which has binding. So far we observe two
        # cases for inputs with binding: 1.) Single value, and 2.) Array For
        # array input, there might be multiple inputBindings, and they might be
        # hidden deep within the data structure. We might want to be careful on
        # handling that. Each array element might have an inputBinding which is
        # used to tell the prefix. However, the array itself can also have
        # inputBinding which tells which position should we put the array at. 
        # For example:
        #  filesB:
        #     type: 
        #       type: array
        #       items: string
        #       inputBinding:
        #           prefix: -B=
        #           separate: false
        #     inputBinding:
        #       position: 2
        # Tells us that filesB array should be in 2nd position, while each element
        # should be prefixed with -B= 
        inputs_with_binding = []
        for e in self.cwl_inputs:
            # Case of single value
            if 'inputBinding' in e:
                inputs_with_binding.append(e)
            else:
                # Case of array of items
                if 'type' in e:
                    if 'type' in e['type']:
                        if e['type']['type'] == 'array':
                            # We put a flag here for further handling
                            e['is_array'] = True
                            e['inputBinding'] = {}
                            if 'inputBinding' in e['type']:
                                inputs_with_binding.append(e)

        # Here we construct a new data structure which contained both elements
        # from arguments and also inputs. As argument elements do not have
        # inputBinding, our existing code might break, so we wrangle their type
        # a little bit append inputBinding to each element. Another thing to
        # note is we add a flag "is_argument: true" at the same level as
        # inputBinding. This will allow us to do different processing between
        # elements from arguments and inputs.
        #
        # We should still be able to differentiate elements form arguments and
        # inputs. For now we assume that arguments will never have 'id' as one
        # of their keys.
        new_arguments = []
        if self.cwl_arguments:
            new_arguments = [{"inputBinding":e} for e in self.cwl_arguments]
            for e in new_arguments:
                e['is_argument'] = True

        # If either input or new_arguments does not have position, we assign
        # priority based on the following rules: 1.) Arguments have higher
        # priority than inputs 2.) Among arguments and inputs themselves, we
        # assign the same arbitrary big value
        cmd_arguments = inputs_with_binding + new_arguments
        for e in cmd_arguments:
            if "inputBinding" in e:
                if "position" in e['inputBinding']:
                    pass
                else:
                    # Here we ensure that arguments have higher priority than
                    # inputs
                    if 'is_argument' in e:
                        assert(e['is_argument'] == True)
                        e['inputBinding']['position'] = 9998
                    else:
                        e['inputBinding']['position'] = 9999
        sorted_inputs = sorted(cmd_arguments, 
                            key=lambda x: x['inputBinding']['position'])
        arguments = [" "]
        for e in sorted_inputs:
            prefix_and_val = ""
            prefix = ""
            arg_val = ""
            separator = " "
            if 'is_argument' not in e:
                if "inputBinding" in e:
                    # Here we get the appropriate value for prefix, separator,
                    # and argument value.

                    # First we check the separator
                    if 'separate' in e['inputBinding']:
                        is_separator = e['inputBinding']['separate']
                        separator = " " if is_separator else ""

                    # We only include e in arguments if it is not optional, else
                    # we define a new variable and ensure "opt_<e>" resolves to
                    # ""
                    if "null" not in e["type"]:
                        # Here we handle an array
                        if 'is_array' in e:
                            prefix = e['type']['inputBinding']['prefix']
                        elif "prefix" in e["inputBinding"]:
                            prefix = e["inputBinding"]["prefix"] 
                        arg_val = "${{{}}}".format(e["id"].replace('-', '_'))
                    else:
                        arg_val = "${{{}}}".format("opt_" + e["id"].replace('-', '_'))
                    
                    # Now we construct
                    prefix_and_val = "{}{}{}".format(prefix, separator, arg_val)
                else:
                    # Input without inputBinding is usually processed with
                    # JavaScript.
                    pass
            else:
                separator = " "
                if "prefix" in e["inputBinding"]:
                    prefix = e["inputBinding"]["prefix"]
                if 'separate' in e['inputBinding']:
                    is_separator = e['inputBinding']['separate']
                    separator = " " if is_separator else "" 
                # Now most of the time the argument value here is in the form of
                # javascript. We need to find a way to handle this, but for now
                # let's just the whole thing as a string. This way it's easy to
                # fine-tune the translation output. 
                prefix_and_val = prefix + separator + self.translate_js(e['inputBinding']['valueFrom'])
            arguments.append(prefix_and_val)
            

        ret_script = " ".join(arguments) if len(arguments) > 1 else None
        # # If script is empty at this point, the CWL file must only contain an expression
        # if not ret_script:
        #     assert('expression' in self.cwl_object.keys())

        if self.cwl_stdout:
            ret_script += self.construct_stdout()


        # If inspect_docker is enabled, we check whether the container has entry
        # point. If it does, we add '#!' and '#' at beginning and end of the
        # script, respectively.
        if args.inspect_docker:
            image_url = self.container_uri
            bash_command = "docker inspect {}".format(image_url)
            process = subprocess.run(bash_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output = r"{}".format(process.stdout.decode('utf8'))
            inspect_result = json.loads(output)
            assert(len(inspect_result) == 1) # It doesn't make sense to have more than one
            assert("ContainerConfig" in inspect_result[0])
            entry_point = inspect_result[0]['ContainerConfig']['Entrypoint'] 
            if entry_point:
                ret_script = "#!" + ret_script
                ret_script = ret_script + " #"

        return ret_script

    def construct_stdout(self):
        redirect_outpath = self.translate_js(self.cwl_stdout)
        return " > " + redirect_outpath

    def handle_js_code_block(self, js_string):
        # For now we catalog the types of JS code blocks, then handle each of
        # them in a case-by-case basis.

        # The first pattern is checking whether a prefix is supplied, 
        # then append the prefix to various file names. 
        
        # For each append of optional prefix, we get the declared javascript
        # variable (group 1), CWL input variable in which we store the prefix
        # (group 2), and the appended str (group 4). 
        pattern = r"(?s)var[\s]+(\w+)[\s]+\=[\s]+inputs\.(\w+)[\s]+\?[\s]+.+(\1)[\s]+\+[\s]+[\'|\"](.+)[\'|\"]" 
        is_appending_opt_prefix = re.search(pattern, js_string)
        ret_str = ""
        if is_appending_opt_prefix:
            prefix_var = is_appending_opt_prefix.group(2)
            appended_str = is_appending_opt_prefix.group(4)

            # We need to check whether CWL input variable is optional.            
            is_var_opt = False
            for e in self.cwl_inputs:
                if (e['id'] == prefix_var):
                    is_var_opt = True
                    break
            opt_str = "opt_" if is_var_opt else ""
            ret_str = "${{{}{} ?: ''}}{}".format(opt_str, prefix_var, appended_str)
            
        return ret_str

    def translate_js(self, js_string):
        # The idea here is to first try to classify the code from 'eagle-eyed'
        # view, where we consider the whole code block, before looking into its
        # tokens. Ideally recognizing a whole code block should be done by an 
        # ML model in the future.
        ret_str = ""
        if self.check_is_block_of_js(js_string):
            ret_str += self.handle_js_code_block(js_string)
        else:            
            # Here we look into the code's fine-grained features (tokens).

            # We observe that we can broadly classify the javascript when binding
            #   as: $([string|call][ + (string|call)]*) e.g.: $("-Xmx" +
            # inputs.java_mem.toString() + "G") $(inputs.bam.nameroot + ".flagstat")
            # $('/usr/local/bin/picard.jar')
            #
            # An exception is the pattern: $(inputs.INPUT.nameroot).pdf, which
            # itself can be treated exactly like $(inputs.INPUT.nameroot + ".pdf")
            # (mentioned above). Thus, for now it suffice to handle the three cases
            # above.

            tokens = re.findall("([^(?!(\$\+))\s]+)", js_string)
            for token in tokens:
                if self.check_is_js_method_call(token):
                    ret_str += self.handle_js_method_call(token)
                elif self.check_is_js_get_attr(token):
                    requested_attr = token.lstrip(".inputs")
                    ret_str += "${{{}}}".format(requested_attr)
                else:
                    ret_str += token.replace('\"', '').replace('\'', '')

        return ret_str

    def handle_js_method_call(self, js_method_call):
        # Translate js function/method calls to its nextflow equivalent
        
        # Method calls should have 3 tokens when splitted by '.'
        if len(js_method_call.split('.')) != 3:
            return ""

        # We should get "inputs.<param_name>.<method_call>". We don't need the
        # 'inputs' part for now.
        _, param_name, method_call = js_method_call.split('.')

        # Here we simply maintain a 'lookup table' for js <--> nextflow calls
        ret_str = None
        if method_call == "nameroot":
            ret_str = "${{{}.baseName}}".format(param_name)
        elif str.lower(method_call) == "tostring":
            ret_str = "${{{}}}".format(param_name)
        elif method_call == "basename":
            ret_str = "${{{}.name}}".format(param_name)
        elif method_call == "dirname":
            # So far this method is used to move to a directory which contains
            # all other input files. 
            # This is useful when our WFE is cwltool, since the staged files 
            # might not be at the same place as initial working directory. 
            # But in the case of nextflow, staged files are always in the 
            # initial working directory, so we don't need to move around. 
            # Thus, we put '.' to indicate we will stay at the same directory. 
            ret_str = ".".format(param_name)
        else:
            # Ideally we should always let ret_str=None when we can't handle the
            # JS. This makes sure the program goes into error rather that
            # running with transient bugs. However this branch might be enabled
            # for debugging/development purpose. 
            ret_str = ""
        return ret_str

    def construct_outputs(self):
        ret_outputs = []
        # for e in self.cwl_outputs: output = { 'id': e['id'], 'glob':
        #     e['outputBinding']['glob'] if 'glob' in e['outputBinding'] else
        #         None, 'type':

        return ret_outputs

    def add_tabs(self, target_str: str, num_tabs: int = 1):
        new_str = ""
        for line in target_str.strip().split('\n'):
            tabs = ("\t" * num_tabs) 
            new_str += "{}{}\n".format(tabs, line)
        
        return new_str

    def find_short_filename_alias(self, short_filename):
        # short_filename: xyz.cwl (UNIX's basename, retaining extension)
        alias_name = None
        for rel_path, alias in self.cwl_dependency_aliases.items():
            rel_path_fn = rel_path.split('/')[-1]
            if rel_path_fn == short_filename:
                alias_name = alias
        
        assert(alias_name)
        return alias_name
    
    def get_short_filename(self):
        return self.cwl_id.split('/')[-1].split('#')[0]

    def print_template_workflow(self, alias_name: str = ""):

        in_template = open(
            "{}/templates/template_workflow.txt".format(os.getcwd()), "r"
        )
        template = Template(in_template.read())

        wf_take_list = [e["id"].replace('-', '_') for e in self.cwl_inputs]

        wf_emit_list = [e["id"].replace('-', '_') for e in self.cwl_outputs]

        wf_take = self.workflow_format_list(wf_take_list)

        # We order process calls based on Networkx's topological generations
        wf_main = ""
        for idx_generation, generation in enumerate(list(nx.topological_generations(self.dag))):
            for idx_node, node in enumerate(generation):
                if node in self.step_files: # not an input / output
                    step_file = self.step_files[node]
                    # Here we find the alias name
                    step_short_fn = step_file.get_short_filename()
                    step_alias = self.find_short_filename_alias(step_short_fn)

                    wf_main_entry = self.add_tabs(step_file.print_template_wf_body(step_alias), num_tabs=2)
                    wf_main += "{}".format(wf_main_entry)
    
        wf_emit = self.workflow_format_list(wf_emit_list)
        subs = {
            "workflow_id": alias_name if alias_name else str.upper(self.cwl_name),
            "take": wf_take,
            "main": wf_main,
            "emit": wf_emit,
        }

        result = template.substitute(subs)

        print(result)

    def print_template_standalone_text(self, alias: str = None):
        in_template = open(
            "{}/templates/template_standalone_test.txt".format(os.path.dirname(os.path.abspath(__file__))), "r"
        )
        template = Template(in_template.read())

        opt_declarations_str = self.construct_workflow_opt_declaration()
        process_name_str = str.upper(alias) if alias else str.upper(self.cwl_name)
        process_args = self.construct_workflow_process_args()

        subs = {
            "opt_declarations": opt_declarations_str,
            "process_name": process_name_str,
            "process_args": process_args
        }

        result = template.substitute(subs)

        print(result)

        return result


    def print_template_wf_body(self, alias_name: str = None):
        wf_body = self.construct_workflow_process_args()
        step_alias_name = str.upper(alias_name) if alias_name else str.upper(self.cwl_name)

        ret_str = "{}(\n\t{}\n)\n".format(step_alias_name, wf_body)

        return ret_str


    def print_template_header(self):
        if self.cwl_class == "Workflow":
            for e in self.cwl_dependencies:
                if not type(e) is dict:
                    lib_name = e.strip().split("/")[-1][:-len(".cwl")]
                    lib_alias = self.cwl_dependency_aliases[e]
                    lib_source = e.replace(".cwl", ".nf")
                    print(
                        "include {{ {} as {} }} from \"{}\"".format(
                            lib_name.upper().replace('.', '_'), lib_alias.upper(), lib_source
                        )
                    )
                else:
                    print(e)

    def construct_input_name(self, input_val):
        return input_val["id"].replace('-', '_') if "null" not in input_val["type"] \
                else ("opt_" + input_val["id"].replace('-', '_'))

    def format_input_list(self, inp_list: List) -> str:
        ret_str = ""
        for e in inp_list:
            e_name = self.construct_input_name(e)
            ret_str += self.get_inp_type(e) + " {}\n\t\t".format(e_name)
            # If input has secondaryFiles, we create new inputs as their placeholder
            if 'secondaryFiles' in e:
                for suffix in e['secondaryFiles']:
                    e_name = self.construct_input_name(e)
                    e_name += "_SECONDARY_{}".format(suffix.replace('^', '').replace('.', '').upper())
                    ret_str += self.get_inp_type(e) + " {}\n\t\t".format(e_name)

        return ret_str

    def check_is_array_of_items(self, inp):
        is_array_of_items = False
        for e in inp['type']:
            if type(e) is dict:
                assert('items' in e and 'type' in e)
                if (e['items'] == 'string' and e['type'] == 'array'):
                    is_array_of_items = True
                    break
        
        return is_array_of_items
                
    def check_is_optional(self, inp):
        is_optional = False
        if "null" in inp['type']:
            is_optional = True

        return is_optional

    def check_is_dir(self, inp):
        return inp['type'] == "Directory"

    def check_is_file(self, inp):
        is_file = False
        if type(inp["type"]) is str:
            is_file = inp["type"] == "File"
        elif type(inp["type"]) is list:
            is_file = "File" in inp["type"]
        elif type(inp["type"]) is dict:
            if "items" in inp["type"]:
                if inp["type"]["items"] == "File":
                    is_file = True

        return is_file

    def check_is_block_of_js(self, js_string):
        # How do we know whether a piece of string is a JS code block? By 'code
        # block' here we mean those code which have multiple lines, instead of a
        # single method call whose examples abound in GDC workflows. There
        # should be a better way to classify such code block, but for now we opt
        # for a very simplistic 'does this include variable declaration?'

        return "var " in js_string # Variable declaration like 'var pfx = ..'

    def check_is_js_method_call(self, js_string):
        # We consider a string as a method call if it either:
        # 1.) Has "self.", or
        # 2.) It consists of 3 tokens separated by '.'
        is_method_call = re.match("([^?!\'w]+)\.([^w]+)\.([^w]+)", js_string)
        if ("self." in js_string) or \
            is_method_call:
            return True
        else:
            return False

    def check_is_js_get_attr(self, js_string):
        # We check whether it has "inputs."
        if "inputs." in js_string:
            return True
        else:
            return False

    def get_inp_type(self, inp):
        inp_type = None
        if self.cwl_class == "Workflow":
            inp_type = ""
        else:  # CommandLineTool or ExpressionTool
            if self.check_is_file(inp) or self.check_is_dir(inp):
                inp_type = "path"
            else:
                inp_type = "val"

        return inp_type

    def get_out_type(self, out):
        out_type = None
        if self.cwl_class == "Workflow":
            out_type = ""
        else:  # CommandLineTool or ExpressionTool
            if "outputBinding" in out:
                if out["type"] == "File" or "*" in out["outputBinding"]["glob"]:
                    out_type = "path"
                else:
                    out_type = "val"
            else: 
                out_type = "val"

        return out_type

    def format_output_list(self, out_list: List) -> str:
        ret_str = ""
        for e in out_list:
            if "outputBinding" in e:
                cwl_glob =  e["outputBinding"]["glob"]
                out_glob = None
                if "$" in cwl_glob:
                    # This is most likely a javascript
                    out_glob = self.translate_js(cwl_glob)
                else: 
                    out_glob = cwl_glob
                ret_str += self.get_out_type(e) + ' "{}"\t, emit: {}\n\t\t'.format(
                out_glob, e["id"].replace('-', '_')
                )
                if 'secondaryFiles' in e:
                    for suffix in e['secondaryFiles']:
                        ret_str += '\t' + self.get_out_type(e) + ' "{}{}"\t, emit: {}_{}\n\t'.format(
                            out_glob, suffix, e["id"].replace('-', '_'), suffix.lstrip('.')
                        )
            else:
                ret_str += self.get_out_type(e) + ' {}\t, emit: {}\n\t\t'.format(
                    e["id"].replace('-', '_'), e["id"].replace('-', '_')
                )

        return ret_str

    def workflow_format_list(self, inp_list: List) -> str:
        ret_str = ""
        for inp in inp_list:
            ret_str += "\n\t\t" + inp

        return ret_str

    def print_script(self):
        # self.print_header()
        self.print_body()

    def collect_step_file_recursive(self):
        all_step_file = {}
        for step in self.cwl_steps:
            if 'run' in step:
                if type(step['run']) is str:
                    new_file_path = '/'.join(self.cwl_path.split('/')[:-1]) + '/' + step['run']
                    # all_step_file.append(new_file_path)
                    new_cwl = CWLStep(new_file_path, self)
                    if new_cwl.cwl_class == "Workflow":
                        print()
                        print("Adding workflow {}".format(new_cwl.cwl_id))
                        all_step_file[step['id']] = new_cwl
                        all_step_file.update(new_cwl.collect_step_file_recursive())
                        print("Done workflow {}".format(new_cwl.cwl_id))
                        print()
                        # new_cwl.collect_step_file_recursive()
                    else:
                        all_step_file[step['id']] = new_cwl
                        print("Added {}".format(new_cwl.cwl_id))
                        alias_name = self.cwl_dependency_aliases[step['run']]
                        # new_cwl.print_template_standalone_text(alias=alias_name)
                        # new_cwl.print_body()
                        # new_cwl.print_template_config()
                else:
                    return {}

        self.step_files = all_step_file

        return all_step_file

    def import_yaml(self, tool_name: str):
        yaml_path = "/home/martinlp/nextflow/wf-translation-suite/wfe-log-parser/{}.yaml".format(tool_name)

        yaml_file = open(yaml_path, 'r')
        self.param_dict[tool_name] = (yaml.safe_load(yaml_file))

    def construct_params_from_yaml(self, tool_names):
        tool_names = [
            "gdc_qc_tool.fastqc",
            "gdc_qc_tool.samtools_flagstat",
            "bam_readgroup_to_json",
            "gdc_qc_tool.samtools_idxstat",
            "gdc_qc_tool.readgroups",
            "gdc_qc_tool.picard",
            "picard_collectrnaseqmetrics",
            "samtools_flagstat",
            "samtools_idxstat",
            "samtools_stats",
            "gdc_qc_tool.samtools_stats",
            "gdc_qc_tool.star"
        ]
        for tool_name in tool_names:
            self.import_yaml(tool_name)

        
def convert_to_nextflow(cwl_path):
    # File Input - This is the only thing you will need to adjust or take in as
    # an input to your function: cwl_file =
    # Path("./gdc_dnaseq.bamfastq_align.workflow.cwl")  # or a plain string
    # works as well
    cwl_file = Path(cwl_path)

    # Import CWL Object
    cwl_obj = load_document_by_uri(cwl_file)

    # Export CWL Object into a built-in typed object
    saved_obj = save(cwl_obj)

    print(str_obj)


class CWLStep(CWLFile):
    def __init__(self, cwl_path, parent: CWLFile) -> None:
        super().__init__(cwl_path)
        self.parent = parent

    def get_source_sink_info(self, input_elem):
        # Here we get the following info for ourselves:
        # 1.) Parent's output as source
        # 2.) Our input as sink

        # First find ourselves among parent's steps
        short_fn = self.get_short_filename()
        alias = self.parent.find_short_filename_alias(short_fn)
        source_sink_dict = nx.get_node_attributes(self.parent.dag, "input_sources")[alias]['in_sink_dict']

        # It is possible that the dict does not contain input_elem 
        # (e.g. input_elem has a default value)
        step_source = None
        emit_source = None
        if input_elem['id'] in source_sink_dict:
            step_source = source_sink_dict[input_elem['id']]['step_source']
            emit_source = source_sink_dict[input_elem['id']]['emit_source']
        # Finally, we check whether it comes from 
        # a step within the workflow, or directly from workflow's "take" 
        if step_source == self.parent.cwl_id.split('#')[-1]:
            if emit_source:
                # Comes directly from workflow's "take"
                return "{}".format(emit_source)
            else:
                return None
        else:
            if (emit_source and step_source):
                return "{}.out.{}".format(step_source, emit_source)
            else:
                assert((not emit_source) and (not step_source))
                return None


    def construct_workflow_arg(self, input_elem):
        arg_name = None
        # Get info from source-sink info from parent's DAG
        elem_parent_perspective = self.get_source_sink_info(input_elem)
        # If doesn't exist from the parent's perspective, 
        # it must be one of (or a combination of) optional, has default value, or 
        # comes from the param file.
        if elem_parent_perspective:
            # print("ELEM PARENT", elem_parent_perspective)
            arg_name = elem_parent_perspective
        else:
            if self.check_is_optional(input_elem):
                arg_name = "opt_" + input_elem['id']
            else:
                arg_name = "params['{}']".format(input_elem['id'])  + (
                                ".path" if (self.check_is_file(input_elem)  or \
                                            self.check_is_dir(input_elem))
                                        else ""
                            )   
        return arg_name


if __name__ == "__main__":
    # print(sys.argv[1])
    cwl_path = sys.argv[1]

    cwl_file = CWLFile(cwl_path)
    if cwl_file.cwl_class == "Workflow":
        # cwl_file.print_steps()
        # cwl_file.print_dependencies()
        cwl_file.print_template_header()
        # cwl_file.print_body()
        all_step_file = cwl_file.collect_step_file_recursive()
        for step_id, step_file in all_step_file.items():
            # alias_name =step.cwl_dependency_aliases[step['run']]
            step_file.print_template_header()
            step_file.print_body()               
    else:
        cwl_file.print_body()

    # Here we add options to: 1.) Print out corresponding nextflow.config for
    # docker usage print(cwl_file.cwl_arguments)
    if args.with_tool_config:
        cwl_file.print_template_config()
    # 2.) Concatenate workflow() for testing purpose
    if args.unit_test:
        cwl_file.print_template_standalone_text()
    
    if args.with_wf_config:
        cwl_file.collect_step_file_recursive()

    if args.show_dag:
        cwl_file.show_dag()
        # print(json.dumps(cwl_file.cwl_object, indent=4))
            
    if args.import_yaml:
        # cwl_file.construct_params_from_yaml()
        cwl_file.construct_workflow_all_params()

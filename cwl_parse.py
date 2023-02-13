from argparse import ArgumentParser
from cwl_utils.parser import load_document_by_uri, save
from pathlib import Path
from string import Template
from typing import List

import logging
import os
import re
import sys


# Logging
logger = logging.getLogger(__name__)
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.DEBUG)

# Argument parsing
parser = ArgumentParser(description="Automated CWL <-> X Translator")
parser.add_argument("cwl_source", help="Source CWL file")
parser.add_argument(
    "--unit_test",
    action="store_true",
    help="Append a Workflow(...) call for stand-alone execution.",
)
parser.add_argument(
    "--with_config",
    action="store_true",
    help="Generate <cwl_source>.nextflow.config for docker usage.",
)
parser.set_defaults(unit_test=False)
parser.set_defaults(with_config=False)
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
        self.cwl_id = self.cwl_object["id"]

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
        if "alpine" in self.container_uri:
            self.container_uri = "ubuntu:latest" 

    def get_name(self):
        self.cwl_name = self.cwl_object["id"].split("#")[-1]

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

    def get_dependencies(self):
        if self.cwl_steps:
            self.cwl_dependencies = [e["run"] for e in self.cwl_steps]

    def get_dependency_aliases(self):
        if self.cwl_steps:
            dep_aliases = {}
            for e in self.cwl_steps:
                dep_aliases[e["run"]] = e["id"]

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
        print(self.cwl_object["id"])

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
                opt_name = "opt_" + e["id"]
                opt_declaration_val = "{}{}${{{}}}".format(prefix, separator, opt_name)

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
                opt_name = "opt_" + e["id"]
                # If input_type is path, we check its name against NO_FILE else
                # if it's val we simply check whether it's empty. Apparently
                # Nextflow complains if >1 optional input use file("NO_FILE"),
                # so we number them.
                param_name = None
                default_val = None
                param_name = "params." + e["id"]
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
            arg_name = "params." + input_elem['id']  + (
                            ".path" if (self.check_is_file(input_elem)  or \
                                        self.check_is_dir(input_elem))
                                    else ""
                        )
        return arg_name

    def construct_workflow_process_args(self):
        # If input is optional, we pass opt_<input_name> as argument, which is
        # defined by construct_workflow_opt_declaration. Else we pass
        # params.<input_name>. 
        arg_lines = []
        for e in self.cwl_inputs:
            arg = self.construct_workflow_arg(e)
            arg_lines.append(arg)
            # If input does not request secondary files,
            # we simply add it to the list of arguments.
            # Else we construct new arguments which takes 
            # the input's value plus the secondary files suffix.
            # Note that these new arguments might not be included 
            # in script construction, so we're adding the logic 
            # to include them directly here.
            if 'secondaryFiles' in e:
                for suffix in e['secondaryFiles']:
                    secondary_arg = self.construct_workflow_arg(e) + " + \"{}\"".format(suffix)
                    arg_lines.append(secondary_arg)
                
        ret_args = ',\n\t\t'.join(arg_lines) 
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
                        arg_val = "${{{}}}".format(e["id"])
                    else:
                        arg_val = "${{{}}}".format("opt_" + e["id"])
                    
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
                
        ret_script += " ".join(arguments)
        if self.cwl_stdout:
            ret_script += self.construct_stdout()

        return ret_script

    def construct_stdout(self):
        redirect_outpath = self.translate_js(self.cwl_stdout)
        return " > " + redirect_outpath

    def handle_js_code_block(self, js_string):
        # For now we catalog the types of JS code blocks, then handle each of
        # them in a case-by-case basis.

        # The first pattern is checking whether a prefix is supplied, 
        # then append the prefix to various file names. 
        
        pattern = r"(?s)var[\s]+(\w+)[\s]+\=[\s]+inputs\.(\w+)[\s]+\?[\s]+.+(\1)[\s]+\+[\s]+[\'|\"](.+)[\'|\"]" 
        is_appending_opt_prefix = re.search(pattern, js_string)
        ret_str = ""
        if is_appending_opt_prefix:
            prefix_var = is_appending_opt_prefix.group(2)
            appended_str = is_appending_opt_prefix.group(4)
            
            is_var_opt = False
            for e in self.cwl_inputs:
                if (e['id'] == prefix_var):
                    is_var_opt = True
                    break
            opt_str = "opt_" if is_var_opt else ""
            ret_str = "${{{}{} ? : ''}}{}".format(opt_str, prefix_var, appended_str)
            
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
        assert(len(js_method_call.split('.')) == 3)

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

    def print_template_workflow(self):

        in_template = open(
            "{}/templates/template_workflow.txt".format(os.getcwd()), "r"
        )
        template = Template(in_template.read())

        wf_take_list = [e["id"] for e in self.cwl_inputs]
        wf_take_bindings = []

        wf_main_list = ["{}()".format(e["id"].upper()) for e in self.cwl_steps]
        wf_main_bindings = []

        wf_emit_list = [e["id"] for e in self.cwl_outputs]
        wf_emit_bindings = []

        wf_take = self.format_list(wf_take_list, wf_take_bindings)
        wf_main = self.format_list(wf_main_list, wf_main_bindings)
        wf_emit = self.format_list(wf_emit_list, wf_emit_bindings)

        subs = {
            "workflow_id": str.upper(self.cwl_name),
            "take": wf_take,
            "main": wf_main,
            "emit": wf_emit,
        }

        result = template.substitute(subs)

        print(result)

    def print_template_standalone_text(self):
        in_template = open(
            "{}/templates/template_standalone_test.txt".format(os.path.dirname(os.path.abspath(__file__))), "r"
        )
        template = Template(in_template.read())

        opt_declarations_str = self.construct_workflow_opt_declaration()
        process_name_str = str.upper(self.cwl_name)
        process_args = self.construct_workflow_process_args()

        subs = {
            "opt_declarations": opt_declarations_str,
            "process_name": process_name_str,
            "process_args": process_args
        }

        result = template.substitute(subs)

        print(result)

    def print_template_header(self):
        for e in self.cwl_dependencies:
            lib_name = e.strip().split("/")[-1].strip(".cwl")
            lib_alias = self.cwl_dependency_aliases[e]
            lib_source = e.replace(".cwl", ".nf")
            print(
                "INCLUDE {{ {} as {} }} from {}".format(
                    lib_name.upper(), lib_alias.upper(), lib_source
                )
            )

    def construct_input_name(self, input_val):
        return input_val["id"] if "null" not in input_val["type"] \
                else ("opt_" + input_val["id"])

    def format_input_list(self, inp_list: List) -> str:
        ret_str = ""
        for e in inp_list:
            e_name = self.construct_input_name(e)
            ret_str += self.get_inp_type(e) + " {}\n\t\t".format(e_name)
            # If input has secondaryFiles, we create new inputs as their placeholder
            if 'secondaryFiles' in e:
                for suffix in e['secondaryFiles']:
                    e_name = self.construct_input_name(e)
                    e_name += "_SECONDARY_{}".format(suffix.lstrip('.').upper())
                    ret_str += self.get_inp_type(e) + " {}\n\t\t".format(e_name)

        return ret_str

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
            if out["type"] == "File" or "*" in out["outputBinding"]["glob"]:
                out_type = "path"
            else:
                out_type = "val"

        return out_type

    def format_output_list(self, out_list: List) -> str:
        ret_str = ""
        for e in out_list:
            cwl_glob =  e["outputBinding"]["glob"]
            out_glob = None
            if "$" in cwl_glob:
                # This is most likely a javascript
                out_glob = self.translate_js(cwl_glob)
            else: 
                out_glob = cwl_glob
            ret_str += self.get_out_type(e) + ' "{}"\t, emit: {}\n\t'.format(
               out_glob, e["id"]
            )
            if 'secondaryFiles' in e:
                for suffix in e['secondaryFiles']:
                    ret_str += '\t' + self.get_out_type(e) + ' "{}{}"\t, emit: {}_{}\n\t'.format(
                        out_glob, suffix, e["id"], suffix.lstrip('.')
                    )

        return ret_str

    def format_list(self, inp_list: List, bindings: List) -> str:
        inp_type = "" if self.cwl_class == "Workflow" else "val "

        ret_str = (
            inp_type + inp_list[0] + "({})".format(bindings[0]) if len(bindings) else ""
        )

        if len(inp_list) > 1:
            assert len(inp_list) == len(bindings)
            for e, b in zip(inp_list[1:], bindings[1:]):
                ret_str += "\n\t\t" + inp_type + e + "({})".format(b)

        return ret_str

    def print_script(self):
        # self.print_header()
        self.print_body()


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

    # Here we add options to: 1.) Print out corresponding nextflow.config for
    # docker usage print(cwl_file.cwl_arguments)
    if args.with_config:
        cwl_file.print_template_config()
    # 2.) Concatenate workflow() for testing purpose
    if args.unit_test:
        cwl_file.print_template_standalone_text()
from pathlib import Path
# from ruamel import yaml
import json
import sys
import re
# import json5

from cwl_utils.parser import load_document_by_uri, save
import ast

class CWLParser:
    def __init__(self, input: str) -> None:
        self.remaining_input = input # Strip first {}
        self.current_token = None
        self.current_chars = None # Workaround for ${ and $(
        self.current_open_list = [] # Stack of 'active' quotes or brackets
        self.outfile = open("results.txt", 'w')
    
    def advance(self):
        self.remaining_input = self.remaining_input[len(self.current_chars):]

    def get_chars(self):
        if self.remaining_input[0] == '$':
            self.current_chars = self.remaining_input[:2] # Take ${ or $(
        else:
            self.current_chars = self.remaining_input[0]
        return self.current_chars
    
    def emit_token(self):
        self.outfile.write(self.current_token + '\n')

    def current_open_pair(self):
        assert(self.current_open_list)
        current_open = self.current_open_list[-1]
        if current_open in "'\"":
            return current_open # Same
        else:
            if current_open == '{':
                return '}'
            elif current_open == '[':
                return ']'
            elif current_open == "${":
                return '}'
            elif current_open == "$(":
                return ')'
            else:
                assert(False)

    def get_current_token(self):
        self.current_token = None
        quotes_and_brackets = ["\"", "'", "[", "{", "}", "]", "${", "$(", ")"]
        if self.get_chars() in quotes_and_brackets:
            print("Got quotes and brackets: {}".format(self.get_chars()))
            if not self.current_open_list:
                self.current_token = self.get_chars()
                self.current_open_list.append(self.get_chars())
            else:
                if self.get_chars() == self.current_open_pair():
                    self.current_token = self.get_chars()
                    self.current_open_list.pop()
                else:
                    self.current_token = self.get_chars()
                    self.current_open_list.append(self.get_chars())
            self.advance()
        else:
            if self.get_chars() in ":,":
                self.advance()
                self.advance() 
            else:
                token = ""
                while not self.get_chars() in quotes_and_brackets:
                    token += self.get_chars()
                    self.advance()
                self.current_token = token
        print("Key open: {}, Current token: {}".format(self.current_open_list, self.current_token))

    def parse(self):
        self.get_current_token()
        if self.current_token:
            self.emit_token()
        
        while self.current_open_list:
            # self.get_key()
            # self.get_val()
            self.get_current_token()
            if self.current_token:
                self.emit_token()
            # self.advance()
            print("Remaining input: {}".format(self.remaining_input))

# File Input - This is the only thing you will need to adjust or take in as an input to your function:
# cwl_file = Path("./gdc_dnaseq.bamfastq_align.workflow.cwl")  # or a plain string works as well
cwl_file = Path("../tools/bwa_record_pe.cwl")

# Import CWL Object
cwl_obj = load_document_by_uri(cwl_file)

# View CWL Object
# print("List of object attributes:\n{}".format("\n".join(map(str, dir(cwl_obj)))))

# Export CWL Object into a built-in typed object
saved_obj = save(cwl_obj)

print(str(saved_obj))

parser = CWLParser(str(saved_obj))
parser.parse()

# Parsing
# str_obj = str(saved_obj)
# with open("outfile.txt", 'w') as outfile:
#     outfile.write(str_obj)


# re_keys = re.compile(r"\'([A-z]*[a-z]*)\':")
# re_vals = re.compile(r"\': \'(.+?)\'")
# # keys = re_keys.findall(r'"\1":', str_obj)
# # vals = re_vals.findall("': \"\\1\"", str_obj)
# keys = re_keys.findall(str_obj)
# with open('keys.txt', 'w') as outfile:
#     for e in keys:
#         outfile.write(e + '\n')

# vals = re_vals.findall(str_obj)
# with open('vals.txt', 'w') as outfile:
#     for e in vals:
#         outfile.write(e + '\n')

# for key, val in zip(keys, vals):
#     print(key, val)

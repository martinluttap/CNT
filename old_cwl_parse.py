import json
import os
import re


class CWLParser:
    def __init__(self, input: str) -> None:
        self.remaining_input = input  # Strip first {}
        self.current_token = None
        self.current_chars = None  # Workaround for ${ and $(
        self.current_open_list = []  # Stack of 'active' quotes or brackets
        self.outfile = open("results.txt", "w")
        self.quotes_and_brackets = ['"', "'", "[", "{", "}", "]", "${", "$(", ")"]

    def advance(self):
        self.remaining_input = self.remaining_input[len(self.current_chars) :]

    def get_chars(self):
        if self.remaining_input[0] == "$":
            self.current_chars = self.remaining_input[:2]  # Take ${ or $(
        else:
            self.current_chars = self.remaining_input[0]
        return self.current_chars

    def emit_token(self):
        self.outfile.write(self.current_token + "\n")

    def current_open_pair(self):
        # assert(self.current_open_list)
        if self.current_open_list:
            current_open = self.current_open_list[-1]
            if current_open in "'\"":
                return current_open  # Same
            else:
                if current_open == "{":
                    return "}"
                elif current_open == "[":
                    return "]"
                elif current_open == "${":
                    return "}"
                elif current_open == "$(":
                    return ")"
                else:
                    assert False
        else:
            return None

    def get_current_token(self):
        self.current_token = None
        if self.get_chars() in self.quotes_and_brackets:
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
                while not self.get_chars() in self.quotes_and_brackets:
                    token += self.get_chars()
                    self.advance()
                self.current_token = token
        print(
            "Key open: {}, Current token: {}".format(
                self.current_open_list, self.current_token
            )
        )

    def get_key(self):
        print("Get Key, curr char ", self.get_chars())
        if self.get_chars() == ",":
            self.advance()
            self.advance()
        if self.get_chars() == "{":
            self.advance()
        if self.get_chars() == "}":
            self.advance()
            self.advance()
            self.advance()
            self.advance()
        if self.get_chars() == "]":
            self.advance()
            self.advance()
            self.advance()
        if self.get_chars() == "[":
            self.advance()
            self.advance()

        self.current_open_list.append(self.get_chars())
        self.advance()
        key = ""
        while not self.get_chars() == self.current_open_pair():
            key += self.get_chars()
            self.advance()
        self.current_open_list.pop()
        self.advance()

        return key

    def get_val(self):
        print("Get val, curr char ", self.get_chars())
        self.advance()  # : to <whitespace>
        self.advance()  # <whitespace> to ??, checked below
        val = ""

        if self.get_chars() == "'":
            self.current_open_list.append(self.get_chars())
            self.advance()
            while not self.get_chars() == self.current_open_pair():
                val += self.get_chars()
                self.advance()
            self.current_open_list.pop()
            self.advance()

        elif self.get_chars() == "[":
            val = self.process_array()

        elif self.get_chars() == "T":
            val = self.process_true()

        elif self.get_chars() == "F":
            val = self.process_false()

        return val

    def process_true(self):
        ret = ""
        while not ret == "True":
            ret += self.get_chars()
            self.advance()

        return ret

    def process_false(self):
        ret = ""
        while not ret == "False":
            ret += self.get_chars()
            self.advance()

        return ret

    def process_array(self):
        self.advance()
        arr = []
        if self.get_chars() == "{":
            while not self.get_chars() == "]":
                key = self.get_key()
                val = self.get_val()
                arr.append((key, val))
        else:
            arr_val = ""
            while not self.get_chars() == "]":
                arr_val += self.get_chars()
                self.advance()
            arr.append(arr_val)
        return arr
        # val = ""
        # self.advance()
        # if self.get_chars() == '{':
        #     val += self.get_key()
        #     val += self.get_val()
        # else:
        #     while self.get_chars() != ']':
        #         val += self.get_chars()
        #         self.advance()
        #     self.advance()
        # return val

    def parse(self):
        # self.get_current_token()
        # if self.current_token:
        #     self.emit_token()
        self.remaining_input = self.remaining_input[1:-1]
        while self.remaining_input:
            key = self.get_key()
            val = self.get_val()
            print(key, val)
            # self.get_current_token()
            # if self.current_token:
            #     self.emit_token()
            # self.advance()
        # print("Remaining input: {}".format(self.remaining_input))


def convert_to_json(cwl_path):
    # File Input - This is the only thing you will need to adjust or take in as an input to your function:
    # cwl_file = Path("./gdc_dnaseq.bamfastq_align.workflow.cwl")  # or a plain string works as well
    cwl_file = Path(cwl_path)

    # Import CWL Object
    cwl_obj = load_document_by_uri(cwl_file)

    # Export CWL Object into a built-in typed object
    saved_obj = save(cwl_obj)

    # Parsing
    str_obj = str(saved_obj)

    # Handling javascript

    # First we mark the javascript part,
    # then swap single-quote to double-quotes and vice versa
    # print(str_obj)
    re_mark_js = re.compile(r"(\"|\')\${(.+?)(\'|\")}")
    groups = re_mark_js.search(str_obj)
    mark_js_str = re_mark_js.sub(lambda x: x.group(0).replace("'", "@'"), str_obj)
    mark_js_str = re_mark_js.sub(lambda x: x.group(0).replace('"', "'"), mark_js_str)
    mark_js_str = re_mark_js.sub(lambda x: x.group(0).replace("@'", '"'), mark_js_str)
  
    # Also swap single quote to double quotes and vice versa
    mark_js_str = mark_js_str.replace("'", "@'")
    mark_js_str = mark_js_str.replace('"', "'")
    mark_js_str = mark_js_str.replace("@'", '"')
    # Suppress weird bug
    mark_js_str = mark_js_str.replace("@\"", "\"")
  
    # If javascript is in the middle of json, its double quotes might not be enclosed
    # re_unenclosed_dq = re.compile(r"(\"\$.+)},")
    # mark_js_str = mark_js_str.replace("'},", "}\",")
    # Handling the rest

    # Change boolean into lowercases
    mark_js_str = mark_js_str.replace(": False", ": false")
    mark_js_str = mark_js_str.replace(": True", ": true")
    # Handling expressionLib
    re_expression_lib = re.compile(r"\"expressionLib\": (\[\".+?\"\])")
    re_file_size_multiplier = re.compile(r", \'.+?\}.+\'")
    if re_expression_lib.search(mark_js_str):
        old_file_size_multiplier = re_file_size_multiplier.findall(
            re_expression_lib.findall(mark_js_str)[0]
        )[0]
        new_file_size_multiplier = (
            old_file_size_multiplier.replace("'", "@'").replace('"', "'").replace("@'", '"')
        )
        mark_js_str = mark_js_str.replace(old_file_size_multiplier, new_file_size_multiplier)
    # Handling valueFrom javascripts
    re_value_from = re.compile(r"\'\$\(.+?\)\'")
    # re_value_from = re.compile(r"\"valueFrom\": \'\$\(.+?\)\'")
    if re_value_from.search(mark_js_str):
        old_value_froms = re_value_from.findall(mark_js_str)
        new_value_froms = []
        for old_vf in old_value_froms:
            old_vf = old_vf.replace("'", "@'")
            old_vf = old_vf.replace('"', "'")
            new_vf = old_vf.replace("@'", '"')
            new_value_froms.append(new_vf)

        for old_vf, new_vf in zip(old_value_froms, new_value_froms):
            mark_js_str = mark_js_str.replace(old_vf, new_vf)

    # print(mark_js_str)
    # print(mark_js_str[500:530])
    final_json = json.loads(mark_js_str)

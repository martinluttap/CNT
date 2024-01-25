from argparse import ArgumentParser

# Argument parsing
parser = ArgumentParser(description="Automated CWL <-> X Translator")
parser.add_argument("cwl_source", help="Source CWL file")

parser.add_argument(
    "--import-yaml",
    action="store_true",
    help="Read YAML import_yaml",
)
parser.add_argument(
    "--inspect-docker",
    action="store_true",
    help="Optimize translation by allowing inspection of docker images. \
        WARNING: This option will pull the docker images and does not clean them up.",
)
parser.add_argument(
    "--unit-test",
    action="store_true",
    help="Append a Workflow(...) call for stand-alone execution.",
)
parser.add_argument(
    "--with-tool-config",
    action="store_true",
    help="Generate <cwl_source>.nextflow.config for docker usage.",
)
parser.add_argument(
    "--show-dag",
    action="store_true",
    help="Print the DAG of workflow.",
)
parser.add_argument(
    "--single-file",
    action="store_true",
    help="Output translation only for the given file",
)
parser.add_argument(
    "--with-wf-config",
    action="store_true",
    help="Generate <cwl_source>.nextflow.config with all tool configs associated with this cwl.",
)
parser.add_argument(
    "--export-to-dir",
    nargs=1,
    metavar=("TARGET_DIR"),
    action="store",
    help="Export per-file translation to designated TARGET_DIR following SOURCE_DIR's structure.",
)

parser.set_defaults(inspect_docker=False)
parser.set_defaults(unit_test=False)
parser.set_defaults(with_tool_config=False)
args = parser.parse_args()

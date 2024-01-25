from cwl.objects import CWLFile
from cli import args
from logger import logger

import csv 
import json
import sys

if __name__ == "__main__":
    cwl_path = sys.argv[1]

    cwl_file = CWLFile(cwl_path=cwl_path,
                       logger=logger,
                       inspect_docker=args.inspect_docker)

    if cwl_file.cwl_class == "Workflow":
        all_step_file = cwl_file.collect_step_file_recursive()
        for step_id, step_file in all_step_file.items():
            step_file.print_template_header()
            step_file.print_params_declaration()
            step_file.print_body()

        # We need to the dependencies after adding all step files
        cwl_file.update_self()

        cwl_file.print_template_header()
        cwl_file.print_params_declaration()
        cwl_file.print_body()
    else:
        cwl_file.print_body()
    
    if args.single_file:
        if cwl_file.cwl_class == "Workflow":
            print(cwl_file.print_template_header())
            print(cwl_file.print_params_declaration())
            print(cwl_file.print_body())
        else:
            print(cwl_file.print_body())


    if args.with_tool_config:
        cwl_file.print_template_config()

    if args.unit_test:
        cwl_file.print_template_standalone_text()
    
    if args.with_wf_config:
        # We handle two things when generaitng config
        # First, we manage docker path and entrypoint 
        # Second, we manage both static and dynamic resource requirements. 
        all_step_file = cwl_file.collect_step_file_recursive()
        
        req_d = {}
        for step_id, step_file in all_step_file.items():
            if (step_file.cwl_class != "Workflow") and step_file.container_uri:
                process_name = str.upper(step_file.get_short_filename()[:-len('.cwl')])
                logger.debug("withName:{}{{\n\tcontainer = '{}'\n}}".format(
                    process_name, 
                              step_file.container_uri)
                )
                req_d[process_name] = step_file.resource_requirements
        # Here we handle the extraction of both static and dynamic resource requirements. 
        # (a) We first handle static requirements simply by extracting the ResourceRequirement field in the CWL
        logger.debug(json.dumps(req_d, indent=4))
        
        with open('requirements.csv', 'w', newline='') as csvfile:
            fieldnames = ['withName', 'coresMin', 'ramMin', 'tmpdirMin', 'outdirMin']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            for proc_name, req in req_d.items():
                try:   
                    row = {'withName': proc_name, 
                        'coresMin': req['coresMin'],
                        'ramMin': req['ramMin'],
                        'tmpdirMin': req['tmpdirMin'],
                        'outdirMin': req['outdirMin']
                        }
                except KeyError as e:
                    logger.debug(f"Key {e} does not exist for {proc_name}. Available keys: {req.keys()}")
                writer.writerow(row)

    if args.show_dag:
        cwl_file.show_dag()
            
    if args.import_yaml:
        cwl_file.construct_workflow_all_params()

    if args.export_to_dir:
        cwl_file.export_to_dir(args.export_to_dir[0])
